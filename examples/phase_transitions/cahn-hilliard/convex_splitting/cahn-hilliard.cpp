/** cahn-hilliard.cpp
 ** Algorithms for 2D Cahn-Hilliard model with convex splitting.
 **
 ** This is an advanced example, with precompiler directives switching between
 ** shared-memory (OpenMP) and message-passing (MPI) parallelism, with the equations
 ** of motion discretized using semi-implicit convex splitting for guaranteed energy
 ** stability, i.e. dE/dt<=0, and with successive over-relaxation to accelerate
 ** convergence.
 **
 ** Questions/comments to trevor.keller@gmail.com (Trevor Keller)
 **/

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include<cfloat>
#include"cahn-hilliard.hpp"

/** This code sets parameters for the evolution of two different energy configurations.
 **
 ** "Vanilla" solves the Cahn-Hilliard equations for the most-simple free energy density
 ** producing spinodal decomposition, F(C)=0.25*C^4-0.5*C^2. To build this case, edit the Makefile
 ** to specify -DVANILLA in the flags parameter, e.g.
 ** 	flags = -O3 -DVANILLA -Wall ...
 **
 ** Otherwise, the Cahn-Hilliard equation is solved for a simple binary phase diagram with parametrized
 ** free energy density, F(C)=-0.5*A*(C-C0)^2 + 0.25*B*(C-C0)^4 + 0.25*Ca*(C-Ca)^4 + 0.25*Cb*(C-Cb)^4.
 ** This case demonstrates stiffer numerics, so proceed with caution. To build, drop -DVANILLA, e.g.
 ** 	flags = -O3 -Wall ...
 **/

#ifdef VANILLA
int edge = 32;
double deltaX = 1.0;
double C0 = 0.0; // system composition
double D  = 1.0; // diffusivity
double K  = 1.0; // gradient energy coefficient
#else
int edge = 200;
double deltaX = 1.0;
double C0 = 0.45; // mean composition
double Ca = 0.05; // alpha-phase solvus composition
double Cb = 0.95; // beta-phase solvus composition
double D  = 2.0/(Cb-Ca);
double K  = 2.0;
#endif

// Max. semi-implicit timestep should be timestep = pow(deltaX,2.0)/(32.0*D*K)
double dt = pow(deltaX,2.0)/(D*K);			// timestep
double CFL = 32.0*dt*D*K/pow(deltaX,4.0);	// to judge improvement w.r.t. explicit discretization

// Numerical constants
double tolerance = 1.0e-14;		// Choose wisely. 1e-10 is the minimum toloerance for which mass is conserved.
								// Tighter tolerance is better, but increases runtime.
unsigned int residual_step = 5;	// number of iterations between residual computations
unsigned int max_iter = 10000;	// don't let the solver stagnate

#ifdef VANILLA
double omega = 1.0;				// relaxation parameter for SOR. omega=1 is stock Gauss-Seidel.
								//                               omega<1 is under-relaxation (slower convergence, not recommended),
								//                               omega>1 is over-relaxation (faster convergence), stagnates above omega=1.1.

// Vanilla energetics
#include<time.h>
template<typename T> inline
double energy_density(const T& C)
{
	// Minima at C = +/- 1, local maximum at C=0
	return 0.25*pow(C,4.0) - 0.5*pow(C,2.0);
}
template<typename T> inline
double full_dfdc(const T& C)
{
	return pow(C,3.0) - C;
}
template<typename T> inline
double contractive_dfdc(const T& C)
{
	return pow(C,3.0);
}
template<typename T> inline
double nonlinear_coeff(const T& C)
{
	return pow(C,2.0);
}
template<typename T> inline
double expansive_dfdc(const T& C)
{
	return -C;
}
#else
double omega = 1.2;				// parametric equations are more forgiving, stagnates above omega=1.67.

// Parametric energetics
double A = 2.0;
double B = A/pow(Ca-C0,2.0); // = 9.8765

// Define the set of quartic coefficients Qi
double QA = B + Ca + Cb;
double QB = 3.0*(B*C0 + Ca*Ca + Cb*Cb);
double QC = 3.0*(B*C0*C0 + pow(Ca,3.0) + pow(Cb,3.0));
double QD = B*pow(C0,3.0) + pow(Ca,4.0) + pow(Cb,4.0);

template<typename T> inline
double energy_density(const T& C)
{
	// Minima at C = 0 and 1
	return -0.5*A*pow(C-C0,2.0) + 0.25*B*pow(C-C0,4.0) + 0.25*Ca*pow(C-Ca,4.0) + 0.25*Cb*pow(C-Cb,4.0);
}
template<typename T> inline
double full_dfdc(const T& C)
{
	return -A*(C-C0) + B*pow(C-C0,3.0) + Ca*pow(C-Ca,3.0) + Cb*pow(C-Cb,3.0);
}
template<typename T> inline
double contractive_dfdc(const T& C)
{
	return QA*pow(C,3.0) + QC*C;
}
template<typename T> inline
double nonlinear_coeff(const T& C)
{
	return QA*pow(C,2.0) + QC;
}
template<typename T> inline
double expansive_dfdc(const T& C)
{
	return -A*(C-C0) - QB*pow(C,2.0) - QD;
}
#endif


namespace MMSP
{

// Discrete Laplacian operator for a specific field
template<int dim, typename T>
double field_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
	double laplacian = 0.0;
	vector<int> s = x;

	const T& y = GRID(x)[field];

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s)[field];
		s[i] -= 2;
		const T& yl = GRID(s)[field];
		s[i] += 1;

		double weight = 1.0 / pow(dx(GRID, i),2.0);
		laplacian += weight * (yh - 2.0 * y + yl);
	}
	return laplacian;
}

// Discrete Laplacian operator missing the central value, for implicit source terms
template<int dim, typename T>
double fringe_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
	double laplacian = 0.0;
	vector<int> s = x;

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s)[field];
		s[i] -= 2;
		const T& yl = GRID(s)[field];
		s[i] += 1;

		double weight = 1.0 / pow(dx(GRID, i),2.0);
		laplacian += weight * (yh + yl);
	}
	return laplacian;
}

template<int dim, typename T> void reportEnergy(const MMSP::grid<dim,vector<T> >& GRID)
{
	int rank=0;
	#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	double dV = 1.0;
	for (int d=0; d<dim; d++)
		dV *= dx(GRID,d);

	if (rank==0)
		std::cout<<"Timestep is "<<dt<<" (explicit Co="<<CFL<<")."<<std::endl;

	double energy = 0.0;
	double mass = 0.0;
	for (int n=0; n<nodes(GRID); n++) {
		vector<int> x = position(GRID,n);
		vector<vector<T> > gradC = grad(GRID, x);
		double magSqGradC = 0.0;
		for (int d=0; d<dim; d++)
			magSqGradC += pow(gradC[d][0],2.0);
		double C = GRID(n)[0];
		energy += (energy_density(C) + 0.5*K*magSqGradC);
		mass += C;
	}
	#ifdef MPI_VERSION
	double localEnergy = energy;
	double localMass = mass;
	MPI_Reduce(&localEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&localMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	#endif
	//if (rank==0)
	//	std::cout<<'0'<<'\t'<<dV*energy<<'\t'<<dV*mass<<std::endl;
}

void generate(int dim, const char* filename)
{
	if (dim==1) {
		edge=1024;
		GRID1D initGrid(2,0,edge); // field 0 is c, field 1 is mu
		for (int n=0; n<nodes(initGrid); n++)
			initGrid(n)[0] = C0 + 0.25*(double(rand())/RAND_MAX) - 0.125; // produces noise with amplitude 0.125 about C=C0

		ghostswap(initGrid); // otherwise, parallel jobs have a "window frame" artifact

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid,n);
			initGrid(n)[1] = full_dfdc(initGrid(n)[0]) - K*field_laplacian(initGrid, x, 0);
		}

		output(initGrid,filename);
		//reportEnergy(initGrid);
	} else if (dim==2) {
		GRID2D initGrid(2,0,2*edge,0,edge); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(initGrid,d) = deltaX;

		#ifdef VANILLA
		for (int n=0; n<nodes(initGrid); n++)
			initGrid(n)[0] = C0 + 0.25*(double(rand())/RAND_MAX) - 0.125; // produces noise with amplitude 0.125 about C=C0
		#else
		// Initial conditions after CHiMaD Phase Field Workshop benchmark problem (October 2015)
		double q[2] = {0.1*sqrt(2.0), 0.1*sqrt(3.0)}; // {0.1*sqrt(2.0), 0.1*sqrt(3.0)} produces stripes oriented 45 degrees to horizontal
		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid,n);
			double wave = x[0]*dx(initGrid,0)*q[0] + x[1]*dx(initGrid,1)*q[1];
			initGrid(n)[0] = C0*(1.0 + 0.1 * std::cos(wave));
		}
		#endif

		ghostswap(initGrid); // otherwise, parallel jobs have a "window frame" artifact

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid,n);
			initGrid(n)[1] = full_dfdc(initGrid(n)[0]) - K*field_laplacian(initGrid, x, 0);
		}

		output(initGrid,filename);
		//reportEnergy(initGrid);
	} else if (dim==3) {
		GRID3D initGrid(2,0,2*edge,0,edge,0,edge/4); // field 0 is c, field 1 is mu
		for (int n=0; n<nodes(initGrid); n++)
			initGrid(n)[0] = C0 + 0.25*(double(rand())/RAND_MAX) - 0.125; // produces noise with amplitude 0.125 about C=C0

		ghostswap(initGrid); // otherwise, parallel jobs have a "window frame" artifact

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid,n);
			initGrid(n)[1] = full_dfdc(initGrid(n)[0]) - K*field_laplacian(initGrid, x, 0);
		}

		output(initGrid,filename);
		//reportEnergy(initGrid);
	}
}

template<int dim, typename T>
void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);   // new values at each point and initial guess for iteration

	newGrid.copy(oldGrid); // deep copy: includes data and ghost cells. Expensive.

	// Make sure the grid spacing is correct.
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
	}

	double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
	double localGridSize = gridSize;
	MPI_Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif

	double lapWeight = 0.0;
	double dV = 1.0;
	for (int d=0; d<dim; d++) {
		lapWeight += 2.0/pow(dx(oldGrid,d),2.0); // dim=2 -> lapWeight = 4/h^2 if dy=dx=h
		dV *= dx(oldGrid,d);
	}

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		double residual=1.0;
		unsigned int iter=0;

		ghostswap(oldGrid);
		ghostswap(newGrid);

		while (iter<max_iter && residual>tolerance) {

			/*  ==== RED-BLACK GAUSS SEIDEL ====
			    Iterate over a checkerboard, updating first red then black tiles.
			    This method eliminates the third "guess" grid, and should converge faster.
			    In 2-D and 3-D, if the sum of indices is even, then the tile is Red; else, Black.
			*/

			for (int color=1; color>-1; color--) {
				// If color==1, skip BLACK tiles, which have Σx[d] odd
				// If color==0, skip RED tiles, which have Σx[d] even
				#ifndef MPI_VERSION
				// OpenMP parallel loop over nodes
				#pragma omp parallel for schedule(dynamic)
				#endif
				for (int n=0; n<nodes(oldGrid); n++) {
					vector<int> x = position(oldGrid,n);
					int x_sum=0;
					for (int d=0; d<dim; d++)
						x_sum += x[d];
					if (x_sum%2 == color)
						continue;

					T cOld = oldGrid(n)[0];
					T cGuess = newGrid(n)[0]; // value from last "guess" iteration
					T uGuess = newGrid(n)[1];

					// A is defined by the last guess, stored in newGrid(n). It is a 2x2 matrix.
					//T A11 = 1.0;
					T A12 = dt*D*lapWeight;
					T A21 = -nonlinear_coeff(newGrid(n)[0]) - K*lapWeight;
					//T A22 = 1.0;

					T detA = 1.0 - (A12*A21); // determinant of A

					// B is defined by the last value, stored in oldGrid(n), and the last guess, stored in newGrid(n). It is a 2x1 column.
					T B1 = cOld + D*dt*fringe_laplacian(newGrid, x, 1);
					T B2 = expansive_dfdc(cOld) - K*fringe_laplacian(newGrid, x, 0);

					// Solve the iteration system AX=B using Cramer's rule
					T cNew = (B1 - B2*A12)/detA; // X1
					T uNew = (B2 - B1*A21)/detA; // X2

					// (Don't) Apply relaxation
					newGrid(n)[0] = omega*cNew + (1.0 - omega)*cGuess;
					newGrid(n)[1] = omega*uNew + (1.0 - omega)*uGuess;

				}
				ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial

			}

			iter++;

			/*  ==== RESIDUAL ====
			    The residual is computed from the original matrix form, Ax=b:
			    any Old term goes into B, while any New term goes in AX. Note that
			    this is not the iteration matrix, it is the original system of equations.
			*/

			if (iter<residual_step || iter%residual_step==0) {
				double normB = 0.0;
				residual = 0.0;
				#ifndef MPI_VERSION
				// OpenMP parallel loop over nodes
				#pragma omp parallel for schedule(dynamic)
				#endif
				for (int n=0; n<nodes(oldGrid); n++) {
					vector<int> x = position(oldGrid,n);
					T lapC = field_laplacian(newGrid, x, 0);
					T lapU = field_laplacian(newGrid, x, 1);

					T cOld = oldGrid(n)[0];
					T cNew = newGrid(n)[0];
					T uNew = newGrid(n)[1];

					T B1 = cOld;
					T B2 = expansive_dfdc(cOld);

					T AX1 = cNew - D*dt*lapU;
					T AX2 = uNew - contractive_dfdc(cNew) + K*lapC;

					// Compute the Error from parts of the solution
					double R1 = B1 - AX1;
					double R2 = B2 - AX2;

					double error = R1*R1 + R2*R2;
					#ifndef MPI_VERSION
					#pragma omp critical
					{
					#endif
						residual += error;
						normB += B1*B1 + B2*B2;
					#ifndef MPI_VERSION
					}
					#endif
				}

				#ifdef MPI_VERSION
				double localResidual=residual;
				MPI_Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				double localNormB=normB;
				MPI_Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				#endif

				residual = sqrt(residual/normB)/(2.0*gridSize);
			}

		}

		if (iter==max_iter) {
			if (rank==0)
				std::cerr<<"    Solver stagnated on step "<<step<<": "<<iter<<" iterations with residual="<<residual<<std::endl;

            MMSP::Abort(-1);
		}

		/*
		double energy = 0.0;
		double mass = 0.0;
		for (int n=0; n<nodes(newGrid); n++) {
			vector<int> x = position(newGrid,n);
			vector<vector<T> > gradC = grad(newGrid, x);
			double gradCsq = 0.0;
			for (int d=0; d<dim; d++)
				gradCsq += pow(gradC[d][0],2.0);
			double C = newGrid(n)[0];
			energy += (energy_density(C) + 0.5*K*gradCsq);
			mass += C;
		}
		#ifdef MPI_VERSION
		double localEnergy = energy;
		double localMass = mass;
		MPI_Reduce(&localEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&localMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		#endif
		if (rank==0)
			std::cout<<iter<<'\t'<<dV*energy<<'\t'<<dV*mass<<std::endl;
		*/

		swap(oldGrid,newGrid);
	}
	if (rank==0)
		std::cout<<std::flush;
}

} // MMSP
#endif

#include"MMSP.main.hpp"
