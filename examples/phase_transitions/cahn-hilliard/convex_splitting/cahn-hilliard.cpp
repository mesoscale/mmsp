/** cahn-hilliard.hpp
 ** Algorithms for 2D Cahn-Hilliard model with semi-implicit discretization by convex splitting.
 ** Questions/comments to trevor.keller@gmail.com (Trevor Keller)
 **
 ** This code sets parameters for the evolution of two different energy configurations.
 **
 ** "Vanilla" solves the Cahn-Hilliard equations for the most-simple free energy density
 ** producing spinodal decomposition, F(C)=0.25*C^4-0.5*C^2. To build this case, edit the Makefile
 ** to specify -DVANILLA in the flags parameter.
 **
 ** Otherwise, the Cahn-Hilliard equation is solved for a simple binary phase diagram with parametrized
 ** free energy density, F(C)=-0.5*A*(C-C0)^2 + 0.25*B*(C-C0)^4 + 0.25*Ca*(C-Ca)^4 + 0.25*Cb*(C-Cb)^4.
 ** This case demonstrates stiffer numerics, so proceed with caution. To build, remove the -DVANILLA flag.
 **/

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include<fstream>
#include"cahn-hilliard.hpp"

#ifdef VANILLA
const double D = 1.0; // diffusivity
const double K = 1.0; // gradient energy coefficient
const double C0 = 0.5;// system composition
#else
const double Ca = 0.05; // alpha-phase solvus composition
const double Cb = 0.95; // beta-phase solvus composition
const double C0 = 0.45; // mean composition
const double D = 2.0/(Cb-Ca);
const double K = 2.0;
#endif

// Spatial constants
const int edge = 128;
const double deltaX = 32.0/double(edge-1); //~0.25 if edge=128;
const double dt = 0.25/pow(D*K,4.0);
const double CFL = 32.0*dt*D*K/pow(deltaX,4.0); // to judge improvement w.r.t. explicit discretization
// Max. semi-implicit timestep should be timestep = pow(deltaX,2.0)/(32.0*D*K)

// Numerical constants
const double tolerance = 1.0e-12;     // Choose wisely. 1e-5 is poor, 1e-8 fair: mass is not conserved. 1e-12 is good, mass is conserved. YMMV depending on coefficients.
const unsigned int max_iter = 10000; // don't let the solver stagnate
const double omega = 1.0;

#ifdef VANILLA
// Vanilla energetics
#include<time.h>
template<typename T> inline
T energy_density(const T& C)
{
    return 0.25*pow(C,4.0) - 0.5*pow(C,2.0);
}
template<typename T> inline
T full_dfdc(const T& C)
{
    return pow(C,3.0) - C;
}
template<typename T> inline
T contractive_dfdc(const T& C)
{
    return pow(C,3.0);
}
template<typename T> inline
T nonlinear_coeff(const T& C)
{
    return pow(C,2.0);
}
template<typename T> inline
T expansive_dfdc(const T& C)
{
    return -C;
}
#else
// Parametric energetics
const double A = 2.0;
const double B = A/pow(Ca-C0,2.0); // = 9.8765

const double AA = B + Ca + Cb;
const double BB = 3.0*(B*C0 + pow(Ca,2.0) + pow(Cb,2.0));
const double CC = 3.0*(B*pow(C0,2.0) + pow(Ca,3.0) + pow(Cb,3.0)) - A;
const double DD = B*pow(C0,3.0) + pow(Ca,4.0) + pow(Cb,4.0) - A*C0;

template<typename T> inline
T energy_density(const T& C)
{
    return -0.5*A*pow(C-C0,2.0) + 0.25*B*pow(C-C0,4.0) + 0.25*Ca*pow(C-Ca,4.0) + 0.25*Cb*pow(C-Cb,4.0);
}
template<typename T> inline
T full_dfdc(const T& C)
{
    return AA*pow(C,3.0) - BB*pow(C,2.0) + CC*C - DD;
}
template<typename T> inline
T contractive_dfdc(const T& C)
{
    return AA*pow(C,3.0) - BB*pow(C,2.0);
}
template<typename T> inline
T nonlinear_coeff(const T& C)
{
    return AA*pow(C,2.0) - BB*C;
}
template<typename T> inline
T expansive_dfdc(const T& C)
{
    return CC*C - DD;
}
#endif


namespace MMSP {

// Define a Laplacian function for a specific field
template<int dim, typename T>
T field_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian = 0.0;
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

// Define a Laplacian function missing the central value, for implicit source terms
template<int dim, typename T>
T fringe_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian(0.0);
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

void generate(int dim, const char* filename)
{
    // Initial conditions after CHiMaD Phase Field Workshop benchmark problem (October 2015)

	if (dim!=2) {
		std::cerr<<"ERROR: Convex splitting discretization is only 2-D, for now."<<std::endl;
		std::exit(-1);
	}

	int rank=0;
	int np=1;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	np = MPI::COMM_WORLD.Get_size();
	#endif

    srand(time(NULL)/(rank+2));

    #ifdef DEBUG
    // Validate free energy functions, splitting, and linearization
    for (int i=0; i<10000/np; i++) {
        double c = 0.5 + 0.6*((double(rand())/RAND_MAX) - 0.5);
        if (fabs(full_dfdc(c) - (contractive_dfdc(c) + expansive_dfdc(c))) > tolerance) {
            std::cerr<<"ERROR: full derivative ft("<<c<<") ≠ fc + fe. Check functions."<<std::endl;
            std::exit(-1);
        }
        if (fabs(contractive_dfdc(c) - nonlinear_coeff(c)*c) > tolerance) {
            std::cerr<<"ERROR: contractive derivative fc("<<c<<") ≠ c*l(c). Check linearization."<<std::endl;
            std::exit(-1);
        }
    }
    #endif

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(2,0,edge,0,edge); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		#ifdef VANILLA
    	for (int n=0; n<nodes(grid); n++)
		    grid(n)[0] = C0*(1.0 + 0.1 * double(rand())/RAND_MAX);
		#else
        const double q[2] = {0.1*sqrt(2.0), 0.1*sqrt(3.0)}; // produces stipes oriented 45 degrees to horizontal
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			double wave = x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1];
			grid(n)[0] = C0*(1.0 + 0.1 * std::cos(wave));
		}
        #endif

		ghostswap(grid); // otherwise, parallel jobs have a "window frame" artifact

		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
		    grid(n)[1] = full_dfdc(grid(n)[0]) - K*field_laplacian(grid, x, 0);
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		output(grid,filename);
		if (rank==0)
			std::cout<<"Timestep is "<<dt<<" (Co="<<CFL<<")."<<std::endl;

        double dV = 1.0;
        for (int d=0; d<dim; d++)
            dV *= dx(grid,d);

		double energy = 0.0;
		double mass = 0.0;
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			MMSP::vector<MMSP::vector<double> > gradC = grad(grid, x);
            double magSqGradC = 0.0;
            for (int d=0; d<dim; d++)
                magSqGradC += pow(gradC[d][0],2.0);
			const double C = grid(x)[0];
			energy += dV*(energy_density(C) + 0.5*K*magSqGradC);
			mass += dV*C;
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		double localEnergy = energy;
		double localMass = mass;
		MPI::COMM_WORLD.Reduce(&localEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0);
		MPI::COMM_WORLD.Reduce(&localMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0);
		#endif
		if (rank==0)
			std::cout<<'0'<<'\t'<<energy<<'\t'<<mass<<std::endl;

        #ifdef DEBUG
        std::ofstream ferr;
        if (rank==0) {
            ferr.open("error.log");
    	    ferr<<"step\titer\tnormBminusAX\tnormB\tresidual\n";
    	    ferr.close();
    	}
    	#endif
	}
}

template<int dim, typename T>
void update(MMSP::grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    ghostswap(oldGrid);

	MMSP::grid<dim,vector<T> > newGrid(oldGrid);   // new values at each point and initial guess for iteration

    newGrid.copy(oldGrid); // deep copy: includes data and ghost cells. Expensive.

	// Make sure the grid spacing is correct. Modify at will.
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
	}

    double gridSize = static_cast<double>(nodes(oldGrid));
	#ifdef MPI_VERSION
    double localGridSize = gridSize;
    MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

    double lapWeight = 0.0;
    double dV = 1.0;
    for (int d=0; d<dim; d++) {
        lapWeight += 2.0/pow(dx(oldGrid,d),2.0); // dim=2 -> lapWeight = 4/h^2 if dy=dx=h
        dV *= dx(oldGrid,d);
    }

    #ifdef DEBUG
    std::ofstream ferr;
    if (rank==0)
        ferr.open("error.log", std::ofstream::app);
    #endif

	for (int step=0; step<steps; step++) {
        double residual=1.0;
        unsigned int iter=0;

        while (iter<max_iter && residual>tolerance) {

            /*  ==== RED-BLACK GAUSS SEIDEL ====
                Iterate over a checkerboard, updating first red then black tiles.
                This method eliminates the third "guess" grid, and should converge faster.
                In 2-D and 3-D, if the sum of indices is even, then the tile is Red; else, Black.
            */

    		for (int color=1; color>-1; color--) {
    		    // If color==1, skip BLACK tiles, which have Σx[d] odd
    		    // If color==0, skip RED tiles, which have Σx[d] even
        		for (int n=0; n<nodes(oldGrid); n++) {
	    		    MMSP::vector<int> x = position(oldGrid,n);
	    	    	int x_sum=0;
	        		for (int d=0; d<dim; d++)
	        		    x_sum += x[d];
    	    		if (x_sum%2 == color)
	    		        continue;

	    	    	const T cOld = oldGrid(n)[0];
	        		const T cLast = newGrid(n)[0];

	        		// A is defined by the last guess, stored in newGrid(x). It is a 2x2 matrix.
    	    		//const double A11 = 1.0;
    	    		const double A12 = dt*D*lapWeight;
	    		    const double A21 = -nonlinear_coeff(cLast) - K*lapWeight;
	    		    //const double A22 = 1.0;

	    	    	// B is defined by the last value, stored in oldGrid(x), and the last guess, stored in newGrid(x). It is a 2x1 column.
	        		const T lapC = fringe_laplacian(newGrid, x, 0); // excludes central term
	        		const T lapU = fringe_laplacian(newGrid, x, 1); // excludes central term

    	    		const T B1 = cOld + D*dt*lapU;
	    		    const T B2 = expansive_dfdc(cOld) - K*lapC;

                    // Solve the iteration system AX=B using Cramer's rule
	        		const double detA = 1.0 - A12*A21; // det|A|
	        		const T cNew = (B1 - B2*A12)/detA; // X1
    	    		const T uNew = (B2 - B1*A21)/detA; // X2

                    // (Don't) Apply relaxation
	    	    	const T uOld = oldGrid(n)[1];
                    newGrid(n)[0] = (1.0 - omega)*cOld + omega*cNew;
                    newGrid(n)[1] = (1.0 - omega)*uOld + omega*uNew;
                }
                ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial
            }

            /*  ==== RESIDUAL ====
                The residual is computed from the original matrix form, Ax=b:
                any Old term goes into B, while any New term goes in AX. Note that
                this is not the iteration matrix, it is the original system of equations.
            */

            if (iter<5 || (iter<100 && iter%10==0) || iter%100==0) {
                double normB = 0.0;
                residual = 0.0;
       	        for (int n=0; n<nodes(oldGrid); n++) {
        	    	MMSP::vector<int> x = position(oldGrid,n);
            	    const T lapC = field_laplacian(newGrid, x, 0);
    	       	    const T lapU = field_laplacian(newGrid, x, 1);

           	     	const T cOld = oldGrid(n)[0];
                    const T cNew = newGrid(n)[0];
                    const T uNew = newGrid(n)[1];

        	  	    const T B1 = cOld;
    	   	        const T B2 = expansive_dfdc(cOld);

                    const T AX1 = cNew - D*dt*lapU;
                    const T AX2 = uNew - contractive_dfdc(cNew) + K*lapC;

                    // Compute the Error from parts of the solution
                    const T R1 = AX1 - B1;
                    const T R2 = AX2 - B2;

                    const double error = pow(R1,2.0) + pow(R2,2.0);
                    residual += error;
                    normB += pow(B1,2.0) + pow(B2,2.0);
            	}

   	        	#ifdef MPI_VERSION
   	            double localResidual=residual;
           	    MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
           	    double localNormB=normB;
           		MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
   	            #endif

   	            #ifdef DEBUG
           	    if (rank==0)
   	        	    ferr<<step<<'\t'<<iter<<'\t'<<residual<<'\t'<<normB<<'\t';
   	            #endif

   	            residual = sqrt(residual/normB)/(2.0*gridSize);

           	    #ifdef DEBUG
   	        	if (rank==0)
  	                ferr<<residual<<std::endl;
                #endif
            }

	    	iter++;
        }

        if (iter==max_iter) {
            if (rank==0)
                std::cerr<<"    Solver stagnated on step "<<step<<": "<<iter<<" iterations with residual="<<residual<<std::endl;
            std::exit(-1);
        } else if (!std::isfinite(residual) || std::isnan(residual)) {
            if (rank==0)
                std::cerr<<"    Numerical error on step "<<step<<" after "<<iter<<" iterations"<<std::endl;
            std::exit(-1);
        }

		double energy = 0.0;
		double mass = 0.0;
		for (int n=0; n<nodes(newGrid); n++) {
			MMSP::vector<int> x = position(newGrid,n);
			MMSP::vector<MMSP::vector<T> > gradC = grad(newGrid, x);
            double gradCsq = 0.0;
            for (int d=0; d<dim; d++)
                gradCsq += pow(gradC[d][0],2.0);
			const T C = newGrid(x)[0];
			energy += dV*(energy_density(C) + 0.5*K*gradCsq);
			mass += dV*C;
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		double localEnergy = energy;
		double localMass = mass;
		MPI::COMM_WORLD.Reduce(&localEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0);
		MPI::COMM_WORLD.Reduce(&localMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0);
		#endif
		if (rank==0)
			std::cout<<iter<<'\t'<<energy<<'\t'<<mass<<std::endl;

    	swap(oldGrid,newGrid);
    	ghostswap(oldGrid);

	}
   	#ifdef DEBUG
	ferr.close();
	#endif
	if (rank==0)
		std::cout<<std::flush;
}

} // MMSP
#endif

#include"MMSP.main.hpp"
