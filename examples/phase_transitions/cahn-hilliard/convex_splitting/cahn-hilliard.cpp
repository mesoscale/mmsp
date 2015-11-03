// cahn-hilliard.hpp
// Algorithms for 2D Cahn-Hilliard model with semi-implicit convex splitting discretization
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"
#include"../energy.hpp"

const int edge = 200;
const double deltaX = 1.0;
const double CFL = 20.0;
const double dt = pow(deltaX, 4)*CFL/(32.0*D*K);
const double tolerance = 1.0e-8; // "Fair" value. 1e-5 is poor, 1e-12 is publication-quality but slow.
const unsigned int max_iter = 5000; // don't let the solver stagnate
const int resfreq=2; // number of iterations per residual computation


// Vanilla Cahn-Hilliard
const double AA = 1.0;
const double BB = 0.0;
const double CC = -1.0;
const double DD = 0.0;

template <typename T>
T energy_density(const T& C)
{
    return 0.25*pow(C,4) - 0.5*pow(C,2);
}
template <typename T>
T full_dfdc(const T& C)
{
    return pow(C,3) - C;
}
template <typename T>
T contractive_dfdc(const T& C)
{
    return pow(C,3);
}
template <typename T>
T expansive_dfdc(const T& C)
{
    return -C;
}

/*
// Parametric Cahn-Hilliard
const double AA = B + Ca + Cb;
const double BB = 3.0*(B*Cm + pow(Ca,2) + pow(Cb,2));
const double CC = 3.0*(B*pow(Cm,2) + pow(Ca,3) + pow(Cb,3)) - A;
const double DD = B*pow(Cm,3) + pow(Ca,4) + pow(Cb,4) - A*Cm;

template <typename T>
T energy_density(const T& C)
{
    return -0.5*A*pow(C-Cm,2) + 0.25*B*pow(C-Cm,4) + 0.25*Ca*pow(C-Ca,4) + 0.25*Cb*pow(C-Cb,4);
}
template <typename T>
T full_dfdc(const T& C)
{
    return AA*pow(C,3) - BB*pow(C,2) + CC*C - DD;
}
template <typename T>
T contractive_dfdc(const T& C)
{
    return AA*pow(C,3) - BB*pow(C,2);
}
template <typename T>
T expansive_dfdc(const T& C)
{
    return CC*C - DD;
}
*/

namespace MMSP {

// Define a Laplacian function for a specific field
template<int dim, typename T>
T field_laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian(0.0);
  vector<int> s = x;

  const T& y = GRID(x)[field];

  for (int i=0; i<dim; i++) {
    s[i] += 1;
    const T& yh = GRID(s)[field];
    s[i] -= 2;
    const T& yl = GRID(s)[field];
    s[i] += 1;

    double weight = 1.0 / pow(dx(GRID, i),2);
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

    double weight = 1.0 / pow(dx(GRID, i),2);
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
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    const double q[2] = {0.1*sqrt(2.0), 0.1*sqrt(3.0)};

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(2,0,edge,0,edge); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			grid(x)[0] = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
		}

		ghostswap(grid); // otherwise, parallel jobs have a "window frame" artifact

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
		    grid(x)[1] = full_dfdc(grid(x)[0]) - K*field_laplacian(grid, x, 0);
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		output(grid,filename);
		if (rank==0)
			std::cout<<"Timestep is "<<dt<<" (Co="<<CFL<<")\nIters\tEnergy\tMass"<<std::endl;
	}
}

template <int dim, typename T>
void update(MMSP::grid<dim,vector<T> >& grid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    ghostswap(grid);

	MMSP::grid<dim,vector<T> > update(grid); // new values at each point and initial guess for iteration
	MMSP::grid<dim,vector<T> > guess(grid); // new guess for iteration

	// Make sure the grid spacing is correct. Modify at will.
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		dx(update,d) = deltaX;
		dx(guess,d) = deltaX;
	}

    double dV = 1.0;
    for (int d=0; d<dim; d++)
        dV *= dx(grid,d);
    double gridSize=1.0*nodes(grid);
	#ifdef MPI_VERSION
    double localGridSize=gridSize;
    MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&localGridSize, &gridSize, 1, MPI_DOUBLE, MPI_SUM);
	#endif

    double lapWeight = 0.0;
    for (int d=0; d<dim; d++)
        lapWeight += 2.0/pow(dx(grid,d),2);

    update.copy(grid); // deep copy: includes data and ghost cells

	for (int step=0; step<steps; step++) {

        double residual=1.0;
        unsigned int iter=0;

        while (iter<max_iter && residual>tolerance) {
            // update stores the old gues, Jacobi lexicographic index l
            // guess  stores the new gues, Jacobi lexicographic index l+1

            // Solve AX=B using Cramer's rule
    		for (int n=0; n<nodes(grid); n++) {
	    		MMSP::vector<int> x = position(grid,n);

	    		const double B1 = grid(x)[0] + dt*D*fringe_laplacian(update, x, 1);
	    		const double B2 = expansive_dfdc(grid(x)[0]) - K*fringe_laplacian(update, x, 0);
	    		const double A12 = dt*D*lapWeight;
	    		const double A21 = -contractive_dfdc(update(x)[0]) - K*lapWeight;

	    		const double denom = 1.0 - A12*A21;
	    		guess(x)[0] =(B1 -  B2*A12)/denom;
	    		guess(x)[1] =(B2 - A21*B1 )/denom;
             }

            swap(update, guess);
            ghostswap(update);

	    	iter++;

    		// Strictly, ||b-Ax|| should be re-computed using the latest guess after each iteration.
    		// But this almost doubles the computational cost, so I'll cheat and infrequently compute it.
            double normB = 0.0;
    		if (iter%resfreq==0) {
                residual = 0.0;
    		    for (int n=0; n<nodes(grid); n++) {
	    	    	MMSP::vector<int> x = position(grid,n);

                    // Compute Σ(dV||B-AX||) / Σ(dV||B||) using L2 vector norms
                    const double R1 = guess(x)[0] - grid(x)[0] - dt*D*field_laplacian(guess, x, 1);
                    const double R2 = guess(x)[1] - contractive_dfdc(guess(x)[0]) - expansive_dfdc(grid(x)[0]) + K*field_laplacian(guess, x, 0);

                    const double normBminusAX = pow(R1,2) + pow(R2,2);

                    residual += normBminusAX;
                    normB += pow(grid(x)[0],2) + pow(expansive_dfdc(grid(x)[0]),2);
    	    	}

    	    	#ifdef MPI_VERSION
    	    	double localResidual=residual;
        		MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
    	    	double localNormB=normB;
        		MPI::COMM_WORLD.Allreduce(&localNormB, &normB, 1, MPI_DOUBLE, MPI_SUM);
    	    	#endif

    	    	residual = sqrt(residual)/(2.0*gridSize*sqrt(normB));
	    	}

	    	#ifdef DEBUG
       		if (rank==0)
   	    	    std::cout<<step<<'\t'<<iter<<'\t'<<normB<<'\t'<<residual<<std::endl;
   		    #endif
        }

		swap(grid,guess);
		ghostswap(grid);

  		#ifndef DEBUG
		double energy = 0.0;
		double mass = 0.0;
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			energy += dx(grid)*dy(grid)*energy_density(update(x)[0]);
			mass += dx(grid)*dy(grid)*update(x)[0];
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		double localEnergy = energy;
		double localMass = mass;
		MPI::COMM_WORLD.Reduce(&localEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0);
		MPI::COMM_WORLD.Reduce(&localMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0);
		#endif
		if (rank==0)
			std::cout<<iter<<'\t'<<energy<<'\t'<<mass<<'\n';
		#endif

	}
	#ifndef DEBUG
	if (rank==0)
		std::cout<<std::flush;
	#endif
}

} // MMSP
#endif

#include"MMSP.main.hpp"
