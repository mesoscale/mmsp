// cahn-hilliard.hpp
// Algorithms for 2D Cahn-Hilliard model with semi-implicit convex splitting discretization
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"
#include"../energy.hpp"

const double AA = 0.25*(B + Ca + Cb);
const double BB = B*Cm + pow(Ca,2) + pow(Cb,2);
const double CC = 0.5*(A - 3.0*B*pow(Cm,2) - 3.0*pow(Ca,3) - 3.0*pow(Cb,3));
const double DD = A*Cm - B*pow(Cm,3) - pow(Ca,4) - pow(Cb,4);
const double EE = 0.5*(A*pow(Cm,2) - 0.5*B*pow(Cm,4) - 0.5*pow(Ca,5) - 0.5*pow(Cb,5));

const double deltaX = 1.0;
const double CFL = 25.0;
const double dt = pow(deltaX, 4)*CFL/(32.0*D*K);
const double tolerance = 1.0e-8; // "Fair" value. 1e-5 is poor, 1e-12 is publication-quality but slow.

double parametric_energy_density(const double& C)
{
    return AA*pow(C,4) - BB*pow(C,3) - CC*pow(C,2) + DD*C - EE;
}

double expansive_dfdc(const double& C)
{
    return DD - 3.0*BB*pow(C,3) - 2.0*CC*pow(C,2);
}

namespace MMSP {

// Define a Laplacian function for a specific field
template<int dim, typename T>
double field_laplacian(const grid<dim, vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian(0.0);
  vector<int> s = x;

  const double& y = GRID(x)[field];

  for (int i=0; i<dim; i++) {
    s[i] += 1;
    const double& yh = GRID(s)[field];
    s[i] -= 2;
    const double& yl = GRID(s)[field];
    s[i] += 1;

    double weight = 1.0 / pow(dx(GRID, i),2);
    laplacian += weight * (yh - 2.0 * y + yl);
  }
  return laplacian;
}

// Define a Laplacian function missing the central value, for implicit source terms
template<int dim, typename T>
double ring_laplacian(const grid<dim, vector<T> >& GRID, const vector<int>& x, const int field)
{
  T laplacian(0.0);
  vector<int> s = x;

  for (int i=0; i<dim; i++) {
    s[i] += 1;
    const double& yh = GRID(s)[field];
    s[i] -= 2;
    const double& yl = GRID(s)[field];
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

    if (AA<0) {
        if (rank==0)
            std::cerr<<"Error: parameter AA is negative. Double-check your constants."<<std::endl;
            std::exit(-1);
    } else if (BB<0) {
        if (rank==0)
            std::cerr<<"Error: parameter BB is negative. Double-check your constants."<<std::endl;
        std::exit(-1);
    } else if (CC<0) {
        if (rank==0)
            std::cerr<<"Error: parameter CC is negative. Double-check your constants."<<std::endl;
        std::exit(-1);
    } // DD and EE don't matter, since their convexity is not in question.

    const double q[2] = {0.1*std::sqrt(2.0), 0.1*std::sqrt(3.0)};

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(2,0,200,0,200); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			grid(x)[0] = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
			grid(x)[1] = 1.0;
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		output(grid,filename);
		if (rank==0)
			std::cout<<"Timestep is "<<dt<<" (Co="<<CFL<<')'<<std::endl;
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

	// Make sure the grid spacing is correct and isotropic
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		dx(update,d) = deltaX;
		dx(guess,d) = deltaX;
	}

	for (int step=0; step<steps; step++) {

        update.copy(grid); // deep copy: includes data and ghost cells

        double residual=1.0;
        unsigned int iter=0;
        while (residual>tolerance) {
            residual = 0.0;
    		for (int n=0; n<nodes(grid); n++) {
	    		MMSP::vector<int> x = position(grid,n);

	    		const double S1 = grid(x)[0] + dt*D*ring_laplacian(update, x, 1);
	    		const double S2 = DD - 3.0*BB*pow(grid(x)[0],2) - 2.0*CC*grid(x)[0] - K*ring_laplacian(update, x, 0);
	    		const double A12 = (2.0*dim)*D/pow(deltaX,2);
	    		const double A21 = -(2.0*dim)*pow(K,2)/pow(deltaX,2) - 4.0*AA*pow(update(x)[0],2);

	    		// Apply Cramer's Rule to the [2x2][2]=[2] matrix equation
	    		const double cramerDenom = 1.0 - A12*A21;
	    		const double cramerX1 = (S1 - S2*A12)/cramerDenom;
	    		const double cramerX2 = (S2 - A21*S1)/cramerDenom;

                // Copy values to guess. Forgive the mismatched C (0-indexed) and matrix (1-indexed) notations.
                guess(x)[0] = cramerX1;
                guess(x)[1] = cramerX2;

                // Compute ||b-Ax||
                const double Ax1 = cramerX1 + A12*cramerX2;
                const double Ax2 = cramerX2 - cramerX1*A21;
                const double normB = std::sqrt(pow(deltaX,2)*(pow(S1,2)+pow(S2,2)));
                const double normBminusAx = std::sqrt(pow(deltaX,2)*pow(S1-Ax1,2) + pow(deltaX,2)*pow(S2-Ax2,2));

                residual += normBminusAx/(normB*nodes(grid));
	    	}
       		#ifdef MPI_VERSION
	        double localResidual=residual;
	        MPI::COMM_WORLD.Barrier();
       		MPI::COMM_WORLD.Allreduce(&localResidual, &residual, 1, MPI_DOUBLE, MPI_SUM);
       		#endif

            swap(update, guess);
            ghostswap(update);

	    	iter++;
        }

		swap(grid,update);
		ghostswap(grid);

  		#ifdef DEBUG
   		if (rank==0)
   		    std::cout<<"Step "<<step<<" converged with residual "<<residual<<" after "<<iter<<" iterations."std::endl;
   		#else
		double energy = 0.0;
		double mass = 0.0;
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			energy += dx(grid)*dy(grid)*parametric_energy_density(update(x)[0]);
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
