// cahn-hilliard.hpp
// Algorithms for 2D Cahn-Hilliard model with semi-implicit convex splitting discretization
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include<fstream>
#include"cahn-hilliard.hpp"

#ifdef VANILLA
// Vanilla Cahn-Hilliard
#include<time.h>
const double D = 1.0;
const double K = 1.0;

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
// Parametric Cahn-Hilliard
const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca + Cb);      // = 0.5
const double A = 2.0;
const double B = A/pow(Ca-Cm,2.0); // = 9.8765
const double D = 2.0/(Cb-Ca);         // = 2.2222
const double K = 2.0;

template<typename T> inline
double energydensity(const T& C)
{
        return -0.5*A*pow(C-Cm,2.0) + 0.25*B*pow(C-Cm,4.0) + 0.25*Ca*pow(C-Ca,4.0) + 0.25*Cb*pow(C-Cb,4.0);
}

template<typename T> inline
double dfdc(const T& C)
{
    return -A*(C-Cm) + B*pow(C-Cm, 3) + Ca*pow(C-Ca, 3) + Cb*pow(C-Cb, 3);
}

const double AA = B + Ca + Cb;
const double BB = 3.0*(B*Cm + pow(Ca,2.0) + pow(Cb,2.0));
const double CC = 3.0*(B*pow(Cm,2.0) + pow(Ca,3.0) + pow(Cb,3.0)) - A;
const double DD = B*pow(Cm,3.0) + pow(Ca,4.0) + pow(Cb,4.0) - A*Cm;

template<typename T> inline
T energy_density(const T& C)
{
    return -0.5*A*pow(C-Cm,2.0) + 0.25*B*pow(C-Cm,4.0) + 0.25*Ca*pow(C-Ca,4.0) + 0.25*Cb*pow(C-Cb,4.0);
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


const int edge = 96;
const double deltaX = 1.0;
const double CFL = 2.0;
const double dt = pow(deltaX,4.0)*CFL/(32.0*D*K);
// Max. semi-implicit timestep should be timestep = pow(deltaX,2.0)/(32.0*D*K)

const double tolerance = 5.0e-10; // Choose wisely. 1e-5 is poor, 1e-8 fair, 1e-12 publication-quality -- may not converge before your deadline.
const unsigned int max_iter = 10000; // don't let the solver stagnate

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
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    #ifdef VANILLA
    srand(time(NULL)/(rank+2));
    #endif

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(2,0,edge,0,edge); // field 0 is c, field 1 is mu
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		#ifdef VANILLA
    	for (int n=0; n<nodes(grid); n++)
		    grid(n)[0] = 0.45 + 0.02 * ((double(rand())/RAND_MAX) - 0.5);
		#else
        const double q[2] = {0.1*sqrt(2.0), 0.1*sqrt(3.0)}; // produces stipes oriented 45 degrees to horizontal
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			grid(n)[0] = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
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
		if (rank==0) {
			std::cout<<"Timestep is "<<dt<<" (Co="<<CFL<<").";
	    #ifdef VANILLA
	        std::cout<<std::endl
	    #else
	     std::cout<<" Run "<<150000*0.005/dt<<" steps to match explicit dataset.\nIters\tEnergy\tMass"<<std::endl;
	    #endif
	    }

        double dV = 1.0;
        for (int d=0; d<dim; d++)
            dV *= dx(grid,d);

		double energy = 0.0;
		double mass = 0.0;
		for (int n=0; n<nodes(grid); n++) {
			MMSP::vector<int> x = position(grid,n);
			const double C = grid(x)[0];
			energy += dV*energy_density(C);
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

    newGrid.copy(oldGrid); // deep copy: includes data and ghost cells

	// Make sure the grid spacing is correct. Modify at will.
	for (int d=0; d<dim; d++) {
		dx(oldGrid,d) = deltaX;
		dx(newGrid,d) = deltaX;
	}

    double gridSize = 1.0*nodes(oldGrid);
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
        ferr.open("error.log", std::ofstream::out | std::ofstream::app);
    #endif

    // From Numerical Recipes: the following Jacobi spectral radius only applies to the Laplace equation, which this code is not.
    // Still, it's the hint of a direction for Chebyshev acceleration of omega.
    //const double jacobi_radius = (cos(2.0*M_PI/(g1(oldGrid,0)-g0(oldGrid,0)))  +  pow(dx(oldGrid,0)/dx(oldGrid,1),2)*cos(2.0*M_PI/(g1(oldGrid,1)-g0(oldGrid,1)))) /
    //                             (1.0+pow(dx(oldGrid,0)/dx(oldGrid,1),2)); // note: PI for Diriclet or Neumann, 2PI for periodic

    const double omega = 0.99; // Successive Over-Relaxation relaxation parameter (w=1 is Gauss-Seidel)

	for (int step=0; step<steps; step++) {
        double residual=1.0;
        unsigned int iter=0;

        while (iter<max_iter && residual>tolerance) {

            /*
                ==== RED-BLACK GAUSS SEIDEL ====

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
	    		        continue; // if odd

	    	    	const double cOld = oldGrid(n)[0];
	        		const double cLast = newGrid(n)[0];

	        		// A is defined by the last guess, stored in newGrid(x). It is a 2x2 matrix.
    	    		const double A11 = 1.0;                                   const double A12 = dt*D*lapWeight;
	    		    const double A21 = -nonlinear_coeff(cLast) - K*lapWeight; const double A22 = 1.0;

	    	    	// B is defined by the last value, stored in oldGrid(x), and the last guess, stored in newGrid(x). It is a 2x1 column.
	        		const double lapC = fringe_laplacian(newGrid, x, 0); // excludes central term
	        		const double lapU = fringe_laplacian(newGrid, x, 1); // excludes central term

    	    		const double B1 = cOld + D*dt*lapU;
	    		    const double B2 = expansive_dfdc(cOld) - K*lapC;

                    // Solve AX=B using Cramer's rule
	        		const double denom = A11*A22 - A12*A21; // det|A|
	        		const double cNew = (A22*B1 - B2*A12)/denom; // X1
    	    		const double uNew = (A11*B2 - B1*A21)/denom; // X2

                    newGrid(n)[0] = cNew;
                    newGrid(n)[1] = uNew;
                }
                ghostswap(newGrid);   // fill in the ghost cells; does nothing in serial
            }

            /*
                ==== RESIDUAL ====
                The residual is computed from the original matrix form, AX=B:
                any Old term goes into B, while any New term goes in AX.
            */

            double normB = 0.0;
            residual = 0.0;
       		for (int n=0; n<nodes(oldGrid); n++) {
    		    MMSP::vector<int> x = position(oldGrid,n);
           		const double lapC = field_laplacian(newGrid, x, 0);
	       		const double lapU = field_laplacian(newGrid, x, 1);

       	    	const double cOld = oldGrid(n)[0];
                const double cNew = newGrid(n)[0];
                const double uOld = oldGrid(n)[1];
                const double uNew = newGrid(n)[1];

        		const double B1 = cOld;
	   		    const double B2 = expansive_dfdc(cOld);

                const double AX1 = cNew - D*dt*lapU;
                const double AX2 = uNew - contractive_dfdc(cNew) + K*lapC;

                // Compute the Error from parts of the solution
                const double R1 = AX1 - B1;
                const double R2 = AX2 - B2;

                const double normBminusAX = pow(R1,2.0)/(2.0*gridSize) + pow(R2,2.0)/(2.0*gridSize);
                residual += normBminusAX;
                normB += pow(B1,2.0) + pow(B2,2.0);

                //newGrid(n)[0] = B1 + omega*(cNew - B1);
                //newGrid(n)[1] = B2 + omega*(uNew - B2);
                newGrid(n)[0] = (1.0 - omega)*cOld + omega*cNew;
                newGrid(n)[1] = (1.0 - omega)*uOld + omega*uNew;
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

   	    	residual = sqrt(residual)/sqrt(2.0*gridSize*normB);

   	    	#ifdef DEBUG
   	    	if (rank==0)
  	    	    ferr<<residual<<std::endl;
  	    	    //ferr<<residual<<"\t\t"<<omega<<std::endl;
            #endif

	    	iter++;
        }

        // Chebyshev acceleration of relaxation parameter
        // omega = (step==0) ? 1.0/(1.0-0.5*pow(jacobi_radius,2)) : 1.0/(1.0-0.25*omega*pow(jacobi_radius,2))  ;

        if (iter==max_iter && rank==0)
                std::cerr<<"    Solver stagnated on step "<<step<<": "<<iter<<" with residual="<<residual<<std::endl;

		oldGrid.copy(newGrid); // want to preserve the values in update for next iteration's first guess -- so don't swap, deep copy

		double energy = 0.0;
		double mass = 0.0;
		for (int n=0; n<nodes(oldGrid); n++) {
			MMSP::vector<int> x = position(oldGrid,n);
			const double C = oldGrid(x)[0];
			energy += dV*energy_density(C);
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
