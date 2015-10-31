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

double expansive_dfdc(const double& C)
{
    return -A*(C-Cm) + B*pow(C-Cm, 3) + Ca*pow(C-Ca, 3) + Cb*pow(C-Cb, 3);
}

namespace MMSP {

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
		MMSP::grid<2,vector<double> > grid(2,0,200,0,200);
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

	MMSP::grid<dim,vector<T> > update(grid);
	MMSP::grid<dim,vector<T> > temp(grid);
	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		dx(update,d) = deltaX;
		dx(temp,d) = deltaX;
	}

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			double c = grid(x)[0];
			temp(x)[0] = dfdc(c) - K*laplacian(grid,x)[0];
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		ghostswap(temp);

		double energy = 0.0;
		double mass = 0.0;
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			update(x)[0] = grid(x)[0] + dt*D*laplacian(temp,x)[0];
			energy += dx(grid)*dy(grid)*energydensity(update(x)[0]);
			mass += dx(grid)*dy(grid)*update(x)[0];
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		double myEnergy = energy;
		double myMass = mass;
		MPI::COMM_WORLD.Reduce(&myEnergy, &energy, 1, MPI_DOUBLE, MPI_SUM, 0);
		MPI::COMM_WORLD.Reduce(&myMass, &mass, 1, MPI_DOUBLE, MPI_SUM, 0);
		#endif
		#ifndef DEBUG
		if (rank==0)
			std::cout<<energy<<'\t'<<mass<<'\n';
		#endif

		swap(grid,update);
		ghostswap(grid);
	}
	#ifndef DEBUG
	if (rank==0)
		std::cout<<std::flush;
	#endif
}

} // MMSP
#endif

#include"MMSP.main.hpp"
