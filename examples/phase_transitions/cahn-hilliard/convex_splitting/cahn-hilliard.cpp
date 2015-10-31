// cahn-hilliard.hpp
// Algorithms for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"
#include"../energy.hpp"

namespace MMSP {

void generate(int dim, const char* filename)
{
	if (dim!=2) {
		std::cerr<<"ERROR: CHiMaD problems are 2-D, only!"<<std::endl;
		std::exit(-1);
	}

	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    const double q[2] = {0.1*std::sqrt(2.0), 0.1*std::sqrt(3.0)};

	if (dim==2) {
		MMSP::grid<2,double> grid(1,0,200,0,200);
		for (int d=0; d<dim; d++)
			dx(grid,d) = deltaX;

		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			grid(x) = 0.45 + 0.01 * std::cos(x[0]*dx(grid,0)*q[0] + x[1]*dx(grid,1)*q[1]);
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
void update(MMSP::grid<dim,T>& grid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

    ghostswap(grid);

	MMSP::grid<dim,T> update(grid);
	MMSP::grid<dim,T> temp(grid);
	// Make sure the grid spacing is correct
	for (int d=0; d<dim; d++) {
		dx(grid,d) = deltaX;
		dx(update,d) = deltaX;
		dx(temp,d) = deltaX;
	}

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			double c = grid(x);
			temp(x) = dfdc(c) - K*laplacian(grid,x);
		}
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		ghostswap(temp);

		double energy = 0.0;
		double mass = 0.0;
		for (int i=0; i<nodes(grid); i++) {
			MMSP::vector<int> x = position(grid,i);
			update(x) = grid(x) + dt*D*laplacian(temp,x);
			energy += dx(grid)*dy(grid)*energydensity(update(x));
			mass += dx(grid)*dy(grid)*update(x);
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
