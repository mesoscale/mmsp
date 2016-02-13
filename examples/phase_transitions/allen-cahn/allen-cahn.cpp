// allen-cahn.cpp
// Algorithms for 2D and 3D Allen-Cahn model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ALLENCAHN_UPDATE
#define ALLENCAHN_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"allen-cahn.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		GRID1D initGrid(1,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(1,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,T>& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    rank = MPI::COMM_WORLD.Get_rank();
    #endif

	grid<dim,T> newGrid(oldGrid);

	double r = 1.0;
	double u = 1.0;
	double K = 1.0;
	double M = 1.0;
	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int i=0; i<nodes(oldGrid); i++) {
			T phi = oldGrid(i);
			newGrid(i) = phi-dt*M*(-r*phi+u*pow(phi,3)-K*laplacian(oldGrid,i));
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
