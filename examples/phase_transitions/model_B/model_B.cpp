// model_B.cpp
// Algorithms for 2D and 3D implementation of "model B"
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MODELB_UPDATE
#define MODELB_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"model_B.hpp"

namespace MMSP{

double gaussian(double ave, double std)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	static double u = 0;
	static double v = 0;
	static bool saved = false;

	if (not saved) {
		start:
		double x = 2.0*double(rand())/double(RAND_MAX)-1.0;
		double y = 2.0*double(rand())/double(RAND_MAX)-1.0;

		double r = x*x+y*y;
		if (r<=0.0 or r>1.0) goto start;
		double d = sqrt(-2.0*log(r)/r);

		u = d*x;
		v = d*y;

		saved = true;
		return ave+u*std;
	}
	else {
		saved = false;
		return ave+v*std;
	}
}

void generate(int dim, const char* filename)
{
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(0,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 0.0;

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(0,0,2*L,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 0.0;

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 0.0;

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,T>& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(oldGrid);

	grid<dim,T> newGrid(oldGrid);
	grid<dim,T> temp(oldGrid);

	T r = 1.0;
	T u = 1.0;
	T K = 1.0;
	T M = 1.0;
	T dt = 0.01;
	T kT = 0.01;
	T dV = 1.0;

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int i=0; i<nodes(oldGrid); i++) {
			T phi = oldGrid(i);
			T noise = gaussian(0.0,sqrt(2.0*kT/(dt*dV)));
			temp(i) = -r*phi+u*pow(phi,3)-K*laplacian(oldGrid,i)+noise;
		}
		ghostswap(temp);

		for (int i=0; i<nodes(oldGrid); i++)
			newGrid(i) = oldGrid(i)+dt*M*laplacian(temp,i);

		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
