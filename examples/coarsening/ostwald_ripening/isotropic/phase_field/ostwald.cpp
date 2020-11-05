// ostwald.cpp
// Ostwald ripening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef OSTWALD_UPDATE
#define OSTWALD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"ostwald.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(2,0,L);

		for (int n=0; n<nodes(initGrid); n++) {
				double r = double(rand())/double(RAND_MAX);
				initGrid(n)[0] = r;
				initGrid(n)[1] = 1.0-r;
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=128;
		GRID2D initGrid(2,0,2*L,0,L);

		for (int n=0; n<nodes(initGrid); n++) {
				double r = double(rand())/double(RAND_MAX);
				initGrid(n)[0] = r;
				initGrid(n)[1] = 1.0-r;
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(2,0,2*L,0,L,0,L/4);

		for (int n=0; n<nodes(initGrid); n++) {
				double r = double(rand())/double(RAND_MAX);
				initGrid(n)[0] = r;
				initGrid(n)[1] = 1.0-r;
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);
	grid<dim,T> wspace(oldGrid,1);

	double dt = 0.001;
	double L = 1.0;
	double D = 1.0;

	double Calpha = 0.05;
	double Cbeta = 0.95;
	double Cmatrix = 0.5*(Calpha+Cbeta);
	double A = 2.0;
	double B = A/pow(Cbeta-Cmatrix,2);
	double gamma = 2.0/pow(Cbeta-Calpha,2);
	double delta = 1.0;
	double epsilon = 3.0;
	double Dalpha = gamma/pow(delta,2);
	double Dbeta = gamma/pow(delta,2);
	double kappa = 2.0;


	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);
		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			double sum = 0.0;
			for (int i=1; i<fields(oldGrid); i++)
				sum += pow(oldGrid(n)[i],2);

			T C = oldGrid(n)[0];
			T lap = laplacian(oldGrid, x, 0);

			wspace(x) = -A*(C-Cmatrix)+B*pow(C-Cmatrix,3)
			               +Dalpha*pow(C-Calpha,3)+Dbeta*pow(C-Cbeta,3)
			               -gamma*(C-Calpha)*sum-kappa*lap;
		}
		ghostswap(wspace);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			T lap = laplacian(wspace,x);
			T C = oldGrid(n)[0];

			newGrid(n)[0] = C+dt*D*lap;

			double sum = 0.0;
			for (int i=1; i<fields(oldGrid); i++)
				sum += pow(oldGrid(n)[i],2);

			vector<T> vlap = laplacian(oldGrid,x);
			for (int i=1; i<fields(oldGrid); i++) {
				T value = oldGrid(n)[i];
				newGrid(n)[i] = value-dt*L*(-gamma*pow(C-Calpha,2)+delta*pow(value,3)
				                            +epsilon*value*(sum-pow(value,2))-kappa*vlap[i]);
			}
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
