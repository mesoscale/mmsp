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
	if (dim==1) {
		GRID1D initGrid(2,0,128);

		int x0 = x0(initGrid);
		int x1 = x1(initGrid);

		for (int x=x0; x<x1; x++)
				double r = double(rand())/double(RAND_MAX);
				initGrid[x][0] = r;
				initGrid[x][1] = 1.0-r;
			}

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(2,0,128,0,128);

		int x0 = x0(initGrid);
		int x1 = x1(initGrid);
		int y0 = y0(initGrid);
		int y1 = y1(initGrid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++) {
				double r = double(rand())/double(RAND_MAX);
				initGrid[x][y][0] = r;
				initGrid[x][y][1] = 1.0-r;
			}

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(2,0,64,0,64,0,64);

		int x0 = x0(initGrid);
		int x1 = x1(initGrid);
		int y0 = y0(initGrid);
		int y1 = y1(initGrid);
		int z0 = z0(initGrid);
		int z1 = z1(initGrid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++) {
					double r = double(rand())/double(RAND_MAX);
					initGrid[x][y][z][0] = r;
					initGrid[x][y][z][1] = 1.0-r;
				}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	grid<dim,vector<T> > newGrid(oldGrid);
	grid<dim,T> wspace(oldGrid,1);

	double dt = 0.01;
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
		for (int n=0; n<nodes(oldGrid); n++) {
			double sum = 0.0;
			for (int i=1; i<fields(oldGrid); i++)
				sum += pow(oldGrid(n)[i],2);

			T C = oldGrid(n)[0];
			T lap = laplacian(oldGrid, n, 0);

			wspace[x][y] = -A*(C-Cmatrix)+B*pow(C-Cmatrix,3)
			               +Dalpha*pow(C-Calpha,3)+Dbeta*pow(C-Cbeta,3)
			               -gamma*(C-Calpha)*sum-kappa*lap;
		}
		ghostswap(wspace);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid, n);
			T lap = laplacian(wspace,x);
			T C = oldGrid(n)[0];

			newGrid(x)[0] = C+dt*D*lap;

			double sum = 0.0;
			for (int i=1; i<fields(oldGrid); i++)
				sum += pow(oldGrid(n)[i],2);

			vector<T> vlap = laplacian(oldGrid,x);
			for (int i=1; i<fields(oldGrid); i++) {
				T value = oldGrid(n)[i];
				newGrid(x)[i] = value-dt*L*(-gamma*pow(C-Calpha,2)+delta*pow(value,3)
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
