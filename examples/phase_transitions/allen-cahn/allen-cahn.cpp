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
	if (dim==1) {
		MMSP::grid<1,double> grid(1,0,128);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,double> grid(1,0,128,0,128);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,double> grid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = 1.0-2.0*double(rand())/double(RAND_MAX);

		MMSP::output(grid,filename);
	}
}

template <int dim, typename T> void update(MMSP::grid<dim,T>& grid, int steps)
{
	MMSP::grid<dim,T> update(grid);

	double r = 1.0;
	double u = 1.0;
	double K = 1.0;
	double M = 1.0;
	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			double phi = grid(i);
			update(i) = phi-dt*M*(-r*phi+u*pow(phi,3)-K*laplacian(grid,i));
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
