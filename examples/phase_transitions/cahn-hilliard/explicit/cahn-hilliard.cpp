// cahn-hilliard.cpp
// Algorithms for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"

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
	MMSP::grid<dim,T> temp(grid);

	double r = 1.0;
	double u = 1.0;
	double K = 1.0;
	double M = 1.0;
	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			double phi = grid(i);
			temp(i) = -r*phi+u*pow(phi,3)-K*laplacian(grid,i);
		}
		ghostswap(temp);

		for (int i=0; i<nodes(grid); i++)
			update(i) = grid(i)+dt*M*laplacian(temp,i);

		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
