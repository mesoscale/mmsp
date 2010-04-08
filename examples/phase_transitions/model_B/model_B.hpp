// model_B.hpp
// Algorithms for 2D and 3D implementation of "model B"
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MODELB_UPDATE
#define MODELB_UPDATE
#include"MMSP.hpp"
#include<cmath>

namespace MMSP{

double gaussian(double ave, double std)
{
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
	if (dim==2) {
		MMSP::grid<2,double> grid(1,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(x) = 0.0;
		}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,double> grid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(x) = 0.0;
		}

		MMSP::output(grid,filename);
	}
}

template <int dim, typename T> void update(MMSP::grid<dim,T>& grid, int steps)
{
	MMSP::grid<dim,T> update(grid);
	MMSP::grid<dim,T> temp(grid);
	MMSP::grid<dim,MMSP::vector<T> > noise(grid,2);

	double r = 1.0;
	double u = 1.0;
	double K = 1.0;
	double M = 1.0;
	double dt = 0.01;
	double kT = 0.02;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			temp(x) = -r*grid(x)+u*pow(grid(x),3)-K*laplacian(grid,x);
			for (int j=0; j<fields(noise); j++)
				noise(x)[j] = gaussian(0.0,sqrt(2.0*kT*M/(dt*volume(grid,x))));
		}
		ghostswap(temp);
		ghostswap(noise);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			update(x) = grid(x)+dt*M*laplacian(temp,x)+dt*div(noise,x);
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
