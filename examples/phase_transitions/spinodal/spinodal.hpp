// spinodal.hpp
// Algorithm for 2D and 3D spinodal decomposition phase field model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef SPINODAL_UPDATE
#define SPINODAL_UPDATE
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

		for (int i=0; i<nodes(grid); i++) grid(i) = 0.0;

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,double> grid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) grid(i) = 0.0;

		MMSP::output(grid,filename);
	}
}

template <int dim> void update(MMSP::grid<dim,double>& grid, int steps)
{
	MMSP::grid<dim,double> update(grid);
	MMSP::grid<dim,double> temp(grid);

	double dt = 0.01;
	double dV = 1.0;
	double epsilon = 0.05;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			double noise = gaussian(0.0,sqrt(epsilon*dt/dV));
			double phi = grid(i);
			temp(i) = -phi+pow(phi,3)-laplacian(grid,i)+noise;
		}
		ghostswap(temp);

		for (int i=0; i<nodes(grid); i++)
			update(i) = grid(i)+dt/(2.0*dV)*laplacian(temp,i);

		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
