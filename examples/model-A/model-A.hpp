// model-A.hpp
// Algorithms for 2D and 3D implementation of "model A"
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MODELA_UPDATE
#define MODELA_UPDATE
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


void update(MMSP::grid<2,double>& grid, int steps)
{
	MMSP::grid<2,double> update(grid);

	double dt = 0.01;
	double M = 1.0;
	double K = 1.0;
	double r = 1.0;
	double u = 1.0;
	double kT = 1.0;
	double dV = dx(grid)*dy(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double value = grid[x][y];
				double lap = (grid[x+1][y]-2.0*grid[x][y]+grid[x-1][y])/(dx(grid)*dx(grid))
				            +(grid[x][y+1]-2.0*grid[x][y]+grid[x][y-1])/(dy(grid)*dy(grid));
				double noise = gaussian(0.0,sqrt(2.0*kT*M/dV/dt));

				update[x][y] = grid[x][y]-dt*M*(-r*value+u*pow(value,3)-K*lap)+dt*noise;
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(MMSP::grid<3,double>& grid, int steps)
{
	MMSP::grid<3,double> update(grid);

	double dt = 0.01;
	double M = 1.0;
	double K = 1.0;
	double r = 1.0;
	double u = 1.0;
	double kT = 1.0;
	double dV = dx(grid)*dy(grid)*dz(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double value = grid[x][y][z];
					double lap = (grid[x+1][y][z]-2.0*grid[x][y][z]+grid[x-1][y][z])/(dx(grid)*dx(grid))
					            +(grid[x][y+1][z]-2.0*grid[x][y][z]+grid[x][y-1][z])/(dy(grid)*dy(grid))
					            +(grid[x][y][z+1]-2.0*grid[x][y][z]+grid[x][y][z-1])/(dz(grid)*dz(grid));
					double noise = gaussian(0.0,sqrt(2.0*kT*M/dV/dt));

					update[x][y][z] = grid[x][y][z]-dt*M*(-r*value+u*pow(value,3)-K*lap)+dt*noise;
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
