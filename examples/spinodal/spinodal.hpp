// spinodal.hpp
// Spinodal decomposition algorithms for 2D and 3D phase field methods
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


void update(MMSP::grid<2,double>& grid, int steps)
{
	MMSP::grid<2,double> update(grid);
	MMSP::grid<2,double> space(grid);
	MMSP::grid<2,vector<double> > noise(grid,2);

	double dt = 0.01;
	double epsilon = 0.05;
	double dV = dx(grid)*dy(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double value = grid[x][y];
				double lap = (grid[x+1][y]-2.0*grid[x][y]+grid[x-1][y])/(dx(grid)*dx(grid))
					        +(grid[x][y+1]-2.0*grid[x][y]+grid[x][y-1])/(dy(grid)*dy(grid));

				space[x][y] = value*(value*value-1.0)-lap;
				noise[x][y][0] = gaussian(0.0,epsilon*dt/dV);
				noise[x][y][1] = gaussian(0.0,epsilon*dt/dV);
			}
		ghostswap(space);
		ghostswap(noise);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double lap1 = (space[x+1][y]-2.0*space[x][y]+space[x-1][y])/(dx(grid)*dx(grid))
					         +(space[x][y+1]-2.0*space[x][y]+space[x][y-1])/(dy(grid)*dy(grid));

				double lap2 = (noise[x+1][y][0]-noise[x-1][y][0])/(2.0*dx(grid))
				             +(noise[x][y+1][1]-noise[x][y-1][1])/(2.0*dy(grid));

				update[x][y] = grid[x][y]+dt/(2.0*dV)*lap1+lap2;
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(MMSP::grid<3,double>& grid, int steps)
{
	MMSP::grid<3,double> update(grid);
	MMSP::grid<3,double> space(grid);
	MMSP::grid<3,vector<double> > noise(grid,3);

	double dt = 0.01;
	double epsilon = 0.05;
	double dV = dx(grid)*dy(grid)*dz(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double value = grid[x][y][z];
					double lap = (grid[x+1][y][z]-2.0*grid[x][y][z]+grid[x-1][y][z])/(dx(grid)*dx(grid))
						        +(grid[x][y+1][z]-2.0*grid[x][y][z]+grid[x][y-1][z])/(dy(grid)*dy(grid))
						        +(grid[x][y][z+1]-2.0*grid[x][y][z]+grid[x][y][z-1])/(dz(grid)*dz(grid));

					space[x][y][z] = value*(value*value-1.0)-lap;
					noise[x][y][z][0] = gaussian(0.0,epsilon*dt/dV);
					noise[x][y][z][1] = gaussian(0.0,epsilon*dt/dV);
					noise[x][y][z][2] = gaussian(0.0,epsilon*dt/dV);
				}
		ghostswap(space);
		ghostswap(noise);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double lap1 = (space[x+1][y][z]-2.0*space[x][y][z]+space[x-1][y][z])/(dx(grid)*dx(grid))
						         +(space[x][y+1][z]-2.0*space[x][y][z]+space[x][y-1][z])/(dy(grid)*dy(grid))
						         +(space[x][y][z+1]-2.0*space[x][y][z]+space[x][y][z-1])/(dz(grid)*dz(grid));

					double lap2 = (noise[x+1][y][z][0]-noise[x-1][y][z][0])/(2.0*dx(grid))
								 +(noise[x][y+1][z][1]-noise[x][y-1][z][1])/(2.0*dy(grid))
								 +(noise[x][y][z+1][2]-noise[x][y][z-1][2])/(2.0*dz(grid));

					update[x][y][z] = grid[x][y][z]+dt/(2.0*dV)*lap1+lap2;
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
