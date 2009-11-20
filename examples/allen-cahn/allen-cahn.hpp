// allen-cahn.hpp
// Algorithms for 2D and 3D Allen-Cahn model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ALLENCAHN_UPDATE
#define ALLENCAHN_UPDATE
#include"MMSP.hpp"
#include<cmath>

namespace MMSP{

void update(MMSP::grid<2,double>& grid, int steps)
{
	MMSP::grid<2,double> update(grid);

	double dt = 0.01;
	double M = 1.0;
	double K = 1.0;
	double r = 1.0;
	double u = 1.0;

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double value = grid[x][y];
				double lap = (grid[x+1][y]-2.0*grid[x][y]+grid[x-1][y])/(dx(grid)*dx(grid))
				            +(grid[x][y+1]-2.0*grid[x][y]+grid[x][y-1])/(dy(grid)*dy(grid));

				update[x][y] = grid[x][y]-dt*M*(-r*value+u*pow(value,3)-K*lap);
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

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double value = grid[x][y][z];
					double lap = (grid[x+1][y][z]-2.0*grid[x][y][z]+grid[x-1][y][z])/(dx(grid)*dx(grid))
					            +(grid[x][y+1][z]-2.0*grid[x][y][z]+grid[x][y-1][z])/(dy(grid)*dy(grid))
					            +(grid[x][y][z+1]-2.0*grid[x][y][z]+grid[x][y][z-1])/(dz(grid)*dz(grid));

					update[x][y] = grid[x][y][z]-dt*M*(-r*value+u*pow(value,3)-K*lap);
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
