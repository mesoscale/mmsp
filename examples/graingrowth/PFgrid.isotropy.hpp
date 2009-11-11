// PFgrid.isotropy.hpp
// Isotropic coarsening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef PFGRID_UPDATE
#define PFGRID_UPDATE
#include"PFgrid.hpp"
#include<cmath>

namespace MMSP{

void update(PFgrid2D& grid, int steps)
{
	PFgrid2D update(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double sum = 0.0;
				for (int i=0; i<fields(grid); i++)
					sum += pow(grid[x][y][i],2);

				for (int i=0; i<fields(grid); i++) {
					double value = grid[x][y][i];
					double lap =
					  (grid[x+1][y][i]-2.0*grid[x][y][i]+grid[x-1][y][i])/(dx(grid)*dx(grid))
					 +(grid[x][y+1][i]-2.0*grid[x][y][i]+grid[x][y-1][i])/(dy(grid)*dy(grid));

					const double dt = 0.01;
					update[x][y][i] = value-dt*(-value-pow(value,3)+2.0*(value*sum-lap));
				}
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(PFgrid3D& grid, int steps)
{
	PFgrid3D update(grid);

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double sum = 0.0;
					for (int i=0; i<fields(grid); i++)
						sum += pow(grid[x][y][z][i],2);

					for (int i=0; i<fields(grid); i++) {
						double value = grid[x][y][z][i];
						double lap =
						  (grid[x+1][y][z][i]-2.0*grid[x][y][z][i]+grid[x-1][y][z][i])/(dx(grid)*dx(grid))
						 +(grid[x][y+1][z][i]-2.0*grid[x][y][z][i]+grid[x][y-1][z][i])/(dy(grid)*dy(grid))
						 +(grid[x][y][z+1][i]-2.0*grid[x][y][z][i]+grid[x][y][z-1][i])/(dz(grid)*dz(grid));

						const double dt = 0.01;
						update[x][y][z][i] = value-dt*(-value-pow(value,3)+2.0*(value*sum-lap));
					}
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
