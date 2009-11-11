// sparsePF.isotropy.hpp
// Isotropic coarsening algorithms for 2D and 3D sparsePF methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef SPARSEPF_UPDATE
#define SPARSEPF_UPDATE
#include"sparsePF.hpp"
#include<cmath>

namespace MMSP{

void update(sparsePF2D& grid, int steps)
{
	for (int step=0; step<steps; step++) {
		sparsePF2D update(grid);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				vector<int> nbors = neighbors(grid,x,y);
				if (length(nbors)==1)
					set(update[x][y],nbors[0]) = grid[x][y][nbors[0]];

				else {
					double sum = 0.0;
					for (int k=0; k<length(grid[x][y]); k++)
						sum += pow(value(grid[x][y],k),2);

					for (int k=0; k<length(nbors); k++) {
						int i = nbors[k];
						double value = grid[x][y][i];
						double lap =
						  (grid[x+1][y][i]-2.0*grid[x][y][i]+grid[x-1][y][i])/(dx(grid)*dx(grid))
						 +(grid[x][y+1][i]-2.0*grid[x][y][i]+grid[x][y-1][i])/(dy(grid)*dy(grid));

						const double dt = 0.01;
						const double epsilon = 1.0e-8;
						double temp = value-dt*(-value-pow(value,3)+2.0*(value*sum-lap));
						if (temp>epsilon) set(update[x][y],i) = temp;
					}
				}
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(sparsePF3D& grid, int steps)
{
	for (int step=0; step<steps; step++) {
		sparsePF3D update(grid);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					vector<int> nbors = neighbors(grid,x,y,z);
					if (length(nbors)==1)
						set(update[x][y][z],nbors[0]) = grid[x][y][z][nbors[0]];

					else {
						double sum = 0.0;
						for (int k=0; k<length(grid[x][y][z]); k++)
							sum += pow(value(grid[x][y][z],k),2);

						for (int k=0; k<length(grid[x][y][z]); k++) {
							int i = index(grid[x][y][z],k);
							double value = grid[x][y][z][i];
							double lap =
							  (grid[x+1][y][z][i]-2.0*grid[x][y][z][i]+grid[x-1][y][z][i])/(dx(grid)*dx(grid))
							 +(grid[x][y+1][z][i]-2.0*grid[x][y][z][i]+grid[x][y-1][z][i])/(dy(grid)*dy(grid))
							 +(grid[x][y][z+1][i]-2.0*grid[x][y][z][i]+grid[x][y][z-1][i])/(dz(grid)*dz(grid));

							const double dt = 0.01;
							const double epsilon = 1.0e-8;
							double temp = value-dt*(-value-pow(value,3)+2.0*(value*sum-lap));
							if (temp>epsilon) set(grid[x][y][z],i) = temp;
						}
					}
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
