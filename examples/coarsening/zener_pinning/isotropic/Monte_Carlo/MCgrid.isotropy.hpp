// MCgrid.isotropy.hpp
// Isotropic coarsening algorithm for 2D and 3D Monte Carlo methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MCGRID_UPDATE
#define MCGRID_UPDATE
#include"MCgrid.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::MCgrid2D grid(128,128);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				grid[x][y] = 1+rand()%100;

		for (int i=0; i<100; i++) {
			int x = x0+rand()%(x1-x0);
			int y = y0+rand()%(y1-y0);
			grid[x][y] = 0;
		}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::MCgrid3D grid(64,64,64);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++)
					grid[x][y][z] = 1+rand()%100;

		for (int i=0; i<100; i++) {
			int x = x0+rand()%(x1-x0);
			int y = y0+rand()%(y1-y0);
			int z = z0+rand()%(z1-z0);
			grid[x][y][z] = 0;
		}

		MMSP::output(grid,filename);
	}
}

void update(MCgrid2D& grid, int steps)
{
	const float kT = 0.50;
	int nx = xlength(grid);
	int ny = ylength(grid);
	int n = nx*ny;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = x0(grid)+rand()%nx;
			int y = y0(grid)+rand()%ny;
			int spin1 = grid[x][y];
			
			if (spin1!=0) {
				vector<int> nbors = neighbors(grid,x,y);
				int spin2 = nbors[rand()%length(nbors)];

				if (spin1!=spin2 and spin2!=0) {
					float dE = -1.0;
					for (int i=-1; i<=1; i++)
						for (int j=-1; j<=1; j++) {
							int spin = grid[x+i][y+j];
							dE += (spin!=spin2)-(spin!=spin1);
						}
					float r = float(rand())/float(RAND_MAX);
					if (dE<=0.0) grid[x][y] = spin2;
					else if (r<exp(-dE/kT)) grid[x][y] = spin2;
				}

				if (h%(ny)==0) ghostswap(grid);
			}
		}
	}
}

void update(MCgrid3D& grid, int steps)
{
	const float kT = 0.75;
	int nx = xlength(grid);
	int ny = ylength(grid);
	int nz = zlength(grid);
	int n = nx*ny*nz;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = x0(grid)+rand()%nx;
			int y = y0(grid)+rand()%ny;
			int z = z0(grid)+rand()%nz;
			int spin1 = grid[x][y][z];
			
			if (spin1!=0) {
				vector<int> nbors = neighbors(grid,x,y,z);
				int spin2 = nbors[rand()%length(nbors)];

				if (spin1!=spin2 and spin2!=0) {
					float dE = -1.0;
					for (int i=-1; i<=1; i++)
						for (int j=-1; j<=1; j++)
							for (int k=-1; k<=1; k++) {
								int spin = grid[x+i][y+j][z+k];
								dE += (spin!=spin2)-(spin!=spin1);
							}
					float r = float(rand())/float(RAND_MAX);
					if (dE<=0.0) grid[x][y][z] = spin2;
					else if (r<exp(-dE/kT)) grid[x][y][z] = spin2;
				}

				if (h%(ny*nz)==0) ghostswap(grid);
			}
		}
	}
}

} // namespace MMSP

#endif
