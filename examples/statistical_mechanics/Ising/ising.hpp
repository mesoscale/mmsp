// ising.hpp
// 2D and 3D Ising model 
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ISING_UPDATE
#define ISING_UPDATE
#include"MMSP.grid.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::grid<2,int> grid(1,0,128,0,128);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				grid[x][y] = rand()%2;

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,int> grid(1,0,64,0,64,0,64);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++)
					grid[x][y][z] = rand()%2;

		MMSP::output(grid,filename);
	}
}

void update(grid<2,int>& grid, int steps)
{
	double J = 1.0;
	double H = 1.0;

	const double kT = 0.50;
	int nx = xlength(grid);
	int ny = ylength(grid);
	int n = nx*ny;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = xmin(grid)+rand()%nx;
			int y = ymin(grid)+rand()%ny;
			int spin1 = grid[x][y];
			int spin2 = -spin1;
			
			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++) {
					int spin = grid[x+i][y+j];
					sum += (spin!=spin2)-(spin!=spin1);
				}
			double dE = -0.5*J*sum-H*(spin2-spin1);

			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid[x][y] = spin2;
			else if (r<exp(-dE/kT)) grid[x][y] = spin2;

			if (h%(ny)==0) ghostswap(grid);
		}
	}
}

void update(grid<3,int>& grid, int steps)
{
	double J = 1.0;
	double H = 1.0;

	const double kT = 0.75;
	int nx = xlength(grid);
	int ny = ylength(grid);
	int nz = zlength(grid);
	int n = nx*ny*nz;

	for (int step=0; step<steps; step++) {
		for (int h=0; h<n; h++) {
			int x = xmin(grid)+rand()%nx;
			int y = ymin(grid)+rand()%ny;
			int z = zmin(grid)+rand()%nz;
			int spin1 = grid[x][y][z];
			int spin2 = -spin1;

			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++)
					for (int k=-1; k<=1; k++) {
						int spin = grid[x+i][y+j][z+k];
						sum += (spin!=spin2)-(spin!=spin1);
					}
			double dE = -0.5*J*sum-H*(spin2-spin1);

			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid[x][y][z] = spin2;
			else if (r<exp(-dE/kT)) grid[x][y][z] = spin2;

			if (h%(ny*nz)==0) ghostswap(grid);
		}
	}
}

} // namespace MMSP

#endif
