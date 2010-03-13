// MCgrid.anisotropy.hpp
// Anisotropic coarsening algorithm for 2D and 3D Monte Carlo methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MCGRID_UPDATE
#define MCGRID_UPDATE
#include"MCgrid.hpp"
#include"anisotropy.hpp"
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
			for (int y=y0; y<y1; y++) {
				if (x<32) {
					if (y<64) grid[x][y] = 2;
					else grid[x][y] = 3;
				}
				else if (x>96) {
					if (y<64) grid[x][y] = 2;
					else grid[x][y] = 3;
				}
				else {
					if (y<32 or y>96) grid[x][y] = 1;
					else grid[x][y] = 0;
				}
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
				for (int z=z0; z<z1; z++) {
					if (x<32) {
						if (y<64) grid[x][y][z] = 2;
						else grid[x][y][z] = 3;
					}
					else if (x>96) {
						if (y<64) grid[x][y][z] = 2;
						else grid[x][y][z] = 3;
					}
					else {
						if (y<32 or y>96) grid[x][y][z] = 1;
						else grid[x][y][z] = 0;
					}
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
			
			vector<int> nbors = neighbors(grid,x,y);
			int spin2 = nbors[rand()%length(nbors)];

			if (spin1!=spin2) {
				float dE = -energy(spin1,spin2);
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++) {
						int spin = grid[x+i][y+j];
						dE += energy(spin,spin2)-energy(spin,spin1);
					}
				float r = float(rand())/float(RAND_MAX);
				float ME = mobility(spin1,spin2)*energy(spin1,spin2);
				if (dE<=0.0 and r<ME) grid[x][y] = spin2;
				if (dE>0.0 and r<ME*exp(-dE/(kT*energy(spin1,spin2)))) grid[x][y] = spin2;
			}

			if (h%(ny/200)==0) ghostswap(grid);
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
			
			vector<int> nbors = neighbors(grid,x,y,z);
			int spin2 = nbors[rand()%length(nbors)];

			if (spin1!=spin2) {
				float dE = -energy(spin1,spin2);
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++)
						for (int k=-1; k<=1; k++) {
							int spin = grid[x+i][y+j][z+k];
							dE += energy(spin,spin2)-energy(spin,spin1);
						}
				float r = float(rand())/float(RAND_MAX);
				float ME = mobility(spin1,spin2)*energy(spin1,spin2);
				if (dE<=0.0 and r<ME) grid[x][y][z] = spin2;
				if (dE>0.0 and r<ME*exp(-dE/(kT*energy(spin1,spin2)))) grid[x][y][z] = spin2;
			}

			if (h%(ny*nz/200)==0) ghostswap(grid);
		}
	}
}

} // namespace MC

#endif
