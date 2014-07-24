// graingrowth.hpp
// Algorithms for 2D and 3D isotropic Monte Carlo grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"graingrowth.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		MMSP::grid<1,int> grid(0,0,128);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = rand()%100;

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,int> grid(0,0,128,0,128);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = rand()%100;

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,int> grid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = rand()%100;

		output(grid,filename);
	}
}

void update(MMSP::grid<2,int>& grid, int steps)
{
	const double kT = 0.50;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node 
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);
			
			// determine neighboring spins
			sparse<bool> neighbors;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++) {
					int spin = grid[x[0]+i][x[1]+j];
					set(neighbors,spin) = true;
				}

			// choose a random neighbor spin
			int spin2 = index(neighbors,rand()%length(neighbors));

			if (spin1!=spin2) {
				// compute energy change
				double dE = -1.0;
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++) {
						int spin = grid[x[0]+i][x[1]+j];
						dE += (spin!=spin2)-(spin!=spin1);
					}

				// attempt a spin flip
				double r = double(rand())/double(RAND_MAX);
				if (dE<=0.0) grid(p) = spin2;
				else if (r<exp(-dE/kT)) grid(p) = spin2;
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

void update(MMSP::grid<3,int>& grid, int steps)
{
	const double kT = 0.75;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node 
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);

			// determine neighboring spins
			sparse<bool> neighbors;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++)
					for (int k=-1; k<=1; k++) {
						int spin = grid[x[0]+i][x[1]+j][x[2]+k];
						set(neighbors,spin) = true;
					}

			// choose a random neighbor spin
			int spin2 = index(neighbors,rand()%length(neighbors));

			if (spin1!=spin2) {
				// compute energy change
				double dE = -1.0;
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++)
						for (int k=-1; k<=1; k++) {
							int spin = grid[x[0]+i][x[1]+j][x[2]+k];
							dE += (spin!=spin2)-(spin!=spin1);
						}

				// attempt a spin flip
				double r = double(rand())/double(RAND_MAX);
				if (dE<=0.0) grid(p) = spin2;
				else if (r<exp(-dE/kT)) grid(p) = spin2;
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
