// graingrowth.cpp
// Algorithms for 2D and 3D Monte Carlo grain growth with Zener pinning
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MCGRID_UPDATE
#define MCGRID_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"zener.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		GRID1D grid(0,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = 64.0-x[0];
			if (d<32.0) grid(i) = 2;
			else grid(i) = 1;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> x = position(grid,i);
			grid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				grid(x) = 0;
				x[d]+=2;
				grid(x) = 0;
				x[d]--;
			}
		}

		output(grid,filename);
	}

	if (dim==2) {
		GRID2D grid(0,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) grid(i) = 2;
			else grid(i) = 1;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> x = position(grid,i);
			grid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				grid(x) = 0;
				x[d]+=2;
				grid(x) = 0;
				x[d]--;
			}
		}

		output(grid,filename);
	}

	if (dim==3) {
		GRID3D grid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) grid(i) = 2;
			else grid(i) = 1;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> x = position(grid,i);
			grid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				grid(x) = 0;
				x[d]+=2;
				grid(x) = 0;
				x[d]--;
			}
		}

		output(grid,filename);
	}
}

template <typename T> void update(MMSP::grid<1,T>& grid, int steps)
{
	const double kT = 0.50;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node
			int p = rand()%nodes(grid);
			int spin1 = grid(p);

			if (spin1!=0) {
				vector<int> x = position(grid,p);
				// determine neighboring spins
				sparse<bool> neighbors;
				set(neighbors,grid(x)) = true;
				for (int d=0; d<1; d++) {
					x[d]--;
					set(neighbors,grid(x)) = true;
					x[d]+=2;
					set(neighbors,grid(x)) = true;
					x[d]--;
				}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
					// compute energy change
					double dE = -1.0;
					for (int i=-1; i<=1; i++) {
						int spin = grid[x[0]+i];
						dE += (spin!=spin2)-(spin!=spin1);
					}

					// attempt a spin flip
					double r = double(rand())/double(RAND_MAX);
					if (dE<=0.0) grid(p) = spin2;
					else if (r<exp(-dE/kT)) grid(p) = spin2;
				}
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

template <typename T> void update(MMSP::grid<2,T>& grid, int steps)
{
	const double kT = 0.50;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node
			int p = rand()%nodes(grid);
			int spin1 = grid(p);

			if (spin1!=0) {
				vector<int> x = position(grid,p);
				// determine neighboring spins
				sparse<bool> neighbors;
				set(neighbors,grid(x)) = true;
				for (int d=0; d<2; d++) {
					x[d]--;
					set(neighbors,grid(x)) = true;
					x[d]+=2;
					set(neighbors,grid(x)) = true;
					x[d]--;
				}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
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
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

template <typename T> void update(MMSP::grid<3,T>& grid, int steps)
{
	const double kT = 0.75;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node
			int p = rand()%nodes(grid);
			int spin1 = grid(p);

			if (spin1!=0) {
				vector<int> x = position(grid,p);
				// determine neighboring spins
				sparse<bool> neighbors;
				set(neighbors,grid(x)) = true;
				for (int d=0; d<3; d++) {
					x[d]--;
					set(neighbors,grid(x)) = true;
					x[d]+=2;
					set(neighbors,grid(x)) = true;
					x[d]--;
				}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
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
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

