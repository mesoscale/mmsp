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

template <int dim, typename T> void update(grid<dim,T>& spinGrid, int steps)
{
	const double kT = (dim==3)?0.75:0.50;
	int gss = int(sqrt(nodes(spinGrid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(spinGrid); h++) {
			// choose a random node
			int p = rand()%nodes(spinGrid);
			int spin1 = spinGrid(p);

			if (spin1!=0) {
				vector<int> x = position(spinGrid,p);
				// determine neighboring spins
				sparse<bool> neighbors;
				set(neighbors,spinGrid(x)) = true;
				for (int d=0; d<1; d++) {
					x[d]--;
					set(neighbors,spinGrid(x)) = true;
					x[d]+=2;
					set(neighbors,spinGrid(x)) = true;
					x[d]--;
				}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
					// compute energy change
					double dE = -1.0;
					if (dim==1) {
						for (int i=-1; i<=1; i++) {
							x[0] += i;
							int spin = spinGrid(x);
							dE += (spin!=spin2)-(spin!=spin1);
							x[0] -= i;
						}
					} else if (dim==2) {
						for (int i=-1; i<=1; i++) {
							x[0] += i;
							for (int j=-1; j<=1; j++) {
								x[1] += j;
								int spin = spinGrid(x);
								dE += (spin!=spin2)-(spin!=spin1);
								x[1] -= j;
							}
							x[0] -= i;
						}
					} else if (dim==3) {
						for (int i=-1; i<=1; i++) {
							x[0] += i;
							for (int j=-1; j<=1; j++) {
								x[1] += j;
								for (int k=-1; k<=1; k++) {
									x[2] += k;
									int spin = spinGrid(x);
									dE += (spin!=spin2)-(spin!=spin1);
									x[2] -= k;
								}
								x[1] -= j;
							}
							x[0] -= i;
						}
					}

					// attempt a spin flip
					double r = double(rand())/double(RAND_MAX);
					if (dE<=0.0) spinGrid(p) = spin2;
					else if (r<exp(-dE/kT)) spinGrid(p) = spin2;
				}
			}
			if (h%gss==0) ghostswap(spinGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

