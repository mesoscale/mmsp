// potts.cpp
// 2D and 3D Potts model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef POTTS_UPDATE
#define POTTS_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"potts.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		GRID1D initGrid(1,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(1,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}
}

template <int dim> void update(grid<dim,int>& spinGrid, int steps)
{
	int Q = 20;
	double J = 1.0;
	double kT = (dim==3)?0.75:0.50;

	int gss = (dim==1)?nodes(spinGrid):int(sqrt(nodes(spinGrid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(spinGrid); h++) {
			// choose a random node
			int p = rand()%nodes(spinGrid);
			vector<int> x = position(spinGrid,p);
			int spin1 = spinGrid(p);

			// choose a new spin randomly
			int spin2 = rand()%Q;

			// compute energy change
			double dE = -1.0;
			if (dim==2) {
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
			dE *= J;

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) spinGrid(p) = spin2;
			else if (r<exp(-dE/kT)) spinGrid(p) = spin2;

			if (h%gss==0) ghostswap(spinGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
