// ising.cpp
// 2D and 3D Ising model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ISING_UPDATE
#define ISING_UPDATE
#include"MMSP.grid.hpp"
#include<cmath>
#include"ising.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		GRID1D initGrid(1,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(1,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

		output(initGrid,filename);
	}
}

template <int dim> void update(grid<dim,int>& spinGrid, int steps)
{
	double J = 2.0;
	double H = 1.0;
	double kT = (dim==2)?0.50:0.75;

	int gss = (dim==1)?nodes(spinGrid):int(sqrt(nodes(spinGrid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(spinGrid); h++) {
			// choose a random site
			int p = rand()%nodes(spinGrid);
			vector<int> x = position(spinGrid,p);
			int spin1 = spinGrid(p);

			// new spin is opposite old spin
			int spin2 = -spin1;

			// compute energy change
			double sum = -1.0;
			if (dim==1) {
				for (int i=-1; i<=1; i++) {
					int spin = spinGrid[x[0]+i];
					sum += (spin!=spin2)-(spin!=spin1);
				}
			} else if (dim==2) {
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++) {
						int spin = spinGrid[x[0]+i][x[1]+j];
						sum += (spin!=spin2)-(spin!=spin1);
					}
			} else if (dim==3) {
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++)
						for (int k=-1; k<=1; k++) {
							int spin = spinGrid[x[0]+i][x[1]+j][x[2]+k];
							sum += (spin!=spin2)-(spin!=spin1);
						}
			}
			double dE = -0.5*J*sum-H*(spin2-spin1);

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
