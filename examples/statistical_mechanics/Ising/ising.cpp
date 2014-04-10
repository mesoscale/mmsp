// ising.hpp
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
	if (dim==2) {
		MMSP::grid<2,int> grid(1,0,128,0,128);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = 2*(rand()%2)-1;

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,int> grid(1,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++)
			grid(i) = 2*(rand()%2)-1;

		output(grid,filename);
	}
}

void update(grid<2,int>& grid, int steps)
{
	double J = 2.0;
	double H = 1.0;
	double kT = 0.50;

	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random site
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);

			// new spin is opposite old spin
			int spin2 = -spin1;

			// compute energy change
			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++) {
					int spin = grid[x[0]+i][x[1]+j];
					sum += (spin!=spin2)-(spin!=spin1);
				}
			double dE = -0.5*J*sum-H*(spin2-spin1);

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid(p) = spin2;
			else if (r<exp(-dE/kT)) grid(p) = spin2;

			if (h%gss==0) ghostswap(grid);
		}
	}
}

void update(grid<3,int>& grid, int steps)
{
	double J = 2.0;
	double H = 1.0;
	double kT = 0.75;

	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random site
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);

			// new spin is opposite old spin
			int spin2 = -spin1;

			// compute energy change
			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++)
					for (int k=-1; k<=1; k++) {
						int spin = grid[x[0]+i][x[1]+j][x[2]+k];
						sum += (spin!=spin2)-(spin!=spin1);
					}
			double dE = -0.5*J*sum-H*(spin2-spin1);

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid(p) = spin2;
			else if (r<exp(-dE/kT)) grid(p) = spin2;

			if (h%gss==0) ghostswap(grid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
