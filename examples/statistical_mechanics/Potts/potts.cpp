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
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(0,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(0,0,2*L,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = rand()%20;

		output(initGrid,filename);
	}
}

template <int dim> void update(grid<dim,int>& spinGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(spinGrid);

	int Q = 20;
	double J = 1.0;
	double kT = (dim==3)?0.75:0.50;

	int gss = (dim==1)?nodes(spinGrid):int(sqrt(nodes(spinGrid)));
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

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
			if (dE<=0.0)
				spinGrid(p) = spin2;
			else if (r<exp(-dE/kT))
				spinGrid(p) = spin2;

			if (h%gss==0)
				ghostswap(spinGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
