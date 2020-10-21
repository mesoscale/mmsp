// ising.cpp
// 2D and 3D Ising model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ISING_UPDATE
#define ISING_UPDATE
#include"MMSP.hpp"
#include"ising.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(0,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(0,0,2*L,0,L);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++)
			initGrid(i) = 2*(rand()%2)-1;

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

	double J = 2.0;
	double H = 1.0;
	double kT = (dim==2)?0.50:0.75;

	int gss = (dim==1)?nodes(spinGrid):int(sqrt(nodes(spinGrid)));
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

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
				for (int i=-1; i<2; i++) {
					x[0] += i;
					int spin = spinGrid(x);
					sum += (spin!=spin2)-(spin!=spin1);
					x[0] -= i;
				}
			} else if (dim==2) {
				for (int i=-1; i<2; i++) {
					x[0] += i;
					for (int j=-1; j<2; j++) {
						x[1] += j;
						int spin = spinGrid(x);
						sum += (spin!=spin2)-(spin!=spin1);
						x[1] -= j;
					}
					x[0] -= i;
				}
			} else if (dim==3) {
				for (int i=-1; i<2; i++) {
					x[0] += i;
					for (int j=-1; j<2; j++) {
						x[1] += j;
						for (int k=-1; k<2; k++) {
							x[2] += k;
							int spin = spinGrid(x);
							sum += (spin!=spin2)-(spin!=spin1);
							x[2] -= k;
						}
						x[1] -= j;
					}
					x[0] -= i;
				}
			}
			double dE = -0.5*J*sum-H*(spin2-spin1);

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
