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
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(0,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = 32-x[0]%64;
			if (d<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			initGrid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				initGrid(x) = 0;
				x[d]+=2;
				initGrid(x) = 0;
				x[d]--;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=128;
		GRID2D initGrid(0,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (d<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			initGrid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				initGrid(x) = 0;
				x[d]+=2;
				initGrid(x) = 0;
				x[d]--;
			}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,L,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2)+pow(32-x[2]%64,2));
			if (d<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			initGrid(x) = 0;
			for (int d=0; d<dim; d++) {
				x[d]--;
				initGrid(x) = 0;
				x[d]+=2;
				initGrid(x) = 0;
				x[d]--;
			}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,T>& spinGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    rank = MPI::COMM_WORLD.Get_rank();
    #endif

	ghostswap(spinGrid);

	const double kT = (dim==3)?0.75:0.50;
	int gss = int(sqrt(nodes(spinGrid)));
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

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
					if (dE<=0.0)
						spinGrid(p) = spin2;
					else if (r<exp(-dE/kT))
						spinGrid(p) = spin2;
				}
			}
			if (h%gss==0)
				ghostswap(spinGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

