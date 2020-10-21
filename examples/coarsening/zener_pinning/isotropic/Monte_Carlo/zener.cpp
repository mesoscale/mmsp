// zener.cpp
// Algorithms for 2D and 3D Monte Carlo grain growth with Zener pinning
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ZENER_UPDATE
#define ZENER_UPDATE
#include<cmath>
#include"MMSP.hpp"
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
			double r = 32-x[0]%64;
			if (r<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				initGrid(p) = 0;
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(0,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double r = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (r<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++)
					initGrid(p) = 0;
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,2*L,0,L,0,L/2);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double r = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (r<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++)
					for (p[2]=x[2]-1; p[2]<=x[2]+1; p[2]++)
						initGrid(p) = 0;
		}

		output(initGrid,filename);
	}
}

template <int dim> void update(grid<dim,int>& mcGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(mcGrid);

	const double kT = (dim==3)?0.75:0.50;
	int gss = int(nodes(mcGrid));
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int h=0; h<nodes(mcGrid); h++) {
			// choose a random node
			int p = rand()%nodes(mcGrid);
			vector<int> x = position(mcGrid,p);
			int spin1 = mcGrid(p);

			if (spin1!=0) {
				// determine neighboring spins
				sparse<bool> neighbors;
				if (dim==1) {
					for (int i=-1; i<2; i++) {
						x[0] += i;
						int spin = mcGrid(x);
						set(neighbors,spin) = true;
						x[0] -= i;
					}
				} else if (dim==2) {
					for (int i=-1; i<2; i++) {
						x[0] += i;
						for (int j=-1; j<2; j++) {
							x[1] += j;
							int spin = mcGrid(x);
							set(neighbors,spin) = true;
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
								int spin = mcGrid(x);
								set(neighbors,spin) = true;
								x[2] -= k;
							}
							x[1] -= j;
						}
						x[0] -= i;
					}
				}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
					// compute energy change
					double dE = -1.0;
					if (dim==1) {
						for (int i=-1; i<2; i++) {
							x[0] += i;
							int spin = mcGrid(x);
							dE += (spin!=spin2)-(spin!=spin1);
							x[0] -= i;
						}
					} else if (dim==2) {
						for (int i=-1; i<2; i++) {
							x[0] += i;
							for (int j=-1; j<2; j++){
								x[1] += j;
								int spin = mcGrid(x);
								dE += (spin!=spin2)-(spin!=spin1);
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
									int spin = mcGrid(x);
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
						mcGrid(p) = spin2;
					else if (r<exp(-dE/kT))
						mcGrid(p) = spin2;
				}
			}
			if (h%gss==0)
				ghostswap(mcGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
