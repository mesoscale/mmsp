// zener.cpp
// Anisotropic coarsening algorithms for 2D and 3D Monte Carlo methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ZENER_UPDATE
#define ZENER_UPDATE
#include"MMSP.hpp"
#include"anisotropy.hpp"
#include"zener.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		GRID1D initGrid(0,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = 64.0-x[0];
			if (d<32.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		for (int j=0; j<50; j++) {
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
			double d = sqrt(pow(L/2-x[0],2)+pow(L/2-x[1],2));
			if (d<L/4.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		for (int j=0; j<50; j++) {
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
		GRID3D initGrid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) initGrid(i) = 2;
			else initGrid(i) = 1;
		}

		for (int j=0; j<50; j++) {
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
    rank = MPI::COMM_WORLD.Get_rank();
    #endif

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
					double dE = -energy(spin1,spin2);
					if (dim==1) {
						for (int i=-1; i<2; i++) {
							x[0] += i;
							int spin = mcGrid(x);
							dE += energy(spin,spin2)-energy(spin,spin1);
							x[0] -= i;
						}
					} else if (dim==2) {
						for (int i=-1; i<2; i++) {
							x[0] += i;
							for (int j=-1; j<2; j++){
								x[1] += j;
								int spin = mcGrid(x);
								dE += energy(spin,spin2)-energy(spin,spin1);
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
									dE += energy(spin,spin2)-energy(spin,spin1);
									x[2] -= k;
								}
								x[1] -= j;
							}
							x[0] -= i;
						}
					}

					// compute boundary energy, mobility
					double E = energy(spin1,spin2);
					double M = mobility(spin1,spin2);

					// attempt a spin flip
					double r = double(rand())/double(RAND_MAX);
					if (dE<=0.0 and r<M*E) mcGrid(p) = spin2;
					if (dE>0.0 and r<M*E*exp(-dE/(E*kT))) mcGrid(p) = spin2;
				}
			}
			if (h%gss==0) ghostswap(mcGrid);
		}
	}
}

} // namespace MC

#endif

#include"MMSP.main.hpp"
