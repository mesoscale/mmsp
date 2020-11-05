// zener.cpp
// Algorithms for 2D and 3D isotropic phase field growth with Zener pinning
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ZENER_UPDATE
#define ZENER_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"zener.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(3,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = 32-x[0]%64;
			if (d<16.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++) {
				initGrid(p)[0] = 1.0;
				initGrid(p)[1] = 0.0;
				initGrid(p)[2] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(3,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (d<16.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++) {
					initGrid(p)[0] = 1.0;
					initGrid(p)[1] = 0.0;
					initGrid(p)[2] = 0.0;
				}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(3,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (d<16.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++)
					for (p[2]=x[2]-1; p[2]<=x[2]+1; p[2]++) {
						initGrid(p)[0] = 1.0;
						initGrid(p)[1] = 0.0;
						initGrid(p)[2] = 0.0;
					}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);

	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int i=0; i<nodes(oldGrid); i++) {
			// compute laplacians
			vector<T> lap = laplacian(oldGrid,i);

			// compute sums of squares
			double sum = 0.0;
			for (int j=0; j<fields(oldGrid); j++) {
				T phi = oldGrid(i)[j];
				sum += phi*phi;
			}

			// compute update values
			for (int j=0; j<fields(oldGrid); j++) {
				T phi = oldGrid(i)[j];
				// particles have zero mobility
				if (j==0) newGrid(i)[j] = phi;
				else newGrid(i)[j] = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[j]));
			}
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

