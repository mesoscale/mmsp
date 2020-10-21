// zener.cpp
// Algorithms for 2D and 3D isotropic sparsePF growth with Zener pinning
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
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
			if (d<16.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++) {
				set(initGrid(p),0) = 1.0;
				set(initGrid(p),1) = 0.0;
				set(initGrid(p),2) = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=128;
		GRID2D initGrid(0,0,2*L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (d<16.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++) {
					set(initGrid(p),0) = 1.0;
					set(initGrid(p),1) = 0.0;
					set(initGrid(p),2) = 0.0;
				}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(0,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32-x[0]%64,2)+pow(32-x[1]%64,2));
			if (d<16.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		int nParticles=std::max(50,(50*8192)/nodes(initGrid));
		for (int j=0; j<nParticles; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> x = position(initGrid,i);
			vector<int> p(x);
			for (p[0]=x[0]-1; p[0]<=x[0]+1; p[0]++)
				for (p[1]=x[1]-1; p[1]<=x[1]+1; p[1]++)
					for (p[2]=x[2]-1; p[2]<=x[2]+1; p[2]++) {
						set(initGrid(p),0) = 1.0;
						set(initGrid(p),1) = 0.0;
						set(initGrid(p),2) = 0.0;
					}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,sparse<T> >& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(oldGrid);

	grid<dim,sparse<T> > newGrid(oldGrid);

	double dt = 0.01;
	double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid,n);

			// determine nonzero fields within
			// the neighborhood of this node
			sparse<bool> neighbors;
			for (int j=0; j<dim; j++)
				for (int k=-1; k<=1; k++) {
					x[j] += k;
					for (int h=0; h<length(oldGrid(x)); h++) {
						int i = index(oldGrid(x),h);
						set(neighbors,i) = 1;
					}
					x[j] -= k;
				}

			// if there is only one nonzero field,
			// then it remains the same
			if (length(neighbors)<2)
				newGrid(n) = oldGrid(n);

			else {
				// compute laplacians
				sparse<T> lap = laplacian(oldGrid,n);

				// compute sums of squares
				double sum = 0.0;
				for (int j=0; j<length(oldGrid(n)); j++) {
					T phi = value(oldGrid(n),j);
					sum += phi*phi;
				}

				// compute update values
				for (int j=0; j<length(neighbors); j++) {
					int i = index(neighbors,j);
					T phi = oldGrid(n)[i];
					// particles have zero mobility
					if (n==0) set(newGrid(n),i) = phi;
					else {
						T value = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[i]));
						if (value>epsilon) set(newGrid(n),i) = value;
					}
				}
			}
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

