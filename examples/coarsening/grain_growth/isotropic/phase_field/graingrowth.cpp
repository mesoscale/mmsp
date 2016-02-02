// graingrowth.cpp
// Algorithms for 2D and 3D isotropic phase field grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"graingrowth.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		GRID1D initGrid(2,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = 64.0-x[0];
			if (d<32.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(2,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(2,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	grid<dim,vector<T> > newGrid(oldGrid);

	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(oldGrid); i++) {
			// compute laplacian
			vector<T> lap = laplacian(oldGrid,i);

			// compute sum of squares
			double sum = 0.0;
			for (int j=0; j<fields(oldGrid); j++) {
				double phi = oldGrid(i)[j];
				sum += phi*phi;
			}

			// compute update values
			for (int j=0; j<fields(oldGrid); j++) {
				T phi = oldGrid(i)[j];
				newGrid(i)[j] = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[j]));
			}
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
