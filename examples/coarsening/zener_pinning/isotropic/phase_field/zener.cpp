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
	if (dim==1) {
		GRID1D initGrid(3,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = 64.0-x[0];
			if (d<32.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++) {
				initGrid[x][0] = 1.0;
				initGrid[x][1] = 0.0;
				initGrid[x][2] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		GRID2D initGrid(3,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++) {
					initGrid[x][y][0] = 1.0;
					initGrid[x][y][1] = 0.0;
					initGrid[x][y][2] = 0.0;
				}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(3,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			initGrid(i)[0] = 0.0;
			initGrid(i)[1] = 0.0;
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) initGrid(i)[2] = 1.0;
			else initGrid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					for (int z=p[1]-1; z<=p[2]+1; z++) {
						initGrid[x][y][z][0] = 1.0;
						initGrid[x][y][z][1] = 0.0;
						initGrid[x][y][z][2] = 0.0;
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

