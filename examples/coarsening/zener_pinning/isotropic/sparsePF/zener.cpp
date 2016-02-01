// graingrowth.hpp
// Algorithms for 2D and 3D isotropic sparsePF grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"zener.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		grid<1,sparse<double> > initGrid(0,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = 64.0-x[0];
			if (d<32.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++) {
				set(initGrid[x],0) = 1.0;
				set(initGrid[x],1) = 0.0;
				set(initGrid[x],2) = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		grid<2,sparse<double> > initGrid(0,0,128,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++) {
					set(initGrid[x][y],0) = 1.0;
					set(initGrid[x][y],1) = 0.0;
					set(initGrid[x][y],2) = 0.0;
				}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		grid<3,sparse<double> > initGrid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) set(initGrid(i),2) = 1.0;
			else set(initGrid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(initGrid);
			vector<int> p = position(initGrid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					for (int z=p[1]-1; z<=p[2]+1; z++) {
						set(initGrid[x][y][z],0) = 1.0;
						set(initGrid[x][y][z],1) = 0.0;
						set(initGrid[x][y][z],2) = 0.0;
					}
		}

		output(initGrid,filename);
	}
}

template <int dim> void update(grid<dim,sparse<double> >& oldGrid, int steps)
{
	double dt = 0.01;
	double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		// update grid must be overwritten each time
		grid<dim,sparse<double> > newGrid(oldGrid);

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
				sparse<double> lap = laplacian(oldGrid,n);

				// compute sums of squares
				double sum = 0.0;
				for (int j=0; j<length(oldGrid(n)); j++) {
					double phi = value(oldGrid(n),j);
					sum += phi*phi;
				}

				// compute update values
				for (int j=0; j<length(neighbors); j++) {
					int i = index(neighbors,j);
					double phi = oldGrid(n)[i];
					// particles have zero mobility
					if (n==0) set(newGrid(n),i) = phi;
					else {
						double value = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[i]));
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

