// graingrowth.hpp
// Algorithms for 2D and 3D isotropic sparsePF grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		MMSP::grid<1,sparse<double> > grid(0,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = 64.0-x[0];
			if (d<32.0) set(grid(i),2) = 1.0;
			else set(grid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++) {
				set(grid[x],0) = 1.0;
				set(grid[x],1) = 0.0;
				set(grid[x],2) = 0.0;
			}
		}

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,sparse<double> > grid(0,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) set(grid(i),2) = 1.0;
			else set(grid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++) {
					set(grid[x][y],0) = 1.0;
					set(grid[x][y],1) = 0.0;
					set(grid[x][y],2) = 0.0;
				}
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,sparse<double> > grid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) set(grid(i),2) = 1.0;
			else set(grid(i),1) = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					for (int z=p[1]-1; z<=p[2]+1; z++) {
						set(grid[x][y][z],0) = 1.0;
						set(grid[x][y][z],1) = 0.0;
						set(grid[x][y][z],2) = 0.0;
					}
		}

		output(grid,filename);
	}
}

template <int dim> void update(MMSP::grid<dim,sparse<double> >& grid, int steps)
{
	double dt = 0.01;
	double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		// update grid must be overwritten each time
		MMSP::grid<dim,sparse<double> > update(grid);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);

			// determine nonzero fields within
			// the neighborhood of this node
			sparse<bool> neighbors;
			for (int j=0; j<dim; j++)
				for (int k=-1; k<=1; k++) {
					x[j] += k;
					for (int h=0; h<length(grid(x)); h++) {
						int index = MMSP::index(grid(x),h);
						set(neighbors,index) = 1;
					}
					x[j] -= k;
				}

			// if there is only one nonzero field,
			// then it remains the same 
			if (length(neighbors)<2)
				update(i) = grid(i);

			else {
				// compute laplacians
				sparse<double> lap = laplacian(grid,i);

				// compute sums of squares
				double sum = 0.0;
				for (int j=0; j<length(grid(i)); j++) {
					double phi = MMSP::value(grid(i),j);
					sum += phi*phi;
				}

				// compute update values
				for (int j=0; j<length(neighbors); j++) {
					int index = MMSP::index(neighbors,j);
					double phi = grid(i)[index];
					// particles have zero mobility
					if (index==0) set(update(i),index) = phi;
					else {
						double value = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[index]));
						if (value>epsilon) set(update(i),index) = value;
					}
				}
			}
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
