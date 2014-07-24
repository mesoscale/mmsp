// graingrowth.hpp
// Algorithms for 2D and 3D isotropic sparsePF grain growth
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
		MMSP::grid<1,sparse<double> > grid(0,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = 64.0-x[0];
			if (d<32.0) set(grid(i),1)= 1.0;
			else set(grid(i),0) = 1.0;
		}

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,sparse<double> > grid(0,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) set(grid(i),1)= 1.0;
			else set(grid(i),0) = 1.0;
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,sparse<double> > grid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) set(grid(i),1)= 1.0;
			else set(grid(i),0) = 1.0;
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

			// if there is only one nonzero
			// field, it remains the same
			if (length(neighbors)<2)
				update(i) = grid(i);

			else {
				// compute laplacian
				sparse<double> lap = laplacian(grid,i);

				// compute sum of squares
				double sum = 0.0;
				for (int j=0; j<length(grid(i)); j++) {
					double phi = MMSP::value(grid(i),j);
					sum += phi*phi;
				}

				// compute update values
				for (int j=0; j<length(neighbors); j++) {
					int index = MMSP::index(neighbors,j);
					double phi = grid(i)[index];
					double value = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[index]));
					if (value>epsilon) set(update(i),index) = value;
				}
			}
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
