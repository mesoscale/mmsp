// graingrowth.hpp
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
	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(2,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(x)[0] = 0.0;
			grid(x)[1] = 0.0;
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) grid(x)[1] = 1.0;
			else grid(x)[0] = 1.0;
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,vector<double> > grid(2,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(x)[0] = 0.0;
			grid(x)[1] = 0.0;
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) grid(x)[1] = 1.0;
			else grid(x)[0] = 1.0;
		}

		output(grid,filename);
	}
}

template <int dim> void update(MMSP::grid<dim,vector<double> >& grid, int steps)
{
	MMSP::grid<dim,vector<double> > update(grid);

	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			// compute laplacian
			vector<double> lap = laplacian(grid,i);

			// compute sum of squares
			double sum = 0.0;
			for (int j=0; j<fields(grid); j++) {
				double phi = grid(i)[j];
				sum += phi*phi;
			}

			// compute update values
			for (int j=0; j<fields(grid); j++) {
				double phi = grid(i)[j];
				update(i)[j] = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[j]));
			}
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
