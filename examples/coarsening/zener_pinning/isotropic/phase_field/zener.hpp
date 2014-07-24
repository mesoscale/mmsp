// zener.hpp
// Algorithms for 2D and 3D isotropic phase field growth with Zener pinning
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ZENER_UPDATE
#define ZENER_UPDATE
#include"MMSP.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		MMSP::grid<1,vector<double> > grid(3,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(i)[0] = 0.0;
			grid(i)[1] = 0.0;
			double d = 64.0-x[0];
			if (d<32.0) grid(i)[2] = 1.0;
			else grid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++) {
				grid[x][0] = 1.0;
				grid[x][1] = 0.0;
				grid[x][2] = 0.0;
			}
		}

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(3,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(i)[0] = 0.0;
			grid(i)[1] = 0.0;
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) grid(i)[2] = 1.0;
			else grid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++) {
					grid[x][y][0] = 1.0;
					grid[x][y][1] = 0.0;
					grid[x][y][2] = 0.0;
				}
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,vector<double> > grid(3,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			grid(i)[0] = 0.0;
			grid(i)[1] = 0.0;
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) grid(i)[2] = 1.0;
			else grid(i)[1] = 1.0;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					for (int z=p[1]-1; z<=p[2]+1; z++) {
						grid[x][y][z][0] = 1.0;
						grid[x][y][z][1] = 0.0;
						grid[x][y][z][2] = 0.0;
					}
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
			// compute laplacians
			vector<double> lap = laplacian(grid,i);

			// compute sums of squares
			double sum = 0.0;
			for (int j=0; j<fields(grid); j++) {
				double phi = grid(i)[j];
				sum += phi*phi;
			}

			// compute update values
			for (int j=0; j<fields(grid); j++) {
				double phi = grid(i)[j];
				// particles have zero mobility
				if (j==0) update(i)[j] = phi;
				else update(i)[j] = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[j]));
			}
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
