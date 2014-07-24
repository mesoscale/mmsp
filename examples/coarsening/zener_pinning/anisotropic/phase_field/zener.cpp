// zener.hpp
// Anisotropic coarsening algorithms for 2D and 3D phase field methods
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
	double width = 8.0;

	for (int step=0; step<steps; step++) {
		for (int i=0; i<nodes(grid); i++) {
			// determine nonzero fields within 
			// the neighborhood of this node
			double S = 0.0;
			vector<int> s(fields(grid),0);
			vector<int> x = position(grid,i);

			for (int h=0; h<fields(grid); h++) {
				for (int j=0; j<dim; j++)
					for (int k=-1; k<=1; k++) {
						x[j] += k;
						if (grid(x)[h]>0.0) {
							s[h] = 1;
							x[j] -= k;
							goto next;
						}
						x[j] -= k;
					}
				next: S += s[h];
			}

			// if only one field is nonzero,
			// then copy this node to update
			if (S<2.0) update(i) = grid(i);

			else {
				// compute laplacian of each field
				vector<double> lap = laplacian(grid,i);

				// compute variational derivatives
				vector<double> dFdp(fields(grid),0.0);
				for (int h=0; h<fields(grid); h++)
					if (s[h]>0.0)
						for (int j=h+1; j<fields(grid); j++)
							if (s[j]>0.0) {
								double gamma = energy(h,j);
								double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
								double w = 4.0*gamma/width;
								dFdp[h] += 0.5*eps*eps*lap[j]+w*grid(i)[j];
								dFdp[j] += 0.5*eps*eps*lap[h]+w*grid(i)[h];
								for (int k=j+1; k<fields(grid); k++)
									if (s[k]>0.0) {
										dFdp[h] += 3.0*grid(i)[j]*grid(i)[k];
										dFdp[j] += 3.0*grid(i)[k]*grid(i)[h];
										dFdp[k] += 3.0*grid(i)[h]*grid(i)[j];
									}
							}

				// compute time derivatives
				vector<double> dpdt(fields(grid),0.0);
				for (int h=0; h<fields(grid); h++)
					if (s[h]>0.0)
						for (int j=h+1; j<fields(grid); j++)
							if (s[j]>0.0) {
								double mu = mobility(h,j);
								// set mobility of particles to zero
								if (h==0 or j==0) mu = 0.0;
								dpdt[h] -= mu*(dFdp[h]-dFdp[j]);
								dpdt[j] -= mu*(dFdp[j]-dFdp[h]);
							}

				// compute update values
				double sum = 0.0;
				for (int h=0; h<fields(grid); h++) {
					double value = grid(i)[h]+dt*(2.0/S)*dpdt[h];
					if (value>1.0) value = 1.0;
					if (value<0.0) value = 0.0;
					update(i)[h] = value;
					sum += value;
				}

				// project onto Gibbs simplex
				double rsum = 0.0;
				if (fabs(sum)>0.0) rsum = 1.0/sum;
				for (int h=0; h<fields(grid); h++)
					update(i)[h] *= rsum;
			}
		}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
