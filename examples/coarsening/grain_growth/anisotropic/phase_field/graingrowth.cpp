// graingrowth.cpp
// Anisotropic coarsening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include"anisotropy.hpp"
#include"graingrowth.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		MMSP::grid<1,vector<double> > grid(4,0,128);

		for (int i=0; i<nodes(grid); i++) {
			for (int h=0; h<fields(grid); h++)
				grid(i)[h] = 0.0;

			vector<int> x = position(grid,i);

			if (x[0]<32)      grid(i)[3] = 1.0;
			else if (x[0]>96) grid(i)[3] = 1.0;
			else              grid(i)[0] = 1.0;
		}

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(4,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			for (int h=0; h<fields(grid); h++)
				grid(i)[h] = 0.0;

			vector<int> x = position(grid,i);

			if (x[0]<32) {
				if (x[1]<64) grid(i)[2] = 1.0;
				else grid(i)[3] = 1.0;
			}
			else if (x[0]>96) {
				if (x[1]<64) grid(i)[2] = 1.0;
				else grid(i)[3] = 1.0;
			}
			else {
				if (x[1]<32 or x[1]>96) grid(i)[1] = 1.0;
				else grid(i)[0] = 1.0;
			}
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<2,vector<double> > grid(4,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			for (int h=0; h<fields(grid); h++)
				grid(i)[h] = 0.0;

			vector<int> x = position(grid,i);

			if (x[0]<16) {
				if (x[1]<32) grid(i)[2] = 1.0;
				else grid(i)[3] = 1.0;
			}
			else if (x[0]>48) {
				if (x[1]<32) grid(i)[2] = 1.0;
				else grid(i)[3] = 1.0;
			}
			else {
				if (x[1]<16 or x[1]>48) grid(i)[1] = 1.0;
				else grid(i)[0] = 1.0;
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
			vector<int> x = position(grid,i);

			// determine nonzero fields within 
			// the neighborhood of this node
			double S = 0.0;
			vector<int> s(fields(grid),0);
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
