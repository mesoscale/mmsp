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
		int L=1024;
		GRID1D initGrid(4,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			for (int h=0; h<fields(initGrid); h++)
				initGrid(i)[h] = 0.0;

			vector<int> x = position(initGrid,i);

			if (x[0]%128 < 32)
				initGrid(i)[3] = 1.0;
			else if (x[0]%128 < 64)
				initGrid(i)[2] = 1.0;
			else if (x[0]%128 < 96)
				initGrid(i)[1] = 1.0;
			else
				initGrid(i)[0] = 1.0;
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(4,0,L,0,L);

		// Divide domain into "unit cells", 128 points on an edge
		for (int i=0; i<nodes(initGrid); i++) {
			for (int h=0; h<fields(initGrid); h++)
				initGrid(i)[h] = 0.0;

			vector<int> x = position(initGrid,i);

			if (x[0]%128 < 32 || x[0]%128 > 96) { // less than 1/4, more than 3/4
				if ((x[1]%128) < 64) // less than 1/2
					initGrid(i)[2] = 1.0;
				else
					initGrid(i)[3] = 1.0;
			} else {
				if (x[1]%128 < 32 || x[1]%128 > 96) // less than 1/4, more than 3/4
					initGrid(i)[1] = 1.0;
				else
					initGrid(i)[0] = 1.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=64;
		GRID3D initGrid(4,0,2*L,0,L,0,L/4);

		// Divide domain into "unit cells", 64 points on an edge
		for (int i=0; i<nodes(initGrid); i++) {
			for (int h=0; h<fields(initGrid); h++)
				initGrid(i)[h] = 0.0;

			vector<int> x = position(initGrid,i);

			if (x[0]%64 < 16 || x[0]%64 > 48) { // less than 1/4, more than 3/4
				if ((x[1]%64) < 32) // less than 1/2
					initGrid(i)[2] = 1.0;
				else
					initGrid(i)[3] = 1.0;
			} else {
				if (x[1]%64 < 16 || x[1]%64 > 48) // less than 1/4, more than 3/4
					initGrid(i)[1] = 1.0;
				else
					initGrid(i)[0] = 1.0;
			}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	ghostswap(oldGrid);

	grid<dim,vector<T> > newGrid(oldGrid);

	double dt = 0.01;
	double width = 8.0;

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step,steps);
		for (int i=0; i<nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid,i);

			// determine nonzero fields within
			// the neighborhood of this node
			T S = 0.0;
			vector<int> s(fields(oldGrid),0);
			for (int h=0; h<fields(oldGrid); h++) {
				for (int j=0; j<dim; j++)
					for (int k=-1; k<=1; k++) {
						x[j] += k;
						if (oldGrid(x)[h]>0.0) {
							s[h] = 1;
							x[j] -= k;
							goto next;
						}
						x[j] -= k;
					}
				next: S += s[h];
			}

			// if only one field is nonzero,
			// then copy this node to newGrid
			if (S<2.0) newGrid(i) = oldGrid(i);

			else {
				// compute laplacian of each field
				vector<T> lap = laplacian(oldGrid,i);

				// compute variational derivatives
				vector<T> dFdp(fields(oldGrid),0.0);
				for (int h=0; h<fields(oldGrid); h++)
					if (s[h]>0.0)
						for (int j=h+1; j<fields(oldGrid); j++)
							if (s[j]>0.0) {
								double gamma = energy(h,j);
								double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
								double w = 4.0*gamma/width;
								dFdp[h] += 0.5*eps*eps*lap[j]+w*oldGrid(i)[j];
								dFdp[j] += 0.5*eps*eps*lap[h]+w*oldGrid(i)[h];
								for (int k=j+1; k<fields(oldGrid); k++)
									if (s[k]>0.0) {
										dFdp[h] += 3.0*oldGrid(i)[j]*oldGrid(i)[k];
										dFdp[j] += 3.0*oldGrid(i)[k]*oldGrid(i)[h];
										dFdp[k] += 3.0*oldGrid(i)[h]*oldGrid(i)[j];
									}
							}

				// compute time derivatives
				vector<T> dpdt(fields(oldGrid),0.0);
				for (int h=0; h<fields(oldGrid); h++)
					if (s[h]>0.0)
						for (int j=h+1; j<fields(oldGrid); j++)
							if (s[j]>0.0) {
								double mu = mobility(h,j);
								dpdt[h] -= mu*(dFdp[h]-dFdp[j]);
								dpdt[j] -= mu*(dFdp[j]-dFdp[h]);
							}

				// compute newGrid values
				T sum = 0.0;
				for (int h=0; h<fields(oldGrid); h++) {
					T value = oldGrid(i)[h]+dt*(2.0/S)*dpdt[h];
					if (value>1.0) value = 1.0;
					if (value<0.0) value = 0.0;
					newGrid(i)[h] = value;
					sum += value;
				}

				// project onto Gibbs simplex
				T rsum = 0.0;
				if (fabs(sum)>0.0) rsum = 1.0/sum;
				for (int h=0; h<fields(oldGrid); h++)
					newGrid(i)[h] *= rsum;
			}
		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
