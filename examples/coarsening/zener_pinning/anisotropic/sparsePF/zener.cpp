// zener.cpp
// Anisotropic coarsening algorithms for 2D and 3D sparse phase field (sparsePF) methods
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
	double width = 8.0;
	double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		grid<dim,sparse<double> > newGrid(oldGrid);

		for (int n=0; n<nodes(oldGrid); n++) {
			// determine nonzero fields within 
			// the neighborhood of this node
			sparse<int> s;
			vector<int> x = position(oldGrid,n);

			for (int j=0; j<dim; j++)
				for (int k=-1; k<=1; k++) {
					x[j] += k;
					for (int h=0; h<length(oldGrid(x)); h++) {
						int i = index(oldGrid(x),h);
						set(s,i) = 1;
					}
					x[j] -= k;
				}
			double S = double(length(s));

			// if only one field is nonzero,
			// then copy this node to newGrid
			if (S<2.0) newGrid(n) = oldGrid(n);

			else {
				// compute laplacian of each field
				sparse<double> lap = laplacian(oldGrid,n);

				// compute variational derivatives
				sparse<double> dFdp;
				for (int h=0; h<length(s); h++) {
					int hindex = index(s,h);
					for (int j=h+1; j<length(s); j++) {
						int jindex = index(s,j);
						double gamma = energy(hindex,jindex);
						double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
						double w = 4.0*gamma/width;
						set(dFdp,hindex) += 0.5*eps*eps*lap[jindex]+w*oldGrid(n)[jindex];
						set(dFdp,jindex) += 0.5*eps*eps*lap[hindex]+w*oldGrid(n)[hindex];
						for (int k=j+1; k<length(s); k++) {
							int kindex = index(s,k);
							set(dFdp,hindex) += 3.0*oldGrid(n)[jindex]*oldGrid(n)[kindex];
							set(dFdp,jindex) += 3.0*oldGrid(n)[kindex]*oldGrid(n)[hindex];
							set(dFdp,kindex) += 3.0*oldGrid(n)[hindex]*oldGrid(n)[jindex];
						}
					}
				}

				// compute time derivatives
				sparse<double> dpdt;
				for (int h=0; h<length(s); h++) {
					int hindex = index(s,h);
					for (int j=h+1; j<length(s); j++) {
						int jindex = index(s,j);
						double mu = mobility(hindex,jindex);
						// set mobility of particles to zero
						if (hindex==0 or jindex==0) mu = 0.0;
						set(dpdt,hindex) -= mu*(dFdp[hindex]-dFdp[jindex]);
						set(dpdt,jindex) -= mu*(dFdp[jindex]-dFdp[hindex]);
					}
				}

				// compute update values
				double sum = 0.0;
				for (int h=0; h<length(s); h++) {
					int i = index(s,h);
					double value = oldGrid(n)[i]+dt*(2.0/S)*dpdt[i];
					if (value>1.0) value = 1.0;
					if (value<0.0) value = 0.0;
					if (value>epsilon) set(newGrid(n),i) = value;
					sum += newGrid(n)[i];
				}

				// project onto Gibbs simplex
				double rsum = 0.0;
				if (fabs(sum)>0.0) rsum = 1.0/sum;
				for (int h=0; h<length(newGrid(n)); h++) {
					int i = index(newGrid(n),h);
					set(newGrid(n),i) *= rsum;
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
