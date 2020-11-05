// heisenberg.cpp
// 2D and 3D heisenberg model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef HEISENBERG_UPDATE
#define HEISENBERG_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"heisenberg.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.
	if (dim==1) {
		int L=1024;
		GRID1D initGrid(3,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			initGrid(i)[0] = cos(psi)*sin(theta);
			initGrid(i)[1] = sin(psi)*sin(theta);
			initGrid(i)[2] = cos(theta);
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(3,0,2*L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			initGrid(i)[0] = cos(psi)*sin(theta);
			initGrid(i)[1] = sin(psi)*sin(theta);
			initGrid(i)[2] = cos(theta);
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=32;
		GRID3D initGrid(3,0,2*L,0,L,0,L/4);

		for (int i=0; i<nodes(initGrid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			initGrid(i)[0] = cos(psi)*sin(theta);
			initGrid(i)[1] = sin(psi)*sin(theta);
			initGrid(i)[2] = cos(theta);
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& spinGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

	ghostswap(spinGrid);

	double J = 1.0;
	double kT = (dim==3)?0.75:0.50;
	double pi = acos(-1.0);

	int gss = (dim==1)?nodes(spinGrid):int(sqrt(nodes(spinGrid)));
	// srand() is called exactly once in MMSP.main.hpp. Do not call it here.

	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int h=0; h<nodes(spinGrid); h++) {
			// choose a random site
			int p = rand()%nodes(spinGrid);
			vector<int> x = position(spinGrid,p);
			vector<T>& s1 = spinGrid(p);

			// choose a random unit vector
			vector<T> s2(3);
			double psi = 2.0*pi*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			s2[0] = cos(psi)*sin(theta);
			s2[1] = sin(psi)*sin(theta);
			s2[2] = cos(theta);

			// compute energy change
			double sum = -1.0;
			if (dim==1) {
				for (int i=-1; i<2; i++) {
					x[0] += i;
					vector<T>& s = spinGrid(x);
					sum += s[0]*(s1[0]-s2[0])+s[1]*(s1[1]-s2[1])+s[2]*(s1[2]-s2[2]);
					x[0] -= i;
				}
			} else if (dim==2) {
				for (int i=-1; i<2; i++) {
					x[0] += i;
					for (int j=-1; j<2; j++) {
						x[1] += j;
						vector<T>& s = spinGrid(x);
						sum += s[0]*(s1[0]-s2[0])+s[1]*(s1[1]-s2[1])+s[2]*(s1[2]-s2[2]);
						x[1] -= j;
					}
					x[0] -= i;
				}
			} else if (dim==3) {
				for (int i=-1; i<2; i++) {
					x[0] += i;
					for (int j=-1; j<2; j++) {
						x[1] += j;
						for (int k=-1; k<2; k++) {
							x[2] += k;
							vector<T>& s = spinGrid(x);
							sum += s[0]*(s1[0]-s2[0])+s[1]*(s1[1]-s2[1])+s[2]*(s1[2]-s2[2]);
							x[2] -= k;
						}
						x[1] -= j;
					}
					x[0] -= i;
				}
			}
			double dE = -J*sum;

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0)
				spinGrid(p) = s2;
			else if (r<exp(-dE/kT))
				spinGrid(p) = s2;

			if (h%gss==0)
				ghostswap(spinGrid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
