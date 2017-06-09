#include <stdio.h>
#include "MMSP.hpp"
using namespace MMSP;

int main(int argc, char* argv[])
{
	Init(argc, argv);

	int nx = 10;
	double dx = 1.0;
	int iterations = 1000;
	double diffusionCoefficient = 1.0;

	// Choice of dt depends on stability criteria.  This
	// choice maintains stability.
	double dt = dx*dx/diffusionCoefficient/4;

	grid<1,scalar<double> > newGrid(0,0,nx);
	grid<1,scalar<double> > oldGrid(0,0,nx);

	for (int x=x0(newGrid); x<x1(newGrid); x++)
		if (x<nx/2) {
			newGrid[x]=1;
			oldGrid[x]=1;
		} else {
			newGrid[x]=0;
			oldGrid[x]=0;
		}

	for (int i=0; i<iterations; i++) {
		for (int x=x0(newGrid); x<x1(newGrid); x++) {
			if (x==0 || x==nx-1) {
			} else {
				double laplacian = (oldGrid[x-1]-2*oldGrid[x]+oldGrid[x+1])/(dx*dx);
				newGrid[x] = (diffusionCoefficient*dt)*laplacian + oldGrid[x];
			}
		}
		swap(oldGrid,newGrid);
	}

	//This prints the results of the grid to cout
	for (int x=x0(oldGrid); x<x1(oldGrid); x++) {
		std::cout<<oldGrid[x]<<std::endl;
	}

	Finalize();

	return 0;
}

