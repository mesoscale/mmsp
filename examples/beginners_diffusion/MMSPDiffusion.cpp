#include"MMSP.hpp"
using namespace MMSP;

int main(int argc, char* argv[])
{
	Init(argc, argv);

	int nx = 10;
	double dx = 1.0;
	int iterations = 100;
	double diffusionCoefficient = 1.0;

	double dt = dx*dx/diffusionCoefficient/4.0;

	grid<1,scalar<double> > newGrid(1,0,nx);
	grid<1,scalar<double> > oldGrid(1,0,nx);

	for (int x=x0(newGrid); x<x1(newGrid); x++)
		if (x<nx/2) {
			newGrid[x]=1;
			oldGrid[x]=1;
		} else {
			newGrid[x]=0;
			oldGrid[x]=0;
		}

	if (x0(newGrid)==g0(newGrid, 0)) {
		b0(newGrid, 0) = Dirichlet;
		b0(oldGrid, 0) = Dirichlet;
	}
	if (x1(newGrid)==g1(newGrid, 0)) {
		b1(newGrid, 0) = Dirichlet;
		b1(oldGrid, 0) = Dirichlet;
	}

	ghostswap(oldGrid);

	for (int k=0; k<iterations; k++) {
		for (int i=0; i<nodes(newGrid); i++) {
			newGrid(i) = (diffusionCoefficient*dt)*laplacian(oldGrid,i) + oldGrid[i];
			oldGrid[x0(newGrid)] = 1.0;
			oldGrid[x1(newGrid)] = 0.0;
		}
		swap(newGrid,oldGrid);
		ghostswap(newGrid);
	}


	for (int x=x0(newGrid); x<x1(newGrid); x++)
		std::cout<<newGrid[x]<<std::endl;

	Finalize();
	return 0;
}

