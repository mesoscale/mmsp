#include"MMSP.hpp"
using namespace MMSP;

int main(int argc, char* argv[])
{
	Init(argc, argv);

	int nx;
	int iterations;
	float dt, diffusionCoefficient, dx;

	nx = 10;
	iterations = 100;
	diffusionCoefficient = 1.0;
	dx = 1.0;
	dt = dx*dx/diffusionCoefficient/4;

	grid<1,scalar<double> > GRID(1,0,nx);
	grid<1,scalar<double> > update(1,0,nx);

	for (int x=x0(GRID); x<x1(GRID); x++)
		if (x<nx/2) {
			GRID[x]=1;
			update[x]=1;
		} else {
			GRID[x]=0;
			update[x]=1;
		}

	b0(GRID,0) = Dirichlet;
	b1(GRID,0) = Dirichlet;
	b0(update,0) = Dirichlet;
	b1(update,0) = Dirichlet;

	ghostswap(GRID);

	for (int k=0; k<iterations; k++) {
		for (int i=0; i<nodes(GRID); i++) {
			update(i)=(diffusionCoefficient*dt/dx/dx)*laplacian(GRID,i)+GRID[i];
			update[x0(GRID)] = 1.0;
			update[x1(GRID)] = 0.0;
		}
		swap(GRID,update);
		ghostswap(GRID);
	};


	for (int x=x0(GRID); x<x1(GRID); x++)
		std::cout<<GRID[x]<<std::endl;

	Finalize();
	return 0;
}

