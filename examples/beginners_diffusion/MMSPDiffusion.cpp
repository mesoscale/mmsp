#include"MMSP.hpp"
using namespace MMSP;

// Add in Laplacian function that is built in to MMSP
// Try using a grid of vectors >> one iterated with this code, one with laplacian


//we start the program off the same way as before, but this time we do not need the offset length variable
int main(int argc, char* argv[])
{
	Init(argc, argv);
	int length;
	int iterate;

	std::cout<<"input couple length"<<std::endl;
	std::cin>>length;
	std::cout<<""<<std::endl;

	std::cout<<"input number of iterations"<<std::endl;
	std::cin>>iterate;
	std::cout<<""<<std::endl;

	grid<1,scalar<double> > GRID(1,0,length);
	grid<1,scalar<double> > update(1,0,length);

	for (int x=x0(GRID); x<x1(GRID); x++)
		if (x<length/2) {
			GRID[x]=1;
			update[x]=1;
		} else {
			GRID[x]=0;
			update[x]=1;
		}

	//now we set the boundary conditions of both grids.  By choosing the Dirichlet conditions, it is nearly identical to the manually set boundaires.
	//the difference is that the first and last nodes of the grid can change, and the theoretical points outside the grid are fixed.
	b0(GRID,0) = Dirichlet;
	b1(GRID,0) = Dirichlet;
	b0(update,0) = Dirichlet;
	b1(update,0) = Dirichlet;

	for (int k=0; k<iterate; k++) {
		for (int i=0; i<nodes(GRID); i++) {
			//we can use MMSP's definition for laplacian instead of hard coding it.
			update(i)=0.5*laplacian(GRID,i)+GRID[i];
		}
		swap(GRID,update);
		ghostswap(GRID);
	};


	for (int x=x0(GRID); x<x1(GRID); x++)
		std::cout<<GRID[x]<<std::endl;

	Finalize();
	return 0;
}

