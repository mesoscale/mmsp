#include"MMSP.hpp"
using namespace MMSP;


// Add in Laplacian funciton that is built in to MMSP
// Try using a grid of vectors >> one iterated with this code, one with laplacian
// Make it more flexible for documentation (get the n-dimension thing down)


int main()
{
	//Here is where all variables are called.  In c++ all variables must be defined before the can be used.
	//The first word is the data type associated with the vairable.
	int length;
	int offlength;
	int iterate;

	//Here we ask for user input, the statement cin>>variable saves whatever is entered by the user in variable
	//We make it so that the length of the diffusion couple and the number of iterations can be changed each time the code is run
	std::cout<<"input couple length"<<std::endl;
	std::cin>>length;
	std::cout<<""<<std::endl;

	std::cout<<"input number of iterations"<<std::endl;
	std::cin>>iterate;
	std::cout<<""<<std::endl;

	//here we define some 1 dimensional grids with float variable types
	grid<1,scalar<float> > GRID(1,0,length);
	grid<1,scalar<float> > GRID2(1,0,length);
	//this value is defined for looping control
	offlength=x1(GRID)-3;

	//this creates two identical grids GRID and GRID2 that are 1 for the first half and 0 for the second. These represent diffusion couples.
	for (int x=x0(GRID); x<x1(GRID); x++)
		if (x<length/2) {
			GRID[x]=1;
			GRID2[x]=1;
		} else {
			GRID[x]=0;
			GRID2[x]=0;
		}

	//This step controls the number of time steps based on the user input from before
	for (int i=0; i<iterate; i++) {
		//Iterate through grid
		for (int x=x0(GRID); x<x1(GRID); x++) {
			//Define fixed boundaries by preventing the first and last nodes of the grid from changing
			if (x==0 || x==length-1) {
			}
			//Take one time step of the discrete Fick's Law with maximum stability criterion (.5)
			//to keep calculations from interfering with each other, the results of computations on GRID are stored in GRID2, then copied back to GRID after the last computation
			else {
				GRID2[x]=0.5*(GRID[x-1]-2*GRID[x]+GRID[x+1])+GRID[x];
				if (x>offlength) {
					swap(GRID,GRID2);
				}
			}
		}
	}
	//This prints the results of the grid to cout
	for (int x=x0(GRID); x<x1(GRID); x++) {
		std::cout<<GRID[x]<<std::endl;
	}
	Finalize();
}

