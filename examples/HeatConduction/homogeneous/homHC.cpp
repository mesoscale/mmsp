// hetHC.cpp
// Algorithms for 2D homogeneous heat conduction
// Questions/comments to thomas.flint@manchester.ac.uk (Tom Flint)

#ifndef HETHC_UPDATE
#define HETHC_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"homHC.hpp"

namespace MMSP{



void generate(int dim, const char* filename)
{
	if (dim==1) {
	std::cerr << "Heat Conduction code is only implemented for 2D." << std::endl;
		MMSP::Abort(-1);
	}

	if (dim==2) {
		int L=256;
		GRID2D initGrid(0,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++)
		//vector<int> x = position(initGrid,i);
		//if(x[0]<=1){
			//initGrid(i)=1.0;
		//	}
		//else{
			initGrid(i) = 0.0;
		//}


		output(initGrid,filename);
	}

	if (dim==3) {
std::cerr << "Heat Conduction code is only implemented for 2D." << std::endl;
		MMSP::Abort(-1);
	}
}

template <int dim, typename T> void update(grid<dim,T>& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    rank = MPI::COMM_WORLD.Get_rank();
    #endif

	ghostswap(oldGrid);


 MMSP::grid<2,double > GRID_Diffusivity(1,0,256,0,256);     //for storing phase data
//dx(GRID_Diffusivity,0)=deltaX;
//dx(GRID_Diffusivity,1)=deltaY;



	grid<dim, double> newGrid(oldGrid);


	double dt = 1e-2;


	for (int step=0; step<steps; step++) {
		if (rank==0)
			print_progress(step, steps);

		for (int i=0; i<nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid,i);
			GRID_Diffusivity(i)=10.0;

			if(x[0]<=1){
			oldGrid(i)=1.0;
			}

		}
		MMSP::output(GRID_Diffusivity, "alpha");


		for (int i=0; i<nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid,i);
			double old = oldGrid(i);

			vector<double> GRID_gradD = gradnew(GRID_Diffusivity,x);
			vector<double> GRID_gradT = gradnew(oldGrid,x);
			
			double lapT=laplacian(oldGrid,i);
			
			double dTdt=(GRID_Diffusivity(i)*lapT);

			newGrid(i) = old+(dt*dTdt);
			
			if(x[0]<=1){
			newGrid(i)=1.0;
			}
			if(x[0]>=250){
			newGrid(i)=oldGrid(i);
			}

		}
		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}

}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
