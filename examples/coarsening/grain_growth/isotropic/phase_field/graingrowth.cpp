// graingrowth.cpp
// Algorithms for 2D and 3D isotropic phase field grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"graingrowth.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		GRID1D initGrid(2,0,128);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = 64.0-x[0];
			if (d<32.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=128;
		GRID2D initGrid(2,0,L,0,L);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(L/2-x[0],2)+pow(L/2-x[1],2));
			if (d<L/4.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		GRID3D initGrid(2,0,64,0,64,0,64);

		for (int i=0; i<nodes(initGrid); i++) {
			vector<int> x = position(initGrid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) {
				initGrid(i)[0] = 0.0;
				initGrid(i)[1] = 1.0;
			} else {
				initGrid(i)[0] = 1.0;
				initGrid(i)[1] = 0.0;
			}
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,vector<T> >& oldGrid, int steps)
{
	int rank=0;
    #ifdef MPI_VERSION
    rank = MPI::COMM_WORLD.Get_rank();
    #endif

	ghostswap(oldGrid);

   	grid<dim,vector<T> > newGrid(oldGrid);

	double dt = 0.01;

	for (int step=0; step<steps; step++) {
		#ifndef UNITTESTING
		if (rank==0)
			print_progress(step, steps);
		#endif
		for (int i=0; i<nodes(oldGrid); i++) {
			// compute laplacian
			vector<T> lap = laplacian(oldGrid,i);

			// compute sum of squares
			double sum = 0.0;
			for (int j=0; j<fields(oldGrid); j++) {
				double phi = oldGrid(i)[j];
				sum += phi*phi;
			}

			// compute update values
			for (int j=0; j<fields(oldGrid); j++) {
				T phi = oldGrid(i)[j];
				newGrid(i)[j] = phi-dt*(-phi-pow(phi,3)+2.0*(phi*sum-lap[j]));
			}
		}
		ghostswap(oldGrid);
		swap(oldGrid,newGrid);
	}
}

} // namespace MMSP

#endif

#ifndef UNITTESTING

#include"MMSP.main.hpp"

#else

/* validation with Google Test suite
 * http://github.com/google/googletest
 */
#include<gtest/gtest.h>

TEST(grid1DTest, Finite) {
	MMSP::vector<double> length(2,0.0);
	MMSP::generate(1, "short1.00000.dat");
	GRID1D myGrid1D("short1.00000.dat");
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		length[0] += myGrid1D(n)[1];
	MMSP::update(myGrid1D, 20000);
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		length[1] += myGrid1D(n)[1];
	myGrid1D.output("short1.20000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<2; i++) {
		double myLen(length[i]);
		MPI::COMM_WORLD.Allreduce(&myLen,&length[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		for (int i=0; i<MMSP::fields(myGrid1D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid1D(n)[i]);
	EXPECT_NEAR(length[1],length[0],1.0);
}

TEST(grid2DTest, Finite) {
	MMSP::vector<double> area(2,0.0);
	MMSP::generate(2, "short2.0000.dat");
	GRID2D myGrid2D("short2.0000.dat");
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		area[0] += myGrid2D(n)[1];
	MMSP::update(myGrid2D, 4000);
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		area[1] += myGrid2D(n)[1];
	myGrid2D.output("short2.4000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<2; i++) {
		double myAre(area[i]);
		MPI::COMM_WORLD.Allreduce(&myAre,&area[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		for (int i=0; i<MMSP::fields(myGrid2D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid2D(n)[i]);
	EXPECT_LE(area[1],area[0]);
}

TEST(grid3DTest, Finite) {
	MMSP::vector<double> volume(2,0.0);
	MMSP::generate(3, "short3.0000.dat");
	GRID3D myGrid3D("short3.0000.dat");
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		volume[0] += myGrid3D(n)[1];
	MMSP::update(myGrid3D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		volume[1] += myGrid3D(n)[1];
	myGrid3D.output("short3.1000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<1; i++) {
		double myVol(volume[i]);
		MPI::COMM_WORLD.Allreduce(&myAre,&volume[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		for (int i=0; i<MMSP::fields(myGrid3D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid3D(n)[i]);
	EXPECT_LE(volume[1],volume[0]);
}

#endif
