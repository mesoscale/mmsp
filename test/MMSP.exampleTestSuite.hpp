/* MMSP.unitTest.hpp
 * Use Google Test to build unit tests for MMSP classes
 * http://github.com/google/googletest
 */

#include"MMSP.hpp"
#include<cmath>
#include<gtest/gtest.h>

// GTest supports templated unit testing. Declare the types of interest.

//using testing::Types;
//typedef Types<char,short,int,long,float,double> datatypes;

TEST(grid1DTest, Finite) {
	MMSP::vector<double> length(3,0.0);
	MMSP::generate(1, "short1.0000.dat");
	GRID1D myGrid1D("short1.0000.dat");
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		length[0] += myGrid1D(n)[1];
	MMSP::update(myGrid1D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		length[1] += myGrid1D(n)[1];
	myGrid1D.output("short1.1000.dat");
	MMSP::update(myGrid1D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		length[2] += myGrid1D(n)[1];
	myGrid1D.output("short1.2000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<3; i++) {
		double myLen(length[i]);
		MPI::COMM_WORLD.Allreduce(&myLen,&length[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid1D); n++)
		for (int i=0; i<MMSP::fields(myGrid1D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid1D(n)[i]);
	for (int i=1; i<3; i++)
		EXPECT_GE(length[i-1],length[i]);
}

TEST(grid2DTest, Finite) {
	MMSP::vector<double> area(3,0.0);
	MMSP::generate(2, "short2.0000.dat");
	GRID2D myGrid2D("short2.0000.dat");
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		area[0] += myGrid2D(n)[1];
	MMSP::update(myGrid2D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		area[1] += myGrid2D(n)[1];
	myGrid2D.output("short2.1000.dat");
	MMSP::update(myGrid2D, 2000);
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		area[2] += myGrid2D(n)[1];
	myGrid2D.output("short2.2000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<3; i++) {
		double myAre(area[i]);
		MPI::COMM_WORLD.Allreduce(&myAre,&area[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid2D); n++)
		for (int i=0; i<MMSP::fields(myGrid2D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid2D(n)[i]);
	for (int i=1; i<3; i++)
		EXPECT_GE(area[i-1],area[i]);
}

TEST(grid3DTest, Finite) {
	MMSP::vector<double> volume(3,0.0);
	MMSP::generate(3, "short3.0000.dat");
	GRID3D myGrid3D("short3.0000.dat");
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		volume[0] += myGrid3D(n)[1];
	MMSP::update(myGrid3D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		volume[1] += myGrid3D(n)[1];
	myGrid3D.output("short3.1000.dat");
	MMSP::update(myGrid3D, 1000);
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		volume[2] += myGrid3D(n)[1];
	myGrid3D.output("short3.2000.dat");
	#ifdef MPI_VERSION
	for (int i=0; i<3; i++) {
		double myVol(volume[i]);
		MPI::COMM_WORLD.Allreduce(&myAre,&volume[i],1,MPI_DOUBLE,MPI_SUM);
	}
	#endif
	for (int n=0; n<MMSP::nodes(myGrid3D); n++)
		for (int i=0; i<MMSP::fields(myGrid3D); i++)
			EXPECT_PRED1(std::isfinite<double>,myGrid3D(n)[i]);
	for (int i=1; i<3; i++)
		EXPECT_GE(volume[i-1],volume[i]);
}
