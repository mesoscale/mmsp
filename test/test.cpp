// test.cpp
// Simple test program using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MMSP.hpp"

int main(int argc, char* argv[])
{
	MMSP::Init(argc,argv);

	#ifndef MPI_VERSION
	std::cout<<"This is a test."<<std::endl;
	#endif

	#ifdef MPI_VERSION
	int id = 0, np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
	std::cout<<"This is a test (processor "<<id<<" of "<<np<<")."<<std::endl;
	#endif

	MMSP::Finalize();
}
