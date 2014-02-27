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
	int id = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	std::cout<<"This is a test (processor "<<id<<" of "<<np<<")."<<std::endl;
	#endif

	MMSP::Finalize();
}
