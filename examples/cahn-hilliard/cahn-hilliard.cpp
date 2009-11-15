// cahn-hilliard.cpp 
// Example program for the Cahn-Hilliard model using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"cahn-hilliard.hpp"

using namespace MMSP;

int main(int argc, char* argv[])
{
	Init(argc,argv);

	// command line error check
	if (argc<3) {
		std::cout<<"Usage: "<<argv[0]<<" inputfile outputfile timesteps\n";
		exit(-1);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open ";
		std::cerr<<argv[1]<<"."<<std::endl;
		exit(-1);
	}

	// read grid data type
	std::string type;
	getline(input,type,'\n');

	// grid type error check
	if (type.substr(0,4)!="grid") {
		std::cerr<<"File input error: file does not contain grid data."<<std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input>>dim;

	input.close();

	// set grain growth type
	bool cahn_hilliard = (type=="grid:double") or (type=="grid:scalar:double");

	// Cahn-Hilliard model
	if (cahn_hilliard) {
		if (dim==2) {
			// input the grid object
			MMSP::grid<2,double> grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}

		if (dim==3) {
			// input the grid object
			MMSP::grid<3,double> grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}
	}

	// unsupported grid type error
	else {
		std::cerr<<"File input error: unsupported grid data type."<<std::endl;
		exit(-1);
	}

	Finalize();
}
