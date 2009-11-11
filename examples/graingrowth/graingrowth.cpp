// graingrowth.cpp 
// Example program for isotropic grain growth using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MCgrid.isotropy.hpp"
#include"PFgrid.isotropy.hpp"
#include"sparsePF.isotropy.hpp"

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
	bool MC = (type=="grid:int") or (type=="grid:unsigned int") or
              (type=="grid:scalar:int") or (type=="grid:scalar:unsigned int") or
              (type=="grid:scalar:long") or (type=="grid:scalar:unsigned long") or
              (type=="grid:scalar:short") or (type=="grid:scalar:unsigned short");

	bool PF = (type=="grid:vector:double");

	bool sparsePF = (type=="grid:sparse:double");


	// grain growth (MCgrid)
	if (MC) {
		if (dim==2) {
			// input the grid object
			MCgrid2D grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}

		if (dim==3) {
			// input the grid object
			MCgrid3D grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}
	}

	// grain growth (PFgrid)
	else if (PF) {
		if (dim==2) {
			// input the grid object
			PFgrid2D grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}

		if (dim==3) {
			// input the grid object
			PFgrid3D grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}
	}

	// grain growth (sparsePF)
	else if (sparsePF) {
		if (dim==2) {
			// input the grid object
			sparsePF2D grid(argv[1]);

			// grain growth simulation
			update(grid,atoi(argv[3]));

			// output the grid to a file
			output(grid,argv[2]);
		}

		if (dim==3) {
			// input the grid object
			sparsePF3D grid(argv[1]);

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
