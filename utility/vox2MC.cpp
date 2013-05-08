// vox2MC.cpp
// Convert MBuilder voxel data to MMSP grid data
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MMSP.hpp"

int main(int argc, char* argv[]) {
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] infile [outfile]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MBuilder voxel data to MMSP grid format.\n";
		std::cout << "Usage: " << argv[0] << " [--help] infile [outfile]\n\n";
		std::cout << "Questions/comments to gruberja@gmail.com (Jason Gruber).\n\n";
		exit(0);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
	}

	// generate output file name
	std::stringstream filename;
	if (argc < 3)
		filename << std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of(".")) << ".MC";
	else
		filename << argv[2];

	// get grid bounds
	int nx, ny, nz;
	input >> nx >> ny >> nz;

	// second and third line are not used
	std::string line;
	getline(input, line);
	getline(input, line);
	getline(input, line);

	// generate MMSP grid
	MMSP::grid<3, int> GRID(1, 0, nx, 0, ny, 0, nz);

	// read in MBuilder grid data
	for (int i = 0; i < MMSP::nodes(GRID); i++) input >> GRID(i);

	// output MMSP grid data
	MMSP::output(GRID, filename.str().c_str());
}

