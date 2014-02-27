// File:    mmsp2xyz.cpp
// Purpose: reads MMSP grid containing sparse floats, tracks specified grain
// Output:  CSV file specifying XYZ+phase
// Depends: MMSP, zlib

// Questions/Comments to kellet@rpi.edu (Trevor Keller)

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <zlib.h>

#include "MMSP.hpp"

int main(int argc, char* argv[]) {
	if ( argc != 4 ) {
		std::cout << "Usage: " << argv[0] << " data.dat grain_id output.csv\n";
		return ( 1 );
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
	}

	const int id(atoi(argv[2]));
	if (id == 0) {
		std::cerr << "Input warning: " << argv[2] << " does not appear to be an integer." << std::endl;
	}

	// read data type
	std::string type;
	getline(input, type, '\n');

	// grid type error check: read line, "grid:sparse:float"
	if (type.substr(0, 4) != "grid") {
		std::cerr << "File input error: file does not contain grid data." << std::endl;
		exit(-1);
	} else if (type.substr(5, 6) != "sparse") {
		std::cerr << "File input error: grid does not contain sparse data." << std::endl;
		exit(-1);
	} else if (type.substr(12, 5) != "float") {
		std::cerr << "File input error: vector data does not contain floats." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;

	std::vector<MMSP::vector<int> > points;
	std::vector<float> weights;

	if (dim == 2) {
		std::cout << "XYZ implies 3D data." << std::endl;
		exit(1);
	} else if (dim == 3) {
		// construct grid object
		MMSP::grid<3, MMSP::sparse<float> > grid(argv[1]);
		for (int d = 0; d < dim; ++d) dx(grid, d) = 0.375;

		// Populate the image from MMSP data.
		MMSP::vector<int> min(3, 500); // origin of grain
		for (int n = 0; n < nodes(grid); ++n) {
			MMSP::vector<int> x = position(grid, n);
			int S = length(grid(n));
			for (int s = 0; s < S; ++s) {
				if (MMSP::index(grid(n), s) == id) {
					points.push_back(x);
					weights.push_back(grid(n)[id]);
					for (int d = 0; d < 3; ++d) if (x[d] < min[d]) min[d] = x[d];
				}
			}
		}
		if (points.size() != weights.size()) {
			std::cout << "Error: XYZ vs weight size mismatch." << std::endl;
			exit(1);
		}
		std::ofstream output(argv[3]);
		if (!output) {
			std::cerr << "File output error: could not create " << argv[3] << ".\n\n";
			exit(-1);
		}
		output << "#X,Y,Z,nx,ny,nz\n";
		for (int i = 0; i < points.size(); ++i) {
			MMSP::vector<MMSP::sparse<float> > normal = MMSP::gradient(grid, points[i]);
			for (int d = 0; d < 3; ++d) output << dx(grid, d)*double(points[i][d] - min[d]) << '\t';
			for (int d = 0; d < 3; ++d) {
				output << dx(grid, d)*double(normal[d][id]);
				if (d < dim - 1) output << '\t';
				else output << '\n';
			}
		}
	} else {
		std::cerr << "Error: " << dim << "-D data is not supported!" << std::endl;
		exit(1);
	}

	return 0;
}
