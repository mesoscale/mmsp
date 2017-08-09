// sPF2PF.cpp
// Convert sparsePF data to phase field data
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
		std::cout << argv[0] << ": convert MMSP sparsePF data format to PFgrid data format.\n";
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
		filename << std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of(".")) << ".PF";
	else
		filename << argv[2];

	// read data type
	std::string type;
	getline(input, type, '\n');

	// grid type error check
	if (type.substr(0, 4) != "grid") {
		std::cerr << "File input error: file does not contain grid data." << std::endl;
		exit(-1);
	}

	// parse data type
	bool bool_type = (type.find("bool") != std::string::npos);
	bool char_type = (type.find("char") != std::string::npos);
	bool unsigned_char_type = (type.find("unsigned char") != std::string::npos);
	bool int_type = (type.find("int") != std::string::npos);
	bool unsigned_int_type = (type.find("unsigned int") != std::string::npos);
	bool long_type = (type.find("long") != std::string::npos);
	bool unsigned_long_type = (type.find("unsigned long") != std::string::npos);
	bool short_type = (type.find("short") != std::string::npos);
	bool unsigned_short_type = (type.find("unsigned short") != std::string::npos);
	bool float_type = (type.find("float") != std::string::npos);
	bool double_type = (type.find("double") != std::string::npos);
	bool long_double_type = (type.find("long double") != std::string::npos);

	bool sparse_type = (type.find("sparse") != std::string::npos);

	if (not bool_type    and
	    not char_type    and  not unsigned_char_type   and
	    not int_type     and  not unsigned_int_type    and
	    not long_type    and  not unsigned_long_type   and
	    not short_type   and  not unsigned_short_type  and
	    not float_type   and
	    not double_type  and  not long_double_type) {
		std::cerr << "File input error: unknown grid data type." << std::endl;
		exit(-1);
	}

	// check for valid MCgrid data
	if (not sparse_type or not double_type) {
		std::cerr << "File input error: data must be of type sparse::double." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;
	if (dim < 1 or dim > 3) {
		std::cerr << "File input error: grid dimension must be 1, 2, or 3." << std::endl;
		exit(-1);
	}

	if (dim == 1) {
		MMSP::grid<1, MMSP::sparse<double> > grid1(argv[1]);
		int x0 = MMSP::x0(grid1);
		int x1 = MMSP::x1(grid1);

		int fields = 0;
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				if (index > fields) fields = index;
			}
		}
		fields += 1;

		MMSP::grid<1, MMSP::vector<double> > grid2(fields, x0, x1);
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			for (int j = 0; j < fields; j++)
				grid2(i)[j] = 0.0;
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				double value = MMSP::value(grid1(i), j);
				grid2(i)[index] = value;
			}
		}
		MMSP::output(grid2, filename.str().c_str());
	}

	if (dim == 2) {
		MMSP::grid<2, MMSP::sparse<double> > grid1(argv[1]);
		int x0 = MMSP::x0(grid1);
		int x1 = MMSP::x1(grid1);
		int y0 = MMSP::y0(grid1);
		int y1 = MMSP::y1(grid1);

		int fields = 0;
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				if (index > fields) fields = index;
			}
		}
		fields += 1;

		MMSP::grid<2, MMSP::vector<double> > grid2(fields, x0, x1, y0, y1);
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			for (int j = 0; j < fields; j++)
				grid2(i)[j] = 0.0;
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				double value = MMSP::value(grid1(i), j);
				grid2(i)[index] = value;
			}
		}
		MMSP::output(grid2, filename.str().c_str());
	}

	if (dim == 3) {
		MMSP::grid<3, MMSP::sparse<double> > grid1(argv[1]);
		int x0 = MMSP::x0(grid1);
		int x1 = MMSP::x1(grid1);
		int y0 = MMSP::y0(grid1);
		int y1 = MMSP::y1(grid1);
		int z0 = MMSP::z0(grid1);
		int z1 = MMSP::z1(grid1);

		int fields = 0;
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				if (index > fields) fields = index;
			}
		}
		fields += 1;

		MMSP::grid<2, MMSP::vector<double> > grid2(fields, x0, x1, y0, y1, z0, z1);
		for (int i = 0; i < MMSP::nodes(grid1); i++) {
			for (int j = 0; j < fields; j++)
				grid2(i)[j] = 0.0;
			int size = MMSP::length(grid1(i));
			for (int j = 0; j < size; j++) {
				int index = MMSP::index(grid1(i), j);
				double value = MMSP::value(grid1(i), j);
				grid2(i)[index] = value;
			}
		}
		MMSP::output(grid2, filename.str().c_str());
	}
}

