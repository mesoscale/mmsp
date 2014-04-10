// MC2vox.cpp
// Convert MMSP grid data to MBuilder voxel data
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MMSP.hpp"
#include<set>

int main(int argc, char* argv[]) {
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] infile [outfile]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP MCgrid data format to PFgrid data format.\n";
		std::cout << "Usage: " << argv[0] << " [--help] infile [outfile]\n";
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

	bool scalar_type = (type.find("scalar") != std::string::npos);
	bool vector_type = (type.find("vector") != std::string::npos);
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
	if (not int_type or vector_type or sparse_type) {
		std::cerr << "File input error: data must be of type int or scalar::int." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;
	if (not dim == 3) {
		std::cerr << "File input error: grid dimension must be 3." << std::endl;
		exit(-1);
	}

	// open output file
	std::ofstream output(filename.str().c_str());
	if (!output) {
		std::cerr << "File output error: could not open " << filename.str() << ".\n\n";
		exit(-1);
	}

	// read in MMSP grid data
	MMSP::grid<3, int> GRID(argv[1]);

	// compute number of grains
	std::set<int> gset;
	for (int i = 0; i < MMSP::nodes(GRID); i++) gset.insert(GRID(i));
	int grains = gset.size();

	// output grid dimensions
	output << MMSP::xlength(GRID) << " ";
	output << MMSP::ylength(GRID) << " ";
	output << MMSP::zlength(GRID) << "\n";
	output << "\'grwXXX\'				52.00 1.000 1.0		" << grains << "\n";
	output << "3.000 0.000 0.000			0\n";

	// output grid data
	int count = 0;
	for (int i = 0; i < MMSP::nodes(GRID); i++) {
		output << GRID(i) << " ";
		count += 1;
		if (count == 20) {
			output << "\n";
			count = 0;
		}
	}
}

