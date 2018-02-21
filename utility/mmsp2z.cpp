// mmsp2z.cpp
// Compress MMSP grid data using zlib
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<sstream>
#include<zlib.h>
#include"MMSP.hpp"

bool fexists(const char *filename) {
	std::ifstream ifile(filename);
	return (bool)ifile;
}

int main(int argc, char* argv[]) {
#ifdef DEBUG
	int checkpoint = 0;
#endif
	// command line error check
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " [--help] infile outfile\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": compress MMSP grid data using zlib.\n";
		std::cout << "Usage: " << argv[0] << " [--help] infile outfile\n\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n\n";
		exit(0);
	}

	// input file error check
	std::istream* input;
	if (std::string(argv[1]) == "--") input = &std::cin;
	else {
		input = new std::ifstream(argv[1]);
		if (!(*input)) {
			std::cerr << "File (*input) error: could not open " << argv[1] << ".\n\n";
			exit(-1);
		}
	}

	// output file error check
	if ( fexists(argv[2]) ) {
		std::cerr << "File output warning: " << argv[2] << " already exists." << std::endl;
		//exit(-1);
	}
	std::ofstream output(argv[2]);
	if (!output) {
		std::cerr << "File output error: could not open ";
		std::cerr << argv[2] << "." << std::endl;
		exit(-1);
	}

	// read data type
	std::string type;
	getline((*input), type, '\n');
	output << type << '\n';

	// read grid dimension
	int dim;
	(*input) >> dim;
	if (dim < 1 or dim > 3) {
		std::cerr << "File input error: grid dimension must be 1, 2, or 3." << std::endl;
		exit(-1);
	}

	output << dim << '\n';

	// read number of fields
	int fields;
	(*input) >> fields;
	output << fields << '\n';

	// read grid sizes
	int x0[3] = {0, 0, 0};
	int x1[3] = {0, 0, 0};
	for (int i = 0; i < dim; i++) {
		(*input) >> x0[i] >> x1[i];
		output << x0[i] << ' ' << x1[i] << '\n';
	}

	// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int i = 0; i < dim; i++) {
		(*input) >> dx[i];
		output << dx[i] << '\n';
	}

	// ignore trailing endlines
	input->ignore(10, '\n');



	// read number of blocks
	int blocks;
	input->read(reinterpret_cast<char*>(&blocks), sizeof(blocks));
	output.write(reinterpret_cast<char*>(&blocks), sizeof(blocks));

	for (int i = 0; i < blocks; i++) {
#ifdef DEBUG
		checkpoint = 0;
		std::cout << "\nProcessing block " << i << ".\n" << std::endl;
#endif
		// read block limits
		int lmin[3] = {0, 0, 0};
		int lmax[3] = {0, 0, 0};
		for (int j = 0; j < dim; j++) {
			input->read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
			output.write(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
			input->read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
			output.write(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
		}
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif

		// read grid data
		unsigned long size, compressedSize;
		input->read(reinterpret_cast<char*>(&size), sizeof(size)); // read raw size
		output.write(reinterpret_cast<char*>(&size), sizeof(size)); // write raw size
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif
		compressedSize = 1.125 * size + 12;
		assert(compressedSize - 12 > size);
		char* raw_buffer = new char[size];
		input->read(raw_buffer, size);
		char* buffer = new char[compressedSize];
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif
		int status;
		const int level = 9;
		status = compress2(reinterpret_cast<unsigned char*>(buffer), &compressedSize, reinterpret_cast<unsigned char*>(raw_buffer), size, level);
		//status = compress2(reinterpret_cast<unsigned char*>(buffer), &compressedSize, reinterpret_cast<unsigned char*>(input), size, level); // this segfaults
		switch( status ) {
		case Z_OK:
			break;
		case Z_MEM_ERROR:
			std::cerr << "Compress: out of memory." << std::endl;
			exit(1);    // quit.
			break;
		case Z_BUF_ERROR:
			std::cerr << "Compress: output buffer wasn't large enough." << std::endl;
			exit(1);    // quit.
			break;
		}
		delete [] raw_buffer;
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif
		output.write(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize)); // write compressed size
		output.write(reinterpret_cast<char*>(buffer), compressedSize); // write compressed data
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif

		// clean up
		delete [] buffer;
#ifdef DEBUG
		++checkpoint;
		std::cout << "Debug:\tcheckpoint " << checkpoint << " reached." << std::endl;
#endif

	}

	if (input != &std::cin) delete input;
#ifdef DEBUG
	std::cout << "\nDebug:\tFinished conversion of " << argv[2] << "." << std::endl;
#endif
	return 0;
}

