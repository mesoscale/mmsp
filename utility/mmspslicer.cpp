// File:    mmspslicer.cpp
// Purpose: reads MMSP grid containing sparse floats, creates image from magnitude.
// Output:  grayscale portable network graphics
// Depends: MMSP, DevIL image library, zlib

// Questions/Comments to trevor.keller@gmail.com (Trevor Keller)

// DevIL usage after http://bobobobo.wordpress.com/2009/03/02/how-to-use-openil-to-generate-and-save-an-image/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <zlib.h>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

#include "MMSP.hpp"
#include "devil_cpp_wrapper.hpp"

using namespace std;

int main(int argc, char* argv[]) {
	if ( argc != 3 ) {
		cout << "Usage: " << argv[0] << " data.dat output.csv\n";
		return ( 1 );
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
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

	if (dim == 2) {
		// Initialize image
		ilInit();

		// construct grid object
		MMSP::grid<2, MMSP::sparse<float> > grid(argv[1]);
		int bytesPerPx = 1; // set to 1 for 8-bit grayscale
		int width = MMSP::g1(grid, 0) - MMSP::g0(grid, 0);
		int height = MMSP::g1(grid, 1) - MMSP::g0(grid, 1);
		unsigned int theSize = height * width * bytesPerPx; // number of char's
		unsigned char* imData = new unsigned char[theSize];

#ifdef DEBUG
		cout << "Grid contains " << MMSP::nodes(grid) << " nodes. Allocated array of " << theSize << " bytes." << endl;
#endif

		// Populate the image from MMSP data.
		MMSP::vector<int> x(2);
		int i = 0;
		for (x[1] = MMSP::y0(grid); x[1] < MMSP::y1(grid); ++x[1]) {
			for (x[0] = MMSP::x0(grid); x[0] < MMSP::x1(grid); ++x[0]) {
				double sum = 0.;
				int s = MMSP::length(grid(x));
				for (int h = 0; h < s; ++h) sum += static_cast<double>(grid(x).value(h) * grid(x).value(h));
				for (int j = 0; j < bytesPerPx; ++j) imData[ i + j ] = 255 * sqrt(sum);
				++i;
			}
		}

		ILenum Error;
		ILuint imageID = ilGenImage() ;
		ilBindImage(imageID);
		ilTexImage(width, height, 1, bytesPerPx, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData);
		Error = ilGetError();
		if (Error != IL_NO_ERROR) cout << "Error making image: " << iluErrorString(Error) << endl;
		ilEnable(IL_FILE_OVERWRITE);
		ilSave( IL_PNG, argv[2] ) ;
		Error = ilGetError();
		if (Error != IL_NO_ERROR) cout << "Error saving image: " << iluErrorString(Error) << endl;

		delete [] imData;
	} else if (dim == 3) {
		// Initialize image
		ilInit();

		// construct grid object
		MMSP::grid<3, MMSP::sparse<float> > grid(argv[1]);
		int bytesPerPx = 1; // set to 1 for 8-bit grayscale
		int width = MMSP::g1(grid, 0) - MMSP::g0(grid, 0);
		int height = MMSP::g1(grid, 1) - MMSP::g0(grid, 1);
		unsigned int theSize = height * width * bytesPerPx; // number of char's
		unsigned char* imData[5];
		for (int i = 0; i < 5; ++i) imData[i] = new unsigned char[theSize];

#ifdef DEBUG
		cout << "Grid contains " << MMSP::nodes(grid) << " nodes. Allocated array of " << theSize << " bytes." << endl;
#endif

		// Populate the image from MMSP data.
		MMSP::vector<int> x(3);
		for (int f = 0; f < 5; ++f) {
			x[2] = MMSP::z0(grid) + int((0.25 + 0.125 * f) * (MMSP::z1(grid) - MMSP::z0(grid)));
			int i = 0;
			for (x[1] = MMSP::y0(grid); x[1] < MMSP::y1(grid); ++x[1]) {
				for (x[0] = MMSP::x0(grid); x[0] < MMSP::x1(grid); ++x[0]) {
					double sum = 0.;
					int s = MMSP::length(grid(x));
					for (int h = 0; h < s; ++h) sum += static_cast<double>(grid(x).value(h) * grid(x).value(h));
					for (int j = 0; j < bytesPerPx; ++j) imData[f][ i + j ] = 255 * sqrt(sum); //(sqrt(sum)>0.45 && sqrt(sum)<0.55)?0:255;
					++i;
				}
			}
			string filename(argv[2]);
			filename.insert(filename.find_last_of('.') + 1, 1, char(49 + f));
			filename.insert(filename.find_last_of('.') + 2, 1, '.');
			ILenum Error;
			ILuint imageID = ilGenImage() ;
			ilBindImage(imageID);
			ilTexImage(width, height, 1, bytesPerPx, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData[f]);
			Error = ilGetError();
			if (Error != IL_NO_ERROR) cout << "Error making image: " << iluErrorString(Error) << endl;
			ilEnable(IL_FILE_OVERWRITE);
			ilSave( IL_PNG, filename.c_str() ) ;
			Error = ilGetError();
			if (Error != IL_NO_ERROR) cout << "Error saving image: " << iluErrorString(Error) << endl;
		}

		for (int i = 0; i < 5; ++i) delete [] imData[i];
	} else {
		cerr << "Error: " << dim << "-D data is not supported!" << endl;
		exit(1);
	}

	return 0;
}
