// mmsp2png.cpp
// Convert MMSP grid data to grayscale PNG image format
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include <iostream>
#include <cstdlib>
#include<zlib.h>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<cmath>
#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>
#include "devil_cpp_wrapper.hpp"

#include"MMSP.hpp"

int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename);

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,MMSP::scalar<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const unsigned int& bufsize, unsigned char* buffer);

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer);

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer);

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] [--mag|--field|--exclude|--contour] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to VTK image data format.\n";
		std::cout << "Usage:    " << argv[0] << " [--help] [--slice=X,n] [--mag|--field=n|--exclude=n|--contour=n,a] [--invert] infile [outfile]\n";
		std::cout << "Examples: " << argv[0] << " --help\n"
		          << "             displays this message.\n";
		std::cout << "          " << argv[0] << " infile\n"
		          << "             writes magnitude of infile values to infile.png, rescaled on [0,255].\n"
		          << "             If infile is 3-D, holds X-coord constant at the midpoint.\n";
		std::cout << "          " << argv[0] << " --slice=X,32 infile outfile\n"
		          << "             writes magnitude of infile values with constant X-coord of 32 to outfile, rescaled on [0,255].\n"
		          << "             Use for 3-D data in conjunction with any of the other flags. Y and Z also work. Defaults to the mid-plane parallel to the specified axis.\n";
		std::cout << "          " << argv[0] << " infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --mag infile.dat \n"
		          << "             writes magnitude of infile.dat values to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --field=0 infile.dat\n"
		          << "             writes field 0 of infile.dat values, only, to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --field=1,2,3 infile.dat\n"
		          << "             writes magnitude of fields 1,2,3 of infile.dat values, only, to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --exclude=1,2,3 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except fields 1,2,3, to infile.png, rescaled on [0,255].\n";
		std::cout << "          " << argv[0] << " --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, with inversions where field 0 equals 0.25,0.50,0.75 to infile.png, rescaled on [0,255].\n"
		          << "             Works in conjunction with --exclude as well.\n";
		std::cout << "          " << argv[0] << " --invert infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, rescaled on [255,0].\n"
		          << "             Use in conjunction with any other flag: inversion is the last operation before output.\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n";
		std::exit(0);
	}

	int datindex = 0; // in typical usage, input filename comes immediately after executable
	int pngindex = 0; // in typical usage, output filename comes immediately after input
	int mode = 0;
	bool invert=false;
	int sliceaxis = -1;
	int slicelevel = -1;
	std::set<int> fieldset; // for --field= or exclude=
	std::set<double> levelset; // for contours

	// Parse command line flags
	for (int i=1; i<argc; i++) {
		std::string flag(argv[i]);
		if (flag == "--mag") {
			mode=1;
		} else if (flag == "--invert") {
			invert=true;
		} else if (flag.substr(0,7) == "--slice") {
			if (flag[7]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=X[|Y|Z][,n].\n"<<std::endl;
				std::exit(-1);
			}
			switch(flag[8]) {
				case 'x': sliceaxis=0;
				          break;
				case 'y': sliceaxis=1;
				          break;
				case 'z': sliceaxis=2;
				          break;
				case 'X': sliceaxis=0;
				          break;
				case 'Y': sliceaxis=1;
				          break;
				case 'Z': sliceaxis=2;
				          break;
			}
			if (sliceaxis==-1) {
				std::cerr<<"Error: axis "<<flag[8]<<" not recognized.\n"<<std::endl;
				std::exit(-1);
			}
			if (flag.length()>10)
				slicelevel = atoi(flag.substr(10,20).c_str());
			else
				std::cout<<"Slicing through the midpoint of "<<flag[8]<<" axis."<<std::endl;
		} else if (flag.substr(0,7) == "--field") {
			if (flag[7]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(2,mode);
			std::string val = flag.substr(8,255);
			char * cstr = new char [val.length()+1];
			std::strcpy (cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(cstr,",");
			}
			delete [] cstr;
		} else if (flag.substr(0,9) == "--exclude") {
			if (flag[9]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(3,mode);
			std::string val = flag.substr(10,255);
			char * cstr = new char [val.length()+1];
			std::strcpy (cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(cstr,",");
			}
			delete [] cstr;
		} else if (flag.substr(0,9) == "--contour") {
			if (flag[9]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n,a[,b,c,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(4,mode);
			std::string val = flag.substr(10,255);
			char * cstr = new char [val.length()+1];
			std::strcpy (cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			fieldset.insert(atoi(ptr));
			ptr = strtok(cstr,",");
			while (ptr != NULL) {
				levelset.insert(atof(ptr));
				ptr = strtok(cstr,",");
			}
			delete [] cstr;
		} else if (flag.find(".png") != std::string::npos) {
			pngindex = i;
		} else {
			datindex = i;
		}
	}

// file open error check
	if (datindex==0) {
		std::cerr << "File input error: no input specified."<<".\n";
		exit(-1);
	}
	std::ifstream input(argv[datindex]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[datindex] << ".\n";
		exit(-1);
	}

// generate output file name
	std::stringstream pngname;
	if (pngindex==0) {
		std::string datname(argv[datindex]);
		int extpos = datname.find_last_of(".");
		if (datname.find_first_of("0123456789",extpos) != std::string::npos)
			pngname << datname << ".png";
		else
			pngname << datname.substr(0, extpos) << ".png";
	} else
		pngname << argv[pngindex];

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

// read grid dimension
	int dim;
	input >> dim;

// read number of fields
	int fields;
	input >> fields;

// read grid sizes
	int x0[3] = {0, 0, 0};
	int x1[3] = {0, 0, 0};
	for (int d=0; d<dim; d++)
		input >> x0[d] >> x1[d];

// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int d=0; d<dim; d++)
		input >> dx[d];

// ignore trailing endlines
	input.ignore(10, '\n');

	int image_size[2] = {1,1};
	int i=0;
	for (int d=0; d<dim; d++) {
		if (d==sliceaxis)
			continue;
		image_size[i] = x1[d]-x0[d];
		i++;
	}
  unsigned int theSize=image_size[0] * image_size[1];
  unsigned char* buffer = new unsigned char[theSize];

		// write grid data
		if (scalar_type or (not vector_type and not sparse_type)) { // must be scalar or built-in
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<bool> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<bool> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<bool> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<float> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<float> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<float> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<long double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<long double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<long double> > GRID(argv[datindex]);
					convert_scalars(GRID, mode, sliceaxis, slicelevel, levelset, theSize, buffer);
				}
			}
		}
		else if (vector_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, sliceaxis, slicelevel, levelset, fieldset, theSize, buffer);
				}
			}
		}

	if (invert)
		for (unsigned int n=0; n<theSize; n++)
			buffer[n] = 255-buffer[n];

	int result = writePNG(image_size[0], image_size[1], 1, buffer, pngname.str().c_str());

	// clean up
	delete [] buffer;

	return result;
}

int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename)
{
  int result = 0;
  // Initialize image
  ilInit();
  ILenum Error;
  ILuint imageID = ilGenImage() ;
  ilBindImage(imageID);
  ilTexImage(w, h, 1, bpp, IL_LUMINANCE, IL_UNSIGNED_BYTE, imData);
  Error = ilGetError();
  if (Error!=IL_NO_ERROR) {
  	std::cout<<"Error making image: "<<iluErrorString(Error)<<std::endl;
  	result = -1;
  }
  ilEnable(IL_FILE_OVERWRITE);
  ilSave( IL_PNG, filename) ;
  Error = ilGetError();
  if (Error!=IL_NO_ERROR) {
  	std::cout<<"Error saving image: "<<iluErrorString(Error)<<std::endl;
  	result = -1;
  }
  return result;
}

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,MMSP::scalar<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const unsigned int& bufsize, unsigned char* buffer)
{
	T min=0;
	T max=1;
	for (int n=0; n<MMSP::nodes(GRID); n++) {
		if (GRID(n)>max)
			max=GRID(n);
		else if (GRID(n)<min)
			min=GRID(n);
	}
	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1, 0);
		for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
			T val = GRID(x);
			assert(n<bufsize);
			buffer[n] = 255*((val-min)/(max-min));
			if (mode==4) //contour
				for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<1.0e-2)
						buffer[n] = 255-buffer[n];
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2, 0);
		for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
			for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
				T val = GRID(x);
				assert(n<bufsize);
				buffer[n] = 255*((val-min)/(max-min));
				if (mode==4) //contour
					for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
						if (std::fabs(val-*it)/std::fabs(*it)<1.0e-2)
							buffer[n] = 255-buffer[n];
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3, 0);
		for (x[2] = MMSP::x0(GRID,2); x[2] < MMSP::x1(GRID,2); x[2]++)
			for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
				for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					T val = GRID(x);
					assert(n<bufsize);
					buffer[n] = 255*((val-min)/(max-min));
					if (mode==4) //contour
						for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
							if (std::fabs(val-*it)/std::fabs(*it)<1.0e-2)
								buffer[n] = 255-buffer[n];
					n++;
				}
	}
}

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer)
{
	T min=0;
	T max=1;
	int included = MMSP::fields(GRID) - fieldset.size();

	for (int n=0; n<MMSP::nodes(GRID); n++) {
		double sum=0.0;
		if (mode<2) { //          --mag
			for (int i=0; i<MMSP::fields(GRID); i++)
				sum += pow(GRID(n)[i],2.0);
		} else if (mode==2) { //  --field
			if (fieldset.size()==1)
				sum = GRID(n)[*fieldset.begin()];
			else
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
					sum += pow(GRID(n)[*it],2.0);
		} else if (mode==3) { //  --exclude
			for (int i=0; i<MMSP::fields(GRID); i++) {
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end()) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			}
		} else if (mode==4) { //  --contour
			// Same as --exclude
			for (int i=0; i<MMSP::fields(GRID); i++) {
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end()) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			}
		}
		if (mode < 2 || fieldset.size()>1 || included>1)
			sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1, 0);
		for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int i=0; i<MMSP::fields(GRID); i++)
					sum += pow(GRID(n)[i],2.0);
			} else if (mode==2) { //  --field
				if (fieldset.size()==1)
					sum = GRID(n)[*fieldset.begin()];
				else
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
						sum += pow(GRID(n)[*it],2.0);
			} else if (mode==3) { //  --exclude
				for (int i=0; i<MMSP::fields(GRID); i++) {
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end()) {
						if (included>1)
							sum += pow(GRID(n)[i],2.0);
						else
							sum = GRID(n)[i];
					}
				}
			} else if (mode==4) { //  --contour
				// Same as --exclude
				for (int i=0; i<MMSP::fields(GRID); i++) {
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end()) {
						if (included>1)
							sum += pow(GRID(n)[i],2.0);
						else
							sum = GRID(n)[i];
					}
				}
			}
			if (mode < 2 || fieldset.size()>1 || included>1)
				sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (mode==4) // --contour
				for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
							buffer[n] = 255-buffer[n];
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2, 0);
		for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
			for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode<2) { //          --mag
					for (int i=0; i<MMSP::fields(GRID); i++)
						sum += pow(GRID(n)[i],2.0);
				} else if (mode==2) { //  --field
					if (fieldset.size()==1) {
						sum = GRID(n)[*fieldset.begin()];
					} else {
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
							sum += pow(GRID(n)[*it],2.0);
					}
				} else if (mode==3) { //  --exclude
					for (int i=0; i<MMSP::fields(GRID); i++) {
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end()) {
							if (included>1)
								sum += pow(GRID(n)[i],2.0);
							else
								sum = GRID(n)[i];
						}
					}
				} else if (mode==4) { //  --contour
					// Same as --exclude
					for (int i=0; i<MMSP::fields(GRID); i++) {
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end()) {
							if (included>1)
								sum += pow(GRID(n)[i],2.0);
							else
								sum = GRID(n)[i];
						}
					}
				}
				if (mode < 2 || fieldset.size()>1 || included>1)
					sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (mode==4) // --contour
					for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
								buffer[n] = 255-buffer[n];
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3, 0);
		for (x[2] = MMSP::x0(GRID,2); x[2] < MMSP::x1(GRID,2); x[2]++)
			for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
				for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode<2) { //          --mag
						for (int i=0; i<MMSP::fields(GRID); i++)
							sum += pow(GRID(n)[i],2.0);
					} else if (mode==2) { //  --field
						if (fieldset.size()==1)
							sum = GRID(n)[*fieldset.begin()];
						else
							for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
								sum += pow(GRID(n)[*it],2.0);
					} else if (mode==3) { //  --exclude
						for (int i=0; i<MMSP::fields(GRID); i++) {
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end()) {
								if (included>1)
									sum += pow(GRID(n)[i],2.0);
								else
									sum = GRID(n)[i];
							}
						}
					} else if (mode==4) { //  --contour
						// Same as --exclude
						for (int i=0; i<MMSP::fields(GRID); i++) {
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end()) {
								if (included>1)
									sum += pow(GRID(n)[i],2.0);
								else
									sum = GRID(n)[i];
							}
						}
					}
					if (mode < 2 || fieldset.size()>1 || included>1)
						sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (mode==4) // --contour
						for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
							for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
								if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
									buffer[n] = 255-buffer[n];
					n++;
				}
	}
}

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                    const int& mode, const int& sliceaxis, const int& slicelevel,
                                                    const std::set<double>& levelset, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer)
{
	T min=0;
	T max=1;

	for (int n=0; n<MMSP::nodes(GRID); n++) {
		double sum=0.0;
		if (mode<2) { //          --mag
			for (int h=0; h<GRID(n).length(); h++)
				sum += pow(GRID(n).value(h),2.0);
		} else if (mode==2) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
					sum += pow(GRID(n)[*it],2.0);
		} else if (mode==3) { //  --exclude
			for (int h=0; h<MMSP::fields(GRID); h++) {
				int i=GRID(n).index(i);
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					sum += pow(GRID(n).value(h),2.0);
			}
		} else if (mode==4) { //  --contour
			// Same as --exclude
			for (int h=0; h<MMSP::fields(GRID); h++) {
				int i=GRID(n).index(i);
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					sum += pow(GRID(n).value(h),2.0);
			}
		}
		sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1, 0);
		for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int h=0; h<GRID(n).length(); h++)
					sum += pow(GRID(n).value(h),2.0);
			} else if (mode==2) { //  --field
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
						sum += pow(GRID(n)[*it],2.0);
			} else if (mode==3) { //  --exclude
				for (int h=0; h<MMSP::fields(GRID); h++) {
					int i=GRID(n).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(n).value(h),2.0);
				}
			} else if (mode==4) { //  --contour
				// Same as --exclude
				for (int h=0; h<MMSP::fields(GRID); h++) {
					int i=GRID(n).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(n).value(h),2.0);
				}
			}
			sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (mode==4) // --contour
				for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
							buffer[n] = 255-buffer[n];
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2, 0);
		for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
			for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode<2) { //          --mag
					for (int h=0; h<GRID(n).length(); h++)
						sum += pow(GRID(n).value(h),2.0);
				} else if (mode==2) { //  --field
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
							sum += pow(GRID(n)[*it],2.0);
				} else if (mode==3) { //  --exclude
					for (int h=0; h<MMSP::fields(GRID); h++) {
						int i=GRID(n).index(i);
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end())
							sum += pow(GRID(n).value(h),2.0);
					}
				} else if (mode==4) { //  --contour
					// Same as --exclude
					for (int h=0; h<MMSP::fields(GRID); h++) {
						int i=GRID(n).index(i);
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end())
							sum += pow(GRID(n).value(h),2.0);
					}
				}
				sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (mode==4) // --contour
					for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
								buffer[n] = 255-buffer[n];
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3, 0);
		for (x[2] = MMSP::x0(GRID,2); x[2] < MMSP::x1(GRID,2); x[2]++)
			for (x[1] = MMSP::x0(GRID,1); x[1] < MMSP::x1(GRID,1); x[1]++)
				for (x[0] = MMSP::x0(GRID,0); x[0] < MMSP::x1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode<2) { //          --mag
						for (int h=0; h<GRID(n).length(); h++)
							sum += pow(GRID(n).value(h),2.0);
					} else if (mode==2) { //  --field
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
							sum += pow(GRID(n)[*it],2.0);
					} else if (mode==3) { //  --exclude
						for (int h=0; h<MMSP::fields(GRID); h++) {
							int i=GRID(n).index(i);
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end())
								sum += pow(GRID(n).value(h),2.0);
						}
					} else if (mode==4) { //  --contour
						// Same as --exclude
						for (int h=0; h<MMSP::fields(GRID); h++) {
							int i=GRID(n).index(i);
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end())
								sum += pow(GRID(n).value(h),2.0);
						}
					}
					sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (mode==4) // --contour
						for (std::set<int>::iterator itf=fieldset.begin(); itf!=fieldset.end(); itf++)
							for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
								if (std::fabs(GRID(x)[*itf]-*itl)/std::fabs(*itl)<1.0e-2)
									buffer[n] = 255-buffer[n];
					n++;
				}
	}
}
