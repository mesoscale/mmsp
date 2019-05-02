// mmsp2png.cpp
// Convert MMSP grid data to grayscale PNG image format
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<iostream>
#include<cstdlib>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<cmath>
#include <png.h>

#include"MMSP.hpp"

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] [--slice] [--zoom] [--contour] [--mag|--field|--exclude] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to PNG image format.\n";
		std::cout << "Usage:    " << argv[0] << " [--help] [--slice=X,n] [--zoom] [--contour=n,a,b] [--mag|--field=n|--exclude=n] [--invert] infile [outfile]\n\n";
		std::cout << "Examples: " << argv[0] << " --help\n"
		          << "             displays this message.\n";
		std::cout << "          " << argv[0] << " infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded on [0,1].\n"
		          << "             If infile is 3-D, holds X-coord constant at the midpoint.\n";
		std::cout << "          " << argv[0] << " --zoom infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded by data min and max exactly.\n"
		          << "             Use in conjunction with any other flag for data with very small values.\n";
		std::cout << "          " << argv[0] << " --zoom=a,b infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded on [a,b].\n"
		          << "             Use in conjunction with any other flag for data with very small values.\n";
		std::cout << "          " << argv[0] << " --slice=X,32 infile outfile\n"
		          << "             writes magnitude of infile values with constant X-coord of 32 to outfile, bounded on [0,1].\n"
		          << "             Use for 3-D data in conjunction with any other flag. Y and Z also work.\n"
		          << "             Defaults to the mid-plane parallel to the specified axis.\n";
		std::cout << "          " << argv[0] << " infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --mag infile.dat \n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [0,1].\n"
		          << "             Use with only one field (--field=n or --exclude all others) to visualize its absolute value.\n";
		std::cout << "          " << argv[0] << " --field=0 infile.dat\n"
		          << "             writes field 0 of infile.dat values, only, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --field=1,2,3 infile.dat\n"
		          << "             writes magnitude of fields 1,2,3 of infile.dat values, only, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=1,2,3 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except fields 1,2,3,to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             makes black pixels where field 0 equals 0.25,0.50,0.75. Use with any other flag.\n";
		std::cout << "          " << argv[0] << " --coninv --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             makes white pixels where field 0 equals 0.25,0.50,0.75. Use with any other flag.\n";
		std::cout << "          " << argv[0] << " --invert infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [255,0].\n"
		          << "             Use in conjunction with any other flag: inversion is the last operation before output.\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n";
		std::exit(0);
	}

	int datindex = 0; // in typical usage, input filename comes immediately after executable
	int pngindex = 0; // in typical usage, output filename comes immediately after input
	int mode = 0;
	bool invert=false; // filter to invert the whole image before writing
	bool coninv=false; // filter to set contour pixels white instead of black
	double contol = 0.06; // contour tolerance (5-7% works OK)
	int sliceaxis = -1;
	int slicelevel = -1;
	double zoomin=0.0, zoomax=1.0;
	int lvlfield = -1; // store the field on which contour are calculated
	std::set<int> fieldset; // for --field or --exclude
	std::set<double> levelset; // for contours

	// Parse command line flags
	for (int i=1; i<argc; i++) {
		std::string flag(argv[i]);
		if (flag == "--mag") {
			mode=std::max(1,mode);
		} else if (flag == "--invert") {
			invert=true;
		} else if (flag.substr(0,6) == "--zoom") {
			if (flag.length()==6 || flag[6]!='=') {
				zoomin = std::sqrt(std::numeric_limits<double>::max());
				zoomax = -zoomin;
			} else {
				std::string val = flag.substr(7,255);
				char* cstr = new char [val.length()+1];
				std::strcpy(cstr, val.c_str());
				char* ptr = strtok(cstr,",");
				zoomin = atof(ptr);
				ptr = strtok(NULL,",");
				zoomax = atof(ptr);
				ptr = strtok(NULL,",");
				assert(ptr==NULL);
				delete [] cstr;
			}
		} else if (flag.substr(0,7) == "--slice") {
			if (flag.length()==7 || flag[7]!='=') {
				std::cerr<<"Error: expected --slice=X[|Y|Z][,n].\n"<<std::endl;
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
		} else if (flag.substr(0,8) == "--coninv") {
			coninv = true;
		} else if (flag.substr(0,8) == "--contol") {
			if (flag.length()==8 || flag[8]!='=') {
				std::cerr<<"Error: expected --contol=a\n"<<std::endl;
				std::exit(-1);
			}
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			contol=(atof(ptr));
			ptr = strtok(NULL,",");
			assert(ptr==NULL);
			delete [] cstr;
		} else if (flag.substr(0,9) == "--contour") {
			if (flag.length()==9 || flag[9]!='=') {
				std::cerr<<"Error: expected --contour=n,a[,b,c,...]\n"<<std::endl;
				std::exit(-1);
			}
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			lvlfield=(atoi(ptr));
			ptr = strtok(NULL,",");
			while (ptr != NULL) {
				levelset.insert(atof(ptr));
				ptr = strtok(NULL,",");
			}
			delete [] cstr;
			if (lvlfield==-1) {
				std::cerr<<"Error: field "<<lvlfield<<" invalid.\n"<<std::endl;
				std::exit(-1);
			}
		} else if (flag.substr(0,7) == "--field") {
			if (flag[7]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(2,mode);
 			std::string val = flag.substr(8,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(NULL,",");
			}
			delete [] cstr;
		} else if (flag.substr(0,9) == "--exclude") {
			if (flag.length()==9 || flag[9]!='=') {
				std::cerr<<"Error: expected --exclude=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(3,mode);
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(NULL,",");
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
	if (dim < 1 or dim > 3) {
		std::cerr << "File input error: grid dimension must be 1, 2, or 3." << std::endl;
		exit(-1);
	}

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

	// if slice was not specified, handle it
	if (dim==3 && sliceaxis<0) {
		sliceaxis=2;
		if (slicelevel<x0[sliceaxis])
			slicelevel = (x1[sliceaxis] - x0[sliceaxis])/2;
	}

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
					MMSP::grid<1,bool> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,bool> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,bool> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,char> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,int> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,long> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,short> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,float> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,float> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,float> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,long double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,long double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,long double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
			else if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,double> GRID(argv[datindex]);
					scalar_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
				}
			}
		}
		else if (vector_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<bool> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<bool> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<bool> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<char> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<int> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<short> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<float> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<float> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<float> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<double> > GRID(argv[datindex]);
					vector_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<bool> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<bool> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<bool> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<char> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<int> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<short> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<float> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<float> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<float> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
			else if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<double> > GRID(argv[datindex]);
					sparse_field_to_png(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
				}
			}
		}

	if (invert)
		for (unsigned int n=0; n<theSize; n++)
			buffer[n] = 255-buffer[n];

	int result = MMSP::writePNG(image_size[0], image_size[1], 1,buffer, pngname.str().c_str());

	// clean up
	delete [] buffer;

	return result;
}

