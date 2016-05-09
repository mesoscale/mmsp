// mmsp2xyz.cpp
// Convert MMSP grid data to XYZ point cloud file format, readable e.g. by MeshLab for surface reconstruction
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<iostream>
#include<cstdlib>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<cmath>

#include"MMSP.hpp"

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol,
                                                    std::ofstream& xyzfil);

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    std::ofstream& xyzfil);

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    std::ofstream& xyzfil);

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] --contour=n,a[,b,...] [--normal] [--mag|--field|--exclude] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to XYZ point cloud format.\n";
		std::cout << "Usage:    " << argv[0] << " [--help] [--contour=n,a,b] [--normal] [--mag|--field=n|--exclude=n] infile [outfile]\n\n";
		std::cout << "Examples: " << argv[0] << " --help\n"
		          << "             displays this message.\n";
		std::cout << "          " << argv[0] << " --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             writes positions where field 0 equals 0.25,0.50,0.75 to infile.xyz.\n";
		std::cout << "          " << argv[0] << " --normal --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             writes positions and normals where field 0 equals 0.25,0.50,0.75 to infile.xyz.\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n";
		std::exit(0);
	}

	int datindex = 0; // in typical usage, input filename comes immediately after executable
	int xyzindex = 0; // in typical usage, output filename comes immediately after input
	int mode = 0;
	bool normal=false;
	double contol = 0.06; // contour tolerance (5-7% works OK)
	int lvlfield = -1; // store the field on which contour are calculated
	std::set<int> fieldset; // for --field or --exclude
	std::set<double> levelset; // for contours

	// Parse command line flags
	for (int i=1; i<argc; i++) {
		std::string flag(argv[i]);
		if (flag == "--mag") {
			mode=std::max(1,mode);
		} else if (flag == "--normal") {
			normal=true;
		} else if (flag.substr(0,8) == "--contol") {
			if (flag.length()==8 || flag[8]!='=') {
				std::cerr<< "Error: expected --contol=a\n"<<std::endl;
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
				std::cerr<< "Error: expected --contour=n,a[,b,c,...]\n"<<std::endl;
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
				std::cerr<< "Error: field "<<lvlfield<< " invalid.\n"<<std::endl;
				std::exit(-1);
			}
		} else if (flag.substr(0,7) == "--field") {
			if (flag[7]!='=') {
				std::cerr<< "Error: expected "<<flag<< "=n[,p,q,...]\n"<<std::endl;
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
				std::cerr<< "Error: expected --exclude=n[,p,q,...]\n"<<std::endl;
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
		} else if (flag.find(".xyz") != std::string::npos) {
			xyzindex = i;
		} else {
			datindex = i;
		}
	}

	// file open error check
	if (datindex==0) {
		std::cerr << "File input error: no input specified."<< ".\n";
		exit(-1);
	}
	std::ifstream input(argv[datindex]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[datindex] << ".\n";
		exit(-1);
	}

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

	// generate output file name
	std::stringstream xyzname;
	if (xyzindex==0) {
		std::string datname(argv[datindex]);
		int extpos = datname.find_last_of(".");
		if (datname.find_first_of("0123456789",extpos) != std::string::npos)
			xyzname << datname << ".xyz";
		else
			xyzname << datname.substr(0, extpos) << ".xyz";
	} else
		xyzname << argv[xyzindex];

	// Open output file
	std::ofstream xyzfil(xyzname.str().c_str());
	xyzfil << "#X\tY\tZ";
	if (normal)
		xyzfil << "nx\tny\tnz";
	xyzfil << '\n';

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

		// write grid data
		if (scalar_type or (not vector_type and not sparse_type)) { // must be scalar or built-in
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,bool> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,bool> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,bool> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned char> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned int> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned long> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned short> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,float> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,float> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,float> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,long double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,long double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,long double> GRID(argv[datindex]);
					convert_scalars(GRID, mode, normal, levelset, lvlfield, contol, xyzfil);
				}
			}
		}
		else if (vector_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<bool> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned char> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned int> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned long> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned short> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<float> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long double> > GRID(argv[datindex]);
					convert_vectors(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<bool> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<float> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long double> > GRID(argv[datindex]);
					convert_sparses(GRID, mode, normal, levelset, lvlfield, contol, fieldset, xyzfil);
				}
			}
		}

	return 0;
}

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, std::ofstream& xyzfil)
{
	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			bool excluded=true;
			T val = GRID(n);
			if (mode==1) //mag
				val = std::abs(val);
			if (levelset.size()!=0) //contour
				for (std::set<double>::iterator it=levelset.begin(); excluded && it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << "\t0\t0";
				if (normal) {
					MMSP::vector<T> norm = MMSP::gradient(GRID,x);
					xyzfil << '\t' << norm[0] << "\t0\t0";
				}
				xyzfil << '\n';
			}
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			bool excluded=true;
			T val = GRID(n);
			if (mode==1) //mag
				val = std::abs(val);
			if (levelset.size()!=0) //contour
				for (std::set<double>::iterator it=levelset.begin(); excluded && it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << "\t0";
				if (normal) {
					MMSP::vector<T> norm = MMSP::gradient(GRID,x);
					xyzfil << '\t' << norm[0] << '\t' << norm[1] << "\t0";
				}
				xyzfil << '\n';
			}
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			bool excluded=true;
			T val = GRID(n);
			if (mode==1) //mag
				val = std::abs(val);
			if (levelset.size()!=0) //contour
				for (std::set<double>::iterator it=levelset.begin(); excluded && it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << '\t' << dx(GRID,2)*x[2];
				if (normal) {
					MMSP::vector<T> norm = MMSP::gradient(GRID,x);
					xyzfil << '\t' << norm[0] << '\t' << norm[1] << '\t' << norm[2];
				}
				xyzfil << std::endl;
			}
		}
	}
}

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    std::ofstream& xyzfil)
{
	int included = (mode==2) ? fieldset.size() : MMSP::fields(GRID)-fieldset.size();

	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode==0) { //         default selection
				for (int i=0; i<MMSP::fields(GRID); i++) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			} else if (mode==1) { //  --mag
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
			}
			if (mode==1 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << "\t0\t0\n";
			}
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode==0) { //         default selection
				for (int i=0; i<MMSP::fields(GRID); i++) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			} else if (mode==1) { //  --mag
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
			}
			if (mode==1 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << "\t0\n";
			}
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode==0) { //         default selection
				for (int i=0; i<MMSP::fields(GRID); i++) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			} else if (mode==1) { //  --mag
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
			}
			if (mode==1 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << '\t' << dx(GRID,2)*x[2] << "\n";
			}
		}
	}
}

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID, const int& mode, const bool& normal,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    std::ofstream& xyzfil)
{
	int included = fieldset.size();

	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int h=0; h<GRID(n).length(); h++)
					sum += pow(GRID(n).value(h),2.0);
			} else if (mode==2) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
					if (included>1)
						sum += pow(GRID(n)[*it],2.0);
					else
						sum = GRID(n)[*it];
				}
			} else if (mode==3) { //  --exclude
				for (int h=0; h<MMSP::fields(GRID); h++) {
					int i=GRID(n).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(n).value(h),2.0);
				}
			}
			if (mode!=2 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << "\t0\t0\n";
			}
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int h=0; h<GRID(n).length(); h++)
					sum += pow(GRID(n).value(h),2.0);
			} else if (mode==2) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
					if (included>1)
						sum += pow(GRID(n)[*it],2.0);
					else
						sum = GRID(n)[*it];
				}
			} else if (mode==3) { //  --exclude
				for (int h=0; h<MMSP::fields(GRID); h++) {
					int i=GRID(n).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(n).value(h),2.0);
				}
			}
			if (mode!=2 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << "\t0\n";
			}
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int h=0; h<GRID(n).length(); h++)
					sum += pow(GRID(n).value(h),2.0);
			} else if (mode==2) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
					if (included>1)
						sum += pow(GRID(n)[*it],2.0);
					else
						sum = GRID(n)[*it];
				}
			} else if (mode==3) { //  --exclude
				for (int h=0; h<MMSP::fields(GRID); h++) {
					int i=GRID(n).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(n).value(h),2.0);
				}
			}
			if (mode!=2 || included!=1)
				sum = std::sqrt(sum);
			bool excluded=true;
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); excluded && itl!=levelset.end(); itl++)
					if (std::fabs(GRID(n)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						excluded=false;
			if (!excluded) {
				MMSP::vector<int> x=MMSP::position(GRID,n);
				xyzfil << dx(GRID,0)*x[0] << '\t' << dx(GRID,1)*x[1] << '\t' << dx(GRID,2)*x[2] << "\n";
			}
		}
	}
}
