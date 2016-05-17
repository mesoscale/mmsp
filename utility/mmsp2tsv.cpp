// mmsp2tsv.cpp
// Convert MMSP grid data to tab-delimited ASCII format (TSV)
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

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID, std::ofstream& tsvfil);

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID, std::ofstream& tsvfil);

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID, std::ofstream& tsvfil);

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] --contour=n,a[,b,...] [--normal] [--mag|--field|--exclude] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to tab-delimited ASCII (TSV) format.\n";
		std::cout << "Usage:    " << argv[0] << " [--help] infile [outfile]\n\n";
		std::cout << "Examples: " << argv[0] << " --help\n"
		          << "             displays this message.\n";
		std::cout << "          " << argv[0] << " input.dat\n"
		          << "             converts grid data from input.dat into input.tsv.\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n";
		std::exit(0);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n";
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
	std::stringstream tsvname;
	if (argc==2) {
		std::string datname(argv[1]);
		int extpos = datname.find_last_of(".");
		if (datname.find_first_of("0123456789",extpos) != std::string::npos)
			tsvname << datname << ".tsv";
		else
			tsvname << datname.substr(0, extpos) << ".tsv";
	} else
		tsvname << argv[2];

	// Open output file
	std::ofstream tsvfil(tsvname.str().c_str());

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
					MMSP::grid<1,bool> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,bool> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,bool> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned char> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned int> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned long> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,unsigned short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,unsigned short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,unsigned short> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,float> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,float> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,float> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,long double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,long double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,long double> GRID(argv[1]);
					convert_scalars(GRID, tsvfil);
				}
			}
		}
		else if (vector_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<bool> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<bool> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<bool> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned char> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned int> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned long> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<unsigned short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<unsigned short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<unsigned short> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<float> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<float> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<float> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::vector<long double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::vector<long double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::vector<long double> > GRID(argv[1]);
					convert_vectors(GRID, tsvfil);
				}
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<bool> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<bool> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<bool> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned char> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned int> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned long> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<unsigned short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<unsigned short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<unsigned short> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<float> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<float> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<float> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::sparse<long double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::sparse<long double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::sparse<long double> > GRID(argv[1]);
					convert_sparses(GRID, tsvfil);
				}
			}
		}

	return 0;
}

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID, std::ofstream& tsvfil)
{
	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			tsvfil << GRID(n) << '\n';
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			tsvfil << GRID(n) << '\n';
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			tsvfil << GRID(n) << '\n';
		}
	}
}

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID, std::ofstream& tsvfil)
{

	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int i=0; i<fields(GRID); i++)
				tsvfil << '\t' << GRID(n)[i];
			tsvfil << '\n';
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int i=0; i<fields(GRID); i++)
				tsvfil << '\t' << GRID(n)[i];
			tsvfil << '\n';
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int i=0; i<fields(GRID); i++)
				tsvfil << '\t' << GRID(n)[i];
			tsvfil << '\n';
		}
	}
}

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID, std::ofstream& tsvfil)
{
	if (dim==1) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int h=0; h<length(GRID(n)); h++)
				tsvfil << '\t' << GRID(n).index(h) << '\t' << GRID(n).value(h);
			tsvfil << '\n';
		}
	} else if (dim==2) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int h=0; h<length(GRID(n)); h++)
				tsvfil << '\t' << GRID(n).index(h) << '\t' << GRID(n).value(h);
			tsvfil << '\n';
		}
	} else if (dim==3) {
		for (int n=0; n<MMSP::nodes(GRID); n++) {
			MMSP::vector<int> x=MMSP::position(GRID,n);
			tsvfil << dx(GRID,0)*x[0];
			for (int d=1; d<dim; d++)
				tsvfil << '\t' << dx(GRID,d)*x[d];
			for (int h=0; h<length(GRID(n)); h++)
				tsvfil << '\t' << GRID(n).index(h) << '\t' << GRID(n).value(h);
			tsvfil << '\n';
		}
	}
}
