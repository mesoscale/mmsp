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

template <int dim, typename T> void scalar_limits(const MMSP::grid<dim,MMSP::scalar<T> >& GRID, T& min, T& max)
{
	min=0;
	max=1;
	for (int n=0; n<MMSP::nodes(GRID); n++) {
		if (GRID(n)>max)
			max=GRID(n);
		else if (GRID(n)<min)
			min=GRID(n);
	}
}
template <int dim, typename T> void vector_limits(const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                  const int& mode, const std::set<int>& fieldset, T& min, T& max)
{
	min=0;
	max=1;
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
				if (it == fieldset.end())
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
			}
		} else if (mode==4) { //  --contour
			// Same as --exclude
			for (int i=0; i<MMSP::fields(GRID); i++) {
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
			}
		}
		if (mode < 2 || fieldset.size()>1 || included>1)
			sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}
}
template <int dim, typename T> void sparse_limits(const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                  const int& mode, const std::set<int>& fieldset, T& min, T& max)
{
	min=0;
	max=1;

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
					sum += pow(GRID(n)[i],2.0);
			}
		} else if (mode==4) { //  --contour
			// Same as --exclude
			for (int h=0; h<MMSP::fields(GRID); h++) {
				int i=GRID(n).index(i);
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					sum += pow(GRID(n)[i],2.0);
			}
		}
		sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}
}

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] [--mag|--max] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to VTK image data format.\n";
		std::cout << "Usage: " << argv[0] << " [--help] [--mag|--max|--field=n] infile [outfile]\n";
		std::cout << "       Select either --mag or --max to flatten vector or sparse data by the specified method.\n";
		std::cout << "Questions/comments to gruberja@gmail.com (Jason Gruber).\n";
		std::exit(0);
	}

	int datindex = 0; // in typical usage, input filename comes immediately after executable
	int pngindex = 0; // in typical usage, output filename comes immediately after input
	int mode = 0;
	bool invert=false;
	int slice = -1;
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
				case 'x': slice=0;
				          break;
				case 'y': slice=1;
				          break;
				case 'z': slice=2;
				          break;
				case 'X': slice=0;
				          break;
				case 'Y': slice=1;
				          break;
				case 'Z': slice=2;
				          break;
			}
			if (slice==-1) {
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
			mode=2;
			std::string val = flag.substr(8,255);
			char* cstr = strtok(val.c_str(),",");
			while (cstr != NULL) {
				fieldset.insert(atoi(cstr));
				cstr = strtok(val.c_str(),",");
			}
		} else if (flag.substr(0,9) == "--exclude") {
			if (flag[9]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=3;
			std::string val = flag.substr(10,255);
			char* cstr = strtok(val.c_str(),",");
			while (cstr != NULL) {
				fieldset.insert(atoi(cstr));
				cstr = strtok(val.c_str(),",");
			}
		} else if (flag.substr(0,9) == "--contour") {
			if (flag[9]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n,a[,b,c,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=4;
			std::string val = flag.substr(10,255);
			char* cstr = strtok(val.c_str(),",");
			fieldset.insert(atoi(cstr));
			cstr = strtok(val.c_str(),",");
			while (cstr != NULL) {
				levelset.insert(atof(cstr));
				cstr = strtok(val.c_str(),",");
			}
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
	if (pngindex==0)
		pngname << std::string(argv[datindex]).substr(0, std::string(argv[datindex]).find_last_of(".")) << ".png";
	else
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
		input >> x0[i] >> x1[i];

// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int d=0; d<dim; d++)
		input >> dx[i];

// ignore trailing endlines
	input.ignore(10, '\n');

	int image_size[2] = {1,1};
	int i=0;
	for (int d=0; d<dim; d++) {
		if (d==slice)
			continue;
		image_size[i] = x1[d]-x0[d];
		i++;
	}
  unsigned int theSize=image_size[0] * image_size[1];
  unsigned char* buffer = new unsigned char[theSize];

		// write grid data
		if (not vector_type and not sparse_type) { // must be scalar or built-in
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<bool> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
						buffer[n] = GRID(x);
						n++;
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<bool> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
							buffer[n] = GRID(x);
							n++;
						}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<bool> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
								buffer[n] = GRID(x);
								n++;
							}
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<char> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
						buffer[n] = GRID(x);
						n++;
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<char> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
							buffer[n] = GRID(x);
							n++;
						}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<char> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
								buffer[n] = GRID(x);
								n++;
							}
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					unsigned int n=0;
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++) {
						buffer[n] = GRID(x);
						n++;
					{
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x)<< " ";
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x)<< " ";
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x)<< " ";
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x)<< " ";
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x)<< " ";
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x)<< " ";
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::scalar<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						buffer[n] = GRID(x)<< " ";
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::scalar<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							buffer[n] = GRID(x)<< " ";
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::scalar<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								buffer[n] = GRID(x)<< " ";
				}
			}
		}

		else if (vector_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::vector<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::vector<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_vector(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::vector<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_vector(output, GRID(x), flatten, field);
				}
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<bool> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned char> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned int> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned long> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned short> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (float_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<float> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(1, 0);
					for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
						scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(2, 0);
					for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
						for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
							scale_sparse(output, GRID(x), flatten, field);
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long double> > GRID(argv[datindex]);
					
					MMSP::vector<int> x(3, 0);
					for (x[2] = x0[2]; x[2] < x1[2]; x[2]++)
						for (x[1] = x0[1]; x[1] < x1[1]; x[1]++)
							for (x[0] = x0[0]; x[0] < x1[0]; x[0]++)
								scale_sparse(output, GRID(x), flatten, field);
				}
			}
		}

	int result = writePNG(image_size[0], image_size[1], 1, buffer, pngname.c_str());

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

