// mmsp2pvd.cpp
// Convert MMSP grid data to ParaView PVD file format
// Note: results in PVD file and a sequence of VTK image data files
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include<sstream>
#include"MMSP.hpp"

template<int dim, typename T> void print_scalars(std::ofstream& fstr, const MMSP::grid<dim,T>& GRID, const int& mode)
{
	if (dim==1) {
		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			if (mode==1) { // --mag
				fstr << std::sqrt(GRID(x)*GRID(x)) << " ";
			} else {
				fstr << GRID(x) << " ";
			}
		}
	} else if (dim==2) {
		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				if (mode==1) { // --mag
					fstr << std::sqrt(GRID(x)*GRID(x)) << " ";
				} else {
					fstr << GRID(x) << " ";
				}
			}
		}
	} else if (dim==3) {
		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					if (mode==1) { // --mag
						fstr << std::sqrt(GRID(x)*GRID(x)) << " ";
					} else {
						fstr << GRID(x) << " ";
					}
				}
			}
		}
	}
}

template<int dim, typename T> void print_vectors(std::ofstream& fstr, const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                const int& mode, const int& field)
{
	if (dim==1) {
		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			const MMSP::vector<T>& v = GRID(x);
			if (mode==1) { // --mag
				double sum = 0.0;
				for (int h = 0; h < v.length(); h++)
					sum += v[h]*v[h];
				fstr << std::sqrt(sum) << " ";
			} else if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < v.length(); h++)
					if (v[h] > v[max])
						max = h;
				fstr << max << " ";
			} else if (mode==3) { // --field
				fstr << v[field] << " ";
			} else {
				for (int h = 0; h < v.length(); h++)
					fstr << v[h] << " ";
			}
		}
	} else if (dim==2) {
		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				const MMSP::vector<T>& v = GRID(x);
				if (mode==1) { // --mag
					double sum = 0.0;
					for (int h = 0; h < v.length(); h++)
						sum += v[h]*v[h];
					fstr << std::sqrt(sum) << " ";
				} else if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < v.length(); h++)
						if (v[h] > v[max])
							max = h;
					fstr << max << " ";
				} else if (mode==3) { // --field
					fstr << v[field] << " ";
				} else {
					for (int h = 0; h < v.length(); h++)
						fstr << v[h] << " ";
				}
			}
		}
	} else if (dim==3) {
		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					const MMSP::vector<T>& v = GRID(x);
					if (mode==1) { // --mag
						double sum = 0.0;
						for (int h = 0; h < v.length(); h++)
							sum += v[h]*v[h];
						fstr << std::sqrt(sum) << " ";
					} else if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < v.length(); h++)
							if (v[h] > v[max])
								max = h;
						fstr << max << " ";
					} else if (mode==3) { // --field
						fstr << v[field] << " ";
					} else {
						for (int h = 0; h < v.length(); h++)
							fstr << v[h] << " ";
					}
				}
			}
		}
	}
}

template<int dim, typename T> void print_sparses(std::ofstream& fstr, const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                const int& mode, const int& field)
{
	if (dim==1) {
		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			const MMSP::sparse<T>& s = GRID(x);
			if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < s.length(); h++)
					if (s.value(h) > s.value(max))
						max = h;
				fstr << s.index(max) << " ";
			} else if (mode==3) { // --field
				fstr << s[field] << " ";
			} else { // --mag is redundant for sparse
				double sum = 0.0;
				for (int h = 0; h < s.length(); h++)
					sum += s.value(h)*s.value(h);
				fstr << std::sqrt(sum) << " ";
			}
		}
	} else if (dim==2) {
		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				const MMSP::sparse<T>& s = GRID(x);
				if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < s.length(); h++)
						if (s.value(h) > s.value(max))
							max = h;
					fstr << s.index(max) << " ";
				} else if (mode==3) { // --field
					fstr << s[field] << " ";
				} else { // --mag is redundant for sparse
					double sum = 0.0;
					for (int h = 0; h < s.length(); h++)
						sum += s.value(h)*s.value(h);
					fstr << std::sqrt(sum) << " ";
				}
			}
		}
	} else if (dim==3) {
		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					const MMSP::sparse<T>& s = GRID(x);
					if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < s.length(); h++)
							if (s.value(h) > s.value(max))
								max = h;
						fstr << s.index(max) << " ";
					} else if (mode==3) { // --field
						fstr << s[field] << " ";
					} else { // --mag is redundant for sparse
						double sum = 0.0;
						for (int h = 0; h < s.length(); h++)
							sum += s.value(h)*s.value(h);
						fstr << std::sqrt(sum) << " ";
					}
				}
			}
		}
	}
}

int main(int argc, char* argv[]) {
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] [--mag|--max|--field=n] infile1 [infile2 [infile3 ... ]] [[--output=]outfile]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to ParaView Data timeseries format.\n";
		std::cout << "Output: PVD index file and time-sequence of VTI data files. Open the PVD (alone) in ParaView.\n";
		std::cout << "Usage: " << argv[0] << " [--help] [--mag|--max|--field=n] [[--output=]outfile] infile1 [infile2 [infile3 ... ]]\n\n";
		std::cout << "       Select either --mag or --max to flatten vector or sparse data by the specified method.\n";
		std::cout << "       Select --field=n with vector or sparse data to extract the value of field n alone.\n";
		std::cout << "Questions/comments to gruberja@gmail.com (Jason Gruber).\n\n";
		exit(0);
	}

	// parse command-line options
	int dat_idx = 1; // in typical usage, filename comes immediately after executable
	int flatten = 0;
	int field = 0;

	if (std::string(argv[1]) == "--mag") {
		flatten=1;
		dat_idx=2;
	} else if (std::string(argv[1]) == "--max") {
		flatten=2;
		dat_idx=2;
	} else if (std::string(argv[1]).substr(0,8) == "--field=") {
		flatten=3;
		dat_idx=2;
		std::string str(argv[1]);
		field = atoi(str.substr(8,12).c_str());
	}

	int pvd_idx=-1; // output file may not be on the command line
	std::ifstream input;
	std::stringstream filename;
	bool flagged=false;
	// Look for --output flag.
	for (int i=dat_idx; i<argc && pvd_idx<0; i++) {
	    if (std::string(argv[i]).substr(0, 9) == "--output=") {
	    	flagged=true;
			filename << std::string(argv[i]).substr(9);
			pvd_idx = i;
			dat_idx = i+1;
	    }
	}
	// Perhaps --output was not specified. Try opening each file, and assign output to the one that fails.
	for (int i=dat_idx; i<argc && pvd_idx<0; i++) {
		input.open(argv[i]);
		if (!input) {
			// Found our PVD name
			pvd_idx = i;
			dat_idx = i+1;
			std::string fstr(argv[i]);
			unsigned int start=0;
			if (fstr.find_last_of("=")!=std::string::npos)
				start=fstr.find_last_of("=");
			filename << fstr.substr(start, fstr.find_last_of(".")) << ".pvd";
		} else {
			input.close();
		}
	}

	if (pvd_idx<0) {// not found on command line
		filename << std::string(argv[dat_idx]).substr(0, std::string(argv[dat_idx]).find_last_of(".")) << ".pvd";
	} else { // found
		if (flagged)
			std::cout<<"PVD name specified as argv["<<pvd_idx<<"]: writing to "<<filename.str()<<".\n";
		else
			std::cout<<"PVD name inferred as argv["<<pvd_idx<<"]: writing to "<<filename.str()<<".\n";
	}

	// file open error check
	std::ofstream pvdfile(filename.str().c_str());
	if (!pvdfile) {
		std::cerr << "File output error: could not open "
		          << filename.str() << ", check permissions." << std::endl;
		exit(-1);
	}

	// determine byte order
	std::string byte_order;
	if (0x01 & static_cast<int>(1)) byte_order = "LittleEndian";
	else byte_order = "BigEndian";

	// write PVD file header
	pvdfile << "<?xml version=\"1.0\"?>\n";
	pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << byte_order << "\">\n";
	pvdfile << "  <Collection>\n";


	// process each file on the command line
	for (int if_idx = dat_idx; if_idx < argc; if_idx++) {
		if (if_idx==pvd_idx)
			continue;
		filename.str(""); // clear stream for new input
		// file open error check
		input.open(argv[if_idx]);
		if (!input) {
			std::cerr << "File input error: could not open " << argv[if_idx] << ".\n\n";
			exit(-1);
		}

		// generate output file name
		filename << std::string(argv[if_idx]).substr(0, std::string(argv[if_idx]).find_last_of(".")) << ".vti";

		// file open error check
		std::ofstream output(filename.str().c_str());
		if (!output) {
			std::cerr << "File output error: could not open ";
			std::cerr << filename.str() << "." << std::endl;
			exit(-1);
		}

		// write PVD file data
		pvdfile << "    <DataSet timestep=\"" << if_idx - dat_idx << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>\n";

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
		for (int i = 0; i < dim; i++)
			input >> x0[i] >> x1[i];

		// read cell spacing
		double dx[3] = {1.0, 1.0, 1.0};
		for (int i = 0; i < dim; i++)
			input >> dx[i];

		// ignore trailing endlines
		input.ignore(10, '\n');

		// output header markup
		output << "<?xml version=\"1.0\"?>\n";
		output << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << byte_order << "\">\n";
		if (dim == 1) {
			output << "  <ImageData WholeExtent=\"" << x0[0] << " " << x1[0] << " 0 0 0 0\"";
			output << "   Origin=\"0 0 0\" Spacing=\"" << dx[0] << " 1 1\">\n";
		}
		if (dim == 2) {
			output << "  <ImageData WholeExtent=\"" << x0[0] << " " << x1[0] << " "
			                                        << x0[1] << " " << x1[1] << " 0 0\"";
			output << "   Origin=\"0 0 0\" Spacing=\"" << dx[0] << " " << dx[1] << " 1\">\n";
		}
		if (dim == 3) {
			output << "  <ImageData WholeExtent=\"" << x0[0] << " " << x1[0] << " "
			                                        << x0[1] << " " << x1[1] << " "
			                                        << x0[2] << " " << x1[2] << "\"";
			output << "   Origin=\"0 0 0\" Spacing=\"" << dx[0] << " " << dx[1] << " " << dx[2] << "\">\n";
		}

		// read number of blocks
	    int blocks;
        input.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));

        for (int i = 0; i < blocks; i++) {
			// read local block limits
			int lmin[3] = {0, 0, 0};
			int lmax[3] = {0, 0, 0};
			for (int j = 0; j < dim; j++) {
				input.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
				input.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
			}
			int blo[3] = {0, 0, 0};
			int bhi[3] = {0, 0, 0};
			// read boundary conditions
			for (int j = 0; j < dim; j++) {
				input.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
				input.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
			}

			// write header markup
			if (dim == 1)
				output << "    <Piece Extent=\"" << lmin[0] << " " << lmax[0] << " 0 0 0 0\">\n";
			if (dim == 2)
				output << "    <Piece Extent=\"" << lmin[0] << " " << lmax[0] << " "
				                                 << lmin[1] << " " << lmax[1] << " 0 0\">\n";
			if (dim == 3)
				output << "    <Piece Extent=\"" << lmin[0] << " " << lmax[0] << " "
				                                 << lmin[1] << " " << lmax[1] << " "
				                                 << lmin[2] << " " << lmax[2] << "\">\n";

			// write cell data markup
			if (scalar_type || flatten>0) {
				output << "      <CellData>\n";
				output << "        <DataArray Name=\"scalar_data\"";
			}

			else if (vector_type) {
				output << "      <CellData>\n";
				output << "        <DataArray Name=\"vector_data\" NumberOfComponents=\"" << fields << "\"";
			}

			else if (sparse_type) {
				output << "      <CellData>\n";
				output << "        <DataArray Name=\"scalar_data\"";
			}

			else { /* built-in data types */
				output << "      <CellData>\n";
				output << "        <DataArray Name=\"scalar_data\"";
			}

			if (flatten==1) // mag
				output << " type=\"Float32\" format=\"ascii\">\n";
			else if (flatten==2) // max
				output << " type=\"Int32\" format=\"ascii\">\n";
			else if (bool_type)
				output << " type=\"UInt8\" format=\"ascii\">\n";
			else if (char_type)
				output << " type=\"Int8\" format=\"ascii\">\n";
			else if (unsigned_char_type)
				output << " type=\"UInt8\" format=\"ascii\">\n";
			else if (int_type)
				output << " type=\"Int32\" format=\"ascii\">\n";
			else if (unsigned_int_type)
				output << " type=\"UInt32\" format=\"ascii\">\n";
			else if (long_type)
				output << " type=\"Int32\" format=\"ascii\">\n";
			else if (unsigned_long_type)
				output << " type=\"UInt32\" format=\"ascii\">\n";
			else if (short_type)
				output << " type=\"Int16\" format=\"ascii\">\n";
			else if (unsigned_short_type)
				output << " type=\"UInt16\" format=\"ascii\">\n";
			else if (float_type)
				output << " type=\"Float32\" format=\"ascii\">\n";
			else if (double_type)
				output << " type=\"Float64\" format=\"ascii\">\n";
			else if (long_double_type)
				output << " type=\"Float128\" format=\"ascii\">\n";


			// read grid data
			unsigned long size, rawSize;
			input.read(reinterpret_cast<char*>(&rawSize), sizeof(rawSize)); // read raw size
			input.read(reinterpret_cast<char*>(&size), sizeof(size)); // read compressed size
			char* compressed_buffer = new char[size];
			input.read(compressed_buffer, size);
			char* buffer = NULL;
			if (size != rawSize) {
				// Decompress data
				buffer = new char[rawSize];
				int status;
				status = uncompress(reinterpret_cast<unsigned char*>(buffer), &rawSize, reinterpret_cast<unsigned char*>(compressed_buffer), size);
				switch(status) {
				case Z_OK:
					break;
				case Z_MEM_ERROR:
					std::cerr << "Uncompress: out of memory." << std::endl;
					exit(1);
					break;
				case Z_BUF_ERROR:
					std::cerr << "Uncompress: output buffer wasn't large enough." << std::endl;
					exit(1);
					break;
				}
				delete [] compressed_buffer;
				compressed_buffer = NULL;
			} else {
				buffer = compressed_buffer;
				compressed_buffer = NULL;
			}

			// write grid data
			if (not vector_type and not sparse_type) { // must be scalar or built-in
				if (bool_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (unsigned_char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
							print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (unsigned_int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (unsigned_long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (unsigned_short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (float_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
				if (long_double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::scalar<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::scalar<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::scalar<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_scalars(output, GRID, flatten);
					}
				}
			}

			else if (vector_type) {
				if (bool_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (unsigned_char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (unsigned_int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (unsigned_long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (unsigned_short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (float_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
				if (long_double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::vector<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::vector<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::vector<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_vectors(output, GRID, flatten, field);
					}
				}
			}

			else if (sparse_type) {
				if (bool_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (unsigned_char_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (unsigned_int_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (unsigned_long_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (unsigned_short_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (float_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<float> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
				if (long_double_type) {
					if (dim == 1) {
						MMSP::grid<1, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 2) {
						MMSP::grid<2, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					} else if (dim == 3) {
						MMSP::grid<3, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
						GRID.from_buffer(buffer);
						print_sparses(output, GRID, flatten, field);
					}
				}
			}

			// clean up
			delete [] buffer;

			// output closing markup
			output << "\n";
			output << "        </DataArray>\n";
			output << "      </CellData>\n";
			output << "    </Piece>\n";
		} // loop for blocks
		output << "  </ImageData>\n";
		output << "</VTKFile>\n";
		input.close();
	} // loop for files

	// output closing PVD markup
	pvdfile << "  </Collection>\n";
	pvdfile << "</VTKFile>\n";
	pvdfile.close();

	return 0;
}

