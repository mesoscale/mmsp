// mmsp2vti.cpp
// Convert MMSP grid data to VTK image data format
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MMSP.hpp"
#include<zlib.h>
#include<sstream>
#include<cmath>
#include<vector>
#include<map>

bool fexists(const char *filename) {
	std::ifstream ifile(filename);
	return ifile;
}

int main(int argc, char* argv[]) {
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] infile [outfile]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to VTK image data format.\n";
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
	#ifdef DEBUG
	std::cout<<"Reading "<<argv[1]<<std::endl;
	#endif

	// generate output file name
	std::stringstream filename;
	if (argc < 3)
		filename << std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of(".")) << ".vti";
	else
		filename << argv[2];

	std::ofstream output(filename.str().c_str());
	if (!output) {
		std::cerr << "File output error: could not open ";
		std::cerr << filename.str() << "." << std::endl;
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input, type, '\n');
	#ifdef DEBUG
	std::cout<<"Grid type is "<<type<<std::endl;
	#endif

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
	if (not sparse_type) {
		std::cerr << "File input error: sparse data expected." << std::endl;
	}

	// read grid dimension
	int dim;
	input >> dim;
	#ifdef DEBUG
	std::cout<<"Grid is "<<dim<<"-dimensional."<<std::endl;
	#endif

	// read number of fields
	int fields;
	input >> fields;
	#ifdef DEBUG
	std::cout<<"Grid has "<<fields<<" fields."<<std::endl;
	#endif

	// read grid sizes
	int g0[3] = {0, 0, 0};
	int g1[3] = {0, 0, 0};
	for (int i = 0; i < dim; i++)
		input >> g0[i] >> g1[i];
	#ifdef DEBUG
	std::cout<<"Grid edge is "<<g1[0] - g0[0]<<std::endl;
	#endif

	// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int i = 0; i < dim; i++)
		input >> dx[i];
	#ifdef DEBUG
	std::cout<<"Grid spacing is "<<dx[0]<<std::endl;
	#endif

	// ignore trailing endlines
	input.ignore(10, '\n');


	// determine byte order
	std::string byte_order;
	if (0x01 & static_cast<int>(1)) byte_order = "LittleEndian";
	else byte_order = "BigEndian";
	#ifdef DEBUG
	std::cout<<"Grid is "<<byte_order<<std::endl;
	#endif

	// output header markup
	output << "<?xml version=\"1.0\"?>\n";
	output << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << byte_order << "\">\n";

	if (dim == 1) {
		output << "  <ImageData WholeExtent=\"" << g0[0] << " " << g1[0] << " 0 0 0 0\"";
		output << "   Origin=\"0 0 0\" Spacing=\"" << dx[0] << " 1 1\">\n";
	}
	if (dim == 2) {
		output << "  <ImageData WholeExtent=\"" << g0[1] << " " << g1[1] << " " << g0[0] << " " << g1[0] << " 0 0\"";
		output << "   Origin=\"0 0 0\" Spacing=\"" << dx[1] << " " << dx[0] << " 1\">\n";
	}
	if (dim == 3) {
		output << "  <ImageData WholeExtent=\"" << g0[2] << " " << g1[2] << " " << g0[1] << " " << g1[1] << " " << g0[0] << " " << g1[0] << "\"";
		output << "   Origin=\"0 0 0\" Spacing=\"" << dx[2] << " " << dx[1] << " " << dx[0] << "\">\n";
	}


	// Estimate number of grains, for color randomization
	int est_grains = 10000;
	if (dim==2) est_grains=static_cast<int>(1.5*float((g1[0]-g0[0])*(g1[1]-g0[1]))/(M_PI*10.*10.)); // average grain is a disk of radius 10
	else if (dim==3) est_grains=static_cast<int>(1.5*float((g1[0]-g0[0])*(g1[1]-g0[1])*(g1[2]-g0[2]))/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
	#ifdef DEBUG
	std::cout<<"Grid contains approx. "<<est_grains<<" grains."<<std::endl;
	#endif
	std::vector<int> colors;
	for (unsigned int i=0; i<est_grains; i++)
		colors.push_back(rand() % est_grains);

	// read number of blocks
	int blocks;
	input.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));

	for (int i = 0; i < blocks; i++) {
		#ifdef DEBUG
		std::cout<<"  Reading block "<<i+1<<" of "<<blocks<<std::endl;
		#endif
		// read block limits
		int lmin[3] = {0, 0, 0};
		int lmax[3] = {0, 0, 0};
		for (int j = 0; j < dim; j++) {
			input.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
			input.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
		}
		#ifdef DEBUG
		std::cout<<"  Block edge is "<<lmax[0] - lmin[0]<<std::endl;
		#endif
		int blo[dim];
    int bhi[dim];
    // read boundary conditions
    for (int j = 0; j < dim; j++) {
      input.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
      input.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
    }

		// write header markup
		if (dim == 1) output << "    <Piece Extent=\"" << lmin[0] << " " << lmax[0] << " 0 0 0 0\">\n";
		if (dim == 2) output << "    <Piece Extent=\"" << lmin[1] << " " << lmax[1] << " " << lmin[0] << " " << lmax[0] << " 0 0\">\n";
		if (dim == 3) output << "    <Piece Extent=\"" << lmin[2] << " " << lmax[2] << " " << lmin[1] << " " << lmax[1] << " " << lmin[0] << " " << lmax[0] << "\">\n";

		// write cell data markup
		if (scalar_type) {
			output << "      <CellData>\n";
			output << "        <DataArray Name=\"scalar_data\"";
		}

		else if (vector_type) {
			output << "      <CellData>\n";
			output << "        <DataArray Name=\"vector_data\" NumberOfComponents=\"" << fields << "\"";
		}

		else if (sparse_type) {
			output << "      <CellData>\n";
			output << "        <DataArray Name=\"grain_id\"";
		}

		else { /* built-in data types */
			output << "      <CellData>\n";
			output << "        <DataArray Name=\"scalar_data\"";
		}

		if (bool_type)
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
		else if (double_type && !sparse_type)
			output << " type=\"Float64\" format=\"ascii\">\n";
		else if (double_type && sparse_type)
			output << " type=\"UInt16\" format=\"ascii\">\n";
		else if (long_double_type)
			output << " type=\"Float128\" format=\"ascii\">\n";

		// read grid data
		unsigned long size, rawSize;
		input.read(reinterpret_cast<char*>(&rawSize), sizeof(rawSize)); // read raw size
		input.read(reinterpret_cast<char*>(&size), sizeof(size)); // read compressed size
		char* compressed_buffer = new char[size];
		input.read(compressed_buffer, size);
		#ifdef DEBUG
		std::cout<<"  Read "<<size<<" B, compressed data."<<std::endl;
		#endif
		char* buffer;
		if (size!=rawSize) {
			// Decompress data
			buffer = new char[rawSize];
			int status;
			status = uncompress(reinterpret_cast<unsigned char*>(buffer), &rawSize, reinterpret_cast<unsigned char*>(compressed_buffer), size);
			switch( status ) {
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
		} else {
			buffer=compressed_buffer;
			compressed_buffer=NULL;
		}

		// write grid data
		if (sparse_type) {
			if (bool_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						bool sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<bool>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						bool sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<bool>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<bool> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						bool sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<bool>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (unsigned_char_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned char> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned char sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned char>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (unsigned_int_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned int sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned int>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (unsigned_long_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned long> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned long sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned long>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			if (unsigned_short_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<unsigned short> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						unsigned short sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<unsigned short>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
			// === FLOAT ===
			if (float_type) {
				#ifdef DEBUG
				std::cout<<"  Writing grain IDs from sparse floats."<<std::endl;
				#endif
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<float> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << GRID(k).grain_id() << " ";
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<float> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k).grain_id()%est_grains] << " ";
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<float> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					#ifdef DEBUG
					std::cout<<"  Opened 3D grid from buffer."<<std::endl;
					#endif
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						#ifdef DEBUG
						assert(GRID(k).grain_id()%est_grains < est_grains);
						#endif
						output << colors[GRID(k).grain_id()%est_grains] << " ";
					}
				}
			}
			// === DOUBLE ===
			if (double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k).grain_id()%est_grains] << " ";
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k).grain_id()%est_grains] << " ";
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k).grain_id()%est_grains] << " ";
				}
			}
			if (long_double_type) {
				if (dim == 1) {
					MMSP::grid<1, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long double sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long double>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 2) {
					MMSP::grid<2, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long double sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long double>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				} else if (dim == 3) {
					MMSP::grid<3, MMSP::sparse<long double> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++) {
						long double sum = 0;
						int nonzero = MMSP::length(GRID(k));
						for (int h = 0; h < nonzero; h++)
							sum += static_cast<long double>(GRID(k).value(h) * GRID(k).value(h));
						output << sum << " ";
					}
				}
			}
		} else if (not scalar_type and not vector_type) {
			if (int_type) {
				if (dim == 1) {
					MMSP::grid<1,MMSP::scalar<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k)%est_grains] << " ";
				} else if (dim == 2) {
					MMSP::grid<2,MMSP::scalar<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k)%est_grains] << " ";
				} else if (dim == 3) {
					MMSP::grid<3,MMSP::scalar<int> > GRID(fields, lmin, lmax);
					GRID.from_buffer(buffer);
					for (int k = 0; k < MMSP::nodes(GRID); k++)
						output << colors[GRID(k)%est_grains] << " ";
				}
			}
		}

		// clean up
		delete [] buffer;

		// write closing markup
		output << "\n";
		output << "        </DataArray>\n";
		output << "      </CellData>\n";
		output << "    </Piece>\n";
	}

	// output closing markup
	output << "  </ImageData>\n";
	output << "</VTKFile>\n";
}

