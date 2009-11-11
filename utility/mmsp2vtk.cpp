// mmsp2vtk.cpp
// Convert MMSP grid data format to VTK image data format
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"MMSP.hpp"

int main(int argc, char* argv[])
{
	// command line error check
	if (argc<3) {
		std::cout<<"Usage: "<<argv[0]<<" inputfile outputfile\n";
		exit(-1);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open ";
		std::cerr<<argv[1]<<"."<<std::endl;
		exit(-1);
	}

	// file open error check
	std::ofstream output(argv[2]);
	if (!output) {
		std::cerr<<"File output error: could not open ";
		std::cerr<<argv[2]<<"."<<std::endl;
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input,type,'\n');

	// grid type error check
	if (type.substr(0,4)!="grid") {
		std::cerr<<"File input error: file does not contain grid data."<<std::endl;
		exit(-1);
	}

	// parse data type
	bool bool_type = (type.find("bool")!=std::string::npos);
	bool char_type = (type.find("char")!=std::string::npos);
	bool unsigned_char_type = (type.find("unsigned char")!=std::string::npos);
	bool int_type = (type.find("int")!=std::string::npos);
	bool unsigned_int_type = (type.find("unsigned int")!=std::string::npos);
	bool long_type = (type.find("long")!=std::string::npos);
	bool unsigned_long_type = (type.find("unsigned long")!=std::string::npos);
	bool short_type = (type.find("short")!=std::string::npos);
	bool unsigned_short_type = (type.find("unsigned short")!=std::string::npos);
	bool float_type = (type.find("float")!=std::string::npos);
	bool double_type = (type.find("double")!=std::string::npos);
	bool long_double_type = (type.find("long double")!=std::string::npos);

	bool scalar_type = (type.find("scalar")!=std::string::npos);
	bool vector_type = (type.find("vector")!=std::string::npos);
	bool sparse_type = (type.find("sparse")!=std::string::npos);

	if (not bool_type    and
        not char_type    and  not unsigned_char_type   and
	    not int_type     and  not unsigned_int_type    and
	    not long_type    and  not unsigned_long_type   and
	    not short_type   and  not unsigned_short_type  and
	    not float_type   and
        not double_type  and  not long_double_type) {
		std::cerr<<"File input error: unknown grid data type."<<std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input>>dim;

	// read number of fields 
	int fields;
	input>>fields;

	// read grid sizes
	int x0[3] = {0,0,0};
	int x1[3] = {0,0,0};
	for (int i=0; i<dim; i++)
		input>>x0[i]>>x1[i];

	// read cell spacing
	float dx[3] = {1.0,1.0,1.0};
	for (int i=0; i<dim; i++)
		input>>dx[i];

	// ignore trailing endlines
	input.ignore(10,'\n');

	// compute slice size
	int slice = 1;
	for (int i=1; i<dim; i++)
		slice *= (x1[i]-x0[i]);

	// compute number of cells
	int cells = 1;
	for (int i=0; i<dim; i++)
		cells *= (x1[i]-x0[i]);


	// determine byte order
	std::string byte_order;
	if (0x01 & static_cast<int>(1)) byte_order = "LittleEndian";
	else byte_order = "BigEndian";

	// output header markup
	output<<"<?xml version=\"1.0\"?>\n";
	output<<"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\""<<byte_order<<"\">\n";
	output<<"  <ImageData WholeExtent=\""<<x0[0]<<" "<<x1[0]<<" "<<x0[1]<<" "<<x1[1]<<" "<<x0[2]<<" "<<x1[2]<<"\"";
	output<<"   Origin=\"0 0 0\" Spacing=\""<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<"\">\n";
	output<<"    <Piece Extent=\""<<x0[0]<<" "<<x1[0]<<" "<<x0[1]<<" "<<x1[1]<<" "<<x0[2]<<" "<<x1[2]<<"\">\n";

	// output cell data markup
	if (scalar_type) {
		output<<"      <CellData Scalars=\"scalar_data\">\n";
		output<<"        <DataArray Name=\"scalar_data\"";
	}

	else if (vector_type) {
		output<<"      <CellData Vectors=\"vector_data\">\n";
		output<<"        <DataArray Name=\"vector_data\" NumberOfComponents=\""<<fields<<"\"";
	}

	else if (sparse_type) {
		output<<"      <CellData Scalars=\"scalar_data\">\n";
		output<<"        <DataArray Name=\"scalar_data\"";
	}

	else /* built-in data types */ {
		output<<"      <CellData Scalars=\"scalar_data\">\n";
		output<<"        <DataArray Name=\"scalar_data\"";
	}

	if (bool_type)
		output<<" type=\"UInt8\" format=\"ascii\">\n";
	else if (char_type)
		output<<" type=\"Int8\" format=\"ascii\">\n";
	else if (unsigned_char_type)
		output<<" type=\"UInt8\" format=\"ascii\">\n";
	else if (int_type)
		output<<" type=\"Int32\" format=\"ascii\">\n";
	else if (unsigned_int_type)
		output<<" type=\"UInt32\" format=\"ascii\">\n";
	else if (long_type)
		output<<" type=\"Int32\" format=\"ascii\">\n";
	else if (unsigned_long_type)
		output<<" type=\"UInt32\" format=\"ascii\">\n";
	else if (short_type)
		output<<" type=\"Int16\" format=\"ascii\">\n";
	else if (unsigned_short_type)
		output<<" type=\"UInt16\" format=\"ascii\">\n";
	else if (float_type)
		output<<" type=\"Float32\" format=\"ascii\">\n";
	else if (double_type)
		output<<" type=\"Float64\" format=\"ascii\">\n";
	else if (long_double_type)
		output<<" type=\"Float128\" format=\"ascii\">\n";

	// output grid data
	for (int x=x0[0]; x<x1[0]; x++) {
		int size;
		input.read(reinterpret_cast<char*>(&size),sizeof(size));
		char* buffer = new char[size];
		input.read(reinterpret_cast<char*>(buffer),size);

		if (scalar_type) {
			if (bool_type) {
				MMSP::grid<1,MMSP::scalar<bool> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (char_type) {
				MMSP::grid<1,MMSP::scalar<char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_char_type) {
				MMSP::grid<1,MMSP::scalar<unsigned char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (int_type) {
				MMSP::grid<1,MMSP::scalar<int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_int_type) {
				MMSP::grid<1,MMSP::scalar<unsigned int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (long_type) {
				MMSP::grid<1,MMSP::scalar<long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_long_type) {
				MMSP::grid<1,MMSP::scalar<unsigned long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (short_type) {
				MMSP::grid<1,MMSP::scalar<short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_short_type) {
				MMSP::grid<1,MMSP::scalar<unsigned short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (float_type) {
				MMSP::grid<1,MMSP::scalar<float> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (double_type) {
				MMSP::grid<1,MMSP::scalar<double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (long_double_type) {
				MMSP::grid<1,MMSP::scalar<long double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
		}

		else if (vector_type) {
			if (bool_type) {
				MMSP::grid<1,MMSP::vector<bool> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (char_type) {
				MMSP::grid<1,MMSP::vector<char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (unsigned_char_type) {
				MMSP::grid<1,MMSP::vector<unsigned char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (int_type) {
				MMSP::grid<1,MMSP::vector<int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (unsigned_int_type) {
				MMSP::grid<1,MMSP::vector<unsigned int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (long_type) {
				MMSP::grid<1,MMSP::vector<long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (unsigned_long_type) {
				MMSP::grid<1,MMSP::vector<unsigned long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (short_type) {
				MMSP::grid<1,MMSP::vector<short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (unsigned_short_type) {
				MMSP::grid<1,MMSP::vector<unsigned short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (float_type) {
				MMSP::grid<1,MMSP::vector<float> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (double_type) {
				MMSP::grid<1,MMSP::vector<double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
			else if (long_double_type) {
				MMSP::grid<1,MMSP::vector<long double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					for (int j=0; j<fields; j++)
						output<<GRID[i][j]<<" ";
			}
		}

		else if (sparse_type) {
			if (bool_type) {
				MMSP::grid<1,MMSP::sparse<bool> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					bool sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<bool>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (char_type) {
				MMSP::grid<1,MMSP::sparse<char> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					char sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<char>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (unsigned_char_type) {
				MMSP::grid<1,MMSP::sparse<unsigned char> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					unsigned char sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<unsigned char>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (int_type) {
				MMSP::grid<1,MMSP::sparse<int> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					int sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<int>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (unsigned_int_type) {
				MMSP::grid<1,MMSP::sparse<unsigned int> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					unsigned int sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<unsigned int>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (long_type) {
				MMSP::grid<1,MMSP::sparse<long> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					long sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<long>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (unsigned_long_type) {
				MMSP::grid<1,MMSP::sparse<unsigned long> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					unsigned long sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<unsigned long>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (short_type) {
				MMSP::grid<1,MMSP::sparse<short> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					short sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<short>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (unsigned_short_type) {
				MMSP::grid<1,MMSP::sparse<unsigned short> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					unsigned short sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<unsigned short>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (float_type) {
				MMSP::grid<1,MMSP::sparse<float> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					float sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<float>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (double_type) {
				MMSP::grid<1,MMSP::sparse<double> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					double sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<double>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
			else if (long_double_type) {
				MMSP::grid<1,MMSP::sparse<long double> > GRID(0,0,slice,0);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++) {
					long double sum = 0;
					int nonzero = MMSP::length(GRID[i]);
					for (int j=0; j<nonzero; j++)
						sum += static_cast<long double>(pow(GRID[i].value(j),2));
					output<<sum<<" ";
				}
			}
		}

		else /* built-in data types */ {
			if (bool_type) {
				MMSP::grid<1,MMSP::scalar<bool> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (char_type) {
				MMSP::grid<1,MMSP::scalar<char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_char_type) {
				MMSP::grid<1,MMSP::scalar<unsigned char> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (int_type) {
				MMSP::grid<1,MMSP::scalar<int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_int_type) {
				MMSP::grid<1,MMSP::scalar<unsigned int> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (long_type) {
				MMSP::grid<1,MMSP::scalar<long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_long_type) {
				MMSP::grid<1,MMSP::scalar<unsigned long> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (short_type) {
				MMSP::grid<1,MMSP::scalar<short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (unsigned_short_type) {
				MMSP::grid<1,MMSP::scalar<unsigned short> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (float_type) {
				MMSP::grid<1,MMSP::scalar<float> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (double_type) {
				MMSP::grid<1,MMSP::scalar<double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
			else if (long_double_type) {
				MMSP::grid<1,MMSP::scalar<long double> > GRID(fields,0,slice);
				GRID.from_buffer(buffer);
				for (int i=0; i<slice; i++)
					output<<GRID[i]<<" ";
			}
		}

		delete [] buffer;
	}

	// output closing markup
	output<<"\n";
	output<<"        </DataArray>\n";
	output<<"      </CellData>\n";
	output<<"    </Piece>\n";
	output<<"  </ImageData>\n";
	output<<"</VTKFile>\n";
}

