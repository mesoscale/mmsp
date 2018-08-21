// mmsp2vti.cpp
// Convert MMSP grid data to XML VTK image data format
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include <cmath>
#include <fstream>
#include <string>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLWriter.h>
#include "MMSP.hpp"

template<int dim, typename T> void print_scalars(std::string filename, const MMSP::grid<dim,T>& GRID, const int mode)
{
	vtkSmartPointer<vtkImageData> vtkData = vtkSmartPointer<vtkImageData>::New();
	vtkData->SetDimensions(MMSP::xlength(GRID, 0), MMSP::xlength(GRID, 1), MMSP::xlength(GRID, 2));
	#if VTK_MAJOR_VERSION <= 5
	vtkData->SetNumberOfScalarComponents(1);
	vtkData->SetScalarTypeToDouble();
	#else
	vtkData->AllocateScalars(VTK_DOUBLE, 1);
	#endif

	if (dim==1) {
		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], 0, 0));
			if (mode==1) { // --mag
				pixel[0] = std::sqrt(GRID(x)*GRID(x));
			} else {
				pixel[0] = GRID(x);
			}
		}
	} else if (dim==2) {
		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], 0));
				if (mode==1) { // --mag
					pixel[0] = std::sqrt(GRID(x)*GRID(x));
				} else {
					pixel[0] = GRID(x);
				}
			}
		}
	} else if (dim==3) {
		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==1) { // --mag
						pixel[0] = std::sqrt(GRID(x)*GRID(x));
					} else {
						pixel[0] = GRID(x);
					}
				}
			}
		}
	}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetCompressorTypeToZLib();
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(vtkData->GetProducerPort());
	#else
	writer->SetInputData(vtkData);
	#endif
	writer->Write();

	vtkData = NULL;
	writer = NULL;
}

template<int dim, typename T> void print_vectors(std::string filename, const MMSP::grid<dim,MMSP::vector<T> >& GRID,
        const int mode, const int field)
{
	vtkSmartPointer<vtkImageData> vtkData = vtkSmartPointer<vtkImageData>::New();
	#if VTK_MAJOR_VERSION <= 5
	vtkData->SetScalarTypeToDouble();
	#endif

	if (dim==1) {
		vtkData->SetDimensions(MMSP::x1(GRID)-MMSP::x0(GRID), 1, 1);
	    vtkData->SetExtent(MMSP::x0(GRID), MMSP::x1(GRID) - 1, 0, 0, 0, 0);
		vtkData->SetSpacing(MMSP::dx(GRID), 1, 1);
		#if VTK_MAJOR_VERSION <= 5
		if (mode==1 || mode==2 || mode==3)
			vtkData->SetNumberOfScalarComponents(1);
		else
			vtkData->SetNumberOfScalarComponents(MMSP::fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vtkData->AllocateScalars(VTK_DOUBLE, 1);
		else
			vtkData->AllocateScalars(VTK_DOUBLE, MMSP::fields(GRID));
		#endif

		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			const MMSP::vector<T>& v = GRID(x);
			double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], 0, 0));
			if (mode==1) { // --mag
				double sum = 0.0;
				for (int h = 0; h < v.length(); h++)
					sum += v[h]*v[h];
				pixel[0] = std::sqrt(sum);
			} else if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < v.length(); h++)
					if (v[h] > v[max])
						max = h;
				pixel[0] = max;
			} else if (mode==3) { // --field
				pixel[0] = v[field];
			} else {
				for (int h = 0; h < v.length(); h++)
					pixel[h] = v[h];
			}
		}
	} else if (dim==2) {
		vtkData->SetDimensions(MMSP::x1(GRID)-MMSP::x0(GRID), MMSP::y1(GRID)-MMSP::y0(GRID), 1);
		vtkData->SetExtent(MMSP::x0(GRID), MMSP::x1(GRID) - 1, MMSP::y0(GRID), MMSP::y1(GRID) - 1, 0, 0);
		vtkData->SetSpacing(MMSP::dx(GRID, 0), MMSP::dx(GRID, 1), 1);
		#if VTK_MAJOR_VERSION <= 5
		if (mode==1 || mode==2 || mode==3)
			vtkData->SetNumberOfScalarComponents(1);
		else
			vtkData->SetNumberOfScalarComponents(MMSP::fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vtkData->AllocateScalars(VTK_DOUBLE, 1);
		else
			vtkData->AllocateScalars(VTK_DOUBLE, MMSP::fields(GRID));
		#endif

		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				const MMSP::vector<T>& v = GRID(x);
				double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], 0));
				if (mode==1) { // --mag
					double sum = 0.0;
					for (int h = 0; h < v.length(); h++)
						sum += v[h]*v[h];
					pixel[0] = std::sqrt(sum);
				} else if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < v.length(); h++)
						if (v[h] > v[max])
							max = h;
					pixel[0] = max;
				} else if (mode==3) { // --field
					pixel[0] = v[field];
				} else {
					for (int h = 0; h < v.length(); h++)
						pixel[h] = v[h];
				}
			}
		}
	} else if (dim==3) {
		vtkData->SetDimensions(MMSP::x1(GRID)-MMSP::x0(GRID), MMSP::y1(GRID)-MMSP::y0(GRID), MMSP::z1(GRID)-MMSP::z0(GRID));
		vtkData->SetExtent(MMSP::x0(GRID), MMSP::x1(GRID) - 1, MMSP::y0(GRID), MMSP::y1(GRID) - 1, MMSP::z0(GRID), MMSP::z1(GRID) - 1);
		vtkData->SetSpacing(MMSP::dx(GRID, 0), MMSP::dx(GRID, 1), MMSP::dx(GRID, 2));
		#if VTK_MAJOR_VERSION <= 5
		if (mode==1 || mode==2 || mode==3)
			vtkData->SetNumberOfScalarComponents(1);
		else
			vtkData->SetNumberOfScalarComponents(MMSP::fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vtkData->AllocateScalars(VTK_DOUBLE, 1);
		else
			vtkData->AllocateScalars(VTK_DOUBLE, MMSP::fields(GRID));
		#endif

		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					const MMSP::vector<T>& v = GRID(x);
					double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==1) { // --mag
						double sum = 0.0;
						for (int h = 0; h < v.length(); h++)
							sum += v[h]*v[h];
						pixel[0] = std::sqrt(sum);
					} else if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < v.length(); h++)
							if (v[h] > v[max])
								max = h;
						pixel[0] = max;
					} else if (mode==3) { // --field
						pixel[0] = v[field];
					} else {
						for (int h = 0; h < v.length(); h++)
							pixel[h] = v[h];
					}
				}
			}
		}
	}
	vtkData->GetPointData()->GetAbstractArray(0)->SetName("vector_data");

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetCompressorTypeToZLib();
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(vtkData->GetProducerPort());
	#else
	writer->SetInputData(vtkData);
	#endif
	writer->Write();

	vtkData = NULL;
	writer = NULL;
}

template<int dim, typename T> void print_sparses(std::string filename, const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
        const int mode, const int field)
{
	vtkSmartPointer<vtkImageData> vtkData = vtkSmartPointer<vtkImageData>::New();
	vtkData->SetDimensions(MMSP::xlength(GRID, 0), MMSP::xlength(GRID, 1), MMSP::xlength(GRID, 2));
	#if VTK_MAJOR_VERSION <= 5
	vtkData->SetNumberOfScalarComponents(1);
	vtkData->SetScalarTypeToDouble();
	#else
	vtkData->AllocateScalars(VTK_DOUBLE, 1);
	#endif

	if (dim==1) {
		MMSP::vector<int> x(1,0);
		for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
			const MMSP::sparse<T>& s = GRID(x);
			double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], 0, 0));
			if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < s.length(); h++)
					if (s.value(h) > s.value(max))
						max = h;
				pixel[0] =  s.index(max);
			} else if (mode==3) { // --field
				pixel[0] = s[field];
			} else { // --mag is redundant for sparse
				double sum = 0.0;
				for (int h = 0; h < s.length(); h++)
					sum += s.value(h)*s.value(h);
				pixel[0] = std::sqrt(sum);
			}
		}
	} else if (dim==2) {
		MMSP::vector<int> x(2,0);
		for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
			for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
				const MMSP::sparse<T>& s = GRID(x);
				double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], 0));
				if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < s.length(); h++)
						if (s.value(h) > s.value(max))
							max = h;
					pixel[0] =  s.index(max);
				} else if (mode==3) { // --field
					pixel[0] = s[field];
				} else { // --mag is redundant for sparse
					double sum = 0.0;
					for (int h = 0; h < s.length(); h++)
						sum += s.value(h)*s.value(h);
					pixel[0] = std::sqrt(sum);
				}
			}
		}
	} else if (dim==3) {
		MMSP::vector<int> x(3,0);
		for (x[2]=MMSP::z0(GRID); x[2]<MMSP::z1(GRID); x[2]++) {
			for (x[1]=MMSP::y0(GRID); x[1]<MMSP::y1(GRID); x[1]++) {
				for (x[0]=MMSP::x0(GRID); x[0]<MMSP::x1(GRID); x[0]++) {
					const MMSP::sparse<T>& s = GRID(x);
					double* pixel = static_cast<double*>(vtkData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < s.length(); h++)
							if (s.value(h) > s.value(max))
								max = h;
						pixel[0] =  s.index(max);
					} else if (mode==3) { // --field
						pixel[0] = s[field];
					} else { // --mag is redundant for sparse
						double sum = 0.0;
						for (int h = 0; h < s.length(); h++)
							sum += s.value(h)*s.value(h);
						pixel[0] = std::sqrt(sum);
					}
				}
			}
		}
	}

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetCompressorTypeToZLib();
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(vtkData->GetProducerPort());
	#else
	writer->SetInputData(vtkData);
	#endif
	writer->Write();

	vtkData = NULL;
	writer = NULL;
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

	int fileindex = 1; // in typical usage, filename comes immediately after executable
	int flatten = 0;
	int field = 0;

	if (std::string(argv[1]) == "--mag") {
		flatten=1;
		fileindex=2;
	} else if (std::string(argv[1]) == "--max") {
		flatten=2;
		fileindex=2;
	} else if (std::string(argv[1]).substr(0,8) == "--field=") {
		flatten=3;
		fileindex=2;
		std::string str(argv[1]);
		field = atoi(str.substr(8,12).c_str());
	}

	// file open error check
	std::ifstream input(argv[fileindex]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[fileindex] << ".\n";
		exit(-1);
	}

	// generate output file name
	std::stringstream filename;
	if (argc < 3 || (flatten>0 && argc<4))
		filename << std::string(argv[fileindex]).substr(0, std::string(argv[fileindex]).find_last_of(".")) << ".vti";
	else
		filename << argv[fileindex+1];

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

	// write grid data
	if (not vector_type and not sparse_type) { // must be scalar or built-in
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<bool> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<bool> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<bool> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<unsigned char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<unsigned char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<unsigned char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<char> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<unsigned int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<unsigned int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<unsigned int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<int> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<unsigned long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<unsigned long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<unsigned long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<long> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<unsigned short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<unsigned short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<unsigned short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<short> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<float> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<float> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<float> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<long double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<long double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<long double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		} else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::scalar<double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::scalar<double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::scalar<double> > GRID(argv[fileindex]);
				print_scalars(filename.str(), GRID, flatten);
			}
		}
	}

	else if (vector_type) {
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<bool> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<bool> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<bool> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<unsigned char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<unsigned char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<unsigned char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<char> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<unsigned int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<unsigned int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<unsigned int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<int> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<unsigned long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<unsigned long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<unsigned long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<long> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<unsigned short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<unsigned short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<unsigned short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<short> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<float> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<float> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<float> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<long double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<long double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<long double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		} else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::vector<double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::vector<double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::vector<double> > GRID(argv[fileindex]);
				print_vectors(filename.str(), GRID, flatten, field);
			}
		}
	}

	else if (sparse_type) {
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<bool> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<bool> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<bool> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<unsigned char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<unsigned char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<unsigned char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<char> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<unsigned int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<unsigned int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<unsigned int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<int> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<unsigned long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<unsigned long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<unsigned long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<long> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<unsigned short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<unsigned short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<unsigned short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<short> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<float> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<float> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<float> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<long double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<long double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<long double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		} else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1, MMSP::sparse<double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 2) {
				MMSP::grid<2, MMSP::sparse<double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			} else if (dim == 3) {
				MMSP::grid<3, MMSP::sparse<double> > GRID(argv[fileindex]);
				print_sparses(filename.str(), GRID, flatten, field);
			}
		}
	}

}
