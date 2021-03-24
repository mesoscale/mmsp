// MMSP.output.cpp
// Convert MMSP grid data to XML VTK image data format
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include <cmath>
#include <fstream>

namespace MMSP {

#ifdef PNG_LIBPNG_VER
int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write PNG in parallel." <<std::endl;
	Abort(-1);
	#endif

	// using libpng
	// After "A simple libpng example program,"
	// http://zarb.org/~gc/html/libpng.html
	// and the libpng manual, http://www.libpng.org/pub/png

	png_byte color_type = PNG_COLOR_TYPE_GRAY;
	// valid choices: PNG_COLOR_TYPE_GRAY       (bit depths 1,2,4, 8, 16)
	//                PNG_COLOR_TYPE_GRAY_ALPHA (bit depths 8, 16)
	//                PNG_COLOR_TYPE_PALETTE    (bit depths 1,2,4, 8)
	//                PNG_COLOR_TYPE_RGB        (bit_depths 8, 16)
	//                PNG_COLOR_TYPE_RGB_ALPHA  (bit_depths 8, 16)
	//                PNG_COLOR_MASK_PALETTE
	//                PNG_COLOR_MASK_COLOR
	//                PNG_COLOR_MASK_ALPHA

	png_byte bit_depth = 8; // valid choices: 1,2,4, 8, 16
	png_structp png_ptr;
	png_infop info_ptr;

	png_bytepp row_pointers = new png_bytep[h];
	for (int j=0; j<h; j++)
		row_pointers[j] = &imData[j*w];

	// Setup PNG file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		std::cerr<<"Error making image: check permissions on "<<filename<<std::endl;
		return (-1);
	}
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		std::cerr<<"Error making image: png_create_write_struct failed."<<std::endl;
		return (-1);
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to init_io."<<std::endl;
		return (-1);
	}
	png_init_io(png_ptr, fp);

	// Write PNG header
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to write header."<<std::endl;
		return (-1);
	}
	png_set_IHDR(png_ptr, info_ptr, w, h,
	                 bit_depth, color_type, PNG_INTERLACE_NONE,
	                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	// Write image
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to write data."<<std::endl;
		return (-1);
	}
	png_write_image(png_ptr, row_pointers);

	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to finish writing."<<std::endl;
		return (-1);
	}
	png_write_end(png_ptr, NULL);

	// Clean up
	delete [] row_pointers;

	fclose(fp);

	return 0;
}

template <int dim, typename T>
void scalar_field_to_png(const grid<dim,T>& GRID,
						 const int& mode, int sliceaxis, int slicelevel,
						 const double& zoomin, const double& zoomax,
						 const bool coninv, const std::set<double>& levelset, const int& lvlfield,
						 const double& contol, const unsigned int& bufsize, unsigned char* buffer)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write PNG in parallel." <<std::endl;
	Abort(-1);
	#endif

	double min=zoomin;
	double max=zoomax;

	for (int n=0; n<nodes(GRID); n++) {
		double val = GRID(n);
		if (mode==1) // mag
			val = std::fabs(val);
		if (val>max)
			max=val;
		else if (val<min)
			min=val;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		vector<int> x(1,0);
		for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
			double val = GRID(x);
			if (mode==1) //mag
				val = std::fabs(val);
			assert(n<bufsize);
			buffer[n] = 255*((val-min)/(max-min));
			if (mode==4) //contour
				for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<contol)
						buffer[n] = 255-buffer[n];
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		vector<int> x(2,0);
		for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
			for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
				double val = GRID(x);
				if (mode==1) //mag
					val = std::fabs(val);
				assert(n<bufsize);
				buffer[n] = 255*((val-min)/(max-min));
				if (mode==4) //contour
					for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
						if (std::fabs(val-*it)/std::fabs(*it)<contol)
							buffer[n] = 255-buffer[n];
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		vector<int> x(3,0);
		for (x[2] = g0(GRID,2); x[2] < g1(GRID,2); x[2]++)
			for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
				for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double val = GRID(x);
					if (mode==1) //mag
						val = std::fabs(val);
					assert(n<bufsize);
					buffer[n] = 255*((val-min)/(max-min));
					if (levelset.size()>0) //contour
						for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
							if (std::fabs(val-*it)/std::fabs(*it)<contol)
								buffer[n] = coninv ? 255 : 0;
					n++;
				}
	}
}

template <int dim, typename T>
void vector_field_to_png(const grid<dim,vector<T> >& GRID,
						 const int& mode, int sliceaxis, int slicelevel,
						 const double& zoomin, const double& zoomax,
						 const bool coninv, const std::set<double>& levelset, const int& lvlfield,
						 const double& contol, const std::set<int>& fieldset,
						 const unsigned int& bufsize, unsigned char* buffer)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write PNG in parallel." <<std::endl;
	Abort(-1);
	#endif

	double min=zoomin;
	double max=zoomax;
	int included = (mode==2) ? fieldset.size() : fields(GRID)-fieldset.size();

	for (int n=0; n<nodes(GRID); n++) {
		double sum=0.0;
		if (mode==0) {        //  no option specified
			for (int i=0; i<fields(GRID); i++)
				if (included>1)
					sum += pow(GRID(n)[i],2.0);
				else
					sum = GRID(n)[i];
		} else if (mode==1) { //  --mag
			for (int i=0; i<fields(GRID); i++)
				sum += pow(GRID(n)[i],2.0);
		} else if (mode==2) { //  --field
			if (included>1)
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
					sum += pow(GRID(n)[*it],2.0);
			else
				sum = GRID(n)[*fieldset.begin()];
		} else if (mode==3) { //  --exclude
			for (int i=0; i<fields(GRID); i++) {
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
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		vector<int> x(1,0);
		for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode==0) { //         default selection
				for (int i=0; i<fields(GRID); i++) {
					if (included>1)
						sum += pow(GRID(x)[i],2.0);
					else
						sum = GRID(x)[i];
				}
			} else if (mode==1) { //  --mag
				for (int i=0; i<fields(GRID); i++)
					sum += pow(GRID(x)[i],2.0);
			} else if (mode==2) { //  --field
				if (fieldset.size()==1)
					sum = GRID(x)[*fieldset.begin()];
				else
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
						sum += pow(GRID(x)[*it],2.0);
			} else if (mode==3) { //  --exclude
				for (int i=0; i<fields(GRID); i++) {
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end()) {
						if (included>1)
							sum += pow(GRID(x)[i],2.0);
						else
							sum = GRID(x)[i];
					}
				}
			}
			if (mode==1 || included!=1)
				sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
					if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						buffer[n] = coninv ? 255 : 0;
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		vector<int> x(2,0);
		for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
			for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode==0) { //         default selection
					for (int i=0; i<fields(GRID); i++) {
						if (included>1)
							sum += pow(GRID(x)[i],2.0);
						else
							sum = GRID(x)[i];
					}
				} else if (mode==1) { //  --mag
					for (int i=0; i<fields(GRID); i++)
						sum += pow(GRID(x)[i],2.0);
				} else if (mode==2) { //  --field
					if (fieldset.size()==1) {
						sum = GRID(x)[*fieldset.begin()];
					} else {
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
							sum += pow(GRID(x)[*it],2.0);
					}
				} else if (mode==3) { //  --exclude
					for (int i=0; i<fields(GRID); i++) {
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end()) {
							if (included>1)
								sum += pow(GRID(x)[i],2.0);
							else
								sum = GRID(x)[i];
						}
					}
				}
				if (mode==1 || included!=1)
					sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (levelset.size()>0) // --contour
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
							buffer[n] = coninv ? 255 : 0;
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		vector<int> x(3,0);
		for (x[2] = g0(GRID,2); x[2] < g1(GRID,2); x[2]++)
			for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
				for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode==0) { //         default selection
						for (int i=0; i<fields(GRID); i++) {
							if (included>1)
								sum += pow(GRID(x)[i],2.0);
							else
								sum = GRID(x)[i];
						}
					} else if (mode==1) { //  --mag
						for (int i=0; i<fields(GRID); i++)
							sum += pow(GRID(x)[i],2.0);
					} else if (mode==2) { //  --field
						if (fieldset.size()==1)
							sum = GRID(x)[*fieldset.begin()];
						else
							for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
								sum += pow(GRID(x)[*it],2.0);
					} else if (mode==3) { //  --exclude
						for (int i=0; i<fields(GRID); i++) {
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end()) {
								if (included>1)
									sum += pow(GRID(x)[i],2.0);
								else
									sum = GRID(x)[i];
							}
						}
					}
					if (mode==1 || included!=1)
						sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (levelset.size()>0) // --contour
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
								buffer[n] = coninv ? 255 : 0;
					n++;
				}
	}
}

template <int dim, typename T>
void sparse_field_to_png(const grid<dim,sparse<T> >& GRID,
						 const int& mode, int sliceaxis, int slicelevel,
						 const double& zoomin, const double& zoomax,
						 const bool coninv, const std::set<double>& levelset, const int& lvlfield,
						 const double& contol, const std::set<int>& fieldset,
						 const unsigned int& bufsize, unsigned char* buffer)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write PNG in parallel." <<std::endl;
	Abort(-1);
	#endif

	double min=zoomin;
	double max=zoomax;
	int included = fieldset.size();

	for (int n=0; n<nodes(GRID); n++) {
		double sum=0.0;
		if (mode<2) { //         --mag or default selection
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
			for (int h=0; h<GRID(n).length(); h++) {
				int i = GRID(n).index(h);
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					sum += pow(GRID(n).value(h),2.0);
			}
		}
		if (included!=1)
			sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		vector<int> x(1,0);
		for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode<2) { //          --mag
				for (int h=0; h<GRID(x).length(); h++)
					sum += pow(GRID(x).value(h),2.0);
			} else if (mode==2) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
					if (included>1)
						sum += pow(GRID(x)[*it],2.0);
					else
						sum = GRID(x)[*it];
				}
			} else if (mode==3) { //  --exclude
				for (int h=0; h<GRID(x).length(); h++) {
					int i = GRID(x).index(h);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(x).value(h),2.0);
				}
			}
			if (mode!=2 || included!=1)
				sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (levelset.size()>0) // --contour
				for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
					if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
						buffer[n] = coninv ? 255 : 0;
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		vector<int> x(2,0);
		for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
			for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode<2) { //          --mag
					for (int h=0; h<GRID(x).length(); h++)
						sum += pow(GRID(x).value(h),2.0);
				} else if (mode==2) { //  --field
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
						if (included>1)
							sum += pow(GRID(x)[*it],2.0);
						else
							sum = GRID(x)[*it];
					}
				} else if (mode==3) { //  --exclude
					for (int h=0; h<GRID(x).length(); h++) {
						int i = GRID(x).index(h);
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end())
							sum += pow(GRID(x).value(h),2.0);
					}
				}
				if (mode!=2 || included!=1)
					sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (levelset.size()>0) // --contour
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
							buffer[n] = coninv ? 255 : 0;
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		vector<int> x(3,0);
		for (x[2] = g0(GRID,2); x[2] < g1(GRID,2); x[2]++)
			for (x[1] = g1(GRID,1)-1; x[1] >= g0(GRID,1); x[1]--)
				for (x[0] = g0(GRID,0); x[0] < g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode<2) { //          --mag
						for (int h=0; h<GRID(x).length(); h++)
							sum += pow(GRID(x).value(h),2.0);
					} else if (mode==2) { //  --field
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
							if (included>1)
								sum += pow(GRID(x)[*it],2.0);
							else
								sum = GRID(x)[*it];
						}
					} else if (mode==3) { //  --exclude
						for (int h=0; h<GRID(x).length(); h++) {
							int i = GRID(x).index(h);
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end())
								sum += pow(GRID(x).value(h),2.0);
						}
					}
					if (mode!=2 || included!=1)
						sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (levelset.size()>0) // --contour
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(GRID(x)[lvlfield]-*itl)/std::fabs(*itl)<contol)
								buffer[n] = coninv ? 255 : 0;
					n++;
				}
	}
}
#endif // PNG

#ifdef VTK_VERSION
template<int dim, typename T>
void scalar_field_to_vtk(std::string filename, const grid<dim,T>& GRID, const int mode)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	Abort(-1);
	#endif

	vtkSmartPointer<vtkImageData> scalarData = vtkSmartPointer<vtkImageData>::New();

	if (dim==1) {
		scalarData->SetDimensions(x1(GRID)-x0(GRID), 1, 1);
		scalarData->SetExtent(x0(GRID), x1(GRID) - 1, 0, 0, 0, 0);
		scalarData->SetSpacing(dx(GRID), 1, 1);
		#if VTK_MAJOR_VERSION <= 5
		scalarData->SetScalarTypeToFloat();
		scalarData->SetNumberOfScalarComponents(1);
		#else
		scalarData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(1,0);
		for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
			float* pixel = static_cast<float*>(scalarData->GetScalarPointer(x[0], 0, 0));
			if (mode==1) { // --mag
				*pixel = std::sqrt(GRID(x)*GRID(x));
			} else {
				*pixel = GRID(x);
			}
		}
	} else if (dim==2) {
		scalarData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       1);
		scalarData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   0, 0);
		scalarData->SetSpacing(dx(GRID, 0), dx(GRID, 1), 1);
		#if VTK_MAJOR_VERSION <= 5
		scalarData->SetScalarTypeToFloat();
		scalarData->SetNumberOfScalarComponents(1);
		#else
		scalarData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(2,0);
		for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
			for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
				float* pixel = static_cast<float*>(scalarData->GetScalarPointer(x[0], x[1], 0));
				if (mode==1) { // --mag
					*pixel = std::sqrt(GRID(x)*GRID(x));
				} else {
					*pixel = GRID(x);
				}
			}
		}
	} else if (dim==3) {
		scalarData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       z1(GRID)-z0(GRID));
		scalarData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   z0(GRID), z1(GRID) - 1);
		scalarData->SetSpacing(dx(GRID, 0), dx(GRID, 1), dx(GRID, 2));
		#if VTK_MAJOR_VERSION <= 5
		scalarData->SetScalarTypeToFloat();
		scalarData->SetNumberOfScalarComponents(1);
		#else
		scalarData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(3,0);
		for (x[2]=z0(GRID); x[2]<z1(GRID); x[2]++) {
			for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
				for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
					float* pixel = static_cast<float*>(scalarData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==1) { // --mag
						*pixel = std::sqrt(GRID(x)*GRID(x));
					} else {
						*pixel = GRID(x);
					}
				}
			}
		}
	}
	scalarData->GetPointData()->GetAbstractArray(0)->SetName("scalar_data");

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(scalarData->GetProducerPort());
	#else
	writer->SetInputData(scalarData);
	#endif
	writer->SetDataModeToBinary();
	writer->SetCompressorTypeToZLib();
	writer->Write();

	scalarData = NULL;
	writer = NULL;
}

template<int dim, typename T>
void vector_field_to_vtk(std::string filename, const grid<dim,vector<T> >& GRID,
						 const int mode, const int field)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	Abort(-1);
	#endif

	vtkSmartPointer<vtkImageData> vectorData = vtkSmartPointer<vtkImageData>::New();

	if (dim==1) {
		vectorData->SetDimensions(x1(GRID)-x0(GRID), 1, 1);
		vectorData->SetExtent(x0(GRID), x1(GRID) - 1, 0, 0, 0, 0);
		vectorData->SetSpacing(dx(GRID), 1, 1);
		#if VTK_MAJOR_VERSION <= 5
		vectorData->SetScalarTypeToFloat();
		if (mode==1 || mode==2 || mode==3)
			vectorData->SetNumberOfScalarComponents(1);
		else
			vectorData->SetNumberOfScalarComponents(fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vectorData->AllocateScalars(VTK_FLOAT, 1);
		else
			vectorData->AllocateScalars(VTK_FLOAT, fields(GRID));
		#endif

		vector<int> x(1,0);
		for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
			const vector<T>& v = GRID(x);
			float* pixel = static_cast<float*>(vectorData->GetScalarPointer(x[0], 0, 0));
			if (mode==1) { // --mag
				double sum = 0.0;
				for (int h = 0; h < v.length(); h++)
					sum += v[h]*v[h];
				*pixel = std::sqrt(sum);
			} else if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < v.length(); h++)
					if (v[h] > v[max])
						max = h;
				*pixel = max;
			} else if (mode==3) { // --field
				*pixel = v[field];
			} else {
				for (int h = 0; h < v.length(); h++)
					pixel[h] = v[h];
			}
		}
	} else if (dim==2) {
		vectorData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       1);
		vectorData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   0, 0);
		vectorData->SetSpacing(dx(GRID, 0), dx(GRID, 1), 1);
		#if VTK_MAJOR_VERSION <= 5
		vectorData->SetScalarTypeToFloat();
		if (mode==1 || mode==2 || mode==3)
			vectorData->SetNumberOfScalarComponents(1);
		else
			vectorData->SetNumberOfScalarComponents(fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vectorData->AllocateScalars(VTK_DOUBLE, 1);
		else
			vectorData->AllocateScalars(VTK_DOUBLE, fields(GRID));
		#endif

		vector<int> x(2,0);
		for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
			for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
				const vector<T>& v = GRID(x);
				double* pixel = static_cast<double*>(vectorData->GetScalarPointer(x[0], x[1], 0));
				if (mode==1) { // --mag
					double sum = 0.0;
					for (int h = 0; h < v.length(); h++)
						sum += v[h]*v[h];
					*pixel = std::sqrt(sum);
				} else if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < v.length(); h++)
						if (v[h] > v[max])
							max = h;
					*pixel = max;
				} else if (mode==3) { // --field
					*pixel = v[field];
				} else {
					for (int h = 0; h < v.length(); h++)
						pixel[h] = v[h];
				}
			}
		}
	} else if (dim==3) {
		vectorData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       z1(GRID)-z0(GRID));
		vectorData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   z0(GRID), z1(GRID) - 1);
		vectorData->SetSpacing(dx(GRID, 0), dx(GRID, 1), dx(GRID, 2));
		#if VTK_MAJOR_VERSION <= 5
		vectorData->SetScalarTypeToFloat();
		if (mode==1 || mode==2 || mode==3)
			vectorData->SetNumberOfScalarComponents(1);
		else
			vectorData->SetNumberOfScalarComponents(fields(GRID));
		#else
		if (mode==1 || mode==2 || mode==3)
			vectorData->AllocateScalars(VTK_FLOAT, 1);
		else
			vectorData->AllocateScalars(VTK_FLOAT, fields(GRID));
		#endif

		vector<int> x(3,0);
		for (x[2]=z0(GRID); x[2]<z1(GRID); x[2]++) {
			for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
				for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
					const vector<T>& v = GRID(x);
					float* pixel = static_cast<float*>(vectorData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==1) { // --mag
						double sum = 0.0;
						for (int h = 0; h < v.length(); h++)
							sum += v[h]*v[h];
						*pixel = std::sqrt(sum);
					} else if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < v.length(); h++)
							if (v[h] > v[max])
								max = h;
						*pixel = max;
					} else if (mode==3) { // --field
						*pixel = v[field];
					} else {
						for (int h = 0; h < v.length(); h++)
							pixel[h] = v[h];
					}
				}
			}
		}
	}

	// vectorData->GetPointData()->GetAbstractArray(0)->SetName("vector_data");

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(vectorData->GetProducerPort());
	#else
	writer->SetInputData(vectorData);
	#endif
	writer->SetDataModeToBinary();
	writer->SetCompressorTypeToZLib();
	writer->Write();

	vectorData = NULL;
	writer = NULL;
}

template<int dim, typename T>
void sparse_field_to_vtk(std::string filename, const grid<dim,sparse<T> >& GRID,
						 const int mode, const int field)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	Abort(-1);
	#endif

	vtkSmartPointer<vtkImageData> sparseData = vtkSmartPointer<vtkImageData>::New();

	if (dim==1) {
		sparseData->SetDimensions(x1(GRID)-x0(GRID), 1, 1);
		sparseData->SetExtent(x0(GRID), x1(GRID) - 1, 0, 0, 0, 0);
		sparseData->SetSpacing(dx(GRID, 0), 1, 1);
		#if VTK_MAJOR_VERSION <= 5
		sparseData->SetNumberOfScalarComponents(1);
		sparseData->SetScalarTypeToFloat();
		#else
		sparseData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(1,0);
		for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
			const sparse<T>& s = GRID(x);
			float* pixel = static_cast<float*>(sparseData->GetScalarPointer(x[0], 0, 0));
			if (mode==2) { // --max
				// Export index of field with greatest magnitude
				int max = 0;
				for (int h = 1; h < s.length(); h++)
					if (s.value(h) > s.value(max))
						max = h;
				*pixel =  s.index(max);
			} else if (mode==3) { // --field
				*pixel = s[field];
			} else { // --mag is redundant for sparse
				double sum = 0.0;
				for (int h = 0; h < s.length(); h++)
					sum += s.value(h)*s.value(h);
				*pixel = std::sqrt(sum);
			}
		}
	} else if (dim==2) {
		sparseData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       1);
		sparseData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   0, 0);
		sparseData->SetSpacing(dx(GRID, 0), dx(GRID, 1), 1);
		#if VTK_MAJOR_VERSION <= 5
		sparseData->SetNumberOfScalarComponents(1);
		sparseData->SetScalarTypeToFloat();
		#else
		sparseData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(2,0);
		for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
			for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
				const sparse<T>& s = GRID(x);
				float* pixel = static_cast<float*>(sparseData->GetScalarPointer(x[0], x[1], 0));
				if (mode==2) { // --max
					// Export index of field with greatest magnitude
					int max = 0;
					for (int h = 1; h < s.length(); h++)
						if (s.value(h) > s.value(max))
							max = h;
					*pixel =  s.index(max);
				} else if (mode==3) { // --field
					*pixel = s[field];
				} else { // --mag is redundant for sparse
					double sum = 0.0;
					for (int h = 0; h < s.length(); h++)
						sum += s.value(h)*s.value(h);
					*pixel = std::sqrt(sum);
				}
			}
		}
	} else if (dim==3) {
		sparseData->SetDimensions(x1(GRID)-x0(GRID),
		                       y1(GRID)-y0(GRID),
		                       z1(GRID)-z0(GRID));
		sparseData->SetExtent(x0(GRID), x1(GRID) - 1,
		                   y0(GRID), y1(GRID) - 1,
		                   z0(GRID), z1(GRID) - 1);
		sparseData->SetSpacing(dx(GRID, 0), dx(GRID, 1), dx(GRID, 2));
		#if VTK_MAJOR_VERSION <= 5
		sparseData->SetNumberOfScalarComponents(1);
		sparseData->SetScalarTypeToFloat();
		#else
		sparseData->AllocateScalars(VTK_FLOAT, 1);
		#endif

		vector<int> x(3,0);
		for (x[2]=z0(GRID); x[2]<z1(GRID); x[2]++) {
			for (x[1]=y0(GRID); x[1]<y1(GRID); x[1]++) {
				for (x[0]=x0(GRID); x[0]<x1(GRID); x[0]++) {
					const sparse<T>& s = GRID(x);
					float* pixel = static_cast<float*>(sparseData->GetScalarPointer(x[0], x[1], x[2]));
					if (mode==2) { // --max
						// Export index of field with greatest magnitude
						int max = 0;
						for (int h = 1; h < s.length(); h++)
							if (s.value(h) > s.value(max))
								max = h;
						*pixel =  s.index(max);
					} else if (mode==3) { // --field
						*pixel = s[field];
					} else { // --mag is redundant for sparse
						double sum = 0.0;
						for (int h = 0; h < s.length(); h++)
							sum += s.value(h)*s.value(h);
						*pixel = std::sqrt(sum);
					}
				}
			}
		}
	}
	sparseData->GetPointData()->GetAbstractArray(0)->SetName("sparse_data");

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(sparseData->GetProducerPort());
	#else
	writer->SetInputData(sparseData);
	#endif
	writer->SetDataModeToBinary();
	writer->SetCompressorTypeToZLib();
	writer->Write();

	sparseData = NULL;
	writer = NULL;
}
#endif // VTK

} // namespace
