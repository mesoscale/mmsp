// MMSP.output.cpp
// Convert MMSP grid data to XML VTK image data format
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include <cmath>
#include <fstream>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLImageDataWriter.h>

namespace MMSP {

template<int dim, typename T> void print_scalars(std::string filename, const grid<dim,T>& GRID, const int mode)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	MMSP::Abort(-1);
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

template<int dim, typename T> void print_vectors(std::string filename, const grid<dim,vector<T> >& GRID,
        const int mode, const int field)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	MMSP::Abort(-1);
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

template<int dim, typename T> void print_sparses(std::string filename, const grid<dim,sparse<T> >& GRID,
        const int mode, const int field)
{
	#ifdef MPI_VERSION
	std::cerr << "Error: cannot write VTK in parallel." <<std::endl;
	MMSP::Abort(-1);
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

} // namespace
