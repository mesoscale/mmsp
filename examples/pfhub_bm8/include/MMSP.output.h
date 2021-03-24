// MMSP.output.hpp
// Declaration of MMSP output functions
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include <cstring>
#include <set>

namespace MMSP {
	#ifdef PNG_LIBPNG_VER
	int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename);

	template <int dim, typename T>
	void scalar_field_to_png(const grid<dim,T>& GRID,
							 const int& mode, int sliceaxis, int slicelevel,
							 const double& zoomin, const double& zoomax,
							 const bool coninv, const std::set<double>& levelset, const int& lvlfield,
							 const double& contol, const unsigned int& bufsize, unsigned char* buffer);

	template <int dim, typename T>
	void vector_field_to_png(const grid<dim,vector<T> >& GRID,
							 const int& mode, int sliceaxis, int slicelevel,
							 const double& zoomin, const double& zoomax,
							 const bool coninv, const std::set<double>& levelset, const int& lvlfield,
							 const double& contol, const std::set<int>& fieldset,
							 const unsigned int& bufsize, unsigned char* buffer);

	template <int dim, typename T>
	void sparse_field_to_png(const grid<dim,sparse<T> >& GRID,
							 const int& mode, int sliceaxis, int slicelevel,
							 const double& zoomin, const double& zoomax,
							 const bool coninv, const std::set<double>& leveOBlset, const int& lvlfield,
							 const double& contol, const std::set<int>& fieldset,
							 const unsigned int& bufsize, unsigned char* buffer);
	#endif // PNG

	#ifdef VTK_VERSION
	template<int dim, typename T>
	void scalar_field_to_vtk(std::string filename, const grid<dim,T>& GRID, const int mode);

	template<int dim, typename T>
	void vector_field_to_vtk(std::string filename, const grid<dim,vector<T> >& GRID,
							 const int mode, const int field);

	template<int dim, typename T>
	void sparse_field_to_vtk(std::string filename, const grid<dim,sparse<T> >& GRID,
							 const int mode, const int field);
	#endif // VTK

} // namespace

#include "MMSP.output.cpp"
