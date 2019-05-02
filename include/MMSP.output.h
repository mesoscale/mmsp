// MMSP.output.hpp
// Declaration of MMSP output functions
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include <cstring>

namespace MMSP {

	template<int dim, typename T> void print_scalars(std::string filename, const grid<dim,T>& GRID, const int mode);
	template<int dim, typename T> void print_vectors(std::string filename, const grid<dim,vector<T> >& GRID,
													 const int mode, const int field);
	template<int dim, typename T> void print_sparses(std::string filename, const grid<dim,sparse<T> >& GRID,
													 const int mode, const int field);
} // namespace

#include "MMSP.output.cpp"
