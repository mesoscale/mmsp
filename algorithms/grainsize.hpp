// grainsize.hpp
// Grain size algorithm for MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINSIZE
#define GRAINSIZE
#include"MMSP.grid.hpp"

namespace MMSP{

MMSP::sparse<int> grainsize(const MMSP::grid<int>& grid)
{
	// Compute the size of each domain within a grid and
	// return a sparse<int> structure containing each size.

	// For integer-valued grids, this amounts to counting
	// the number of occurances of each value.  We allow
	// the indices to take any integer value, possibly not
	// consecutive, so we use the sparse data structure 
	// which in most cases should save both time and memory.  

	// Note that there are no checks to determine whether
	// the domains are simply connected; users who wish
	// to compute the sizes of connected components should
	// apply a burn algorithm to the grid before using
	// this function.  

	sparse<int> grainsize;

	// accumulate local domain sizes
	for (int i=0; i<nodes(grid); i++)
		MMSP::set(grainsize,grid(i)) += 1;

	// sum up domain sizes over all partitions
	grainsize = MMSP::global(grainsize,"sum");

	return grainsize;
}

MMSP::sparse<int> grainsize(const MMSP::grid<MMSP::vector<double> >& grid)
{
	// Compute the size of each domain within a grid and
	// return a sparse<int> structure containing each size.

	// For vector-valued grids, we choose that index which
	// has the largest corresponding value to represent the
	// domain id.  Here, we use the sparse data structure 
	// primarily to retain consistency between all instances
	// of "grainsize." 

	// Note that there are no checks to determine whether
	// the domains are simply connected; users who wish
	// to compute the sizes of connected components should
	// apply a burn algorithm to the grid before using
	// this function.  

	sparse<int> grainsize;

	// accumulate local domain sizes
	for (int i=0; i<nodes(grid); i++) {
		int index = 0;
		double value = grid(i)[0];
		for (int j=1; j<MMSP::length(grid(i)); j++)
			if (grid(i)[j]>value) {
				index = j;
				value = grid(i)[j];
			}
		MMSP::set(grainsize,index) += 1;
	}

	// sum up domain sizes over all partitions
	grainsize = MMSP::global(grainsize,"sum");

	return grainsize;
}

MMSP::sparse<int> grainsize(const MMSP::grid<MMSP::sparse<double> >& grid)
{
	// Compute the size of each domain within a grid and
	// return a sparse<int> structure containing each size.

	// For grids containing sparse data, we choose the index
	// with the largest corresponding value to represent the
	// domain id.  We allow the indices to take any integer
	// value, possibly not consecutive, so we use the sparse
	// data structure which in most cases should save both
	// time and memory. 

	// Note that there are no checks to determine whether
	// the domains are simply connected; users who wish
	// to compute the sizes of connected components should
	// apply a burn algorithm to the grid before using
	// this function.  

	sparse<int> grainsize;

	// accumulate local domain sizes
	for (int i=0; i<nodes(grid); i++) {
		if (MMSP::length(grid(i))>0) {
			int index = 0;
			double value = 0.0;
			for (int j=0; j<MMSP::length(grid(i)); j++)
				if (MMSP::value(grid(i),j)>value) {
					index = MMSP::index(grid(i),j);
					value = MMSP::value(grid(i),j);
				}
			MMSP::set(grainsize,index) += 1;
		}
	}

	// sum up domain sizes over all partitions
	grainsize = MMSP::global(grainsize,"sum");

	return grainsize;
}

} // namespace MMSP

#endif
