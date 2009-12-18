// PFgrid.hpp
// Class definitions for rectilinear phase field grids
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef PFGRID
#define PFGRID
#include"MMSP.grid.hpp"

namespace MMSP{

class PFgrid1D : public grid<1,vector<double> > {
public:
	// constructors
	PFgrid1D(int FIELDS, int x, int GHOSTS=1)
		: grid<1,vector<double> >(FIELDS,0,x,GHOSTS) {}
	PFgrid1D(int FIELDS, int x0, int x1, int GHOSTS=1)
		: grid<1,vector<double> >(FIELDS,x0,x1,GHOSTS) {}
	PFgrid1D(const PFgrid1D& GRID, int FIELDS)
		: grid<1,vector<double> >(GRID,FIELDS) {}
	PFgrid1D(const char* filename, int GHOSTS=1)
		: grid<1,vector<double> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x) const
	{
		MMSP::vector<int> neighbors;
		const PFgrid1D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int index=0; index<length(grid[x+i]); index++) {
				for (int h=0; h<length(neighbors); h++)
					if (index==neighbors[h]) goto skip;
				resize(neighbors,length(neighbors)+1);
				neighbors[length(neighbors)-1] = index;
				skip: continue;
			}
	return neighbors;
	}
};

MMSP::vector<int> neighbors(const PFgrid1D& GRID, int x) {return GRID.neighbors(x);}

template <> void copy(PFgrid1D& GRID1, const PFgrid1D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(PFgrid1D& GRID1, PFgrid1D& GRID2) {GRID1.swap(GRID2);}


class PFgrid2D : public grid<2,vector<double> > {
public:
	// constructors
	PFgrid2D(int FIELDS, int x, int y, int GHOSTS=1)
		: grid<2,vector<double> >(FIELDS,0,x,0,y,GHOSTS) {}
	PFgrid2D(int FIELDS, int x0, int x1, int y0, int y1, int GHOSTS=1)
		: grid<2,vector<double> >(FIELDS,x0,x1,y0,y1,GHOSTS) {}
	PFgrid2D(const PFgrid2D& GRID, int FIELDS)
		: grid<2,vector<double> >(GRID,FIELDS) {}
	PFgrid2D(const char* filename, int GHOSTS=1)
		: grid<2,vector<double> >(filename,GHOSTS) {} 

	// utility functions
	MMSP::vector<int> neighbors(int x, int y) const
	{
		MMSP::vector<int> neighbors;
		const PFgrid2D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++)
				for (int index=0; index<length(grid[x+i][y+j]); index++) {
					for (int h=0; h<length(neighbors); h++)
						if (index==neighbors[h]) goto skip;
					resize(neighbors,length(neighbors)+1);
					neighbors[length(neighbors)-1] = index;
					skip: continue;
				}
		return neighbors;
	}
};

MMSP::vector<int> neighbors(const PFgrid2D& GRID, int x, int y) {return GRID.neighbors(x,y);}

template <> void copy(PFgrid2D& GRID1, const PFgrid2D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(PFgrid2D& GRID1, PFgrid2D& GRID2) {GRID1.swap(GRID2);}


class PFgrid3D : public grid<3,vector<double> > {
public:
	// constructors
	PFgrid3D(int FIELDS, int x, int y, int z, int GHOSTS=1)
		: grid<3,vector<double> >(FIELDS,0,x,0,y,0,z,GHOSTS) {}
	PFgrid3D(int FIELDS, int x0, int x1, int y0, int y1, int z0, int z1, int GHOSTS=1)
		: grid<3,vector<double> >(FIELDS,x0,x1,y0,y1,z0,z1,GHOSTS) {}
	PFgrid3D(const PFgrid3D& GRID, int FIELDS)
		: grid<3,vector<double> >(GRID,FIELDS) {}
	PFgrid3D(const char* filename, int GHOSTS=1)
		: grid<3,vector<double> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x, int y, int z) const
	{
		MMSP::vector<int> neighbors;
		const PFgrid3D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++)
				for (int k=-1; k<=1; k++)
					for (int index=0; index<length(grid[x+i][y+j][z+k]); index++) {
						for (int h=0; h<length(neighbors); h++)
							if (index==neighbors[h]) goto skip;
						resize(neighbors,length(neighbors)+1);
						neighbors[length(neighbors)-1] = index;
						skip: continue;
					}
		return neighbors;
	}
};

MMSP::vector<int> neighbors(const PFgrid3D& GRID, int x, int y, int z) {return GRID.neighbors(x,y,z);}

template <> void copy(PFgrid3D& GRID1, const PFgrid3D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(PFgrid3D& GRID1, PFgrid3D& GRID2) {GRID1.swap(GRID2);}

} // namespace MMSP

#endif
