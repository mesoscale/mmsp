// MCgrid.hpp
// Class definitions for rectilinear Monte Carlo grids
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MCGRID
#define MCGRID
#include"MMSP.grid.hpp"

namespace MMSP{

class MCgrid1D : public grid<1,scalar<int> > {
public:
	// constructors
	MCgrid1D(int x, int GHOSTS=1)
		: grid<1,scalar<int> >(1,0,x,GHOSTS) {}
	MCgrid1D(int x0, int x1, int GHOSTS=1)
		: grid<1,scalar<int> >(1,x0,x1,GHOSTS) {}
	MCgrid1D(const MCgrid1D& GRID, int FIELDS)
		: grid<1,scalar<int> >(GRID,FIELDS) {}
	MCgrid1D(const char* filename, int GHOSTS=1)
		: grid<1,scalar<int> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x) 
	{
		const MCgrid1D& grid = *this;
		MMSP::vector<int> neighbors;
		for (int i=-1; i<=1; i++) {
			int index = grid[x+i];
			for (int h=0; h<length(neighbors); h++)
				if (index==neighbors[h]) goto skip;
			resize(neighbors,length(neighbors)+1);
			neighbors[length(neighbors)-1] = index;
			skip: continue;
		}
		return neighbors;
	}

	friend MMSP::vector<int> neighbors(MCgrid1D& GRID, int x)
		{return GRID.neighbors(x);}
};

template <> void copy(MCgrid1D& GRID1, const MCgrid1D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(MCgrid1D& GRID1, MCgrid1D& GRID2) {GRID1.swap(GRID2);}


class MCgrid2D : public grid<2,scalar<int> > {
public:
	// constructors
	MCgrid2D(int x, int y, int GHOSTS=1)
		: grid<2,scalar<int> >(1,0,x,0,y,GHOSTS) {}
	MCgrid2D(int x0, int x1, int y0, int y1, int GHOSTS=1)
		: grid<2,scalar<int> >(1,x0,x1,y0,y1,GHOSTS) {}
	MCgrid2D(const MCgrid2D& GRID, int FIELDS)
		: grid<2,scalar<int> >(GRID,FIELDS) {}
	MCgrid2D(const char* filename, int GHOSTS=1)
		: grid<2,scalar<int> >(filename,GHOSTS) {} 

	// utility functions
	MMSP::vector<int> neighbors(int x, int y) 
	{
		const MCgrid2D& grid = *this;
		MMSP::vector<int> neighbors;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++) {
				int index = grid[x+i][y+j];
				for (int h=0; h<length(neighbors); h++)
					if (index==neighbors[h]) goto skip;
				resize(neighbors,length(neighbors)+1);
				neighbors[length(neighbors)-1] = index;
				skip: continue;
			}
		return neighbors;
	}

	friend MMSP::vector<int> neighbors(MCgrid2D& GRID, int x, int y)
		{return GRID.neighbors(x,y);}
};

template <> void copy(MCgrid2D& GRID1, const MCgrid2D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(MCgrid2D& GRID1, MCgrid2D& GRID2) {GRID1.swap(GRID2);}


class MCgrid3D : public grid<3,scalar<int> > {
public:
	// constructors
	MCgrid3D(int x, int y, int z, int GHOSTS=1)
		: grid<3,scalar<int> >(1,0,x,0,y,0,z,GHOSTS) {}
	MCgrid3D(int x0, int x1, int y0, int y1, int z0, int z1, int GHOSTS=1)
		: grid<3,scalar<int> >(1,x0,x1,y0,y1,z0,z1,GHOSTS) {}
	MCgrid3D(const MCgrid3D& GRID, int FIELDS)
		: grid<3,scalar<int> >(GRID,FIELDS) {}
	MCgrid3D(const char* filename, int GHOSTS=1)
		: grid<3,scalar<int> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x, int y, int z) 
	{
		const MCgrid3D& grid = *this;
		MMSP::vector<int> neighbors;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++)
				for (int k=-1; k<=1; k++) {
					int index = grid[x+i][y+j][z+k];
					for (int h=0; h<length(neighbors); h++)
						if (index==neighbors[h]) goto skip;
					resize(neighbors,length(neighbors)+1);
					neighbors[length(neighbors)-1] = index;
					skip: continue;
				}
		return neighbors;
	}

	friend MMSP::vector<int> neighbors(MCgrid3D& GRID, int x, int y, int z)
		{return GRID.neighbors(x,y,z);}
};

template <> void copy(MCgrid3D& GRID1, const MCgrid3D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(MCgrid3D& GRID1, MCgrid3D& GRID2) {GRID1.swap(GRID2);}

} // namespace MMSP

#endif
