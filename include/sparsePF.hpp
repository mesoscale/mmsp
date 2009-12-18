// sparsePF.hpp
// Class definitions for rectilinear sparse phase field grids
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef SPARSEPF 
#define SPARSEPF 
#include"MMSP.grid.hpp"

namespace MMSP{

class sparsePF1D : public grid<1,sparse<double> > {
public:
	// constructors
	sparsePF1D(int x, int GHOSTS=1)
		: grid<1,sparse<double> >(0,0,x,GHOSTS) {}
	sparsePF1D(int x0, int x1, int GHOSTS=1)
		: grid<1,sparse<double> >(0,x0,x1,GHOSTS) {}
	sparsePF1D(const sparsePF1D& GRID, int FIELDS)
		: grid<1,sparse<double> >(GRID,FIELDS) {}
	sparsePF1D(const char* filename, int GHOSTS=1)
		: grid<1,sparse<double> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x) const
	{
		MMSP::vector<int> neighbors;
		const sparsePF1D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int n=0; n<length(grid[x+i]); n++) {
				int ndex = index(grid[x+i],n);
				for (int h=0; h<length(neighbors); h++)
					if (ndex==neighbors[h]) goto skip;
				resize(neighbors,length(neighbors)+1);
				neighbors[length(neighbors)-1] = ndex;
				skip: continue;
			}
		return neighbors;
	}
};

MMSP::vector<int> neighbors(const sparsePF1D& GRID, int x) {return GRID.neighbors(x);}

template <> void copy(sparsePF1D& GRID1, const sparsePF1D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(sparsePF1D& GRID1, sparsePF1D& GRID2) {GRID1.swap(GRID2);}


class sparsePF2D : public grid<2,sparse<double> > {
public:
	// constructors
	sparsePF2D(int x, int y, int GHOSTS=1)
		: grid<2,sparse<double> >(0,0,x,0,y,GHOSTS) {}
	sparsePF2D(int x0, int x1, int y0, int y1, int GHOSTS=1)
		: grid<2,sparse<double> >(0,x0,x1,y0,y1,GHOSTS) {}
	sparsePF2D(const sparsePF2D& GRID, int FIELDS)
		: grid<2,sparse<double> >(GRID,FIELDS) {}
	sparsePF2D(const char* filename, int GHOSTS=1)
		: grid<2,sparse<double> >(filename,GHOSTS) {} 

	// utility functions
	MMSP::vector<int> neighbors(int x, int y) const
	{
		MMSP::vector<int> neighbors;
		const sparsePF2D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++)
				for (int n=0; n<length(grid[x+i][y+j]); n++) {
					int ndex = index(grid[x+i][y+j],n);
					for (int h=0; h<length(neighbors); h++)
						if (ndex==neighbors[h]) goto skip;
					resize(neighbors,length(neighbors)+1);
					neighbors[length(neighbors)-1] = ndex;
					skip: continue;
				}
		return neighbors;
	}
};

MMSP::vector<int> neighbors(const sparsePF2D& GRID, int x, int y) {return GRID.neighbors(x,y);}

template <> void copy(sparsePF2D& GRID1, const sparsePF2D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(sparsePF2D& GRID1, sparsePF2D& GRID2) {GRID1.swap(GRID2);}


class sparsePF3D : public grid<3,sparse<double> > {
public:
	// constructors
	sparsePF3D(int x, int y, int z, int GHOSTS=1)
		: grid<3,sparse<double> >(0,0,x,0,y,0,z,GHOSTS) {}
	sparsePF3D(int x0, int x1, int y0, int y1, int z0, int z1, int GHOSTS=1)
		: grid<3,sparse<double> >(0,x0,x1,y0,y1,z0,z1,GHOSTS) {}
	sparsePF3D(const sparsePF3D& GRID, int FIELDS)
		: grid<3,sparse<double> >(GRID,FIELDS) {}
	sparsePF3D(const char* filename, int GHOSTS=1)
		: grid<3,sparse<double> >(filename,GHOSTS) {}

	// utility functions
	MMSP::vector<int> neighbors(int x, int y, int z) const
	{
		MMSP::vector<int> neighbors;
		const sparsePF3D& grid = *this;
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++)
				for (int k=-1; k<=1; k++)
					for (int n=0; n<length(grid[x+i][y+j][z+k]); n++) {
						int ndex = index(grid[x+i][y+j][z+k],n);
						for (int h=0; h<length(neighbors); h++)
							if (ndex==neighbors[h]) goto skip;
						resize(neighbors,length(neighbors)+1);
						neighbors[length(neighbors)-1] = ndex;
						skip: continue;
					}
		return neighbors;
	}
};

MMSP::vector<int> neighbors(const sparsePF3D& GRID, int x, int y, int z) {return GRID.neighbors(x,y,z);}

template <> void copy(sparsePF3D& GRID1, const sparsePF3D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(sparsePF3D& GRID1, sparsePF3D& GRID2) {GRID1.swap(GRID2);}

} // namespace MMSP

#endif
