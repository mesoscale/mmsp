// CAgrid.hpp
// Class definitions for rectilinear cellular automata grids
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAGRID
#define CAGRID
#include"MMSP.grid.hpp"

namespace MMSP{

class CAgrid1D : public grid<1,scalar<int> > {
public:
	// constructors
	CAgrid1D(int x, int GHOSTS=1)
		: grid<1,scalar<int> >(1,0,x,GHOSTS) {}
	CAgrid1D(int x0, int x1, int GHOSTS=1)
		: grid<1,scalar<int> >(1,x0,x1,GHOSTS) {}
	CAgrid1D(const CAgrid1D& GRID, int FIELDS)
		: grid<1,scalar<int> >(GRID,FIELDS) {}
	CAgrid1D(const char* filename, int GHOSTS=1)
		: grid<1,scalar<int> >(filename,GHOSTS) {}
};

template <> void copy(CAgrid1D& GRID1, const CAgrid1D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(CAgrid1D& GRID1, CAgrid1D& GRID2) {GRID1.swap(GRID2);}


class CAgrid2D : public grid<2,scalar<int> > {
public:
	// constructors
	CAgrid2D(int x, int y, int GHOSTS=1)
		: grid<2,scalar<int> >(1,0,x,0,y,GHOSTS) {}
	CAgrid2D(int x0, int x1, int y0, int y1, int GHOSTS=1)
		: grid<2,scalar<int> >(1,x0,x1,y0,y1,GHOSTS) {}
	CAgrid2D(const CAgrid2D& GRID, int FIELDS)
		: grid<2,scalar<int> >(GRID,FIELDS) {}
	CAgrid2D(const char* filename, int GHOSTS=1)
		: grid<2,scalar<int> >(filename,GHOSTS) {} 
};

template <> void copy(CAgrid2D& GRID1, const CAgrid2D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(CAgrid2D& GRID1, CAgrid2D& GRID2) {GRID1.swap(GRID2);}


class CAgrid3D : public grid<3,scalar<int> > {
public:
	// constructors
	CAgrid3D(int x, int y, int z, int GHOSTS=1)
		: grid<3,scalar<int> >(1,0,x,0,y,0,z,GHOSTS) {}
	CAgrid3D(int x0, int x1, int y0, int y1, int z0, int z1, int GHOSTS=1)
		: grid<3,scalar<int> >(1,x0,x1,y0,y1,z0,z1,GHOSTS) {}
	CAgrid3D(const CAgrid3D& GRID, int FIELDS)
		: grid<3,scalar<int> >(GRID,FIELDS) {}
	CAgrid3D(const char* filename, int GHOSTS=1)
		: grid<3,scalar<int> >(filename,GHOSTS) {}
};

template <> void copy(CAgrid3D& GRID1, const CAgrid3D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(CAgrid3D& GRID1, CAgrid3D& GRID2) {GRID1.swap(GRID2);}

} // namespace MMSP

#endif
