// FDgrid.hpp
// Class definitions for rectilinear finite difference grids
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef FDGRID
#define FDGRID
#include"MMSP.grid.hpp"

namespace MMSP{

class FDgrid1D : public grid<1,scalar<double> > {
public:
	// constructors
	FDgrid1D(int x, int GHOSTS=1)
		: grid<1,scalar<double> >(1,0,x,GHOSTS) {}
	FDgrid1D(int x0, int x1, int GHOSTS=1)
		: grid<1,scalar<double> >(1,x0,x1,GHOSTS) {}
	FDgrid1D(const char* filename, int GHOSTS=1)
		: grid<1,scalar<double> >(filename,GHOSTS) {}
};

template <> void copy(FDgrid1D& GRID1, const FDgrid1D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(FDgrid1D& GRID1, FDgrid1D& GRID2) {GRID1.swap(GRID2);}


class FDgrid2D : public grid<2,scalar<double> > {
public:
	// constructors
	FDgrid2D(int x, int y, int GHOSTS=1)
		: grid<2,scalar<double> >(1,0,x,0,y,GHOSTS) {}
	FDgrid2D(int x0, int x1, int y0, int y1, int GHOSTS=1)
		: grid<2,scalar<double> >(1,x0,x1,y0,y1,GHOSTS) {}
	FDgrid2D(const char* filename, int GHOSTS=1)
		: grid<2,scalar<double> >(filename,GHOSTS) {} 
};

template <> void copy(FDgrid2D& GRID1, const FDgrid2D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(FDgrid2D& GRID1, FDgrid2D& GRID2) {GRID1.swap(GRID2);}


class FDgrid3D : public grid<3,scalar<double> > {
public:
	// constructors
	FDgrid3D(int x, int y, int z, int GHOSTS=1)
		: grid<3,scalar<double> >(1,0,x,0,y,0,z,GHOSTS) {}
	FDgrid3D(int x0, int x1, int y0, int y1, int z0, int z1, int GHOSTS=1)
		: grid<3,scalar<double> >(1,x0,x1,y0,y1,z0,z1,GHOSTS) {}
	FDgrid3D(const char* filename, int GHOSTS=1)
		: grid<3,scalar<double> >(filename,GHOSTS) {}
};

template <> void copy(FDgrid3D& GRID1, const FDgrid3D& GRID2) {GRID1.copy(GRID2);}
template <> void swap(FDgrid3D& GRID1, FDgrid3D& GRID2) {GRID1.swap(GRID2);}

} // namespace MMSP

#endif
