// poisson.hpp
// smooth() and defect() functions for multigrid solution of Poisson equation
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef POISSON
#define POISSON

#include"multigrid.hpp"

namespace MMSP {

template <typename T>
void smooth(MMSP::grid<2,T>& u, const MMSP::grid<2,T>& f, int stride, int iterations=1)
{
	// red-black Gauss-Seidel iteration
	int s = stride;
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double x1 = MMSP::x1(u);
	double y1 = MMSP::y1(u);
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double w   = 1.0/(2.0/dx2+2.0/dy2);
	double wx  = w/dx2;
	double wy  = w/dy2;

	for (int i=0; i<iterations; i++) {
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				if ((x+y)%2==0)
					u[x][y] = -w*f[x][y]+wx*(u[x-s][y]+u[x+s][y])
					                    +wy*(u[x][y-s]+u[x][y+s]);
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				if ((x+y)%2==1)
					u[x][y] = -w*f[x][y]+wx*(u[x-s][y]+u[x+s][y])
					                    +wy*(u[x][y-s]+u[x][y+s]);
	}
}

template <typename T>
void smooth(MMSP::grid<3,T>& u, const MMSP::grid<3,T>& f, int stride, int iterations=1)
{
	// red-black Gauss-Seidel iteration
	int s = stride;
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double dz = s*MMSP::dz(u);
	double x1 = MMSP::x1(u);
	double y1 = MMSP::y1(u);
	double z1 = MMSP::z1(u);
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dz2 = dz*dz;
	double w   = 1.0/(2.0/dx2+2.0/dy2+2.0/dz2);
	double wx  = w/dx2;
	double wy  = w/dy2;
	double wz  = w/dz2;

	for (int i=0; i<iterations; i++) {
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				for (int z=s; z<z1-1; z+=s)
					if ((x+y+z)%2==0)
						u[x][y][z] = -w*f[x][y][z]+wx*(u[x-s][y][z]+u[x+s][y][z])
												  +wy*(u[x][y-s][z]+u[x][y+s][z])
												  +wz*(u[x][y][z-s]+u[x][y][z+s]);
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				for (int z=s; z<z1-1; z+=s)
					if ((x+y+z)%2!=0)
						u[x][y][z] = -w*f[x][y][z]+wx*(u[x-s][y][z]+u[x+s][y][z])
												  +wy*(u[x][y-s][z]+u[x][y+s][z])
												  +wz*(u[x][y][z-s]+u[x][y][z+s]);
	}
}

template <typename T>
void defect(const MMSP::grid<2,T>& u, const MMSP::grid<2,T>& f, MMSP::grid<2,T>& d, int stride)
{
	// compute defect for Poisson equation lap(u) = f
	int s = stride;
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double x1 = MMSP::x1(u);
	double y1 = MMSP::y1(u);
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double wx  = 1.0/dx2;
	double wy  = 1.0/dy2;

	for (int x=s; x<x1-1; x+=s)
		for (int y=s; y<y1-1; y+=s)
			d[x][y] = f[x][y]-wx*(u[x-s][y]-2.0*u[x][y]+u[x+s][y])
			                 -wy*(u[x][y-s]-2.0*u[x][y]+u[x][y+s]);

	for (int x=s; x<x1-1; x+=s) {
		d[x][0] = 0.0;
		d[x][y1-1] = 0.0;
	}
	for (int y=s; y<y1-1; y+=s) {
		d[0][y] = 0.0;
		d[x1-1][y] = 0.0;
	}
}

template <typename T>
void defect(const MMSP::grid<3,T>& u, const MMSP::grid<3,T>& f, MMSP::grid<3,T>& d, int stride)
{
	// compute defect for Poisson equation lap(u) = f
	int s = stride;
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double dz = s*MMSP::dz(u);
	double x1 = MMSP::x1(u);
	double y1 = MMSP::y1(u);
	double z1 = MMSP::z1(u);
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dz2 = dz*dz;
	double wx  = 1.0/dx2;
	double wy  = 1.0/dy2;
	double wz  = 1.0/dz2;

	for (int x=s; x<x1-1; x+=s)
		for (int y=s; y<y1-1; y+=s)
			for (int z=s; z<z1-1; z+=s)
				d[x][y][z] = f[x][y][z]-wx*(u[x-s][y][z]-2.0*u[x][y][z]+u[x+s][y][z])
			                           -wy*(u[x][y-s][z]-2.0*u[x][y][z]+u[x][y+s][z])
			                           -wz*(u[x][y][z-s]-2.0*u[x][y][z]+u[x][y][z+s]);

	for (int x=s; x<x1-1; x+=s)
		for (int y=s; y<y1-1; y+=s) {
			d[x][y][0] = 0.0;
			d[x][y][z1-1] = 0.0;
		}
	for (int y=s; y<y1-1; y+=s)
		for (int z=s; z<z1-1; z+=s) {
			d[0][y][z] = 0.0;
			d[x1-1][y][z] = 0.0;
		}
	for (int z=s; z<z1-1; z+=s)
		for (int x=s; x<x1-1; x+=s) {
			d[x][0][z] = 0.0;
			d[x][y1-1][z] = 0.0;
		}
}

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::grid<2,double> grid(1,0,129,0,129);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++) {
				double X = double(x)/129.0;
				double Y = double(y)/129.0;
				grid[x][y] = exp(X*Y);
			}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,double> grid(1,0,65,0,65,0,65);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++) {
					double X = double(x)/129.0;
					double Y = double(y)/129.0;
					double Z = double(z)/129.0;
					grid[x][y][z] = exp(X*Y*Z);
				}

		MMSP::output(grid,filename);
	}
}

void update(MMSP::grid<2,double>& grid, int steps)
{
	const MMSP::grid<2,double>& f = grid;

	MMSP::grid<2,double> u(f);
	int x0 = MMSP::x0(u);
	int x1 = MMSP::x1(u);
	int y0 = MMSP::y0(u);
	int y1 = MMSP::y1(u);

	for (int x=x0; x<x1; x++)
		for (int y=y0; y<y1; y++)
			u[x][y] = 0.0;

	FMG(u,f,1,2,2);

	MMSP::swap(grid,u);
}

void update(MMSP::grid<3,double>& grid, int steps)
{
	const MMSP::grid<3,double>& f = grid;

	MMSP::grid<3,double> u(f);
	int x0 = MMSP::x0(u);
	int x1 = MMSP::x1(u);
	int y0 = MMSP::y0(u);
	int y1 = MMSP::y1(u);
	int z0 = MMSP::z0(u);
	int z1 = MMSP::z1(u);

	for (int x=x0; x<x1; x++)
		for (int y=y0; y<y1; y++)
			for (int z=z0; z<z1; z++)
				u[x][y][z] = 0.0;

	FMG(u,f,1,2,2);

	MMSP::swap(grid,u);
}

} // namespace MMSP

#endif
