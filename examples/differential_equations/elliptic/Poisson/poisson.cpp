// poisson.cpp
// smooth() and defect() functions for multigrid solution of Poisson equation
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef POISSON
#define POISSON

#include"multigrid.hpp"
#include"poisson.hpp"

namespace MMSP {

template <typename T>
void smooth(grid<1,T>& u, const grid<1,T>& f, int stride, int iterations=1)
{
	// red-black Gauss-Seidel iteration
	int s  = stride;
	int x1 = MMSP::x1(u);
	double dx = s*MMSP::dx(u);
	double dx2 = dx*dx;
	double w   = 1.0/(2.0/dx2);
	double wx  = w/dx2;

	for (int i=0; i<iterations; i++) {
		for (int x=s; x<x1-1; x+=s)
			if (x%2==0)
				u[x] = -w*f[x]+wx*(u[x-s]+u[x+s]);

		for (int x=s; x<x1-1; x+=s)
			if (x%2==1)
				u[x] = -w*f[x]+wx*(u[x-s]+u[x+s]);
	}
}

template <typename T>
void smooth(grid<2,T>& u, const grid<2,T>& f, int stride, int iterations=1)
{
	// red-black Gauss-Seidel iteration
	int s  = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
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
void smooth(grid<3,T>& u, const grid<3,T>& f, int stride, int iterations=1)
{
	// red-black Gauss-Seidel iteration
	int s  = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	int z1 = MMSP::z1(u);
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double dz = s*MMSP::dz(u);
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
void defect(const grid<1,T>& u, const grid<1,T>& f, grid<1,T>& d, int stride)
{
	// compute defect for Poisson equation lap(u) = f
	int s  = stride;
	int x1 = MMSP::x1(u);
	double dx = s*MMSP::dx(u);
	double dx2 = dx*dx;
	double wx  = 1.0/dx2;

	for (int x=s; x<x1-1; x+=s)
		d[x] = f[x]-wx*(u[x-s]-2.0*u[x]+u[x+s]);

	for (int x=s; x<x1-1; x+=s) {
		d[0] = 0.0;
		d[x1-1] = 0.0;
	}
}

template <typename T>
void defect(const grid<2,T>& u, const grid<2,T>& f, grid<2,T>& d, int stride)
{
	// compute defect for Poisson equation lap(u) = f
	int s  = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
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
void defect(const grid<3,T>& u, const grid<3,T>& f, grid<3,T>& d, int stride)
{
	// compute defect for Poisson equation lap(u) = f
	int s  = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	int z1 = MMSP::z1(u);
	double dx = s*MMSP::dx(u);
	double dy = s*MMSP::dy(u);
	double dz = s*MMSP::dz(u);
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
	if (dim==1) {
		grid<1,double> initGrid(1,0,2049);

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			double X = double(x[0])/128.0;
			initGrid(n) = exp(X);
		}

		output(initGrid,filename);
	}

	if (dim==2) {
		int L=513;
		grid<2,double> initGrid(1,0,L,0,L);

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			double X = double(x[0])/L;
			double Y = double(x[1])/L;
			initGrid(n) = exp(X*Y);
		}

		output(initGrid,filename);
	}

	if (dim==3) {
		int L=65;
		grid<3,double> initGrid(1,0,L,0,L,0,L);

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			double X = double(x[0])/L;
			double Y = double(x[1])/L;
			double Z = double(x[2])/L;
			initGrid(n) = exp(X*Y*Z);
		}

		output(initGrid,filename);
	}
}

template <int dim, typename T> void update(grid<dim,T>& multiGrid, int steps)
{
	const grid<dim,T>& f = multiGrid;

	grid<dim,T> u(f);

	for (int n=0; n<nodes(u); n++)
		u(n) = 0.0;

	FMG(u,f,1,2,2);

	swap(multiGrid,u);
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
