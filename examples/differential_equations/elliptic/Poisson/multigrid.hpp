// multigrid.hpp
// multigrid functionality for MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MULTIGRID
#define MULTIGRID

#include"MMSP.grid.hpp"

namespace MMSP {

template <typename T> void coarsen(MMSP::grid<2,T>& u, int stride, std::string method="full-weighting")
{
	// coarsen from stride s to stride 2s
	int s = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);

	if (method=="injection") return;

	else if (method=="full-weighting")
		for (int x=2*s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s)
				u[x][y] = 0.2500*u[x][y]
				         +0.1250*(u[x+s][y]+u[x-s][y]+u[x][y+s]+u[x][y-s])
				         +0.0625*(u[x+s][y+s]+u[x-s][y+s]+u[x+s][y-s]+u[x-s][y-s]);
}

template <typename T> void coarsen(MMSP::grid<3,T>& u, int stride, std::string method="full-weighting")
{
	// coarsen from stride s to stride 2s
	int s = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	int z1 = MMSP::z1(u);

	if (method=="injection") return;

	else if (method=="full-weighting")
		for (int x=2*s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s)
				for (int z=2*s; z<z1-1; z+=2*s)
					u[x][y][z] = 0.125000*u[x][y][z]
					            +0.062500*(u[x+s][y][z]+u[x-s][y][z]+u[x][y+s][z]+u[x][y-s][z]+u[x][y][z+s]+u[x][y][z-s])
					            +0.031250*(u[x+s][y+s][z]+u[x-s][y+s][z]+u[x][y+s][z+s]+u[x][y-s][z+s]+u[x+s][y][z+s]+u[x+s][y][z-s])
					            +0.031250*(u[x+s][y-s][z]+u[x-s][y-s][z]+u[x][y+s][z-s]+u[x][y-s][z-s]+u[x-s][y][z+s]+u[x-s][y][z-s])
					            +0.015625*(u[x+s][y+s][z+s]+u[x+s][y+s][z-s]+u[x+s][y-s][z+s]+u[x-s][y+s][z+s])
					            +0.015625*(u[x+s][y-s][z-s]+u[x-s][y+s][z-s]+u[x-s][y-s][z+s]+u[x-s][y-s][z-s]);
}

template <typename T> void refine(MMSP::grid<2,T>& u, int stride, std::string method="linear")
{
	// refine from stride 2s to stride s
	int s = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);

	if (method=="linear") {
		for (int x=s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s)
				u[x][y] = 0.5*(u[x+s][y]+u[x-s][y]);
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=2*s)
				u[x][y] = 0.5*(u[x][y+s]+u[x][y-s]);
	}

	else if (method=="cubic") {
		for (int x=s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s) {
				if (x==s)
					u[x][y] = 5.0/16.0*(u[x-s][y]-u[x+3*s][y])+15.0/16.0*u[x+s][y]+1.0/16.0*u[x+5*s][y];
				else if (x==x1-1-s)
					u[x][y] = 1.0/16.0*u[x-5*s][y]+15.0/16.0*u[x-s][y]+5.0/16.0*(u[x+s][y]-u[x-3*s][y]);
				else
					u[x][y] = 9.0/16.0*(u[x-s][y]+u[x+s][y])-1.0/16.0*(u[x-3*s][y]+u[x+3*s][y]);
			}
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=2*s) {
				if (y==s)
					u[x][y] = 5.0/16.0*(u[x][y-s]-u[x][y+3*s])+15.0/16.0*u[x][y+s]+1.0/16.0*u[x][y+5*s];
				else if (y==y1-1-s)
					u[x][y] = 1.0/16.0*u[x][y-5*s]+15.0/16.0*u[x][y-s]+5.0/16.0*(u[x][y+s]-u[x][y-3*s]);
				else
					u[x][y] = 9.0/16.0*(u[x][y-s]+u[x][y+s])-1.0/16.0*(u[x][y-3*s]+u[x][y+3*s]);
			}
	}
}
template <typename T> void refine(MMSP::grid<3,T>& u, int stride, std::string method="linear")
{
	// refine from stride 2s to stride s
	int s = stride;
	int x1 = MMSP::x1(u);
	int y1 = MMSP::y1(u);
	int z1 = MMSP::z1(u);

	if (method=="linear") {
		for (int x=s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s)
				for (int z=2*s; z<z1-1; z+=2*s)
					u[x][y][z] = 0.5*(u[x+s][y][z]+u[x-s][y][z]);
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=2*s)
				for (int z=2*s; z<z1-1; z+=2*s)
					u[x][y][z] = 0.5*(u[x][y+s][z]+u[x][y-s][z]);
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				for (int z=s; z<z1-1; z+=2*s)
					u[x][y][z] = 0.5*(u[x][y][z+s]+u[x][y][z-s]);
	}

	else if (method=="cubic") {
		for (int x=s; x<x1-1; x+=2*s)
			for (int y=2*s; y<y1-1; y+=2*s)
				for (int z=2*s; z<z1-1; z+=2*s) {
					if (x==s)
						u[x][y][z] = 5.0/16.0*(u[x-s][y][z]-u[x+3*s][y][z])+15.0/16.0*u[x+s][y][z]+1.0/16.0*u[x+5*s][y][z];
					else if (x==x1-1-s)
						u[x][y][z] = 1.0/16.0*u[x-5*s][y][z]+15.0/16.0*u[x-s][y][z]+5.0/16.0*(u[x+s][y][z]-u[x-3*s][y][z]);
					else
						u[x][y][z] = 9.0/16.0*(u[x-s][y][z]+u[x+s][y][z])-1.0/16.0*(u[x-3*s][y][z]+u[x+3*s][y][z]);
				}
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=2*s)
				for (int z=2*s; z<z1-1; z+=2*s) {
					if (y==s)
						u[x][y][z] = 5.0/16.0*(u[x][y-s][z]-u[x][y+3*s][z])+15.0/16.0*u[x][y+s][z]+1.0/16.0*u[x][y+5*s][z];
					else if (y==y1-1-s)
						u[x][y][z] = 1.0/16.0*u[x][y-5*s][z]+15.0/16.0*u[x][y-s][z]+5.0/16.0*(u[x][y+s][z]-u[x][y-3*s][z]);
					else
						u[x][y][z] = 9.0/16.0*(u[x][y-s][z]+u[x][y+s][z])-1.0/16.0*(u[x][y-3*s][z]+u[x][y+3*s][z]);
				}
		for (int x=s; x<x1-1; x+=s)
			for (int y=s; y<y1-1; y+=s)
				for (int z=s; z<z1-1; z+=2*s) {
					if (z==s)
						u[x][y][z] = 5.0/16.0*(u[x][y][z-s]-u[x][y][z+3*s])+15.0/16.0*u[x][y][z+s]+1.0/16.0*u[x][y][z+5*s];
					else if (z==z1-1-s)
						u[x][y][z] = 1.0/16.0*u[x][y][z-5*s]+15.0/16.0*u[x][y][z-s]+5.0/16.0*(u[x][y][z+s]-u[x][y][z-3*s]);
					else
						u[x][y][z] = 9.0/16.0*(u[x][y][z-s]+u[x][y][z+s])-1.0/16.0*(u[x][y][z-3*s]+u[x][y][z+3*s]);
				}
	}
}

template <int dim, typename T> 
void MG(MMSP::grid<dim,T>& u, const MMSP::grid<dim,T>& f, int stride, int gamma=1, int nu1=2, int nu2=2)
{
	// standard multigrid cycle 
	int s = stride;

	// solve if at coarsest level
	bool solve = false;
	for (int i=0; i<dim; i++)
		if (2*s+1==x1(u,i))
			solve = true;
	if (solve)
		return smooth(u,f,s,40);

	// presmooth
	smooth(u,f,s,nu1);

	// compute defect
	MMSP::grid<dim,T> d(u);
	defect(u,f,d,s);

	// restrict defect
	coarsen(d,s);

	// solve for correction
	MMSP::grid<dim,T> v(u);
	v = static_cast<T>(0.0);

	// multigrid iteration
	for (int k=0; k<gamma; k++)
		MG(v,d,2*s,gamma,nu1,nu2);

	// interpolate correction
	refine(v,s);

	// compute corrected approximation
	u += v;

	// postsmooth
	smooth(u,f,s,nu2);
}

template <int dim, typename T> 
void FMG(MMSP::grid<dim,T>& u, const MMSP::grid<dim,T>& f, int gamma=1, int nu1=2, int nu2=2)
{
	// solve at coarsest level
	int s = (x1(u,0)-1)/4;
	for (int i=1; i<dim; i++)
		if ((x1(u,i)-1)/4<s)
			s = (x1(u,i)-1)/4;
	smooth(u,f,s,40);

	// solve at finer levels
	for (s=s/2; s>0; s/=2) {
		refine(u,s,"cubic");
		MG(u,f,s,gamma,nu1,nu2);
	}
}

} // namespace MMSP

#endif
