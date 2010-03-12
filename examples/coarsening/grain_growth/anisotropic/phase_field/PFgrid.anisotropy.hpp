// PFgrid.anisotropy.hpp
// Anisotropic coarsening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef PFGRID_UPDATE
#define PFGRID_UPDATE
#include"PFgrid.hpp"
#include"anisotropy.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::PFgrid2D grid(2,128,128);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++) {
				grid[x][y][0] = 0.0;
				grid[x][y][1] = 0.0;
				double d = sqrt(pow(64.0-x,2)+pow(64.0-y,2));
				if (d<32.0) grid[x][y][1] = 1.0;
				else grid[x][y][0] = 1.0;
			}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::PFgrid3D grid(2,64,64,64);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++) {
					grid[x][y][z][0] = 0.0;
					grid[x][y][z][1] = 0.0;
					double d = sqrt(pow(32.0-x,2)+pow(32.0-y,2)+pow(32.0-z,2));
					if (d<16.0) grid[x][y][z][1] = 1.0;
					else grid[x][y][z][0] = 1.0;
				}

		MMSP::output(grid,filename);
	}
}

void update(PFgrid2D& grid, int steps)
{
	PFgrid2D update(grid);
	const double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {

				double S = 0.0;
				vector<double> s(fields(grid),0.0);
				for (int i=0; i<fields(grid); i++) {
					for (int h=-1; h<=1; h++)
						for (int k=-1; k<=1; k++)
							if (grid[x+h][y+k][i]>epsilon) {
								s[i] = 1.0;
								goto next;
							}
					next: S += s[i];
				}

				vector<double> lap(fields(grid),0.0);
				for (int i=0; i<fields(grid); i++)
					if (s[i]>0.0)
						lap[i] = (grid[x+1][y][i]-2.0*grid[x][y][i]+grid[x-1][y][i])/(dx(grid)*dx(grid))
						        +(grid[x][y+1][i]-2.0*grid[x][y][i]+grid[x][y-1][i])/(dy(grid)*dy(grid));

				vector<double> dFdp(fields(grid),0.0);
				for (int i=0; i<fields(grid); i++)
					if (s[i]>0.0)
						for (int j=0; j<fields(grid); j++)
							if (s[j]>0.0) {
								double gamma = energy(i,j);
								double width = 8.0;
								double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
								double w = 4.0*gamma/width;
								dFdp[i] += 0.5*eps*eps*lap[j]+w*grid[x][y][j];
								dFdp[j] += 0.5*eps*eps*lap[i]+w*grid[x][y][i];
								for (int k=0; k<fields(grid); k++)
									if (s[k]>0.0) {
										dFdp[i] += 3.0*grid[x][y][j]*grid[x][y][k];
										dFdp[j] += 3.0*grid[x][y][k]*grid[x][y][i];
										dFdp[k] += 3.0*grid[x][y][i]*grid[x][y][j];
									}
							}

				vector<double> dpdt(fields(grid),0.0);
				for (int i=0; i<fields(grid); i++)
					if (s[i]>0.0)
						for (int j=0; j<fields(grid); j++)
							if (s[j]>0.0) {
								double mu = mobility(i,j);
								dpdt[i] -= mu*(dFdp[i]-dFdp[j]);
								dpdt[j] -= mu*(dFdp[j]-dFdp[i]);
							}

				for (int i=0; i<fields(grid); i++) {
					double dt = 0.02;
					update[x][y][i] = grid[x][y][i]+dt*(2.0/S)*dpdt[i];
					if (update[x][y][i]<epsilon) update[x][y][i] = 0.0;
					else if (update[x][y][i]>1.0) update[x][y][i] = 1.0;
				}
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(PFgrid3D& grid, int steps)
{
	PFgrid3D update(grid);
	const double epsilon = 1.0e-8;

	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {

					double S = 0.0;
					vector<double> s(fields(grid),0.0);
					for (int i=0; i<fields(grid); i++) {
						for (int h=-1; h<=1; h++)
							for (int k=-1; k<=1; k++)
								for (int l=-1; l<=1; l++)
									if (grid[x+h][y+k][z+l][i]>epsilon) {
										s[i] = 1.0;
										goto next;
									}
						next: S += s[i];
					}

					vector<double> lap(fields(grid),0.0);
					for (int i=0; i<fields(grid); i++)
						if (s[i]>0.0)
							lap[i] = (grid[x+1][y][z][i]-2.0*grid[x][y][z][i]+grid[x-1][y][z][i])/(dx(grid)*dx(grid))
									+(grid[x][y+1][z][i]-2.0*grid[x][y][z][i]+grid[x][y-1][z][i])/(dy(grid)*dy(grid))
									+(grid[x][y][z+1][i]-2.0*grid[x][y][z][i]+grid[x][y][z-1][i])/(dz(grid)*dz(grid));

					vector<double> dFdp(fields(grid),0.0);
					for (int i=0; i<fields(grid); i++)
						if (s[i]>0.0)
							for (int j=0; j<fields(grid); j++)
								if (s[j]>0.0) {
									double gamma = energy(i,j);
									double width = 8.0;
									double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
									double w = 4.0*gamma/width;
									dFdp[i] += 0.5*eps*eps*lap[j]+w*grid[x][y][z][j];
									dFdp[j] += 0.5*eps*eps*lap[i]+w*grid[x][y][z][i];
									for (int k=0; k<fields(grid); k++)
										if (s[k]>0.0) {
											dFdp[i] += 3.0*grid[x][y][z][j]*grid[x][y][z][k];
											dFdp[j] += 3.0*grid[x][y][z][k]*grid[x][y][z][i];
											dFdp[k] += 3.0*grid[x][y][z][i]*grid[x][y][z][j];
										}
								}

					vector<double> dpdt(fields(grid),0.0);
					for (int i=0; i<fields(grid); i++)
						if (s[i]>0.0)
							for (int j=0; j<fields(grid); j++)
								if (s[j]>0.0) {
									double mu = mobility(i,j);
									dpdt[i] -= mu*(dFdp[i]-dFdp[j]);
									dpdt[j] -= mu*(dFdp[j]-dFdp[i]);
								}

					for (int i=0; i<fields(grid); i++) {
						double dt = 0.02;
						update[x][y][z][i] = grid[x][y][z][i]+dt*(2.0/S)*dpdt[i];
						if (update[x][y][z][i]<epsilon) update[x][y][z][i] = 0.0;
						else if (update[x][y][z][i]>1.0) update[x][y][z][i] = 1.0;
					}
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
