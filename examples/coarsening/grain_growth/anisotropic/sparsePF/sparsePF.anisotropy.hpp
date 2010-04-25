// sparsePF.anisotropy.hpp
// Anisotropic coarsening algorithms for 2D and 3D sparse phase field (sparsePF) methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef SPARSEPF_UPDATE
#define SPARSEPF_UPDATE
#include"sparsePF.hpp"
#include"anisotropy.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::sparsePF2D grid(128,128);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++) {
				if (x<32) {
					if (y<64) MMSP::set(grid[x][y],2) = 1.0;
					else MMSP::set(grid[x][y],3) = 1.0;
				}
				else if (x>96) {
					if (y<64) MMSP::set(grid[x][y],2) = 1.0;
					else MMSP::set(grid[x][y],3) = 1.0;
				}
				else {
					if (y<32 or y>96) MMSP::set(grid[x][y],1) = 1.0;
					else MMSP::set(grid[x][y],0) = 1.0;
				}
			}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::sparsePF3D grid(64,64,64);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++) {
					if (x<16) {
						if (y<32) MMSP::set(grid[x][y][z],2) = 1.0;
						else MMSP::set(grid[x][y][z],3) = 1.0;
					}
					else if (x>48) {
						if (y<32) MMSP::set(grid[x][y][z],2) = 1.0;
						else MMSP::set(grid[x][y][z],3) = 1.0;
					}
					else {
						if (y<16 or y>48) MMSP::set(grid[x][y][z],1) = 1.0;
						else MMSP::set(grid[x][y][z],0) = 1.0;
					}
				}

		MMSP::output(grid,filename);
	}
}

void update(sparsePF2D& grid, int steps)
{
	for (int step=0; step<steps; step++) {
		sparsePF2D update(grid);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				vector<int> nbors = neighbors(grid,x,y);
				if (length(nbors)==1)
					set(update[x][y],nbors[0]) = grid[x][y][nbors[0]];

				else {
					double S = static_cast<double>(length(nbors));

					sparse<double> lap;
					for (int i=0; i<length(nbors); i++) {
						int index = nbors[i];
						set(lap,index) = (grid[x+1][y][index]-2.0*grid[x][y][index]+grid[x-1][y][index])/(dx(grid)*dx(grid))
						                +(grid[x][y+1][index]-2.0*grid[x][y][index]+grid[x][y-1][index])/(dy(grid)*dy(grid));
					}

					sparse<double> dFdp;
					for (int i=0; i<length(nbors); i++)
						for (int j=i+1; j<length(nbors); j++) {
							int index1 = nbors[i];
							int index2 = nbors[j];
							double gamma = energy(index1,index2);
							const double width = 8.0;
							double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
							double w = 4.0*gamma/width;
							set(dFdp,index1) += 0.5*eps*eps*lap[index2]+w*grid[x][y][index2];
							set(dFdp,index2) += 0.5*eps*eps*lap[index1]+w*grid[x][y][index1];
							for (int k=j+1; k<length(nbors); k++) {
								int index3 = nbors[k];
								set(dFdp,index1) += 6.0*grid[x][y][index2]*grid[x][y][index3];
								set(dFdp,index2) += 6.0*grid[x][y][index3]*grid[x][y][index1];
								set(dFdp,index3) += 6.0*grid[x][y][index1]*grid[x][y][index2];
							}
						}

					sparse<double> dpdt;
					for (int i=0; i<length(nbors); i++)
						for (int j=i+1; j<length(nbors); j++) {
							int index1 = nbors[i];
							int index2 = nbors[j];
							double mu = mobility(index1,index2);
							set(dpdt,index1) -= mu*(dFdp[index1]-dFdp[index2]);
							set(dpdt,index2) -= mu*(dFdp[index2]-dFdp[index1]);
						}

					for (int i=0; i<length(nbors); i++) {
						const double dt = 0.02;
						const double epsilon = 1.0e-8;
						int index = nbors[i];
						double value = grid[x][y][index]+dt*(2.0/S)*dpdt[index];
						if (value>1.0) value = 1.0;
						if (value>epsilon) set(update[x][y],index) = value;
					}
				}
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(sparsePF3D& grid, int steps)
{
	for (int step=0; step<steps; step++) {
		sparsePF3D update(grid);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					vector<int> nbors = neighbors(grid,x,y,z);
					if (length(nbors)==1)
						set(update[x][y][z],nbors[0]) = grid[x][y][z][nbors[0]];

					else {
						double S = static_cast<double>(length(nbors));

						sparse<double> lap;
						for (int i=0; i<length(nbors); i++) {
							int index = nbors[i];
							set(lap,index) = (grid[x+1][y][z][index]-2.0*grid[x][y][z][index]+grid[x-1][y][z][index])/(dx(grid)*dx(grid))
											+(grid[x][y+1][z][index]-2.0*grid[x][y][z][index]+grid[x][y-1][z][index])/(dy(grid)*dy(grid))
											+(grid[x][y][z+1][index]-2.0*grid[x][y][z][index]+grid[x][y][z-1][index])/(dy(grid)*dy(grid));
						}

						sparse<double> dFdp;
						for (int i=0; i<length(nbors); i++)
							for (int j=i+1; j<length(nbors); j++) {
								int index1 = nbors[i];
								int index2 = nbors[j];
								double gamma = energy(index1,index2);
								const double width = 8.0;
								double eps = 4.0/acos(-1.0)*sqrt(0.5*gamma*width);
								double w = 4.0*gamma/width;
								set(dFdp,index1) += 0.5*eps*eps*lap[index2]+w*grid[x][y][z][index2];
								set(dFdp,index2) += 0.5*eps*eps*lap[index1]+w*grid[x][y][z][index1];
								for (int k=j+1; k<length(nbors); k++) {
									int index3 = nbors[k];
									set(dFdp,index1) += 3.0*grid[x][y][z][index2]*grid[x][y][z][index3];
									set(dFdp,index2) += 3.0*grid[x][y][z][index3]*grid[x][y][z][index1];
									set(dFdp,index3) += 3.0*grid[x][y][z][index1]*grid[x][y][z][index2];
								}
							}

						sparse<double> dpdt;
						for (int i=0; i<length(nbors); i++)
							for (int j=i+1; j<length(nbors); j++) {
								int index1 = nbors[i];
								int index2 = nbors[j];
								double mu = mobility(index1,index2);
								set(dpdt,index1) -= mu*(dFdp[index1]-dFdp[index2]);
								set(dpdt,index2) -= mu*(dFdp[index2]-dFdp[index1]);
							}

						for (int i=0; i<length(nbors); i++) {
							const double dt = 0.02;
							const double epsilon = 1.0e-8;
							int index = nbors[i];
							double value = grid[x][y][z][index]+dt*(2.0/S)*dpdt[index];
							if (value>1.0) value = 1.0;
							if (value>epsilon) set(update[x][y][z],index) = value;
						}
					}
				}
		swap(grid,update);
		ghostswap(grid);
	}
}

} // namespace MMSP

#endif
