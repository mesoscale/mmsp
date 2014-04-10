// ostwald.hpp
// Ostwald ripening algorithms for 2D and 3D phase field methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef OSTWALD_UPDATE
#define OSTWALD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"ostwald.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::grid<2,MMSP::vector<double> > grid(2,0,128,0,128);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++) {
				double r = double(rand())/double(RAND_MAX);
				grid[x][y][0] = r;
				grid[x][y][1] = 1.0-r;
			}

		MMSP::output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,MMSP::vector<double> > grid(2,0,64,0,64,0,64);
		int x0 = MMSP::x0(grid);
		int x1 = MMSP::x1(grid);
		int y0 = MMSP::y0(grid);
		int y1 = MMSP::y1(grid);
		int z0 = MMSP::z0(grid);
		int z1 = MMSP::z1(grid);

		for (int x=x0; x<x1; x++)
			for (int y=y0; y<y1; y++)
				for (int z=z0; z<z1; z++) {
					double r = double(rand())/double(RAND_MAX);
					grid[x][y][z][0] = r;
					grid[x][y][z][1] = 1.0-r;
				}

		MMSP::output(grid,filename);
	}
}

void update(MMSP::grid<2,MMSP::vector<double> >& grid, int steps)
{
	MMSP::grid<2,MMSP::vector<double> > update(grid);
	MMSP::grid<2,double> wspace(grid,1);

	double dt = 0.01;
	double L = 1.0;
	double D = 1.0;

	double Calpha = 0.05;
	double Cbeta = 0.95;
	double Cmatrix = 0.5*(Calpha+Cbeta);
	double A = 2.0;
	double B = A/pow(Cbeta-Cmatrix,2);
	double gamma = 2.0/pow(Cbeta-Calpha,2);
	double delta = 1.0;
	double epsilon = 3.0;
	double Dalpha = gamma/pow(delta,2);
	double Dbeta = gamma/pow(delta,2);
	double kappa = 2.0;


	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				double sum = 0.0;
				for (int i=1; i<fields(grid); i++)
					sum += pow(grid[x][y][i],2);

				double C = grid[x][y][0];
				double lap = 
					  (grid[x+1][y][0]-2.0*grid[x][y][0]+grid[x-1][y][0])/(dx(grid)*dx(grid))
					 +(grid[x][y+1][0]-2.0*grid[x][y][0]+grid[x][y-1][0])/(dy(grid)*dy(grid));

				wspace[x][y] = -A*(C-Cmatrix)+B*pow(C-Cmatrix,3)
				               +Dalpha*pow(C-Calpha,3)+Dbeta*pow(C-Cbeta,3)
				               -gamma*(C-Calpha)*sum-kappa*lap;
			}
		ghostswap(wspace);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++) {
				if (1) {
					double lap = 
					  (wspace[x+1][y]-2.0*wspace[x][y]+wspace[x-1][y])/(dx(grid)*dx(grid))
					 +(wspace[x][y+1]-2.0*wspace[x][y]+wspace[x][y-1])/(dy(grid)*dy(grid));

					update[x][y][0] = grid[x][y][0]+dt*D*lap;
				}

				double sum = 0.0;
				for (int i=1; i<fields(grid); i++)
					sum += pow(grid[x][y][i],2);

				for (int i=1; i<fields(grid); i++) {
					double C = grid[x][y][0];
					double value = grid[x][y][i];
					double lap = 
						  (grid[x+1][y][i]-2.0*grid[x][y][i]+grid[x-1][y][i])/(dx(grid)*dx(grid))
						 +(grid[x][y+1][i]-2.0*grid[x][y][i]+grid[x][y-1][i])/(dy(grid)*dy(grid));

					update[x][y][i] = grid[x][y][i]-dt*L*(-gamma*pow(C-Calpha,2)+delta*pow(value,3)
					                                      +epsilon*value*(sum-pow(value,2))-kappa*lap);
				}
			}
		swap(grid,update);
		ghostswap(grid);
	}
}

void update(MMSP::grid<3,MMSP::vector<double> >& grid, int steps)
{
	MMSP::grid<3,MMSP::vector<double> > update(grid);
	MMSP::grid<3,double> wspace(grid,1);

	double dt = 0.01;
	double L = 1.0;
	double D = 1.0;

	double Calpha = 0.05;
	double Cbeta = 0.95;
	double Cmatrix = 0.5*(Calpha+Cbeta);
	double A = 2.0;
	double B = A/pow(Cbeta-Cmatrix,2);
	double gamma = 2.0/pow(Cbeta-Calpha,2);
	double delta = 1.0;
	double epsilon = 3.0;
	double Dalpha = gamma/pow(delta,2);
	double Dbeta = gamma/pow(delta,2);
	double kappa = 2.0;


	for (int step=0; step<steps; step++) {
		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					double sum = 0.0;
					for (int i=1; i<fields(grid); i++)
						sum += pow(grid[x][y][z][i],2);

					double C = grid[x][y][z][0];
					double lap = 
						  (grid[x+1][y][z][0]-2.0*grid[x][y][z][0]+grid[x-1][y][z][0])/(dx(grid)*dx(grid))
						 +(grid[x][y+1][z][0]-2.0*grid[x][y][z][0]+grid[x][y-1][z][0])/(dy(grid)*dy(grid))
						 +(grid[x][y][z+1][0]-2.0*grid[x][y][z][0]+grid[x][y][z-1][0])/(dz(grid)*dz(grid));

					wspace[x][y][z] = -A*(C-Cmatrix)+B*pow(C-Cmatrix,3)
					                  +Dalpha*pow(C-Calpha,3)+Dbeta*pow(C-Cbeta,3)
								      -gamma*(C-Calpha)*sum-kappa*lap;
				}
		ghostswap(wspace);

		for (int x=x0(grid); x<x1(grid); x++)
			for (int y=y0(grid); y<y1(grid); y++)
				for (int z=z0(grid); z<z1(grid); z++) {
					if (1) {
						double lap = 
						  (wspace[x+1][y][z]-2.0*wspace[x][y][z]+wspace[x-1][y][z])/(dx(grid)*dx(grid))
						 +(wspace[x][y+1][z]-2.0*wspace[x][y][z]+wspace[x][y-1][z])/(dy(grid)*dy(grid))
						 +(wspace[x][y][z+1]-2.0*wspace[x][y][z]+wspace[x][y][z-1])/(dz(grid)*dz(grid));

						update[x][y][z][0] = grid[x][y][z][0]+dt*D*lap;
					}

					double sum = 0.0;
					for (int i=1; i<fields(grid); i++)
						sum += pow(grid[x][y][z][i],2);

					for (int i=1; i<fields(grid); i++) {
						double C = grid[x][y][z][0];
						double value = grid[x][y][z][i];
						double lap = 
							  (grid[x+1][y][z][i]-2.0*grid[x][y][z][i]+grid[x-1][y][z][i])/(dx(grid)*dx(grid))
							 +(grid[x][y+1][z][i]-2.0*grid[x][y][z][i]+grid[x][y-1][z][i])/(dy(grid)*dy(grid))
							 +(grid[x][y][z+1][i]-2.0*grid[x][y][z][i]+grid[x][y][z-1][i])/(dz(grid)*dz(grid));

						update[x][y][z][i] = grid[x][y][z][i]-dt*L*(-gamma*pow(C-Calpha,2)+delta*pow(value,3)
																	+epsilon*value*(sum-pow(value,2))-kappa*lap);
					}
				}
		swap(grid,update);
		ghostswap(grid);
	}

}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
