// heisenberg.cpp
// 2D and 3D Heisenberg model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef HEISENBERG_UPDATE
#define HEISENBERG_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"heisenberg.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==1) {
		MMSP::grid<1,vector<double> > grid(3,0,128);

		for (int i=0; i<nodes(grid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			grid(i)[0] = cos(psi)*sin(theta);
			grid(i)[1] = sin(psi)*sin(theta);
			grid(i)[2] = cos(theta);
		}

		output(grid,filename);
	}

	if (dim==2) {
		MMSP::grid<2,vector<double> > grid(3,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			grid(i)[0] = cos(psi)*sin(theta);
			grid(i)[1] = sin(psi)*sin(theta);
			grid(i)[2] = cos(theta);
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,vector<double> > grid(3,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			double psi = 2.0*acos(-1.0)*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			grid(i)[0] = cos(psi)*sin(theta);
			grid(i)[1] = sin(psi)*sin(theta);
			grid(i)[2] = cos(theta);
		}

		output(grid,filename);
	}
}

void update(grid<2,vector<double> >& grid, int steps)
{
	double J = 1.0;
	double kT = 0.50;
	double pi = acos(-1.0);

	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random site
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			vector<double>& s1 = grid(p);

			// choose a random unit vector 
			vector<double> s2(3);
			double psi = 2.0*pi*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			s2[0] = cos(psi)*sin(theta);
			s2[1] = sin(psi)*sin(theta);
			s2[2] = cos(theta);

			// compute energy change
			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++) {
					vector<double>& s = grid[x[0]+i][x[1]+j];
					sum += s[0]*(s1[0]-s2[0])+s[1]*(s1[1]-s2[1])+s[2]*(s1[2]-s2[2]);
				}
			double dE = -J*sum;

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid(p) = s2;
			else if (r<exp(-dE/kT)) grid(p) = s2;

			if (h%gss==0) ghostswap(grid);
		}
	}
}

void update(grid<3,vector<double> >& grid, int steps)
{
	double J = 1.0;
	double kT = 0.75;
	double pi = acos(-1.0);

	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random site
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			vector<double>& s1 = grid(p);

			// choose a random unit vector
			vector<double> s2(3);
			double psi = 2.0*pi*double(rand())/double(RAND_MAX);
			double theta = acos(1.0-2.0*double(rand())/double(RAND_MAX));
			s2[0] = cos(psi)*sin(theta);
			s2[1] = sin(psi)*sin(theta);
			s2[2] = cos(theta);

			// compute energy change
			double sum = -1.0;
			for (int i=-1; i<=1; i++)
				for (int j=-1; j<=1; j++)
					for (int k=-1; k<=1; k++) {
						vector<double>& s = grid[x[0]+i][x[1]+j][x[2]+k];
						sum += s[0]*(s1[0]-s2[0])+s[1]*(s1[1]-s2[1])+s[2]*(s1[2]-s2[2]);
					}
			double dE = -J*sum;

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) grid(p) = s2;
			else if (r<exp(-dE/kT)) grid(p) = s2;

			if (h%gss==0) ghostswap(grid);
		}
	}
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"
