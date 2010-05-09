// zener.hpp
// Anisotropic coarsening algorithms for 2D and 3D Monte Carlo methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ZENER_UPDATE
#define ZENER_UPDATE 
#include"MMSP.hpp"
#include"anisotropy.hpp"
#include<cmath>

namespace MMSP{

void generate(int dim, const char* filename)
{
	if (dim==2) {
		MMSP::grid<2,int> grid(0,0,128,0,128);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(64.0-x[0],2)+pow(64.0-x[1],2));
			if (d<32.0) grid(i) = 2;
			else grid(i) = 1;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					grid[x][y] = 0;
		}

		output(grid,filename);
	}

	if (dim==3) {
		MMSP::grid<3,int> grid(0,0,64,0,64,0,64);

		for (int i=0; i<nodes(grid); i++) {
			vector<int> x = position(grid,i);
			double d = sqrt(pow(32.0-x[0],2)+pow(32.0-x[1],2)+pow(32.0-x[2],2));
			if (d<16.0) grid(i) = 2;
			else grid(i) = 1;
		}

		for (int j=0; j<50; j++) {
			int i = rand()%nodes(grid);
			vector<int> p = position(grid,i);
			for (int x=p[0]-1; x<=p[0]+1; x++)
				for (int y=p[1]-1; y<=p[1]+1; y++)
					for (int z=p[2]-1; z<=p[2]+1; z++)
						grid[x][y][z] = 0;
		}

		output(grid,filename);
	}
}

void update(MMSP::grid<2,int>& grid, int steps)
{
	const double kT = 0.50;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node 
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);

			if (spin1!=0) {
				// determine neighboring spins
				sparse<bool> neighbors;
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++) {
						int spin = grid[x[0]+i][x[1]+j];
						set(neighbors,spin) = true;
					}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
					// compute energy change
					double dE = -E(spin1,spin2);
					for (int i=-1; i<=1; i++)
						for (int j=-1; j<=1; j++)
							for (int k=-1; k<=1; k++) {
								int spin = grid[x[0]+i][x[1]+j][x[2]+k];
								dE += E(spin,spin2)-E(spin,spin1);
							}

					// compute boundary energy, mobility
					double E = energy(spin1,spin2);
					double M = mobility(spin1,spin2);

					// attempt a spin flip
					double r = double(rand())/double(RAND_MAX);
					if (dE<=0.0 and r<M*E) grid(p) = spin2;
					if (dE>0.0 and r<M*E*exp(-dE/(E*kT))) grid(p) = spin2;
				}
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

void update(MMSP::grid<3,int>& grid, int steps)
{
	const double kT = 0.75;
	int gss = int(sqrt(nodes(grid)));

	for (int step=0; step<steps; step++) {
		for (int h=0; h<nodes(grid); h++) {
			// choose a random node 
			int p = rand()%nodes(grid);
			vector<int> x = position(grid,p);
			int spin1 = grid(p);

			if (spin1!=0) {
				// determine neighboring spins
				sparse<bool> neighbors;
				for (int i=-1; i<=1; i++)
					for (int j=-1; j<=1; j++)
						for (int k=-1; k<=1; k++) {
							int spin = grid[x[0]+i][x[1]+j][x[2]+k];
							set(neighbors,spin) = true;
						}

				// choose a random neighbor spin
				int spin2 = index(neighbors,rand()%length(neighbors));

				if (spin1!=spin2 and spin2!=0) {
					// compute energy change
					double dE = -E(spin1,spin2);
					for (int i=-1; i<=1; i++)
						for (int j=-1; j<=1; j++)
							for (int k=-1; k<=1; k++) {
								int spin = grid[x[0]+i][x[1]+j][x[2]+k];
								dE += E(spin,spin2)-E(spin,spin1);
							}

					// compute boundary energy, mobility
					double E = energy(spin1,spin2);
					double M = mobility(spin1,spin2);

					// attempt a spin flip
					double r = double(rand())/double(RAND_MAX);
					if (dE<=0.0 and r<M*E) grid(p) = spin2;
					if (dE>0.0 and r<M*E*exp(-dE/(E*kT))) grid(p) = spin2;
				}
			}
			if (h%gss==0) ghostswap(grid);
		}
	}
}

} // namespace MC

#endif
