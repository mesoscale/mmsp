#include <stdio.h>
#include "MMSP.hpp"
using namespace MMSP;

int main(int argc, char* argv[])
{
  Init(argc, argv);

  int nx;
  int iterations;
  float dt, diffusionCoefficient, dx;

  nx = 10;
  iterations = 100;
  diffusionCoefficient = 1.0;
  dx = 1.0;
  dt = dx*dx/diffusionCoefficient/4;
	  
  grid<1,scalar<double> > GRID(1,0,nx);
  grid<1,scalar<double> > GRID2(1,0,nx);

  for (int x=x0(GRID); x<x1(GRID); x++)
    if (x<nx/2) {
      GRID[x]=1;
      GRID2[x]=1;
    } else {
      GRID[x]=0;
      GRID2[x]=0;
    }

  for (int i=0; i<iterations; i++) {
    for (int x=x0(GRID); x<x1(GRID); x++) {
      if (x==0 || x==nx-1) {
      }
      else {
	GRID2[x]=(diffusionCoefficient*dt/dx/dx)*(GRID[x-1]-2*GRID[x]+GRID[x+1])+GRID[x];
      }
    }
    swap(GRID,GRID2);
  }

  //This prints the results of the grid to cout
  for (int x=x0(GRID); x<x1(GRID); x++) {
    std::cout<<GRID[x]<<std::endl;
  }

  Finalize();

  return 0;
}

