// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"CAgrid.hpp"
#include"MCgrid.hpp"
#include"PFgrid.hpp"
#include"sparsePF.hpp"
#include<cmath>

int main(int argc, char* argv[])
{
	MMSP::Init(argc,argv);

	// circle geometry (MCgrid)
	if (1) {
		// create the grid object
		MMSP::MCgrid2D grid(0,100,0,100);

		// generate a circle
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++) {
			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2));

			// set cell values
			if (d<=25.0) grid[x][y] = 1;
			else grid[x][y] = 0;
		}

		// output the grid to a file
		output(grid,"circle.MC");
	}

	// circle geometry (PFgrid)
	if (1) {
		// create the grid object
		MMSP::PFgrid2D grid(2,0,100,0,100);

		// generate a circle
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++) {
			// set default field values
			grid[x][y][0] = 0.0;
			grid[x][y][1] = 0.0;

			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2));

			// set field values
			if (d<=25.0) grid[x][y][1] = 1.0;
			else grid[x][y][0] = 1.0;
		}

		// output the grid to a file
		output(grid,"circle.PF");
	}

	// circle geometry (sparsePF)
	if (1) {
		// create the grid object
		MMSP::sparsePF2D grid(0,100,0,100);

		// generate a circle
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++) {
			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2));

			// set field values
			if (d<=25.0) MMSP::set(grid[x][y],1) = 1.0;
			else MMSP::set(grid[x][y],0) = 1.0;
		}

		// output the grid to a file
		output(grid,"circle.sPF");
	}

	// sphere geometry (MCgrid)
	if (1) {
		// create the grid object
		MMSP::MCgrid3D grid(0,100,0,100,0,100);

		// generate a sphere
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++)
		for (int z=0; z<100; z++) {
			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2)+pow(49.5-z,2));

			// set cell values
			if (d<=25.0) grid[x][y][z] = 1;
			else grid[x][y][z] = 0;
		}

		// output the grid to a file
		output(grid,"sphere.MC");
	}

	// sphere geometry (PFgrid)
	if (1) {
		// create the grid object
		MMSP::PFgrid3D grid(2,0,100,0,100,0,100);

		// generate a sphere
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++)
		for (int z=0; z<100; z++) {
			// set default field values
			grid[x][y][z][0] = 0.0;
			grid[x][y][z][1] = 0.0;

			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2)+pow(49.5-z,2));

			// set field values
			if (d<=25.0) grid[x][y][z][1] = 1.0;
			else grid[x][y][z][0] = 1.0;
		}

		// output the grid to a file
		output(grid,"sphere.PF");
	}

	// sphere geometry (sparsePF)
	if (1) {
		// create the grid object
		MMSP::sparsePF3D grid(0,100,0,100,0,100);

		// generate a sphere
		for (int x=0; x<100; x++)
		for (int y=0; y<100; y++)
		for (int z=0; z<100; z++) {
			// compute distance from center of grid
			double d = sqrt(pow(49.5-x,2)+pow(49.5-y,2)+pow(49.5-z,2));

			// set field values
			if (d<=25.0) MMSP::set(grid[x][y][z],1) = 1.0;
			else MMSP::set(grid[x][y][z],0) = 1.0;
		}

		// output the grid to a file
		output(grid,"sphere.sPF");
	}

	MMSP::Finalize();
}
