#include"MMSP.hpp"
#include<math.h>
using namespace MMSP;


void main()
{
  // Here we set the critical parameters of the simulation.
  int number_points = 100;
  float dx = 0.1;
  float conc_lhs = 1.0;
  float conc_rhs = 0.0;
  float simulation_time = 1.0;
  float diffusion_coefficient = 1.0;

  // Things that get computed at runtime are done here.
  // The timestep is determined by considering the Courant condition.
  //
  float dt;
  dt = (dx*dx)/(2*diffusion_coefficient);
  //
  // The length of the domain is computed from the grid spacing and the
  // number of points used.
  //
  float length;
  length = dx*number_points;
  
  int timesteps;
  timesteps = simulation_time/dx;
  // Does this automagially truncate the fractional part?

  // Create the data structures to hold the concentration values.  Two grids
  // are needed for storing the "new" and the "old" values of concentration.
  //
  // Short explanation of the data structure here.
  //
  MMSP::grid<1,scalar<float>> concentration_new(1,x_lower_bound,x_upper_bound);
  MMSP::grid<1,scalar<float>> concentration_old(1,x_lower_bound,x_upper_bound);
  
  // We use two functions that are part of the grid class for accesing information
  // about the grid.  x0() gives us the global lower bound of the grid passed as
  // an argument.  x1() gives us the global upper bound of the grid passed as an
  // argument.
  //
  // In this part of the code we are setting the initial conditions.
  //
  for (int x = MMSP::x0(concentration_new); x < MMSP::x1(concentration_new); x++)
    if (x < x_upper_bound/2) {
      concentration_new[x]=conc_lhs;
      concentration_old[x]=conc_lhs;
    }
    else {
      concentration_new[x]=conc_rhs;
      concentration_old[x]=conc_rhs;
    }

  // Run the simulation for i timesteps.
  for (int n = 0; n < timesteps; n++) {

    // Iterate through the spatial coordinate on the i gridpoints.
    // Here we are using the functions x0 and x1 defined as part of the 
    // grid class.
    //
    // Note that we leave the first and last grid points alone.  These are
    // our boundary grid points, or ghost cells that will be updated after
    // each timestep.
    //
    // We compute the coefficient for the discrete equation.  We need only
    // compute this once.
    //
    coefficient = diffusion_coefficient*dt
    for (int i = x0(concentration_old)+1; i < x1(concentration_old)-1; i++) {
      // We compute the laplacian on the grid points here.  Later this will be
      // replaced by a function in MMSP.
      laplacian = (concentration_old[x-1]-2*concentration_old[x]+concentration_old[x+1])
      concentration_new[x]=coefficient*laplacian + concentration_old[x];
      // This closes the spatial iterating loop.
    }
    // Set the boundary conditions here.  In this example the boundary conditions
    // are no-flux.
    concentration_new[x0(concentration_new)] = concentration_new[x0(concentration_new)+1]
    concentration_new[x0(concentration_new)-1] = concentration_new[x0(concentration_new)-2]
    // Swap old and new concentrations for the next timestep.
    swap(concentration_old, concentration_new);
    // This closes the timestep loop.
  }
  // This prints the results of the grid to cout.  We print the old value here
  // because of the swap up above.
  for (int i = x0(concentration_old)+1; i < x1(concentration_old)-1; i++) {
    std::cout<<concentration_old[x]<< "\n";
    // Finishes the output loop.
  }
  Finalize();
  // This bracket ends main()
}

