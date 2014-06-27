// Uncomment the following definition if you have gnuplot installed on 
// your system and you've want the gnuplot_i library linked into
// this code.  This will permit interactive display of your results.
// If you don't have gnuplot or gnuplot_i just leave that line commented
// out.
//
// Add a line here about what the makefile or compilation command would 
// look like.  Given the complexity of MMSP it might be OK to assume that 
// folks can, in general, compile code.
//

#define GNUPLOT_I

// On systems where you can apt-get packages (or similar) make sure that
// your installation of gnuplot is the x11 compatible version.  On a debian
// system you would: apt-get install gnuplot-x11

#ifdef GNUPLOT_I
#include "gnuplot_i.hpp"
#endif

#include "../../include/MMSP.hpp"
#include <math.h>
using namespace MMSP;

int main()
{
  // Here we set the critical parameters of the simulation.
  int number_points = 100;
  float dx = 0.1;
  float conc_lhs = 1.0;
  float conc_rhs = 0.0;
  float simulation_time = 1.0;
  float const diffusion_coefficient = 1.0;
  float coefficient;
  float raw_laplacian;

  // Here we set the grid parameters.
  int x_lower_bound = 0;
  int x_upper_bound = number_points-1;

#ifdef GNUPLOT_I
  std::vector<double> data;
  Gnuplot conc_plot;
  char t;
#endif

  // Things that get computed at runtime are done here.
  // The timestep is determined by the Courant condition.
  //
  float dt = (dx*dx)/(2*diffusion_coefficient);
  //
  // The length of the domain is computed from the grid spacing and the
  // number of points used.
  //
  float length = dx*number_points;
  
  int timesteps = simulation_time/dx;

  // Create the data structures to hold the concentration values.  Two 
  // grids are needed for storing the "new" and the "old" values of
  // concentration.
  //
  // ToDo:  Short explanation of the data structure here.
  //
  MMSP::grid< 1, scalar <float> > concentration_new(0,x_lower_bound,x_upper_bound);
  MMSP::grid< 1, scalar <float> > concentration_old(0,x_lower_bound,x_upper_bound);
  
  // We use two functions that are part of the grid class for accessing
  // information about the grid.  x0() gives us the global lower bound 
  // of the grid passed as an argument.  x1() gives us the global upper
  // bound of the grid passed as an argument.
  //
  // In this part of the code we are setting the initial conditions.
  //
  for (int x = MMSP::x0(concentration_new); x < MMSP::x1(concentration_new); x++) {
    if (x < x_upper_bound/2) {
      concentration_new[x]=conc_lhs;
      concentration_old[x]=conc_lhs;
    }
    else {
      concentration_new[x]=conc_rhs;
      concentration_old[x]=conc_rhs;
    }
  }
  // Run the simulation for n timesteps.
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
    coefficient = diffusion_coefficient*dt;
    for (int i = x0(concentration_old)+1; i < x1(concentration_old)-1; i++) {
      // We compute the laplacian on the grid points here.  Later this 
      // will be replaced by a function in MMSP.
      raw_laplacian = (concentration_old[i-1]-2*concentration_old[i]+concentration_old[i+1]);
      concentration_new[i]=coefficient*raw_laplacian + concentration_old[i];
      // This closes the spatial iterating loop.
    }
    // Set the boundary conditions here.  In this example the boundary conditions
    // are no-flux.
    concentration_new[x0(concentration_new)] = concentration_new[x0(concentration_new)+1];
    concentration_new[x0(concentration_new)-1] = concentration_new[x0(concentration_new)-2];
    // Swap old and new concentrations for the next timestep.
    swap(concentration_old, concentration_new);
    std::cout << "Timestep: " << n << " \n";
  }

    // This is where we print the output.
#ifdef GNUPLOT_I
    // If you have gnuplot_i on your system, then do these lines.
  for (int i = x0(concentration_old); i < x1(concentration_old); i++)
    {
      data.push_back((double)concentration_old[i]);
    }
  conc_plot.reset_plot();
  conc_plot.cmd("set nolabel");
  conc_plot.cmd("set yrange[-0.1:1.1]");
  conc_plot.plot_x(data, "Concentration");
  std::cout << "Press a key then press enter: ";
  std::cin >> t;
  //
#else
  // ELSE print the data to the terminal.
  for (int i = x0(concentration_old); i < x1(concentration_old)-1; i++) {
    // print the data.
    std::cout << concentration_old[i] << '\n';
  }
  //
#endif

  Finalize();
  // return success.
  return 0;
  // This bracket ends main()
}

