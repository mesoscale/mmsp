Mesoscale Microstructure Simulation Project
====

The goal of the Mesoscale Microstructure Simulation Project (MMSP) is to provide a simple,
consistent, and extensible programming interface for all grid and mesh based microstructure
evolution methods. Simple means that the package has a very small learning curve, and for
most routine simulations, only a minimal amount of code must be written. By consistent we
mean, for example, that code for two-dimensional simulations is nearly identical to that
for three-dimensional simulations, single processor programs are easily parallelized, and
fundamentally different methods like Monte Carlo or phase field have the same look and feel.
Finally, extensible means that it's straightforward to add new grid types or physical behaviors
to the package. Other considerations include efficiency and portability (MMSP is written
entirely in ISO compliant c++). For more details, see the documentation.

MMSP is nothing more than a collection of c++ header files that declare a number of grid objects
(classes) and define how most of their methods (member functions) are implemented.

*Some things MMSP provides include:*

 * A simple, extensible programming interface
 * Computational grids of arbitrary dimension
 * Parallel implementations using MPI
 * Automatic, optimal parallel mesh topologies
 * Utility programs for grid visualization
 * Monte Carlo methods
 * Cellular automata methods
 * Phase field methods (conventional)
 * Phase field methods (sparsePF)
 * General finite difference PDE solvers
 * 22+ example problems that run in 2D and 3D, single and parallel

*Typical MMSP applications include:*

 * Grain growth and coarsening
 * Precipitation reactions
 * Crystal growth and solidification
 * Lattice based kinetic Monte Carlo
 * Statistical mechanics: Ising model, classical Heisenberg model, etc.
 * Spinodal decomposition and other second order transformations

*MMSP requires:*

 * Minimal programming experience
 * An ISO compliant c++ compiler (e.g. GCC 2.95 or later)
 * zlib libraries for data compression (e.g. zlib 1.2.7)
 * MPI libraries are required if compiling parallel programs (e.g. Open MPI )

*Documentation*

The MMSP manual is a work in progress. It is currently the only source for detailed documentation about MMSP.

*Contact us*

The administrators for the MMSP source code are Jason Gruber (gruberja@gmail.com) and Trevor Keller (trevor.keller@gmail.com). Please do not hesitate to send questions or comments.


This work is supported in part by the US NSF under award #1056704 through the Metals and Metallic Nanostructures Program, Division of Materials Research. 
