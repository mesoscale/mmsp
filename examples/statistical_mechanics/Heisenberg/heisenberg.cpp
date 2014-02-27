// heisenberg.cpp 
// Classical Heisenberg model using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"heisenberg.hpp"

std::string PROGRAM = "heisenberg";
std::string MESSAGE = "Classical Heisenberg model using MMSP";

typedef MMSP::grid<2,MMSP::vector<double> > GRID2D;
typedef MMSP::grid<3,MMSP::vector<double> > GRID3D;

#include"MMSP.main.hpp"

