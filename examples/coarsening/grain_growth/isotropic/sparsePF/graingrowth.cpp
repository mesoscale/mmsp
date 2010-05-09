// graingrowth.cpp 
// Isotropic sparse phase field (sparsePF) grain growth example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"graingrowth.hpp"

std::string PROGRAM = "graingrowth";
std::string MESSAGE = "Isotropic sparse phase field (sparsePF) grain growth example code";

typedef MMSP::grid<2,MMSP::sparse<double> > GRID2D;
typedef MMSP::grid<3,MMSP::sparse<double> > GRID3D;

#include"MMSP.main.hpp"

