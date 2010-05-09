// graingrowth.cpp 
// Anisotropic phase field grain growth example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"graingrowth.hpp"

std::string PROGRAM = "graingrowth";
std::string MESSAGE = "Anisotropic phase field grain growth example code";

typedef MMSP::grid<2,MMSP::vector<double> > GRID2D;
typedef MMSP::grid<3,MMSP::vector<double> > GRID3D;

#include"MMSP.main.hpp"

