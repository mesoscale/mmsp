// zener.cpp
// Isotropic phase field grain growth with Zener pinning example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"zener.hpp"

std::string PROGRAM = "zener";
std::string MESSAGE = "Isotropic phase field grain growth with Zener pinning example code";

typedef MMSP::grid<1,MMSP::vector<double> > GRID1D;
typedef MMSP::grid<2,MMSP::vector<double> > GRID2D;
typedef MMSP::grid<3,MMSP::vector<double> > GRID3D;

#include"MMSP.main.hpp"

