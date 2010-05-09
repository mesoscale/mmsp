// zener.cpp 
// Anisotropic sparse phase field (sparsePF) Zener pinning example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"zener.hpp"

std::string PROGRAM = "zener";
std::string MESSAGE = "Anisotropic sparse phase field (sparsePF) Zener pinning example code";

typedef MMSP::grid<2,MMSP::sparse<double> > GRID2D;
typedef MMSP::grid<3,MMSP::sparse<double> > GRID3D;

#include"MMSP.main.hpp"

