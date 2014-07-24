// zener.cpp
// isotropic sparse phase field (sparsePF) model for Zener pinning using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"zener.hpp"

std::string PROGRAM = "zener";
std::string MESSAGE = "isotropic sparse phase field (sparsePF) model for Zener pinning using MMSP";

typedef MMSP::grid<1,MMSP::sparse<double> > GRID1D;
typedef MMSP::grid<2,MMSP::sparse<double> > GRID2D;
typedef MMSP::grid<3,MMSP::sparse<double> > GRID3D;

#include"MMSP.main.hpp"

