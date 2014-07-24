// zener.cpp
// Isotropic Monte Carlo grain growth with Zener pinning example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"zener.hpp"

std::string PROGRAM = "zener";
std::string MESSAGE = "Isotropic Monte Carlo grain growth with Zener pinning example code";

typedef MMSP::grid<1,int> GRID1D;
typedef MMSP::grid<2,int> GRID2D;
typedef MMSP::grid<3,int> GRID3D;

#include"MMSP.main.hpp"

