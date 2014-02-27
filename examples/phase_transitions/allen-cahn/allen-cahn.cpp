// allen-cahn.cpp 
// Example program for the Allen-Cahn model using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"allen-cahn.hpp"

std::string PROGRAM = "allen-cahn";
std::string MESSAGE = "Allen-Cahn model using MMSP";

typedef MMSP::grid<2,double> GRID2D;
typedef MMSP::grid<3,double> GRID3D;

#include"MMSP.main.hpp"

