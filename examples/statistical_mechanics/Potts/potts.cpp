// potts.cpp
// Q-state Potts model using MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include"potts.hpp"

std::string PROGRAM = "potts";
std::string MESSAGE = "Q-state Potts model using MMSP";

typedef MMSP::grid<2,int> GRID2D;
typedef MMSP::grid<3,int> GRID3D;

#include"MMSP.main.hpp"

