// energy.hpp
// Energy functions for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_ENERGY
#define CAHNHILLIARD_ENERGY
#include<cmath>

const double deltaX = 1.0;
const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca + Cb);
const double A = 2.0;
const double B = A/((Ca-Cm)*(Ca-Cm));
const double D = 2.0/(Cb-Ca);
const double K = 2.0;
const double CFL = 0.25;
const double dt = std::pow(deltaX, 4)*CFL/(32.0*D*K);

double energydensity(const double& C)
{
	return -0.5*A*pow(C-Cm,2) + 0.25*B*pow(C-Cm,4) + 0.25*Ca*pow(C-Ca,4) + 0.25*Cb*pow(C-Cb,4);
}

double dfdc(const double& C)
{
    return -A*(C-Cm) + B*pow(C-Cm, 3) + Ca*pow(C-Ca, 3) + Cb*pow(C-Cb, 3);
}

#endif
