// energy.hpp
// Energy parameters and density function for 2D Cahn-Hilliard model
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#ifndef CAHNHILLIARD_ENERGY
#define CAHNHILLIARD_ENERGY
#include<cmath>

const double Ca = 0.05;
const double Cb = 0.95;
const double Cm = 0.5*(Ca + Cb);      // = 0.5
const double A = 2.0;
const double B = A/((Ca-Cm)*(Ca-Cm)); // = 9.8765

template<typename T>
double energydensity(const T& C)
{
	return -0.5*A*pow(C-Cm,2) + 0.25*B*pow(C-Cm,4) + 0.25*Ca*pow(C-Ca,4) + 0.25*Cb*pow(C-Cb,4);
}

template<typename T>
double dfdc(const T& C)
{
    return -A*(C-Cm) + B*pow(C-Cm, 3) + Ca*pow(C-Ca, 3) + Cb*pow(C-Cb, 3);
}

#endif
