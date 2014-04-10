// anisotropy.hpp
// Anisotropic energy and mobility functions
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ANISOTROPY
#define ANISOTROPY
#include<map>

// global energy and mobility storage
namespace anisotropy{
	std::map<int,std::map<int,float> > energy_table;
	std::map<int,std::map<int,float> > mobility_table;
}

template <typename T> T min(const T& a, const T& b) {return (a<b?a:b);}
template <typename T> T max(const T& a, const T& b) {return (a>b?a:b);}

float energy(int i, int j)
{
	using namespace anisotropy; 

	// trivial case: no boundary
	if (i==j) return 0.0;

	// use computed value, if possible
	int a = min(i,j);
	int b = max(i,j);
	float energy = energy_table[a][b];
	if (energy==0.0) {
		// compute energy here...
		energy = 1.0;
		energy_table[a][b] = energy;
	}
	return energy;
}

float mobility(int i, int j)
{
	using namespace anisotropy; 

	// trivial case: no boundary
	if (i==j) return 0.0;

	// use computed value, if possible
	int a = min(i,j);
	int b = max(i,j);
	float mobility = mobility_table[a][b];
	if (mobility==0.0) {
		// compute mobility here...
		mobility = 1.0;
		mobility_table[a][b] = mobility;
	}
	return mobility;
}

#endif
