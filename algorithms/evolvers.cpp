// evolvers.cpp
// Iplementations of methods to evolve MMSP grids
#ifndef _EVOLVERS_CPP_
#define _EVOLVERS_CPP_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <limits>
#include <cassert>

#include "evolvers.h"

void print_progress(const int step, const int steps, const int iterations)
{
	 char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout << "•] "
							<<std::setw(2)<<std::right<<deltat/3600<<"h:"
							<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
							<<std::setw(2)<<std::right<<deltat%60<<"s"
							<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}

template<int dim, typename T>
void interface_field_evolves(const T& dt, const T& width, const T& gamma, const T& epsilon, const T& w, const T& mu,
														 const MMSP::grid<dim,MMSP::sparse<T> >& former, MMSP::grid<dim,MMSP::sparse<T> >& latter, const int node)
{
	MMSP::vector<int> x = position(former, node);

	// determine nonzero fields within
	// the neighborhood of this node
	// (2 adjacent voxels along each cardinal direction)
	MMSP::sparse<int> s;
	for (int j = 0; j < dim; j++)
		for (int k = -1; k <= 1; k++) {
			x[j] += k;
			for (int h = 0; h < MMSP::length(former(x)); h++) {
				int index = MMSP::index(former(x), h);
				set(s, index) = 1;
			}
			x[j] -= k;
		}
	T S = T(length(s));

	// if only one field is nonzero,
	// then copy this node to latter grid
	if (S < 2.0) latter(node) = former(node);
	else {
		// compute laplacian of each field
		MMSP::sparse<T> lap = MMSP::laplacian(former, node);

		// compute variational derivatives
		MMSP::sparse<T> dFdp;
		for (int h = 0; h < MMSP::length(s); h++) {
			int hindex = MMSP::index(s, h);
			for (int j = h + 1; j < MMSP::length(s); j++) {
				int jindex = MMSP::index(s, j);
				// latter dFdp_h and dFdp_j, so the inner loop can be over j>h instead of j≠h
				set(dFdp, hindex) += 0.5 * epsilon * epsilon * lap[jindex] + w * former(node)[jindex];
				set(dFdp, jindex) += 0.5 * epsilon * epsilon * lap[hindex] + w * former(node)[hindex];
			}
		}

		// compute time derivatives
		MMSP::sparse<T> dpdt;
		for (int h = 0; h < MMSP::length(s); h++) {
			int hindex = MMSP::index(s, h);
			for (int j = h + 1; j < MMSP::length(s); j++) {
				int jindex = MMSP::index(s, j);
				set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
				set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
			}
		}

		// compute latter values
		T sum = 0.0;
		for (int h = 0; h < MMSP::length(s); h++) {
			int index = MMSP::index(s, h);
			T value = former(node)[index] + dt * (2.0 / S) * dpdt[index]; // Extraneous factor of 2?
			if (value > 1.0) value = 1.0;
			if (value < 0.0) value = 0.0;
			if (value > 1.0e-8) MMSP::set(latter(node), index) = value;
			sum += latter(node)[index];
		}

		// project onto Gibbs simplex (enforce Σφ=1)
		T rsum = 0.0;
		if (fabs(sum) > 0.0) rsum = 1.0 / sum;
		for (int h = 0; h < MMSP::length(latter(node)); h++) {
			int index = MMSP::index(latter(node), h);
			MMSP::set(latter(node), index) *= rsum;
		}
	}
}

#endif
