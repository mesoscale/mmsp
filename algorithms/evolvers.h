// evolvers.hpp
// Definitions of methods to evolve MMSP grids
#ifndef _EVOLVERS_H_
#define _EVOLVERS_H_
void print_progress(const int step, const int steps, const int iterations);

template<int dim, typename T>
void interface_field_evolves(const T& dt, const T& width, const T& gamma, const T& epsilon, const T& w, const T& mu,
														 const MMSP::grid<dim,MMSP::sparse<T> >& former, MMSP::grid<dim,MMSP::sparse<T> >& latter, const int node);

#endif // _EVOLVERS_H_
