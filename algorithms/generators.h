// generators.hpp
// Definitions of methods to initialize MMSP grids
#ifndef _GENERATORS_H_
#define _GENERATORS_H_

template <int dim, typename T>
double radius(const MMSP::vector<T>& a, const MMSP::vector<T>& b);

template <int dim, typename T>
int radiussq(const MMSP::vector<T>& a, const MMSP::vector<T>& b);

int seeds_to_buffer(const std::vector<MMSP::vector<int> >& vp, int* &q);

void seeds_from_buffer(std::vector<MMSP::vector<int> >& vp, int* &q, const int& size);

void seedswap(std::vector<std::vector<MMSP::vector<int> > >& seeds);

template<int dim, typename T>
void exact_voronoi(MMSP::grid<dim,MMSP::sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds);

template<int dim, typename T>
void approximate_voronoi(MMSP::grid<dim,MMSP::sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds);

template<int dim, typename T>
void tessellate(MMSP::grid<dim,MMSP::sparse<T> >& grid, const int& nseeds);

template<int dim>
void seeds_from_poisson_process(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const int& nseeds, std::vector<std::vector<MMSP::vector<int> > >& seeds);

template<int dim>
void seeds_to_file(const int g0[dim], const int g1[dim], const std::vector<std::vector<MMSP::vector<int> > >& seeds, const char* filename);

template<int dim>
void seeds_from_file(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const char* seedfilename, std::vector<std::vector<MMSP::vector<int> > >& seeds);

template<int dim>
void honeycomb_seeds(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const int a, std::vector<std::vector<MMSP::vector<int> > >& seeds);

void print_progress(const int i, const int N);

#endif // _GENERATORS_H_
