// File contains code to initialize an MMSP simulation with a Voronoi tessellation
// of the simulation domain. Depends on Fast Marching code written supplied in part by
// Dr. Barb Cutler, Computer Science Department, Rensselaer Polytechnic Institute,
// and completed by Trevor Keller.

#ifndef _TESSELLATE_HPP_
#define _TESSELLATE_HPP_

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <ctime>
#include <limits>
#include <cassert>
#include "MersenneTwister.h"
#include "point.hpp"
#include "rdtsc.h"

// MMSP boundary conditions -- copied from MMSP.utility.hpp
enum {
    mirror    = 0,
    Neumann   = 1,
    periodic  = 2,
    parallel  = 3,
    Dirichlet = 4
};

template <int dim, typename T> double radius(const MMSP::vector<T>& a, const MMSP::vector<T>& b)
{
	double radius=0.;
	for (int d=0; d<dim; ++d) radius+=double(b[d]-a[d])*double(b[d]-a[d]);
	return sqrt(radius);
}

template <int dim, typename T> double radius(const MMSP::vector<T>& a, const Point<T>& b)
{
	double radius=0.;
	for (int d=0; d<dim; ++d) radius+=double(b[d]-a[d])*double(b[d]-a[d]);
	return sqrt(radius);
}

template <int dim, typename T> double radius(const Point<T>& a, const Point<T>& b)
{
	double radius=0.;
	for (int d=0; d<dim; ++d) radius+=double(b[d]-a[d])*double(b[d]-a[d]);
	return sqrt(radius);
}

double radius(int x0 = 0, int x1 = 0, int y0 = 0, int y1 = 0, int z0 = 0, int z1 = 0)
{
	return sqrt( pow(double(x0 - x1), 2) + pow(double(y0 - y1), 2) + pow(double(z0 - z1), 2) );
}

template <int dim, class T>
MMSP::vector<int> getPosition(const Point<T>& p)
{
	MMSP::vector<T> answer(dim, 0); // makes new vector of length dim, each entry is 0
	answer[0] = p.x;
	answer[1] = p.y;
	if (dim == 3) answer[2] = p.z;
	return answer;
}

class DistanceVoxel
{
public:
	DistanceVoxel() : x( -1 ), y( -1 ), z( 0 ), distance( std::numeric_limits<double>::max() ) {}
	// accessor
	int getX() const { return x; }
	int getY() const { return y; }
	int getZ() const { return z; }
	double getValue() const { return distance; }
	unsigned int getID() const { return id; }
	// modifier
	void setX( int _x ) { x = _x; }
	void setY( int _y ) { y = _y; }
	void setZ( int _z ) { z = _z; }
	void setValue( double v ) { distance = v; }
	void setID( unsigned int _i ) { id = _i; }
private:
	// REPRESENTATION
	int x;           // a distance voxel
	int y;           // knows its position in the
	int z;           // image and which
	unsigned int id; // seed voxel it's closest to
	double distance; // and how far away that is
};

#include "priority_queue.h"

namespace MMSP
{

#ifdef MPI_VERSION

#ifdef PHASEFIELD
// Voronoi tessellation for MMSP::Grid<dim,MMSP::sparse<T>>

template<int dim, typename T>
struct exact_voronoi_thread_para {
	MMSP::grid<dim,sparse<T> >* grid;
	std::vector<std::vector<Point<int> > >* seeds;
	std::set<unsigned int>* neigh;
	unsigned long nstart;
	unsigned long nend;
};

template<int dim, typename T>
void * exact_voronoi_threads_helper( void* s )
{
	exact_voronoi_thread_para<dim,T>* ss = ( exact_voronoi_thread_para<dim,T>* ) s ;

	for (unsigned long n=ss->nstart; n < ss->nend; ++n) {
		const MMSP::vector<int> x=position(*(ss->grid),n);
		double min_distance=std::numeric_limits<double>::max();
		int min_identity=-1;

		for (std::set<unsigned int>::const_iterator i=(*ss->neigh).begin(); i != (*ss->neigh).end(); i++) {
			int identity=-1;
			unsigned int rank=*i;
			for (unsigned int j=0; j<rank; j++) identity+=(*(ss->seeds))[j].size();
			for (unsigned int s=0; s<(*(ss->seeds))[rank].size(); ++s) {
				++identity;
				Point<int> seed=(*(ss->seeds))[rank][s];
				double distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check coordinates across periodic boundary
/*
				for (int d=0; d<dim; d++)
					check_boundary(seed[d], x0((*(ss->grid)),d), x1(*(ss->grid),d), b0(*(ss->grid),d), b1(*(ss->grid),d));
				if (seed==(*(ss->seeds))[rank][s]) continue;
				distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
*/
			}
		}
		set((*(ss->grid))(n), min_identity) = 1.;
	}

	pthread_exit(0);
	return NULL;
} // exact_voronoi

template<int dim, typename T>
unsigned long exact_voronoi_threads(MMSP::grid<dim,sparse<T> >& grid, std::vector<std::vector<Point<int> > >& seeds, const int& nthreads)
{
	// Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
	int id=MPI::COMM_WORLD.Get_rank();
	unsigned int np=MPI::COMM_WORLD.Get_size();

	// Determine neighborhood of seeds to scan
	// based on determination of n0, n1 in MMSP.grid.hpp
	std::set<unsigned int> neighbors;
	neighbors.insert(id);
	for (int d=0; d<dim; d++) {
		neighbors.insert(N0(grid,d));
		neighbors.insert(N1(grid,d));
	}
	for (int d=0; d<dim; d++) {
		int Nid=N0(grid,d);
		int pos[dim];
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
		Nid=N1(grid,d);
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
	}

	pthread_t* p_threads = new pthread_t[nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);
	exact_voronoi_thread_para<dim,T>* voronoi_para = new exact_voronoi_thread_para<dim,T>[nthreads];

	const unsigned long nincr = nodes(grid)/nthreads;
	unsigned long ns = 0;

	unsigned long timer=rdtsc();
	for (int i=0; i<nthreads; i++) {
		voronoi_para[i].nstart=ns;
		ns+=nincr;
		voronoi_para[i].nend=ns;

		voronoi_para[i].grid = &grid;
		voronoi_para[i].seeds = &seeds;
		voronoi_para[i].neigh = &neighbors;

		pthread_create(&p_threads[i], &attr, exact_voronoi_threads_helper<dim,T>, (void*) &voronoi_para[i] );
	}

	for (int i=0; i!= nthreads ; i++)
		pthread_join(p_threads[i], NULL);

	delete [] p_threads;
	delete [] voronoi_para;
	return rdtsc() - timer;
}

#else
// Voronoi tessellation for MMSP::grid<dim,T>

template<int dim, typename T>
struct exact_voronoi_thread_para {
	MMSP::grid<dim,T>* grid;
	std::vector<std::vector<Point<int> > >* seeds;
	std::set<unsigned int>* neigh;
	unsigned long nstart;
	unsigned long nend;
};

template<int dim, typename T>
void * exact_voronoi_threads_helper( void* s )
{
	exact_voronoi_thread_para<dim,T>* ss = ( exact_voronoi_thread_para<dim,T>* ) s ;

	for (unsigned long n=ss->nstart; n < ss->nend; ++n) {
		const MMSP::vector<int> x=position(*(ss->grid),n);
		double min_distance=std::numeric_limits<double>::max();
		int min_identity=-1;

		for (std::set<unsigned int>::const_iterator i=(*ss->neigh).begin(); i != (*ss->neigh).end(); i++) {
			int identity=-1;
			unsigned int rank=*i;
			for (unsigned int j=0; j<rank; j++) identity+=(*(ss->seeds))[j].size();
			for (unsigned int s=0; s<(*(ss->seeds))[rank].size(); ++s) {
				++identity;
				Point<int> seed=(*(ss->seeds))[rank][s];
				double distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check coordinates across periodic boundary
/*
				for (int d=0; d<dim; d++)
					check_boundary(seed[d], x0((*(ss->grid)),d), x1(*(ss->grid),d), b0(*(ss->grid),d), b1(*(ss->grid),d));
				if (seed==(*(ss->seeds))[rank][s]) continue;
				distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
*/
			}
		}
		(*(ss->grid))(n) = static_cast<T>(min_identity);
	}

	pthread_exit(0);
	return NULL;
} // exact_voronoi

template<int dim, typename T>
unsigned long exact_voronoi_threads(MMSP::grid<dim,T>& grid, std::vector<std::vector<Point<int> > >& seeds, const int& nthreads)
{
	// Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
	int id=MPI::COMM_WORLD.Get_rank();
	unsigned int np=MPI::COMM_WORLD.Get_size();

	// Determine neighborhood of seeds to scan
	// based on determination of n0, n1 in MMSP.grid.hpp
	std::set<unsigned int> neighbors;
	neighbors.insert(id);
	for (int d=0; d<dim; d++) {
		neighbors.insert(N0(grid,d));
		neighbors.insert(N1(grid,d));
	}
	for (int d=0; d<dim; d++) {
		int Nid=N0(grid,d);
		int pos[dim];
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
		Nid=N1(grid,d);
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
	}

	pthread_t* p_threads = new pthread_t[nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);
	exact_voronoi_thread_para<dim,T>* voronoi_para = new exact_voronoi_thread_para<dim,T>[nthreads];

	const unsigned long nincr = nodes(grid)/nthreads;
	unsigned long ns = 0;

	unsigned long timer = rdtsc();
	for (int i=0; i<nthreads; i++) {
		voronoi_para[i].nstart=ns;
		ns+=nincr;
		voronoi_para[i].nend=ns;

		voronoi_para[i].grid = &grid;
		voronoi_para[i].seeds = &seeds;
		voronoi_para[i].neigh = &neighbors;

		pthread_create(&p_threads[i], &attr, exact_voronoi_threads_helper<dim,T>, (void*) &voronoi_para[i] );
	}

	for (int i=0; i!= nthreads ; i++)
		pthread_join(p_threads[i], NULL);

	delete [] p_threads;
	delete [] voronoi_para;
	return rdtsc() - timer;
}

#endif

template<int dim, typename T>
unsigned long exact_voronoi(MMSP::grid<dim, sparse<T> >& grid, const std::vector<std::vector<Point<int> > >& seeds)
{
	// Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
	int id=MPI::COMM_WORLD.Get_rank();
	int np=MPI::COMM_WORLD.Get_size();

	// Determine neighborhood of seeds to scan
	// based on determination of n0, n1 in MMSP.grid.hpp
	std::set<unsigned int> neighbors;
	neighbors.insert(id);
	for (int d=0; d<dim; d++) {
		// Add first-nearest neighbors
		neighbors.insert(N0(grid, d));
		neighbors.insert(N1(grid, d));
	}
	for (int d=0; d<dim; d++) {
		// Add second-nearest neighbors
		int Nid=N0(grid,d);
		int pos[dim];
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
		Nid=N1(grid,d);
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
	}

	unsigned long timer=rdtsc();
	for (unsigned long n=0; n<nodes(grid); ++n) {
		const MMSP::vector<int> x=position(grid,n);
		double min_distance=std::numeric_limits<double>::max();
		int min_identity=-1;

		for (std::set<unsigned int>::const_iterator i=neighbors.begin(); i!=neighbors.end(); i++) {
			int identity=-1;
			unsigned int rank=*i;
			for (unsigned int j=0; j<rank; j++) identity+=seeds[j].size();
			for (unsigned int s=0; s<seeds[rank].size(); ++s) {
				++identity;
				Point<int> seed=seeds[rank][s];
				double distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check coordinates across periodic boundary
/*
				for (int d=0; d<dim; d++) check_boundary(seed[d], x0(grid,d), x1(grid,d), b0(grid,d), b1(grid,d));
				if (seed==seeds[rank][s]) continue;
				distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
*/
			}
		}
		set(grid(n), min_identity) = 1.;
	}
	return rdtsc() - timer;
} // exact_voronoi

template<int dim, typename T>
unsigned long exact_voronoi(MMSP::grid<dim,T>& grid, const std::vector<std::vector<Point<int> > >& seeds)
{
	// Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
	int id=MPI::COMM_WORLD.Get_rank();
	int np=MPI::COMM_WORLD.Get_size();

	// Determine neighborhood of seeds to scan
	// based on determination of n0, n1 in MMSP.grid.hpp
	std::set<unsigned int> neighbors;
	neighbors.insert(id);
	for (int d=0; d<dim; d++) {
		// Add first-nearest neighbors
		neighbors.insert(N0(grid, d));
		neighbors.insert(N1(grid, d));
	}
	for (int d=0; d<dim; d++) {
		// Add second-nearest neighbors
		int Nid=N0(grid,d);
		int pos[dim];
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
		Nid=N1(grid,d);
		for (int i=0; i<dim; i++)
			pos[i]=Nid/sp(grid,i);
		for (int i=0; i<dim; i++) {
			if (i==d) continue; // exclude 3rd-nearest
			int snpos[dim];
			for (int j=0; j<dim; j++)
				snpos[j]=pos[j];
			unsigned int snid=0;
			snpos[i] = (pos[i]-1+P1(grid,i))%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
			snid=0;
			snpos[i] = (pos[i]+1)%P1(grid,i);
			for (int j=0; j<dim; j++)
				snid+=sp(grid,j)*snpos[j];
			if (snid<np) neighbors.insert(snid);
		}
	}

	unsigned long timer=rdtsc();
	for (unsigned long n=0; n<nodes(grid); ++n) {
		const MMSP::vector<int> x=position(grid,n);
		double min_distance=std::numeric_limits<double>::max();
		T min_identity=-1;

		for (std::set<unsigned int>::const_iterator i=neighbors.begin(); i!=neighbors.end(); i++) {
			T identity=-1;
			unsigned int rank=*i;
			for (unsigned int j=0; j<rank; j++) identity+=seeds[j].size();
			for (unsigned int s=0; s<seeds[rank].size(); ++s) {
				++identity;
				Point<int> seed=seeds[rank][s];
				double distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check coordinates across periodic boundary
/*
				for (int d=0; d<dim; d++) check_boundary(seed[d], x0(grid,d), x1(grid,d), b0(grid,d), b1(grid,d));
				if (seed==seeds[rank][s]) continue;
				distance=radius<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
*/
			}
		}
		grid(n) = reinterpret_cast<int>(min_identity);
	}
	return rdtsc() - timer;
} // exact_voronoi
#endif

template<int dim>
void propagate_distance( const DistanceVoxel* core_voxel, MMSP::grid<dim, DistanceVoxel>& grid, DistanceVoxel_PriorityQueue& queue)
{
	if (dim == 2) {
		const int x = core_voxel->getX();
		const int y = core_voxel->getY();
		MMSP::vector<int> position = getPosition<dim, int>(Point<int>(x, y));
		const double core_distance = (grid(position)).getValue();
		const unsigned int core_id = (grid(position)).getID();
		// Loop over 9 neighboring Voxels
		for ( int i = x - 1; i <= x + 1; ++i ) {
			for ( int j = y - 1; j <= y + 1; ++j ) {
				if ( i == x && j == y) continue; // know thyself
				MMSP::vector<int> p = getPosition<dim, int>(Point<int>(i, j));
				// Take care of periodic boundary conditions
				for (int d = 0; d < dim; ++d) check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
				if ((p[0] < x0(grid, 0)) || (p[0] >= x1(grid, 0))) continue;
				if ((p[1] < x0(grid, 1)) || (p[1] >= x1(grid, 1))) continue;
				double distance = core_distance + radius(x, i, y, j); // using pre-check distance
				DistanceVoxel* Voxel = &( grid( p ) );
				if ( Voxel->getValue() > distance ) {
					Voxel->setValue( distance );
					Voxel->setID( core_id );
					if ( !queue.in_heap( Voxel ) ) queue.push( Voxel ); // Add new Voxel to heap
					else queue.update_position( Voxel ); // Update position of Voxel in the heap based on new value
				}
			}
		}
	} else if (dim == 3) {
		const int x = core_voxel->getX();
		const int y = core_voxel->getY();
		const int z = core_voxel->getZ();
		MMSP::vector<int> position = getPosition<dim, int>(Point<int>(x, y, z));
		const double core_distance = (grid(position)).getValue();
		const unsigned int core_id = (grid(position)).getID();
		// Loop over 27 neighboring Voxels
		for ( int i = x - 1; i <= x + 1; ++i ) {
			for ( int j = y - 1; j <= y + 1; ++j ) {
				for (int k = z - 1; k <= z + 1; ++k) {
					if ( i == x && j == y && k == z) continue; // know thyself
					MMSP::vector<int> p = getPosition<dim, int>(Point<int>(i, j, k));
					// Take care of periodic boundary conditions
					for (int d = 0; d < dim; ++d) check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
					if ((p[0] < x0(grid, 0)) || (p[0] >= x1(grid, 0))) continue;
					if ((p[1] < x0(grid, 1)) || (p[1] >= x1(grid, 1))) continue;
					if ((p[2] < x0(grid, 2)) || (p[2] >= x1(grid, 2))) continue;
					double distance = core_distance + radius(x, i, y, j, z, k); // using pre-check distance
					DistanceVoxel* Voxel = &( grid(p) );
					if ( Voxel->getValue() > distance ) {
						Voxel->setValue( distance );
						Voxel->setID( core_id );
						if ( !queue.in_heap( Voxel ) ) queue.push( Voxel ); // Add new Voxel to heap
						else queue.update_position( Voxel ); // Update position of Voxel in the heap based on new value
					}
				}
			}
		}
	}
} // propagate distance

// Given a populated vector<Point<int>> and empty array of T[size], populate the array
int seeds_to_buffer(const std::vector<Point<int> >& vp, int*& q)
{
	int size = 0;
	// q should already point to an array of T[size]
	if (q == NULL) {
		std::cerr << "\nError in seeds_to_buffer: send_buffer not initialized." << std::endl;
		exit(1);
	}
	int* p = q;
	for (unsigned int i = 0; i < vp.size(); ++i) {
		size += 3;
		*p = vp[i].x;
		++p;
		*p = vp[i].y;
		++p;
		*p = vp[i].z;
		++p;
	}
	return size; // number of Points created
}

// Given an empty vector<Point<int>> and a populated array of T[size], populate the vector
void seeds_from_buffer(std::vector<Point<int> >& vp, int*& q, const int& size)
{
	for (int* p = q + 2; p < q + size; p += 3) {
		int x = *(p - 2);
		int y = *(p - 1);
		int z = *p;
		vp.push_back(Point<int>(x, y, z));
	}
}


template<int dim, typename T>
unsigned long approximate_voronoi(MMSP::grid<dim, sparse<T> >& grid, const std::vector<std::vector<Point<int> > >& seeds)
{
	// Implements a fast marching algorithm to generate the distance map
	// Based on code written by Barb Cutler, RPI Comp. Sci. Dept., for CSCI-1200.
	int id = 0;
	#ifdef MPI_VERSION
	id = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	// Perform the tessellation, using fast-marching fanciness
	unsigned long timer=0;
	if (dim == 2) {
		int min[2];
		int max[2];
		for (int i = 0; i < dim; ++i) {
			min[i] = x0(grid, i);
			max[i] = x1(grid, i);
		}
		// Constructor requires global limits of the grid. Magical!
		MMSP::grid<2, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1));
		DistanceVoxel value;
		value.setValue( std::numeric_limits<double>::max() );
		value.setID( 0 );

		// Initialize distance_grid with very-large distances
		for (int i = 0; i < nodes(distance_grid); ++i) {
			MMSP::vector<int> pos = grid.position(i);
			for (int d = 0; d < dim; ++d) {
				assert(pos[d] >= x0(distance_grid, d));
				assert(pos[d] < x1(distance_grid, d));
			}
			value.setX( pos[0] );
			value.setY( pos[1] );
			value.setZ( pos[2] );
			distance_grid(i) = value;
		}

		// create the voxel Heap
		DistanceVoxel_PriorityQueue queue;

		int nseeds = 0;
		for ( int i = 0; i < id; ++i ) nseeds += seeds[i].size();

		// Enqueue this node's seeds
		for ( int i = 0; i < seeds[id].size(); ++i ) {
			MMSP::vector<int> pos = getPosition<dim, int>(seeds[id][i]);
			DistanceVoxel* p = &( distance_grid(pos) );
			for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
			p->setValue( 0. );
			p->setID( nseeds + i );
			// Propagate distance from each seed to its neighbors. Start adding to the Heap.
			propagate_distance( p, distance_grid, queue );
		}

		// Fast-march the local seeds
		unsigned long qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		#ifdef MPI_VERSION
		// Copy ghost voxels from adjacent ranks
		ghostswap(distance_grid);

		// Propagate ghost distances
		int n0[3];
		int n1[3];
		for (int d = 0; d < dim; ++d) {
			n0[d] = N0(grid, d); // Identify processor below in this dimension
			n1[d] = N1(grid, d); // Identify processor above in this dimension
			// Parallel boundary condition
			if (n0[d] != id) {
				// Ghosts from below
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				int scan_axis = (d == 0) ? 1 : 0;
				for (x[scan_axis] = x0(grid, scan_axis); x[scan_axis] < x1(grid, scan_axis); ++x[scan_axis]) {
					MMSP::vector<int> g = x;
					--g[d]; // ghost point next to x
					DistanceVoxel* px = &( distance_grid(x) );
					DistanceVoxel* pg = &( distance_grid(g) );
					double r = pg->getValue() + 1;
					if (r < px->getValue() ) {
						px->setValue( pg->getValue() + 1 );
						px->setID( pg->getID() );
						propagate_distance( px, distance_grid, queue );
					}
					++x[scan_axis];
				}
			}
			if (n1[d] != id) {
				// Ghosts from above
				MMSP::vector<int> x = position(distance_grid, nodes(distance_grid) - 1); // real point at top corner
				int scan_axis = (d == 0) ? 1 : 0;
				for (x[scan_axis] = x0(grid, scan_axis); x[scan_axis] < x1(grid, scan_axis); ++x[scan_axis]) {
					MMSP::vector<int> g = x;
					++g[d]; // ghost point next to x
					DistanceVoxel* px = &( distance_grid(x) );
					DistanceVoxel* pg = &( distance_grid(g) );
					double r = pg->getValue() + 1;
					if (r < px->getValue() ) {
						px->setValue( pg->getValue() + 1 );
						px->setID( pg->getID() );
						propagate_distance( px, distance_grid, queue );
					}
				}
			}
		}
		#endif

		// Fast-march the ghost points
		qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		// Copy result from distance_grid to phase-field grid
		for (int i = 0; i < nodes(distance_grid); ++i)
			MMSP::set(grid(i),distance_grid(i).getID()) = 1.;
	} else if (dim == 3) {
		int min[3];
		int max[3];
		for (int i = 0; i < dim; ++i) {
			min[i] = x0(grid, i);
			max[i] = x1(grid, i);
		}
		// Constructor requires global limits of the grid. Magical!
		MMSP::grid<3, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1), g0(grid, 2), g1(grid, 2));
		DistanceVoxel value;
		value.setValue( std::numeric_limits<double>::max() );
		value.setID( 0 );
		for (int i = 0; i < nodes(distance_grid); ++i) {
			MMSP::vector<int> pos = grid.position(i);
			for (int d = 0; d < dim; ++d) {
				assert(pos[d] >= x0(distance_grid, d));
				assert(pos[d] < x1(distance_grid, d));
			}
			value.setX( pos[0] );
			value.setY( pos[1] );
			value.setZ( pos[2] );
			distance_grid(i) = value;
		}
		// create the voxel Heap
		DistanceVoxel_PriorityQueue queue;

		int nseeds = 0;
		for ( int i = 0; i < id; ++i ) nseeds += seeds[i].size();

		// Start queue with this node's seeds
		for ( int i = 0; i < seeds[id].size(); ++i ) {
			MMSP::vector<int> pos = getPosition<dim, int>(seeds[id][i]);
			DistanceVoxel* p = &( distance_grid(pos) );
			for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
			p->setValue( 0. );
			p->setID( nseeds + i );
			// Propagate distance from each seed to its neighbors. Start adding to the Heap.
			propagate_distance( p, distance_grid, queue );
		}

		// Fast-march the local seeds
		unsigned long qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		#ifdef MPI_VERSION
		// Copy ghost voxels from adjacent ranks
		ghostswap(distance_grid);

		// Propagate ghost distances
		int n0[3];
		int n1[3];
		for (int d = 0; d < dim; ++d) {
			n0[d] = N0(grid, d); // Identify processor below in this dimension
			n1[d] = N1(grid, d); // Identify processor above in this dimension
			// Parallel boundary condition
			if (n0[d] != id) {
				// Ghosts from below
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				MMSP::vector<int> scan_axes(dim - 1, 0);
				for (int e = 0; e < dim - 1; ++e) scan_axes[e] = (d == 0 || e == d) ? e + 1 : e;
				for (x[scan_axes[0]] = x0(grid, scan_axes[0]); x[scan_axes[0]] < x1(grid, scan_axes[0]); ++x[scan_axes[0]]) {
					for (x[scan_axes[1]] = x0(grid, scan_axes[1]); x[scan_axes[1]] < x1(grid, scan_axes[1]); ++x[scan_axes[1]]) {
						MMSP::vector<int> g = x;
						--g[d]; // ghost point next to x
						DistanceVoxel* px = &( distance_grid(x) );
						DistanceVoxel* pg = &( distance_grid(g) );
						double r = pg->getValue() + 1;
						if (r < px->getValue() ) {
							px->setValue( pg->getValue() + 1 );
							px->setID( pg->getID() );
							propagate_distance( px, distance_grid, queue );
						}
					}
				}
			}
			if (n1[d] != id) {
				// Ghosts from above
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				MMSP::vector<int> scan_axes(dim - 1, 0);
				for (int e = 0; e < dim - 1; ++e) scan_axes[e] = (d == 0 || e == d) ? e + 1 : e;
				for (x[scan_axes[0]] = x0(grid, scan_axes[0]); x[scan_axes[0]] < x1(grid, scan_axes[0]); ++x[scan_axes[0]]) {
					for (x[scan_axes[1]] = x0(grid, scan_axes[1]); x[scan_axes[1]] < x1(grid, scan_axes[1]); ++x[scan_axes[1]]) {
						MMSP::vector<int> g = x;
						++g[d]; // ghost point next to x
						DistanceVoxel* px = &( distance_grid(x) );
						DistanceVoxel* pg = &( distance_grid(g) );
						double r = pg->getValue() + 1;
						if (r < px->getValue() ) {
							px->setValue( pg->getValue() + 1 );
							px->setID( pg->getID() );
							propagate_distance( px, distance_grid, queue );
						}
					}
				}
			}
		}
		#endif

		// Fast-march the ghost points
		qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		// Copy result from distance_grid to phase-field grid
		for (int i = 0; i < nodes(distance_grid); ++i)
			MMSP::set(grid(i),distance_grid(i).getID()) = 1.;
	} else {
		std::cerr << "Error: Invalid dimension (" << dim << ") in tessellation." << std::endl;
		std::exit(1);
	}
	return timer;
} // approximate_voronoi

template<int dim, typename T>
unsigned long approximate_voronoi(MMSP::grid<dim,T>& grid, const std::vector<std::vector<Point<int> > >& seeds)
{
	// Implements a fast marching algorithm to generate the distance map
	// Based on code written by Barb Cutler, RPI Comp. Sci. Dept., for CSCI-1200.
	int id = 0;
	#ifdef MPI_VERSION
	id = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	// Perform the tessellation, using fast-marching fanciness
	unsigned long timer=0;
	if (dim == 2) {
		int min[2];
		int max[2];
		for (int i = 0; i < dim; ++i) {
			min[i] = x0(grid, i);
			max[i] = x1(grid, i);
		}
		// Constructor requires global limits of the grid. Magical!
		MMSP::grid<2, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1));
		DistanceVoxel value;
		value.setValue( std::numeric_limits<double>::max() );
		value.setID( 0 );

		// Initialize distance_grid with very-large distances
		for (int i = 0; i < nodes(distance_grid); ++i) {
			MMSP::vector<int> pos = grid.position(i);
			for (int d = 0; d < dim; ++d) {
				assert(pos[d] >= x0(distance_grid, d));
				assert(pos[d] < x1(distance_grid, d));
			}
			value.setX( pos[0] );
			value.setY( pos[1] );
			value.setZ( pos[2] );
			distance_grid(i) = value;
		}

		// create the voxel Heap
		DistanceVoxel_PriorityQueue queue;

		int nseeds = 0;
		for ( int i = 0; i < id; ++i ) nseeds += seeds[i].size();

		// Enqueue this node's seeds
		for ( int i = 0; i < seeds[id].size(); ++i ) {
			MMSP::vector<int> pos = getPosition<dim, int>(seeds[id][i]);
			DistanceVoxel* p = &( distance_grid(pos) );
			for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
			p->setValue( 0. );
			p->setID( nseeds + i );
			// Propagate distance from each seed to its neighbors. Start adding to the Heap.
			propagate_distance( p, distance_grid, queue );
		}

		// Fast-march the local seeds
		unsigned long qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		#ifdef MPI_VERSION
		// Copy ghost voxels from adjacent ranks
		ghostswap(distance_grid);

		// Propagate ghost distances
		int n0[3];
		int n1[3];
		for (int d = 0; d < dim; ++d) {
			n0[d] = N0(grid, d); // Identify processor below in this dimension
			n1[d] = N1(grid, d); // Identify processor above in this dimension
			// Parallel boundary condition
			if (n0[d] != id) {
				// Ghosts from below
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				int scan_axis = (d == 0) ? 1 : 0;
				for (x[scan_axis] = x0(grid, scan_axis); x[scan_axis] < x1(grid, scan_axis); ++x[scan_axis]) {
					MMSP::vector<int> g = x;
					--g[d]; // ghost point next to x
					DistanceVoxel* px = &( distance_grid(x) );
					DistanceVoxel* pg = &( distance_grid(g) );
					double r = pg->getValue() + 1;
					if (r < px->getValue() ) {
						px->setValue( pg->getValue() + 1 );
						px->setID( pg->getID() );
						propagate_distance( px, distance_grid, queue );
					}
					++x[scan_axis];
				}
			}
			if (n1[d] != id) {
				// Ghosts from above
				MMSP::vector<int> x = position(distance_grid, nodes(distance_grid) - 1); // real point at top corner
				int scan_axis = (d == 0) ? 1 : 0;
				for (x[scan_axis] = x0(grid, scan_axis); x[scan_axis] < x1(grid, scan_axis); ++x[scan_axis]) {
					MMSP::vector<int> g = x;
					++g[d]; // ghost point next to x
					DistanceVoxel* px = &( distance_grid(x) );
					DistanceVoxel* pg = &( distance_grid(g) );
					double r = pg->getValue() + 1;
					if (r < px->getValue() ) {
						px->setValue( pg->getValue() + 1 );
						px->setID( pg->getID() );
						propagate_distance( px, distance_grid, queue );
					}
				}
			}
		}
		#endif

		// Fast-march the ghost points
		qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		// Copy result from distance_grid to phase-field grid
		for (int i = 0; i < nodes(distance_grid); ++i)
			grid(i) = distance_grid(i).getID();

	} else if (dim == 3) {
		int min[3];
		int max[3];
		for (int i = 0; i < dim; ++i) {
			min[i] = x0(grid, i);
			max[i] = x1(grid, i);
		}
		// Constructor requires global limits of the grid. Magical!
		MMSP::grid<3, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1), g0(grid, 2), g1(grid, 2));
		DistanceVoxel value;
		value.setValue( std::numeric_limits<double>::max() );
		value.setID( 0 );
		for (int i = 0; i < nodes(distance_grid); ++i) {
			MMSP::vector<int> pos = grid.position(i);
			for (int d = 0; d < dim; ++d) {
				assert(pos[d] >= x0(distance_grid, d));
				assert(pos[d] < x1(distance_grid, d));
			}
			value.setX( pos[0] );
			value.setY( pos[1] );
			value.setZ( pos[2] );
			distance_grid(i) = value;
		}
		// create the voxel Heap
		DistanceVoxel_PriorityQueue queue;

		int nseeds = 0;
		for ( int i = 0; i < id; ++i ) nseeds += seeds[i].size();

		// Start queue with this node's seeds
		for ( int i = 0; i < seeds[id].size(); ++i ) {
			MMSP::vector<int> pos = getPosition<dim, int>(seeds[id][i]);
			DistanceVoxel* p = &( distance_grid(pos) );
			for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
			p->setValue( 0. );
			p->setID( nseeds + i );
			// Propagate distance from each seed to its neighbors. Start adding to the Heap.
			propagate_distance( p, distance_grid, queue );
		}

		// Fast-march the local seeds
		unsigned long qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		#ifdef MPI_VERSION
		// Copy ghost voxels from adjacent ranks
		ghostswap(distance_grid);

		// Propagate ghost distances
		int n0[3];
		int n1[3];
		for (int d = 0; d < dim; ++d) {
			n0[d] = N0(grid, d); // Identify processor below in this dimension
			n1[d] = N1(grid, d); // Identify processor above in this dimension
			// Parallel boundary condition
			if (n0[d] != id) {
				// Ghosts from below
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				MMSP::vector<int> scan_axes(dim - 1, 0);
				for (int e = 0; e < dim - 1; ++e) scan_axes[e] = (d == 0 || e == d) ? e + 1 : e;
				for (x[scan_axes[0]] = x0(grid, scan_axes[0]); x[scan_axes[0]] < x1(grid, scan_axes[0]); ++x[scan_axes[0]]) {
					for (x[scan_axes[1]] = x0(grid, scan_axes[1]); x[scan_axes[1]] < x1(grid, scan_axes[1]); ++x[scan_axes[1]]) {
						MMSP::vector<int> g = x;
						--g[d]; // ghost point next to x
						DistanceVoxel* px = &( distance_grid(x) );
						DistanceVoxel* pg = &( distance_grid(g) );
						double r = pg->getValue() + 1;
						if (r < px->getValue() ) {
							px->setValue( pg->getValue() + 1 );
							px->setID( pg->getID() );
							propagate_distance( px, distance_grid, queue );
						}
					}
				}
			}
			if (n1[d] != id) {
				// Ghosts from above
				MMSP::vector<int> x = position(distance_grid, 0); // real point at lower corner
				MMSP::vector<int> scan_axes(dim - 1, 0);
				for (int e = 0; e < dim - 1; ++e) scan_axes[e] = (d == 0 || e == d) ? e + 1 : e;
				for (x[scan_axes[0]] = x0(grid, scan_axes[0]); x[scan_axes[0]] < x1(grid, scan_axes[0]); ++x[scan_axes[0]]) {
					for (x[scan_axes[1]] = x0(grid, scan_axes[1]); x[scan_axes[1]] < x1(grid, scan_axes[1]); ++x[scan_axes[1]]) {
						MMSP::vector<int> g = x;
						++g[d]; // ghost point next to x
						DistanceVoxel* px = &( distance_grid(x) );
						DistanceVoxel* pg = &( distance_grid(g) );
						double r = pg->getValue() + 1;
						if (r < px->getValue() ) {
							px->setValue( pg->getValue() + 1 );
							px->setID( pg->getID() );
							propagate_distance( px, distance_grid, queue );
						}
					}
				}
			}
		}
		#endif

		// Fast-march the ghost points
		qtimer=rdtsc();
		while ( !queue.empty() ) {
			const DistanceVoxel* p = queue.top();
			queue.pop();
			propagate_distance( p, distance_grid, queue );
		}
		timer+=rdtsc()-qtimer;

		// Copy result from distance_grid to phase-field grid
		for (int i = 0; i < nodes(distance_grid); ++i)
			grid(i) = distance_grid(i).getID();

	} else {
		std::cerr << "Error: Invalid dimension (" << dim << ") in tessellation." << std::endl;
		std::exit(1);
	}
	return timer;
} // approximate_voronoi


template<int dim, typename T>
unsigned long tessellate(MMSP::grid<dim,T>& grid, const int& nseeds, const int& nthreads)
{
	int id=0;
	unsigned int np=1;
	#ifdef MPI_VERSION
	id=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	#endif
	unsigned long int pseudorand_seed = time(NULL);
	#ifndef SILENT
	if (id == 0) std::cout << "Master seed is " << std::setw(10) << std::right << pseudorand_seed << ". <---- Record this value!" << std::endl;
	#endif
	#ifdef MPI_VERSION
	pseudorand_seed = pseudorand_seed / (id + 1);
	#endif
	MTRand pseudorand_number( pseudorand_seed );
	std::vector<Point<int> > local_seeds; // blank for now
	std::vector<std::vector<Point<int> > > seeds;
	while (seeds.size() <= np) seeds.push_back(local_seeds); // avoid a segfault

	// Generate the seeds
	if (dim == 2) {
		int x = 0, y = 0;
		//if (id == 0) seed_points << "x,y\n";
		for (int i = 0; i < nseeds; ++i) {
			x = x0(grid, 0) + pseudorand_number.randInt( x1(grid, 0) - x0(grid, 0) - 1 );
			y = x0(grid, 1) + pseudorand_number.randInt( x1(grid, 1) - x0(grid, 1) - 1 );
			bool dupe = false;
			for (unsigned int j = 0; j < seeds[id].size(); ++j) {
				// No duplicates!
				if ((seeds[id][j].x == x) && (seeds[id][j].y == y)) {
					--i;
					dupe = true;
				}
				if (dupe) break; // stop scanning other seeds
			}
			if (dupe) continue; // don't add this seed
			local_seeds.push_back( Point<int>(x, y, 0) );
		}
	} else if (dim == 3) {
		int x = 0, y = 0, z = 0;
		//if (id == 0) seed_points << "x,y,z\n";
		for (int i = 0; i < nseeds; ++i) {
			x = x0(grid, 0) + pseudorand_number.randInt( x1(grid, 0) - x0(grid, 0) - 1 );
			y = x0(grid, 1) + pseudorand_number.randInt( x1(grid, 1) - x0(grid, 1) - 1 );
			z = x0(grid, 2) + pseudorand_number.randInt( x1(grid, 2) - x0(grid, 2) - 1 );
			bool dupe = false;
			for (unsigned int j = 0; j < seeds[id].size(); ++j) {
				// No duplicates!
				if ((seeds[id][j].x == x) && (seeds[id][j].y == y) && (seeds[id][j].z == z)) {
					--i;
					dupe = true;
				}
				if (dupe) break; // stop scanning other seeds
			}
			if (dupe) continue; // don't add this seed
			local_seeds.push_back( Point<int>(x, y, z) );
		}
	} else {
		std::cerr << "Error: Invalid dimension (" << dim << ") in tessellation." << std::endl;
		std::exit(1);
	}


	#ifndef MPI_VERSION
	seeds[id].insert(seeds[id].end(), local_seeds.begin(), local_seeds.end());
	#else
	// Exchange seeds between all processors
	// Gather  of seeds per processor
	int send_size=(local_seeds.size()) * (sizeof(Point<int>) / sizeof(int)); // number of integers
	int* send_buffer = new int[send_size]; // number of integers: should be 3*seeds[id].size()
	send_size = MMSP::seeds_to_buffer(local_seeds, send_buffer);
	int* seed_sizes = new int[np];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgather(&send_size, 1, MPI_INT, seed_sizes, 1, MPI_INT);
	int total_size=0;
	for (unsigned int i=0; i<np; ++i) total_size+=seed_sizes[i];
	int* offsets = new int[np];
	offsets[0]=0;
	for (unsigned int i=1; i<np; ++i) offsets[i]=seed_sizes[i-1]+offsets[i-1];
	int* seed_block = new int[total_size];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgatherv(send_buffer, send_size, MPI_INT, seed_block, seed_sizes, offsets, MPI_INT);
	delete [] send_buffer;
	send_buffer=NULL;

	for (unsigned int i=0; i<np; ++i) {
		int* p=seed_block+offsets[i];
		MMSP::seeds_from_buffer(seeds[i], p, seed_sizes[i]);
	}
	delete [] seed_sizes;
	seed_sizes=NULL;
	delete [] offsets;
	offsets=NULL;
	delete [] seed_block;
	seed_block=NULL;
	int vote=1;
	int total_procs=0;
	MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
	#ifndef SILENT
	if (id==0) std::cout<<"Synchronized seeds on "<<total_procs<<" ranks."<<std::endl;
	#endif
	#endif

	// Perform the actual tessellation
	unsigned long timer=0;
	#ifndef MPI_VERSION
	timer=approximate_voronoi<dim,T>(grid, seeds);
	#else
	timer=exact_voronoi_threads<dim,T>(grid, seeds, nthreads);
	MPI::COMM_WORLD.Barrier();
	total_procs=0;
	MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
	#ifndef SILENT
	if (id==0) std::cout<<"Tessellated the domain on "<<total_procs<<" ranks."<<std::endl;
	#endif
	#endif
	return timer;
} // tessellate

template<int dim, typename T>
unsigned long tessellate(MMSP::grid<dim, MMSP::sparse<T> >& grid, const int& nseeds, const int& nthreads)
{
	int id=0;
	int np=1;
	#ifdef MPI_VERSION
	id=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();
	#endif
	unsigned long int pseudorand_seed = time(NULL);
	#ifndef SILENT
	if (id == 0) std::cout << "Master seed is " << std::setw(10) << std::right << pseudorand_seed << ". <---- Record this value!" << std::endl;
	#endif
	#ifdef MPI_VERSION
	pseudorand_seed = pseudorand_seed / (id + 1);
	#endif
	MTRand pseudorand_number( pseudorand_seed );
	std::vector<Point<int> > local_seeds; // blank for now
	std::vector<std::vector<Point<int> > > seeds;
	while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault

	// Generate the seeds
	if (dim == 2) {
		int x = 0, y = 0;
		//if (id == 0) seed_points << "x,y\n";
		for (int i = 0; i < nseeds; ++i) {
			x = x0(grid, 0) + pseudorand_number.randInt( x1(grid, 0) - x0(grid, 0) - 1 );
			y = x0(grid, 1) + pseudorand_number.randInt( x1(grid, 1) - x0(grid, 1) - 1 );
			bool dupe = false;
			for (unsigned int j = 0; j < seeds[id].size(); ++j) {
				// No duplicates!
				if ((seeds[id][j].x == x) && (seeds[id][j].y == y)) {
					--i;
					dupe = true;
				}
				if (dupe) break; // stop scanning other seeds
			}
			if (dupe) continue; // don't add this seed
			local_seeds.push_back( Point<int>(x, y, 0) );
		}
	} else if (dim == 3) {
		int x = 0, y = 0, z = 0;
		//if (id == 0) seed_points << "x,y,z\n";
		for (int i = 0; i < nseeds; ++i) {
			x = x0(grid, 0) + pseudorand_number.randInt( x1(grid, 0) - x0(grid, 0) - 1 );
			y = x0(grid, 1) + pseudorand_number.randInt( x1(grid, 1) - x0(grid, 1) - 1 );
			z = x0(grid, 2) + pseudorand_number.randInt( x1(grid, 2) - x0(grid, 2) - 1 );
			bool dupe = false;
			for (unsigned int j = 0; j < seeds[id].size(); ++j) {
				// No duplicates!
				if ((seeds[id][j].x == x) && (seeds[id][j].y == y) && (seeds[id][j].z == z)) {
					--i;
					dupe = true;
				}
				if (dupe) break; // stop scanning other seeds
			}
			if (dupe) continue; // don't add this seed
			local_seeds.push_back( Point<int>(x, y, z) );
		}
	} else {
		std::cerr << "Error: Invalid dimension (" << dim << ") in tessellation." << std::endl;
		std::exit(1);
	}


	#ifndef MPI_VERSION
	seeds[id].insert(seeds[id].end(), local_seeds.begin(), local_seeds.end());
	#else
	// Exchange seeds between all processors
	// Gather  of seeds per processor
	int send_size=(local_seeds.size()) * (sizeof(Point<int>) / sizeof(int)); // number of integers
	int* send_buffer = new int[send_size]; // number of integers: should be 3*seeds[id].size()
	send_size = MMSP::seeds_to_buffer(local_seeds, send_buffer);
	int* seed_sizes = new int[np];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgather(&send_size, 1, MPI_INT, seed_sizes, 1, MPI_INT);
	int total_size=0;
	for (int i=0; i<np; ++i) total_size+=seed_sizes[i];
	int* offsets = new int[np];
	offsets[0]=0;
	for (int i=1; i<np; ++i) offsets[i]=seed_sizes[i-1]+offsets[i-1];
	int* seed_block = new int[total_size];
	MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allgatherv(send_buffer, send_size, MPI_INT, seed_block, seed_sizes, offsets, MPI_INT);
	delete [] send_buffer;
	send_buffer=NULL;

	for (int i=0; i<np; ++i) {
		int* p=seed_block+offsets[i];
		MMSP::seeds_from_buffer(seeds[i], p, seed_sizes[i]);
	}
	delete [] seed_sizes;
	seed_sizes=NULL;
	delete [] offsets;
	offsets=NULL;
	delete [] seed_block;
	seed_block=NULL;
	int vote=1;
	int total_procs=0;
	MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
	#ifndef SILENT
	if (id==0) std::cout<<"Synchronized seeds on "<<total_procs<<" ranks."<<std::endl;
	#endif
	#endif

	// Perform the actual tessellation
	unsigned long timer=0;
	#ifndef MPI_VERSION
	timer=approximate_voronoi<dim,T>(grid, seeds);
	#else
	timer=exact_voronoi_threads<dim,T>(grid, seeds, nthreads);
	MPI::COMM_WORLD.Barrier();
	total_procs=0;
	MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
	#ifndef SILENT
	if (id==0) std::cout<<"Tessellated the domain on "<<total_procs<<" ranks."<<std::endl;
	#endif
	#endif
	return timer;
} // tessellate

} // namespace
#endif

// Formatted using astyle:
//  astyle --style=linux --indent-col1-comments --indent=tab --indent-preprocessor --pad-header --align-pointer=type --keep-one-line-blocks --suffix=none

