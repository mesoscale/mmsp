// generators.cpp
// Implementations of methods to initialize MMSP grids
#ifndef _GENERATORS_CPP_
#define _GENERATORS_CPP_
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <limits>
#include <cassert>
#include "MersenneTwister.h"
#include "distancevoxel.h"
#include "priority_queue.h"

#include "generators.h"

// MMSP boundary conditions -- copied from MMSP.utility.hpp
enum {
  mirror    = 0,
  Neumann   = 1,
  periodic  = 2,
  parallel  = 3,
  Dirichlet = 4
};

template <int dim, typename T> double radius(const MMSP::vector<T>& a, const MMSP::vector<T>& b) {
  double radius=0.0;
	for (int d=0; d<dim; ++d) radius+=std::pow(double(b[d]-a[d]), 2.0);
  return sqrt(radius);
}

template <int dim, typename T> int radiussq(const MMSP::vector<T>& a, const MMSP::vector<T>& b) {
  int radiussq=0;
	for (int d=0; d<dim; ++d) radiussq+=(b[d]-a[d])*(b[d]-a[d]);
  return radiussq;
}

// Given a populated vector<MMSP::vector<int>> and empty array of T[size], populate the array
int seeds_to_buffer(const std::vector<MMSP::vector<int> >& vp, int* &q) {
  int size = 0;
  // q should already point to an array of T[size]
  if (q == NULL) {
    std::cerr << "\nError in seeds_to_buffer: send_buffer not initialized." << std::endl;
    exit(1);
  }
  int* p = q;
  for (unsigned int i=0; i<vp.size(); ++i) {
  	for (int d=0; d<3; ++d) {
    	++size;
    	*p = vp[i][d];
    	++p;
	  }
	}
  return size; // number of seeds created
}

// Given an empty vector<MMSP::vector<int>> and a populated array of T[size], populate the vector
void seeds_from_buffer(std::vector<MMSP::vector<int> >& vp, int* &q, const int& size) {
  if (q == NULL) {
    std::cerr << "\nError in seeds_from_buffer: recv_buffer not initialized." << std::endl;
    exit(1);
  }
  MMSP::vector<int> v(3);
  for (int* p = q + 2; p < q + size; p += 3) {
    v[0] = *(p - 2);
    v[1] = *(p - 1);
    v[2] = *p;
    vp.push_back(v);
  }
}

MMSP::vector<int> getPosition(const DistanceVoxel& dv) {
	MMSP::vector<int> x(3);
	for (int d=0; d<3; ++d) x[d]=dv[d];
	return x;
}


template<int dim, typename T>
void exact_voronoi(MMSP::grid<dim, MMSP::sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds) {
	int id=0;
	#ifdef MPI_VERSION
	id=MPI::COMM_WORLD.Get_rank();
	#endif
  // Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
  // seeds must contain every seed from every rank.
  #ifdef DEBUG
  unsigned long tstart=time(NULL);
  #endif
	if (id==0) print_progress(0,nodes(grid));
	for (int n=0; n<nodes(grid); ++n) {
		const MMSP::vector<int> x=position(grid,n);
		int min_distance=std::numeric_limits<int>::max();
		int identity=-1, min_identity=identity;

		for (std::vector<std::vector<MMSP::vector<int> > >::const_iterator rank=seeds.begin(); rank!=seeds.end(); ++rank) {
			for (std::vector<MMSP::vector<int> >::const_iterator s=rank->begin(); s!=rank->end(); ++s) {
				++identity;
				MMSP::vector<int> seed=*s;
				int distance=radiussq<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check position across parallel and periodic boundaries to ensure continuity across node boundaries
				for (int d=0; d<dim; d++) MMSP::check_boundary(seed[d], x0(grid,d), x1(grid,d), b0(grid,d), b1(grid,d));
				if (seed==*s) continue;
				distance=radiussq<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
			}
		}
		if (id==0) print_progress(n+1,nodes(grid));
		MMSP::set(grid(n), min_identity) = 1.;
	}
	#ifdef DEBUG
	#ifdef MPI_VERSION
	if (id==0)
	#endif
	std::cout<<"Completed exact tessellation ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	#endif
} // exact_voronoi

template<int dim>
void propagate_distance( const DistanceVoxel* core_voxel, MMSP::grid<dim, DistanceVoxel>& grid, DistanceVoxel_PriorityQueue& queue)
{
  if (dim == 2) {
    MMSP::vector<int> position = getPosition(*core_voxel);
    const int core_distance = (grid(position)).getValue();
    const unsigned int core_id = (grid(position)).getID();
    MMSP::vector<int> x(position);
    // Loop over 9 neighboring Voxels
    for (x[0]=position[0]-1; x[0]<position[0]+2; ++x[0]) {
    	for (x[1]=position[1]-1; x[1]<position[1]+2; ++x[1]) {
        if (x==position) continue; // know thyself
        MMSP::vector<int> p(x);
        // Take care of periodic boundary conditions
        for (int d=0; d<dim; ++d) MMSP::check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
        if ((p[0] < x0(grid, 0)) || (p[0] >= x1(grid, 0))) continue;
        if ((p[1] < x0(grid, 1)) || (p[1] >= x1(grid, 1))) continue;
        int distance = core_distance + radiussq<2,int>(position, p); // using pre-check distance
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
  else if (dim == 3) {
    MMSP::vector<int> position = getPosition(*core_voxel);
    const int core_distance = (grid(position)).getValue();
    const unsigned int core_id = (grid(position)).getID();
    // Loop over 27 neighboring Voxels
    MMSP::vector<int> x(position);
    for (x[0]=position[0]-1; x[0]<position[0]+2; ++x[0]) {
	    for (x[1]=position[1]-1; x[1]<position[1]+2; ++x[1]) {
  		  for (x[2]=position[2]-1; x[2]<position[2]+2; ++x[2]) {
          if (x==position) continue; // know thyself
          MMSP::vector<int> p(x);
          // Take care of periodic boundary conditions
          for (int d=0; d<dim; ++d) MMSP::check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
          if ((p[0] < x0(grid, 0)) || (p[0] >= x1(grid, 0))) continue;
          if ((p[1] < x0(grid, 1)) || (p[1] >= x1(grid, 1))) continue;
          if ((p[2] < x0(grid, 2)) || (p[2] >= x1(grid, 2))) continue;
          int distance = core_distance + radiussq<3,int>(position, p); // using pre-check distance
          DistanceVoxel* Voxel = &( grid(p) );
          if ( Voxel->getValue() > distance )
          {
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

template<int dim, typename T>
void approximate_voronoi(MMSP::grid<dim, MMSP::sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds) {
  // Implements a fast marching algorithm to generate the distance map
  // Based on code written by Barb Cutler, RPI Comp. Sci. Dept., for CSCI-1200.
  #ifdef DEBUG
  unsigned long tstart=time(NULL);
  #endif
  int id = 0;
	#ifdef MPI_VERSION
  id = MPI::COMM_WORLD.Get_rank();
	#endif
  // Perform the tessellation, using fast-marching fanciness
  if (dim == 2) {
    // Constructor requires global limits of the grid. Magical!
    MMSP::grid<2, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1));
    DistanceVoxel value;
    value.setValue( std::numeric_limits<int>::max() );
    value.setID( 0 );

    // Initialize distance_grid with very-large distances
    for (long i=0; i<nodes(distance_grid); ++i) {
      MMSP::vector<int> pos = grid.position(i);
      for (int d = 0; d < dim; ++d)
        value[d]=pos[d];
      distance_grid(i) = value;
    }

    // create the voxel Heap
    DistanceVoxel_PriorityQueue queue;

    int nseeds = 0;
    for (int i=0; i<id; ++i) nseeds += seeds[i].size();

    // Enqueue this node's seeds
    for (int i=0; i<int(seeds[id].size()); ++i) {
      MMSP::vector<int> pos = seeds[id][i];
      DistanceVoxel* p = &( distance_grid(pos) );
      for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
      p->setValue( 0. );
      p->setID( nseeds + i );
      // Propagate distance from each seed to its neighbors. Start adding to the Heap.
      propagate_distance( p, distance_grid, queue );
    }

    // Fast-march the local seeds
    while( !queue.empty() ) {
      const DistanceVoxel* p = queue.top();
      queue.pop();
      propagate_distance( p, distance_grid, queue );
    }

		#ifdef MPI_VERSION
    // Copy ghost voxels from adjacent ranks
    MMSP::ghostswap(distance_grid);

    // Propagate ghost distances
    int n0[3];
    int n1[3];
    for (int d = 0; d < dim; ++d) {
      n0[d] = N0(grid, d); // Identify processor below in this dimension
      n1[d] = N1(grid, d); // Identify processor above in this dimension
      // Parallel boundary condition
      if (n0[d] != id) { // Ghosts from below
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
        }
      }
      if (n1[d] != id) { // Ghosts from above
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
    while( !queue.empty() ) {
      const DistanceVoxel* p = queue.top();
      queue.pop();
      propagate_distance( p, distance_grid, queue );
    }

    // Copy result from distance_grid to phase-field grid
    for (int i=0; i<nodes(distance_grid); ++i) MMSP::set(grid(i), (distance_grid(i)).getID()) = 1.;
  } else if (dim == 3) {
    // Constructor requires global limits of the grid. Magical!
    MMSP::grid<3, DistanceVoxel> distance_grid(1, g0(grid, 0), g1(grid, 0), g0(grid, 1), g1(grid, 1), g0(grid, 2), g1(grid, 2));
    DistanceVoxel value;
    value.setValue( std::numeric_limits<int>::max() );
    value.setID( 0 );
    for (int i=0; i<nodes(distance_grid); ++i) {
      MMSP::vector<int> pos = position(grid, i);
      for (int d = 0; d < dim; ++d) {
        assert(pos[d] >= x0(distance_grid, d));
        assert(pos[d] < x1(distance_grid, d));
        value[d]=pos[d];
      }
      distance_grid(i) = value;
    }
    // create the voxel Heap
    DistanceVoxel_PriorityQueue queue;

    int nseeds = 0;
    for ( int i=0; i<id; ++i ) nseeds += seeds[i].size();

    // Start queue with this node's seeds
    for ( int i=0; i<int(seeds[id].size()); ++i ) {
      MMSP::vector<int> pos = seeds[id][i];
      DistanceVoxel* p = &( distance_grid(pos) );
      for (int j = 0; j < dim; ++j) assert((pos[j] < x1(grid, j)) && (pos[j] >= x0(grid, j)));
      p->setValue( 0. );
      p->setID( nseeds + i );
      // Propagate distance from each seed to its neighbors. Start adding to the Heap.
      propagate_distance( p, distance_grid, queue );
    }

    // Fast-march the local seeds
    while( !queue.empty() ) {
      const DistanceVoxel* p = queue.top();
      queue.pop();
      propagate_distance( p, distance_grid, queue );
    }

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
    while( !queue.empty() ) {
      const DistanceVoxel* p = queue.top();
      queue.pop();
      propagate_distance( p, distance_grid, queue );
    }

    // Copy result from distance_grid to phase-field grid
    for (int i=0; i<nodes(distance_grid); ++i) MMSP::set(grid(i), (distance_grid(i)).getID()) = 1.;
  }
  else
  {
    std::cerr << "Error: Invalid dimension (" << dim << ") in tessellation." << std::endl;
    std::exit(1);
  }
	#ifdef DEBUG
	#ifdef MPI_VERSION
	if (MPI::COMM_WORLD.Get_rank()==0)
	#endif
	std::cout<<"Completed approximate tessellation ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	#endif
} // approximate_voronoi

void seedswap(std::vector<std::vector<MMSP::vector<int> > >& seeds)
{
	#ifdef MPI_VERSION
  int id=MPI::COMM_WORLD.Get_rank();
  int np=MPI::COMM_WORLD.Get_size();
  // Exchange seeds between all processors
  int send_size=3*seeds[id].size(); // number of integers
  int* send_buffer = new int[send_size]; // number of integers
  send_size = seeds_to_buffer(seeds[id], send_buffer);
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
	delete [] send_buffer; send_buffer=NULL;

	for (int i=0; i<np; ++i) {
		int* p=seed_block+offsets[i];
		seeds_from_buffer(seeds[i], p, seed_sizes[i]);
	}
	delete [] seed_sizes; seed_sizes=NULL;
	delete [] offsets; offsets=NULL;
	delete [] seed_block; seed_block=NULL;
  int vote=1;
  int total_procs=0;
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
  #ifdef DEBUG
  if (id==0) std::cout<<"Synchronized "<<np*seeds[id].size()<<" seeds on "<<total_procs<<" ranks."<<std::endl;
  #endif
	#endif
}

template<int dim>
void seeds_from_poisson_process(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const int& nseeds, std::vector<std::vector<MMSP::vector<int> > >& seeds)
{
  int id=0;
  int np=1;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
	#endif
  std::vector<MMSP::vector<int> > local_seeds; // blank for now
  seeds.clear();
  while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault

  unsigned long int pseudorand_seed = time(NULL);
  if (id == 0) std::cout << "Master seed is " << std::setw(10) << std::right << pseudorand_seed << ". <--- Record this value!" << std::endl;
	#ifdef MPI_VERSION
  pseudorand_seed = pseudorand_seed / (id + 1);
	#endif
  MTRand pseudorand_number( pseudorand_seed );

  // Generate the seeds. Assume repeats don't happen with the Mersenne Twister.
	MMSP::vector<int> seed(3,0.0);
  for (int i=0; i<nseeds; ++i) {
  	for (int d=0; d<dim; ++d)
      seed[d] = g0[d]+pseudorand_number.randInt(g1[d]-g0[d]-1);
    local_seeds.push_back(seed);
  }

	seeds[id].insert(seeds[id].end(), local_seeds.begin(), local_seeds.end());
	seedswap(seeds);
}

template<int dim>
void seeds_to_file(const int g0[dim], const int g1[dim], const std::vector<std::vector<MMSP::vector<int> > >& seeds, const char* filename)
{
  int id=0;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
	#endif
	if (id==0) {
		std::ofstream output(filename);
		output<<dim<<'\n';
		output.precision(6);
		for (std::vector<std::vector<MMSP::vector<int> > >::const_iterator i=seeds.begin(); i!=seeds.end(); ++i)
			for (std::vector<MMSP::vector<int> >::const_iterator j=i->begin(); j!=i->end(); ++j) {
				for (int d=0; d<dim; ++d) output<<std::fixed<<double((*j)[d]-g0[d])/(g1[d]-g0[d])<<'\t';
				output<<'\n';
			}
		output.close();
	}
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Barrier();
	#endif
}

template<int dim>
void seeds_from_file(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const char* seedfilename, std::vector<std::vector<MMSP::vector<int> > >& seeds)
{
  int id=0;
  int np=1;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
	#endif
	std::ifstream input(seedfilename);
	if (!input) {
		std::cerr<<"\nError: "<<seedfilename<<" does not exist.\n"<<std::endl;
		std::exit(-1);
	}
	int seed_dim;
	input >> seed_dim;
	if (dim!=seed_dim) {
		std::cerr<<"\nError: "<<seedfilename<<" is "<<seed_dim<<"-dimensional.\n"<<std::endl;
		std::exit(-1);
	}
  std::vector<MMSP::vector<int> > local_seeds; // blank for now
  seeds.clear();
  while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault
	std::string sentence;
	MMSP::vector<int> seed(3,0.0);
	double s;
	while(input.good()) {
		for (int d=0; d<dim; ++d) {
			input >> s;
			if (s>1.0001) {
				std::cerr<<"\nError: Seeds in "<<seedfilename<<" need to be normalized (found s="<<s<<").\n"<<std::endl;
				std::exit(-1);
			}
			// Seeds in file should be normalized (unit line, square, or cube).
			seed[d]=g0[d]+(g1[d]-g0[d]-1)*s;
		}
		bool isLocal=true;
		for (int d=0; isLocal && d<dim; ++d)
			if (seed[d]<x0[d] || seed[d]>=x1[d]) isLocal=false;
		if (isLocal) seeds[id].push_back(seed);
	}
	seedswap(seeds);
}

template<int dim>
void honeycomb_seeds(const int x0[dim], const int x1[dim], const int g0[dim], const int g1[dim], const int a, std::vector<std::vector<MMSP::vector<int> > >& seeds)
{
  int id=0;
  int np=1;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
	#endif
  std::vector<MMSP::vector<int> > local_seeds; // blank for now
  seeds.clear();
  while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault

  if (id==0) std::cout<<"Making "<<double((g1[0]-g0[0]))/a<<" honeycombs along x."<<std::endl;

	if (dim==2) {
		MMSP::vector<int> seed(3,0.0); // seeds must be 3D for seedswap to work properly
		bool offset=(x0[0]+a)%(a/2)>0?true:false;
		for (seed[1]=(x0[1]+a)%a; seed[1]<x1[1]; seed[1]+=a) {
			for (seed[0]=(x0[0]+a)%a+int(offset)*a/2; seed[0]<x1[0]; seed[0]+=a)
				seeds[id].push_back(seed);
			offset=!offset;
		}
	} else {
		std::cerr<<"\nError: "<<dim<<"-dimensional honeycomb is not implemented yet.\n"<<std::endl;
		std::exit(-1);
	}
	seedswap(seeds);
}

template<int dim, typename T>
void tessellate(MMSP::grid<dim, MMSP::sparse<T> >& grid, const int& nseeds) {
	int* lmin = new int[dim];
	int* lmax = new int[dim];
	int* gmin = new int[dim];
	int* gmax = new int[dim];
	for (int d=0; d<dim; ++d){
		lmin[d]=x0(grid,d);
		lmax[d]=x1(grid,d);
		gmin[d]=g0(grid,d);
		gmax[d]=g1(grid,d);
	}
	std::vector<std::vector<MMSP::vector<int> > > seeds;
	seeds_from_poisson_process<dim>(lmin,lmax,gmin,gmax,nseeds,seeds);
  delete [] lmin;
  delete [] lmax;
  delete [] gmin;
  delete [] gmax;
	// Perform the actual tessellation
	#ifndef MPI_VERSION
	approximate_voronoi<dim,T>(grid, seeds);
	#else
  exact_voronoi<dim,T>(grid, seeds);
	MPI::COMM_WORLD.Barrier();
	#endif
} // tessellate

void print_progress(const int i, const int N) {
  char* timestring;
  static unsigned long tstart;
  struct tm* timeinfo;

  if (i==0) {
    tstart = time(NULL);
    std::time_t rawtime;
    std::time( &rawtime );
    timeinfo = std::localtime( &rawtime );
    timestring = std::asctime(timeinfo);
    timestring[std::strlen(timestring)-1] = '\0';
    std::cout<<timestring<<" ["<<std::flush;
  } else if (i==N) {
    unsigned long deltat = time(NULL)-tstart;
    std::cout<<"•] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
                    <<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
                    <<std::setw(2)<<std::right<<deltat%60<<"s"
                    <<'.'<<std::endl;
  } else if ((20 * i) % N == 0) std::cout<<"• "<<std::flush;
}

#endif // _GENERATORS_
