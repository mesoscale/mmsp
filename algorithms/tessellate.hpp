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
#include <cmath>
#include <ctime>
#include <limits>
#include <cassert>
#include "MersenneTwister.h"

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

#ifndef MPI_VERSION
// MMSP boundary conditions -- copied from MMSP.utility.hpp
enum {
  mirror    = 0,
  Neumann   = 1,
  periodic  = 2,
  parallel  = 3,
  Dirichlet = 4
};

class DistanceVoxel {
public:
  DistanceVoxel() : x( -1 ), y( -1 ), z( 0 ), distance( std::numeric_limits<int>::max() ) {}
  // accessor
  int getX() const { return x; }
  int getY() const { return y; }
  int getZ() const { return z; }
  int getValue() const { return distance; }
  unsigned int getID() const { return id; }
  int& operator [](int i) {
    if (i==0) return this->x;
    else if (i==1) return this->y;
    return this->z;
  }
  const int& operator [](int i) const {
    if (i==0) return this->x;
    else if (i==1) return this->y;
    return this->z;
  }
  // modifier
  void setX( int _x ) { x = _x; }
  void setY( int _y ) { y = _y; }
  void setZ( int _z ) { z = _z; }
  void setValue( int v ) { distance = v; }
  void setID( unsigned int _i ) { id = _i; }
private:
  // REPRESENTATION
  int x;           // a distance voxel
  int y;           // knows its position in the
  int z;           // image and which
  unsigned int id; // seed voxel it's closest to
  int distance; // and how far away that is
};

MMSP::vector<int> getPosition(const DistanceVoxel& dv) {
	MMSP::vector<int> x(3);
	for (int d=0; d<3; ++d) x[d]=dv[d];
	return x;
}

#include "priority_queue.h"
#endif

namespace MMSP {

#ifdef MPI_VERSION
template<int dim, typename T>
void exact_voronoi(MMSP::grid<dim, sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds) {
	int id=MPI::COMM_WORLD.Get_rank();
  // Exact Voronoi tessellation from seeds, based on Euclidean distance function. Runtime is O(Nseeds*L*W*H).
  // seeds must contain every seed from every rank.
  #ifdef DEBUG
  unsigned long tstart=time(NULL);
  #endif
	for (int n=0; n<nodes(grid); ++n) {
		if (id==0) print_progress(n,nodes(grid));
		const MMSP::vector<int> x=position(grid,n);
		int min_distance=std::numeric_limits<int>::max();
		int identity=-1, min_identity=identity;

		for (std::vector<std::vector<MMSP::vector<int> > >::const_iterator rank=seeds.begin(); rank!=seeds.end(); ++rank) {
			for (std::vector<MMSP::vector<int> >::const_iterator s=rank->begin(); s!=rank->end(); ++s) {
				++identity;
				vector<int> seed=*s;
				int distance=radiussq<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
				// Check position across parallel and periodic boundaries to ensure continuity across node boundaries
				for (int d=0; d<dim; d++) check_boundary(seed[d], x0(grid,d), x1(grid,d), b0(grid,d), b1(grid,d));
				if (seed==*s) continue;
				distance=radiussq<dim,int>(x,seed);
				if (distance<min_distance) {
					min_distance=distance;
					min_identity=identity;
				}
			}
		}
		set(grid(n), min_identity) = 1.;
	}
	#ifdef DEBUG
	#ifdef MPI_VERSION
	if (id==0)
	#endif
	std::cout<<"Completed exact tessellation ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	#endif
} // exact_voronoi

#else

template<int dim>
void propagate_distance( const DistanceVoxel* core_voxel, MMSP::grid<dim, DistanceVoxel>& grid, DistanceVoxel_PriorityQueue& queue)
{
  if (dim == 2) {
    MMSP::vector<int> position = getPosition(*core_voxel);
    const int core_distance = (grid(position)).getValue();
    const unsigned int core_id = (grid(position)).getID();
    MMSP::vector<int> p(position);
    // Loop over 9 neighboring Voxels
    for (p[0]=position[0]-1; p[0]<position[0]+2; ++p[0]) {
    	for (p[1]=position[1]-1; p[1]<position[1]+2; ++p[1]) {
        if (p==position) continue; // know thyself
        // Take care of periodic boundary conditions
        for (int d=0; d<dim; ++d) check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
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
    MMSP::vector<int> p(position);
    for (p[0]=position[0]-1; p[0]<position[0]+2; ++p[0]) {
	    for (p[1]=position[1]-1; p[1]<position[1]+2; ++p[1]) {
  		  for (p[2]=position[2]-1; p[2]<position[2]+2; ++p[2]) {
          if (p==position) continue; // know thyself
          // Take care of periodic boundary conditions
          for (int d=0; d<dim; ++d) check_boundary(p[d], x0(grid, d), x1(grid, d), b0(grid, d), b1(grid, d));
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
void approximate_voronoi(MMSP::grid<dim, sparse<T> >& grid, const std::vector<std::vector<MMSP::vector<int> > >& seeds) {
  // Implements a fast marching algorithm to generate the distance map
  // Based on code written by Barb Cutler, RPI Comp. Sci. Dept., for CSCI-1200.
  #ifdef DEBUG
  unsigned long tstart=time(NULL);
  #endif
  int id = 0;
	#ifdef MPI_VERSION
  id = MPI::COMM_WORLD.Get_rank();
  int np = MPI::COMM_WORLD.Get_size();
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
    ghostswap(distance_grid);

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
    for (int i=0; i<nodes(distance_grid); ++i) set(grid(i), (distance_grid(i)).getID()) = 1.;
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
    for (int i=0; i<nodes(distance_grid); ++i) set(grid(i), (distance_grid(i)).getID()) = 1.;
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
#endif

template<int dim, typename T>
void tessellate(const std::vector<MMSP::vector<int> >& local_seeds, std::vector<std::vector<MMSP::vector<int> > >& seeds, MMSP::grid<dim, sparse<T> >& grid) {
  int id=0;
  int np=1;
	#ifdef MPI_VERSION
	id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
	#endif

	#ifndef MPI_VERSION
	seeds[0].insert(seeds[0].end(), local_seeds.begin(), local_seeds.end());
	#else
  // Exchange seeds between all processors
  int send_size=3*local_seeds.size(); // number of integers
  int* send_buffer = new int[send_size]; // number of integers
  send_size = seeds_to_buffer(local_seeds, send_buffer);
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
  if (id==0) std::cout<<"Synchronized "<<total_size/3<<" seeds on "<<total_procs<<" ranks."<<std::endl;
	#endif

  // Perform the actual tessellation
	#ifndef MPI_VERSION
	approximate_voronoi<dim,T>(grid, seeds);
	#else
  exact_voronoi<dim,T>(grid, seeds);
	MPI::COMM_WORLD.Barrier();
  total_procs=0;
  MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
	#endif
} // tessellate


template<int dim, typename T>
void seeds_from_prng(const MMSP::grid<dim, sparse<T> >& grid, const int& nseeds, std::vector<MMSP::vector<int> >& local_seeds) {
  int id=0;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
	#endif

  unsigned long int pseudorand_seed = time(NULL);
  if (id == 0) std::cout << "Master seed is " << std::setw(10) << std::right << pseudorand_seed << ". <--- Record this value!" << std::endl;
	#ifdef MPI_VERSION
  pseudorand_seed = pseudorand_seed / (id + 1);
	#endif
  MTRand pseudorand_number( pseudorand_seed );

  // Generate the seeds. Assume repeats don't happen with the Mersenne Twister.
	MMSP::vector<int> s(3,0.0);
  for (int i=0; i<nseeds; ++i) {
  	for (int d=0; d<dim; ++d)
   		#ifdef MPI_VERSION
      s[d] = g0(grid,d)+pseudorand_number.randInt(g1(grid,d)-g0(grid,d)-1);
   		#else
      s[d] = x0(grid,d)+pseudorand_number.randInt(x1(grid,d)-x0(grid,d)-1);
			#endif
    local_seeds.push_back(s);
  }
}

template<int dim, typename T>
void seeds_from_file(const MMSP::grid<dim, sparse<T> >& grid, const std::string& seed_filename, std::vector<MMSP::vector<int> >& local_seeds) {
  int id=0;
	#ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
	#endif

  if (id == 0) std::cout << "Reading seeds from " << seed_filename << std::endl;
  // Read in the seeds. Assume repeats don't happen.
  std::ifstream sfile(seed_filename.c_str());
  if (!sfile) {
  	if (id==0) std::cerr<<"ERROR: "<<seed_filename<<" not found!"<<std::endl;
  	std::exit(-1);
  }
	MMSP::vector<int> s(3,0);
	std::string line;
	while (getline(sfile,line)) {
		std::stringstream ss(line);
		double x;
  	for (int d=0; d<dim; ++d) {
  		ss>>x;
  		s[d] = g0(grid,d)+int(double(g1(grid,d)-g0(grid,d))*x);
		}
		// Check if this seed is local
		#ifndef MPI_VERSION
   	local_seeds.push_back(s);
   	#else
   	bool isLocal=true;
   	for (int d=0; d<dim && isLocal; d++)
   		if (s[d]<x0(grid,d) || s[d]>=x1(grid,d))
   			isLocal=false;
   	if (isLocal)
   		local_seeds.push_back(s);
   	#endif
  }
  sfile.close();

}


// Randomly generate seeds, using the Mersenne Twister.
template<int dim, typename T>
void tessellate(MMSP::grid<dim, sparse<T> >& grid, const int& nseeds) {
  int np=1;
	#ifdef MPI_VERSION
  np=MPI::COMM_WORLD.Get_size();
	#endif
  std::vector<MMSP::vector<int> > local_seeds; // blank for now
  std::vector<std::vector<MMSP::vector<int> > > seeds;
  while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault

	seeds_from_prng(grid, nseeds, local_seeds);
	tessellate(local_seeds, seeds, grid);

} // random tessellation


// Tessellate seeds specified in a file.
// Seed file should contain tab-delimited columns x	y	z with no header row.
template<int dim, typename T>
void tessellate(MMSP::grid<dim, sparse<T> >& grid, const std::string& seed_filename) {
  int np=1;
	#ifdef MPI_VERSION
  np=MPI::COMM_WORLD.Get_size();
	#endif
  std::vector<MMSP::vector<int> > local_seeds; // blank for now
  std::vector<std::vector<MMSP::vector<int> > > seeds;
  while (int(seeds.size()) <= np) seeds.push_back(local_seeds); // avoid a segfault

		seeds_from_file(grid, seed_filename, local_seeds);
		tessellate(local_seeds, seeds, grid);

} // prescribed tessellation

} // namespace

#endif
