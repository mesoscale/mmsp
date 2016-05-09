// MMSP.grid.hpp
// MMSP grid class definition and implementation
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_GRID
#define MMSP_GRID
#include<iostream>
#include<cstdlib>
#include<cstdarg>
#include<sstream>
#include<cmath>
#include<iomanip>
#include<limits>
#include<cassert>
#include <sys/statvfs.h>
#ifndef RAW
#include<zlib.h>
#endif

#include"MMSP.utility.hpp"
#include"MMSP.scalar.hpp"
#include"MMSP.vector.hpp"
#include"MMSP.sparse.hpp"

namespace MMSP
{

// declaration of grid class
template <int dim, typename T> class grid;

// grid utility functions
template <int dim, typename T> int nodes(const grid<dim, T>& GRID)
{
	return nodes(GRID);
}
template <int dim, typename T> int fields(const grid<dim, T>& GRID)
{
	return fields(GRID);
}
template <int dim, typename T> int ghosts(const grid<dim, T>& GRID)
{
	return ghosts(GRID);
}
template <int dim, typename T> int g0(const grid<dim, T>& GRID, int i)
{
	return g0(GRID, i);
}
template <int dim, typename T> int g1(const grid<dim, T>& GRID, int i)
{
	return g1(GRID, i);
}
template <int dim, typename T> int b0(const grid<dim, T>& GRID, int i)
{
	return b0(GRID, i);
}
template <int dim, typename T> int b1(const grid<dim, T>& GRID, int i)
{
	return b1(GRID, i);
}
template <int dim, typename T> int& b0(grid<dim, T>& GRID, int i)
{
	return b0(GRID, i);
}
template <int dim, typename T> int& b1(grid<dim, T>& GRID, int i)
{
	return b1(GRID, i);
}

// grid utility functions (all directions)
template <int dim, typename T> int x0(const grid<dim, T>& GRID, int i)
{
	return x0(GRID, i);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID, int i)
{
	return x1(GRID, i);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID, int i)
{
	return xmin(GRID, i);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID, int i)
{
	return xmax(GRID, i);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID, int i)
{
	return xlength(GRID, i);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID, int i)
{
	return dx(GRID, i);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID, int i)
{
	return dx(GRID, i);
}

// grid utility functions (x direction)
template <int dim, typename T> int x0(const grid<dim, T>& GRID)
{
	return x0(GRID);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID)
{
	return x1(GRID);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID)
{
	return xmin(GRID);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID)
{
	return xmax(GRID);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID)
{
	return xlength(GRID);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID)
{
	return dx(GRID);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID)
{
	return dx(GRID);
}

// grid utility functions (y direction)
template <int dim, typename T> int y0(const grid<dim, T>& GRID)
{
	return y0(GRID);
}
template <int dim, typename T> int y1(const grid<dim, T>& GRID)
{
	return y1(GRID);
}
template <int dim, typename T> int ymin(const grid<dim, T>& GRID)
{
	return ymin(GRID);
}
template <int dim, typename T> int ymax(const grid<dim, T>& GRID)
{
	return ymax(GRID);
}
template <int dim, typename T> int ylength(const grid<dim, T>& GRID)
{
	return ylength(GRID);
}
template <int dim, typename T> double dy(const grid<dim, T>& GRID)
{
	return dy(GRID);
}
template <int dim, typename T> double& dy(grid<dim, T>& GRID)
{
	return dy(GRID);
}

// grid utility functions (z direction)
template <int dim, typename T> int z0(const grid<dim, T>& GRID)
{
	return z0(GRID);
}
template <int dim, typename T> int z1(const grid<dim, T>& GRID)
{
	return z1(GRID);
}
template <int dim, typename T> int zmin(const grid<dim, T>& GRID)
{
	return zmin(GRID);
}
template <int dim, typename T> int zmax(const grid<dim, T>& GRID)
{
	return zmax(GRID);
}
template <int dim, typename T> int zlength(const grid<dim, T>& GRID)
{
	return zlength(GRID);
}
template <int dim, typename T> double dz(const grid<dim, T>& GRID)
{
	return dz(GRID);
}
template <int dim, typename T> double& dz(grid<dim, T>& GRID)
{
	return dz(GRID);
}


// instantiation of grid class
template <int dim, typename T>
class grid
{
public:
	// constructors
	grid(int FIELDS, int min[dim], int max[dim], int GHOSTS = 1, bool SINGLE = false)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = min[i];
			g1[i] = max[i];
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = GHOSTS;
		#endif

		// setup grid properties
		setup(SINGLE);
	}

	grid(int FIELDS, ...)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		va_list list;
		va_start(list, FIELDS);
		for (int i=0; i<dim; i++) {
			g0[i] = va_arg(list, int);
			g1[i] = va_arg(list, int);
		}

		va_end(list);

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = 1;
		#endif

		// setup grid properties
		setup();
	}

	grid(const grid& GRID)
	{
		// set number of fields
		fields = MMSP::fields(GRID);

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = MMSP::g0(GRID, i);
			g1[i] = MMSP::g1(GRID, i);
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = MMSP::ghosts(GRID);
		#endif

		// setup grid properties
		setup();

		// replace defaults
		for (int i=0; i<dim; i++) {
			b0[i] = MMSP::b0(GRID, i);
			b1[i] = MMSP::b1(GRID, i);
			dx[i] = MMSP::dx(GRID, i);
		}
	}

	template <typename U>
	grid(const grid<dim, U>& GRID, int FIELDS)
	{
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i=0; i<dim; i++) {
			g0[i] = MMSP::g0(GRID, i);
			g1[i] = MMSP::g1(GRID, i);
		}

		// set number of ghosts
		ghosts = 0;

		#ifdef MPI_VERSION
		ghosts = MMSP::ghosts(GRID);
		#endif

		// setup grid properties
		setup();

		// replace defaults
		for (int i=0; i<dim; i++) {
			b0[i] = MMSP::b0(GRID, i);
			b1[i] = MMSP::b1(GRID, i);
			dx[i] = MMSP::dx(GRID, i);
		}
	}

	grid(const char* filename, int GHOSTS = 1)
	{
		// initialize data
		data=NULL;

		// read data from file
		input(filename, GHOSTS);
	}

	void setup(bool SINGLE = false)
	{
		// setup default grid parameters
		for (int i=0; i<dim; i++) {
			x0[i] = g0[i];
			x1[i] = g1[i];

			b0[i] = periodic;
			b1[i] = periodic;

			dx[i] = 1.0;

			p0[i] = 0;
			p1[i] = 1;
			sp[i] = 1;

			n0[i] = 0;
			n1[i] = 0;
		}

		#ifdef MPI_VERSION
		// get global grid data, set neighbor processors
		unsigned int rank = MPI::COMM_WORLD.Get_rank();
		unsigned int np = MPI::COMM_WORLD.Get_size();

		// if bool SINGLE is set to true,
		// we generate grid data only on proc 0
		if (not SINGLE) {
			// compute integral factors of "np"
			int nfac = 0;
			for (unsigned int i=1; i<=np; i++)
				if ((np / i)*i == np) nfac += 1;

			int* factors = new int[nfac];
			nfac = 0;
			for (unsigned int i=1; i<=np; i++)
				if ((np / i)*i == np) {
					factors[nfac] = i;
					nfac += 1;
				}

			// compute global slice areas
			int area[dim];
			for (int i=0; i<dim; i++) {
				area[i] = 1;
				for (int j=0; j<dim; j++)
					if (i != j) area[i] *= (g1[j] - g0[j]);
			}

			// initialize optimal ghost area
			int minarea = -1;

			// compute all combinations of "dim" factors
			int ncom = 1;
			for (int i=0; i<dim; i++)
				ncom *= nfac;

			for (int i=0; i<ncom; i++) {
				int combo[dim];
				int total = i;
				for (int j=0; j<dim; j++) {
					int slice = 1;
					for (int k = j + 1; k < dim; k++)
						slice *= nfac;
					combo[j] = total / slice;
					total -= combo[j] * slice;
				}

				// compute the product of "dim" factors
				unsigned int product = 1;
				for (int j=0; j<dim; j++)
					product *= factors[combo[j]];

				// if product equals "np", compute ghost area
				if (product == np) {
					int test = 0;
					for (int j=0; j<dim; j++)
						test += area[j] * (1 + factors[combo[j]]);
					// choose optimal (smallest) ghost area
					if (test < minarea or minarea < 0) {
						minarea = test;
						for (int k=0; k<dim; k++)
							p1[k] = factors[combo[k]];
					}
				}
			}

			// clean up
			if (factors!=NULL) {
				delete [] factors;
				factors=NULL;
			}

			// compute slice sizes
			for (int i=0; i<dim; i++) {
				sp[i] = 1;
				for (int j = i + 1; j < dim; j++)
					sp[i] *= p1[j];
			}

			// determine local grid limits
			int pos[dim];

			int total = rank;
			for (int i=0; i<dim; i++) {
				// compute position
				pos[i] = total / sp[i];
				total -= pos[i] * sp[i];

				int nom = (g1[i] - g0[i]) / p1[i];
				int rem = (g1[i] - g0[i]) - p1[i] * nom;

				// set limits to "nominal" values
				x0[i] = g0[i] + pos[i] * nom;
				x1[i] = x0[i] + nom;

				// adjust limits according to "remainder"
				if (rem > 0) {
					if (pos[i] < rem) {
						x0[i] += pos[i];
						x1[i] = x0[i] + nom + 1;
					} else if (pos[i] == rem) {
						x0[i] += pos[i];
						x1[i] = x0[i] + nom;
					} else {
						x0[i] += rem;
						x1[i] = x0[i] + nom;
					}
				}
			}

			// determine neighbor processors
			for (int i=0; i<dim; i++) {
				int npos[dim];
				for (int j=0; j<dim; j++)
					npos[j] = pos[j];

				// set neighbor below
				n0[i] = 0;
				npos[i] = (pos[i] - 1 + p1[i]) % p1[i];
				for (int j=0; j<dim; j++)
					n0[i] += sp[j] * npos[j];

				// set neighbor above
				n1[i] = 0;
				npos[i] = (pos[i] + 1) % p1[i];
				for (int j=0; j<dim; j++)
					n1[i] += sp[j] * npos[j];
			}

			// adjust boundary conditions
			for (int i=0; i<dim; i++) {
				if (x0[i] != g0[i]) b0[i] = parallel;
				if (x1[i] != g1[i]) b1[i] = parallel;
			}
		}
		#endif

		// compute slice sizes
		for (int i=0; i<dim; i++) {
			s0[i] = x0[i] - ghosts;
			s1[i] = x1[i] + ghosts;
		}
		for (int i=0; i<dim; i++) {
			sx[i] = 1;
			xx[i] = 1;
			for (int j = i + 1; j < dim; j++) {
				sx[i] *= (s1[j] - s0[j]);
				xx[i] *= (x1[j] - x0[j]);
			}
		}

		// compute number of cells, nodes
		cells = sx[0] * (s1[0] - s0[0]);
		nodes = xx[0] * (x1[0] - x0[0]);

		// resize data structures
		data = new T[cells];
		for (int i=0; i<cells; i++)
			resize(data[i], fields);

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}

	// destructor
	~grid()
	{
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
	}


	// assignment operators
	template <typename U> grid& operator=(const U& value)
	{
		for (int i=0; i<cells; i++)
			data[i] = static_cast<T>(value);
	}

	template <typename U> grid& operator=(const grid<dim, U>& GRID)
	{
		for (int i=0; i<cells; i++)
			data[i] = static_cast<T>(GRID.data[i]);
	}

	template <typename U> grid& operator+=(const grid<dim, U>& GRID)
	{
		for (int i=0; i<cells; i++)
			data[i] += static_cast<T>(GRID.data[i]);
	}

	template <typename U> grid& operator-=(const grid<dim, U>& GRID)
	{
		for (int i=0; i<cells; i++)
			data[i] -= static_cast<T>(GRID.data[i]);
	}


	// subscript operators
	target < dim - 1, 0, T > operator [](int x) const
	{
		check_boundary(x, x0[0], x1[0], b0[0], b1[0]);
		return target < dim - 1, 0, T > (data + (x - s0[0]) * sx[0], s0, sx, x0, x1, b0, b1);
	}

	T& operator()(MMSP::vector<int> x) const
	{
		T* p = data;
		for (int i=0; i<dim; i++) {
			check_boundary(x[i], x0[i], x1[i], b0[i], b1[i]);
			p += (x[i] - s0[i]) * sx[i];
		}
		return *p;
	}

	T& operator()(int i) const
	{
		#ifdef MPI_VERSION
		int x[dim];
		for (int j=0; j<dim; j++) {
			int n = i / xx[j];
			x[j] = n + x0[j];
			i -= n * xx[j];
		}
		i = 0;
		for (int j=0; j<dim; j++)
			i += (x[j] - s0[j]) * sx[j];
		#endif

		return data[i];
	}


	// position utility function
	MMSP::vector<int> position(int i) const
	{
		MMSP::vector<int> x(dim);
		for (int j=0; j<dim; j++) {
			int n = i / xx[j];
			x[j] = n + x0[j];
			i -= n * xx[j];
		}
		return x;
	}


	// parallelization
	void ghostswap()
	{
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		for (int i=0; i<dim; i++) {
			if (1) {
				// send to processor above and receive from processor below
				int send_proc = n1[i];
				int recv_proc = n0[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j] - ghosts;
					send_max[j] = x1[j] + ghosts;
					recv_min[j] = x0[j] - ghosts;
					recv_max[j] = x1[j] + ghosts;
				}

				send_min[i] = x1[i] - ghosts;
				send_max[i] = x1[i];
				recv_min[i] = x0[i] - ghosts;
				recv_max[i] = x0[i];

				unsigned long send_size = this->buffer_size(send_min, send_max);
				unsigned long recv_size = 0;

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 100); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 100); // receive number of ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				send_size = this->to_buffer(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 200); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 200); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				this->from_buffer(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer
				delete [] send_buffer;
				send_buffer=NULL;
				delete [] recv_buffer;
				recv_buffer=NULL;
			}

			if (1) {
				// send to processor below and receive from processor above
				int send_proc = n0[i];
				int recv_proc = n1[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j] - ghosts;
					send_max[j] = x1[j] + ghosts;
					recv_min[j] = x0[j] - ghosts;
					recv_max[j] = x1[j] + ghosts;
				}

				send_min[i] = x0[i];
				send_max[i] = x0[i] + ghosts;
				recv_min[i] = x1[i];
				recv_max[i] = x1[i] + ghosts;

				unsigned long send_size = this->buffer_size(send_min, send_max);
				unsigned long recv_size = 0;

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 300); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 300); // receive number of incoming ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				send_size = this->to_buffer(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 400); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 400); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();

				this->from_buffer(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer

				delete [] send_buffer;
				send_buffer=NULL;
				delete [] recv_buffer;
				recv_buffer=NULL;
			}
		}
		MPI::COMM_WORLD.Barrier();
		#endif
	}

	void ghostswap(const int sublattice)   // ghostswap for Monte Carlo communiation reduction.
	{
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		for (int i=0; i<dim; i++) {
			if (1) {
				// send to processor above and receive from processor below
				int send_proc = n1[i];
				int recv_proc = n0[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j] - ghosts;
					send_max[j] = x1[j] + ghosts;
					recv_min[j] = x0[j] - ghosts;
					recv_max[j] = x1[j] + ghosts;
				}

				send_min[i] = x1[i] - ghosts; //ghosts = 1; x0 x1 excluding ghost
				send_max[i] = x1[i];
				recv_min[i] = x0[i] - ghosts;
				recv_max[i] = x0[i];

				if (dim==2) {
					if (sublattice==0) {
						send_min[i] = ((x1[i] - ghosts)%2==0? x1[i]-ghosts:x1[i]); // only send even number
						recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i]-ghosts:x0[i]); // only recv even number
						send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
					}

					if (sublattice==3) {
						send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i]-ghosts:x1[i]); // only send odd number
						recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i]-ghosts:x0[i]); // only recv odd number
						send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
					}

					if (sublattice==1) {
						if (i==0) {
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i]-ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i]-ghosts:x0[i]); // only recv even number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) {
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i]-ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i]-ghosts:x0[i]); // only recv odd number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						}
					}

					if (sublattice==2) {
						if (i==0) {
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i]-ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i]-ghosts:x0[i]); // only recv odd number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) {
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i]-ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i]-ghosts:x0[i]); // only recv even number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						}
					}
				} else if (dim==3) {
					if (sublattice==0) { // 0,0,0
						send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
						recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
						for (int l=0; l<dim; l++) {
							if (l!=i) {
								send_min[l] = (send_min[l]%2==0?send_min[l]:send_min[l]+1); //the first element of that 1 layer strip should be even
								recv_min[l] = (recv_min[l]%2==0?recv_min[l]:recv_min[l]+1); //the first element of that 1 layer strip should be even
							}
						}
					}

					else if (sublattice==1) { //0,0,1
						if (i==0) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==2) { //0,1,0
						if (i==0) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==3) { //0,1,1
						if (i==0) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==4) { //1,0,0
						if (i==0) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==5) { //1,0,1
						if (i==0) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==6) { //1,1,0
						if (i==0) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv even
							send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i] - ghosts:x1[i]); // only send odd number
							recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i] - ghosts:x0[i]); // only recv odd number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv odd
							send_min[i] = ((x1[i] - ghosts)%2==0? x1[i] - ghosts:x1[i]); // only send even number
							recv_min[i] = ((x0[i] - ghosts)%2==0? x0[i] - ghosts:x0[i]); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==7) {
						send_min[i] = ((x1[i] - ghosts)%2!=0? x1[i]-ghosts:x1[i]); // only send odd number
						recv_min[i] = ((x0[i] - ghosts)%2!=0? x0[i]-ghosts:x0[i]); // only recv odd number
						for (int l=0; l<dim; l++) {
							if (l!=i) {
								send_min[l] = (send_min[l]%2!=0?send_min[l]:send_min[l]+1); //the first element of that 1 layer strip should be odd
								recv_min[l] = (recv_min[l]%2!=0?recv_min[l]:recv_min[l]+1); //the first element of that 1 layer strip should be odd
							}
						}
					}
				}//dim == 3

				unsigned long send_size = this->buffer_size_save(send_min, send_max);
				unsigned long recv_size = 0;

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 100); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 100); // receive number of ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				send_size = this->to_buffer_save(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 200); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 200); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();

				this->from_buffer_save(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer

				delete [] send_buffer;
				send_buffer=NULL;
				delete [] recv_buffer;
				recv_buffer=NULL;
			}

			if (1) {
				// send to processor below and receive from processor above
				int send_proc = n0[i];
				int recv_proc = n1[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j=0; j<dim; j++) {
					send_min[j] = x0[j] - ghosts;
					send_max[j] = x1[j] + ghosts;
					recv_min[j] = x0[j] - ghosts;
					recv_max[j] = x1[j] + ghosts;
				}

				send_min[i] = x0[i];
				send_max[i] = x0[i] + ghosts;
				recv_min[i] = x1[i];
				recv_max[i] = x1[i] + ghosts;

				if (dim==2) {
					if (sublattice==0) {
						send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
						recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
						send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
					}

					if (sublattice==3) {
						send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
						recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
						send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
					}

					if (sublattice==1) {
						if (i==0) {
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) {
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						}
					}

					if (sublattice==2) {
						if (i==0) {
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2==0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2==0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) {
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[abs(1-i)] = (send_min[abs(1-i)]%2!=0?send_min[abs(1-i)]:send_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
							recv_min[abs(1-i)] = (recv_min[abs(1-i)]%2!=0?recv_min[abs(1-i)]:recv_min[abs(1-i)]+1); //the first element of that 1 layer strip should be odd
						}
					}
				} else if (dim==3) {
					if (sublattice==0) { // 0,0,0
						send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
						recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
						for (int l=0; l<dim; l++) {
							if (l!=i) {
								send_min[l] = (send_min[l]%2==0?send_min[l]:send_min[l]+1); //the first element of that 1 layer strip should be even
								recv_min[l] = (recv_min[l]%2==0?recv_min[l]:recv_min[l]+1); //the first element of that 1 layer strip should be even
							}
						}
					}

					else if (sublattice==1) { //0,0,1
						if (i==0) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==2) { //0,1,0
						if (i==0) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==3) { //0,1,1
						if (i==0) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv even
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2==0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be even
							recv_min[0] = (recv_min[0]%2==0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be even
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==4) { //1,0,0
						if (i==0) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==5) { //1,0,1
						if (i==0) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==1) { //only send/recv even
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2!=0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be odd
							recv_min[2] = (recv_min[2]%2!=0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be odd
						} else if (i==2) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2==0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be even
							recv_min[1] = (recv_min[1]%2==0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be even
						}
					}

					else if (sublattice==6) { //1,1,0
						if (i==0) { //only send/recv odd
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==1) { //only send/recv even
							send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
							recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[2] = (send_min[2]%2==0?send_min[2]:send_min[2]+1); //the first element of that 1 layer strip should be even
							recv_min[2] = (recv_min[2]%2==0?recv_min[2]:recv_min[2]+1); //the first element of that 1 layer strip should be even
						} else if (i==2) { //only send/recv odd
							send_min[i] = (x0[i]%2==0? x0[i]:x0[i]+ghosts); // only send even number
							recv_min[i] = (x1[i]%2==0? x1[i]:x1[i]+ghosts); // only recv even number
							send_min[0] = (send_min[0]%2!=0?send_min[0]:send_min[0]+1); //the first element of that 1 layer strip should be odd
							recv_min[0] = (recv_min[0]%2!=0?recv_min[0]:recv_min[0]+1); //the first element of that 1 layer strip should be odd
							send_min[1] = (send_min[1]%2!=0?send_min[1]:send_min[1]+1); //the first element of that 1 layer strip should be odd
							recv_min[1] = (recv_min[1]%2!=0?recv_min[1]:recv_min[1]+1); //the first element of that 1 layer strip should be odd
						}
					}

					else if (sublattice==7) { //1,1,1
						send_min[i] = (x0[i]%2!=0? x0[i]:x0[i]+ghosts); // only send odd number
						recv_min[i] = (x1[i]%2!=0? x1[i]:x1[i]+ghosts); // only recv odd number
						for (int l=0; l<dim; l++) {
							if (l!=i) {
								send_min[l] = (send_min[l]%2!=0?send_min[l]:send_min[l]+1); //the first element of that 1 layer strip should be odd
								recv_min[l] = (recv_min[l]%2!=0?recv_min[l]:recv_min[l]+1); //the first element of that 1 layer strip should be odd
							}
						}
					}
				}//dim == 3

				unsigned long send_size = this->buffer_size_save(send_min, send_max);
				unsigned long recv_size = 0;

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 300); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 300); // receive number of incoming ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];

				send_size = this->to_buffer_save(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 400); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 400); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();

				this->from_buffer_save(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer

				delete [] send_buffer;
				send_buffer=NULL;
				delete [] recv_buffer;
				recv_buffer=NULL;
			}
		}
		MPI::COMM_WORLD.Barrier();
		#endif
	}

	unsigned long buffer_size_save(const int min[dim], const int max[dim]) const
	{
		return buffer_size_save(data, 0, min, max);
	}
	unsigned long buffer_size_save(T* p, int i, const int min[dim], const int max[dim]) const
	{
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x+=2)
				size += MMSP::buffer_size(*(p + (x - s0[i]) * sx[i]));

		else
			for (int x = min[i]; x < max[i]; x+=2)
				size += buffer_size_save(p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long to_buffer_save(char* buffer, const int min[dim], const int max[dim]) const
	{
		return to_buffer_save(buffer, data, 0, min, max);
	}

	unsigned long to_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
	{
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x+=2)
				size += MMSP::to_buffer(*(p + (x - s0[i]) * sx[i]), buffer + size);
		else
			for (int x = min[i]; x < max[i]; x+=2)
				size += to_buffer_save(buffer + size, p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long from_buffer_save(char* buffer, const int min[dim], const int max[dim])
	{
		return from_buffer_save(buffer, data, 0, min, max);
	}

	unsigned long from_buffer_save(char* buffer, T* p, int i, const int min[dim], const int max[dim])
	{
		unsigned long size = 0;
		if (i == dim - 1) {
			for (int x = min[i]; x < max[i]; x+=2)
				size += MMSP::from_buffer(*(p + (x - s0[i]) * sx[i]), buffer + size);
		} else {
			for (int x = min[i]; x < max[i]; x+=2)
				size += from_buffer_save(buffer + size, p + (x - s0[i]) * sx[i], i + 1, min, max);
		}
		return size;
	}
	// buffer I/O
	unsigned long buffer_size() const
	{
		return buffer_size(x0, x1);
	}

	unsigned long buffer_size(const int min[dim], const int max[dim]) const
	{
		return buffer_size(data, 0, min, max);
	}

	unsigned long buffer_size(T* p, int i, const int min[dim], const int max[dim]) const
	{
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x++)
				size += MMSP::buffer_size(*(p + (x - s0[i]) * sx[i]));
		else
			for (int x = min[i]; x < max[i]; x++)
				size += buffer_size(p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long to_buffer(char* buffer) const
	{
		return to_buffer(buffer, x0, x1);
	}

	unsigned long to_buffer(char* buffer, const int min[dim], const int max[dim]) const
	{
		return to_buffer(buffer, data, 0, min, max);
	}

	unsigned long to_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const
	{
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x++)
				size += MMSP::to_buffer(*(p + (x - s0[i]) * sx[i]), buffer + size);
		else
			for (int x = min[i]; x < max[i]; x++)
				size += to_buffer(buffer + size, p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long from_buffer(char* buffer)
	{
		return from_buffer(buffer, x0, x1);
	}

	unsigned long from_buffer(char* buffer, const int min[dim], const int max[dim])
	{
		return from_buffer(buffer, data, 0, min, max);
	}

	unsigned long from_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim])
	{
		unsigned long size = 0;
		if (i == dim - 1) {
			for (int x = min[i]; x < max[i]; x++)
				size += MMSP::from_buffer(*(p + (x - s0[i]) * sx[i]), buffer + size);
		} else {
			for (int x = min[i]; x < max[i]; x++)
				size += from_buffer(buffer + size, p + (x - s0[i]) * sx[i], i + 1, min, max);
		}
		return size;
	}

	// file I/O
	void input(const char* filename, int GHOSTS = 1, int SINGLE = false)
	{
		// file open error check
		std::ifstream input(filename);
		if (!input) {
			std::cerr << "File input error: could not open ";
			std::cerr << filename << "." << std::endl;
			exit(-1);
		}

		// grid data type error check
		std::string type;
		getline(input, type, '\n');
		if (type != name(*this)) {
			std::cerr << "File read error: wrong data type (" << type << ")." << std::endl;
			exit(-2);
		}

		// dimension error check
		int dimen;
		input >> dimen;
		if (dimen != dim) {
			std::cerr << "File read error: wrong dimension (" << dimen << ")." << std::endl;
			exit(-3);
		}

		// read number of fields
		input >> fields;

		// read grid size
		for (int i=0; i<dim; i++)
			input >> g0[i] >> g1[i];

		// set number of ghosts
		ghosts = GHOSTS;
		#ifndef MPI_VERSION
		ghosts = 0;
		#endif

		// setup grid parameters
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		setup(SINGLE);

		// read cell spacing
		for (int i=0; i<dim; i++)
			input >> dx[i];

		// ignore trailing endlines
		input.ignore(10, '\n');

		// input grid data
		read(input, GHOSTS);

		// ghostswap if necessary
		if (not SINGLE) ghostswap();
		input.close();
	}

	void read(std::ifstream& file, int GHOSTS = 1)
	{
		/*
			Reads each block of the input file and constructs a
			temporary grid, "GRID". Then compares the spatial extents of
			the block (lmin, lmax) to the local grid (this->x0, this->x1).
			Data is copied from the region where both overlap, if it exists.
			Then the next block is read.
		*/
		// read number of blocks
		int blocks;
		file.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));

		#ifdef GRIDDEBUG
		int actual_read=0;
		unsigned long data_read=0;
		#endif

		// for each block...
		for (int i=0; i<blocks; i++) {
			int lmin[dim];
			int lmax[dim];
			// read block limits
			for (int j=0; j<dim; j++) {
				file.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
				file.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
			}
			int blo[dim];
			int bhi[dim];
			// read boundary conditions
			for (int j=0; j<dim; j++) {
				file.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
				file.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
			}

			// read block data
			unsigned long size_on_disk, size_in_mem;
			file.read(reinterpret_cast<char*>(&size_in_mem), sizeof(size_in_mem)); // read grid data size
			file.read(reinterpret_cast<char*>(&size_on_disk), sizeof(size_on_disk)); // read compressed size

			// find overlapping region
			int min[dim];
			int max[dim];
			bool overlap = true;
			for (int j=0; j<dim; j++) {
				min[j] = (lmin[j] < x0[j]) ? x0[j] : lmin[j];
				max[j] = (lmax[j] > x1[j]) ? x1[j] : lmax[j];
				if (min[j] >= max[j]) overlap = false;
			}

			if (overlap) {
				#ifdef GRIDDEBUG
				++actual_read;
				data_read+=size_on_disk;
				#endif
				char* buffer = new char[size_on_disk];
				file.read(buffer, size_on_disk);
				grid<dim, T> GRID(fields, lmin, lmax, 0, true);
				if (size_in_mem!=size_on_disk) {
					#ifdef RAW
					std::cerr<<"Unable to uncompress data: compiled without zlib."<<std::endl;
					exit(1);
					#else
					// Uncompress data
					char* raw = new char[size_in_mem];
					int status = uncompress(reinterpret_cast<Bytef*>(raw), &size_in_mem, reinterpret_cast<Bytef*>(buffer), size_on_disk);
					switch ( status ) {
					case Z_OK:
						break;
					case Z_MEM_ERROR:
						std::cerr << "Uncompress: out of memory." << std::endl;
						exit(1);    // quit.
						break;
					case Z_BUF_ERROR:
						std::cerr << "Uncompress: output buffer wasn't large enough." << std::endl;
						exit(1);    // quit.
						break;
					}
					GRID.from_buffer(raw);
					delete [] raw;
					raw=NULL;
					#endif
				} else GRID.from_buffer(buffer);
				delete [] buffer;
				buffer=NULL;

				// copy block data that overlaps
				unsigned long size = GRID.buffer_size(min, max);
				buffer = new char[size];
				GRID.to_buffer(buffer, min, max);
				this->from_buffer(buffer, min, max);
				delete [] buffer;
				buffer=NULL;

				// set boundary conditions from file
				for (int j=0; j<dim; j++) {
					if (x0[j]==lmin[j]) b0[j]=blo[j];
					if (x1[j]==lmax[j]) b1[j]=bhi[j];
				}
			} else {
				// No overlap.
				file.seekg(size_on_disk, file.cur);
			}
		}

		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}

	void output(const char* filename) const
	{
		#ifndef MPI_VERSION
		unsigned int np=1;
		// file open error check
		std::ofstream output(filename);
		if (!output) {
			std::cerr << "File output error: could not open ";
			std::cerr << filename << "." << std::endl;
			exit(-1);
		}

		std::stringstream outstr;

		// get grid data type
		std::string type = name(*this);
		outstr << type << '\n';

		// get grid dimension
		outstr << dim << '\n';

		// get number of fields
		outstr << fields << '\n';

		// get grid size
		for (int i=0; i<dim; i++) outstr << g0[i] << " " << g1[i] << '\n';

		// get cell spacing
		for (int i=0; i<dim; i++) outstr << dx[i] << '\n';

		// Write file header to file
		output.write(outstr.str().c_str(), outstr.str().size());
		// Write number of blocks (processors) to file
		output.write(reinterpret_cast<const char*>(&np), sizeof(np));

		// get grid data to write
		char* buffer;
		unsigned long size_of_buffer = this->write_buffer(buffer);
		// output grid data
		output.write(buffer, size_of_buffer);
		delete [] buffer;
		buffer=NULL;


		#else
		/* MPI-IO write to disk */

		// MPI C-style functions require char*, not const char*, filename
		char fname[FILENAME_MAX] = {}; // initializes array with null chars
		// C-style strings are null-terminated ('\0') by definition
		for (unsigned int i=0; i<FILENAME_MAX && filename[i]!='\0'; i++)
			fname[i]=filename[i];
		MPI::COMM_WORLD.Barrier();
		unsigned int rank = MPI::COMM_WORLD.Get_rank();
		unsigned int np = MPI::COMM_WORLD.Get_size();
		MPI_Request request;
		MPI_Status status;
		int mpi_err = 0;

		// Read filesystem block size (using statvfs). Default to 4096 B.
		struct statvfs buf;
		const unsigned long blocksize = (statvfs(".", &buf) == -1)?4096:buf.f_bsize;

		if (blocksize<1048576) {
			// Standard MPI-IO: every rank writes to disk
			#ifdef BGQ
			if (rank==0) std::cout<<"Bug: using stock MPI-IO, instead of BGQ IO!"<<std::endl;
			#endif
			MPI_File output;
			mpi_err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY|MPI_MODE_EXCL|MPI_MODE_CREATE, MPI_INFO_NULL, &output);
			if (mpi_err != MPI_SUCCESS) {
				if (rank==0)
					MPI_File_delete(fname,MPI_INFO_NULL);
				mpi_err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY|MPI_MODE_EXCL|MPI_MODE_CREATE, MPI_INFO_NULL, &output);
				assert(mpi_err==MPI_SUCCESS);
			}
			if (!output) {
				std::cerr << "File output error: could not open " << fname << "." << std::endl;
				exit(-1);
			}

			// Generate MMSP header from rank 0
			unsigned long header_offset=0;
			if (rank == 0) {
				std::stringstream outstr;
				// get grid data type
				std::string type = name(*this);

				outstr << type << '\n';
				outstr << dim << '\n';
				outstr << fields << '\n';

				for (int i=0; i<dim; i++) outstr << g0[i] << " " << g1[i] << '\n'; // global grid dimensions
				for (int i=0; i<dim; i++) outstr << dx[i] << '\n'; // grid spacing

				// Write MMSP header to file
				header_offset=outstr.str().size();
				char* header = new char[header_offset];
				for (unsigned int i=0; i<header_offset; i++)
					header[i] = outstr.str()[i];
				MPI_File_sync(output);
				MPI_File_iwrite_at(output,0,header, header_offset, MPI_CHAR, &request);
				MPI_Wait(&request, &status);
				MPI_File_sync(output);
				// Write number of blocks (processors) to file
				MPI_File_iwrite_at(output,header_offset,reinterpret_cast<char*>(&np), sizeof(np), MPI_CHAR, &request);
				MPI_Wait(&request, &status);
				MPI_File_sync(output);
				header_offset+=sizeof(np);
				delete [] header;
				header=NULL;
			}
			MPI::COMM_WORLD.Barrier();
			MPI_File_sync(output);
			MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0

			// get grid data to write
			char* buffer=NULL;
			unsigned long size_of_buffer = this->write_buffer(buffer);
			assert(buffer!=NULL);

			// Compute file offsets based on buffer sizes
			unsigned long* datasizes = new unsigned long[np];
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Allgather(&size_of_buffer, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);

			// Pre-allocate disk space
			unsigned long filesize=0;
			for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
			MPI::COMM_WORLD.Barrier();
			MPI_File_preallocate(output, filesize);

			unsigned long* offsets = new unsigned long[np];
			offsets[0]=header_offset;
			for (unsigned int n=1; n<np; ++n) {
				assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
				offsets[n]=offsets[n-1]+datasizes[n-1];
			}
			#ifdef GRIDDEBUG
			assert(datasizes[rank]==size_of_buffer);
			if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
			#endif

			// Write buffer to disk
			MPI_File_sync(output);
			MPI::COMM_WORLD.Barrier();
			MPI_File_iwrite_at(output,offsets[rank],buffer,datasizes[rank],MPI_CHAR,&request);
			MPI_Wait(&request, &status);
			#ifdef GRIDDEBUG
			int error, write_errors=0;
			MPI_Get_count(&status, MPI_INT, &error);
			if (error!=0)
				std::cerr<<"  Error on Rank "<<rank<<": "<<MPI::Get_error_class(error)<<std::endl;
			MPI::COMM_WORLD.Allreduce(&error, &write_errors, 1, MPI_INT, MPI_SUM);
			if (rank==0) std::cout<<"  Write finished on "<<np-write_errors<<'/'<<np<<" ranks."<<std::endl;
			assert(write_errors==0);
			#endif
			delete [] buffer;
			buffer=NULL;

			MPI::COMM_WORLD.Barrier();
			MPI_File_sync(output);
			// Make sure everything's written before closing the file.
			MPI_Offset actual_size;
			MPI_File_get_size(output,&actual_size);
			if (rank==0) {
				#ifdef GRIDDEBUG
				std::cout<<fname<<" should be "<<offsets[np-1]+datasizes[np-1]<<" B;";
				std::cout<<" wrote "<<actual_size<<" B to disk."<<std::endl;
				#endif
				assert(offsets[np-1]+datasizes[np-1] == static_cast<unsigned long>(actual_size));
			}
			MPI::COMM_WORLD.Barrier();
			MPI_File_close(&output);
			delete [] offsets;
			offsets=NULL;
			delete [] datasizes;
			datasizes=NULL;

		} else {

			// Custom MPI-IO for large-block filesystems, such as IBM GPFS:
			// Two-stage aggregate+write IO, one rank per block at MOST
			// create buffer pointers
			unsigned long* datasizes = NULL;
			unsigned long* offsets = NULL;
			unsigned long* aoffsets = NULL;
			unsigned long* misalignments = NULL;
			char* databuffer=NULL;
			char* headbuffer=NULL;
			char* filebuffer=NULL;
			unsigned int* writeranks=NULL;
			MPI_Request* recvrequests = NULL;
			MPI_Status* recvstatuses = NULL;
			int mpi_err = 0;

			// get grid data to write
			const unsigned long size_of_buffer = write_buffer(databuffer);
			assert(databuffer!=NULL);
			// Generate MMSP header from rank 0
			unsigned long header_offset=0;
			if (rank==0) {
				// get grid data type
				std::string type = name(*this);

				std::stringstream outstr;
				outstr << type << '\n';
				outstr << dim << '\n';
				outstr << MMSP::fields(*this) << '\n';

				for (int i=0; i<dim; i++) outstr << MMSP::g0(*this,i) << " " << MMSP::g1(*this,i) << '\n'; // global grid dimensions
				for (int i=0; i<dim; i++) outstr << MMSP::dx(*this,i) << '\n'; // grid spacing

				// Write MMSP header to buffer
				header_offset=outstr.str().size();
				headbuffer = new char[header_offset+sizeof(rank)];
				memcpy(headbuffer, outstr.str().c_str(), header_offset);
				memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
				header_offset+=sizeof(rank);
			}
			MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0
			#ifdef GRIDDEBUG
			if (rank==0) std::cout<<"Prepared file header."<<std::endl;
			#endif
			MPI::COMM_WORLD.Barrier();

			// Compute file offsets based on buffer sizes
			datasizes = new unsigned long[np];
			MPI::COMM_WORLD.Allgather(&size_of_buffer, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);
			#ifdef GRIDDEBUG
			if (rank==0) std::cout<<"Synchronized data sizes."<<std::endl;
			#endif

			// Determine disk space requirement
			unsigned long filesize=header_offset;
			for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
			MPI::COMM_WORLD.Barrier();

			offsets = new unsigned long[np];
			offsets[0]=header_offset;
			for (unsigned int n=1; n<np; ++n) {
				assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
				offsets[n]=offsets[n-1]+datasizes[n-1];
			}
			offsets[0]=0;
			#ifdef GRIDDEBUG
			assert(datasizes[rank]==size_of_buffer);
			if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
			#endif

			// Calculate number of  writers & write size
			unsigned long blocks = filesize/blocksize;
			while (blocks*blocksize<filesize)	++blocks;
			const unsigned int nwriters = (blocks>np)?np:blocks;
			const unsigned long writesize=blocksize*(blocks/nwriters);
			assert(writesize % blocksize==0);
			const unsigned long excessblocks=blocks % nwriters;
			bool isWriter=false;
			#ifdef GRIDDEBUG
			if (rank==0) std::cout<<"  Preparing "<<nwriters<<" aggregator/writers; writesize is "<<writesize<<" B, with "<<excessblocks<<" excess blocks."<<std::endl;
			#endif

			// Scan to determine which ranks are writers
			writeranks = new unsigned int[nwriters+1];
			aoffsets = new unsigned long[nwriters];
			writeranks[nwriters]=np-1; // generalization for last writer's benefit
			unsigned int temprank=0;
			for (unsigned int w=0; w<nwriters; w++) {
				unsigned long ws=(w<=excessblocks)?writesize+blocksize:writesize;
				// file offset of the w^th writer
				aoffsets[w]=(w>0)?ws+aoffsets[w-1]:0;
				while ((aoffsets[w] > offsets[temprank]+datasizes[temprank]) && temprank<np)
					temprank++;
				writeranks[w]=temprank;
				if (rank==temprank)
					isWriter=true;
				temprank++;
			}

			// Determine which rank to send data to
			unsigned int prevwriter=nwriters, nextwriter=0;
			if (rank==0) {
				prevwriter=0;
			} else {
				while (writeranks[prevwriter]>=rank)
					--prevwriter;
			}
			if (rank>=writeranks[nwriters-1]) {
				nextwriter=nwriters;
			} else {
				while (writeranks[nextwriter]<=rank)
					++nextwriter;
			}

			unsigned long ws = writesize;
			if (nextwriter<=excessblocks)
				ws+=blocksize;
			if (rank>=writeranks[nwriters-1])
				ws=filesize-aoffsets[nwriters-1]; // last block may be only partially-filled

			unsigned long deficiency=0;
			if (rank>0) {
				unsigned long prevws = (prevwriter>=excessblocks)?writesize:writesize+blocksize;
				deficiency = aoffsets[prevwriter]+prevws - offsets[rank];
				if (deficiency>datasizes[rank])
					deficiency=datasizes[rank];
			}
			// Collect block misalignments
			misalignments = new unsigned long[np];
			MPI::COMM_WORLD.Barrier();
			MPI::COMM_WORLD.Allgather(&deficiency, 1, MPI_UNSIGNED_LONG, misalignments, 1, MPI_UNSIGNED_LONG);

			#ifdef GRIDDEBUG
			if (datasizes[rank]-deficiency>ws)
				std::fprintf(stderr, "Error on Rank %u, alignment: buffered %lu B > writesize %lu B.\n", rank, datasizes[rank]-deficiency, ws);
			#endif

			// Accumulate data
			const unsigned int silentranks=writeranks[nextwriter]-rank; // number of MPI ranks between this rank and the next writer
			MPI_Request sendrequest;
			MPI::COMM_WORLD.Barrier();
			if (isWriter) {
				// This rank is a writer.
				assert(misalignments[rank] < datasizes[rank]);
				#ifdef GRIDDEBUG
				if (rank>0 && writeranks[prevwriter+1]!=rank)
					std::fprintf(stderr, "Error on Rank %u, writer ID: %u != %u\n", rank, writeranks[prevwriter+1], rank);
				#endif

				// Copy local data into filebuffer
				filebuffer = new char[ws];
				char* p = filebuffer;
				if (rank==0) {
					memcpy(p, headbuffer, header_offset);
					p+=header_offset;
				}
				#ifdef GRIDDEBUG
				if (datasizes[rank]-misalignments[rank]>ws)
					std::fprintf(stderr, "Error on Rank %u, memcpy: %lu B > %lu B\n", rank, datasizes[rank]-misalignments[rank], ws);
				#endif
				char* q=databuffer+misalignments[rank];
				memcpy(p, q, datasizes[rank]-misalignments[rank]);
				p+=datasizes[rank]-misalignments[rank];

				// Recv remote data into filebuffer
				if (silentranks>0) {
					recvrequests = new MPI_Request[silentranks];
					recvstatuses = new MPI_Status[silentranks];
				}
				for (unsigned int i=0; i<silentranks && rank+i+1<np; i++) {
					unsigned int recv_proc = rank+i+1;
					assert(recv_proc!=rank && recv_proc<np);
					#ifdef GRIDDEBUG
					if (recv_proc<rank || recv_proc>np)
						std::fprintf(stderr, "Error on Rank %u, receiving: recv_proc=%i\n", rank, recv_proc);
					#endif
					unsigned long recv_size = misalignments[recv_proc];
					if (recv_size==0) continue;
					#ifdef GRIDDEBUG
					if (p+recv_size>filebuffer+ws)
						std::fprintf(stderr, "Error on Rank %u, receiving from %i: %lu B > %lu B\n", rank, recv_proc, p-filebuffer, ws-recv_size);
					#endif
					MPI_Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc, MPI::COMM_WORLD, &recvrequests[i]);
					p+=recv_size;
				}
				#ifdef GRIDDEBUG
				if (p-filebuffer!=int(ws))
					std::fprintf(stderr, "Error on Rank %u, total received: %i B != %lu B\n", rank, int(p-filebuffer), ws);
				#endif
				if (rank>0 && misalignments[rank]>0) {
					q=databuffer;
					assert(writeranks[prevwriter]<rank);
					MPI_Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
				}
			}
			if (misalignments[rank] >= datasizes[rank]) {
				assert(writeranks[prevwriter]<rank && writeranks[prevwriter]<np);
				MPI_Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
			}
			if (recvrequests != NULL)
				MPI_Waitall(silentranks, recvrequests, recvstatuses);
			if (rank>0) MPI_Wait(&sendrequest, &status);
			MPI::COMM_WORLD.Barrier();

			// file open error check
			#ifdef GRIDDEBUG
			if (rank==0) std::cout<<"  Opening "<<std::string(fname)<<" for output."<<std::endl;
			#endif
			MPI_File output;
			mpi_err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY|MPI_MODE_EXCL|MPI_MODE_CREATE, MPI_INFO_NULL, &output);
			if (mpi_err != MPI_SUCCESS) {
				if (rank==0)
					MPI_File_delete(fname,MPI_INFO_NULL);
				mpi_err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY|MPI_MODE_EXCL|MPI_MODE_CREATE, MPI_INFO_NULL, &output);
				assert(mpi_err==MPI_SUCCESS);
			}
			if (!output) {
				if (rank==0) std::cerr << "File output error: could not open " << fname << "." << std::endl;
				if (rank==0) std::cerr << "                   If it already exists, delete it and try again." << std::endl;
				exit(-1);
			}
			if (mpi_err != MPI_SUCCESS) {
				char error_string[256];
				int length_of_error_string=256;
				MPI_Error_string(mpi_err, error_string, &length_of_error_string);
				fprintf(stderr, "%3d: %s\n", rank, error_string);
			}

			// Write to disk
			if (filebuffer!=NULL) {
				unsigned int w=0;
				while (writeranks[w]!=rank) ++w;
				assert(w<nwriters);
				if (w==nwriters-1)
					assert(filesize-aoffsets[w]==ws);
				mpi_err = MPI_File_iwrite_at(output, aoffsets[w], filebuffer, ws, MPI_CHAR, &request);
				MPI_Wait(&request, &status);
				if (mpi_err != MPI_SUCCESS) {
					char error_string[256];
					int length_of_error_string=256;
					MPI_Error_string(mpi_err, error_string, &length_of_error_string);
					fprintf(stderr, "%3d: %s\n", rank, error_string);
				}
			} else {
				ws = 0; // not a writer
			}

			MPI::COMM_WORLD.Barrier();
			MPI_File_close(&output);
			if (recvrequests!=NULL) {
				delete [] recvrequests;
				recvrequests=NULL;
				delete [] recvstatuses;
				recvstatuses=NULL;
			}
			delete [] misalignments;
			misalignments = NULL;
			delete [] writeranks;
			writeranks=NULL;
			delete [] offsets;
			offsets=NULL;
			delete [] aoffsets;
			aoffsets=NULL;
			delete [] datasizes;
			datasizes=NULL;
			delete [] databuffer;
			databuffer=NULL;
			if (filebuffer!=NULL) {
				delete [] filebuffer;
				filebuffer=NULL;
			}
		}
		#endif
	}

	unsigned long write_buffer(char*& buf) const
	{
		// Find out how big the dataset is
		unsigned long size_in_mem = this->buffer_size();
		unsigned long size_on_disk = 1.125 * size_in_mem + 12;
		#ifdef RAW
		size_on_disk=size_in_mem;
		#endif

		// Figure out the block extents
		unsigned long header_size = 0;
		for (int j=0; j<dim; j++) {
			header_size += static_cast<unsigned long>(sizeof(x0[j]));
			header_size += static_cast<unsigned long>(sizeof(x1[j]));
			header_size += static_cast<unsigned long>(sizeof(b0[j]));
			header_size += static_cast<unsigned long>(sizeof(b1[j]));
		}
		// Make a buffer to hold all the data
		const unsigned long size_of_buffer = header_size + static_cast<unsigned long>(sizeof(size_in_mem))
		                     + size_on_disk + static_cast<unsigned long>(sizeof(size_on_disk));
		buf = new char[size_of_buffer];
		char* dst = buf;
		unsigned long increment=0; // number of bytes to copy

		// Write local limits
		for (int j=0; j<dim; j++) {
			increment = sizeof(x0[j]);
			memcpy(dst, reinterpret_cast<const char*>(&x0[j]), increment);
			dst += increment;
			increment = sizeof(x1[j]);
			memcpy(dst, reinterpret_cast<const char*>(&x1[j]), increment);
			dst += increment;
		}

		// Write local boundary conditions
		for (int j=0; j<dim; j++) {
			increment = sizeof(b0[j]);
			memcpy(dst, reinterpret_cast<const char*>(&b0[j]), increment);
			dst += increment;
			increment = sizeof(b1[j]);
			memcpy(dst, reinterpret_cast<const char*>(&b1[j]), increment);
			dst += increment;
		}

		// Write the size of the raw data block
		increment = sizeof(size_in_mem);
		memcpy(dst, reinterpret_cast<const char*>(&size_in_mem), increment);
		dst += increment;

		// Write the size of the compressed block
		#ifndef RAW
		char* q(dst); // save current location: need to re-write this value later
		#endif
		increment = sizeof(size_on_disk);
		memcpy(dst, reinterpret_cast<const char*>(&size_on_disk), increment);
		dst += increment;

		// Read the data block from grid private data
		#ifdef RAW
		size_on_disk=this->to_buffer(dst);
		#else
		char* raw = new char[size_in_mem];
		size_in_mem = this->to_buffer(raw);

		// Compress the data block to the buffer
		int level=9; // highest compression level (slowest speed)
		int status = compress2(reinterpret_cast<Bytef*>(dst), &size_on_disk, reinterpret_cast<Bytef*>(raw), size_in_mem, level);
		switch (status) {
		case Z_OK:
			break;
		case Z_MEM_ERROR:
			std::cerr << "Compress: out of memory." << std::endl;
			exit(1);
			break;
		case Z_BUF_ERROR:
			std::cerr << "Compress: output buffer wasn't large enough." << std::endl;
			exit(1);
			break;
		}
		assert(size_on_disk<=size_of_buffer); // Abort if data was lost in compression process. Duplicate of Z_BUF_ERROR check.
		dst=NULL;
		delete [] raw;
		raw=NULL;

		// Re-write the size of the (compressed) data block
		increment = sizeof(size_on_disk);
		memcpy(q, reinterpret_cast<const char*>(&size_on_disk), increment);
		#endif

		// Return total size of the populated buffer
		return header_size + static_cast<unsigned long>(sizeof(size_in_mem)) + static_cast<unsigned long>(sizeof(size_on_disk)) + size_on_disk;
	}

	// grid utility functions
	friend int nodes(const grid& GRID)
	{
		return GRID.nodes;
	}
	friend int fields(const grid& GRID)
	{
		return GRID.fields;
	}
	friend int ghosts(const grid& GRID)
	{
		return GRID.ghosts;
	}
	friend int g0(const grid& GRID, int i)
	{
		return GRID.g0[i];
	}
	friend int g1(const grid& GRID, int i)
	{
		return GRID.g1[i];
	}
	friend int b0(const grid& GRID, int i)
	{
		return GRID.b0[i];
	}
	friend int b1(const grid& GRID, int i)
	{
		return GRID.b1[i];
	}
	friend int& b0(grid& GRID, int i)
	{
		return GRID.b0[i];
	}
	friend int& b1(grid& GRID, int i)
	{
		return GRID.b1[i];
	}

	// grid utility functions (all directions)
	friend int x0(const grid& GRID, int i)
	{
		return GRID.x0[i];
	}
	friend int x1(const grid& GRID, int i)
	{
		return GRID.x1[i];
	}
	friend int xmin(const grid& GRID, int i)
	{
		return GRID.x0[i];
	}
	friend int xmax(const grid& GRID, int i)
	{
		return GRID.x1[i];
	}
	friend int xlength(const grid& GRID, int i)
	{
		return GRID.x1[i] - GRID.x0[i];
	}
	friend double dx(const grid& GRID, int i)
	{
		return GRID.dx[i];
	}
	friend double& dx(grid& GRID, int i)
	{
		return GRID.dx[i];
	}
	friend int N0(const grid& GRID, int i)
	{
		return GRID.n0[i];
	}
	friend int N1(const grid& GRID, int i)
	{
		return GRID.n1[i];
	}
	friend int P0(const grid& GRID, int i)
	{
		return GRID.p0[i];
	}
	friend int P1(const grid& GRID, int i)
	{
		return GRID.p1[i];
	}
	friend int sp(const grid& GRID, int i)
	{
		return GRID.sp[i];
	}

	// grid utility functions (x direction)
	friend int x0(const grid& GRID)
	{
		return GRID.x0[0];
	}
	friend int x1(const grid& GRID)
	{
		return GRID.x1[0];
	}
	friend int xmin(const grid& GRID)
	{
		return GRID.x0[0];
	}
	friend int xmax(const grid& GRID)
	{
		return GRID.x1[0];
	}
	friend int xlength(const grid& GRID)
	{
		return GRID.x1[0] - GRID.x0[0];
	}
	friend double dx(const grid& GRID)
	{
		return GRID.dx[0];
	}
	friend double& dx(grid& GRID)
	{
		return GRID.dx[0];
	}

	// grid utility functions (y direction)
	friend int y0(const grid& GRID)
	{
		return GRID.x0[1];
	}
	friend int y1(const grid& GRID)
	{
		return GRID.x1[1];
	}
	friend int ymin(const grid& GRID)
	{
		return GRID.x0[1];
	}
	friend int ymax(const grid& GRID)
	{
		return GRID.x1[1];
	}
	friend int ylength(const grid& GRID)
	{
		return GRID.x1[1] - GRID.x0[1];
	}
	friend double dy(const grid& GRID)
	{
		return GRID.dx[1];
	}
	friend double& dy(grid& GRID)
	{
		return GRID.dx[1];
	}

	// grid utility functions (z direction)
	friend int z0(const grid& GRID)
	{
		return GRID.x0[2];
	}
	friend int z1(const grid& GRID)
	{
		return GRID.x1[2];
	}
	friend int zmin(const grid& GRID)
	{
		return GRID.x0[2];
	}
	friend int zmax(const grid& GRID)
	{
		return GRID.x1[2];
	}
	friend int zlength(const grid& GRID)
	{
		return GRID.x1[2] - GRID.x0[2];
	}
	friend double dz(const grid& GRID)
	{
		return GRID.dx[2];
	}
	friend double& dz(grid& GRID)
	{
		return GRID.dx[2];
	}


	// utility functions
	void swap(grid& GRID)
	{
		// swap grid data
		T* DATA = data;
		data = GRID.data;
		GRID.data = DATA;

		// swap number of nodes
		int NODES = nodes;
		nodes = GRID.nodes;
		GRID.nodes = NODES;

		// swap number of cells
		int CELLS = cells;
		cells = GRID.cells;
		GRID.cells = CELLS;

		// swap number of fields
		int FIELDS = fields;
		fields = GRID.fields;
		GRID.fields = FIELDS;

		// swap number of ghosts
		int GHOSTS = ghosts;
		ghosts = GRID.ghosts;
		GRID.ghosts = GHOSTS;

		// swap grid parameters
		for (int i=0; i<dim; i++) {
			int G0 = g0[i];
			g0[i] = GRID.g0[i];
			GRID.g0[i] = G0;
			int G1 = g1[i];
			g1[i] = GRID.g1[i];
			GRID.g1[i] = G1;
			int X0 = x0[i];
			x0[i] = GRID.x0[i];
			GRID.x0[i] = X0;
			int X1 = x1[i];
			x1[i] = GRID.x1[i];
			GRID.x1[i] = X1;
			int XX = xx[i];
			xx[i] = GRID.xx[i];
			GRID.xx[i] = XX;
			int S0 = s0[i];
			s0[i] = GRID.s0[i];
			GRID.s0[i] = S0;
			int S1 = s1[i];
			s1[i] = GRID.s1[i];
			GRID.s1[i] = S1;
			int SX = sx[i];
			sx[i] = GRID.sx[i];
			GRID.sx[i] = SX;
			int B0 = b0[i];
			b0[i] = GRID.b0[i];
			GRID.b0[i] = B0;
			int B1 = b1[i];
			b1[i] = GRID.b1[i];
			GRID.b1[i] = B1;
			double DX = dx[i];
			dx[i] = GRID.dx[i];
			GRID.dx[i] = DX;
			int P0 = p0[i];
			p0[i] = GRID.p0[i];
			GRID.p0[i] = P0;
			int P1 = p1[i];
			p1[i] = GRID.p1[i];
			GRID.p1[i] = P1;
			int SP = sp[i];
			sp[i] = GRID.sp[i];
			GRID.sp[i] = SP;
			int N0 = n0[i];
			n0[i] = GRID.n0[i];
			GRID.n0[i] = N0;
			int N1 = n1[i];
			n1[i] = GRID.n1[i];
			GRID.n1[i] = N1;
		}
	}

	void copy(const grid& GRID)
	{
		// initialize data
		if (data != NULL) {
			delete [] data;
			data=NULL;
		}

		// copy number of nodes
		nodes = GRID.nodes;

		// copy number of cells
		cells = GRID.cells;

		// copy number of fields
		fields = GRID.fields;

		// copy number of ghosts
		ghosts = GRID.ghosts;

		// copy grid parameters
		for (int i=0; i<dim; i++) {
			g0[i] = GRID.g0[i];
			g1[i] = GRID.g1[i];
			x0[i] = GRID.x0[i];
			x1[i] = GRID.x1[i];
			xx[i] = GRID.xx[i];
			s0[i] = GRID.s0[i];
			s1[i] = GRID.s1[i];
			sx[i] = GRID.sx[i];
			b0[i] = GRID.b0[i];
			b1[i] = GRID.b1[i];
			dx[i] = GRID.dx[i];
			p0[i] = GRID.p0[i];
			p1[i] = GRID.p1[i];
			sp[i] = GRID.sp[i];
			n0[i] = GRID.n0[i];
			n1[i] = GRID.n1[i];
		}

		// resize data structures
		data = new T[cells];
		for (int i=0; i<cells; i++)
			resize(data[i], fields);

		// copy grid data
		for (int i=0; i<cells; i++) {
			unsigned long size_of_buffer = MMSP::buffer_size(data[i]);
			char* buffer = new char[size_of_buffer];
			MMSP::to_buffer(GRID.data[i], buffer);
			MMSP::from_buffer(data[i], buffer);
			delete [] buffer;
			buffer=NULL;
		}
	}


protected:
	T* data;        // local grid data

	int nodes;			// number of nodes (excluding ghosts)
	int cells;			// number of nodes (including ghosts)
	int fields;			// number of fields
	int ghosts;			// ghost cell depth

	int g0[dim];    // global lower coordinate limit (excluding ghosts)
	int g1[dim];    // global upper coordinate limit (excluding ghosts)

	int x0[dim];    // local lower coordinate limit (excluding ghosts)
	int x1[dim];    // local upper coordinate limit (excluding ghosts)
	int xx[dim];    // local cells in slice (excluding ghosts)

	int s0[dim];    // local lower coordinate limit (including ghosts)
	int s1[dim];    // local upper coordinate limit (including ghosts)
	int sx[dim];    // local cells in slice (including ghosts)

	int b0[dim];    // boundary condition at x0
	int b1[dim];    // boundary condition at x1

	double dx[dim]; // global cell spacing

	int p0[dim];
	int p1[dim];
	int sp[dim];    // global processors in slice

	int n0[dim];    // neighbor processor at x0
	int n1[dim];    // neighbor processor at x1
};


// math functions
template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, const vector<int>& x)
{
	T laplacian = 0.0;
	MMSP::vector<int> s = x;
	const T& y = GRID(x);

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s);
		s[i] -= 2;
		const T& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
		laplacian += weight * (yh - 2.0 * y + yl);
	}
	return laplacian;
}

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, const vector<int>& x)
{
	int N = fields(GRID);
	vector<T> laplacian(N, 0.0);
	vector<int> s = x;

	const vector<T>& y = GRID(x);

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const vector<T>& yh = GRID(s);
		s[i] -= 2;
		const vector<T>& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
		for (int j=0; j<N; j++)
			laplacian[j] += weight * (yh[j] - 2.0 * y[j] + yl[j]);
	}
	return laplacian;
}

template<int dim, typename T> T laplacian(const grid<dim,vector<T> >& GRID, const vector<int>& x, const int field)
{
    double laplacian = 0.0;
    vector<int> s = x;

    const T& y = GRID(x)[field];

    for (int i=0; i<dim; i++) {
        s[i] += 1;
        const T& yh = GRID(s)[field];
        s[i] -= 2;
        const T& yl = GRID(s)[field];
        s[i] += 1;

        double weight = 1.0 / pow(dx(GRID, i),2.0);
        laplacian += weight * (yh - 2.0 * y + yl);
    }
    return laplacian;
}

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, const vector<int>& x)
{
	sparse<T> laplacian;
	vector<int> s = x;

	const sparse<T>& y = GRID(x);

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const sparse<T>& yh = GRID(s);
		s[i] -= 2;
		const sparse<T>& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
		laplacian += weight * (yh - 2.0 * y + yl);
	}
	return laplacian;
}

template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, int i)
{
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, int i)
{
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> T laplacian(const grid<dim, vector<T> >& GRID, int i, int f)
{
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x, f);
}

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, int i)
{
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> vector<T> gradient(const grid<dim, T>& GRID, const vector<int>& x)
{
	vector<T> gradient(dim);
	vector<int> s = x;

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s);
		s[i] -= 2;
		const T& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (2.0 * dx(GRID, i));
		gradient[i] = weight * (yh - yl);
	}
	return gradient;
}

template <int dim, typename T> vector<T> grad(const grid<dim, T>& GRID, const vector<int>& x)
{
	return gradient(GRID, x);
}

template <int dim, typename T> T divergence(const grid<dim, T>& GRID, const vector<int>& x)
{
	T divergence = 0.0;
	vector<int> s = x;

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = GRID(s);
		s[i] -= 2;
		const T& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (2.0 * dx(GRID, i));
		divergence += weight * (yh - yl);
	}
	return divergence;
}

template <int dim, typename T> vector<T> divergence(const grid<dim, vector<T> >& GRID, const vector<int>& x)
{
	vector<T> divergence(dim, 0.0);
	vector<int> s = x;

	int N = length(GRID(x));

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const vector<T>& yh = GRID(s);
		s[i] -= 2;
		const vector<T>& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (2.0 * dx(GRID, i));
		for (int j=0; j<N; j++)
			divergence[j] += weight * (yh[j] - yl[j]);
	}
	return divergence;
}

template <int dim, typename T> T div(const grid<dim, T>& GRID, const vector<int>& x)
{
	return divergence(GRID, x);
}

template <int dim, typename T> vector<T> div(const grid<dim, vector<T> >& GRID, const vector<int>& x)
{
	return divergence(GRID, x);
}

// position utility function
template <int dim, typename T> MMSP::vector<int> position(const grid<dim, T>& GRID, int index)
{
	return GRID.position(index);
}

// parallelization
template <int dim, typename T> void ghostswap(grid<dim, T>& GRID)
{
	GRID.ghostswap();
}

template <int dim, typename T> void ghostswap(grid<dim, T>& GRID, const int sublattice)
{
	GRID.ghostswap(sublattice);
}
// buffer I/O functions
template <int dim, typename T> unsigned long buffer_size(const grid<dim, T>& GRID)
{
	return GRID.buffer_size();
}
template <int dim, typename T> unsigned long to_buffer(const grid<dim, T>& GRID, char* buffer)
{
	return GRID.to_buffer(buffer);
}
template <int dim, typename T> unsigned long from_buffer(grid<dim, T>& GRID, const char* buffer)
{
	return GRID.from_buffer(buffer);
}

// file I/O functions
template <int dim, typename T> void read(grid<dim, T>& GRID, std::ifstream& file)
{
	GRID.read(file);
}
template <int dim, typename T> void write(const grid<dim, T>& GRID, std::ifstream& file)
{
	GRID.write(file);
}
template <int dim, typename T> void input(grid<dim, T>& GRID, const char* filename, int GHOSTS = 1, int SINGLE = false)
{
	GRID.input(filename, GHOSTS, SINGLE);
}
template <int dim, typename T> void output(const grid<dim, T>& GRID, const char* filename)
{
	GRID.output(filename);
}
template <int dim, typename T> unsigned long write_buffer(const grid<dim, T>& GRID, char*& buf)
{
	return GRID.write_buffer(buf);
}

// utility functions
template <int dim, typename T> int length(const grid<dim, T>& GRID)
{
	return nodes(GRID);
}
template <int dim, typename T> void resize(int n, grid<dim, T>& GRID) {}
template <int dim, typename T> void swap(grid<dim, T>& GRID1, grid<dim, T>& GRID2)
{
	GRID1.swap(GRID2);
}
template <int dim, typename T> void copy(grid<dim, T>& GRID1, grid<dim, T>& GRID2)
{
	GRID2.copy(GRID1);
}
template <int dim, typename T> std::string name(const grid<dim, T>& GRID)
{
	return std::string("grid:") + name(T());
}

} // namespace MMSP

#endif
