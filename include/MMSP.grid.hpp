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
#include<zlib.h>

#include"MMSP.utility.hpp"
#include"MMSP.scalar.hpp"
#include"MMSP.vector.hpp"
#include"MMSP.sparse.hpp"

namespace MMSP {

// declaration of grid class
template <int dim, typename T> class grid;

// grid utility functions
template <int dim, typename T> int nodes(const grid<dim, T>& GRID) {
	return nodes(GRID);
}
template <int dim, typename T> int fields(const grid<dim, T>& GRID) {
	return fields(GRID);
}
template <int dim, typename T> int ghosts(const grid<dim, T>& GRID) {
	return ghosts(GRID);
}
template <int dim, typename T> int g0(const grid<dim, T>& GRID, int i) {
	return g0(GRID, i);
}
template <int dim, typename T> int g1(const grid<dim, T>& GRID, int i) {
	return g1(GRID, i);
}
template <int dim, typename T> int b0(const grid<dim, T>& GRID, int i) {
	return b0(GRID, i);
}
template <int dim, typename T> int b1(const grid<dim, T>& GRID, int i) {
	return b1(GRID, i);
}
template <int dim, typename T> int& b0(grid<dim, T>& GRID, int i) {
	return b0(GRID, i);
}
template <int dim, typename T> int& b1(grid<dim, T>& GRID, int i) {
	return b1(GRID, i);
}

// grid utility functions (all directions)
template <int dim, typename T> int x0(const grid<dim, T>& GRID, int i) {
	return x0(GRID, i);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID, int i) {
	return x1(GRID, i);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID, int i) {
	return xmin(GRID, i);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID, int i) {
	return xmax(GRID, i);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID, int i) {
	return xlength(GRID, i);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID, int i) {
	return dx(GRID, i);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID, int i) {
	return dx(GRID, i);
}

// grid utility functions (x direction)
template <int dim, typename T> int x0(const grid<dim, T>& GRID) {
	return x0(GRID);
}
template <int dim, typename T> int x1(const grid<dim, T>& GRID) {
	return x1(GRID);
}
template <int dim, typename T> int xmin(const grid<dim, T>& GRID) {
	return xmin(GRID);
}
template <int dim, typename T> int xmax(const grid<dim, T>& GRID) {
	return xmax(GRID);
}
template <int dim, typename T> int xlength(const grid<dim, T>& GRID) {
	return xlength(GRID);
}
template <int dim, typename T> double dx(const grid<dim, T>& GRID) {
	return dx(GRID);
}
template <int dim, typename T> double& dx(grid<dim, T>& GRID) {
	return dx(GRID);
}

// grid utility functions (y direction)
template <int dim, typename T> int y0(const grid<dim, T>& GRID) {
	return y0(GRID);
}
template <int dim, typename T> int y1(const grid<dim, T>& GRID) {
	return y1(GRID);
}
template <int dim, typename T> int ymin(const grid<dim, T>& GRID) {
	return ymin(GRID);
}
template <int dim, typename T> int ymax(const grid<dim, T>& GRID) {
	return ymax(GRID);
}
template <int dim, typename T> int ylength(const grid<dim, T>& GRID) {
	return ylength(GRID);
}
template <int dim, typename T> double dy(const grid<dim, T>& GRID) {
	return dy(GRID);
}
template <int dim, typename T> double& dy(grid<dim, T>& GRID) {
	return dy(GRID);
}

// grid utility functions (z direction)
template <int dim, typename T> int z0(const grid<dim, T>& GRID) {
	return z0(GRID);
}
template <int dim, typename T> int z1(const grid<dim, T>& GRID) {
	return z1(GRID);
}
template <int dim, typename T> int zmin(const grid<dim, T>& GRID) {
	return zmin(GRID);
}
template <int dim, typename T> int zmax(const grid<dim, T>& GRID) {
	return zmax(GRID);
}
template <int dim, typename T> int zlength(const grid<dim, T>& GRID) {
	return zlength(GRID);
}
template <int dim, typename T> double dz(const grid<dim, T>& GRID) {
	return dz(GRID);
}
template <int dim, typename T> double& dz(grid<dim, T>& GRID) {
	return dz(GRID);
}


// instantiation of grid class
template <int dim, typename T>
class grid {
public:
	// constructors
	grid(int FIELDS, int min[dim], int max[dim], int GHOSTS = 1, bool SINGLE = false) {
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i = 0; i < dim; i++) {
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

	grid(int FIELDS, ...) {
		// set number of fields
		fields = FIELDS;

		// read function arguments
		va_list list;
		va_start(list, FIELDS);
		for (int i = 0; i < dim; i++) {
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

	grid(const grid& GRID) {
		// set number of fields
		fields = MMSP::fields(GRID);

		// read function arguments
		for (int i = 0; i < dim; i++) {
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
		for (int i = 0; i < dim; i++) {
			b0[i] = MMSP::b0(GRID, i);
			b1[i] = MMSP::b1(GRID, i);
			dx[i] = MMSP::dx(GRID, i);
		}
	}

	template <typename U>
	grid(const grid<dim, U>& GRID, int FIELDS) {
		// set number of fields
		fields = FIELDS;

		// read function arguments
		for (int i = 0; i < dim; i++) {
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
		for (int i = 0; i < dim; i++) {
			b0[i] = MMSP::b0(GRID, i);
			b1[i] = MMSP::b1(GRID, i);
			dx[i] = MMSP::dx(GRID, i);
		}
	}

	grid(const char* filename, int GHOSTS = 1) {
		// initialize data
		data = NULL;

		// read data from file
		input(filename, GHOSTS);
	}

	void setup(bool SINGLE = false) {
		// setup default grid parameters
		for (int i = 0; i < dim; i++) {
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
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		// if bool SINGLE is set to true,
		// we generate grid data only on proc 0
		if (not SINGLE) {
			// compute integral factors of "np"
			int nfac = 0;
			for (int i = 1; i <= np; i++)
				if ((np / i)*i == np) nfac += 1;

			int* factors = new int[nfac];
			nfac = 0;
			for (int i = 1; i <= np; i++)
				if ((np / i)*i == np) {
					factors[nfac] = i;
					nfac += 1;
				}

			// compute global slice areas
			int area[dim];
			for (int i = 0; i < dim; i++) {
				area[i] = 1;
				for (int j = 0; j < dim; j++)
					if (i != j) area[i] *= (g1[j] - g0[j]);
			}

			// initialize optimal ghost area
			int minarea = -1;

			// compute all combinations of "dim" factors
			int ncom = 1;
			for (int i = 0; i < dim; i++)
				ncom *= nfac;

			for (int i = 0; i < ncom; i++) {
				int combo[dim];
				int total = i;
				for (int j = 0; j < dim; j++) {
					int slice = 1;
					for (int k = j + 1; k < dim; k++)
						slice *= nfac;
					combo[j] = total / slice;
					total -= combo[j] * slice;
				}

				// compute the product of "dim" factors
				int product = 1;
				for (int j = 0; j < dim; j++)
					product *= factors[combo[j]];

				// if product equals "np", compute ghost area
				if (product == np) {
					int test = 0;
					for (int j = 0; j < dim; j++)
						test += area[j] * (1 + factors[combo[j]]);
					// choose optimal (smallest) ghost area
					if (test < minarea or minarea < 0) {
						minarea = test;
						for (int k = 0; k < dim; k++)
							p1[k] = factors[combo[k]];
					}
				}
			}

			// clean up
			delete [] factors;

			// compute slice sizes
			for (int i = 0; i < dim; i++) {
				sp[i] = 1;
				for (int j = i + 1; j < dim; j++)
					sp[i] *= p1[j];
			}

			// determine local grid limits
			int pos[dim];

			int total = id;
			for (int i = 0; i < dim; i++) {
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
			for (int i = 0; i < dim; i++) {
				int npos[dim];
				for (int j = 0; j < dim; j++)
					npos[j] = pos[j];

				// set neighbor below
				n0[i] = 0;
				npos[i] = (pos[i] - 1 + p1[i]) % p1[i];
				for (int j = 0; j < dim; j++)
					n0[i] += sp[j] * npos[j];

				// set neighbor above
				n1[i] = 0;
				npos[i] = (pos[i] + 1) % p1[i];
				for (int j = 0; j < dim; j++)
					n1[i] += sp[j] * npos[j];
			}

			// adjust boundary conditions
			for (int i = 0; i < dim; i++) {
				if (x0[i] != g0[i]) b0[i] = parallel;
				if (x1[i] != g1[i]) b1[i] = parallel;
			}
		}
#endif

		// compute slice sizes
		for (int i = 0; i < dim; i++) {
			s0[i] = x0[i] - ghosts;
			s1[i] = x1[i] + ghosts;
		}
		for (int i = 0; i < dim; i++) {
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
		for (int i = 0; i < cells; i++)
			resize(data[i], fields);

#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
#endif
	}


	// destructor
	~grid() {
		delete [] data;
	}


	// assignment operators
	template <typename U> grid& operator=(const U& value) {
		for (int i = 0; i < cells; i++)
			data[i] = static_cast<T>(value);
	}

	template <typename U> grid& operator=(const grid<dim, U>& GRID) {
		for (int i = 0; i < cells; i++)
			data[i] = static_cast<T>(GRID.data[i]);
	}

	template <typename U> grid& operator+=(const grid<dim, U>& GRID) {
		for (int i = 0; i < cells; i++)
			data[i] += static_cast<T>(GRID.data[i]);
	}

	template <typename U> grid& operator-=(const grid<dim, U>& GRID) {
		for (int i = 0; i < cells; i++)
			data[i] -= static_cast<T>(GRID.data[i]);
	}


	// subscript operators
	target < dim - 1, 0, T > operator [](int x) const {
		check_boundary(x, x0[0], x1[0], b0[0], b1[0]);
		return target < dim - 1, 0, T > (data + (x - s0[0]) * sx[0], s0, sx, x0, x1, b0, b1);
	}

	T& operator()(MMSP::vector<int> x) const {
		T* p = data;
		for (int i = 0; i < dim; i++) {
			check_boundary(x[i], x0[i], x1[i], b0[i], b1[i]);
			p += (x[i] - s0[i]) * sx[i];
		}
		return *p;
	}

	T& operator()(int i) const {
#ifdef MPI_VERSION
		int x[dim];
		for (int j = 0; j < dim; j++) {
			int n = i / xx[j];
			x[j] = n + x0[j];
			i -= n * xx[j];
		}
		i = 0;
		for (int j = 0; j < dim; j++)
			i += (x[j] - s0[j]) * sx[j];
#endif

		return data[i];
	}


	// position utility function
	MMSP::vector<int> position(int i) const {
		MMSP::vector<int> x(dim);
		for (int j = 0; j < dim; j++) {
			int n = i / xx[j];
			x[j] = n + x0[j];
			i -= n * xx[j];
		}
		return x;
	}


	// parallelization
	void ghostswap() {
#ifdef MPI_VERSION
		int np = MPI::COMM_WORLD.Get_size();
		int id = MPI::COMM_WORLD.Get_rank();
		for (int i = 0; i < dim; i++) {
			if (1) {
				// send to processor above and receive from processor below
				int send_proc = n1[i];
				int recv_proc = n0[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j = 0; j < dim; j++) {
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
				if (send_size > std::numeric_limits<int>::max()) {
					std::cout << "Error on rank " << id;
					std::cout << ": send_size is " << send_size;
					std::cout << ", exceeds limit by " << send_size - (unsigned long)std::numeric_limits<int>::max() << std::endl;
					exit(-1);
				}

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 100); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 100); // receive number of ghosts
				MPI::Request::Waitall(2, requests);
				//MPI::COMM_WORLD.Barrier(); // shouldn't be necessary
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				unsigned long buff_size = this->to_buffer(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 200); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 200); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				this->from_buffer(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer
				delete [] send_buffer;
				delete [] recv_buffer;
			}

			if (1) {
				// send to processor below and receive from processor above
				int send_proc = n0[i];
				int recv_proc = n1[i];

				int send_min[dim], send_max[dim];
				int recv_min[dim], recv_max[dim];
				for (int j = 0; j < dim; j++) {
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
				if (send_size > std::numeric_limits<int>::max()) {
					std::cout << "Error on rank " << id;
					std::cout << ": send_size is " << send_size;
					std::cout << ", exceeds limit by " << send_size - (unsigned long)std::numeric_limits<int>::max() << std::endl;
					exit(-1);
				}

				// Small data transfer: blocking Sendrecv should scale -- but don't plan on it.
				MPI::Request requests[2];
				requests[0] = MPI::COMM_WORLD.Isend(&send_size, 1, MPI_UNSIGNED_LONG, send_proc, 300); // send number of ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(&recv_size, 1, MPI_UNSIGNED_LONG, recv_proc, 300); // receive number of incoming ghosts
				MPI::Request::Waitall(2, requests);
				char* send_buffer = new char[send_size];
				char* recv_buffer = new char[recv_size];
				unsigned long buff_size = this->to_buffer(send_buffer, send_min, send_max);

				// Large data transfer requires non-blocking communication
				requests[0] = MPI::COMM_WORLD.Isend(send_buffer, send_size, MPI_CHAR, send_proc, 400); // send ghosts
				requests[1] = MPI::COMM_WORLD.Irecv(recv_buffer, recv_size, MPI_CHAR, recv_proc, 400); // receive ghosts
				MPI::Request::Waitall(2, requests);
				MPI::COMM_WORLD.Barrier();
				this->from_buffer(recv_buffer, recv_min, recv_max); // populate ghost cell data from buffer
				delete [] send_buffer;
				delete [] recv_buffer;
			}
		}
		MPI::COMM_WORLD.Barrier();
#endif
	}


	// buffer I/O
	unsigned long buffer_size() const {
		return buffer_size(x0, x1);
	}

	unsigned long buffer_size(const int min[dim], const int max[dim]) const {
		return buffer_size(data, 0, min, max);
	}

	unsigned long buffer_size(T* p, int i, const int min[dim], const int max[dim]) const {
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x++)
				size += MMSP::buffer_size(*(p + (x - s0[i]) * sx[i]));
		else
			for (int x = min[i]; x < max[i]; x++)
				size += buffer_size(p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long to_buffer(char* buffer) const {
		return to_buffer(buffer, x0, x1);
	}

	unsigned long to_buffer(char* buffer, const int min[dim], const int max[dim]) const {
		return to_buffer(buffer, data, 0, min, max);
	}

	unsigned long to_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) const {
		unsigned long size = 0;
		if (i == dim - 1)
			for (int x = min[i]; x < max[i]; x++)
				size += MMSP::to_buffer(*(p + (x - s0[i]) * sx[i]), buffer + size);
		else
			for (int x = min[i]; x < max[i]; x++)
				size += to_buffer(buffer + size, p + (x - s0[i]) * sx[i], i + 1, min, max);
		return size;
	}

	unsigned long from_buffer(char* buffer) {
		return from_buffer(buffer, x0, x1);
	}

	unsigned long from_buffer(char* buffer, const int min[dim], const int max[dim]) {
		return from_buffer(buffer, data, 0, min, max);
	}

	unsigned long from_buffer(char* buffer, T* p, int i, const int min[dim], const int max[dim]) {
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
	void input(const char* filename, int GHOSTS = 1, int SINGLE = false) {
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
		for (int i = 0; i < dim; i++)
			input >> g0[i] >> g1[i];

		// set number of ghosts
		ghosts = GHOSTS;
#ifndef MPI_VERSION
		ghosts = 0;
#endif

		// setup grid parameters
		delete [] data;
		setup(SINGLE);

		// read cell spacing
		for (int i = 0; i < dim; i++)
			input >> dx[i];

		// ignore trailing endlines
		input.ignore(10, '\n');

		// input grid data
		read(input, GHOSTS);

		// ghostswap if necessary
		if (not SINGLE) ghostswap();
	}

	void output(const char* filename) const {
		int id = 0;
		int np = 1;
#ifdef MPI_VERSION
		id = MPI::COMM_WORLD.Get_rank();
		np = MPI::COMM_WORLD.Get_size();
		std::cout << "\nWARNING: Use output_mpi(filename) instead of output(filename).\n" << std::endl;
#endif
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
		for (int i = 0; i < dim; i++) outstr << g0[i] << " " << g1[i] << '\n';

		// get cell spacing
		for (int i = 0; i < dim; i++) outstr << dx[i] << '\n';

		// Write file header to file
		output.write(outstr.str().c_str(), outstr.str().size());
		// Write number of blocks (processors) to file
		output.write(reinterpret_cast<const char*>(&np), sizeof(np));
#ifdef DEBUG
		std::cout << "Header complete. Writing data.\n";
#endif

		// get grid data to write
		char* buffer;
		unsigned long size = write_buffer(buffer);
		// output grid data
		output.write(buffer, size);
		delete [] buffer;
		buffer = NULL;
	}

	void read(std::ifstream& file, int GHOSTS = 1) {
		// read number of blocks
		int blocks;
		file.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));

		// for each block...
		for (int i = 0; i < blocks; i++) {
			int lmin[dim];
			int lmax[dim];
			// read block limits
			for (int j = 0; j < dim; j++) {
				file.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
				file.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
			}
			int blo[dim];
			int bhi[dim];
			// read boundary conditions
			for (int j = 0; j < dim; j++) {
				file.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
				file.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
			}

			// read block data
			unsigned long size, rawSize;
			file.read(reinterpret_cast<char*>(&rawSize), sizeof(size)); // read raw size
			file.read(reinterpret_cast<char*>(&size), sizeof(size)); // read compressed size
			char* buffer = new char[size];
			file.read(buffer, size);
			grid<dim, T> GRID(fields, lmin, lmax, 0, true);
			// Decompress data
			char* raw = new char[rawSize];
			int status;
			status = uncompress(reinterpret_cast<unsigned char*>(raw), &rawSize, reinterpret_cast<unsigned char*>(buffer), size);
			switch( status ) {
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

			delete [] buffer;
			GRID.from_buffer(raw);
			delete [] raw;

			// find overlapping region
			int min[dim];
			int max[dim];
			bool overlap = true;
			for (int j = 0; j < dim; j++) {
				min[j] = (lmin[j] < x0[j] ? x0[j] : lmin[j]); // if lmin[j]<x0[j] then min[j]=x0[j], else min[j] = lmin[j]
				max[j] = (lmax[j] > x1[j] ? x1[j] : lmax[j]);
				if (min[j] >= max[j]) overlap = false;
			}

			// copy block data that overlaps
			if (overlap) {
				unsigned long size = GRID.buffer_size(min, max);
				char* buffer = new char[size];
				GRID.to_buffer(buffer, min, max);
				this->from_buffer(buffer, min, max);
				delete [] buffer;
			}

			// set boundary conditions from file
			for (int j = 0; j < dim; j++) {
				if (x0[j]==g0[j]) b0[j]=blo[j];
				if (x1[j]==g1[j]) b1[j]=bhi[j];
			}
		}

#ifdef MPI_VERSION
		int id = MPI::COMM_WORLD.Get_rank();
		MPI::COMM_WORLD.Barrier();
#endif
	}

	void output_mpi(const char* filename) const {
#ifndef MPI_VERSION
		std::cerr << "Function output_mpi() requires MPI." << std::endl;
		exit(-1);
#else
		int id = MPI::COMM_WORLD.Get_rank();
		int np = MPI::COMM_WORLD.Get_size();

		// file open error check
		MPI::File output = MPI::File::Open(MPI::COMM_WORLD, filename, MPI::MODE_CREATE | MPI::MODE_WRONLY, MPI::INFO_NULL);
		if (!output) {
			std::cerr << "File output error: could not open " << filename << "." << std::endl;
			exit(-1);
		}

		if (id == 0) {
			std::stringstream outstr;
			// get grid data type
			std::string type = name(*this);
			outstr << type << '\n';

			// get grid dimension
			outstr << dim << '\n';

			// get number of fields
			outstr << fields << '\n';

			// get grid size
			for (int i = 0; i < dim; i++) outstr << g0[i] << " " << g1[i] << '\n';

			// get cell spacing
			for (int i = 0; i < dim; i++) outstr << dx[i] << '\n';

			// Write file header to file
			MPI::Request requests[2];
			requests[0] = output.Iwrite_shared(outstr.str().c_str(), outstr.str().size(), MPI_CHAR);

			// Write number of blocks (processors) to file
			requests[1] = output.Iwrite_shared(reinterpret_cast<const char*>(&np), sizeof(np), MPI_CHAR);

			MPI::Request::Waitall(2, requests);
			#ifdef PDEBUG
			std::cout<<"\nWrote MMSP header on rank "<<id<<" to "<<filename<<std::endl;
			#endif
		}
		MPI::COMM_WORLD.Barrier();
		output.Sync();
		MPI::COMM_WORLD.Barrier();
		output.Sync();


		// get grid data to write -- compressed with zlib
		char* buffer;
		unsigned long size = write_buffer(buffer);
		#ifdef PDEBUG
		int vote=1;
		int total_procs=0;
		MPI::COMM_WORLD.Barrier();
		MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
		if (id==0) std::cout<<"\nCompressed MMSP grid data on "<<total_procs<<" ranks."<<std::endl;
		#endif

		// Write the buffer to disk!
		MPI::Status status;
		MPI::Request request;

		output.Sync();
		MPI::COMM_WORLD.Barrier();
		output.Sync();

		//output.Write_ordered(buffer,size,MPI_CHAR,status);
		request = output.Iwrite_shared(buffer, size, MPI_CHAR);

		// After non-blocking shared write, wait for file IO to complete
		request.Wait();
		MPI::COMM_WORLD.Barrier();
		#ifdef PDEBUG
		total_procs=0;
		MPI::COMM_WORLD.Allreduce(&vote, &total_procs, 1, MPI_INT, MPI_SUM);
		if (id==0) std::cout<<"\nWrote MMSP grid data from "<<total_procs<<" ranks to "<<filename<<std::endl;
		#endif
		delete [] buffer;
		buffer = NULL;

		// Make sure everything's written before closing the file.
		output.Sync();
		MPI::COMM_WORLD.Barrier();
		output.Sync();
		output.Close();

		// Make sure everything's written before closing the file.
		#ifdef PDEBUG
		if (id==0) std::cout<<"\nFinished writing "<<filename<<'\n'<<std::endl;
		#endif
		std::cout << std::flush;
		MPI::COMM_WORLD.Barrier();
#endif
	}

	unsigned long write_buffer(char* &buf) const {
		int id = 0;
		int np = 1;
#ifdef MPI_VERSION
		id = MPI::COMM_WORLD.Get_rank();
		np = MPI::COMM_WORLD.Get_size();
#endif

		// Find out how big the dataset is
		unsigned long data_size = this->buffer_size();
		unsigned long compressed_size = 1.125 * data_size + 12;

		// write number of blocks
		int blocks = np;

		// Figure out the block extents
		unsigned long header_size = 0;
		for (int j = 0; j < dim; j++) {
			header_size += static_cast<unsigned long>(sizeof(x0[j]));
			header_size += static_cast<unsigned long>(sizeof(x1[j]));
			header_size += static_cast<unsigned long>(sizeof(b0[j]));
			header_size += static_cast<unsigned long>(sizeof(b1[j]));
		}
		// Make a buffer to hold all the data
		unsigned long size = header_size + static_cast<unsigned long>(sizeof(data_size))
		                     + static_cast<unsigned long>(sizeof(compressed_size)) + compressed_size;
		buf = new char[size];
		for (unsigned long i = 0; i < size; ++i) buf[i] = 0;
		char* p = buf;
		int increment = 0; // number of bytes to copy (1 char per byte)

		// Write local limits
		for (int j = 0; j < dim; j++) {
			increment = sizeof(x0[j]);
			memcpy(p, reinterpret_cast<const char*>(&x0[j]), increment);
			p += increment;
			increment = sizeof(x1[j]);
			memcpy(p, reinterpret_cast<const char*>(&x1[j]), increment);
			p += increment;
		}

		// Write local boundary conditions
		for (int j = 0; j < dim; j++) {
			increment = sizeof(b0[j]);
			memcpy(p, reinterpret_cast<const char*>(&b0[j]), increment);
			p += increment;
			increment = sizeof(b1[j]);
			memcpy(p, reinterpret_cast<const char*>(&b1[j]), increment);
			p += increment;
		}

		// Write the size of the raw data block
		increment = sizeof(data_size);
		memcpy(p, reinterpret_cast<const char*>(&data_size), increment);
		p += increment;
		char* q(p); // save current location: need to re-write this value later

		// Write the size of the compressed block
		increment = sizeof(compressed_size);
		memcpy(p, reinterpret_cast<const char*>(&compressed_size), increment);
		p += increment;

		// Read the data block
		char* raw = new char[data_size];
		unsigned long xfer_size = this->to_buffer(raw);
		if (xfer_size > size - header_size) {
			std::cout << "ERROR: Rank" << std::setw(2) << std::right << id << " buffered " << xfer_size << " B instead of " << data_size << " B." << std::endl;
			exit(1);
		}

		// Compress the data block to the buffer
		int status;
		int level = 9; // highest compression level (slowest speed)
		status = compress2(reinterpret_cast<Bytef*>(p), &compressed_size, reinterpret_cast<Bytef*>(raw), data_size, level);
		switch( status ) {
		case Z_OK:
			break;

		case Z_MEM_ERROR:
			std::cerr << "Compress: out of memory." << std::endl;
			exit(1);    // quit.
			break;

		case Z_BUF_ERROR:
			std::cerr << "Compress: output buffer wasn't large enough." << std::endl;
			exit(1);    // quit.
			break;
		}

		// Re-write the size of the (compressed) data block
		increment = sizeof(compressed_size);
		memcpy(q, reinterpret_cast<const char*>(&compressed_size), increment);

		// Update size of the write-buffer
		size = header_size + static_cast<unsigned long>(sizeof(data_size)) + static_cast<unsigned long>(sizeof(compressed_size)) + compressed_size;

		// Cleanup
		delete [] raw;
		return size;
	}

	// grid utility functions
	friend int nodes(const grid& GRID) {
		return GRID.nodes;
	}
	friend int fields(const grid& GRID) {
		return GRID.fields;
	}
	friend int ghosts(const grid& GRID) {
		return GRID.ghosts;
	}
	friend int g0(const grid& GRID, int i) {
		return GRID.g0[i];
	}
	friend int g1(const grid& GRID, int i) {
		return GRID.g1[i];
	}
	friend int b0(const grid& GRID, int i) {
		return GRID.b0[i];
	}
	friend int b1(const grid& GRID, int i) {
		return GRID.b1[i];
	}
	friend int& b0(grid& GRID, int i) {
		return GRID.b0[i];
	}
	friend int& b1(grid& GRID, int i) {
		return GRID.b1[i];
	}

	// grid utility functions (all directions)
	friend int x0(const grid& GRID, int i) {
		return GRID.x0[i];
	}
	friend int x1(const grid& GRID, int i) {
		return GRID.x1[i];
	}
	friend int xmin(const grid& GRID, int i) {
		return GRID.x0[i];
	}
	friend int xmax(const grid& GRID, int i) {
		return GRID.x1[i];
	}
	friend int xlength(const grid& GRID, int i) {
		return GRID.x1[i] - GRID.x0[i];
	}
	friend double dx(const grid& GRID, int i) {
		return GRID.dx[i];
	}
	friend double& dx(grid& GRID, int i) {
		return GRID.dx[i];
	}
	friend int N0(const grid& GRID, int i) {
		return GRID.n0[i];
	}
	friend int N1(const grid& GRID, int i) {
		return GRID.n1[i];
	}

	// grid utility functions (x direction)
	friend int x0(const grid& GRID) {
		return GRID.x0[0];
	}
	friend int x1(const grid& GRID) {
		return GRID.x1[0];
	}
	friend int xmin(const grid& GRID) {
		return GRID.x0[0];
	}
	friend int xmax(const grid& GRID) {
		return GRID.x1[0];
	}
	friend int xlength(const grid& GRID) {
		return GRID.x1[0] - GRID.x0[0];
	}
	friend double dx(const grid& GRID) {
		return GRID.dx[0];
	}
	friend double& dx(grid& GRID) {
		return GRID.dx[0];
	}

	// grid utility functions (y direction)
	friend int y0(const grid& GRID) {
		return GRID.x0[1];
	}
	friend int y1(const grid& GRID) {
		return GRID.x1[1];
	}
	friend int ymin(const grid& GRID) {
		return GRID.x0[1];
	}
	friend int ymax(const grid& GRID) {
		return GRID.x1[1];
	}
	friend int ylength(const grid& GRID) {
		return GRID.x1[1] - GRID.x0[1];
	}
	friend double dy(const grid& GRID) {
		return GRID.dx[1];
	}
	friend double& dy(grid& GRID) {
		return GRID.dx[1];
	}

	// grid utility functions (z direction)
	friend int z0(const grid& GRID) {
		return GRID.x0[2];
	}
	friend int z1(const grid& GRID) {
		return GRID.x1[2];
	}
	friend int zmin(const grid& GRID) {
		return GRID.x0[2];
	}
	friend int zmax(const grid& GRID) {
		return GRID.x1[2];
	}
	friend int zlength(const grid& GRID) {
		return GRID.x1[2] - GRID.x0[2];
	}
	friend double dz(const grid& GRID) {
		return GRID.dx[2];
	}
	friend double& dz(grid& GRID) {
		return GRID.dx[2];
	}


	// utility functions
	void swap(grid& GRID) {
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
		for (int i = 0; i < dim; i++) {
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

	void copy(const grid& GRID) {
		// initialize data
		if (data != NULL) delete [] data;

		// copy number of nodes
		nodes = GRID.nodes;

		// copy number of cells
		cells = GRID.cells;

		// copy number of fields
		fields = GRID.fields;

		// copy number of ghosts
		ghosts = GRID.ghosts;

		// copy grid parameters
		for (int i = 0; i < dim; i++) {
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
		for (int i = 0; i < cells; i++)
			resize(data[i], fields);

		// copy grid data
		for (int i = 0; i < cells; i++) {
			int size = MMSP::buffer_size(data[i]);
			char* buffer = new char[size];
			MMSP::to_buffer(GRID.data[i], buffer);
			MMSP::from_buffer(data[i], buffer);
			delete [] buffer;
		}
	}


protected:
	T* data;        // local grid data

	int nodes;      // number of nodes (excluding ghosts)
	int cells;      // number of nodes (including ghosts)
	int fields;     // number of fields
	int ghosts;     // ghost cell depth

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
template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, const vector<int>& x) {
	T laplacian = 0.0;
	MMSP::vector<int> s = x;
	const T& y = GRID(x);

	for (int i = 0; i < dim; i++) {
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

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, const vector<int>& x) {
	int N = fields(GRID);
	vector<T> laplacian(N, 0.0);
	vector<int> s = x;

	const vector<T>& y = GRID(x);

	for (int i = 0; i < dim; i++) {
		s[i] += 1;
		const vector<T>& yh = GRID(s);
		s[i] -= 2;
		const vector<T>& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
		for (int j = 0; j < N; j++)
			laplacian[j] += weight * (yh[j] - 2.0 * y[j] + yl[j]);
	}
	return laplacian;
}

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, const vector<int>& x) {
	sparse<T> laplacian;
	vector<int> s = x;

	const sparse<T>& y = GRID(x);

	for (int i = 0; i < dim; i++) {
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

template <int dim, typename T> T laplacian(const grid<dim, T>& GRID, int i) {
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> vector<T> laplacian(const grid<dim, vector<T> >& GRID, int i) {
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> sparse<T> laplacian(const grid<dim, sparse<T> >& GRID, int i) {
	vector<int> x = GRID.position(i);
	return laplacian(GRID, x);
}

template <int dim, typename T> vector<T> gradient(const grid<dim, T>& GRID, const vector<int>& x) {
	vector<T> gradient(dim);
	vector<int> s = x;

	for (int i = 0; i < dim; i++) {
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

template <int dim, typename T> vector<T> grad(const grid<dim, T>& GRID, const vector<int>& x) {
	return gradient(GRID, x);
}

template <int dim, typename T> T divergence(const grid<dim, T>& GRID, const vector<int>& x) {
	T divergence = 0.0;
	vector<int> s = x;

	for (int i = 0; i < dim; i++) {
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

template <int dim, typename T> vector<T> divergence(const grid<dim, vector<T> >& GRID, const vector<int>& x) {
	vector<T> divergence(dim, 0.0);
	vector<int> s = x;

	int N = length(GRID(x));

	for (int i = 0; i < dim; i++) {
		s[i] += 1;
		const vector<T>& yh = GRID(s);
		s[i] -= 2;
		const vector<T>& yl = GRID(s);
		s[i] += 1;

		double weight = 1.0 / (2.0 * dx(GRID, i));
		for (int j = 0; j < N; j++)
			divergence[j] += weight * (yh[j] - yl[j]);
	}
	return divergence;
}

template <int dim, typename T> T div(const grid<dim, T>& GRID, const vector<int>& x) {
	return divergence(GRID, x);
}

template <int dim, typename T> vector<T> div(const grid<dim, vector<T> >& GRID, const vector<int>& x) {
	return divergence(GRID, x);
}

// position utility function
template <int dim, typename T> MMSP::vector<int> position(const grid<dim, T>& GRID, int index) {
	return GRID.position(index);
}

// parallelization
template <int dim, typename T> void ghostswap(grid<dim, T>& GRID) {
	GRID.ghostswap();
}

// buffer I/O functions
template <int dim, typename T> unsigned long buffer_size(const grid<dim, T>& GRID) {
	return GRID.buffer_size();
}
template <int dim, typename T> unsigned long to_buffer(const grid<dim, T>& GRID, char* buffer) {
	return GRID.to_buffer(buffer);
}
template <int dim, typename T> unsigned long from_buffer(grid<dim, T>& GRID, const char* buffer) {
	return GRID.from_buffer(buffer);
}

// file I/O functions
template <int dim, typename T> void read(grid<dim, T>& GRID, std::ifstream& file) {
	GRID.read(file);
}
template <int dim, typename T> void write(const grid<dim, T>& GRID, std::ifstream& file) {
	GRID.write(file);
}
template <int dim, typename T> void input(grid<dim, T>& GRID, const char* filename, int GHOSTS = 1, int SINGLE = false) {
	GRID.input(filename, GHOSTS, SINGLE);
}
template <int dim, typename T> void output(const grid<dim, T>& GRID, const char* filename) {
	GRID.output(filename);
}
template <int dim, typename T> unsigned long write_buffer(const grid<dim, T>& GRID, char* &buf) {
	GRID.write_buffer(buf);
}
#ifdef MPI_VERSION
template <int dim, typename T> void write_mpi(const grid<dim, T>& GRID, MPI::File& file) {
	GRID.write_mpi(file);
}
template <int dim, typename T> void output_mpi(const grid<dim, T>& GRID, const char* filename) {
	GRID.output_mpi(filename);
}
#endif

// utility functions
template <int dim, typename T> int length(const grid<dim, T>& GRID) {
	return nodes(GRID);
}
template <int dim, typename T> void resize(int n, grid<dim, T>& GRID) {}
template <int dim, typename T> void swap(grid<dim, T>& GRID1, grid<dim, T>& GRID2) {
	GRID1.swap(GRID2);
}
template <int dim, typename T> void copy(grid<dim, T>& GRID1, grid<dim, T>& GRID2) {
	GRID2.copy(GRID1);
}
template <int dim, typename T> std::string name(const grid<dim, T>& GRID) {
	return std::string("grid:") + name(T());
}

} // namespace MMSP
#endif
