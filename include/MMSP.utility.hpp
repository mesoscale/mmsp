// MMSP.utility.hpp
// Utility functions and classes for MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_UTILITY
#define MMSP_UTILITY
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<ciso646>
#include<ctime>

namespace MMSP {

// MMSP Init function
void Init(int argc, char* argv[]) {
#ifdef MPI_VERSION
	MPI::Init(argc, argv);
#endif
}

// MMSP Finalize function
void Finalize() {
#ifdef MPI_VERSION
	MPI::Finalize();
#endif
}


// MMSP boundary conditions
enum {
	mirror    = 0,
	Neumann   = 1,
	periodic  = 2,
	parallel  = 3,
	Dirichlet = 4
};

// check_boundary: a utility function that adjusts coordinates
// based on limiting coordinates and boundary conditions
void check_boundary(int& x, int x0, int x1, int b0, int b1) {
	if (x < x0) {
		if (b0 == Neumann or b0 == Dirichlet) x = x0;
#ifndef MPI_VERSION
		else if (b0 == periodic) x = x1 - (x0 - x);
#endif
		else if (b0 == mirror) x = 2 * x0 - x;
	} else if (x >= x1) {
		if (b1 == Neumann or b1 == Dirichlet) x = (x1 - 1);
#ifndef MPI_VERSION
		else if (b1 == periodic) x = x0 + (x - x1);
#endif
		else if (b1 == mirror) x = 2 * (x1 - 1) - x;
	}
}


// target class: a utility class that makes
// grids with array-style subscripting possible
template <int dim, int index, typename T>
class target {
public:
	target(T* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}
	target < dim - 1, index + 1, T > operator[](int x) {
		int i = index + 1;
		check_boundary(x, x0[i], x1[i], b0[i], b1[i]);
		return target < dim - 1, index + 1, T > (data + (x - s0[i]) * sx[i], s0, sx, x0, x1, b0, b1);
	}
	const target < dim - 1, index + 1, T > operator[](int x) const {
		int i = index + 1;
		check_boundary(x, x0[i], x1[i], b0[i], b1[i]);
		return target < dim - 1, index + 1, T > (data + (x - s0[i]) * sx[i], s0, sx, x0, x1, b0, b1);
	}
	T* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

// target class: dim = 1 specialization for all data types
template <int index, typename T>
class target<1, index, T> {
public:
	target(T* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}
	T& operator[](int x) {
		int i = index + 1;
		check_boundary(x, x0[i], x1[i], b0[i], b1[i]);
		return *(data + (x - s0[i]) * sx[i]);
	}
	const T& operator[](int x) const {
		int i = index + 1;
		check_boundary(x, x0[i], x1[i], b0[i], b1[i]);
		return *(data + (x - s0[i]) * sx[i]);
	}
	T* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

// target class: dim = 0 specialization for primitive data types
// this class template must be overloaded for new classes, but
// this is only necessary if grids with dim = 1 are used
template <int index, typename T>
class target<0, index, T> {
public:
	target(T* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}
	template <typename U> T& operator=(const U& value) {
		*data = value;
		return *data;
	}
	template <typename U> const T& operator=(const U& value) const {
		*data = value;
		return *data;
	}
	operator T&() {
		return *data;
	}
	operator const T&() const {
		return *data;
	}
	T* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};


// buffer I/O functions for primitive and POD data types
// overload as friend functions within new non-POD class definitions
template <typename T> unsigned long buffer_size(const T& value) {
	return sizeof(value);
}
template <typename T> unsigned long to_buffer(const T& value, char* buffer) {
	memcpy(buffer, &value, sizeof(value));
	return sizeof(value);
}
template <typename T> unsigned long from_buffer(T& value, const char* buffer) {
	memcpy(&value, buffer, sizeof(value));
	return sizeof(value);
}

// file I/O functions for primitive and POD data types
// overload as friend functions within new non-POD class definitions
template <typename T> void read(T& value, std::ifstream& file) {
	file.read(reinterpret_cast<char*>(&value), sizeof(value));
}
template <typename T> void write(const T& value, std::ofstream& file) {
	file.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

// utility functions for primitive and POD data types
// overload as friend functions within new non-POD class definitions
template <typename T> int length(const T& value) {
	return 1;
}
template <typename T> void resize(T& value, int n) {
	return;
}
template <typename T> void copy(T& value, const T& copy) {
	value = copy;
}
template <typename T> void swap(T& value, T& swap) {
	T temp = value;
	value = swap;
	swap = temp;
}

// status function for type identification
std::string name(const bool& value) {
	return "bool";
}
std::string name(const char& value) {
	return "char";
}
std::string name(const unsigned char& value) {
	return "unsigned char";
}
std::string name(const int& value) {
	return "int";
}
std::string name(const unsigned int& value) {
	return "unsigned int";
}
std::string name(const long& value) {
	return "long";
}
std::string name(const unsigned long& value) {
	return "unsigned long";
}
std::string name(const short& value) {
	return "short";
}
std::string name(const unsigned short& value) {
	return "unsigned short";
}
std::string name(const float& value) {
	return "float";
}
std::string name(const double& value) {
	return "double";
}
std::string name(const long double& value) {
	return "long double";
}

// mathematical operations
template <typename T> T max(const T& x, const T& y) {
	return x > y ? x : y;
}
template <typename T> T min(const T& x, const T& y) {
	return x < y ? x : y;
}

// global reducing function
template <typename T> T global(T& value, const char* operation) {
	// initialize global value
	T global = value;

#ifdef MPI_VERSION
	int rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();

	if (rank == 0) {
		// receive local values
		for (int i = 1; i < np; i++) {
			T temp;
			int size;
			MPI::COMM_WORLD.Recv(&size, 1, MPI_INT, i, 100);
			char* buffer = new char[size];
			MPI::COMM_WORLD.Recv(buffer, size, MPI_CHAR, i, 200);
			from_buffer(temp, buffer);
			if (buffer != NULL) {
				delete [] buffer;
				buffer=NULL;
			}

			// perform operation
			if (std::string(operation)=="add" or std::string(operation)=="sum")
				global += temp;
			else if (std::string(operation)=="min" or std::string(operation)=="minimum")
				global = min(global, temp);
			else if (std::string(operation)=="max" or std::string(operation)=="maximum")
				global = max(global, temp);
		}

		// send global value
		for (int i = 1; i < np; i++) {
			int size = buffer_size(global);
			MPI::COMM_WORLD.Send(&size, 1, MPI_INT, i, 300);
			char* buffer = new char[size];
			to_buffer(global, buffer);
			MPI::COMM_WORLD.Send(buffer, size, MPI_CHAR, i, 400);
			if (buffer != NULL) {
				delete [] buffer;
				buffer=NULL;
			}
		}
	}

	else {
		// send local value
		int size = buffer_size(value);
		MPI::COMM_WORLD.Send(&size, 1, MPI_INT, 0, 100);
		char* buffer = new char[size];
		to_buffer(value, buffer);
		MPI::COMM_WORLD.Send(buffer, size, MPI_CHAR, 0, 200);
		if (buffer != NULL) {
			delete [] buffer;
			buffer=NULL;
		}

		// receive global value
		MPI::COMM_WORLD.Recv(&size, 1, MPI_INT, 0, 300);
		buffer = new char[size];
		MPI::COMM_WORLD.Recv(buffer, size, MPI_CHAR, 0, 400);
		from_buffer(global, buffer);
		if (buffer != NULL) {
			delete [] buffer;
			buffer=NULL;
		}
	}

	MPI::COMM_WORLD.Barrier();
#endif

	return global;
}

} // namespace MMSP

/*
	Prints timestamps and a 20-point progress bar to stdout.
	Call once inside the update function (or equivalent).

	for (int step=0; step<steps; step++) {
		print_progress(step, steps);
		...
		for (int n=0; n<nodes(grid); n++) {
			...
		}
	}
*/
void print_progress(const int step, const int steps) {
  char* timestring;
  static unsigned long tstart;
  struct tm* timeinfo;
  static int iterations = 0;

  if (step==0) {
    tstart = time(NULL);
    std::time_t rawtime;
    std::time( &rawtime );
    timeinfo = std::localtime( &rawtime );
    timestring = std::asctime(timeinfo);
    timestring[std::strlen(timestring)-1] = '\0';
    std::cout<<"No. "<<1+iterations/steps<<":\t"<<timestring<<" ["<<std::flush;
  } else if (step==steps-1) {
    unsigned long deltat = time(NULL)-tstart;
    std::cout<<"•] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
                    <<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
                    <<std::setw(2)<<std::right<<deltat%60<<"s"<<std::endl;
  } else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
  iterations++;
}
#endif
