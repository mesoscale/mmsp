// MMSP.utility.h
// Utility function and class declarations for MMSP
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_UTILITY
#define MMSP_UTILITY

#include<cstring>
#include<fstream>

namespace MMSP {

// MMSP Init function
void Init(int argc, char* argv[]);

// MMSP Finalize function
void Finalize();

// MMSP Abort function
void Abort(int err);

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
void check_boundary(int& x, int x0, int x1, int b0, int b1);


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
template <typename T> T global(T& value, const char* operation);

// simple progress bar for the terminal
void print_progress(const int step, const int steps);


} // namespace MMSP

#include "MMSP.utility.cpp"

#endif
