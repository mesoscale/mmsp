// MMSP.vector.h
// Class definition for the MMSP vector data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_VECTOR
#define MMSP_VECTOR

#include<cassert>
#include"MMSP.utility.h"

namespace MMSP {

template <typename T> class vector {
public:
	// constructors / destructor
	vector();

	vector(const vector& v);

	template <typename U> vector(const vector<U>& v);

	vector(int N);

	template <typename U> vector(int N, const U& value);

	~vector();

	// data access operators
	T& operator[](int i) {
		assert(i < size);
		return data[i];
	}
	const T& operator[](int i) const {
		assert(i < size);
		return data[i];
	}

	// assignment operator
	vector& operator=(const T& value);

	vector& operator=(const vector& v);

	template <typename U> vector& operator=(const U& value);

	template <typename U> vector& operator=(const vector<U>& v);

	// buffer I/O functions
	int buffer_size() const;

	int to_buffer(char* buffer) const;

	int from_buffer(const char* buffer);

	// file I/O functions
	void write(std::ofstream& file) const;

	void read(std::ifstream& file);

	// utility functions
	int length() const;

	void resize(int N);

	void copy(const vector& v);

	void swap(vector& v);

	template <typename U> void append(const U& value);

	template <typename U> void append(const vector<U>& v);
private:
	// object data
	T* data;
	int size;
};

// buffer I/O functions
template <typename T> int buffer_size(const vector<T>& v) {
	return v.buffer_size();
}
template <typename T> int to_buffer(const vector<T>& v, char* buffer) {
	return v.to_buffer(buffer);
}
template <typename T> int from_buffer(vector<T>& v, const char* buffer) {
	return v.from_buffer(buffer);
}

// file I/O functions
template <typename T> void write(const vector<T>& v, std::ofstream& file) {
	return v.write(file);
}
template <typename T> void read(vector<T>& v, std::ifstream& file) {
	return v.read(file);
}

// utility functions
template <typename T> int length(const vector<T>& v) {
	return v.length();
}
template <typename T> void resize(vector<T>& v, int n) {
	v.resize(n);
}
template <typename T> void copy(vector<T>& v, const vector<T>& w) {
	v.copy(w);
}
template <typename T> void swap(vector<T>& v, vector<T>& w) {
	v.swap(w);
}
template <typename T, typename U> void append(vector<T>& v, const U& value) {
	v.append(value);
}
template <typename T, typename U> void append(vector<T>& v, const vector<U>& w) {
	v.append(w);
}
template <typename T> std::string name(const vector<T>& s) {
	return std::string("vector:") + name(T());
}

// mathematical operators
template <typename T> vector<T> min(const vector<T>& x, const vector<T>& y) {
	int N = x.length();
	vector<T> z(N);
	for (int i = 0; i < N; i++) z[i] = min(x[i], y[i]);
	return z;
}

template <typename T> vector<T> max(const vector<T>& x, const vector<T>& y) {
	int N = x.length();
	vector<T> z(N);
	for (int i = 0; i < N; i++) z[i] = max(x[i], y[i]);
	return z;
}

template <typename T, typename U> vector<T>& operator+=(vector<T>& x, const vector<U>& y) {
	int N = x.length();
	for (int i = 0; i < N; i++) x[i] += y[i];
	return x;
}
template <typename T, typename U> vector<T> operator+(const vector<T>& x, const vector<U>& y) {
	int N = x.length();
	vector<T> z(N);
	for (int i = 0; i < N; i++) z[i] = x[i] + y[i];
	return z;
}
template <typename T, typename U> vector<T>& operator-=(vector<T>& x, const vector<U>& y) {
	int N = x.length();
	for (int i = 0; i < N; i++) x[i] -= y[i];
	return x;
}
template <typename T, typename U> vector<T> operator-(const vector<T>& x, const vector<U>& y) {
	int N = x.length();
	vector<T> z(N);
	for (int i = 0; i < N; i++) z[i] = x[i] - y[i];
	return z;
}
template <typename T, typename U> vector<T>& operator*=(vector<T>& x, const U& value) {
	int N = x.length();
	for (int i = 0; i < N; i++) x[i] *= static_cast<T>(value);
	return x;
}
template <typename T, typename U> vector<T> operator*(const U& value, const vector<T>& x) {
	int N = x.length();
	vector<T> z(N);
	for (int i = 0; i < N; i++) z[i] = static_cast<T>(value) * x[i];
	return z;
}
template <typename T, typename U> vector<T> operator*(const vector<T>& x, const U& value) {
	return value*x;
}
template <typename T> T operator*(const vector<T>& x, const vector<T>& y) {
	int N = (x.length()>y.length())?y.length():x.length();
	T z(0);
	for (int i = 0; i < N; i++) z += x[i] * y[i];
	return z;
}


// target class: dim = 0 specialization for vector class
template <int ind, typename T>
class target<0, ind, vector<T> > {
public:
	// constructors
	target(vector<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}

	// data access operators
	operator vector<T>&() {
		return *data;
	}
	operator const vector<T>&() const {
		return *data;
	}
	T& operator[](int i) {
		return data->operator[](i);
	}
	const T& operator[](int i) const {
		return data->operator[](i);
	}

	// assignment operators
	vector<T>& operator=(const T& value) const {
		return data->operator=(value);
	}
	vector<T>& operator=(const vector<T>& v) const {
		return data->operator=(v);
	}
	template <typename U> vector<T>& operator=(const U& value) const {
		return data->operator=(value);
	}
	template <typename U> vector<T>& operator=(const vector<U>& v) const {
		return data->operator=(v);
	}

	// buffer I/O functions
	int buffer_size() const {
		return data->buffer_size();
	}
	int to_buffer(char* buffer) const {
		return data->to_buffer(buffer);
	}
	int from_buffer(const char* buffer) const {
		return data->from_buffer(buffer);
	}

	// file I/O functions
	void write(std::ofstream& file) const {
		data->write(file);
	}
	void read(std::ifstream& file) const {
		data->read(file);
	}

	// utility functions
	int length() const {
		return data->length();
	}
	int resize(int n) const {
		return data->resize(n);
	}
	void copy(const target& t) const {
		data->copy(t->data);
	}
	void swap(const target& t) const {
		data->swap(t->data);
	}
	template <typename U> void append(const U& value) const {
		data->append(value);
	}
	template <typename U> void append(const vector<U>& v) const {
		data->append(v);
	}

	// object data
	vector<T>* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

// buffer I/O functions
template <int ind, typename T> int buffer_size(const target<0, ind, vector<T> >& v) {
	return v.buffer_size();
}
template <int ind, typename T> int to_buffer(const target<0, ind, vector<T> >& v, char* buffer) {
	return v.to_buffer(buffer);
}
template <int ind, typename T> int from_buffer(const target<0, ind, vector<T> >& v, const char* buffer) {
	return v.from_buffer(buffer);
}

// file I/O functions
template <int ind, typename T> void write(const target<0, ind, vector<T> >& v, std::ofstream& file) {
	return v.write(file);
}
template <int ind, typename T> void read(const target<0, ind, vector<T> >& v, std::ifstream& file) {
	return v.read(file);
}

// utility functions
template <int ind, typename T> int length(const target<0, ind, vector<T> >& v) {
	return v.length();
}
template <int ind, typename T> void resize(const target<0, ind, vector<T> >& v, int n) {
	v.resize(n);
}
template <int ind, typename T> void copy(const target<0, ind, vector<T> >& v, const target<0, ind, vector<T> >& w) {
	v.copy(w);
}
template <int ind, typename T> void swap(const target<0, ind, vector<T> >& v, const target<0, ind, vector<T> >& w) {
	v.swap(w);
}
template <int ind, typename T, typename U> void append(const target<0, ind, vector<T> >& v, const U& value) {
	v.append(value);
}
template <int ind, typename T, typename U> void append(const target<0, ind, vector<T> >& v, const vector<U>& w) {
	v.append(w);
}
template <int ind, typename T> std::string name(const target<0, ind, vector<T> >& s) {
	return std::string("vector:") + name(T());
}

// mathematical operators
template <int ind, typename T>
vector<T> min(const target<0, ind, vector<T> >& x, const target<0, ind, vector<T> >& y) {
	return min(*(x.data), *(y.data));
}
template <int ind, typename T>
vector<T> max(const target<0, ind, vector<T> >& x, const target<0, ind, vector<T> >& y) {
	return max(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
vector<T>& operator+=(target<0, ind, vector<T> >& x, const target<0, ind, vector<U> >& y) {
	return operator+=(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
vector<T> operator+(const target<0, ind, vector<T> >& x, const target<0, ind, vector<U> >& y) {
	return operator+(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
vector<T>& operator-=(target<0, ind, vector<T> >& x, const target<0, ind, vector<U> >& y) {
	return operator-=(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
vector<T> operator-(const target<0, ind, vector<T> >& x, const target<0, ind, vector<U> >& y) {
	return operator-(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
vector<T>& operator*=(target<0, ind, vector<T> >& x, const U& value) {
	return operator*=(*(x.data), value);
}
template <int ind, typename T, typename U>
vector<T> operator*(const U& value, const target<0, ind, vector<T> >& x) {
	return operator*(value, *(x.data));
}
template <typename T> bool operator==(const vector<T>& a, const vector<T>& b);


} // namespace MMSP

#include "MMSP.vector.cpp"

#endif
