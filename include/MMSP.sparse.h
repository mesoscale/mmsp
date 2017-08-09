// MMSP.sparse.h
// Class definition for the MMSP sparse data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_SPARSE
#define MMSP_SPARSE

#include "MMSP.utility.h"

namespace MMSP {

template <typename T>
struct item {
	int index;
	T value;
};

template <typename T>
class sparse {
public:
	// constructors / destructor
	sparse();
	sparse(const sparse& x);
	~sparse();

	// assignment operator
	sparse& operator=(const sparse& x);
	template <typename U> sparse& operator=(const sparse<U>& x);

	// data access operators
	T& set(int index);

	T operator[](int index) const;

	unsigned int grain_id() const;
	double getMagPhi() const;

	int index(const int& i) const;
	T value(const int& i) const;

	// buffer I/O functions
	int buffer_size() const;
	int to_buffer(char* buffer) const;
	int from_buffer(const char* buffer);

	// file I/O functions
	void write(std::ofstream& file) const;
	void read(std::ifstream& file);

	// utility functions
	int length() const;
	void resize(int n);
	void copy(const sparse& s);
	void swap(sparse& s);

private:
	// object data
	item<T>* data;
	int size;
};

// buffer I/O functions
template <typename T> int buffer_size(const sparse<T>& s) {
	return s.buffer_size();
}
template <typename T> int to_buffer(const sparse<T>& s, char* buffer) {
	return s.to_buffer(buffer);
}
template <typename T> int from_buffer(sparse<T>& s, const char* buffer) {
	return s.from_buffer(buffer);
}

// file I/O functions
template <typename T> void write(const sparse<T>& s, std::ofstream& file) {
	return s.write(file);
}
template <typename T> void read(sparse<T>& s, std::ifstream& file) {
	return s.read(file);
}

// utility functions
template <typename T> int length(const sparse<T>& s) {
	return s.length();
}
template <typename T> void resize(sparse<T>& s, int n) {
	s.resize(n);
}
template <typename T> void copy(sparse<T>& s, const sparse<T>& t) {
	s.copy(t);
}
template <typename T> void swap(sparse<T>& s, sparse<T>& t) {
	s.swap(t);
}

template <typename T> T& set(sparse<T>& s, int index) {
	return s.set(index);
}
template <typename T> int index(const sparse<T>& s, int i) {
	return s.index(i);
}
template <typename T> T value(const sparse<T>& s, int i) {
	return s.value(i);
}
template <typename T> std::string name(const sparse<T>& s) {
	return std::string("sparse:") + name(T());
}

// numerical operators
template <typename T> sparse<T> max(const sparse<T>& x, const sparse<T>& y);

template <typename T> sparse<T> min(const sparse<T>& x, const sparse<T>& y);

template <typename T, typename U> sparse<T>& operator+=(sparse<T>& x, const sparse<U>& y);

template <typename T, typename U> sparse<T> operator+(const sparse<T>& x, const sparse<U>& y);

template <typename T, typename U> sparse<T>& operator-=(sparse<T>& x, const sparse<U>& y);

template <typename T, typename U> sparse<T> operator-(const sparse<T>& x, const sparse<U>& y);

template <typename T, typename U> sparse<T>& operator*=(sparse<T>& x, const U& value);

template <typename T, typename U> sparse<T> operator*(const U& value, const sparse<T>& x);


// target class: dim = 0 specialization for sparse class
template <int ind, typename T>
class target<0, ind, sparse<T> > {
public:
	target(sparse<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1);

	operator sparse<T>&() {
		return *data;
	}
	operator const sparse<T>&() const {
		return *data;
	}
	T operator[](int i) {
		return data->operator[](i);
	}
	const T operator[](int i) const {
		return data->operator[](i);
	}

	sparse<T>& operator=(const sparse<T>& s) const {
		return data->operator=(s);
	}
	template <typename U> sparse<T>& operator=(const sparse<U>& s) const {
		return data->operator=(s);
	}

	int buffer_size() const {
		return data->buffer_size();
	}
	int to_buffer(char* buffer) const {
		return data->to_buffer(buffer);
	}
	int from_buffer(const char* buffer) const {
		return data->from_buffer(buffer);
	}

	void write(std::ofstream& file) const {
		data->write(file);
	}
	void read(std::ifstream& file) const {
		data->read(file);
	}

	int length() const {
		return data->length();
	}
	int resize(int n) const {
		return data->resize(n);
	}
	void copy(const target& t) const {
		return data->copy(t->data);
	}
	void swap(const target& t) const {
		return data->swap(t->data);
	}

	T& set(int index) const {
		return data->set(index);
	}
	int index(int i) const {
		return data->index(i);
	}
	T value(int i) const {
		return data->value(i);
	}

	sparse<T>* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

template <int ind, typename T> T& set(const target<0, ind, sparse<T> >& s, int index) {
	return s.set(index);
}
template <int ind, typename T> int index(const target<0, ind, sparse<T> >& s, int i) {
	return s.index(i);
}
template <int ind, typename T> T value(const target<0, ind, sparse<T> >& s, int i) {
	return s.value(i);
}

// buffer I/O functions
template <int ind, typename T> int buffer_size(const target<0, ind, sparse<T> >& s) {
	return s.buffer_size();
}
template <int ind, typename T> int to_buffer(const target<0, ind, sparse<T> >& s, char* buffer) {
	return s.to_buffer(buffer);
}
template <int ind, typename T> int from_buffer(const target<0, ind, sparse<T> >& s, const char* buffer) {
	return s.from_buffer(buffer);
}

// file I/O functions
template <int ind, typename T> void write(const target<0, ind, sparse<T> >& s, std::ofstream& file) {
	return s.write(file);
}
template <int ind, typename T> void read(const target<0, ind, sparse<T> >& s, std::ifstream& file) {
	return s.read(file);
}

// utility functions
template <int ind, typename T> int length(const target<0, ind, sparse<T> >& s) {
	return s.length();
}
template <int ind, typename T> void resize(const target<0, ind, sparse<T> >& s, int n) {
	s.resize(n);
}
template <int ind, typename T> void copy(const target<0, ind, sparse<T> >& s, const target<0, ind, sparse<T> >& t) {
	s.copy(t);
}
template <int ind, typename T> void swap(const target<0, ind, sparse<T> >& s, const target<0, ind, sparse<T> >& t) {
	s.swap(t);
}
template <int ind, typename T> std::string name(const target<0, ind, sparse<T> >& s) {
	return std::string("sparse:") + name(T());
}

// numerical operators
template <int ind, typename T>
sparse<T>& min(target<0, ind, sparse<T> > x, const target<0, ind, sparse<T> >& y) {
	return min(*(x.data), *(y.data));
}
template <int ind, typename T>
sparse<T>& max(target<0, ind, sparse<T> > x, const target<0, ind, sparse<T> >& y) {
	return max(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
sparse<T>& operator+=(target<0, ind, sparse<T> > x, const target<0, ind, sparse<U> >& y) {
	return operator+=(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
sparse<T> operator+(const target<0, ind, sparse<T> > x, const target<0, ind, sparse<U> >& y) {
	return operator+(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
sparse<T>& operator-=(target<0, ind, sparse<T> > x, const target<0, ind, sparse<U> >& y) {
	return operator-=(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
sparse<T> operator-(const target<0, ind, sparse<T> > x, const target<0, ind, sparse<U> >& y) {
	return operator-(*(x.data), *(y.data));
}
template <int ind, typename T, typename U>
sparse<T>& operator*=(target<0, ind, sparse<T> > x, const U& value) {
	return operator*=(*(x.data), value);
}
template <int ind, typename T, typename U>
sparse<T> operator*(const U& value, const target<0, ind, sparse<T> >& x) {
	return operator*(value, *(x.data));
}
template <typename T> bool operator==(const sparse<T>& a, const sparse<T>& b);


} // namespace MMSP

#include "MMSP.sparse.cpp"

#endif
