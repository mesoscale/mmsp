// MMSP.sparse.hpp
// Class definition for the MMSP sparse data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_SPARSE
#define MMSP_SPARSE
#include <cmath>
#include <limits>
#include <cassert>
#include"MMSP.utility.hpp"

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
	sparse() {
		size = 0;
		data = NULL;
	}
	sparse(const sparse& x) {
		size = x.length();
		data = new item<T>[size];
		memcpy(data, x.data, size * sizeof(item<T>));
	}
	~sparse() {
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
	}

	// assignment operator
	sparse& operator=(const sparse& x) {
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		size = x.length();
		data = new item<T>[size];
		//memcpy(data, x.data, size * sizeof(item<T>));
		for (int i = 0; i < size; i++) {
			data[i].index = x.index(i);
			data[i].value = x.value(i);
		}
		return *this;
	}
	template <typename U> sparse& operator=(const sparse<U>& x) {
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		size = x.length();
		data = new item<T>[size];
		for (int i = 0; i < size; i++) {
			data[i].index = x.index(i);
			data[i].value = static_cast<T>(x.value(i));
		}
		return *this;
	}

	// data access operators
	T& set(int index) {
		for (int i = 0; i < size; i++)
			if (data[i].index == index)
				return data[i].value;

		item<T>* temp = new item<T>[size];
		memcpy(temp, data, size * sizeof(item<T>));
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		data = new item<T>[size + 1];
		memcpy(data, temp, size * sizeof(item<T>));
		if (temp!=NULL){
			delete [] temp;
			temp=NULL;
		}
		size += 1;
		data[size - 1].index = index;
		data[size - 1].value = static_cast<T>(0);
		return data[size - 1].value;
	}

	T operator[](int index) const {
		for (int i = 0; i < size; i++)
			if (data[i].index == index)
				return data[i].value;
		return static_cast<T>(0);
	}

	unsigned int grain_id() const {
		unsigned int max_index = 0;
		T max_value = -1.0;

		for (int i = 0; i < size; i++) {
			if (data[i].value > max_value) {
				max_index = data[i].index;
				max_value = data[i].value;
			}
		}
		assert(max_index < std::numeric_limits<unsigned int>::max());
		return max_index;
	}

	double getMagPhi() const {
		double sum = 0.0;

		for (int i = 0; i < size; i++) {
			double phi = data[i].value;
			sum += phi * phi;
		}

		return sqrt(sum);
	}

	int index(int i) const {
		return data[i].index;
	}
	T value(int i) const {
		return data[i].value;
	}

	// buffer I/O functions
	int buffer_size() const {
		return sizeof(size) + size * sizeof(item<T>);
	}
	int to_buffer(char* buffer) const {
		memcpy(buffer, &size, sizeof(size));
		memcpy(buffer + sizeof(size), data, size * sizeof(item<T>));
		return sizeof(size) + size * sizeof(item<T>);
	}
	int from_buffer(const char* buffer) {
		memcpy(&size, buffer, sizeof(size));
		if (data!=NULL)	{
			delete [] data;
			data=NULL;
		}
		data = new item<T>[size];
		memcpy(data, buffer + sizeof(size), size * sizeof(item<T>));
		return sizeof(size) + size * sizeof(item<T>);
	}

	// file I/O functions
	void write(std::ofstream& file) const {
		file.write(reinterpret_cast<const char*>(&size), sizeof(size));
		file.write(reinterpret_cast<const char*>(data), size * sizeof(item<T>));
	}
	void read(std::ifstream& file) {
		file.read(reinterpret_cast<char*>(&size), sizeof(size));
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		data = new item<T>[size];
		file.read(reinterpret_cast<char*>(data), size * sizeof(item<T>));
	}

	// utility functions
	int length() const {
		return size;
	}
	void resize(int n) {}
	void copy(const sparse& s) {
		size = s.size;
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		data = new item<T>[size];
		memcpy(data, s.data, size * sizeof(item<T>));
	}
	void swap(sparse& s) {
		item<T>* t = data;
		data = s.data;
		s.data = t;
		int* l = size;
		size = s.size;
		s.size = l;
	}

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
template <typename T> sparse<T> max(const sparse<T>& x, const sparse<T>& y) {
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++) {
		T val = x.value(i);
		int ind = x.index(i);
		z.set(ind) = MMSP::max(val, y[ind]);
	}
	int N2 = y.length();
	for (int i = 0; i < N2; i++) {
		T val = y.value(i);
		int ind = y.index(i);
		z.set(ind) = MMSP::max(val, x[ind]);
	}
	return z;
}

template <typename T> sparse<T> min(const sparse<T>& x, const sparse<T>& y) {
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++) {
		T val = x.value(i);
		int ind = x.index(i);
		z.set(ind) = MMSP::min(val, y[ind]);
	}
	int N2 = y.length();
	for (int i = 0; i < N2; i++) {
		T val = y.value(i);
		int ind = y.index(i);
		z.set(ind) = MMSP::min(val, x[ind]);
	}
	return z;
}

template <typename T, typename U> sparse<T>& operator+=(sparse<T>& x, const sparse<U>& y) {
	int N = y.length();
	for (int i = 0; i < N; i++)
		x.set(y.index(i)) += static_cast<T>(y.value(i));
	return x;
}
template <typename T, typename U> sparse<T> operator+(const sparse<T>& x, const sparse<U>& y) {
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++)
		z.set(x.index(i)) += static_cast<T>(x.value(i));
	int N2 = y.length();
	for (int i = 0; i < N2; i++)
		z.set(y.index(i)) += static_cast<T>(y.value(i));
	return z;
}
template <typename T, typename U> sparse<T>& operator-=(sparse<T>& x, const sparse<U>& y) {
	int N = y.length();
	for (int i = 0; i < N; i++)
		x.set(y.index(i)) -= static_cast<T>(y.value(i));
	return x;
}
template <typename T, typename U> sparse<T> operator-(const sparse<T>& x, const sparse<U>& y) {
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++)
		z.set(x.index(i)) += static_cast<T>(x.value(i));
	int N2 = y.length();
	for (int i = 0; i < N2; i++)
		z.set(y.index(i)) -= static_cast<T>(y.value(i));
	return z;
}
template <typename T, typename U> sparse<T>& operator*=(sparse<T>& x, const U& value) {
	int N = x.length();
	for (int i = 0; i < N; i++)
		x.set(x.index(i)) *= static_cast<T>(value);
	return x;
}
template <typename T, typename U> sparse<T> operator*(const U& value, const sparse<T>& x) {
	sparse<T> z;
	int N = x.length();
	for (int i = 0; i < N; i++)
		z.set(x.index(i)) = static_cast<T>(value) * x.value(i);
	return z;
}


// target class: dim = 0 specialization for sparse class
template <int ind, typename T>
class target<0, ind, sparse<T> > {
public:
	target(sparse<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}

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
template <typename T> bool operator==(const sparse<T>& a, const sparse<T>& b) {
  int N=a.length();
  if (N != b.length()) return false;
  for (int i=0; i<N; ++i) {
  	int indexA = a.index(i);
  	bool found=false;
  	bool match=false;
  	for (int j=0; j<N && !found; ++j) {
  		int indexB = b.index(j);
  		if (indexA==indexB) {
  			found=true;
  			match = (a.value(i) == b.value(j));
  		}
  	}
  	if (!found || (found && !match)) return false;
	}
  return true;
}

} // namespace MMSP

#endif
