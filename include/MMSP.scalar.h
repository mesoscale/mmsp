// MMSP.scalar.h
// Class definition for the MMSP scalar data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_SCALAR
#define MMSP_SCALAR
#include"MMSP.utility.h"

namespace MMSP {

template <typename T> class scalar {
public:
	// constructors
	scalar() {}
	scalar(const T& value) {
		data = value;
	}
	scalar(const scalar& s) {
		data = s.data;
	}
	template <typename U> scalar(const U& value) {
		data = static_cast<T>(value);
	}
	template <typename U> scalar(const scalar<U>& s) {
		data = static_cast<T>(s);
	}

	// data access operators
	operator T&() {
		return data;
	}
	operator const T&() const {
		return data;
	}

	// assignment operators
	scalar& operator=(const T& value) {
		data = value;
		return *this;
	}
	scalar& operator=(const scalar& s) {
		data = s.data;
		return *this;
	}
	template <typename U> scalar& operator=(const U& value) {
		data = static_cast<T>(value);
		return *this;
	}
	template <typename U> scalar& operator=(const scalar<U>& s) {
		data = static_cast<T>(s);
		return *this;
	}

	// buffer I/O functions
	int buffer_size() const {
		return sizeof(T);
	}
	int to_buffer(char* buffer) const {
		memcpy(buffer, &data, sizeof(T));
		return sizeof(T);
	}
	int from_buffer(const char* buffer) {
		memcpy(&data, buffer, sizeof(T));
		return sizeof(T);
	}

	// file I/O functions
	void write(std::ofstream& file) const {
		file.write(reinterpret_cast<const char*>(&data), sizeof(T));
	}
	void read(std::ifstream& file) {
		file.read(reinterpret_cast<char*>(&data), sizeof(T));
	}

	// utility functions
	int length() const {
		return 1;
	}
	void resize(int n) {}
	void copy(const scalar& s) {
		memcpy(data, s.data, sizeof(T));
	}
	void swap(scalar& s) {
		T temp = data;
		data = s.data;
		s.data = temp;
	}

private:
	// object data
	T data;
};

// buffer I/O functions
template <typename T> int buffer_size(const scalar<T>& s) {
	return s.buffer_size();
}
template <typename T> int to_buffer(const scalar<T>& s, char* buffer) {
	return s.to_buffer(buffer);
}
template <typename T> int from_buffer(scalar<T>& s, const char* buffer) {
	return s.from_buffer(buffer);
}

// file I/O functions
template <typename T> void write(const scalar<T>& s, std::ofstream& file) {
	return s.write(file);
}
template <typename T> void read(scalar<T>& s, std::ifstream& file) {
	return s.read(file);
}

// utility functions
template <typename T> int length(const scalar<T>& s) {
	return s.length();
}
template <typename T> void resize(scalar<T>& s, int n) {
	s.resize(n);
}
template <typename T> void copy(scalar<T>& s, const scalar<T>& t) {
	s.copy(t);
}
template <typename T> void swap(scalar<T>& s, scalar<T>& t) {
	s.swap(t);
}
template <typename T> std::string name(const scalar<T>& s) {
	return std::string("scalar:") + name(T());
}


// target class: dim = 0 specialization for scalar class
template <int ind, typename T>
class target<0, ind, scalar<T> > {
public:
	// constructor
	target(scalar<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1) {
		data = DATA;
		s0 = S0;
		sx = SX;
		x0 = X0;
		x1 = X1;
		b0 = B0;
		b1 = B1;
	}

	// data access operators
	operator T&() {
		return *data;
	}
	operator const T&() const {
		return *data;
	}

	// assignment operators
	scalar<T>& operator=(const T& value) const {
		return data->operator=(value);
	}
	scalar<T>& operator=(const scalar<T>& s) const {
		return data->operator=(s);
	}
	template <typename U> scalar<T>& operator=(const U& value) const {
		return data->operator=(value);
	}
	template <typename U> scalar<T>& operator=(const scalar<U>& s) const {
		return data->operator=(s);
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

	// object data
	scalar<T>* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

// buffer I/O functions
template <int ind, typename T> int buffer_size(const target<0, ind, scalar<T> >& s) {
	return s.buffer_size();
}
template <int ind, typename T> int to_buffer(const target<0, ind, scalar<T> >& s, char* buffer) {
	return s.to_buffer(buffer);
}
template <int ind, typename T> int from_buffer(const target<0, ind, scalar<T> >& s, const char* buffer) {
	return s.from_buffer(buffer);
}

// file I/O functions
template <int ind, typename T> void write(const target<0, ind, scalar<T> >& s, std::ofstream& file) {
	return s.write(file);
}
template <int ind, typename T> void read(const target<0, ind, scalar<T> >& s, std::ifstream& file) {
	return s.read(file);
}

// utility functions
template <int ind, typename T> int length(const target<0, ind, scalar<T> >& s) {
	return s.length();
}
template <int ind, typename T> void resize(const target<0, ind, scalar<T> >& s, int n) {
	s.resize(n);
}
template <int ind, typename T> void copy(const target<0, ind, scalar<T> >& s, const target<0, ind, scalar<T> >& t) {
	s.copy(t);
}
template <int ind, typename T> void swap(const target<0, ind, scalar<T> >& s, const target<0, ind, scalar<T> >& t) {
	s.swap(t);
}
template <int ind, typename T> std::string name(const target<0, ind, scalar<T> >& s) {
	return std::string("scalar:") + name(T());
}

} // namespace MMSP

#endif
