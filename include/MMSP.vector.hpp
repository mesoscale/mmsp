// MMSP.vector.hpp
// Class definition for the MMSP vector data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef MMSP_VECTOR
#define MMSP_VECTOR
#include"MMSP.utility.hpp"

namespace MMSP{

template <typename T> class vector{
public:
	// constructors / destructor
	vector() {data=NULL; size=0;}
	~vector() {delete [] data;}

	// data access operators
	T& operator[](int i) {return data[i];}
	const T& operator[](int i) const {return data[i];}

	// assignment operator
	vector& operator=(const vector& v) {copy(v); return *this;}

	// buffer I/O functions
	int buffer_size() const
		{return size*sizeof(T);}
	int to_buffer(char* buffer) const
		{memcpy(buffer,data,size*sizeof(T)); return size*sizeof(T);}
	int from_buffer(const char* buffer)
		{memcpy(data,buffer,size*sizeof(T)); return size*sizeof(T);}

	// file I/O functions
	void write(std::ofstream& file) const
		{file.write(reinterpret_cast<const char*>(data),size*sizeof(T));}
	void read(std::ifstream& file)
		{file.read(reinterpret_cast<char*>(data),size*sizeof(T));}

	// utility functions
	int length() const
		{return size;}
	void resize(int n) {
		if (n>size) {
			T* temp = new T[n];
			memcpy(temp,data,size*sizeof(T));
			delete [] data;
			size = n;
			data = new T[size];
			memcpy(data,temp,n*sizeof(T));
			delete [] temp;
		}
		else if (n<size) {
			T* temp = new T[n];
			memcpy(temp,data,n*sizeof(T));
			delete [] data;
			size = n;
			data = new T[size];
			memcpy(data,temp,size*sizeof(T));
			delete [] temp;
		}
	}
	void copy(const vector& v) {
		size = v.size;
		delete [] data;
		data = new T[size];
		memcpy(data,v.data,size*sizeof(T));
	}
	void swap(vector& v) {
		T* t = data;
		data = v.data;
		v.data = t;
		int* s = size;
		size = v.size;
		v.size = s;
	}

private:
	// object data
	T* data;
	int size;
};

// buffer I/O functions
template <typename T> int buffer_size(const vector<T>& v) {return v.buffer_size();}
template <typename T> int to_buffer(const vector<T>& v, char* buffer) {return v.to_buffer(buffer);}
template <typename T> int from_buffer(vector<T>& v, const char* buffer) {return v.from_buffer(buffer);}

// file I/O functions
template <typename T> void write(const vector<T>& v, std::ofstream& file) {return v.write(file);}
template <typename T> void read(vector<T>& v, std::ifstream& file) {return v.read(file);}

// utility functions
template <typename T> int length(const vector<T>& v) {return v.length();}
template <typename T> void resize(vector<T>& v, int n) {v.resize(n);}
template <typename T> void copy(vector<T>& v, const vector<T>& w) {v.copy(w);}
template <typename T> void swap(vector<T>& v, vector<T>& w) {v.swap(w);}
template <typename T> std::string name(const vector<T>& s) {return std::string("vector:")+name(T());}


// target class: dim = 0 specialization for vector class
template <int ind, typename T>
class target<0,ind,vector<T> >{
public:
	target(vector<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1)
		{data=DATA; s0=S0; sx=SX; x0=X0; x1=X1; b0=B0; b1=B1;}

	operator vector<T>&() {return *data;}
	operator const vector<T>&() const {return *data;}

	T& operator[](int i) {return data->operator[](i);}
	const T& operator[](int i) const {return data->operator[](i);}

	int buffer_size() const {return data->buffer_size();}
	int to_buffer(char* buffer) const {return data->to_buffer(buffer);}
	int from_buffer(const char* buffer) const {return data->from_buffer(buffer);}

	void write(std::ofstream& file) const {data->write(file);}
	void read(std::ifstream& file) const {data->read(file);}

	int length() const {return data->length();}
	int resize(int n) const {return data->resize(n);}
	void copy(const target& t) const {return data->copy(t->data);}
	void swap(const target& t) const {return data->swap(t->data);}

	vector<T>* data;
	const int* s0;
	const int* sx;
	const int* x0;
	const int* x1;
	const int* b0;
	const int* b1;
};

// buffer I/O functions
template <int ind, typename T> int buffer_size(const target<0,ind,vector<T> >& v) {return v.buffer_size();}
template <int ind, typename T> int to_buffer(const target<0,ind,vector<T> >& v, char* buffer) {return v.to_buffer(buffer);}
template <int ind, typename T> int from_buffer(const target<0,ind,vector<T> >& v, const char* buffer) {return v.from_buffer(buffer);}

// file I/O functions
template <int ind, typename T> void write(const target<0,ind,vector<T> >& v, std::ofstream& file) {return v.write(file);}
template <int ind, typename T> void read(const target<0,ind,vector<T> >& v, std::ifstream& file) {return v.read(file);}

// utility functions
template <int ind, typename T> int length(const target<0,ind,vector<T> >& v) {return v.length();}
template <int ind, typename T> void resize(const target<0,ind,vector<T> >& v, int n) {v.resize(n);}
template <int ind, typename T> void copy(const target<0,ind,vector<T> >& v, const target<0,ind,vector<T> >& w) {v.copy(w);}
template <int ind, typename T> void swap(const target<0,ind,vector<T> >& v, const target<0,ind,vector<T> >& w) {v.swap(w);}
template <int ind, typename T> std::string name(const target<0,ind,vector<T> >& s) {return std::string("vector:")+name(T());}

} // namespace MMSP

#endif
