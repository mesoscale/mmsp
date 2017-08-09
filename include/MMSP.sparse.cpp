// MMSP.sparse.h
// Class implementation for the MMSP sparse data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#include <cmath>
#include <limits>
#include <cassert>

namespace MMSP
{

// constructors / destructor
template<typename T> sparse<T>::sparse()
{
	size = 0;
	data = NULL;
}

template<typename T> sparse<T>::sparse(const sparse& x)
{
	size = x.length();
	data = new item<T>[size];
	memcpy(data, x.data, size * sizeof(item<T>));
}

template<typename T> sparse<T>::~sparse()
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
}

// assignment operator
template<typename T> sparse<T>& sparse<T>::operator=(const sparse& x)
{
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

template<typename T> template <typename U> sparse<T>& sparse<T>::operator=(const sparse<U>& x)
{
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
template<typename T> T& sparse<T>::set(int index)
{
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
	if (temp!=NULL) {
		delete [] temp;
		temp=NULL;
	}
	size += 1;
	data[size - 1].index = index;
	data[size - 1].value = static_cast<T>(0);
	return data[size - 1].value;
}

template<typename T> T sparse<T>::operator[](int index) const
{
	for (int i = 0; i < size; i++)
		if (data[i].index == index)
			return data[i].value;
	return static_cast<T>(0);
}

template <typename T> unsigned int sparse<T>::grain_id() const
{
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

template <typename T> double sparse<T>::getMagPhi() const
{
	double sum = 0.0;

	for (int i = 0; i < size; i++) {
		double phi = data[i].value;
		sum += phi * phi;
	}

	return sqrt(sum);
}

template <typename T> int sparse<T>::index(const int& i) const
{
	assert(i < size);
	return data[i].index;
}

template <typename T> T sparse<T>::value(const int& i) const
{
	assert(i < size);
	return data[i].value;
}

// buffer I/O functions
template <typename T> int sparse<T>::buffer_size() const
{
	return sizeof(size) + size * sizeof(item<T>);
}

template <typename T> int sparse<T>::to_buffer(char* buffer) const
{
	memcpy(buffer, &size, sizeof(size));
	memcpy(buffer + sizeof(size), data, size * sizeof(item<T>));
	return sizeof(size) + size * sizeof(item<T>);
}

template <typename T> int sparse<T>::from_buffer(const char* buffer)
{
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
template <typename T> void sparse<T>::write(std::ofstream& file) const
{
	file.write(reinterpret_cast<const char*>(&size), sizeof(size));
	file.write(reinterpret_cast<const char*>(data), size * sizeof(item<T>));
}

template <typename T> void sparse<T>::read(std::ifstream& file)
{
	file.read(reinterpret_cast<char*>(&size), sizeof(size));
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	data = new item<T>[size];
	file.read(reinterpret_cast<char*>(data), size * sizeof(item<T>));
}

// utility functions
template <typename T> int sparse<T>::length() const
{
	return size;
}

template <typename T> void sparse<T>::resize(int n) {}

template <typename T> void sparse<T>::copy(const sparse& s)
{
	size = s.size;
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	data = new item<T>[size];
	memcpy(data, s.data, size * sizeof(item<T>));
}

template <typename T> void sparse<T>::swap(sparse& s)
{
	item<T>* t = data;
	data = s.data;
	s.data = t;
	int l = size;
	size = s.size;
	s.size = l;
}

// numerical operators
template <typename T> sparse<T> max(const sparse<T>& x, const sparse<T>& y)
{
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

template <typename T> sparse<T> min(const sparse<T>& x, const sparse<T>& y)
{
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

template <typename T, typename U> sparse<T>& operator+=(sparse<T>& x, const sparse<U>& y)
{
	int N = y.length();
	for (int i = 0; i < N; i++)
		x.set(y.index(i)) += static_cast<T>(y.value(i));
	return x;
}

template <typename T, typename U> sparse<T> operator+(const sparse<T>& x, const sparse<U>& y)
{
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++)
		z.set(x.index(i)) += static_cast<T>(x.value(i));
	int N2 = y.length();
	for (int i = 0; i < N2; i++)
		z.set(y.index(i)) += static_cast<T>(y.value(i));
	return z;
}

template <typename T, typename U> sparse<T>& operator-=(sparse<T>& x, const sparse<U>& y)
{
	int N = y.length();
	for (int i = 0; i < N; i++)
		x.set(y.index(i)) -= static_cast<T>(y.value(i));
	return x;
}

template <typename T, typename U> sparse<T> operator-(const sparse<T>& x, const sparse<U>& y)
{
	sparse<T> z;
	int N1 = x.length();
	for (int i = 0; i < N1; i++)
		z.set(x.index(i)) += static_cast<T>(x.value(i));
	int N2 = y.length();
	for (int i = 0; i < N2; i++)
		z.set(y.index(i)) -= static_cast<T>(y.value(i));
	return z;
}

template <typename T, typename U> sparse<T>& operator*=(sparse<T>& x, const U& value)
{
	int N = x.length();
	for (int i = 0; i < N; i++)
		x.set(x.index(i)) *= static_cast<T>(value);
	return x;
}

template <typename T, typename U> sparse<T> operator*(const U& value, const sparse<T>& x)
{
	sparse<T> z;
	int N = x.length();
	for (int i = 0; i < N; i++)
		z.set(x.index(i)) = static_cast<T>(value) * x.value(i);
	return z;
}

template <int ind, typename T> target<0,ind,sparse<T> >::target(sparse<T>* DATA, const int* S0, const int* SX, const int* X0, const int* X1, const int* B0, const int* B1)
{
	data = DATA;
	s0 = S0;
	sx = SX;
	x0 = X0;
	x1 = X1;
	b0 = B0;
	b1 = B1;
}

template <typename T> bool operator==(const sparse<T>& a, const sparse<T>& b)
{
	int N=a.length();
	if (N != b.length()) return false;
	for (int i = 0; i < N; i++) {
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
