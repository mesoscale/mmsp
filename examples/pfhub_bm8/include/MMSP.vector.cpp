// MMSP.vector.cpp
// Class implementation for the MMSP vector data structure
// Questions/comments to gruberja@gmail.com (Jason Gruber)

namespace MMSP
{

//vector constructors and destructor
template<typename T> vector<T>::vector()
{
	size = 0;
	data = NULL;
}

template<typename T> vector<T>::vector(const vector& v)
{
	size = v.length();
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(v[i]);
}

template <typename T> template<typename U> vector<T>::vector(const vector<U>& v)
{
	size = v.length();
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(v[i]);
}

template<typename T> vector<T>::vector(int N)
{
	size = N;
	data = new T[size];
}

template <typename T> template<typename U> vector<T>::vector(int N, const U& value)
{
	size = N;
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(value);
}

template <typename T> vector<T>::~vector()
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
}

// assignment operators
template<typename T> vector<T>& vector<T>::operator=(const T& value)
{
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(value);
	return *this;
}

template<typename T> vector<T>& vector<T>::operator=(const vector& v)
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	size = v.length();
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(v[i]);
	return *this;
}

template <typename T> template <typename U> vector<T>& vector<T>::operator=(const U& value)
{
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(value);
	return *this;
}

template <typename T> template <typename U> vector<T>& vector<T>::operator=(const vector<U>& v)
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	size = v.length();
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = static_cast<T>(v[i]);
	return *this;
}

// buffer I/O functions
template<typename T> int vector<T>::buffer_size() const
{
	return sizeof(size) + size * sizeof(T);
}

template<typename T> int vector<T>::to_buffer(char* buffer) const
{
	memcpy(buffer, &size, sizeof(size));
	memcpy(buffer + sizeof(size), data, size * sizeof(T));
	return sizeof(size) + size * sizeof(T);
}

template<typename T> int vector<T>::from_buffer(const char* buffer)
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	memcpy(&size, buffer, sizeof(size));
	data = new T[size];
	memcpy(data, buffer + sizeof(size), size * sizeof(T));
	return sizeof(size) + size * sizeof(T);
}

// file I/O functions
template<typename T> void vector<T>::write(std::ofstream& file) const
{
	file.write(reinterpret_cast<const char*>(data), size * sizeof(T));
}

template<typename T> void vector<T>::read(std::ifstream& file)
{
	file.read(reinterpret_cast<char*>(data), size * sizeof(T));
}

// utility functions
template<typename T> int vector<T>::length() const
{
	return size;
}

template<typename T> void vector<T>::resize(int N)
{
	if (size == 0) {
		size = N;
		data = new T[size];
	} else if (N > size) {
		T* temp = new T[N];
		memcpy(temp, data, size * sizeof(T));
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		size = N;
		data = new T[size];
		memcpy(data, temp, size * sizeof(T));
		if (temp!=NULL) {
			delete [] temp;
			temp=NULL;
		}
	} else if (N < size) {
		T* temp = new T[N];
		memcpy(temp, data, N * sizeof(T));
		if (data!=NULL) {
			delete [] data;
			data=NULL;
		}
		size = N;
		data = new T[size];
		memcpy(data, temp, N * sizeof(T));
		if (temp!=NULL) {
			delete [] temp;
			temp=NULL;
		}
	}
}

template<typename T> void vector<T>::copy(const vector& v)
{
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	size = v.size;
	data = new T[size];
	memcpy(data, v.data, size * sizeof(T));
}

template<typename T> void vector<T>::swap(vector& v)
{
	T* temp = data;
	data = v.data;
	v.data = temp;
	int s = size;
	size = v.size;
	v.size = s;
}

template <typename T> template <typename U> void vector<T>::append(const U& value)
{
	T* temp = new T[size];
	memcpy(temp, data, size * sizeof(T));
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	size += 1;
	data = new T[size];
	memcpy(data, temp, size * sizeof(T));
	if (temp!=NULL) {
		delete [] temp;
		temp=NULL;
	}
	data[size - 1] = static_cast<T>(value);
}

template <typename T> template <typename U> void vector<T>::append(const vector<U>& v)
{
	int N = size;
	T* temp = new T[size];
	memcpy(temp, data, size * sizeof(T));
	if (data!=NULL) {
		delete [] data;
		data=NULL;
	}
	size += v.length();
	data = new T[size];
	memcpy(data, temp, size * sizeof(T));
	if (temp!=NULL) {
		delete [] temp;
		temp=NULL;
	}
	for (int i = N; i < size; i++)
		data[i] = static_cast<T>(v[i - N]);
}

template <typename T> bool operator==(const vector<T>& a, const vector<T>& b)
{
	int N=a.length();
	if (N != b.length()) return false;
	for (int i=0; i<N; ++i)
		if (a[i]!=b[i]) return false;
	return true;
}


} // namespace MMSP
