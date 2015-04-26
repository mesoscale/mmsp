#ifndef _Point_h_
#define _Point_h_
#include <iostream>
#include <cmath>

template <class T> class Point {
public:
	// Representation
	T x, y, z;
	// Constructors
	Point<T>(T a = 0, T b = 0, T c = 0): x(a), y(b), z(c) {}
	Point<T>(const Point<T>& old): x(old.x), y(old.y), z(old.z) {}
	template <class U>
	Point<T>(const Point<U>& old): x(old.x), y(old.y), z(old.z) {}
	//Accessors
	T& operator [](int i) {
		if (i==1) return this->y;
		else if (i==2) return this->z;
		return this->x;
	}
	const T& operator [](int i) const {
		if (i==1) return this->y;
		else if (i==2) return this->z;
		return this->x;
	}
	// Modifiers
	Point& operator+=(const Point<T>& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return (*this);
	}
	Point& operator-=(const Point<T>& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return (*this);
	}
	template <class U>
	Point<T>& operator*=(U rhs) {
		T n(rhs);
		this->x = x * n;
		this->y = y * n;
		this->z = z * n;
		return (*this);
	}
	template <class U>
	Point<T>& operator/=(U rhs) {
		T n(rhs);
		this->x = x / n;
		this->y = y / n;
		this->z = z / n;
		return (*this);
	}
};

// Arithmetic operators: Not member functions, but helpful utilities.
template <class T>
Point<T> operator+(const Point<T>& p, const Point<T>& q) {
	return Point<T>(p.x + q.x, p.y + q.y, p.z + q.z);
}

template <class T>
Point<T> operator-(const Point<T>& p, const Point<T>& q) {
	return Point<T>(p.x - q.x, p.y - q.y, p.z - q.z);
}

template <class T>
Point<T> operator*(const Point<T>& p, const Point<T>& q) {
	return Point<T>(p.x * q.x, p.y * q.y, p.z * q.z);
}

template <class U, class T>
Point<U> operator*(const U& f, const Point<T>& p) {
	return Point<U>(f * p.x, f * p.y, f * p.z);
}

template <class T>
Point<T> operator/(const Point<T>& p, const Point<T>& q) {
	return Point<T>(p.x / q.x, p.y / q.y, p.z / q.z);
}

template <class T, class U>
Point<T> operator/(const Point<T>& p, const U& f) {
	return Point<T>(p.x / f, p.y / f, p.z / f);
}


// Point output operator
template <class T> std::ostream& operator<<(std::ostream& out, const Point<T>& p) {
	out << '(' << p.x << ',' << p.y << ',' << p.z << ')';
	return out;
}

// Distance between points
template <class T>
double distance_function( Point<T> p, Point<T> q) {
	return sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
}

template <class T> bool operator==(const Point<T>& a, const Point<T>&b) {
	return (a.x==b.x) && (a.y==b.y) && (a.z==b.z);
}

#endif
