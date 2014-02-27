// Unit tests for point operators

#include <iostream>
#include "point.hpp"

int main() {
	Point<float> a(1.125, 1.125, 1.125);
	Point<float> b(2.125 * a);
	float m = 2.;

	std::cout << b << " + " << a << " = " << b + a << std::endl;
	std::cout << b << " - " << a << " = " << b - a << std::endl;
	std::cout << b << " * " << a << " = " << b*a << std::endl;
	std::cout << a << " / " << b << " = " << a / b << std::endl;
	std::cout << b << " / " << m << " = " << b / m << std::endl;
	std::cout << m << " * " << b << " = " << m*b << std::endl;
	std::cout << "Distance between " << b << " and " << a << " is " << distance_function(a, b) << "." << std::endl;
	return 0;
}
