#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>

#define PI 3.141592653589793238463

class Polynomial
{
    // Coefficient Matrix for Newton Basis
	std::vector<std::vector<double>> coefficient_matrix;

	// Sample points
	std::vector<double> knots_x;
	std::vector<double> knots_y;

	// Degree
	int degree;


  public:

	// Constructor
    Polynomial()
    : coefficient_matrix{}
	, knots_x{}
	, knots_y{}
    {}

    /*Polynomial(int n)
    : a_(std::max(1, n+1), int())
    {}

    Polynomial(std::initializer_list<int> coeffs)
    : a_{coeffs}
    {}*/

    
	// Operations
	// Horner
	double operator()(double x, double y) const
	{
	std::vector<double> coefficients_b(degree);
	for (int j = 0; j <= degree-1; j++) {
		double aux = coefficient_matrix[degree-1][j];
		for (int i = degree-1; i >0 ; i--) {
			aux = aux*(x - knots_x[i-1]) + coefficient_matrix[i-1][j];
		}
		coefficients_b[j] = aux;
	}
	double val = coefficients_b[degree-1];
		for (int j = degree - 1; j > 0; j--) {
			val = val*(y - knots_y[j-1]) + coefficients_b[j-1];
		}
	return val;
	}

	// Interpolation
	void chebyshev_knot(int lo_bound_x,int up_bound_x, int lo_bound_y, int up_bound_y)
	{
		knots_x.resize(degree+1);
		for (int i = 0; i <= degree; i++) {
			knots_x[i] = (lo_bound_x+up_bound_x)/2 + ((up_bound_x-lo_bound_x)/2)*cos(((2 * i) + 1)*PI / (2 * (degree+1)));
		}

		knots_y.resize(degree+1);
		for (int i = 0; i <= degree; i++) {
			knots_y[i] = (lo_bound_y+up_bound_y)/2 + ((up_bound_y-lo_bound_y)/2)*cos(((2 * i) + 1)*PI / (2 * (degree+1)));
		}
	}
	
	// Divided differences
	void divided_differences()
	{
		
	}
};

int main()
{
	Polynomial p;
	p.chebyshev_knot(0,1,0,1);
	p(0,0);


    std::cout << "Alle Tests erfolgreich!\n";
}


