// Test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>

std::vector<double> chebyshev_knot(int order,int lo_bound, int up_bound);
double evaluate_polynomial(double x, double y, std::vector<std::vector<double>> coefficient_matrix, std::vector<double> knots_x, std::vector<double> knots_y);

int main()
{
	int const order = 1;

	std::vector<std::vector<double>> matrix;
	matrix.resize(2, std::vector<double>(2, 0));
	matrix[0][0] = 1;
	matrix[1][0] = 1;
	matrix[0][1] = 1;
	matrix[1][1] = 1;

	std::vector<double> knots(order+1);
	knots = chebyshev_knot(order,-1,1);
	for (int i = 0; i <= order; i++) {
		std::cout << knots[i] << "         ";
	}




	double solution;
	solution = evaluate_polynomial(0, 0, matrix, knots, knots);
	std::cout << solution;

    return 0;
}

std::vector<double> chebyshev_knot(int order,int lo_bound,int up_bound)
{
	const double pi = 3.141592653589793238463;
	std::vector<double> knots(order+1);
	for (int i = 0; i <= order; i++) {
		knots[i] = (lo_bound+up_bound)/2 + ((up_bound-lo_bound)/2)*cos(((2 * i) + 1)*pi / (2 * (order+1)));
	}
	return knots;
}

double evaluate_polynomial(double x, double y, std::vector<std::vector<double>> coefficient_matrix, std::vector<double> knots_x, std::vector<double> knots_y)
{
	int n = knots_x.size();
	std::vector<double> coefficients_b(n);
	for (int j = 0; j <= n-1; j++) {
		double aux = coefficient_matrix[n-1][j];
		for (int i = n-1; i >0 ; i--) {
			//std::cout << "x: " << x << "x_i: " << knots_x[i] << "coeff: " << coefficient_matrix[i][j];
			aux = aux*(x - knots_x[i-1]) + coefficient_matrix[i][j];
			std::cout << "aux "  << aux << "      ";
		}
		coefficients_b[j] = aux;
		std::cout << "b[j] " << aux << "     ";
	}
	double val = coefficients_b[n-1];
	for (int j = n - 1; j > 0; j--) {
		val = val*(y - knots_y[j-1]) + coefficients_b[j];
	}
	return val;
}