#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <limits>
#include <type_traits>
#include <iostream>

#include "../include/math_constexpr.h"
#include "../include/constexpr_array.h"
#include "../include/quadrature_constexpr.h"
#include "../include/polynomial.h"

int main() {
    constexpr unsigned int order = 3;
    Polynomial<order, long double,Quadrature> poly;
   /* Quadrature<order, long double> quad;
    quad.compute_quadrature_points(order+1,1,1);
    quad.compute_quadrature_weights(quad.knots,0,0);

    std::cout << quad.knots[0] << "    " << quad.knots[1] << "    " << quad.knots[2] << "    " << quad.knots[3] << "    " << std::endl << std::endl;

    for (int i=0;i<order+1;i++) {
            std::cout << poly.eval_lagrange(0,quad.knots[i],quad.knots) << "       " << poly.eval_lagrange(1,quad.knots[i],quad.knots) << "       "
                      << poly.eval_lagrange(2,quad.knots[i],quad.knots) << "       " << poly.eval_lagrange(3,quad.knots[i],quad.knots) << std::endl;
        }
    std::cout << std::endl;

    for (int i=0;i<order+1;i++) {
            std::cout << poly.eval_1st_derivative(0,quad.knots[i],quad.knots) << "       " << poly.eval_1st_derivative(1,quad.knots[i],quad.knots) << "       "
                      << poly.eval_1st_derivative(2,quad.knots[i],quad.knots) << "       " << poly.eval_1st_derivative(3,quad.knots[i],quad.knots) << std::endl;
        }
    std::cout << std::endl;

    for (int i=0;i<order+1;i++) {
            std::cout << poly.eval_2nd_derivative(0,quad.knots[i],quad.knots) << "       " << poly.eval_2nd_derivative(1,quad.knots[i],quad.knots) << "       "
                      << poly.eval_2nd_derivative(2,quad.knots[i],quad.knots) << "       " << poly.eval_2nd_derivative(3,quad.knots[i],quad.knots) << std::endl;
        }*/

    return 0;
}
