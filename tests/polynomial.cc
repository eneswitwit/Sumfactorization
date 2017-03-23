#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <limits>
#include <type_traits>
#include <iostream>

#include "../include/constexpr_math.h"
#include "../include/constexpr_array.h"
#include "../include/constexpr_quadrature.h"
#include "../include/polynomialbasis/constexpr_lagrange.h"

constexpr unsigned int order = 3;
Lagrange<order, long double,Quadrature> lagr;
Quadrature<long double,order> quad;

long double *ptr;
long double knot = quad.knots_[1];

int main() {

    std::cout << lagr.eval_lagrange<&knot,0>();

    /*for (int i=0;i<order+1;i++) {
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
