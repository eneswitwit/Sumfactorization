#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <limits>
#include <cassert>
#include <type_traits>
#include <iostream>

#include "../../include/constexpr_math.h"
#include "../../include/constexpr_array.h"
#include "../../include/constexpr_quadrature.h"
#include "../../include/polynomialbasis/constexpr_lagrange.h"


// Epsilon for asserting
constexpr long double eps = 0.00000001;

// Order for Interpolation
constexpr static size_t order = 5;


int main() {

    /** TEST CONSTEXPR **/
    constexpr Lagrange<long double, order, Quadrature> lagr;
    /** lagr object is known at compile time **/


    /** TEST EVAL LAGRANGE **/

    // If we evaluate at node points we expect 1
    for (unsigned int i = 0; i < order + 1; i++) {
        assert( lagr.eval_lagrange(lagr.knots_[i], i) <= 1 + eps);
        assert( lagr.eval_lagrange(lagr.knots_[i], i) >= 1 - eps);
    }

    // If we evalute the i-th lagrange polynomial at node j-th node point we except 0
    for (unsigned int j = 0; j < order + 1 ; j++) {
        for (unsigned int i = 0; i < order + 1; i++) {
            if ( i != j ) {
                assert( lagr.eval_lagrange(lagr.knots_[i], j) <= 0 + eps);
                assert( lagr.eval_lagrange(lagr.knots_[i], j) >= 0 - eps);
            }
        }
    }

    /** TEST EVAL GRADIENT LAGRANGE **/


    return 0;
}
