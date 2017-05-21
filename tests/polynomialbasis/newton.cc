#include <cassert> 
#include "../../include/polynomialbasis/constexpr_newton.h"
#include <iostream>


// Epsilon for asserting
constexpr long double eps = 0.001;

// Order for Interpolation
constexpr static size_t order = 3;


int main() {

    /** TEST CONSTEXPR **/
    constexpr Newton<long double, order, Quadrature> newton;
    /** newton object is known at compile time **/


    /** TEST EVAL**/

    // 1st Newton Basis
    assert( abs(newton.eval(newton.knots_[0],0) - 1) <= eps );
    assert( abs(newton.eval(newton.knots_[1],0) - 1) <= eps );
    assert( abs(newton.eval(newton.knots_[2],0) - 1) <= eps );
    assert( abs(newton.eval(newton.knots_[3],0) - 1) <= eps );
    // 2nd Newton Basis
    assert( abs(newton.eval(newton.knots_[0],1) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[1],1) - 0.276393) <= eps );
    assert( abs(newton.eval(newton.knots_[2],1) - 0.723607) <= eps );
    assert( abs(newton.eval(newton.knots_[3],1) - 1) <= eps );
    // 3rd Newton Basis
    assert( abs(newton.eval(newton.knots_[0],2) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[1],2) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[2],2) - 0.323607) <= eps );
    assert( abs(newton.eval(newton.knots_[3],2) - 0.723607) <= eps );
    // 4th Newton Basis
    assert( abs(newton.eval(newton.knots_[0],3) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[1],3) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[2],3) - 0) <= eps );
    assert( abs(newton.eval(newton.knots_[3],3) - 0.2) <= eps );

    /** TEST EVAL GRADIENT Newton **/

    // 1st Newton Basis
    assert( abs(newton.eval_1st_derivative(newton.knots_[0],0) - 0) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[1],0) - 0) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[2],0) - 0) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[3],0) - 0) <= eps );
    // 2nd Newton Basis
    assert( abs(newton.eval_1st_derivative(newton.knots_[0],1) - 1) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[1],1) - 1) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[2],1) - 1) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[3],1) - 1) <= eps );
    // 3rd Newton Basis
    assert( abs(newton.eval_1st_derivative(newton.knots_[0],2) + 0.276393) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[1],2) - 0.276393) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[2],2) - 1.17082) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[3],2) - 1.72361) <= eps );
    // 4th Newton Basis
    assert( abs(newton.eval_1st_derivative(newton.knots_[0],3) - 0.2) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[1],3) + 0.123607) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[2],3) - 0.323606) <= eps );
    assert( abs(newton.eval_1st_derivative(newton.knots_[3],3) - 1.2) <= eps );

    /** TEST EVAL LAPLACIAN **/

    // 1st Newton Basis
    assert( abs(newton.eval_2nd_derivative(newton.knots_[0],0)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[1],0)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[2],0)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[3],0)) <= eps );
    // 2nd Newton Basis
    assert( abs(newton.eval_2nd_derivative(newton.knots_[0],1)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[1],1)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[2],1)) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[3],1)) <= eps );
    // 3rd Newton Basis
    assert( abs(newton.eval_2nd_derivative(newton.knots_[0],2) - 2) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[1],2) - 2) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[2],2) - 2) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[3],2) - 2) <= eps );
    // 4rd Newton Basis
    assert( abs(newton.eval_2nd_derivative(newton.knots_[0],3) + 2) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[1],3) + 0.341643) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[2],3) - 2.34164) <= eps );
    assert( abs(newton.eval_2nd_derivative(newton.knots_[3],3) - 4) <= eps );

    return 0;
}
