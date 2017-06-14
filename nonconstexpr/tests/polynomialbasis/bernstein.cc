#include <cassert> 
#include "../../include/polynomialbasis/constexpr_bernstein.h"
#include <iostream>


// Epsilon for asserting
constexpr long double eps = 0.001;

// Order for Interpolation
constexpr static size_t order = 3;


int main() {

    /** TEST CONSTEXPR **/
    constexpr Bernstein<long double, order, Quadrature> bernst;
    /** bernst object is known at compile time **/

    /** TEST EVAL**/

    // 1st bernst Basis
    assert( abs(bernst.eval(bernst.knots_[0],0) - 1) <= eps );
    assert( abs(bernst.eval(bernst.knots_[1],0) - 0.378886) <= eps );
    assert( abs(bernst.eval(bernst.knots_[2],0) - 0.0211145) <= eps );
    assert( abs(bernst.eval(bernst.knots_[3],0) - 0) <= eps );
    // 2nd bernst Basis
    assert( abs(bernst.eval(bernst.knots_[0],1) - 0) <= eps );
    assert( abs(bernst.eval(bernst.knots_[1],1) - 0.434164) <= eps );
    assert( abs(bernst.eval(bernst.knots_[2],1) - 0.165836) <= eps );
    assert( abs(bernst.eval(bernst.knots_[3],1) - 0) <= eps );
    // 3rd bernst Basis
    assert( abs(bernst.eval(bernst.knots_[0],2) - 0) <= eps );
    assert( abs(bernst.eval(bernst.knots_[1],2) - 0.165836) <= eps );
    assert( abs(bernst.eval(bernst.knots_[2],2) - 0.434164) <= eps );
    assert( abs(bernst.eval(bernst.knots_[3],2) - 0) <= eps );
    // 4th bernst Basis
    assert( abs(bernst.eval(bernst.knots_[0],3) - 0) <= eps );
    assert( abs(bernst.eval(bernst.knots_[1],3) - 0.0211145) <= eps );
    assert( abs(bernst.eval(bernst.knots_[2],3) - 0.378886) <= eps );
    assert( abs(bernst.eval(bernst.knots_[3],3) - 1) <= eps );

    /** TEST EVAL GRADIENT bernst **/

    // 1st bernst Basis
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[0],0) + 3) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[1],0) + 1.57082) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[2],0) + 0.229179) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[3],0) - 0) <= eps );
    // 2nd bernst Basis
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[0],1) - 3) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[1],1) - 0.370822) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[2],1) + 0.97082) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[3],1) - 0) <= eps );
    // 3rd bernst Basis
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[0],2) - 0) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[1],2) - 0.97082) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[2],2) + 0.370822) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[3],2) + 3) <= eps );
    // 4th bernst Basis
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[0],3) - 0) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[1],3) - 0.229179) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[2],3) - 1.57082) <= eps );
    assert( abs(bernst.eval_1st_derivative(bernst.knots_[3],3) - 3) <= eps );

    /** TEST EVAL LAPLACIAN **/

    // 1st bernst Basis
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[0],0) - 6) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[1],0) - 4.34164) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[2],0) - 1.65836) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[3],0)) <= eps );
    // 2nd bernst Basis
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[0],1) + 12) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[1],1) + 7.02493) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[2],1) - 1.02493) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[3],1) - 6) <= eps );
    // 3rd bernst Basis
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[0],2) - 6) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[1],2) - 1.02493) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[2],2) + 7.02493) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[3],2) + 12) <= eps );
    // 4rd bernst Basis
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[0],3)) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[1],3) - 1.65836) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[2],3) - 4.34164) <= eps );
    assert( abs(bernst.eval_2nd_derivative(bernst.knots_[3],3) - 6) <= eps );

    return 0;
}
