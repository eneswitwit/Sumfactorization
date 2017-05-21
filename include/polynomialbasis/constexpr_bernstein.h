#ifndef __CONSTEXPR_Bernstein_H__
#define __CONSTEXPR_Bernstein_H__

#include <cstdint>
#include "../constexpr_quadrature.h"
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <iostream>


// A class for handling polynomials with Bernstein basis.

/**
 * This class will be able to evaluate a Bernstein basis function at a given point.
 * We will have a number of knots depending on our order with a basis function each corresponding to a certain knot.
 * The desired bilinear forms also require first and second derivatives of the basis function.
 *
 * The template parameter 'Number' defines the number type. 'order' determines how many knots/basis functions we will have.
 * We also pass our quadrature  class, we will use it to compute the points from the Gauss-Lobatto quadrature rule.
 */

template <typename Number, std::size_t order, template<typename, std::size_t> class Quadrature>
class Bernstein {
public:

    const constexpr_array <Number, order+1> knots_;

    constexpr  Bernstein() : knots_(compute_knots()) {}

    constexpr constexpr_array <Number, order + 1> compute_knots() const {
        constexpr Quadrature<Number,order> quad;
        return quad.compute_quadrature_points();
    }

    /**
     * The i-th Bernstein-polynomial is evaluated using the standard formula N_i(x)=\prod_{j=0}^{i-1} (x-x_j) \
     */
    constexpr Number eval(const Number & x, const unsigned int i) const {
        unsigned int nCk = binom_coeff(order,i);
        return  nCk*std::pow(x,i)*std::pow(1-x,order-i);
    }

    /**
     * We compute the derivative using the product rule. N_i^'(x)= \sum_{j=0}^{i-1} \prod_{k\neq j} (x-x_k)
     */

    constexpr Number eval_1st_derivative(const Number & x, const unsigned int i) const {
        if (i==0)
            return (order)*std::pow(1-x,order-1)*(-1);
        if (i==order)
           return order*std::pow(x,order-1);

        Number nCk = binom_coeff(order,i);
        return  nCk*(i*std::pow(x,i-1)*std::pow(1-x,order-i)-(order-i)*std::pow(x,i)*std::pow(1-x,order-i-1));
    }

    /**
     * We apply the product rule again to each summand. P_i^{''}(x)= \sum_{j=0}^{i-1} \sum_{k\neq j} \prod_{l\neq j,k} (x-x_l)
     */

    constexpr Number eval_2nd_derivative(const Number & x, const unsigned int i) const {
        if (i==0)
            return order*(order-1)*std::pow(1-x,order-2);
        if (i==order)
            return order*(order-1)*std::pow(x,order-2);
        if (i==1)
            return order*(order-1)*std::pow(1-x,order-3)*(order*x-2);
        if (i==order-1)
            return order*(order-1)*std::pow(x,order-3)*(order*(1-x)-2);
        Number nCk = binom_coeff(order,i);
        return nCk*(+i*(i-1)*std::pow(x,i-2)*std::pow(1-x,order-i)
                    -2*i*(order-i)*std::pow(x,i-1)*std::pow(1-x,order-i-1)
                    +(order-i)*(order-i-1)*std::pow(x,i)*std::pow(1-x,order-i-2));
    }

    // Returns value of Binomial Coefficient C(n, k)
    constexpr unsigned int binom_coeff(int n, int k) const {
        unsigned int res = 1;
        unsigned int k1 = 0;

        // Since C(n, k) = C(n, n-k)
        if ( k > n - k )
            k1 = n - k;
        else k1=k;

        // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
        for (unsigned int i = 0; i < k1; ++i)
        {
            res *= (n - i);
            res /= (i + 1);
        }

        return res;
    }

};

#endif
