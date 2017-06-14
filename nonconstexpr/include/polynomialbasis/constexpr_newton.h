#ifndef __CONSTEXPR_Newton_H__
#define __CONSTEXPR_Newton_H__

#include <cstdint>
#include "../constexpr_quadrature.h"

// A class for handling polynomials with newton basis.

/**
 * This class will be able to evaluate a newton basis function at a given point.
 * We will have a number of knots depending on our order with a basis function each corresponding to a certain knot.
 * The desired bilinear forms also require first and second derivatives of the basis function.
 *
 * The template parameter 'Number' defines the number type. 'order' determines how many knots/basis functions we will have.
 * We also pass our quadrature  class, we will use it to compute the points from the Gauss-Lobatto quadrature rule.
 */

template <typename Number, std::size_t order, template<typename, std::size_t> class Quadrature>
class Newton {
public:

    const constexpr_array <Number, order+1> knots_;

    constexpr  Newton() : knots_(compute_knots()) {}

    constexpr constexpr_array <Number, order + 1> compute_knots() const {
        constexpr Quadrature<Number,order> quad;
        return quad.compute_quadrature_points();
    }

    /**
     * The i-th Newton-polynomial is evaluated using the standard formula N_i(x)=\prod_{j=0}^{i-1} (x-x_j) \
     */
    constexpr Number eval(const Number & x, const unsigned int i) const {
        Number val = 1.;
        if (i==0)
            return val;

        for (unsigned int j = 0; j < i; j++) {
            val*=(x-knots_[j]);
        }
        return val;
    }

    /**
     * We compute the derivative using the product rule. N_i^'(x)= \sum_{j=0}^{i-1} \prod_{k\neq j} (x-x_k)
     */

    constexpr Number eval_1st_derivative(const Number & x, const unsigned int i) const {
        Number sum=0.;
        if (i==0)
            return sum;
        for (unsigned int j = 0; j < i; j++) {
            Number prod=1.;
            for (unsigned int k = 0; k < i; k++) {
                if (k!=j)
                    prod*=(x-knots_[k]);
            }
            sum+=prod;
        }
        return sum;
    }

    /**
     * We apply the product rule again to each summand. P_i^{''}(x)= \sum_{j=0}^{i-1} \sum_{k\neq j} \prod_{l\neq j,k} (x-x_l)
     */

    constexpr Number eval_2nd_derivative(const Number & x, const unsigned int i) const {
        Number sum=0;
        for (unsigned int j=0; j<i ;j++) {
            Number sum1=0;
            for (unsigned int k=0; k<i; k++) {
                if (k!=j) {
                    Number prod=1.;
                    for (unsigned int l=0; l<i; l++) {
                        if (l!=k && l!=j)
                            prod*=(x-knots_[l]);
                    }
                    sum1+=prod;
                }
            }
            sum+=sum1;
        }
        return sum;
    }



};

#endif
