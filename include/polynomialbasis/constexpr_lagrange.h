#ifndef __CONSTEXPR_LAGRANGE_H__
#define __CONSTEXPR_LAGRANGE_H__

#include <cstdint>
#include "../constexpr_quadrature.h"

// A class for handling polynomials with lagrangian basis.

template <typename Number, std::size_t order, template<typename, std::size_t> class Quadrature>
class Lagrange {
public:
 
   const constexpr_array <Number, order+1> knots_;

   constexpr  Lagrange() : knots_(compute_knots()) {}

    constexpr constexpr_array <Number, order + 1> compute_knots() const {
        constexpr Quadrature<Number,order> quad;
        return quad.compute_quadrature_points();
    }


    constexpr Number eval_lagrange(const Number & x, const unsigned int i) const {
        Number val = 1.;
        for (unsigned int j = 0; j <= order; j++) {
                if (i != j) {
                        val *= (x - knots_[j]) / (knots_[i] - knots_[j]);
                    }
            }
        return val;
    }

    constexpr Number eval_1st_derivative(const Number & x, const unsigned int i) const {
        Number sum=0;
        for (unsigned int m=0;m<order+1;m++) {
                if (m!=i) {
                        Number product=1;
                        for (unsigned int j=0;j<order+1;j++) {
                                if ((j!=i) && (j!=m)) {
                                        product*=(x-knots_[j])/(knots_[i]-knots_[j]);
                                    }
                            }
                        sum+=product/(knots_[i]-knots_[m]);
                    }
            }
        return sum;
    }

    constexpr Number eval_2nd_derivative(const Number & x, const unsigned int i) const {
        Number sum=0;
        for (unsigned int k=0;k<order+1;k++) {
                if (k!=i) {
                        Number sum1=0;
                        for (unsigned int m=0;m<order+1;m++) {
                                if ((m!=i)&&(m!=k)) {
                                        Number product=1;
                                        for (unsigned int j=0;j<order+1;j++) {
                                                if ((j!=i) && (j!=m) && (j!=k)) {
                                                        product*=(x-knots_[j])/(knots_[i]-knots_[j]);
                                                    }
                                            }
                                        sum1+=product/(knots_[i]-knots_[m]);
                                    }
                            }
                        sum+=sum1/(knots_[i]-knots_[k]);
                    }
            }
        return sum;
    }



};

#endif
