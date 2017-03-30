#ifndef __CONSTEXPR_QUADRATURE_H__
#define __CONSTEXPR_QUADRATURE_H__

#include <limits>
#include <cmath>
#include <array>

#include "constexpr_array.h"
#include "constexpr_math.h"

// Template class Quadrature for computing knots and weights of the quadrature.
template<typename Number, std::size_t order >
class Quadrature {
private:

public:
    const constexpr_array < Number, order + 1 > knots_;
    const constexpr_array < Number, order + 1 > weights_;
    constexpr Quadrature():
        knots_ ( compute_quadrature_points()),
        weights_ ( compute_quadrature_weights()){}

    constexpr const Number & operator[](std::size_t n) const {
        return knots_[n];
    }


    // Template constexpr function for calculating the input as constexpr for the constexpr sinus
    constexpr Number sin_input(int i, int m) const {
        long double PI = 3.14159265358979323846;
        long double PI_HALF = 1.57079632679489661923;
        return (((2 * i + 1) / (2 * m) * PI) + PI_HALF);
    }


    // Template constexpr function for computing the Jacobi Polynomial
    template<std::size_t n>
    constexpr Number JacobiP(Number x, int alpha, int beta) const
    {
        constexpr_array < Number, n+1 > p;
        p[0] = 1.0L;
        if (order == 0) return p[0];
        p[1] = ((alpha + beta + 2) * x + (alpha - beta)) / 2;
        if (order == 1) return p[1];

        for (unsigned int i = 1; i <= (n - 1); ++i)
        {
            const int v  = 2 * i + alpha + beta;
            const int a1 = 2 * (i + 1) * (i + alpha + beta + 1) * v;
            const int a2 = (v + 1) * (alpha * alpha - beta * beta);
            const int a3 = v * (v + 1) * (v + 2);
            const int a4 = 2 * (i + alpha) * (i + beta) * (v + 2);

            p[i + 1] = static_cast<Number>( (a2 + a3 * x) * p[i] - a4 * p[i - 1]) / a1;
        }
        return p[n];
    }

    /**
     * The following functions are taken directly from the Deal.II library. We had to modify them slightly as we want them to work as constexpr functions.
     */


    constexpr constexpr_array < Number, order + 1 > compute_quadrature_points() const
    {
        constexpr_array < Number, order + 1 > knots;
        constexpr unsigned int m = order - 1;
        for (unsigned int i = 0; i < m; ++i)
        {
            knots[i + 1] = - sin( sin_input(i, m) ) ;
        }
        for (unsigned int k = 0; k < m; ++k)
        {
            knots[k+1] = compute_kth_entry(k+1, knots);
        }

        knots[0] = -1.L;
        knots[order] = +1.L;

        for (unsigned int j = 0; j < order + 1; j++) {
            knots[j] = transform_kth_entry(j, knots);
        }

        return knots;

    }

    constexpr constexpr_array<Number, order + 1> compute_quadrature_weights() const {

        constexpr_array < Number, order + 1 > weights;

        for (unsigned int j = 0; j < order+1; j++) {
          weights[j] = inv_transform_kth_entry(j,knots_);
        }

        Number s = 0.L;

        const Number factor = 2 * gamma<order+1>() * gamma<order+1>() / ((order)*gamma<order+1>()*gamma<order+2>());
        for (unsigned int i=0; i<order+1; ++i)
          {
            s = JacobiP<order>(weights[i], 0, 0);
            weights[i] = 0.5*factor/(s*s);
          }

        return weights;
    }

    template<std::size_t n>
    constexpr Number gamma() const
    {
      constexpr_array < Number, n-2 > result;
      result[0] = n - 1;
      for (unsigned int i=2; i<n-1; i++)
        result[i-1] = result[i-2]*i;
      return result[n-3];
    }

    constexpr Number compute_kth_entry(const unsigned int k, const constexpr_array < Number, order + 1 > & knots) const
    {
        constexpr long double
        long_double_eps = static_cast<Number>(std::numeric_limits<long double>::epsilon()),
        double_eps      = static_cast<Number>(std::numeric_limits<double>::epsilon());
        constexpr Number runtime_one = 1.0;
        constexpr Number tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );
        constexpr unsigned int m = order - 1;
        Number s{1.};
        Number J_x{1.};
        Number f{1.};
        Number delta{1.};
        Number r = knots[k];
        if (k > 1) {
            r = (r + knots[k - 1]) / 2;
        }

        do
        {
            s = 0.;
            for (unsigned int i = 1; i < k; ++i) {
                s += 1. / (r - knots[i]);
            }

            J_x   =  0.5 * (m + 3) * JacobiP<m-1>(r, 2, 2);
            f     = JacobiP<m>(r, 1, 1);
            delta = f / (f * s - J_x);
            r += delta;
        }
        while (fabs<Number>(delta) >= tolerance);
        return r;
    }

    /**
     * The following functions will transform the entries of a quadrature point. The points are computed for the interval [-1,1] but we want to use the interval [0,1].
     * The inverse transformation is needed for the computation of the quadrature weights.
     */

    constexpr Number transform_kth_entry(const unsigned int k, const constexpr_array < Number, order + 1 > & knots) const
    {
        Number val = knots[k];
        val *= 0.5L;
        val += 0.5L;
        return val;
    }

    constexpr Number inv_transform_kth_entry(const unsigned int k, const constexpr_array < Number, order + 1 > & knots) const
    {
        Number val = knots[k];
        val -= 0.5L;
        val *= 2.L;
        return val;
    }

};



#endif
