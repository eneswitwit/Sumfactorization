#ifndef __QUADRATURE_CONSTEXPR_H__
#define __QUADRATURE_CONSTEXPR_H__



// Template class Quadrature for computing knots and weights of the quadrature.
template<typename y_type, size_t order >
class Quadrature {
private:
    const constexpr_array < y_type, order + 1 > knots_;
public:
    constexpr Quadrature():
        knots_ ( compute_quadrature_points()) {}

    constexpr const y_type & operator[](size_t n) const {
        return knots_[n];
    }


    // Template constexpr function for calculating the input as constexpr for the constexpr sinus
    constexpr y_type sin_input(int i, int m) const {
        long double PI = 3.14159265358979323846;
        long double PI_HALF = 1.57079632679489661923;
        return (((2 * i + 1) / (2 * m) * PI) + PI_HALF);
    }


    // Template constexpr function for computing the Jacobi Polynomial
    template<size_t n>
    constexpr y_type JacobiP(y_type x, int alpha, int beta) const
    {
        constexpr_array < y_type, n+3 > p;
        p[0] = 1.0L;
        if (order == 0) return p[0];
        p[1] = ((alpha + beta + 2) * x + (alpha - beta)) / 2;
        if (order == 1) return p[1];

        for (unsigned int i = 1; i <= (order - 1); ++i)
        {
            const int v  = 2 * i + alpha + beta;
            const int a1 = 2 * (i + 1) * (i + alpha + beta + 1) * v;
            const int a2 = (v + 1) * (alpha * alpha - beta * beta);
            const int a3 = v * (v + 1) * (v + 2);
            const int a4 = 2 * (i + alpha) * (i + beta) * (v + 2);

            p[i + 1] = static_cast<y_type>( (a2 + a3 * x) * p[i] - a4 * p[i - 1]) / a1;
        }
        return p[n];
    }


    constexpr constexpr_array < y_type, order + 1 > compute_quadrature_points() const
    {
        constexpr_array < y_type, order + 1 > knots;
        constexpr unsigned int m = order - 1;
        for (unsigned int i = 0; i < m; ++i)
        {
            knots[i + 1] = -math::sin( sin_input(i, m) ) ;
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

    constexpr y_type compute_kth_entry(const unsigned int k, const constexpr_array < y_type, order + 1 > & knots) const
    {
        constexpr long double
        long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
        double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
        constexpr y_type runtime_one = 1.0;
        constexpr y_type tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );
        constexpr unsigned int m = order - 1;
        y_type s{1.};
        y_type J_x{1.};
        y_type f{1.};
        y_type delta{1.};
        y_type r = knots[k];
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
        while (math::fabs<y_type>(delta) >= tolerance);
        return r;
    }

    constexpr y_type transform_kth_entry(const unsigned int k, const constexpr_array < y_type, order + 1 > & knots) const
    {
        y_type val = knots[k];
        val *= 0.5L;
        val += 0.5L;
        return val;
    }

};



#endif
