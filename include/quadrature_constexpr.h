#ifndef __QUADRATURE_CONSTEXPR_H__
#define __QUADRATURE_CONSTEXPR_H__



constexpr int alpha = 1;
constexpr int beta = 1;
constexpr long double PI = 3.14159265358979323846;
constexpr long double PI_HALF = 1.57079632679489661923;


// Template function for calculating the input as constexpr for sinus for the compute quadrature points function
template<typename y_type>
constexpr y_type sin_input(int i, int m) {
    return (((2 * i + 1) / (2 * m) * PI) - PI_HALF);
}

// Template function for calculating fabs as constexpr
template<typename y_type>
constexpr y_type fabs(y_type val) {
    if (val < 0) {
        return -val;
    }
    else {
        return val;
    }
}


// Template class Quadrature for computing knots and weights of the quadrature.
template<typename y_type, constexpr int order >
class Quadrature {
private:
    const constexpr_array < y_type, order + 1 > knots_;
public:
    constexpr Quadrature(const constexpr_array < y_type, order > & quad_sin_vec):
        knots_ ( compute_quadrature_points(quad_sin_vec)) {}


    constexpr constexpr_array < y_type, order + 1 > compute_quadrature_points(const constexpr_array < y_type, order > & quad_sin_vec) const
    {
        constexpr_array < y_type, order + 1 > knots;
        constexpr unsigned int m = order - 1;
        for (unsigned int i = 0; i < m; ++i)
        {
            knots[i + 1] = quad_sin_vec[i] ;
        }

        for (unsigned int k = 1; k < m; ++k)
        {
            knots[k] = compute_kth_entry(k, knots);
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
        y_type s, J_x, f, delta;
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

            J_x   =  0.5 * (m + 3) * JacobiP(r, 2, 2, m - 1);
            f     = JacobiP(r, alpha, beta, m);
            delta = f / (f * s - J_x);
            r += delta;
        }
        while (fabs<y_type>(delta) >= tolerance);
        return r;
    }

    constexpr y_type transform_kth_entry(const unsigned int k, const constexpr_array < y_type, order + 1 > & knots) const
    {
        y_type val = knots[k];
        val *= 0.5L;
        val += 0.5L;
        return val;
    }


    constexpr double JacobiP(double r, double alpha, double beta, double m) const
    {
        return 1;
    }




};



#endif