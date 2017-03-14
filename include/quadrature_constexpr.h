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
public:
    constexpr constexpr_array < y_type, order + 1 > knots_;
    constexpr Quadrature(constexpr_array < y_type, order > quad_sin_vec) {
        knots_ = compute_quadrature_points(quad_sin_vec);
        //knots = compute_quadrature_points(1);
    }

    constexpr constexpr_array < y_type, order + 1 > compute_quadrature_points(constexpr_array < y_type, order >  quad_sin_vec)
    {
        constexpr_array < y_type, order + 1 > knots;
        constexpr unsigned int m = order - 1;
        constexpr long double
        long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
        double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
        constexpr y_type runtime_one = 1.0;
        constexpr y_type tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );
        for (unsigned int i = 0; i < m; ++i)
        {
            knots[i + 1] = quad_sin_vec[i] ;
        }
        y_type s, J_x, f, delta;

        for (unsigned int k = 1; k < m; ++k)
        {
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
            while (fabs(delta) >= tolerance);

            knots[k] = r;
        }
        knots[0] = -1.L;
        knots[order] = +1.L;
        for (unsigned int j = 0; j < order + 1; j++) {
            knots[j] *= 0.5;
            knots[j] += 0.5;
        }
        return knots;

    }

    /*constexpr std::array < y_type, order + 1 > compute_quadrature_points(int dummy)
    {
        std::array < y_type, order + 1 > knots_constexpr {1.};
        constexpr unsigned int m = order - 1;
        constexpr long double
        long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
        double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
        constexpr y_type runtime_one = 1.0;
        constexpr y_type tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );

        for (unsigned int i = 0; i < m; ++i)
        {
            knots_constexpr[i + 1] = - std::cos( (y_type) (2 * i + 1) / (2 * m) * PI );
        }
        y_type s, J_x, f, delta;

        for (unsigned int k = 1; k < m; ++k)
        {
            y_type r = knots_constexpr[k];
            if (k > 1)
                r = (r + knots_constexpr[k - 1]) / 2;

            do
            {
                s = 0.;
                for (unsigned int i = 1; i < k; ++i)
                    s += 1. / (r - knots_constexpr[i]);

                J_x   =  0.5 * (alpha + beta + m + 1) * JacobiP(r, alpha + 1, beta + 1, m - 1);
                f     = JacobiP(r, alpha, beta, m);
                delta = f / (f * s - J_x);
                r += delta;
            }
            while (std::fabs(delta) >= tolerance);

            knots_constexpr[k] = r;
        }
        knots_constexpr[0] = -1.L;
        knots_constexpr[order] = +1.L;
        for (unsigned int j = 0; j < order + 1; j++) {
            knots_constexpr[j] *= 0.5;
            knots_constexpr[j] += 0.5;
        }

        return knots_constexpr;

    }*/

    constexpr double JacobiP(double r, double alpha, double beta, double m) {
        return 1;
    }




};



#endif