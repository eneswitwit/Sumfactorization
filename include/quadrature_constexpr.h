#ifndef __QUADRATURE_CONSTEXPR_H__
#define __QUADRATURE_CONSTEXPR_H__


#define PI 3.14159265359

constexpr int alpha = 1;
constexpr int beta = 1;

template <typename y_type, size_t size>
constexpr std::array<y_type, size> create_vector()
{
    std::array<y_type, size> arr{1.};
    for (unsigned int i = 0; i < size; ++i)
        arr[i] = i;
    return arr;
}

template <typename y_type, size_t size>
constexpr std::array<std::array<y_type, size>, size> create_array()
{
    std::array<std::array<y_type, size>, size> arr{1.};
    for (unsigned int i = 0; i < size; ++i)
        for (unsigned int j = 0; j < size; ++j)
            arr[i][j] = 1;
    return arr;
}

// Template class Quadrature for computing knots and weights of the quadrature.
template<int order, typename y_type>
class Quadrature {
public:
    constexpr std::array < y_type, order + 1 > knots;

    constexpr Quadrature(){
        //compute_quadrature_points();
        knots = compute_quadrature_points(1);
    }

    constexpr void compute_quadrature_points()
    {
        constexpr unsigned int m = order - 1;
        constexpr long double
        long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
        double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
        constexpr y_type runtime_one = 1.0;
        constexpr y_type tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );

        for (unsigned int i = 0; i < m; ++i)
        {
            knots[i + 1] = - std::cos( (y_type) (2 * i + 1) / (2 * m) * PI );
        }
        y_type s, J_x, f, delta;

        for (unsigned int k = 1; k < m; ++k)
        {
            y_type r = knots[k];
            if (k > 1)
                r = (r + knots[k - 1]) / 2;

            do
            {
                s = 0.;
                for (unsigned int i = 1; i < k; ++i)
                    s += 1. / (r - knots[i]);

                J_x   =  0.5 * (alpha + beta + m + 1) * JacobiP(r, alpha + 1, beta + 1, m - 1);
                f     = JacobiP(r, alpha, beta, m);
                delta = f / (f * s - J_x);
                r += delta;
            }
            while (std::fabs(delta) >= tolerance);

            knots[k] = r;
        }
        knots[0] = -1.L;
        knots[order] = +1.L;
        for (unsigned int j = 0; j < order + 1; j++) {
            knots[j] *= 0.5;
            knots[j] += 0.5;
        }

    }

    constexpr std::array<y_type, order+1> compute_quadrature_points(int dummy)
    {
        std::array<y_type, order+1> knots_constexpr {1.};
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

    }

    constexpr double JacobiP(double r,double alpha,double beta,double m){
        return 1;
    }




};



#endif