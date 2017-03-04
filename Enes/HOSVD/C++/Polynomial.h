#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

template <int order, typename y_type>
class Polynomial {
public:
    std::vector < y_type > knots;
    std::vector < y_type > weights;

    // Constructor
    constexpr Polynomial( std::vector < y_type > knots_)
    {
        knots = knots_ ;
        weights = compute_quadrature_weights(knots_, 0, 0) ;
    };

    y_type eval_lagr(int i, y_type x) const {
        y_type val = 1.;
        for (int j = 0; j <= order; j++) {
            if (i != j) {
                val *= (x - knots[j]) / (knots[i] - knots[j]);
            }
        }
        return val;
    }

};

#endif