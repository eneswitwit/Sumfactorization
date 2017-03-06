#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

template <int order, typename y_type>
class Polynomial {
public:

    constexpr y_type eval_lagr(int i, y_type x, const std::array < y_type, order + 1 > & knots) const {
        y_type val = 1.;
        for (int j = 0; j <= order; j++) {
            if (i != j) {
                val *= (x - knots[j]) / (knots[i] - knots[j]);
            }
        }
        return val;
    }

    constexpr y_type approx_gradient(int i, y_type x, std::array < y_type, order + 1 > & knots) const {
        y_type h = 0;
        h = sqrt(std::nextafter(h, h));
        if (x != 0) {
            h *= x;
        }
        y_type val = (eval_lagr(i, x + h, & knots) - eval_lagr(i, x + h, & knots)) / (2 * h);
    }

};

#endif
