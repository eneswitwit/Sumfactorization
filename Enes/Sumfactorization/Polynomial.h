#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

template <int order, typename y_type>
class Integrate {
private:
    std::array < y_type, order + 1 > knots;
    std::array < y_type, order + 1 > weights;

public:
    // Constructor
    constexpr Polynomial(
        const std::array < y_type, order + 1 > &knots_ = {1.},
        const std::array < y_type, order + 1 > &weights_ = {1.})
        : knots(knots_),
          weights(weights_)
    {};

    constexpr y_type eval_lagr(int i, y_type x) const {
        y_type val = 1.;
        for (int j = 0; j <= order; j++) {
            if (i != j) {
                val *= (x - knots[j]) / (knots[i] - knots[j]);
            }
        }
        return val;
    }

};

endif