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
    constexpr y_type eval_gradient(int i, y_type x, std::array < y_type, order + 1 > & knots) const {
        y_type sum=0;
        if(std::find(knots.begin(), knots.end(), x) != knots.end()) {
                // evaluate gradient on a knot
                for (unsigned int m=0;m<order+1;m++) {
                        if (m!=i) {
                                y_type product=1;
                                for (unsigned int j=0;j<order+1;j++) {
                                        if ((j!=i) && (j!=m)) {
                                                product*=(x-knots[j])/(knots[i]-knots[j]);
                                            }
                                    }
                                sum+=(knots[i]-knots[m])*product;
                            }
                    }
                return sum;
            } else {
                // evaluate gradient in between knots
                for (unsigned int j=0;j<order+1;j++) {
                        if (j!=i) {
                                sum+=1/(x-knots[j]);
                            }
                    }
                return sum*eval_lagr(i,x,& knots);
            }

    }

};

#endif
