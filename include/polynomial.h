#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

template <int order, typename y_type>
class Polynomial {
public:

    y_type eval_lagrange(int i, y_type x, const std::array < y_type, order + 1 > & knots) const {
        y_type val = 1.;
        for (int j = 0; j <= order; j++) {
                if (i != j) {
                        val *= (x - knots[j]) / (knots[i] - knots[j]);
                    }
            }
        return val;
    }

    y_type eval_1st_derivative(int i, y_type x, std::array < y_type, order + 1 > & knots) const {
        y_type sum=0;
        for (unsigned int m=0;m<order+1;m++) {
                if (m!=i) {
                        y_type product=1;
                        for (unsigned int j=0;j<order+1;j++) {
                                if ((j!=i) && (j!=m)) {
                                        product*=(x-knots[j])/(knots[i]-knots[j]);
                                    }
                            }
                        sum+=product/(knots[i]-knots[m]);
                    }
            }
        return sum;
    }

    y_type eval_2nd_derivative(int i, y_type x, std::array < y_type, order + 1 > & knots) const {
        y_type sum=0;
        for (unsigned int k=0;k<order+1;k++) {
                if (k!=i) {
                        y_type sum1=0;
                        for (unsigned int m=0;m<order+1;m++) {
                                if ((m!=i)&&(m!=k)) {
                                        y_type product=1;
                                        for (unsigned int j=0;j<order+1;j++) {
                                                if ((j!=i) && (j!=m) && (j!=k)) {
                                                        product*=(x-knots[j])/(knots[i]-knots[j]);
                                                    }
                                            }
                                        sum1+=product/(knots[i]-knots[m]);
                                    }
                            }
                        sum+=sum1/(knots[i]-knots[k]);
                    }
            }
        return sum;
    }



};

#endif
