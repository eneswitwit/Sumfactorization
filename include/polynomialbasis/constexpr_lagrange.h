#ifndef __CONSTEXPR_LAGRANGE_H__
#define __CONSTEXPR_LAGRANGE_H__


// A class for handling polynomials with lagrangian basis.

template <typename y_type, size_t order, template<typename, size_t> class Quadrature>
class Lagrange {
public:
 
   const constexpr_array <y_type, order+1> knots_;

   constexpr  Lagrange() : knots_(compute_knots()) {}

    constexpr constexpr_array <y_type, order + 1> compute_knots() const {
        constexpr Quadrature<y_type,order> quad;
        return quad.compute_quadrature_points();
    }


    constexpr y_type eval_lagrange(const y_type & x, const unsigned int i) const {
        y_type val = 1.;
        for (unsigned int j = 0; j <= order; j++) {
                if (i != j) {
                        val *= (x - knots_[j]) / (knots_[i] - knots_[j]);
                    }
            }
        return val;
    }

    constexpr y_type eval_1st_derivative(const y_type & x, const unsigned int i) const {
        y_type sum=0;
        for (unsigned int m=0;m<order+1;m++) {
                if (m!=i) {
                        y_type product=1;
                        for (unsigned int j=0;j<order+1;j++) {
                                if ((j!=i) && (j!=m)) {
                                        product*=(x-knots_[j])/(knots_[i]-knots_[j]);
                                    }
                            }
                        sum+=product/(knots_[i]-knots_[m]);
                    }
            }
        return sum;
    }

    constexpr y_type eval_2nd_derivative(const y_type & x, const unsigned int i) const {
        y_type sum=0;
        for (unsigned int k=0;k<order+1;k++) {
                if (k!=i) {
                        y_type sum1=0;
                        for (unsigned int m=0;m<order+1;m++) {
                                if ((m!=i)&&(m!=k)) {
                                        y_type product=1;
                                        for (unsigned int j=0;j<order+1;j++) {
                                                if ((j!=i) && (j!=m) && (j!=k)) {
                                                        product*=(x-knots_[j])/(knots_[i]-knots_[j]);
                                                    }
                                            }
                                        sum1+=product/(knots_[i]-knots_[m]);
                                    }
                            }
                        sum+=sum1/(knots_[i]-knots_[k]);
                    }
            }
        return sum;
    }



};

#endif
