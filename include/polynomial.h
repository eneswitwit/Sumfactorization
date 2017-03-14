#ifndef __POLYNOMIAL_H__
#define __POLYNOMIAL_H__

template <int order, typename y_type>
class Polynomial {
public:

    std::vector <y_type, order+1> knots;

    Polynomial() {
        compute_interpolation_points();
    }

    void compute_interpolation_points()
    {
      const unsigned int m = order - 1;
      const long double
      long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
      double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
      const y_type runtime_one = 1.0;
      const y_type tolerance = (runtime_one + long_double_eps != runtime_one ? std::max (double_eps / 100, long_double_eps * 5) : double_eps * 5 );

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

    y_type eval_lagrange(int i, y_type x) const {
        y_type val = 1.;
        for (int j = 0; j <= order; j++) {
                if (i != j) {
                        val *= (x - knots[j]) / (knots[i] - knots[j]);
                    }
            }
        return val;
    }

    y_type eval_1st_derivative(int i, y_type x) const {
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

    y_type eval_2nd_derivative(int i, y_type x) const {
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
