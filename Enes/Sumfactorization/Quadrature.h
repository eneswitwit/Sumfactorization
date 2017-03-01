#ifndef __QUADRATURE_H__
#define __QUADRATURE_H__

// Helper function for transforming vector to array
template<typename y_type, int order>
std::array<y_type, order> vec_to_arr(std::vector<y_type> vec) {
  std::array<y_type, order> arr;
  for (unsigned int i = 0; i < order ; ++i) {
    arr[i] = vec[i];
  }
  return arr;
}

template<int order, typename y_type>
class Quadrature {
public:
  std::array < y_type, order + 1 > knots;
  std::array < y_type, order + 1 > weights;

  void compute_quadrature_points(const unsigned int q, const int alpha, const int beta)
  {
    const unsigned int m = q - 2; // no. of inner points
    const double PI = 3.14159265359;
    std::vector<y_type> x(m);
    const long double
    long_double_eps = static_cast<y_type>(std::numeric_limits<long double>::epsilon()),
    double_eps      = static_cast<y_type>(std::numeric_limits<double>::epsilon());
    volatile y_type runtime_one = 1.0;
    const y_type tolerance
      = (runtime_one + long_double_eps != runtime_one
         ?
         std::max (double_eps / 100, long_double_eps * 5)
         :
         double_eps * 5
        );
    y_type JacobiP(const y_type x, const int alpha, const int beta, const unsigned int n);
    for (unsigned int i = 0; i < m; ++i)
      x[i] = - std::cos( (y_type) (2 * i + 1) / (2 * m) * PI );
    y_type s, J_x, f, delta;

    for (unsigned int k = 0; k < m; ++k)
    {
      y_type r = x[k];
      if (k > 0)
        r = (r + x[k - 1]) / 2;

      do
      {
        s = 0.;
        for (unsigned int i = 0; i < k; ++i)
          s += 1. / (r - x[i]);

        J_x   =  0.5 * (alpha + beta + m + 1) * JacobiP(r, alpha + 1, beta + 1, m - 1);
        f     = JacobiP(r, alpha, beta, m);
        delta = f / (f * s - J_x);
        r += delta;
      }
      while (std::fabs(delta) >= tolerance);

      x[k] = r;
    }
    x.insert(x.begin(), -1.L);
    x.push_back(+1.L);
    for (unsigned int j = 0; j < x.size(); j++) {
      x[j] *= 0.5;
      x[j] += 0.5;
    }

    knots = vec_to_arr < y_type, order + 1 > (x);
  }



  void compute_quadrature_weights(std::array<y_type,order+1> x, const int alpha, const int beta)
  {
    for (unsigned int j = 0; j < x.size(); j++) {
      x[j] *= 2;
      x[j] += -1;
    }
    y_type gamma(const unsigned int n);
    y_type JacobiP(const y_type x, const int alpha, const int beta, const unsigned int n);
    const unsigned int q = x.size();
    std::vector<y_type> w(q);
    y_type s = 0.L;

    const y_type factor = std::pow(2., alpha + beta + 1) *
                               gamma(alpha + q) *
                               gamma(beta + q) /
                               ((q - 1) * gamma(q) * gamma(alpha + beta + q + 1));
    for (unsigned int i = 0; i < q; ++i)
    {
      s = JacobiP(x[i], alpha, beta, q - 1);
      w[i] = factor / (s * s);
    }
    w[0]   *= (beta + 1);
    w[q - 1] *= (alpha + 1);

    for (unsigned int i = 0; i < w.size(); ++i)
    {
      w[i] = w[i] * 0.5;
    }

    weights = vec_to_arr<y_type,order+1>(w);
  }


  y_type JacobiP(const y_type x, const int alpha, const int beta, const unsigned int n)
  {
    // the Jacobi polynomial is evaluated
    // using a recursion formula.
    std::vector<y_type> p(n + 1);

    // initial values P_0(x), P_1(x):
    p[0] = 1.0L;
    if (n == 0) return p[0];
    p[1] = ((alpha + beta + 2) * x + (alpha - beta)) / 2;
    if (n == 1) return p[1];

    for (unsigned int i = 1; i <= (n - 1); ++i)
    {
      const int v  = 2 * i + alpha + beta;
      const int a1 = 2 * (i + 1) * (i + alpha + beta + 1) * v;
      const int a2 = (v + 1) * (alpha * alpha - beta * beta);
      const int a3 = v * (v + 1) * (v + 2);
      const int a4 = 2 * (i + alpha) * (i + beta) * (v + 2);

      p[i + 1] = static_cast<y_type>( (a2 + a3 * x) * p[i] - a4 * p[i - 1]) / a1;
    } // for
    return p[n];
  }

  y_type gamma(const unsigned int n)
  {
    y_type result = n - 1;
    for (int i = n - 2; i > 1; --i)
      result *= i;
    return result;
  }


};
#endif
