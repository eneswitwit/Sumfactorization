#ifndef __QUADRATURE_NONCLASS_H__
#define __QUADRATURE_NONCLASS_H__

#define PI 3.14159265359

// Helper function for transforming vector to array
template<typename y_type, int order>
std::array<y_type, order> vec_to_arr(std::vector<y_type> vec) {
  std::array<y_type, order> arr;
  for (unsigned int i = 0; i < order ; ++i) {
    arr[i] = vec[i];
  }
  return arr;
}

constexpr int alpha = 1;
constexpr int beta = 1;

template<size_t order, typename y_type>
constexpr std::array < y_type, order + 1 > compute_quadrature_points()
{
  std::array < y_type, order + 1 > knots{1.};
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
  return knots;
}

void compute_quadrature_weights(std::array < y_type, order + 1 > & knots)
{
  const int alpha = 1;
  const int beta = 1;

  for (unsigned int j = 0; j < knots.size(); j++) {
    knots[j] *= 2;
    knots[j] += -1;
  }
  std::vector<y_type> w(order + 1);
  y_type s = 0.L;

  const y_type factor = std::pow(2., alpha + beta + 1) *
                        gamma(alpha + order + 1) *
                        gamma(beta + order + 1) /
                        ((order + 1 - 1) * gamma(order + 1) * gamma(alpha + beta + order + 1 + 1));
  for (unsigned int i = 0; i < order + 1; ++i)
  {
    s = JacobiP(knots[i], alpha, beta, order + 1 - 1);
    w[i] = factor / (s * s);
  }
  w[0]   *= (beta + 1);
  w[order] *= (alpha + 1);

  for (unsigned int i = 0; i < w.size(); ++i)
  {
    w[i] = w[i] * 0.5;
  }

  weights = vec_to_arr < y_type, order + 1 > (w);
}

y_type JacobiP(const y_type x, const int alpha, const int beta, const unsigned int n)
{
  std::vector<y_type> p(n + 1);
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
  }
  return p[n];
}

y_type gamma(const unsigned int n)
{
  y_type result = n - 1;
  for (int i = n - 2; i > 1; --i)
    result *= i;
  return result;
}


#endif