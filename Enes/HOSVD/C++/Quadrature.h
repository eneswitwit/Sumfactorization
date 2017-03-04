#ifndef __QUADRATURE_H__
#define __QUADRATURE_H__

std::vector<long double> compute_quadrature_points(const unsigned int q,
    const int alpha,
    const int beta)
{
  const unsigned int m = q - 2; // no. of inner points
  const double PI = 3.14159265359;
  std::vector<long double> x(m);

  // compute quadrature points with
  // a Newton algorithm.

  // Set tolerance. See class QGauss
  // for detailed explanation.
  const long double
  long_double_eps = static_cast<long double>(std::numeric_limits<long double>::epsilon()),
  double_eps      = static_cast<long double>(std::numeric_limits<double>::epsilon());

  // check whether long double is
  // more accurate than double, and
  // set tolerances accordingly
  volatile long double runtime_one = 1.0;
  const long double tolerance
    = (runtime_one + long_double_eps != runtime_one
       ?
       std::max (double_eps / 100, long_double_eps * 5)
       :
       double_eps * 5
      );

  // The following implementation
  // follows closely the one given in
  // the appendix of the book by
  // Karniadakis and Sherwin:
  // Spectral/hp element methods for
  // computational fluid dynamics
  // (Oxford University Press, 2005)

  long double JacobiP(const long double x, const int alpha, const int beta, const unsigned int n);
  // we take the zeros of the Chebyshev
  // polynomial (alpha=beta=-0.5) as
  // initial values:
  for (unsigned int i = 0; i < m; ++i)
    x[i] = - std::cos( (long double) (2 * i + 1) / (2 * m) * PI );

  long double s, J_x, f, delta;

  for (unsigned int k = 0; k < m; ++k)
  {
    long double r = x[k];
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
  } // for

  // add boundary points:
  x.insert(x.begin(), -1.L);
  x.push_back(+1.L);
  for (unsigned int j = 0; j < x.size(); j++) {
    x[j] *= 0.5;
    x[j] += 0.5;
  }

  return x;
}

std::vector<long double> compute_quadrature_weights(std::vector<long double> x, const int alpha, const int beta)
{
  for (unsigned int j = 0; j < x.size(); j++) {
    x[j] *= 2;
    x[j] += -1;
  }
  long double gamma(const unsigned int n);
  long double JacobiP(const long double x, const int alpha, const int beta, const unsigned int n);
  const unsigned int q = x.size();
  std::vector<long double> w(q);
  long double s = 0.L;

  const long double factor = std::pow(2., alpha + beta + 1) *
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

  return w;
}


long double JacobiP(const long double x, const int alpha, const int beta, const unsigned int n)
{
  // the Jacobi polynomial is evaluated
  // using a recursion formula.
  std::vector<long double> p(n + 1);

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

    p[i + 1] = static_cast<long double>( (a2 + a3 * x) * p[i] - a4 * p[i - 1]) / a1;
  } // for
  return p[n];
}

long double gamma(const unsigned int n)
{
  long double result = n - 1;
  for (int i = n - 2; i > 1; --i)
    result *= i;
  return result;
}

#endif