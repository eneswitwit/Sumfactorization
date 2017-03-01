#include <iostream>
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include "Quadrature.h"
#include "Polynomial.h"

using namespace std;

template <size_t size>
constexpr std::array<double, size> create_vector()
{
  std::array<double, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    arr[i] = i;
  return arr;
}

template <size_t size>
constexpr std::array<std::array<double, size>, size> create_array()
{
  std::array<std::array<double, size>, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      arr[i][j] = i + j;
  return arr;
}


template <int order, typename y_type, template<int order_p, typename y_type_p> class Polynomial >
class Integrate {

public:


  constexpr int integrate_lagrange(std::array < std::array < y_type, order + 1 >, order + 1 > &y, const std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
    Polynomial <order, y_type> lagr(compute_quadrature_points(order, 0, 0));
    int counter = 0;
    for (int k = 0; k <= order; k++)
    {
      for (int l = 0; l <= order; l++)
      {
        y_type val_qx = 0;
        for (int qx = 0; qx <= order; qx++)
        {
          y_type val_i = 0;
          for (int i = 0; i <= order; i++)
          {
            y_type val_qy = 0;
            for (int qy = 0; qy <= order; qy++)
            {
              y_type val_j = 0;
              for (int j = 0; j <= order; j++)
              {
                val_j += lagr.eval_lagr(j, lagr.knots[qy]) * u[i][j];
              }
              val_qy += val_j * lagr.weights[qy] * lagr.eval_lagr(l, lagr.knots[qy]);
              counter++;
              counter++;
            }
            val_i += val_qy * lagr.eval_lagr(i, lagr.knots[qx]);
            counter++;
          }
          val_qx += val_i * lagr.weights[qx] * lagr.eval_lagr(k, lagr.knots[qx]);
          counter++;
          counter++;
        }
        y[k][l] = val_qx;
      }
    }
    return counter;
  };
};

int main()
{
  /*constexpr int order = 3;
  constexpr std::array < double, order + 1 > weights = create_vector < order + 1 > ();
  constexpr std::array < double, order + 1 > knots = create_vector < order + 1 > ();
  constexpr std::array < std::array < double, order + 1 >, order + 1 > u = create_array < order + 1 > ();
  std::cout << u[order][order] << std::endl;
  static_assert(u[order][order] == 2 * order);

  constexpr Polynomial<order, double> p(u, knots, weights);
  constexpr int count = p.integrate();*/

  constexpr int order = 3;
  Integrate<3, double, Polynomial > lagrange;
  constexpr std::array < std::array < double, order + 1 >, order + 1 > u = create_array < order + 1 > ();
  constexpr std::array < std::array < double, order + 1 >, order + 1 > y = create_array < order + 1 > ();
  lagrange.integrate_lagrange(&y,&u);
  return count;
}
