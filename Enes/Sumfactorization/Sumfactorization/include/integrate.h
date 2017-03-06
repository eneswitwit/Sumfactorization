#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

template <int order, typename y_type, template<int, typename> class Polynomial, template<int, typename> class Quadrature >
class Integrate {
public:

  constexpr int integrate_lagrange(std::array < std::array < y_type, order + 1 >, order + 1 > &y, const std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
    Quadrature<order, y_type> quad = {1.};
    quad.compute_quadrature_points(order + 1, 1, 1);
    quad.compute_quadrature_weights(quad.knots, 0, 0);

    std::array < y_type, order + 1 > knots = quad.knots;
    std::array < y_type, order + 1 > weights = quad.weights;

    Polynomial <order, y_type> lagr;
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
                val_j += lagr.eval_lagr(j, quad.knots[qy], quad.knots) * u[i][j];
              }
              val_qy += val_j * quad.weights[qy] * lagr.eval_lagr(l, quad.knots[qy], quad.knots);
              counter++;
              counter++;
            }
            val_i += val_qy * lagr.eval_lagr(i, quad.knots[qx], quad.knots);
            counter++;
          }
          val_qx += val_i * quad.weights[qx] * lagr.eval_lagr(k, quad.knots[qx], quad.knots);
          counter++;
          counter++;
        }
        y[k][l] = val_qx;
      }
    }
    return counter;
  };
};

#endif
