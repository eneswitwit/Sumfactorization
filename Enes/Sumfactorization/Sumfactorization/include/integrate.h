#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

template <int order, typename y_type, template<int, typename> class Polynomial, template<int, typename> class Quadrature >
class Integrate {
public:

  constexpr int vmult_mass(std::array < std::array < y_type, order + 1 >, order + 1 > &y, const std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
    Quadrature<order, y_type> quad = {1.};
    quad.compute_quadrature_points(order + 1, 1, 1);
    quad.compute_quadrature_weights(quad.knots, 0, 0);

    std::array < y_type, order + 1 > knots = quad.knots;
    std::array < y_type, order + 1 > weights = quad.weights;

    Polynomial <order, y_type> lagr;

    // evaluate all basis functions in 1D at all quadrature points.
    // row index: index of basis function
    // column index: index of quadrature point

    std::array < std::array < y_type, order + 1 >, order + 1 > eval_poly;

    for (unsigned int i=0;i<order+1;i++) {
            for (unsigned int j=0;j<order+1;j++) {
                    eval_poly[i][j]=lagr.eval_lagrange(i, quad.knots[j], quad.knots);
                }
        }

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
                                                    val_j += eval_poly[j][qy] * u[i][j];
                                                }
                                            val_qy += val_j * quad.weights[qy] * eval_poly[l][qy];
                                            counter++;
                                            counter++;
                                        }
                                    val_i += val_qy * eval_poly[i][qx];
                                    counter++;
                                }
                            val_qx += val_i * quad.weights[qx] * evalpoly[k][qx];
                            counter++;
                            counter++;
                        }
                    y[k][l] = val_qx;
                }
        }
    return counter;
  }

  constexpr int vmult_laplacian(std::array < std::array < y_type, order + 1 >, order + 1 > &y, const std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
      Quadrature<order, y_type> quad = {1.};
      quad.compute_quadrature_points(order + 1, 1, 1);
      quad.compute_quadrature_weights(quad.knots, 0, 0);

      std::array < y_type, order + 1 > knots = quad.knots;
      std::array < y_type, order + 1 > weights = quad.weights;

      Polynomial <order, y_type> lagr;

      // evaluate all basis functions and their derivatives in 1D at all quadrature points.
      // row index: index of basis function
      // column index: index of quadrature point
      // third index: stores values for 0, stores values of gradient for 1

      // dimension will be a template parameter later on
      unsigned int dim=2;

      std::array < std::array < std::array < y_type, order + 1 >, order + 1 >, 2 > eval_poly;

      for (unsigned int i=0;i<order+1;i++) {
              for (unsigned int j=0;j<order+1;j++) {
                      eval_poly[i][j][0]=lagr.eval_lagrange(i, quad.knots[j], quad.knots);
                      eval_poly[i][j][1]=lagr.eval_gradient(i, quad.knots[j], quad.knots);
                  }
          }
      y_type val_qx, val_qy, val_i, val_j;

      int counter = 0;
      for (int k = 0; k <= order; k++)
          {
              for (int l = 0; l <= order; l++)
                  {
                      y[k][l]=0;
                      for (unsigned int d=0;d<dim;d++) {
                              val_qx = 0;
                              for (int qx = 0; qx <= order; qx++)
                                  {
                                      val_i = 0;
                                      for (int i = 0; i <= order; i++)
                                          {
                                              val_qy = 0;
                                              for (int qy = 0; qy <= order; qy++)
                                                  {
                                                      val_j = 0;
                                                      for (int j = 0; j <= order; j++)
                                                          {
                                                              val_j += eval_poly[j][qy][dim-1-d] * u[i][j];
                                                          }
                                                      val_qy += val_j * quad.weights[qy] * eval_poly[l][qy][dim-1-d];
                                                      counter++;
                                                      counter++;
                                                  }
                                              val_i += val_qy * eval_poly[i][qx][d];
                                              counter++;
                                          }
                                      val_qx += val_i * quad.weights[qx] * evalpoly[k][qx][d];
                                      counter++;
                                      counter++;
                                  }
                              y[k][l] += val_qx;
                          }
                  }
          }
      return counter;
  }

};

#endif
