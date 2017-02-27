#include <iostream>
#include <vector>

std::vector<double> generate_knots(double order)
{
  std::vector<double> knots(order+1);
  for (int i=0; i<=order; i++) {
      knots[i]=i/(order);
    }
  return knots;
}

double eval_lagr(int i, double x, int order) {
  std::vector<double> knots(order);
  knots=generate_knots(order);
  double val;
      for (int j=0;j<=order;j++) {
          if (i!=j) {
          val*=(x-knots[j])/(knots[i]-knots[j]);
        }
      }
  return val;
}

std::vector<std::vector<double>> integrate(std::vector<std::vector<double>> u, int order) {
  std::vector<double> knots(order);
  knots=generate_knots(order);
  std::vector<double> weights(order);
  weights=generate_knots(order);
  std::vector<std::vector<double>> y(order+1,std::vector<double> (order+1));
  int counter=0;
  for (int k=0; k<=order; k++)
    {
      for (int l=0; l<=order; l++)
        {
          double val_qx=0;
          for (int qx=0; qx<=order; qx++)
            {
              double val_i=0;
              for (int i=0; i<=order; i++)
                {
                  double val_qy=0;
                  for (int qy=0; qy<=order; qy++)
                    {
                      double val_j=0;
                      for (int j=0; j<=order; j++)
                        {
                          val_j+=eval_lagr(j,knots[qy],order)*u[i][j];
                        }
                      val_qy+=val_j*weights[qy]*eval_lagr(l,knots[qy],order);
                      counter++;
                      counter++;
                    }
                  val_i+=val_qy*eval_lagr(i,knots[qx],order);
                  counter++;
                }
              val_qx+=val_i*weights[qx]*eval_lagr(k,knots[qx],order);
              counter++;
              counter++;
            }
          y[k][l]=val_qx;
        }
    }
  std::cout << counter << std::endl;
  return y;
}

int main()
{
  int order=3;
  std::vector<std::vector<double>> u(order+1, std::vector<double>(order+1,1));
  std::vector<std::vector<double>> val=integrate(u,order);
  return 0;
}
