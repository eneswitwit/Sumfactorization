#include <QCoreApplication>
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

std::vector<double> eval_lagr(double x, int order) {
  std::vector<double> knots(order);
  knots=generate_knots(order);
  std::vector<double> val(order+1,1);
  for (int i=0;i<=order;i++) {
      for (int j=0;j<=order;j++) {
          if (i!=j) {
          val[i]*=(x-knots[j])/(knots[i]-knots[j]);
        }
      }
    }
  return val;
}

double test_with(std::vector<std::vector<double>> u, int k, int l, int order) {
  std::vector<double> knots(order);
  knots=generate_knots(order);
  std::vector<double> weights(3);
  weights[0]=1.0/3.0;
  weights[1]=4.0/3.0;
  weights[2]=1.0/3.0;
  std::vector<double> sf_eval_x(order+1);
  std::vector<double> sf_eval_y(order+1);
  double val=0;
  for (int qx=0; qx <= order; qx++) {
      for (int qy=0; qy <= order; qy++) {
          sf_eval_x=eval_lagr(knots[qx],order);
          sf_eval_y=eval_lagr(knots[qy],order);
          double sum=0;
          for (int i=0;i <= order ; i++) {
              for (int j=0;j <= order ; j++) {
                  sum+=sf_eval_x[i]*sf_eval_x[k]*sf_eval_y[j]*sf_eval_y[l]*u[i][j];
                }
            }
          val+=sum*weights[qx]*weights[qy];
          std::cout << qx << "  " << qy << "   "  <<val<< std::endl;
        }
    }
  return val;
}

int main()
{
  int order=2;
  std::vector<std::vector<double>> u(order+1, std::vector<double>(order+1,1));
  double val=test_with(u,0,1,order);
  std::cout << val;
  return 0;
}
