#include <iostream>
#include <array>

template <size_t size>
constexpr std::array<double, size> create_vector()
{ 
  std::array<double, size> arr{1.};
  for (unsigned int i=0; i<size; ++i)
    arr[i] = i;
  return arr;
}

template <size_t size>
constexpr std::array<std::array<double, size>, size> create_array()
{ 
  std::array<std::array<double, size>, size> arr{1.};
  for (unsigned int i=0; i<size; ++i)
    for (unsigned int j=0; j<size; ++j)
      arr[i][j] = i+j;
  return arr;
}


template <int order, typename y_type >
class Polynomial {
private:
    std::array < y_type, order + 1 > knots;
    std::array < y_type, order + 1 > weights;
    std::array < std::array < y_type, order + 1 >, order + 1 > u;
    std::array < std::array < y_type, order + 1 >, order + 1 > y={1.};


public:
    // Constructor
    constexpr Polynomial(const std::array < std::array < y_type, order + 1 >, order + 1 > &u_,
                        const std::array < y_type, order + 1 > &knots_={1.},
                        const std::array < y_type, order + 1 > &weights_={1.})
      :knots(knots_),
  weights(weights_),
  u(u_)
    {};    

    constexpr y_type eval_lagr(int i, y_type x) const {
        y_type val=0.;
        for (int j = 0; j <= order; j++) {
            if (i != j) {
                val *= (x - knots[j]) / (knots[i] - knots[j]);
            }
        }
        return val;
    };

    constexpr int integrate() const {
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
                                val_j += eval_lagr(j, knots[qy]) * u[i][j];
                            }
                            val_qy += val_j * weights[qy] * eval_lagr(l, knots[qy]);
                            counter++;
                            counter++;
                        }
                        val_i += val_qy * eval_lagr(i, knots[qx]);
                        counter++;
                    }
                    val_qx += val_i * weights[qx] * eval_lagr(k, knots[qx]);
                    counter++;
                    counter++;
                }
                //y[k][l] = val_qx;
            }
        }
        //std::cout << counter << std::endl;
        return counter;
    };
};

int main()
{
  constexpr int order=3;
  constexpr std::array <double, order + 1 > weights=create_vector<order+1>();
  constexpr std::array <double, order + 1 > knots=create_vector<order+1>();   
  constexpr std::array<std::array<double, order+1>, order+1> u = create_array<order+1>();  
  std::cout << u[order][order] << std::endl;
  static_assert(u[order][order]==2*order);  

  constexpr Polynomial<order, double> p(u, knots, weights);
  constexpr int count = p.integrate();
  return count;  
}
