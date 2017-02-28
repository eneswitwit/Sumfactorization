#include <iostream>
#include <array>
#include <array>

template <int order, typename y_type >
class Polynomial {
private:
    std::array < y_type, order + 1 > knots;
    std::array < y_type, order + 1 > weights;
    std::array < std::array < y_type, order + 1 >, order + 1 > u;
    std::array < std::array < y_type, order + 1 >, order + 1 > y;


public:

    // Constructor

    Polynomial() {};

    Polynomial(std::array < std::array < y_type, order + 1 >, order + 1 > uu)
    {
        u = uu;
    };

    Polynomial(std::array < y_type, order + 1 > knots_ , std::array < y_type, order + 1 > weights_ )
    {
        knots = knots_;
        weights = weights_;
    };


    // Sum factorization

    void generate_knots()
    {
        for (int i = 0; i <= order; i++) {
            knots[i] = i / (order);
        }
        return knots;
    };

    y_type eval_lagr(int i, y_type x) {
        y_type val;
        for (int j = 0; j <= order; j++) {
            if (i != j) {
                val *= (x - knots[j]) / (knots[i] - knots[j]);
            }
        }
        return val;
    };

    void integrate() {
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
                                val_j += eval_lagr(j, knots[qy], order) * u[i][j];
                            }
                            val_qy += val_j * weights[qy] * eval_lagr(l, knots[qy], order);
                            counter++;
                            counter++;
                        }
                        val_i += val_qy * eval_lagr(i, knots[qx], order);
                        counter++;
                    }
                    val_qx += val_i * weights[qx] * eval_lagr(k, knots[qx], order);
                    counter++;
                    counter++;
                }
                y[k][l] = val_qx;
            }
        }
        std::cout << counter << std::endl;
        return y;
    };
};

int main()
{
    Polynomial<3, double> p;
    return 0;
}
