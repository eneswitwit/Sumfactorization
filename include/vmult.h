#ifndef __VMULT_H__
#define __VMULT_H__

template <int order, typename y_type, template<int, typename> class Polynomial, template<int, typename> class Quadrature >
class VMULT {
public:
    constexpr void vmult_mass(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
        // Quadrature Properties
        Quadrature<order, y_type> quad = {1.};
        quad.compute_quadrature_points(order + 1, 1, 1);
        quad.compute_quadrature_weights(quad.knots, 0, 0);
        std::array < y_type, order + 1 > knots = quad.knots;
        std::array < y_type, order + 1 > weights = quad.weights;

        // Polynomial Basis
        Polynomial <order, y_type> poly;

        // Build Matrix N
        std::array < std::array < y_type, order + 1 >, order + 1 > N;
        for (unsigned int i = 0; i < order + 1; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                N[i][j] = poly.eval_lagrange(j, knots[i], knots);
            }
        }

        // Build Matrix N_w
        std::array < std::array < y_type, order + 1 >, order + 1 > NW;
        for (unsigned int i = 0; i < order + 1; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NW[j][i] = weights[i] * poly.eval_lagrange(j, knots[i], knots);
            }
        }

        // Calculate N_w^T * N
        std::array < std::array < y_type, order + 1 >, order + 1 > C_1;
        std::array < std::array < y_type, order + 1 >, order + 1 > C_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NW, N, C_1);
        multiply_matrices<order+1,order+1,order+1,order+1,order+1,order+1,y_type>(C_1,u,C_2);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (C_2, C_1, y);

    }



};

#endif