#ifndef __VMULT_H__
#define __VMULT_H__

template <int order, typename y_type, template<int, typename> class Polynomial, template<int, typename> class Quadrature >
class VMULT {
public:
    constexpr void vmult_mass(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
        // Quadrature Properties
        Quadrature<order, y_type> quad;
        std::array < y_type, order + 1 > knots = quad.knots;
        std::array < y_type, order + 1 > weights  = quad.weights;
        //std::array < y_type, order + 1 > knots = quad.knots;
        //std::array < y_type, order + 1 > weights = quad.weights;

        // Polynomial Basis
        Polynomial <order, y_type> poly;

        // Build Matrix N and NW
        std::array < std::array < y_type, order + 1 >, order + 1 > N;
        std::array < std::array < y_type, order + 1 >, order + 1 > NW_transposed;
        for (unsigned int i = 0; i < order + 1; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                        N[i][j] = poly.eval_lagrange(j, knots[i], knots);
                        NW_transposed[j][i] = weights[i] * N[i][j];
                    }
            }


        // Calculate N_w^T * N
        std::array < std::array < y_type, order + 1 >, order + 1 > C_1;
        std::array < std::array < y_type, order + 1 >, order + 1 > C_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NW_transposed, N, C_1);
        multiply_matrices<order+1,order+1,order+1,order+1,order+1,order+1,y_type>(C_1,u,C_2);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (C_2, C_1, y,1);

    }

    constexpr void vmult_gradient(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
        Quadrature<order, y_type> quad;
        Quadrature<order,y_type> interpolation;

        // Polynomial Basis
        Polynomial <order, y_type> poly;

        // Build Matrix N and NW
        std::array < std::array < y_type, order + 1 >, order + 1 > N;
        std::array < std::array < y_type, order + 1 >, order + 1 > NW_transposed;
        std::array < std::array < y_type, order + 1 >, order + 1 > Ndx;
        std::array < std::array < y_type, order + 1 >, order + 1 > NWdx_transposed;
        for (unsigned int i = 0; i < order + 1; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                        N[i][j] = poly.eval_lagrange(j, quad.knots[i], interpolation.knots);
                        NW_transposed[j][i] = quad.weights[i] * N[i][j];
                        Ndx[i][j] = poly.eval_1st_derivative(j,quad.knots[i], interpolation.knots);
                        NWdx_transposed[j][i] = quad.weights[i] * Ndx[i][j];
                    }
            }
        std::array < std::array < y_type, order + 1 >, order + 1 > P;
        std::array < std::array < y_type, order + 1 >, order + 1 > Pdx;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NW_transposed, N, P);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NWdx_transposed, Ndx, Pdx);

        std::array < std::array < y_type, order + 1 >, order + 1 > S_1;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Pdx, u, y);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (y, P, S_1, 1);

        std::array < std::array < y_type, order + 1 >, order + 1 > S_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (P, u, y);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (y, Pdx, S_2, 1);

        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (S_1, S_2, y);
    }

    constexpr void vmult_laplacian(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) const {
        Quadrature<order, y_type> quad;
        Quadrature<order,y_type> interpolation;

        // Polynomial Basis
        Polynomial <order, y_type> poly;

        // Build Matrix N and NW
        std::array < std::array < y_type, order + 1 >, order + 1 > N;
        std::array < std::array < y_type, order + 1 >, order + 1 > NW_transposed;
        std::array < std::array < y_type, order + 1 >, order + 1 > Ndx;
        std::array < std::array < y_type, order + 1 >, order + 1 > NWdx_transposed;
        for (unsigned int i = 0; i < order + 1; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                        N[i][j] = poly.eval_lagrange(j, quad.knots[i], interpolation.knots);
                        NW_transposed[j][i] = quad.weights[i] * N[i][j];
                        Ndx[i][j] = poly.eval_2nd_derivative(j,quad.knots[i], interpolation.knots);
                    }
            }
        std::array < std::array < y_type, order + 1 >, order + 1 > P;
        std::array < std::array < y_type, order + 1 >, order + 1 > Pdx;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NW_transposed, N, P);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (NW_transposed, Ndx, Pdx);

        std::array < std::array < y_type, order + 1 >, order + 1 > S_1;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Pdx, u, y);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (y, P, S_1, 1);

        std::array < std::array < y_type, order + 1 >, order + 1 > S_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (P, u, y);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (y, Pdx, S_2, 1);

        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (S_1, S_2, y);
    }

};

#endif
