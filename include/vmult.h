#ifndef __VMULT_H__
#define __VMULT_H__

template <typename y_type, size_t order, size_t q_order, int c, template<typename, size_t> class Quadrature, template<typename, size_t, template<typename, size_t> class Quadrature > class Polynomial>
class VMULT {
public:

    const Quadrature < y_type, q_order + 1 > quad;
    const Polynomial <y_type, order, Quadrature> poly;

    const constexpr_array < constexpr_array < y_type, q_order + 1 >, order + 1 > N_Transposed;
    const constexpr_array < constexpr_array < y_type, order + 1 >, q_order + 1 > NW;
    const constexpr_array < constexpr_array < y_type, order + 1 >, order + 1 > N_Product;

    const constexpr_array < constexpr_array < y_type, q_order + 1 >, order + 1 > Ndx_Transposed;
    const constexpr_array < constexpr_array < y_type, order + 1 >, q_order + 1 > NWdx;
    const constexpr_array < constexpr_array < y_type, order + 1 >, order + 1 > Ndx_Product;

    const constexpr_array < constexpr_array < y_type, q_order + 1 >, order + 1 > Ndxdx_Transposed;
    const constexpr_array < constexpr_array < y_type, order + 1 >, order + 1 > Ndxdx_Product;

    constexpr VMULT() :
        quad(),
        poly(),
        N_Transposed(compute_basis_matrix_NT()),
        NW(),
        N_Product(),
        Ndx_Transposed(),
        NWdx(),
        Ndx_Product(),
        Ndxdx_Transposed(),
        Ndxdx_Product() {}

    constexpr constexpr_array < constexpr_array < y_type, q_order + 1 >, order + 1 > compute_basis_matrix_NT() {
        constexpr_array < constexpr_array < y_type, q_order + 1 >, order + 1 > NT;
        for (unsigned int i = 0; i < q_order+1 ; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NT[i][j] = poly.eval_lagrange(quad[i],j);
            }
        }
        return NT;

    }

    /*
    void compute_basis_matrix_NW() {
        for (unsigned int i = 0; i < q_order; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                N_Transposed[i][j] = poly.eval_lagrange(j, quad.knots[i]);
                NW[j][i] = quad.weights[i] * N_Transposed[i][j];
            }
        }
        multiply_matrices < order + 1, q_order, q_order, order + 1, order + 1, order + 1, y_type > (NW, N_Transposed, N_Product);

        for (int i = 0; i < order + 1; i++) {
            for (int j = 0; j < order + 1; j++) {
                std::cout << "N_Transposed[" << i << "][" << j << "] = " << N_Transposed[i][j] << std::endl;
            }
        }

    }


    
    void compute_for_gradient() {
        for (unsigned int i = 0; i < q_order; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                        Ndx_Transposed[i][j] = poly.eval_1st_derivative(j, quad.knots[i]);
                        NWdx[j][i] = quad.weights[i] * Ndx_Transposed[i][j];
                    }
            }
        multiply_matrices < order + 1, q_order, q_order, order + 1, order + 1, order + 1, y_type > (NWdx, Ndx_Transposed, Ndx_Product);
    }

    void compute_for_laplacian() {
        for (unsigned int i = 0; i < q_order; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                        Ndxdx_Transposed[i][j] = poly.eval_2nd_derivative(j, quad.knots[i]);
                    }
            }
        multiply_matrices < order + 1, q_order, q_order, order + 1, order + 1, order + 1, y_type > (NW, Ndxdx_Transposed, Ndxdx_Product);
    }

    void vmult_mass(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) {

        if (c==0) {
                compute_basis_matrix();
            }

        // Initialize matrices for calculations
        std::array < std::array < y_type, order + 1 >, order + 1 > Sub;
        // Calculate (N_w * N^T) * U
        multiply_matrices<order+1,order+1,order+1,order+1,order+1,order+1,y_type>(N_Product,u,Sub);
        // Calculate ((N_w * N^T) * U) * (N_w * N^T)^T
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sub, N_Product, y, 1);

    }

    void vmult_gradient(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) {
        if (c==0) {
                compute_basis_matrix();
                compute_for_gradient();
            }
        if ((c==1)||(c==3)) {
                compute_for_gradient();
            }

        std::array < std::array < y_type, order + 1 >, order + 1 > Sub;

        std::array < std::array < y_type, order + 1 >, order + 1 > Sum_1;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Ndx_Product, u, Sub);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sub, N_Product, Sum_1, 1);

        std::array < std::array < y_type, order + 1 >, order + 1 > Sum_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (N_Product, u, Sub);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sub, Ndx_Product, Sum_2, 1);

        add_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sum_1, Sum_2, y);
    }

    constexpr void vmult_laplacian(std::array < std::array < y_type, order + 1 >, order + 1 > &y, std::array < std::array < y_type, order + 1 >, order + 1 > &u) {
        if (c==0) {
                compute_basis_matrix();
            }
        if (c<3) {
                compute_for_laplacian();
            }

        std::array < std::array < y_type, order + 1 >, order + 1 > Sub;

        std::array < std::array < y_type, order + 1 >, order + 1 > Sum_1;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Ndxdx_Product, u, Sub);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sub, N_Product, Sum_1, 1);

        std::array < std::array < y_type, order + 1 >, order + 1 > Sum_2;
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (N_Product, u, Sub);
        multiply_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sub, Ndxdx_Product, Sum_2, 1);

        add_matrices < order + 1, order + 1, order + 1, order + 1, order + 1, order + 1, y_type > (Sum_1, Sum_2, y);
    }*/

};

#endif

