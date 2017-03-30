#ifndef __VMULT_H__
#define __VMULT_H__

#include "constexpr_matrix.h"
#include "constexpr_quadrature.h"

/**
 * The methods used in this class are further explained in the document "Sumfactorization". The general principle is to exploit the tensor-product structure.
 * Our implementation is directly derived from the teachlet "Efficient evaluation of weak forms in discontinuous Galerkin methods"
 *
 * The template parameter 'order' and 'q_order' correspond to the order of the polynomials and the order of the quadrature rule.
 * We want to utilize two classes for evaluating polynomials and quadrature rules.
 */


template <typename Number, size_t order, size_t q_order, template<typename, size_t> class Quadrature, template<typename, size_t, template<typename, size_t> class Quadrature_ > class Polynomial>
class VMULT {
public:

    const Quadrature < Number, q_order > quad;
    const Polynomial <Number, order, Quadrature> poly;

    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NT_;
    const constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NW_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NP_;

    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXT_;
    const constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NWDX_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXP_;

    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXDXT_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXDXP_;

    constexpr VMULT() :
        quad(),
        poly(),
        NT_(compute_basis_matrix_NT()),
        NW_(compute_basis_matrix_NW()),
        NP_(compute_basis_matrix_NP()),
        NDXT_(compute_gradient_NDXT()),
        NWDX_(compute_gradient_NWDX()),
        NDXP_(compute_gradient_NDXP()),
        NDXDXT_(compute_laplacian_NDXDXT()),
        NDXDXP_(compute_laplacian_NDXDXP()) {}


    /** BASIS MATRICES **/

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_basis_matrix_NT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NT;
        for (unsigned int i = 0; i < q_order + 1 ; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NT[i][j] = poly.eval_lagrange(quad.knots_[i], j);
            }
        }
        return NT;

    }


    constexpr constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > compute_basis_matrix_NW() {
        constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NW;
        for (unsigned int i = 0; i < q_order +1 ; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NW[j][i] = quad.weights_[i] * NT_[i][j];
            }
        }

        return NW;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_basis_matrix_NP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NP;
        NP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NW_, NT_);
        return NP;
    }





    /** GRADIENT **/
    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_gradient_NDXT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXT;
        for (unsigned int i = 0; i < q_order + 1 ; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NDXT[i][j] = poly.eval_1st_derivative(quad.knots_[i], j);
            }
        }
        return NDXT;

    }


    constexpr constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > compute_gradient_NWDX() {
        constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NWDX;
        for (unsigned int i = 0; i < q_order + 1; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NWDX[j][i] = quad.weights_[i] * NDXT_[i][j];
            }
        }

        return NWDX;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_gradient_NDXP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXP;
        NDXP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NWDX_, NDXT_);
        return NDXP;
    }




    /** LAPLACIAN **/
    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_laplacian_NDXDXT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXDXT;
        for (unsigned int i = 0; i < q_order + 1; i++) {
            for (unsigned int j = 0; j < order + 1; j++) {
                NDXDXT[i][j] = -1. * poly.eval_2nd_derivative(quad.knots_[i], j);
            }
        }
        return NDXDXT;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_laplacian_NDXDXP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXDXP;
        NDXDXP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NW_, NDXDXT_);
        return NDXDXP;
    }


    constexpr void mass(constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &y, constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &u) const {
        // Initialize matrices for calculations
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sub;
        // Calculate (N_w * N^T) * U
        Sub = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number > (NP_, u);
        // Calculate ((N_w * N^T) * U) * (N_w * N^T)^T
        y = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number , 1 > (Sub, NP_);

    }

    constexpr void gradient(constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &y, constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &u) const {

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sub;

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sum_1;
        Sub = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number > (NDXP_, u);
        Sum_1 = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number , 1 > (Sub, NP_);

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sum_2;
        Sub = multiply_matrices < order + 1, order + 1, order + 1, order + 1,  Number > (NP_, u);
        Sum_2 = multiply_matrices < order + 1, order + 1, order + 1, order + 1,  Number > (Sub, NDXP_);

        y = add_matrices < order + 1, order + 1, order + 1, order + 1,  Number > (Sum_1, Sum_2);
    }


    constexpr void laplacian(constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &y, constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > &u) const {

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sub;

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sum_1;
        Sub = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number > (NDXDXP_, u);
        Sum_1 = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number, 1 > (Sub, NP_);

        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > Sum_2;
        Sub = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number > (NP_, u);
        Sum_2 = multiply_matrices < order + 1, order + 1, order + 1, order + 1, Number , 1 > (Sub, NDXDXP_);

        y = add_matrices < order + 1, order + 1, order + 1, order + 1 , Number > (Sum_1, Sum_2);
    }

};

#endif

