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
 *
 * The boolean template parameter 'usememory' will determine whether we precompute the shape functions at the quadrature points and store them or compute the values from scratch when needed.
 * The default value is 1, this means we will precompute and store the values by default.
 */


template <typename Number, size_t order, size_t q_order, template<typename, size_t> class Quadrature, template<typename, size_t, template<typename, size_t> class Quadrature_ > class Polynomial, bool usememory=1>
class VMULT {
public:

    const Quadrature < Number, q_order > quad;
    const Polynomial <Number, order, Quadrature> poly;
    //Number mass_counter_;

    /**
     * The following matrices are the same as in 2D, since the matrices themselves correspond to only one dimension.
     */
    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NT_;
    const constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NW_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NP_;

    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXT_;
    const constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NWDX_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXP_;

    const constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXDXT_;
    const constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXDXP_;

    VMULT() :
        quad(),
        poly(),
        NT_(compute_basis_matrix_NT()),
        NW_(compute_basis_matrix_NW()),
        NP_(compute_basis_matrix_NP()),
        NDXT_(compute_gradient_NDXT()),
        NWDX_(compute_gradient_NWDX()),
        NDXP_(compute_gradient_NDXP()),
        NDXDXT_(compute_laplacian_NDXDXT()),
        NDXDXP_(compute_laplacian_NDXDXP()){}


    /** BASIS MATRICES **/

    /**
     * The documentation explains why the matrices are defined the way they are.
     */

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_basis_matrix_NT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NT;
        if (usememory==1) {
            for (unsigned int i = 0; i < q_order + 1 ; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                    NT[i][j] = poly.eval(quad.knots_[i], j);
                }
            }
        }
        return NT;

    }


    constexpr constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > compute_basis_matrix_NW() {
        constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NW;
        if (usememory==1) {
            for (unsigned int i = 0; i < q_order +1 ; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                    NW[j][i] = quad.weights_[i] * NT_[i][j];
                }
            }
        }
        return NW;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_basis_matrix_NP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NP;
        if (usememory==1) {
            NP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NW_, NT_);
        }
        return NP;
    }





    /** GRADIENT **/
    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_gradient_NDXT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXT;
        if (usememory==1) {
            for (unsigned int i = 0; i < q_order + 1 ; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                    NDXT[i][j] = poly.eval_1st_derivative(quad.knots_[i], j);
                }
            }
        }
        return NDXT;

    }


    constexpr constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > compute_gradient_NWDX() {
        constexpr_array < constexpr_array < Number, q_order + 1 >, order + 1 > NWDX;
        if (usememory==1) {
            for (unsigned int i = 0; i < q_order + 1; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                    NWDX[j][i] = quad.weights_[i] * NDXT_[i][j];
                }
            }
        }
        return NWDX;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_gradient_NDXP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXP;
        if (usememory==1) {
            NDXP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NWDX_, NDXT_);
        }
        return NDXP;
    }




    /** LAPLACIAN **/
    constexpr constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > compute_laplacian_NDXDXT() {
        constexpr_array < constexpr_array < Number, order + 1 >, q_order + 1 > NDXDXT;
        if (usememory==1) {
            for (unsigned int i = 0; i < q_order + 1; i++) {
                for (unsigned int j = 0; j < order + 1; j++) {
                    NDXDXT[i][j] = poly.eval_2nd_derivative(quad.knots_[i], j);
                }
            }
        }
        return NDXDXT;
    }

    constexpr constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > compute_laplacian_NDXDXP() {
        constexpr_array < constexpr_array < Number, order + 1 >, order + 1 > NDXDXP;
        if (usememory==1) {
            NDXDXP = multiply_matrices < order + 1, q_order + 1, q_order + 1, order + 1, Number > (NW_, NDXDXT_);
        }
        return NDXDXP;
    }

    /**
     * The function has two inputs: a refference to the solution y and a refference to the coefficient vector u.
     *
     * The explanation on the formulas can be seen in a document.
     *
     * In each loop (consisting of 4 for-loops) we have an additional if-statement.
     * It depends of the boolean template parameter 'usememory', and will either use the precomputed values or compute them from scratch.
     */

    constexpr void mass(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &y, constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &u) const {
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> W;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> M;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Z;
        //unsigned int mass_counter=0;

        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    W[k][i][j]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
          //                  mass_counter++;
                            W[k][i][j]+=u[i][j][l]*NP_[k][l];
                        }
                        else {
                            Number NP_lk = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                NP_lk+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],k);
            //                    mass_counter+=2;
                            }
                            W[k][i][j]+=u[i][j][l]*NP_lk;
              //              mass_counter++;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    M[j][k][i]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                //            mass_counter++;
                            M[j][k][i]+=W[k][i][l]*NP_[j][l];
                        }
                        else {
                            Number NP_lj = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                NP_lj+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],j);
                  //              mass_counter+=2;
                            }
                    //        mass_counter++;
                            M[j][k][i]+=W[k][i][l]*NP_lj;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    Z[i][j][k]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                      //      mass_counter++;
                            Z[i][j][k]+=M[j][k][l]*NP_[i][l];
                        }
                        else {
                            Number NP_li = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                        //        mass_counter+=2;
                                NP_li+=quad.weights_[q]*poly.eval(quad.knots_[q],i)*poly.eval(quad.knots_[q],l);
                            }
                          //  mass_counter++;
                            Z[i][j][k]+=M[j][k][l]*NP_li;
                        }
                    }
                }
            }
        }
        y=Z;
    }

    constexpr void gradient(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &y, constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &u) const {

        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> W;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> M;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Z;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Wd;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Md1;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Md2;

        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    W[k][i][j]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                            W[k][i][j]+=u[i][j][l]*(NP_[k][l]);
                            Wd[k][i][j]+=u[i][j][l]*(NDXP_[k][l]);
                        }
                        else {
                            Number N_lk = 0;
                            Number Ndx_lk = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_lk+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],k);
                                Ndx_lk+=quad.weights_[q]*poly.eval_1st_derivative(quad.knots_[q],l)*poly.eval_1st_derivative(quad.knots_[q],k);
                            }
                            W[k][i][j]+=u[i][j][l]*N_lk;
                            Wd[k][i][j]+=u[i][j][l]*Ndx_lk;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    M[j][k][i]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                            M[j][k][i]+=W[k][i][l]*(NP_[j][l]);
                            Md1[j][k][i]+=W[k][i][l]*(NDXP_[j][l]);
                            Md2[j][k][i]+=Wd[k][i][l]*(NP_[j][l]);
                        }
                        else {
                            Number N_lj = 0;
                            Number Ndx_lj = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_lj+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],j);
                                Ndx_lj+=quad.weights_[q]*poly.eval_1st_derivative(quad.knots_[q],l)*poly.eval_1st_derivative(quad.knots_[q],j);
                            }
                            M[j][k][i]+=W[k][i][l]*N_lj;
                            Md1[j][k][i]+=W[k][i][l]*Ndx_lj;
                            Md2[j][k][i]+=Wd[k][i][l]*N_lj;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    Z[i][j][k]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1)
                            Z[i][j][k]+=(M[j][k][l]*NDXP_[i][l]+Md1[j][k][l]*NP_[i][l]+Md2[j][k][l]*NP_[i][l]);
                        else {
                            Number N_li = 0;
                            Number Ndx_li = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_li+=quad.weights_[q]*poly.eval(quad.knots_[q],i)*poly.eval(quad.knots_[q],l);
                                Ndx_li+=quad.weights_[q]*poly.eval_1st_derivative(quad.knots_[q],l)*poly.eval_1st_derivative(quad.knots_[q],i);
                            }
                            Z[i][j][k]+=(M[j][k][l]*Ndx_li + Md1[j][k][l]*N_li + Md2[j][k][l]*N_li);
                        }
                    }
                }
            }
        }
        y=Z;
    }


    constexpr void laplacian(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &y, constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> &u) const {

        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> W;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> M;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Z;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Wd;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Md1;
        constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> Md2;

        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    W[k][i][j]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                            Wd[k][i][j]+=u[i][j][l]*NDXDXP_[k][l];
                            W[k][i][j]+=u[i][j][l]*NP_[k][l];
                        }
                        else {
                            Number N_lk = 0;
                            Number Ndxdx_lk = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_lk+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],k);
                                Ndxdx_lk+=quad.weights_[q]*poly.eval_2nd_derivative(quad.knots_[q],l)*poly.eval(quad.knots_[q],k);
                            }
                            W[k][i][j]+=u[i][j][l]*N_lk;
                            Wd[k][i][j]+=u[i][j][l]*Ndxdx_lk;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    M[j][k][i]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1) {
                            M[j][k][i]+=W[k][i][l]*NP_[j][l];
                            Md1[j][k][i]+=W[k][i][l]*NDXDXP_[j][l];
                            Md2[j][k][i]+=Wd[k][i][l]*NP_[j][l];
                        }
                        else {
                            Number N_lj = 0;
                            Number Ndxdx_lj = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_lj+=quad.weights_[q]*poly.eval(quad.knots_[q],l)*poly.eval(quad.knots_[q],j);
                                Ndxdx_lj+=quad.weights_[q]*poly.eval_2nd_derivative(quad.knots_[q],l)*poly.eval(quad.knots_[q],j);
                            }
                            M[j][k][i]+=W[k][i][l]*N_lj;
                            Md1[j][k][i]+=W[k][i][l]*Ndxdx_lj;
                            Md2[j][k][i]+=Wd[k][i][l]*N_lj;
                        }
                    }
                }
            }
        }
        for (unsigned int k=0; k < order+1; k++) {
            for (unsigned int i=0 ;i < order+1; i++) {
                for (unsigned int j=0 ;j < order+1; j++) {
                    Z[i][j][k]=0;
                    for (unsigned int l=0;l < order+1; l++) {
                        if (usememory==1)
                            Z[i][j][k]-=(M[j][k][l]*NDXDXP_[i][l] + Md1[j][k][l]*NP_[i][l] + Md2[j][k][l]*NP_[i][l]);
                        else {
                            Number N_li = 0;
                            Number Ndxdx_li = 0;
                            for (unsigned int q=0; q<q_order+1;q++) {
                                N_li+=quad.weights_[q]*poly.eval(quad.knots_[q],i)*poly.eval(quad.knots_[q],l);
                                Ndxdx_li+=quad.weights_[q]*poly.eval_2nd_derivative(quad.knots_[q],l)*poly.eval(quad.knots_[q],i);
                            }
                            Z[i][j][k]-=(M[j][k][l]*Ndxdx_li + Md1[j][k][l]*N_li + Md2[j][k][l]*N_li);
                        }
                    }
                }
            }
        }
        y=Z;
    }

};

#endif

