#ifndef __LA_OPERATIONS_H__
#define __LA_OPERATIONS_H__

#include "constexpr_array.h"

/**
  * This function implements a basic matrix multiplication with the added option to transpose the right hand side.
  * This way we don't need to actually transpose a matrix and compute the desired matrix directly.
  */
template <int A_rows, int A_columns,int B_rows, int B_columns,typename Number, bool transpose_B=0>
constexpr constexpr_array<constexpr_array<Number, B_columns>, A_rows> multiply_matrices(const constexpr_array<constexpr_array<Number, A_columns>, A_rows> &A,const constexpr_array<constexpr_array<Number, B_columns>, B_rows> &B) {

    constexpr_array<constexpr_array<Number, B_columns>, A_rows> C;

    if (transpose_B==0) {
        for (unsigned int i=0;i<B_columns;i++) {
            for (unsigned int j=0;j<A_rows;j++) {
                C[j][i]=0;
                for (unsigned int k=0;k<B_rows;k++) {
                    C[j][i]+=B[k][i]*A[j][k];
                }
            }
        }
    }
    else {
        for (unsigned int i=0;i<B_rows;i++) {
            for (unsigned int j=0;j<A_rows;j++) {
                C[j][i]=0;
                for (unsigned int k=0;k<B_columns;k++) {
                    C[j][i]+=B[i][k]*A[j][k];
                }
            }
        }
    }
    return C;

}
/**
 * Basic matrix addition
 * */
template <int A_rows, int A_columns,int B_rows, int B_columns,typename Number>
constexpr constexpr_array<constexpr_array<Number, A_rows>, A_columns> add_matrices(const constexpr_array<constexpr_array<Number, A_rows>, A_columns> &A,const constexpr_array<constexpr_array<Number, B_rows>, B_columns> &B) {
    // Check if matrix dimensions match for addition
    static_assert ((A_columns==B_columns)&&(A_rows==B_rows), "Matrices not compatible for addition");

    constexpr_array<constexpr_array<Number, A_rows>, A_columns> C;

        for (unsigned int i=0;i<A_columns;i++) {
            for (unsigned int j=0;j<A_rows;j++) {
                    C[i][j]=A[i][j]+B[i][j];
            }
        }

        return C;

}

#endif
