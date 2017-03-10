#ifndef __LA_OPERATIONS_H__
#define __LA_OPERATIONS_H__

template <int A_rows, int A_columns,int B_rows, int B_columns,int C_rows, int C_columns,typename y_type>
void multiply_matrices(std::array<std::array<y_type, A_rows>, A_columns> &A,std::array<std::array<y_type, B_rows>, B_columns> &B,std::array<std::array<y_type, C_rows>, C_columns> &C, bool transpose_B=0) {
    // Check if matrix dimensions match for multiplication
    static_assert(A_columns == B_rows, "Matrices not compatible for multiplication!");

    // Check if reference to solution has the correct dimension
    static_assert((A_rows==C_rows)&&(B_columns==C_columns), "Matrices not compatible for multiplication!");

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

}

template <int A_rows, int A_columns,int B_rows, int B_columns,int C_rows, int C_columns,typename y_type>
void add_matrices(std::array<std::array<y_type, A_rows>, A_columns> &A,std::array<std::array<y_type, B_rows>, B_columns> &B,std::array<std::array<y_type, C_rows>, C_columns> &C) {
    // Check if matrix dimensions match for addition
    static_assert ((A_columns==B_columns)&&(A_rows==B_rows), "Matrices not compatible for addition");

    // Check if reference to solution has the correct dimension
    static_assert ((A_columns==C_columns)&&(A_rows==C_rows), "Dimension of solution doesn't match");

        for (unsigned int i=0;i<B_columns;i++) {
            for (unsigned int j=0;j<A_rows;j++) {
                    C[i][j]=A[i][j]+B[i][j];
            }
        }

}

#endif
