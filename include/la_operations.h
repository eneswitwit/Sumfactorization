#ifndef __LA_OPERATIONS_H__
#define __LA_OPERATIONS_H__

template <int A_rows, int A_columns,int B_rows, int B_columns,int C_rows, int C_columns,typename y_type>
void multiply_matrices(std::array<std::array<y_type, A_rows>, A_columns> &A,std::array<std::array<y_type, B_rows>, B_columns> &B,std::array<std::array<y_type, C_rows>, C_columns> &C, bool transpose_B=0) {
    // Check if matrix dimensions match for multiplication
    if (A_columns!=B_rows) {
        std::cout << "Matrices not compatible for multiplication!" << std::endl;
        return;
    }

    // Check if reference to solution has the correct dimension
    if ((A_rows!=C_rows)||(B_columns!=C_columns)) {
        std::cout << "Dimension of solution and reference to solution does not match!";
        return;
    }

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

template<int A_rows, int A_columns, int v_entries, int u_entries, typename y_type>
constexpr void multiply_vector(std::array<std::array<y_type, A_rows>, A_columns> &A, std::array<y_type, v_entries> &v, std::array<y_type, u_entries> &u) {
    if (A_columns!=v_entries) {
        std::cout << "Dimension mismatch, can't multiply vector with matrix" << std::endl;
        return;
    }

    if (A_rows!=u_entries) {
        std::cout << "Dimension of solution and reference to solution does not match!" << std::endl;
        return;
    }
    for (unsigned int i=0;i<A_rows;i++) {
        u[i]=0;
        for (unsigned int j=0;j<A_columns;j++) {
            u[j]+=A[i][j]*v[j];
        }
    }
}

#endif
