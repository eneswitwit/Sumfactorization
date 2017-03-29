#include <cassert>
#include <cmath>
#include <constexpr_matrix.h>

int main() {
    constexpr int n=2;
    constexpr_array <constexpr_array <long double , n >,n+1> B;
    constexpr_array <constexpr_array <long double , n+1 >,n> A;
    constexpr_array <constexpr_array <long double , n >,n> C;
    constexpr_array <constexpr_array <long double , n >,n> C_hard;

    for (int i=0;i<n;i++) {
        for (int j=0;j<n+1;j++) {
            A[i][j]=i+j;
            B[j][i]=i-2-j;
        }
    }

    C=multiply_matrices<n,n+1,n+1,n,long double>(A,B);


    C_hard[0][0]=-11;
    C_hard[0][1]=-8;
    C_hard[1][0]=-20;
    C_hard[1][1]=-14;

    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            assert(std::abs(C[i][j]-C_hard[i][j])<0.00000000000001);
        }
    }

    return 0;
}
