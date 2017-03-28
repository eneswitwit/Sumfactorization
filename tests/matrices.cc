#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <limits>
#include <type_traits>
#include <iostream>

#include "../include/constexpr_math.h"
#include "../include/constexpr_array.h"
#include "../include/constexpr_matrix.h"

int main() {
    constexpr int n=2;
    constexpr_array <constexpr_array <long double , n+1 >,n> B;
    constexpr_array <constexpr_array <long double , n >,n+1> A;
    constexpr_array <constexpr_array <long double , n >,n> C;

    for (int i=0;i<n;i++) {
        for (int j=0;j<n+3;j++) {
            A[i][j]=i+j;
            B[j][i]=i-2-j;
        }
    }

    C=multiply_matrices<n,n+1,n+1,n,long double>(A,B);

    for (int i=0;i<n;i++) {
        for (int j=0;j<n+1;j++) {
            std::cout << A[i][j] << "  ";
        }
        std::cout << std::endl;
    } std::cout << std::endl;

    for (int i=0;i<n+1;i++) {
        for (int j=0;j<n;j++) {
            std::cout << B[i][j] << "  ";
        }
        std::cout << std::endl;
    } std::cout << std::endl;

    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            std::cout << C[i][j] << "  ";
        }
        std::cout << std::endl;
    } std::cout << std::endl;

    return 0;
}
