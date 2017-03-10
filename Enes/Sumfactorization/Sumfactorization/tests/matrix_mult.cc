#include <iostream>
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include "../include/la_operations.h"

template <typename y_type, size_t size>
constexpr std::array<y_type, size> create_vector()
{
  std::array<y_type, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    arr[i] = std::rand() & 5;
  return arr;
}

template <typename y_type, size_t size>
constexpr std::array<std::array<y_type, size>, size> create_array()
{
  std::array<std::array<y_type, size>, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      arr[i][j] = std::rand() % 5;
  return arr;
}

int main() {
    constexpr int order=2;
    std::array < std::array < long double, order + 1 >, order + 1 > A = create_array < long double, order + 1 > ();
    std::array < std::array < long double, order + 1 >, order + 1 > B = create_array < long double, order + 1 > ();
    std::array < std::array < long double, order + 1 >, order + 1 > C = create_array < long double, order + 1 > ();
    multiply_matrices<order+1,order+1,order+1,order+1,order+1,order+1,long double>(A,B,C,1);

    for (int i=0;i<3;i++) {
        std::cout << A[i][0] << "  " << A[i][1] << "  " << A[i][2] << "             " << B[i][0] << "  " << B[i][1] << "  " << B[i][2] << std::endl;
    }
    std::cout << std::endl;

    for (int i=0;i<3;i++) {
        std::cout << C[i][0] << "  " << C[i][1] << "  " << C[i][2] << std::endl;
    }
    std::cout << std::endl;

}
