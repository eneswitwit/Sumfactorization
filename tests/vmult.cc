#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <vector>
#include <limits>
#include <iostream>

#include "../include/constexpr_math.h"
#include "../include/constexpr_array.h"
#include "../include/la_operations.h"
#include "../include/constexpr_quadrature.h"
#include "../include/polynomialbasis/constexpr_lagrange.h"
#include "../include/vmult.h"

using namespace std;

// Hardcode solution
template <typename y_type, size_t order>
array < array < y_type, order + 1 >, order + 1 > lagrange_nodes(array < array < y_type, order + 1 >, order + 1 > u, array < y_type, order + 1 > q_weights) {
  for (unsigned int i = 0; i < q_weights.size(); i++) {
    for (unsigned int j = 0; j < q_weights.size(); j++) {
      u[i][j] *= q_weights[i] * q_weights[j];
    }
  }
  return u;
}

template <typename y_type, size_t size>
constexpr_array<y_type, size> create_vector()
{
  constexpr_array<y_type, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    arr[i] = i;
  return arr;
}

template <typename y_type, size_t size>
constexpr_array<constexpr_array<y_type, size>, size> create_array()
{
  constexpr_array<constexpr_array<y_type, size>, size> arr;
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      arr[i][j] = 1;
  return arr;
}


int main()
{
  // Initialize VMULT
  constexpr size_t order = 3;
  constexpr size_t q_order = 3;
  VMULT<double, order,q_order,1, Quadrature , Lagrange > vmult;

  // Compute VMULT Mass Matrix
  constexpr_array < constexpr_array < long double, order + 1 >, order + 1 > u_2 = create_array < long double, order + 1 > ();
  constexpr_array < constexpr_array < long double, order + 1 >, order + 1 > y_2 = create_array < long double, order + 1 > ();
  //vmult.vmult_mass(y_2, u_2);

  // Output resulting vector y.
  for (unsigned int i = 0; i < order + 1; i++) {
    for (unsigned int j = 0; j < order + 1; j++) {
      std::cout << "y[" << i  << "," << j << "] = " << y_2[i][j] << endl;
    }
  }


  /** TESTING **/

  /*
  // Hardcode solution
  Quadrature<order, long double> quad;

  for (int i = 0; i < order + 1; i++) {
    std::cout << "quad.knots[" << i << "] = " << quad.knots[i] << std::endl;
  }

  array < array < long double, order + 1 >, order + 1 > y_hard = lagrange_nodes<long double, order>(u_2, quad.weights);

  // Testing for correctness
  for (unsigned int i = 0; i < order + 1; i++) {
    for (unsigned int j = 0; j < order + 1; j++) {
      assert(y_2[i][j] == y_hard[i][j]);
    }
  }

  // Testing for correct complexity

  cout << "Testing was succesful." << endl;
  */

  return 0;
}
