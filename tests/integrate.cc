#include <iostream>
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>
#include "../include/la_operations.h"
#include "../include/quadrature.h"
#include "../include/quadrature_nonclass.h"
#include "../include/polynomial.h"
#include "../include/vmult.h"
#include "../include/integrate.h"

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
constexpr std::array<y_type, size> create_vector()
{
  std::array<y_type, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    arr[i] = i;
  return arr;
}

template <typename y_type, size_t size>
constexpr std::array<std::array<y_type, size>, size> create_array()
{
  std::array<std::array<y_type, size>, size> arr{1.};
  for (unsigned int i = 0; i < size; ++i)
    for (unsigned int j = 0; j < size; ++j)
      arr[i][j] = 1;
  return arr;
}


int main()
{
  // Calculate with old method
  constexpr unsigned int order = 10;
  constexpr unsigned int q_order = order+1;
  //Integrate<order, long double, Polynomial , Quadrature> lagrange;
  //std::array < std::array < long double, order + 1 >, order + 1 > u = create_array < long double, order + 1 > ();
  //std::array < std::array < long double, order + 1 >, order + 1 > y = create_array < long double, order + 1 > ();
  //unsigned int count = lagrange.vmult_mass(y, u);

  // Calculate with new method
  VMULT<order,q_order,5, long double, Polynomial , Quadrature> vmult;
  std::array < std::array < long double, order + 1 >, order + 1 > u_2 = create_array < long double, order + 1 > ();
  std::array < std::array < long double, order + 1 >, order + 1 > y_2 = create_array < long double, order + 1 > ();
  vmult.vmult_mass(y_2, u_2);

  // Output resulting vector y.
  for (unsigned int i = 0; i < order + 1; i++) {
    for (unsigned int j = 0; j < order + 1; j++) {
      std::cout << "y[" << i  << "," << j << "] = " << y_2[i][j] << endl;
    }
  }

  // Hardcode solution
  Quadrature<order, long double> quad;
  constexpr std::array < y_type, order + 1 > knots = quad.compute_quadrature_points();
  std::array < y_type, order + 1 > weights  = quad.compute_quadrature_weights(knots);
  array < array < long double, order + 1 >, order + 1 > y_hard = lagrange_nodes<long double, order>(u_2, weights);

  // Testing for correctness
  for (unsigned int i = 0; i < order + 1; i++) {
    for (unsigned int j = 0; j < order + 1; j++) {
      assert(y_2[i][j] == y_hard[i][j]);
    }
  }

  // Testing for correct complexity

  cout << "Testing was succesful." << endl;
  
  return 0;
}
