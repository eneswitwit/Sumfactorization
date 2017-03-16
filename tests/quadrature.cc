#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <type_traits>

#include "../include/constexpr_sin.h"
#include "../include/constexpr_array.h"
#include "../include/quadrature_constexpr.h"

template<typename y_type, size_t order>
constexpr constexpr_array< y_type , order> quadrature_sin_vector() {
    constexpr_array < y_type, order> quad_sin_vec;
    int m = order-1;
    for (unsigned int i = 0; i < order; i++) {
        quad_sin_vec[i] = -math::sin( sin_input<y_type>(i, m) ) ;
    }
    return quad_sin_vec;
};

constexpr static int order = 4;
constexpr static constexpr_array< double , order> quad_sin_vec = quadrature_sin_vector<double,order>();


int main() {
    constexpr Quadrature<double, order> quad(quad_sin_vec);

    // Hardcode inputs for quadrature
    constexpr double val = math::sin(0.0);

    return 0;
}