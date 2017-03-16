#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <type_traits>

#include "../include/constexpr_sin.h"
#include "../include/constexpr_array.h"
#include "../include/quadrature_constexpr.h"



constexpr static int order = 4;


int main() {
    constexpr Quadrature<double, order> quad;

    // Hardcode inputs for quadrature
    constexpr double val = math::sin(0.0);

    return 0;
}