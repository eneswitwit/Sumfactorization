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
#include "../include/constexpr_quadrature.h"



constexpr static int order = 3;


int main() {
        //std::vector<long double> vec = compute_quadrature_points(order+1,1,1);
	for (int i = 0; i < order+1 ; i++)
	{
        //	std::cout << "vec[" << i << "] = " << vec[i] << std::endl;
	}

	constexpr Quadrature<long double, order> quad;
        for (int i = 0; i < order+1; i++)
	{
                std::cout << "quad[" << i << "] = " << quad.knots_[i] << "              " << quad.weights_[i] << std::endl;
	}

	return 0;
}
