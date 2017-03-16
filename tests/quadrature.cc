#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <cstdint>
#include <array>
#include <type_traits>
#include <iostream>

#include "../include/constexpr_sin.h"
#include "../include/constexpr_array.h"
#include "../include/quadrature_constexpr.h"



constexpr static int order = 4;


int main() {
	constexpr Quadrature<double, order> quad;
	for (int i = 0; i < order + 1; i++)
	{
		std::cout << "quad[" << i << "] = " << quad[i] << std::endl;
	}

	return 0;
}