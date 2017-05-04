#include <array>
#include <cassert>
#include <iostream>
#include <chrono>
#include "../include/constexpr_vmult_3D.h"
#include "../include/constexpr_quadrature.h"
#include "../include/polynomialbasis/constexpr_lagrange.h"

using namespace std;

template <typename Number, size_t size>
constexpr_array<Number, size> create_vector()
{
    constexpr_array<Number, size> arr{1.};
    for (unsigned int i = 0; i < size; ++i)
        arr[i] = i;
    return arr;
}

template <typename Number, size_t size>
constexpr_array<constexpr_array<constexpr_array<Number, size>, size>,size> create_array()
{
    constexpr_array<constexpr_array<constexpr_array<Number, size>, size>, size> arr;
    for (unsigned int i = 0; i < size; ++i)
        for (unsigned int j = 0; j < size; ++j)
            for (unsigned int k = 0; k < size; ++k)
                arr[i][j][k] = 10^i;
    return arr;
}


int main()
{
    //for (unsigned int i=3;i<10;i++) {
        constexpr size_t order = 16;
        constexpr size_t q_order = 16;

        constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> u_1 = create_array < long double, order + 1 > ();
        constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_mass_usememory = create_array < long double, order + 1 > ();
        constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_mass_nomemory = create_array < long double, order + 1 > ();

        std::chrono::steady_clock::time_point begin_usememory = std::chrono::steady_clock::now();
        constexpr VMULT<long double, order, q_order, Quadrature , Lagrange> vmult_usememory;
        vmult_usememory.mass(y_mass_usememory, u_1);
        vmult_usememory.gradient(y_mass_nomemory, u_1);
        vmult_usememory.laplacian(y_mass_nomemory, u_1);
        std::chrono::steady_clock::time_point end_usememory= std::chrono::steady_clock::now();

        std::chrono::steady_clock::time_point begin_nomemory = std::chrono::steady_clock::now();
        constexpr VMULT<long double, order, q_order, Quadrature , Lagrange, 0 > vmult_nomemory;
        vmult_nomemory.mass(y_mass_nomemory, u_1);
        vmult_nomemory.gradient(y_mass_nomemory, u_1);
        vmult_nomemory.laplacian(y_mass_nomemory, u_1);
        std::chrono::steady_clock::time_point end_nomemory= std::chrono::steady_clock::now();

        std::cout << " | " << order + 1 << " || " << std::chrono::duration_cast<std::chrono::milliseconds> (end_usememory - begin_usememory).count()
                  << " ms || " << std::chrono::duration_cast<std::chrono::milliseconds> (end_nomemory - begin_nomemory).count()<< " ms" << std::endl;
    //}

    return 0;
}
