#include <array>
#include <cassert>
#include <iostream>
#include "../include/constexpr_vmult_3D.h"
#include "../include/constexpr_quadrature.h"
#include "../include/polynomialbasis/constexpr_lagrange.h"

using namespace std;

// Hardcode solution
template <typename Number, size_t order, size_t q_order>
constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> test_mass(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> u) {

    Lagrange<long double,order,Quadrature> poly;
    Quadrature<long double,q_order> quad;

    constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> y;

    for (unsigned int l = 0; l < order+1; l++) {
        for (unsigned int m = 0; m < order+1; m++) {
            for (unsigned int n = 0; n < order+1; n++) {
                y[l][m][n]=0;
                for (unsigned int i = 0; i < order+1; i++) {
                    for (unsigned int j = 0; j < order+1; j++) {
                        for (unsigned int k = 0; k < order+1; k++) {

                            Number sumx=0;
                            for (unsigned int qx=0;qx<q_order+1;qx++) {
                                sumx+=quad.weights_[qx]*poly.eval_lagrange(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l);
                            }

                            Number sumy=0;
                            for (unsigned int qy=0;qy<q_order+1;qy++) {
                                sumy+=quad.weights_[qy]*poly.eval_lagrange(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m);
                            }

                            Number sumz=0;
                            for (unsigned int qz=0;qz<q_order+1;qz++) {
                                sumz+=quad.weights_[qz]*poly.eval_lagrange(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n);
                            }

                            y[l][m][n]+=u[i][j][k]*sumy*sumx*sumz;
                        }
                    }
                }

            }
        }
    }
    return y;
}

template <typename Number, size_t order, size_t q_order>
constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order +1> test_gradient(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order +1> u) {

    Lagrange<long double,order,Quadrature> poly;
    Quadrature<long double,q_order> quad;

    constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order +1> y;

    for (unsigned int l = 0; l < order+1; l++) {
        for (unsigned int m = 0; m < order+1; m++) {
            for (unsigned int n = 0; n < order+1; n++) {
                y[l][m][n]=0;
                for (unsigned int i = 0; i < order+1; i++) {
                    for (unsigned int j = 0; j < order+1; j++) {
                        for (unsigned int k = 0; k < order+1; k++) {
                            Number sumx=0;
                            for (unsigned int qx=0;qx<q_order+1;qx++) {
                                sumx+=quad.weights_[qx]*(poly.eval_1st_derivative(quad.knots_[qx],i)*poly.eval_1st_derivative(quad.knots_[qx],l)
                                                         +poly.eval_lagrange(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l)
                                                         +poly.eval_lagrange(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l));
                            }

                            Number sumy=0;
                            for (unsigned int qy=0;qy<q_order+1;qy++) {
                                sumy+=quad.weights_[qy]*(poly.eval_1st_derivative(quad.knots_[qy],j)*poly.eval_1st_derivative(quad.knots_[qy],m)
                                                         +poly.eval_lagrange(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m)
                                                         +poly.eval_lagrange(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m));
                            }

                            Number sumz=0;
                            for (unsigned int qz=0;qz<q_order+1;qz++) {
                                sumz+=quad.weights_[qz]*(poly.eval_1st_derivative(quad.knots_[qz],k)*poly.eval_1st_derivative(quad.knots_[qz],n)
                                                         +poly.eval_lagrange(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n)
                                                         +poly.eval_lagrange(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n));
                            }

                            y[l][m][n]+=u[i][j][k]*sumy*sumx*sumz;
                        }
                    }
                }
            }
        }
    }
    return y;
}

template <typename Number, size_t order, size_t q_order>
constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> test_laplace(constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order + 1> u) {

    Lagrange<long double,order,Quadrature> poly;
    Quadrature<long double,q_order> quad;

    constexpr_array < constexpr_array < constexpr_array < Number, order + 1 >, order + 1 >, order +1> y;

    for (unsigned int l = 0; l < order+1; l++) {
        for (unsigned int m = 0; m < order+1; m++) {
            for (unsigned int n = 0; n < order+1; n++) {
                y[l][m][n]=0;
                for (unsigned int i = 0; i < order+1; i++) {
                    for (unsigned int j = 0; j < order+1; j++) {
                        for (unsigned int k = 0; k < order+1; k++) {
                            Number sumx=0;
                            for (unsigned int qx=0;qx<q_order+1;qx++) {
                                sumx+=quad.weights_[qx]*(+poly.eval_2nd_derivative(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l)
                                                         +poly.eval_lagrange(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l)
                                                         +poly.eval_lagrange(quad.knots_[qx],i)*poly.eval_lagrange(quad.knots_[qx],l));
                            }

                            Number sumy=0;
                            for (unsigned int qy=0;qy<q_order+1;qy++) {
                                sumy+=quad.weights_[qy]*(+poly.eval_2nd_derivative(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m)
                                                         +poly.eval_lagrange(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m)
                                                         +poly.eval_lagrange(quad.knots_[qy],j)*poly.eval_lagrange(quad.knots_[qy],m));
                            }

                            Number sumz=0;
                            for (unsigned int qz=0;qz<q_order+1;qz++) {
                                sumz+=quad.weights_[qz]*(+poly.eval_2nd_derivative(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n)
                                                         +poly.eval_lagrange(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n)
                                                         +poly.eval_lagrange(quad.knots_[qz],k)*poly.eval_lagrange(quad.knots_[qz],n));
                            }

                            y[l][m][n]-=u[i][j][k]*sumy*sumx*sumz;
                        }
                    }
                }
            }
        }
    }
    return y;
}

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
    // Initialize VMULT
    constexpr size_t order = 3;
    constexpr size_t q_order = 3;
    constexpr VMULT<long double, order, q_order, Quadrature , Lagrange > vmult;

    // Compute VMULT Mass Matrix
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> u_1 = create_array < long double, order + 1 > ();
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_mass = create_array < long double, order + 1 > ();
    vmult.mass(y_mass, u_1);

    // Compute VMULT Gradient
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_gradient = create_array < long double, order + 1 > ();
    vmult.gradient(y_gradient, u_1);

    // Compute VMULT Laplacian
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_laplace = create_array < long double, order + 1 > ();
    vmult.laplacian(y_laplace, u_1);




    /** TESTING **/
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_mass_hard = test_mass<long double, order,q_order>(u_1);
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_gradient_hard = test_gradient<long double, order,q_order>(u_1);
    constexpr_array < constexpr_array < constexpr_array < long double, order + 1 >, order + 1 >, order + 1> y_laplace_hard = test_laplace<long double, order,q_order>(u_1);

    const long double eps=0.00000000001;

    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            //std::cout << y_mass[i][j] << "     " <<  y_mass_hard[i][j] << std::endl;
        }
    }

    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            //std::cout << y_gradient[i][j] << "     " <<  y_gradient_hard[i][j] << std::endl;
        }
    }

    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            //std::cout << y_laplace[i][j] << "     " <<  y_laplace_hard[i][j] << std::endl;
        }
    }

    // Test vmult.mass
    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            for (unsigned int k = 0; k < order + 1; k++) {
                assert(abs(y_mass[i][j][k]-y_mass_hard[i][j][k])<=eps);
            }
        }
    }
    std::cout << "VMULT_3D.mass     test successful" << std::endl;

    // Test vmult.gradient
    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            for (unsigned int k = 0; k < order + 1; k++) {
                assert(abs(y_gradient[i][j][k]-y_gradient_hard[i][j][k])<=eps);
            }
        }
    }
    std::cout << "VMULT_3D.gradient test successful" << std::endl;

    // Test vmult.laplace
    for (unsigned int i = 0; i < order + 1; i++) {
        for (unsigned int j = 0; j < order + 1; j++) {
            for (unsigned int k = 0; k < order + 1; k++) {
                assert(abs(y_laplace[i][j][k]-y_laplace_hard[i][j][k])<=eps);
            }
        }
    }
    std::cout << "VMULT_3D.laplace  test successful" << std::endl;

    return 0;
}
