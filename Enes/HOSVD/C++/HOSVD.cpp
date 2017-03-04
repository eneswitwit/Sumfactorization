#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <limits>
#include <cmath>

#include "Quadrature.h"
#include "Polynomial.h"

// Tensor Library
#include "Tensor/tensors/tensor1.h"
#include "Tensor/tensors/tensor2.h"
#include "Tensor/tensors/tensor3.h"
#include "Tensor/tensors/tensor4.h"
#include "Tensor/tensors/operators.h"

// Eigen Library for using SVD
#include "Eigen/SVD"

using namespace std;
using namespace TensorsLib;
using namespace Eigen;

// Lagrange Basis Evaluation for specific x
double eval_lagr(int i, double x, array<double, 2> nodes )
{
    double val = 1.;
    for (int j = 0; j <= 1; j++) {
        if (i != j) {
            val *= (x - nodes[j]) / (nodes[i] - nodes[j]);
        }
    }
    return val;
};

// Compute Entries of Mass Matrix Tensor
template<typename y_type, int order>
Tensor4 < y_type , order + 1 > computeTensor(array < y_type, order + 1 > nodes, array < y_type, order + 1 > weights, array < y_type, order + 1 > vertex) {
    Tensor4 < y_type , order + 1 > M;
    // Calculate Tensor Entries
    for (int i1 = 0 ; i1 < order + 1 ; ++i1 ) {
        for (int i2 = 0 ; i2 < order + 1 ; ++i2 ) {
            for (int j1 = 0 ; j1 < order + 1 ; ++j1 ) {
                for (int j2 = 0 ; j2 < order + 1 ; ++j2) {
                    double sum = 0;
                    for (int k = 0; k < order + 1 ; ++k) {
                        double sum_inner = 0;
                        for (int l = 0; l < order + 1; ++l) {
                            sum_inner += eval_lagr(i2, nodes[l], vertex) * eval_lagr(j2, nodes[l], vertex) * weights[l];
                        }
                        sum += sum_inner * (weights[k] * eval_lagr(i1, nodes[k], vertex) * eval_lagr(j1, nodes[k], vertex));
                    }
                    M[i1][i2][j1][j2] = sum;
                }
            }
        }
    }
    return M;
};

// Print Tensor
template<typename y_type, int order>
void print_Tensor4(Tensor4 < y_type, order + 1 > & TENSOR) {
    cout << "\n" << "TENSOR[:][:][0][0]" << endl;
    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            cout << TENSOR[i][j][0][0] << " " ;
        }
        cout << endl;
    }
    cout << "\n" << "TENSOR[:][:][1][0]" << endl;
    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            cout << TENSOR[i][j][1][0] << " " ;
        }
        cout << endl;
    }
    cout << "\n" << "TENSOR[:][:][0][1]" << endl;
    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            cout << TENSOR[i][j][0][1] << " " ;
        }
        cout << endl;
    }
    cout << "\n" << "TENSOR[:][:][1][1]" << endl;
    for (int i = 0; i < order + 1; ++i) {
        for (int j = 0; j < order + 1; ++j) {
            cout << TENSOR[i][j][1][1] << " " ;
        }
        cout << endl;
    }
};



// Transform number of 2D function to numbers of the corresponding 1D functions
std::array<int, 2> transform(int k, int order) {
    int n = order + 1;
    // Find row
    int row = k / n;
    // Find Column
    int column = k - row * n;

    std::array<int, 2> index = {row, column};
    return index;
};

// Transform number of 1D functions to numbers of the corresponding 2D function
int transform(int row, int column, int order) {
    int k = row * (order + 1) + column;
    return k;
};

// Transform number of 4 1D functions to number number of the 2 corresponding 2D functions.
std::array<int, 2> transform(int i1, int i2, int j1, int j2, int order) {
    int k1 = transform(i1, i2, order);
    int k2 = transform(j1, j2, order);
    std::array<int, 2> k = {k1, k2};
    return k;
};



// Transform Mass Tensor to Mass Matrix
template<typename y_type, int order>
Tensor2 < y_type, (order + 1)*(order + 1) > tensorToMatrix(Tensor4 < y_type, order + 1 > & M) {
    Tensor2 < y_type, (order + 1)*(order + 1) > MM;
    for (int i1 = 0; i1 < order + 1; ++i1) {
        for (int i2 = 0; i2 < order + 1; ++i2) {
            for (int j1 = 0; j1 < order + 1; ++j1) {
                for (int j2 = 0; j2 < order + 1; ++j2) {
                    std::array<int, 2> k = transform(i1, i2, j1, j2, order);
                    MM[k[0]][k[1]] = M[i1][i2][j1][j2];
                }
            }
        }
    }
    return MM;
};

// Calculate J_k (Kolda and Bader Page 460)
int J_(int k, int n) {
    int Jk = 1;
    for (int m = 1; m < k; ++m) {
        if (m != n) {
            Jk *= 2;
        }
    }
    return Jk;
}

template<typename y_type, int order>
void mode_(int n, Tensor4 < y_type, order + 1 > & TENSOR , array<array<y_type, 2>, 8> & MATRIX) {
    for (int i1 = 1; i1 < order + 2 ; ++i1) {
        for (int i2 = 1; i2 < order + 2 ; ++i2) {
            for (int i3 = 1; i3 < order + 2 ; ++i3) {
                for (int i4 = 1; i4 < order + 2 ; ++i4) {
                    // Calculate j (Kolda and Bader Page 460)
                    int j = 1;
                    if (n != 1) {
                        j += (i1 - 1) * J_(1, n);
                    }
                    if (n != 2) {
                        j += (i2 - 1) * J_(2, n);
                    }
                    if (n != 3) {
                        j += (i3 - 1) * J_(3, n);
                    }
                    if (n != 4) {
                        j += (i4 - 1) * J_(4, n);
                    }

                    // Map Tensor entries to Matrix elements
                    if (n == 1) {
                        MATRIX[i1 - 1][j - 1] = TENSOR[i1 - 1][i2 - 1][i3 - 1][i4 - 1];
                    }
                    if (n == 2) {
                        MATRIX[i2 - 1][j - 1] = TENSOR[i1 - 1][i2 - 1][i3 - 1][i4 - 1];
                    }
                    if (n == 3) {
                        MATRIX[i3 - 1][j - 1] = TENSOR[i1 - 1][i2 - 1][i3 - 1][i4 - 1];
                    }
                    if (n == 4) {
                        MATRIX[i4 - 1][j - 1] = TENSOR[i1 - 1][i2 - 1][i3 - 1][i4 - 1];
                    }


                }
            }
        }

    }

}

template<typename y_type, int order, int col>
Tensor4 < y_type, order + 1 > n_mode_product(int n, Tensor4 < y_type, order + 1 > & TENSOR , array < array < y_type, order + 1 >, col > & MATRIX) {
    Tensor4 < y_type, order + 1 > TENSOR_RESULT;
    if (n == 1) {
        for (int j = 0; j < order + 1 ; ++j) {
            for (int i2 = 0; i2 < order + 1 ; ++i2) {
                for (int i3 = 0; i3 < order + 1 ; ++i3) {
                    for (int i4 = 0; i4 < order + 1 ; ++i4) {
                        double sum = 0;
                        for (int in = 0; in < col; ++in) {
                            sum += TENSOR[in][i2][i3][i4] * MATRIX[j][in];
                        }
                        TENSOR_RESULT[j][i2][i3][i4] = sum;
                    }
                }
            }
        }
    }
    if (n == 2) {
        for (int i1 = 0; i1 < order + 1 ; ++i1) {
            for (int j = 0; j < order + 1 ; ++j) {
                for (int i3 = 0; i3 < order + 1 ; ++i3) {
                    for (int i4 = 0; i4 < order + 1 ; ++i4) {
                        double sum = 0;
                        for (int in = 0; in < col; ++in) {
                            sum += TENSOR[i1][in][i3][i4] * MATRIX[j][in];
                        }
                        TENSOR_RESULT[i1][j][i3][i4] = sum;
                    }
                }
            }
        }
    }
    if (n == 3) {
        for (int i1 = 0; i1 < order + 1 ; ++i1) {
            for (int i2 = 0; i2 < order + 1 ; ++i2) {
                for (int j = 0; j < order + 1 ; ++j) {
                    for (int i4 = 0; i4 < order + 1 ; ++i4) {
                        double sum = 0;
                        for (int in = 0; in < col; ++in) {
                            sum += TENSOR[i1][i2][in][i4] * MATRIX[j][in];
                        }
                        TENSOR_RESULT[i1][i2][j][i4] = sum;
                    }
                }
            }
        }
    }
    if (n == 4) {
        for (int i1 = 0; i1 < order + 1 ; ++i1) {
            for (int i2 = 0; i2 < order + 1 ; ++i2) {
                for (int i3 = 0; i3 < order + 1 ; ++i3) {
                    for (int j = 0; j < order + 1 ; ++j) {
                        double sum = 0;
                        for (int in = 0; in < col; ++in) {
                            sum += TENSOR[i1][i2][i3][in] * MATRIX[j][in];
                        }
                        TENSOR_RESULT[j][i2][i3][j] = sum;
                    }
                }
            }
        }
    }
    return TENSOR_RESULT;
};


int main() {

    // Hardcoded Nodes and Weight for Quadrature
    array<double, 2> vertex = {0, 1};

    array<double, 2> nodes = {0.78868, 0.21132};

    array<double, 2> weights = {0.5 , 0.5};

    // Initialize Tensor
    Tensor4<double, 2> M;

    // Calculate Tensor Entries
    M = computeTensor<double, 1>(nodes, weights, vertex);

    // Transform Tensor to Matrix
    Tensor2<double, 4> MM;
    MM = tensorToMatrix<double, 1>(M);

    // Mode-k flattening
    array<array<double, 2>, 8> MATRIX;
    mode_<double, 1>(1, M, MATRIX);

    cout << endl << "Tensor" << "\n" << M << endl << "\n\n";
    cout << endl << "Tensor" << "\n" ;
    print_Tensor4<double,1>(M);
    cout << "Mass Matrix derived from Tensor" << "\n" << MM << endl << "\n" ;

    // Print Mode-k flattening
    cout << "Mode-1 flattening" << endl;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 8; ++j) {
            cout << MATRIX[i][j] << " " ;
        }
        cout << endl;
    }
    cout << "\n";

    // Test n-Mode Product
    array<array<double, 2>, 2> MATRIX2;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            MATRIX2[i][j] = 1;
        }
    }

    Tensor4<double, 2> M_PROD;
    M_PROD = n_mode_product<double, 1, 2>(1, M, MATRIX2);
    cout << "N Mode Product:" << endl << M_PROD << endl;



// Test SVD
    /*MatrixXf m = MatrixXf::Random(3, 2);
    cout << "Here is the matrix m:" << endl << m << endl;
    JacobiSVD<MatrixXf> svd(m, ComputeFullU | ComputeFullV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
    Vector3f rhs(1, 0, 0);
    cout << "Now consider this rhs vector:" << endl << rhs << endl;
    cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
    cout << "Consider U * U.transpose:" << endl << svd.matrixU()*svd.matrixU().transpose() << endl;
    */
}