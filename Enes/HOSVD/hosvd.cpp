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

int main() {

    // Hardcoded Nodes and Weight for Quadrature
    array<double, 2> vertex = {0, 1};

    array<double, 2> nodes = {0.78868, 0.21132};

    array<double, 2> weights = {0.5 , 0.5};

    // Initialize Tensor
    Tensor4<double, 2> M;

    // Calculate Tensor Entries
    for (int i1 = 0 ; i1 < 2 ; ++i1 ) {
        for (int i2 = 0 ; i2 < 2 ; ++i2 ) {
            for (int j1 = 0 ; j1 < 2 ; ++j1 ) {
                for (int j2 = 0 ; j2 < 2 ; ++j2) {
                    double sum = 0;
                    for (int k = 0; k < 2 ; ++k) { // Sum for Quadrature Points
                        for (int l = 0; l < 2; ++l) { // Sum for Quadrature Points
                            int sum_l = 0;
                            sum_l = weights[l] * eval_lagr(i2, nodes[l], vertex) * eval_lagr(j2, nodes[l], vertex);
                        }
                        sum += weights[k] * eval_lagr(i1, nodes[k], vertex) * eval_lagr(j1, nodes[k], vertex);
                    }
                    M[i1, i2, j1, j2] = sum;
                }
            }
        }
    }

    cout << M;

    // Compute Mass Matrix
    Tensor2<double, 2> MM;
    MM[1,1] = M[1,1,1,1];
    MM[1,2] = M[1,1,2,1];
    MM[2,1] = M[2,1,1,1];
    MM[2,2] = M[2,1,2,1];

    cout << MM;

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