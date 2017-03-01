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
                    for (int k = 0; k < order + 1 ; ++k) { // Sum for Quadrature Points
                        for (int l = 0; l < order + 1; ++l) { // Sum for Quadrature Points
                            int sum_l = 0;
                            sum_l = weights[l] * eval_lagr(i2, nodes[l], vertex) * eval_lagr(j2, nodes[l], vertex);
                        }
                        sum += weights[k] * eval_lagr(i1, nodes[k], vertex) * eval_lagr(j1, nodes[k], vertex);
                    }
                    M[i1][i2][j1][j2] = sum;
                }
            }
        }
    }
    return M;
};

// Transform number of 2D function to numbers of the corresponding 1D functions
std::array<int,2> transform(int k, int order){
    int n = order+1;
    // Find row 
    int row = k/n;
    // Find Column
    int column = k-row*n;

    std::array<int,2> index={row,column};
    return index;
};

// Transform number of 1D functions to numbers of the corresponding 2D function
int transform(int row,int column, int order){
    int k = row*(order+1)+column;
    return k;
};

std::array<int,2> transform(int i1,int i2, int j1, int j2, int order){
    int k1 = transform(i1,i2,order);
    int k2 = transform(j1,j2,order);
    std::array<int,2> k = {k1,k2};
    return k;
};



// Transform Mass Matrix Tensor to Mass Matrix
template<typename y_type, int order>
Tensor2<y_type,(order+1)*(order+1)> tensorToMatrix(Tensor4<y_type,order+1> & M){
    Tensor2<y_type,(order+1)*(order+1)> MM;
    for (int i1=0;i1<order+1;++i1){
        for (int i2=0;i2<order+1;++i2){
            for (int j1=0;j1<order+1;++j1){
                for (int j2=0;j2<order+1;++j2){
                    std::array<int,2> k=transform(i1,i2,j1,j2,order);
                    MM[k[0]][k[1]] = M[i1][i2][j1][j2];
                }
            }
        }
    }
    return MM;
};


int main() {

    // Hardcoded Nodes and Weight for Quadrature
    array<double, 2> vertex = {0, 1};

    array<double, 2> nodes = {0.78868, 0.21132};

    array<double, 2> weights = {0.5 , 0.5};

    // Initialize Tensor
    Tensor4<double, 2> M;

    // Calculate Tensor Entries
    M = computeTensor<double, 1>(nodes,weights,vertex);

    // Transform Tensor to Matrix
    Tensor2<double, 4> MM;
    MM = tensorToMatrix<double, 1>(M);

    cout << "Tensor" << "\n" << M << endl << "\n\n";
    cout << "Mass Matrix derived from Tensor" << "\n" << MM << endl;


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