#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>

#define PI 3.141592653589793238463

// Initialize Coefficient matrix
std::vector<std::vector<double>> Coefficient={};

// Initialize data_set
std::vector<std::vector<double>> X={
    {1,1},
    {1,0},
    {0,1},
    {0,0}
};
std::vector<double> y={
    1,
    1,
    1,
    0
};

// Define direction vector from datapoint i to datapoint j
double u(int i, int j, std::vector<double> point){
    double x_1 = (X[j][0]-X[i][0])*(point[0]-X[i][0]);
    double x_2 = (X[j][1]-X[i][1])*(point[1]-X[i][1]);
    return x_1+x_2;
}

// Define evaluation function
double f(std::vector<double> coefficients, std::vector<double> point){
    double val = coefficients[0];
    double product = 1;
    for(int i=1; i < coefficients.size() ; i++){
        for(int j=0; j < i; j++){
            product = product * u(j,i,point);
        }
        val = val + coefficients[i] * product;
    }
}

// Define function for calculating new coefficient depenending on the previous coefficients
std::vector<double> compute_coefficients(std::vector<std::vector<double>> X_data, std::vector<double> y_data){
    std::vector<double> coefficients = {};
    // First coefficient is the first value 
    coefficients.push_back(y[0]);
    for(int i=1; i < y_data.size(); i++){
        double c_i = y[i] - f(coefficients,X_data[i]);
        for(int j=0; j < i; j++){
            c_i = c_i / u(j,i,X_data[i]);
        }
        coefficients.push_back(c_i);
    }
    return coefficients;
}


int main(){
    std::vector<double> Coefficient_vector = compute_coefficients(X,y);
    for(int i=0;i<Coefficient_vector.size();i++){
        std::cout << Coefficient_vector[i] << "\n" ;
    }

    std::cout << "\n" << f(Coefficient_vector,{1,1}) << "\n" ;
    return 0;
}