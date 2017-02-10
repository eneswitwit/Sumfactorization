#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>

#define PI 3.141592653589793238463

// Initialize Coefficient matrix
std::vector<std::vector<double>> Coefficient={};

// Initialize data_set
std::vector<std::vector<double>> data_set={{0,0,1},{0,0.5,0},{0,1,0},{0.5,0,0},{1,0,0},{0.5,0.5,0},{1,0.5,0},
{0.5,1,0},{1,1,0}};

// Define w_{i,j}(x,y)
double w(int i, int j, double x, double y){
    double val = 1;
    for(unsigned k=0; k<i; k++){
        val=val*(x-data_set[k][0]);
    }
    for(unsigned k=0; k<j; k++){
        val=val*(y-data_set[k][1]);
    }
    return val;
}

//Define evaluation function
double p(std::vector<double> coefficients, double x, double y){
    double val = coefficients[0];
    for(int i=0;i<coefficients.size();i++){
        val = val + coefficients[i]*w(i,i,x,y);
    }
    return val;
}

//Define function computing coefficients
std::vector<double> compute_coefficients(std::vector<std::vector<double>> data){
    std::vector<double> coefficients = {};
    coefficients.push_back(data[0][2]);
    for(int i=1;i<data.size();i++){
        double c_i = (data[i][2]-p(coefficients,data[i][0],data[i][1])) / w(i,i,data[i][0],data[i][1]) ;
        coefficients.push_back(c_i);
    }
    return coefficients;
}

int main(){
    std::vector<double> Coefficient_vector = compute_coefficients(data_set);
    for(int i=0;i<Coefficient_vector.size();i++){
        std::cout << Coefficient_vector[i] << "\n" ;
    }

    std::cout << "\n" << p(Coefficient_vector,1,1) << "\n" ;
    return 0;

}