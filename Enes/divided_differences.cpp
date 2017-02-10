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
}

int main(){


    return 0;
}