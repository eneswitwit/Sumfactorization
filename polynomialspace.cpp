#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <numeric>



class Polynomialspace{

    // Coefficient vector for newton basis
    std::vector<double> coefficients;

    // Data
    std::vector<std::vector<double>> X;
    std::vector<double> y;

    public:
    // Constructors
    Polynomialspace(std::vector<std::vector<double>> XX, std::vector<double> yy){
        X=XX;
        y=yy;
        compute_coefficients();
    }

    Polynomialspace(std::vector<double> Coefficients, std::vector<std::vector<double>> XX, std::vector<double> yy){
        coefficients=Coefficients;
        X=XX;
        y=yy;
    }
    // Define direction vector from datapoint i to datapoint j
    double u(int i, int j, std::vector<double> point) const {
        double x_1 = (X[j][0]-X[i][0])*(point[0]-X[i][0]);
        double x_2 = (X[j][1]-X[i][1])*(point[1]-X[i][1]);
        return x_1+x_2;
    }
    // Evaluating function
    double operator()(double x,double y) const {
        double val = coefficients[0];
        double product = 1;
        for(int i=1; i < coefficients.size() ; i++){
            for(int j=0; j < i; j++){
                product = product * u(j,i,{x,y});
            }
            val = val + coefficients[i] * product;
        }
        return val;
    }
    // Define function for computing coefficients for newton basis
    void compute_coefficients(){
        coefficients = {};
        // First coefficient is the first value 
        coefficients.push_back(y[0]);
        for(int i=1; i < y.size(); i++){
            Polynomialspace f(coefficients,X,y);
            double c_i = y[i] - f(X[i][0],X[i][1]);
            for(int j=0; j < i; j++){
                c_i = c_i / u(j,i,X[i]);
            }
            coefficients.push_back(c_i);
        }
    }

    // Standard operations
    bool operator==(Polynomialspace const & other) const
    {
        return coefficients == other.coefficients;
    }

    // Read or write for coefficient with number i
    int & operator[](int i)
    {
        return coefficients[i];
    }
    int operator[](int i) const
    {
        return coefficients[i];
    }


};

int main(){

    // Test data
    std::vector<std::vector<double>> X={{1,1},{1,0},{0,1},{0,0}};
    std::vector<double> y={1,1,1,0};

    Polynomialspace p(X,y);
    // Test whether divided difference algorithm works correctly
    for(int i=0;i<X.size();i++){
        assert(p(X[i][0],X[i][1]) == y[i]);
    }

    return 0;
}