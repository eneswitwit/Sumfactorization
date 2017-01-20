#include <iostream>
#include <vector>
using namespace std;

class Polynomial{

    public:

    int degree;
    int nodes_counter;
    vector<double> coefficients[0];
    vector<double> sample_points[0][0];
    vector<double> sample_values[0];

    // Constructor
    Polynomial(int deg){
        degree = deg;
    }  

    void sample_points (){}

    double evaluation(){}

    void divided_differences(){}


};



int main(){
    Polynomial poly_degree_2(2);
    std::cout << poly_degree_2.getDegree() ;
    return 0;
}