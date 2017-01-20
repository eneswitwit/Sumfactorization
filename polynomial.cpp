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
    vector<double> newtonMatrix[0][0];

    // Constructor
    Polynomial(int deg){
        degree = deg;
    }  

    void sample_points (){}

    double evaluation(){}

    void divided_differences(){}

    double product(double x_val, int x,double y_val, int y){
        if(x>0){
            return (x_val - sample_points[x][1])*product(x_val,x-1,y_val,y);
        }
        else{
            return 1;
        }
    }

    void compute_newtonMatrix(double x_val, double y_val){
        for(int i=0;i<nodes_counter+1;++i){
            for(int j=0;j<nodes_counter+1;++j){
                newtonMatrix[i,j]= product(x_val,j,y_val,i)
            }
        }
    }

};



int main(){
    Polynomial poly_degree_2(2);
    std::cout << poly_degree_2.getDegree() ;
    return 0;
}