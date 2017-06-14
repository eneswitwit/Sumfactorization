#include "../include/constexpr_quadrature.h"
#include <iostream>
#include <cassert>



constexpr static int order = 9;

/**
 * For this test, we will check every entry with a hard coded solution.
 * We gathered the rule for 10 points from: http://keisan.casio.com/exec/system/1280801905
 * A transformation is needed, since the rule is written for the interval [-1,1], our rule works on [0,1]
 */

int main() {
    const long double eps=0.00000000001;
    constexpr Quadrature<long double, order> quad;

    assert(abs(quad.knots_[0]*2-1-(-1))<eps);
    assert(abs(quad.knots_[1]*2-1-(-0.9195339081664588138289))<eps);
    assert(abs(quad.knots_[2]*2-1-(-0.7387738651055050750031))<eps);
    assert(abs(quad.knots_[3]*2-1-(-0.4779249498104444956612))<eps);
    assert(abs(quad.knots_[4]*2-1-(-0.1652789576663870246262))<eps);
    assert(abs(quad.knots_[5]*2-1-(0.1652789576663870246262))<eps);
    assert(abs(quad.knots_[6]*2-1-(0.4779249498104444956612))<eps);
    assert(abs(quad.knots_[7]*2-1-(0.7387738651055050750031))<eps);
    assert(abs(quad.knots_[8]*2-1-(0.9195339081664588138289))<eps);
    assert(abs(quad.knots_[9]*2-1-(1))<eps);

    std::cout << "quad.points   test successful" << std::endl;

    assert(abs(quad.weights_[0]*2-(0.02222222222222222222222))<eps);
    assert(abs(quad.weights_[1]*2-(0.1333059908510701111262))<eps);
    assert(abs(quad.weights_[2]*2-(0.2248893420631264521195))<eps);
    assert(abs(quad.weights_[3]*2-(0.292042683679683757876))<eps);
    assert(abs(quad.weights_[4]*2-(0.327539761183897456657))<eps);
    assert(abs(quad.weights_[5]*2-(0.327539761183897456657))<eps);
    assert(abs(quad.weights_[6]*2-(0.292042683679683757876))<eps);
    assert(abs(quad.weights_[7]*2-(0.2248893420631264521195))<eps);
    assert(abs(quad.weights_[8]*2-(0.1333059908510701111262))<eps);
    assert(abs(quad.weights_[9]*2-(0.02222222222222222222222))<eps);

    std::cout << "quad.weights  test successful" << std::endl;

    return 0;
}
