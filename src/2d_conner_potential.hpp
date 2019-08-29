#ifndef _2D_CONNER_POTENTIAL_HPP
#define _2D_CONNER_POTENTIAL_HPP

#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"


namespace {
    using std::vector;
    using std::complex;

    double param_A = 0.05;
    double param_B = 4.0;
    double param_C = 0.15;
    double param_D = 2.0;
    double param_E = 4.0;
    double param_W = 0.5;

    void output_potential_param() {
        ioer::info("# 2D conner potential [with phi = W * (x + y)] parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " C = ", param_C,
                    " D = ", param_D,
                    " E = ", param_E,
                    " W = ", param_W);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 6, 
                "set_potenial_params: potential paramter vector size must be >= 6");
        param_A = params[0];
        param_B = params[1];
        param_C = params[2];
        param_D = params[3];
        param_E = params[4];
        param_W = params[5];
    }

    double cal_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        return param_W * (x + y);
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const double barrier = param_A * pow(param_E, 2) / (pow(x - y, 2) + pow(param_E, 2));

        vector< complex<double> > H(2,2);
        H[0+0*2] = tanh(x-param_B) - tanh(x+param_B) + param_D + tanh(y) + barrier;
        H[1+1*2] = tanh(y-param_B) - tanh(y+param_B) + param_D + tanh(x) + barrier;
        H[0+1*2] = param_C * eip;
        H[1+0*2] = conj(H[0+1*2]);
        return H;
    }

};

#endif
