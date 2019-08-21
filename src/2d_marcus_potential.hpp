#ifndef _2D_MARCUS_POTENTIAL_HPP
#define _2D_MARCUS_POTENTIAL_HPP

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
    using std::pow;
    using std::complex;

    const double omega = 4.375e-5;
    const double Er = 2.39e-2;
    const double M = pow(0.5 * omega * omega  * Er, 0.5);

    double param_A = omega * omega;
    double param_B = M;
    double param_C = 2.5e-5;
    double param_D = 0.018;

    /*
     * H11 = 0.5 * A * x**2 + B * x
     * H22 = 0.5 * A * x**2 - B * x - D
     * H12 = C
     */

    void output_potential_param() {
        ioer::info("# 2D Marcus potential [Landry's Paper: How to recover ...] parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " C = ", param_C,
                    " D = ", param_D);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 4, 
                "set_potenial_params: potential paramter vector size must be >= 4");
        param_A = params[0];
        param_B = params[1];
        param_C = params[2];
        param_D = params[3];
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];

        vector< complex<double> > H(4);
        H[0+0*2] = 0.5 * param_A * pow(x, 2) + param_B * x;
        H[1+1*2] = 0.5 * param_A * pow(x, 2) - param_B * x - param_D;
        H[0+1*2] = param_C;
        H[1+0*2] = conj(H[0+1*2]);
        return H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        vector< complex<double> > nablaHx(4);

        nablaHx[0+0*2] = param_A * x + param_B;
        nablaHx[1+1*2] = param_A * x - param_B;
        nablaHx[0+1*2] = 0.0;
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        return nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        vector< complex<double> > nablaHy(4, 0.0);
        return nablaHy;
    }

    void cal_info_nume(const vector<double>& r, 
            vector< complex<double> >& Fx, vector< complex<double> >& Fy, 
            vector< complex<double> >& dcx, vector< complex<double> >& dcy, 
            vector<double>& eva, vector< complex<double> >& lastevt 
            )
    {
        vector< complex<double> > evt;
        matrixop::hdiag(cal_H(r), eva, evt);
        // correct phase
        if (not lastevt.empty()) {
            auto tmp = matrixop::matCmat(lastevt, evt, 2);
            for (int j = 0; j < 2; ++j) {
                complex<double> eip = tmp[j+j*2] / abs(tmp[j+j*2]);
                for (int k = 0; k < 2; ++k) {
                    evt[k+j*2] /= eip;
                }
            }
        }
        // F, dc
        dcx = matrixop::matCmatmat(evt, cal_nablaHx(r), evt, 2, 2);
        dcy = matrixop::matCmatmat(evt, cal_nablaHy(r), evt, 2, 2);
        Fx.assign(4, 0.0);
        Fy.assign(4, 0.0);
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                Fx[j+k*2] = -dcx[j+k*2];
                Fy[j+k*2] = -dcy[j+k*2];
                if (j == k) {
                    dcx[j+k*2] = 0.0;
                    dcy[j+k*2] = 0.0;
                }
                else {
                    dcx[j+k*2] /= (eva[k] - eva[j]);
                    dcy[j+k*2] /= (eva[k] - eva[j]);
                }
            }
        }
        lastevt = move(evt);
    }
};

#endif
