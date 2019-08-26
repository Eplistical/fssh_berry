#ifndef _2D_HELIX_POTENTIAL_HPP
#define _2D_HELIX_POTENTIAL_HPP

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

    double param_A = 0.1;
    double param_B = 0.2;
    double param_C = 0.001;
    double param_D = 0.0;
    double param_Wx = 0.0;
    double param_Wy = 0.0;

    /*
     * H00 = 0.5 * A * x**2 + 0.5 * (y - B)**2 + 0.5 * D
     * H11 = 0.5 * (x - B)**2 + 0.5 * A * y**2 - 0.5 * D
     * H01 = C * exp(i * (Wx * x + Wy * y))
     */

    void output_potential_param() {
        ioer::info("# 2D Helix potential [with phi = Wx * x + Wy * y] parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " C = ", param_C,
                    " D = ", param_D,
                    " Wx = ", param_Wx,
                    " Wy = ", param_Wy);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 6, 
                "set_potenial_params: potential paramter vector size must be >= 6");
        param_A = params[0];
        param_B = params[1];
        param_C = params[2];
        param_D = params[3];
        param_Wx = params[4];
        param_Wy = params[5];
    }

    double cal_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        return param_Wx * x + param_Wy * y;
    }

    vector<double> cal_der_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        vector<double> der_phi(r.size(), 0.0);
        der_phi[0] = param_Wx;
        der_phi[1] = param_Wy;
        return der_phi;
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > H(4);
        H[0+0*2] = 0.5 * (param_A * pow(x, 2) + pow(y - param_B, 2) + param_D);
        H[1+1*2] = 0.5 * (pow(x - param_B, 2) + param_A * pow(y, 2) - param_D);
        H[0+1*2] = param_C * eip;
        H[1+0*2] = conj(H[0+1*2]);
        return H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const vector<double> der_phi = cal_der_phi(r);
        vector< complex<double> > nablaHx(4);

        nablaHx[0+0*2] = param_A * x;
        nablaHx[1+1*2] = x - param_B;
        nablaHx[0+1*2] = matrixop::IMAGIZ * param_C * eip * der_phi[0];
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        return nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const vector<double> der_phi = cal_der_phi(r);
        vector< complex<double> > nablaHy(4);

        nablaHy[0+0*2] = y - param_B;
        nablaHy[1+1*2] = param_A * y;
        nablaHy[0+1*2] = matrixop::IMAGIZ * param_C * eip * der_phi[1];
        nablaHy[1+0*2] = conj(nablaHy[0+1*2]);

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
