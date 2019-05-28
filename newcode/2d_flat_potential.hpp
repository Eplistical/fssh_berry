#ifndef _2D_FLAT_POTENTIAL_HPP
#define _2D_FLAT_POTENTIAL_HPP

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
#include "boost/math/special_functions/erf.hpp"

namespace {
    using std::vector;
    using std::complex;

    double param_A = 0.10;
    double param_B = 3.0;
    double param_W = 0.0;

    void output_potential_param() {
        ioer::info("# 2D flat potential parameters: ", 
                    " A = ", param_A,
                    " B = ", param_B,
                    " W = ", param_W);
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 3, 
                "set_potenial_params: potential paramter vector size must be >= 3");
        param_A = params[0];
        param_B = params[1];
        param_W = params[2];
    }

    double cal_theta(double x) {
        return 0.5 * M_PI * (erf(param_B * x) + 1);
    }

    double cal_der_theta(double x) {
        return sqrt(M_PI) * param_B * exp(-param_B * param_B * x * x);
    }

    double cal_der2_theta(double x) {
        return -2.0 * sqrt(M_PI) * pow(param_B, 3) * x * exp(-param_B * param_B * x * x);
    }

    double cal_phi(double y) {
        return param_W * y;
    }

    double cal_der_phi(double y) {
        return param_W;
    }

    double cal_der2_phi(double y) {
        return 0.0;
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double theta = cal_theta(r[0]);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
        vector< complex<double> > H {
            -cos(theta), sin(theta) * conj(eip), sin(theta) * eip, cos(theta)
        };
        return param_A * H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double theta = cal_theta(r[0]);
        const double der_theta = cal_der_theta(r[0]);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
        vector< complex<double> > nablaHx {
            sin(theta), cos(theta) * conj(eip), cos(theta) * eip, -sin(theta)
        };
        return param_A * der_theta * nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        const double theta = cal_theta(r[0]);
        const double der_theta = cal_der_theta(r[0]);
        const double der_phi = cal_der_phi(r[1]);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
        vector< complex<double> > nablaHy {
            0.0, -sin(theta) * conj(eip), sin(theta) * eip, 0.0
        };
        return matrixop::IMAGIZ * param_A * der_phi * nablaHy;
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

    void cal_info_anal(const vector<double>& r, 
            vector< complex<double> >& Fx, vector< complex<double> >& Fy, 
            vector< complex<double> >& dcx, vector< complex<double> >& dcy, 
            vector<double>& eva
            )
    {
        const double theta = cal_theta(r[0]);
        const double der_theta = cal_der_theta(r[0]);
        const double der_phi= cal_der_phi(r[1]);
        const double CC = cos(0.5 * theta);
        const double SS = sin(0.5 * theta);

        // eva
        eva.assign(2, 0.0);
        eva[0] = -param_A;
        eva[1] = param_A;

        // dc
        dcx.assign(4, 0.0);
        dcx[0+0*2] = 0.0;
        dcx[0+1*2] = 0.5 * der_theta;
        dcx[1+0*2] = -0.5 * der_theta;
        dcx[1+1*2] = 0.0;

        dcy.assign(4, 0.0);
        dcy[0+0*2] = matrixop::IMAGIZ * der_phi * CC * CC;
        dcy[0+1*2] = matrixop::IMAGIZ * der_phi * SS * CC;
        dcy[1+0*2] = matrixop::IMAGIZ * der_phi * SS * CC;
        dcy[1+1*2] = matrixop::IMAGIZ * der_phi * SS * SS;

        // F
        Fx.assign(4, 0.0);
        Fx[0+1*2] = dcx[0+1*2] * (eva[0] - eva[1]);
        Fx[1+0*2] = dcx[1+0*2] * (eva[1] - eva[0]);

        Fy.assign(4, 0.0);
        Fy[0+1*2] = dcy[0+1*2] * (eva[0] - eva[1]);
        Fy[1+0*2] = dcy[1+0*2] * (eva[1] - eva[0]);
    }
};

#endif
