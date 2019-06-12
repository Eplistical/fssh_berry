#ifndef _2D_FLAT_PHIXPLUSY_POTENTIAL_HPP
#define _2D_FLAT_PHIXPLUSY_POTENTIAL_HPP

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
        ioer::info("# 2D flat potential [with phi = W * (x + y)] parameters: ", 
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

    double cal_theta(const vector<double>& r) {
        const double x = r[0];
        return 0.5 * M_PI * (erf(param_B * x) + 1);
    }

    vector<double> cal_der_theta(const vector<double>& r) {
        const double x = r[0];
        vector<double> der_theta(r.size(), 0.0);
        der_theta[0] = sqrt(M_PI) * param_B * exp(-param_B * param_B * x * x);
        return der_theta;
    }

    double cal_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        return param_W * (x + y);
    }

    vector<double> cal_der_phi(const vector<double>& r) {
        vector<double> der_phi(r.size(), 0.0);
        der_phi[0] = param_W;
        der_phi[1] = param_W;
        return der_phi;
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double theta = cal_theta(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        vector< complex<double> > H {
            -cos(theta), sin(theta) * conj(eip), sin(theta) * eip, cos(theta)
        };
        return param_A * H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double theta = cal_theta(r);
        const vector<double> der_theta = cal_der_theta(r);
        const vector<double> der_phi = cal_der_phi(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > nablaHx(4, 0.0);
        nablaHx[0+0*2] = param_A * sin(theta) * der_theta[0];
        nablaHx[0+1*2] = param_A * eip * (cos(theta) * der_theta[0] + matrixop::IMAGIZ * sin(theta) * der_phi[0]);
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        nablaHx[1+1*2] = param_A * -sin(theta) * der_theta[0];
        return nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        const double theta = cal_theta(r);
        const vector<double> der_theta = cal_der_theta(r);
        const vector<double> der_phi = cal_der_phi(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));

        vector< complex<double> > nablaHy(4, 0.0);
        nablaHy[0+0*2] = param_A * sin(theta) * der_theta[1];
        nablaHy[0+1*2] = param_A * eip * (cos(theta) * der_theta[1] + matrixop::IMAGIZ * sin(theta) * der_phi[1]);
        nablaHy[1+0*2] = conj(nablaHy[0+1*2]);
        nablaHy[1+1*2] = param_A * -sin(theta) * der_theta[1];
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
