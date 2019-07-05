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
        ioer::info("# 2D flat potential [with phi = W * y] parameters: ", 
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

    vector<double> cal_derder_theta(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        vector<double> derder_theta(2 * 2, 0.0);
        derder_theta[0+0*2] = -2.0 * sqrt(M_PI) * pow(param_B, 3) * x * exp(-param_B * param_B * x * x); // d2theta / dxdx
        derder_theta[0+1*2] = 0.0; // d2theta / dxdy
        derder_theta[1+0*2] = 0.0; // d2theta / dydx 
        derder_theta[1+1*2] = 0.0; // d2theta / dydy
        return derder_theta;
    }

    double cal_phi(const vector<double>& r) {
        const double y = r[1];
        return param_W * y;
    }

    vector<double> cal_der_phi(const vector<double>& r) {
        vector<double> der_phi(r.size(), 0.0);
        der_phi[1] = param_W;
        return der_phi;
    }

    vector<double> cal_derder_phi(const vector<double>& r) {
        const double x = r[0];
        const double y = r[1];
        vector<double> derder_phi(2 * 2, 0.0);
        derder_phi[0+0*2] = 0.0; // d2phi / dxdx
        derder_phi[0+1*2] = 0.0; // d2phi / dxdy
        derder_phi[1+0*2] = 0.0; // d2phi / dydx 
        derder_phi[1+1*2] = 0.0; // d2phi / dydy
        return derder_phi;
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

    void cal_info_anal(const vector<double>& r, 
            vector< complex<double> >& Fx, vector< complex<double> >& Fy, 
            vector< complex<double> >& dcx, vector< complex<double> >& dcy, 
            vector< complex<double> >& dx_dcx, vector< complex<double> >& dy_dcy, 
            vector<double>& eva, vector< complex<double> >& lastevt
            )
    {
        const double theta = cal_theta(r);
        const double phi = cal_phi(r);
        const vector<double> der_theta = cal_der_theta(r);
        const vector<double> der_phi= cal_der_phi(r);
        const vector<double> derder_theta = cal_derder_theta(r);
        const vector<double> derder_phi= cal_derder_phi(r);
        const double CC = cos(0.5 * theta);
        const double SS = sin(0.5 * theta);

        // eva
        eva.assign(2, 0.0);
        eva[0] = -param_A;
        eva[1] = param_A;

        // dc
        dcx.assign(4, 0.0);
        dcx[0+0*2] = matrixop::IMAGIZ * der_phi[0] * CC * CC;
        dcx[0+1*2] = 0.5 * der_theta[0] + matrixop::IMAGIZ * der_phi[0] * SS * CC;
        dcx[1+0*2] = -conj(dcx[0+1*2]);
        dcx[1+1*2] = matrixop::IMAGIZ * der_phi[0] * SS * SS;

        dcy.assign(4, 0.0);
        dcy[0+0*2] = matrixop::IMAGIZ * der_phi[1] * CC * CC;
        dcy[0+1*2] = 0.5 * der_theta[1] + matrixop::IMAGIZ * der_phi[1] * SS * CC;
        dcy[1+0*2] = -conj(dcy[0+1*2]);
        dcy[1+1*2] = matrixop::IMAGIZ * der_phi[1] * SS * SS;

        // F
        Fx.assign(4, 0.0);
        Fx[0+1*2] = dcx[0+1*2] * (eva[0] - eva[1]);
        Fx[1+0*2] = dcx[1+0*2] * (eva[1] - eva[0]);

        Fy.assign(4, 0.0);
        Fy[0+1*2] = dcy[0+1*2] * (eva[0] - eva[1]);
        Fy[1+0*2] = dcy[1+0*2] * (eva[1] - eva[0]);

        // dx_dcx, dy_dcy
        dx_dcx.assign(4, 0.0);
        dx_dcx[0+0*2] = matrixop::IMAGIZ * CC * (CC * derder_phi[0+0*2] - SS * der_phi[0] * der_theta[0]);
        dx_dcx[1+1*2] = matrixop::IMAGIZ * SS * (SS * derder_phi[0+0*2] + CC * der_phi[0] * der_theta[0]);
        dx_dcx[0+1*2] = 0.5 * (derder_theta[0+0*2] 
                    + matrixop::IMAGIZ * sin(theta) * derder_phi[0+0*2]
                    + matrixop::IMAGIZ * cos(theta) * der_phi[0] * der_theta[0]
                    );
        dx_dcx[1+0*2] = -conj(dx_dcx[0+1*2]);

        dy_dcy.assign(4, 0.0);
        dy_dcy[0+0*2] = matrixop::IMAGIZ * CC * (CC * derder_phi[1+1*2] - SS * der_phi[1] * der_theta[1]);
        dy_dcy[1+1*2] = matrixop::IMAGIZ * SS * (SS * derder_phi[1+1*2] + CC * der_phi[1] * der_theta[1]);
        dy_dcy[0+1*2] = 0.5 * (derder_theta[1+1*2] 
                    + matrixop::IMAGIZ * sin(theta) * derder_phi[1+1*2]
                    + matrixop::IMAGIZ * cos(theta) * der_theta[1] * der_phi[1]
                    );
        dy_dcy[1+0*2] = -conj(dy_dcy[0+1*2]);
    }
};

#endif
