#ifndef _2D_YANZE_POTENTIAL_HPP
#define _2D_YANZE_POTENTIAL_HPP

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

    double param_A = 0.2;
    double param_C = 0.02;
    double param_r = 1.0;
    double param_R = 2.0;
    double param_eps = 5.0;
    double param_forbidden = 10.0;
    double param_alpha = 200.0;
    double param_W = 0.0;

    void output_potential_param() {
        ioer::info("# 2D Yanze circular potential [with phi = W * sqrt((x+param_R/2)^2 + (y+param_R/2)^2)] parameters: ", 
                    " A = ", param_A,
                    " C = ", param_C,
                    " r = ", param_r,
                    " R = ", param_R,
                    " eps = ", param_eps,
                    " forbidden = ", param_forbidden,
                    " alpha = ", param_alpha,
                    " W = ", param_W
                    );
    }

    void set_potenial_params(const std::vector<double>& params) {
        misc::crasher::confirm(params.size() >= 8, 
                "set_potenial_params: potential paramter vector size must be >= 8");
        param_A = params[0];
        param_C = params[1];
        param_r = params[2];
        param_R = params[3];
        param_eps = params[4];
        param_forbidden = params[5];
        param_alpha = params[6];
        param_W = params[7];
    }

    double cal_m_r(const vector<double>& r) {
        return sqrt( pow((r[0] + 0.5 * param_R), 2) + pow((r[1] + 0.5 * param_R), 2) );
    }

    vector<double> cal_der_m_r(const vector<double>& r) {
        return vector<double> { r[0] + 0.5 * param_R, r[1] + 0.5 * param_R } / cal_m_r(r);
    }

    double cal_m_theta(const vector<double>& r) {
        return atan2(r[1] + 0.5 * param_R, r[0] + 0.5 * param_R);
    }

    vector<double> cal_der_m_theta(const vector<double>& r) {
        const double x1 = r[0] + 0.5 * param_R;
        const double y1 = r[1] + 0.5 * param_R;
        return vector<double> { -y1, x1 } / (x1*x1 + y1*y1);
    }

    double cal_phi(const vector<double>& r) {
        return param_W * cal_m_r(r);
    }

    vector<double> cal_der_phi(const vector<double>& r) {
        return param_W * cal_der_m_r(r);
    }

    vector< complex<double> > cal_H(const vector<double>& r) {
        const double m_r = cal_m_r(r);
        const double m_theta = cal_m_theta(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        vector< complex<double> > H(4);
        H[0+0*2] = param_A * (tanh((m_theta - M_PI/4) * param_eps) + 1.0) 
                    + param_forbidden * (2.0 + tanh((m_r - param_R) * param_alpha) - tanh((m_r - param_r) * param_alpha)); 
        H[1+1*2] = param_A * (tanh(-(m_theta - M_PI/4) * param_eps) + 1.0) 
                    + param_forbidden * (2.0 + tanh((m_r - param_R) * param_alpha) - tanh((m_r - param_r) * param_alpha)); 
        H[0+1*2] = param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2)) * eip
                    * 0.5 * (tanh((m_r - param_r) * param_alpha) - tanh((m_r - param_R) * param_alpha));
        H[1+0*2] = conj(H[0+1*2]);
        return H;
    }

    vector< complex<double> > cal_nablaHx(const vector<double>& r) {
        const double m_r = cal_m_r(r);
        const double m_theta = cal_m_theta(r);
        const vector<double> der_m_r = cal_der_m_r(r);
        const vector<double> der_m_theta = cal_der_m_theta(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const vector<double> der_phi = cal_der_phi(r);

        vector< complex<double> > nablaHx(4);
        auto der_tanh = [](double x) { 
            return (1.0 - pow(tanh(x), 2)); 
        };
        nablaHx[0+0*2] = param_A * der_tanh((m_theta - M_PI/4) * param_eps) * param_eps * der_m_theta[0]
                            + param_forbidden * ( 
                                    der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[0] 
                                    - der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[0] 
                                  );
        nablaHx[1+1*2] = param_A * der_tanh((m_theta - M_PI/4) * -param_eps) * -param_eps * der_m_theta[0]
                            + param_forbidden * ( 
                                    der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[0] 
                                    - der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[0] 
                                  );
        // (e^{i*phi} * f)' = eip * (i*phi'*f + f')
        nablaHx[0+1*2] = eip * (
                matrixop::IMAGIZ * der_phi[0] 
                    * param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2)) 
                    * 0.5 * (tanh((m_r - param_r) * param_alpha) - tanh((m_r - param_R) * param_alpha))
                + param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2)) * (-2.0) * param_eps * (m_theta - M_PI/4) * der_m_theta[0]
                    * 0.5 * (tanh((m_r - param_r) * param_alpha) - tanh((m_r - param_R) * param_alpha))
                + param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2))
                    * 0.5 * (der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[0] * - der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[0])
                );
        nablaHx[1+0*2] = conj(nablaHx[0+1*2]);
        return nablaHx;
    }

    vector< complex<double> > cal_nablaHy(const vector<double>& r) {
        const double m_r = cal_m_r(r);
        const double m_theta = cal_m_theta(r);
        const vector<double> der_m_r = cal_der_m_r(r);
        const vector<double> der_m_theta = cal_der_m_theta(r);
        const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r));
        const vector<double> der_phi = cal_der_phi(r);

        vector< complex<double> > nablaHy(4);
        // (tanh(f))' = tanh'(f) * f' 
        auto der_tanh = [](double x) { 
            return (1.0 - pow(tanh(x), 2)); 
        };
        nablaHy[0+0*2] = param_A * der_tanh((m_theta - M_PI/4) * param_eps) * param_eps * der_m_theta[1]
                            + param_forbidden * ( 
                                    der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[1] 
                                    - der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[1] 
                                  );
        nablaHy[1+1*2] = param_A * der_tanh((m_theta - M_PI/4) * -param_eps) * -param_eps * der_m_theta[1]
                            + param_forbidden * ( 
                                    der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[1] 
                                    - der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[1] 
                                  );
        // (e^{i*phi} * f)' = eip * (i*phi'*f + f')
        nablaHy[0+1*2] = eip * (
                matrixop::IMAGIZ * der_phi[1] 
                    * param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2)) 
                    * 0.5 * (tanh((m_r - param_r) * param_alpha) - tanh((m_r - param_R) * param_alpha))
                + param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2)) * (-2.0) * param_eps * (m_theta - M_PI/4) * der_m_theta[1]
                    * 0.5 * (tanh((m_r - param_r) * param_alpha) - tanh((m_r - param_R) * param_alpha))
                + param_C * exp(-param_eps * pow((m_theta - M_PI/4), 2))
                    * 0.5 * (der_tanh((m_r - param_r) * param_alpha) * param_alpha * der_m_r[1] * - der_tanh((m_r - param_R) * param_alpha) * param_alpha * der_m_r[1])
                );
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
