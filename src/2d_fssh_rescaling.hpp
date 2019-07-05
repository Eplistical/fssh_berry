#ifndef _2D_FSSH_RESCALING_HPP
#define _2D_FSSH_RESCALING_HPP


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
    using std::string;

    vector<double> get_rescaling_direction(const string& rescaling_alg, 
            const vector<double>& r,
            const vector<double>& v,
            const vector< complex<double> >& c,
            const int from, const int to,
            const vector< complex<double> >& dcx,
            const vector< complex<double> >& dcy,
            const vector< complex<double> >& dx_dcx,
            const vector< complex<double> >& dy_dcy,
            const vector< complex<double> >& Fx,
            const vector< complex<double> >& Fy,
            const vector<double>& eva
            ) {
        vector<double> n(2);

        if (rescaling_alg == "x") {
            // x direction
            n[0] = 1.0;
            n[1] = 0.0;
        }
        else if (rescaling_alg == "p") {
            // p direction
            n[0] = v[0];
            n[1] = v[1];
        }
        else if (rescaling_alg == "m1") {
            // Method #1 : real part of Berry force
            n[0] = (dcx[from+to*2] * (v[0] * dcx[to+from*2] + v[1] * dcy[to+from*2])).real();
            n[1] = (dcy[from+to*2] * (v[0] * dcx[to+from*2] + v[1] * dcy[to+from*2])).real();
        }
        else if (rescaling_alg == "m2") {
            // Method #2
            const vector<double> dcR { dcx[from+to*2].real(), dcy[from+to*2].real() };
            const vector<double> dcI { dcx[from+to*2].imag(), dcy[from+to*2].imag() };
            const double diff_norm2 = norm2(dcR) - norm2(dcI);
            const double twice_eta0 = std::atan(-2 * (dcR[0] * dcI[0] + dcR[1] * dcI[1]) / diff_norm2);
            double eta;
            if (cos(twice_eta0) * diff_norm2 > 0.0) {
                eta = 0.5 * twice_eta0;
            }
            else {
                eta = 0.5 * twice_eta0 + 0.5 * M_PI;
            }
            const complex<double> eieta = exp(matrixop::IMAGIZ * eta);
            n[0] = (eieta * dcx[from+to*2]).real();
            n[1] = (eieta * dcy[from+to*2]).real();
        }
        else if (rescaling_alg == "a1") {
            /*
             *  1->2:
             *  n = Re(g) = hbar / |d'*d| * Re[(\sum_beta d(d21^beta) / (dr^beta)) * d12)
             *
             */
            n[0] = ((dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcx[from+to*2]).real();
            n[1] = ((dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcy[from+to*2]).real();
        }
        else if (rescaling_alg == "a2") {
            /*
             *  1->2:
             *  n = Im(g) = hbar / |d'*d| * Re[(\sum_beta d(d21^beta) / (dr^beta)) * d12)
             *
             */
            n[0] = ((dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcx[from+to*2]).imag();
            n[1] = ((dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcy[from+to*2]).imag();
        }
        else if (rescaling_alg == "a3") {
            /*
             *  1->2:
             *  n = Re(g') = hbar / |d'*d| * Re[(\sum_beta conj(d(d21^beta) / (dr^beta))) * d12)
             *
             */
            n[0] = (conj(dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcx[from+to*2]).real();
            n[1] = (conj(dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcy[from+to*2]).real();
        }
        else if (rescaling_alg == "a4") {
            /*
             *  1->2:
             *  n = Im(conj(g)) = hbar / |d'*d| * Im[(\sum_beta conj(d(d21^beta) / (dr^beta))) * d12)
             *
             */
            n[0] = (conj(dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcx[from+to*2]).imag();
            n[1] = (conj(dx_dcx[to+from*2] + dy_dcy[to+from*2]) * dcy[from+to*2]).imag();
        }
        else {
            misc::crasher::confirm(false, "hopper: Invalid rescaling_alg");
        }

        return n;
    }

};

#endif
