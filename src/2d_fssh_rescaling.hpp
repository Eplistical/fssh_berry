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
             * additional ansatz 1
             *  1->2:
             *  hbar/(E2-E1) * Re[d12 * (d21 \dot (hbar * p / m + (F22 - F11) / (||d12*d21||)))]
             *  ~ Re[d12 * (d21 \dot (hbar * p / m + (F22 - F11)/ (||d12*d21||)))]
             */
            const double normd = sqrt(pow(abs(dcx[to+from*2]), 2) + pow(abs(dcy[to+from*2]), 2));
            const double DeltaFx = Fx[to+to*2].real() - Fx[from+from*2].real();
            const double DeltaFy = Fy[to+to*2].real() - Fy[from+from*2].real();
            n[0] = (dcx[to+from*2] * (dcx[from+to*2] * (v[0] + DeltaFx / normd) + dcy[from+to*2] * (v[1] + DeltaFy / normd))).real();
            n[1] = (dcy[to+from*2] * (dcx[from+to*2] * (v[0] + DeltaFx / normd) + dcy[from+to*2] * (v[1] + DeltaFy / normd))).real();
            
        }
        else if (rescaling_alg == "a2") {
            /*
             * additional ansatz 2
             *  1->2:
             *  hbar/(E2-E1) * Re[d12 * (d21 \dot (hbar * p / m - (F22 - F11)/ (||d12*d21||)))]
             *  ~ Re[d12 * (d21 \dot (hbar * p / m - (F22 - F11)/ (||d12*d21||)))]
             */
            const double normd = sqrt(pow(abs(dcx[to+from*2]), 2) + pow(abs(dcy[to+from*2]), 2));
            const double DeltaFx = Fx[to+to*2].real() - Fx[from+from*2].real();
            const double DeltaFy = Fy[to+to*2].real() - Fy[from+from*2].real();
            n[0] = (dcx[to+from*2] * (dcx[from+to*2] * (v[0] - DeltaFx / normd) + dcy[from+to*2] * (v[1] - DeltaFy / normd))).real();
            n[1] = (dcy[to+from*2] * (dcx[from+to*2] * (v[0] - DeltaFx / normd) + dcy[from+to*2] * (v[1] - DeltaFy / normd))).real();
        }
        else {
            misc::crasher::confirm(false, "hopper: Invalid rescaling_alg");
        }

        return n;
    }

};

#endif
