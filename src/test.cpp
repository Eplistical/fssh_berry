#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/fmtstring.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/matrixop.hpp"

//#include "2d_flat_potential_phixplusy.hpp"
//#include "2d_tully1_potential.hpp"
#include "2d_helix_potential.hpp"
//#include "2d_marcus_potential.hpp"


double cal_crossing_E0(const double x) {
    // get lower adiabat energy for crossing line
    // crossig line: y = x
    const double y = x;
    return 0.5 * (param_A * pow(x, 2) + pow(y - param_B, 2) + param_D) - param_C;
}

int main(int argc, char** argv) {
    output_potential_param();

    const int N = 400;
    const double xmin = -0.1, xmax = 0.5;
    vector<double> xarr = linspace(xmin, xmax, N);
    vector<double> yarr = linspace(xmin, xmax, N);

    vector<double> eva(2);

    for (int i(0); i < N; ++i) {
        const double x = xarr[i];
        const double y = x;
        vector<double> r { x, y };

        auto H = cal_H(r);
        matrixop::eigh(H, eva);

        ioer::tabout(
                x, 
                H[0+0*2].real(), 
                0.5 * 0.1 * x * x + 0.5 * (x - 0.2) * (x - 0.2),
                eva[0], 
                cal_crossing_E0(x)
                //eva[1]
                );
    }
    return 0;
}
