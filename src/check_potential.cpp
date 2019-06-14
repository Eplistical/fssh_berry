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
#include "2d_tully1_potential.hpp"

int main(int argc, char** argv) {
    output_potential_param();

    const int N = 200;
    vector<double> xarr = linspace(-8.0, 8.0, N);
    vector<double> yarr = linspace(-8.0, 8.0, N);

    vector<double> eva(2);

    ioer::tabout("# x", "y", "E0", "E1");

    for (int i(0); i < N; ++i) {
        const double x = xarr[i];
        for (int j(0); j < N; ++j) {
            const double y = yarr[j];
            vector<double> r { x, y };

            auto H = cal_H(r);
            matrixop::eigh(H, eva);

            ioer::tabout(x, y, eva[0], eva[1]);
        }
    }
    return 0;
}
