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

int main(int argc, char** argv) {
    output_potential_param();

    const int N = 400;
    const double xmin = -10.0, xmax = 10.0;
    vector<double> xarr = linspace(xmin, xmax, N);
    vector<double> yarr = linspace(xmin, xmax, N);

    vector<double> eva(2);

    ioer::tabout("# x", "y", "E0", "E1", "H00", "H11");

    for (int i(0); i < N; ++i) {
        const double x = xarr[i];
        for (int j(0); j < N; ++j) {
            const double y = yarr[j];
            vector<double> r { x, y };

            auto H = cal_H(r);
            matrixop::eigh(H, eva);

            ioer::tabout(x, y, eva[0], eva[1], H[0+0*2].real(), H[1+1*2].real());
        }
    }
    return 0;
}
