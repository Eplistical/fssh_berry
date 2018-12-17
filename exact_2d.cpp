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
#include "misc/fft.hpp"
#include "misc/matrixop.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;

double L = 32;
int M = 512;
double mass = 1000.0;

double A = 0.1;
double B = 3.0;
double W = 0.3;

double xI = -3.0;
double yI = 0.0;
double sigmax = 1.0;
double sigmay = 1.0;
double kxI = 30.0;
double kyI = 0.0;

double xwall_left = -99999999.9;
double xwall_right = 99999999.9;

double init_s = 1.0;
int Nstep = 2000;
double dt = 0.1;

int output_step = 100;

string output_mod = "init_s";

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("L", po::value<double>(&L), "grid para, the grid is [-L/2, L/2]")
        ("M", po::value<int>(&M), "grid number")
        ("mass", po::value<double>(&mass), "mass")
        ("A", po::value<double>(&A), "potential para")
        ("B", po::value<double>(&B), "potential para")
        ("W", po::value<double>(&W), "potential para")
        ("init_x", po::value<double>(&xI), "init x")
        ("init_y", po::value<double>(&yI), "init y")
        ("sigma_x", po::value<double>(&sigmax), "init sigma x")
        ("sigma_y", po::value<double>(&sigmay), "init sigma y")
        ("init_px", po::value<double>(&kxI), "init px")
        ("init_py", po::value<double>(&kyI), "init py")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("xwall_left", po::value<double>(&xwall_left), " left boundary x value to check end")
        ("xwall_right", po::value<double>(&xwall_right), " right boundary x value to check end")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("dt", po::value<double>(&dt), "time step")
        ("output_step", po::value<int>(&output_step), "output step")
        ("output_mod", po::value<string>(&output_mod), "output mode, can be init_s or init_px")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

vector< complex<double> > cal_H(const double x, const double y) {
    vector< complex<double> > H(4);
    double theta = 0.5 * M_PI * (erf(B * x) + 1.0);
    double phi = W * y;
    complex<double> eip = exp(matrixop::IMAGIZ * phi);
    H[0+0*2] = -cos(theta);
    H[0+1*2] = sin(theta) * eip;
    H[1+0*2] = conj(H[0+1*2]);
    H[1+1*2] = cos(theta);
    return A * H;
}

vector< complex<double> > myfftshift(const vector< complex<double> >& Z) 
{
    vector< complex<double> > rst(M * M);
    for (int k = 0; k < M/2; ++k) {
        for (int j = 0; j < M/2; ++j) {
            rst[k+j*M] = Z[(k+M/2) + (j+M/2)*M];
            rst[(k+M/2)+(j+M/2)*M] = Z[k + j*M];
            rst[(k+M/2)+j*M] = Z[k + (j+M/2)*M];
            rst[k+(j+M/2)*M] = Z[(k+M/2) + j*M];
        }
    }
    return rst;
}

bool check_end( const vector< complex<double> >& psi0, const vector< complex<double> >& psi1, 
                const vector<double>& xarr, const vector<double>& yarr)
{
    // if the WF is outside the wall, return true to stop the program
    double n0_outside = 0.0;
    double n1_outside = 0.0;
    double n0_total = 0.0;
    double n1_total = 0.0;
    for (int k = 0; k < M; ++k) {
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            double y = yarr[k];
            n0_total += pow(abs(psi0[k+j*M]), 2);
            n1_total += pow(abs(psi1[k+j*M]), 2);
            if (x < xwall_left or x > xwall_right or y < xwall_left or y > xwall_right) {
                n0_outside += pow(abs(psi0[k+j*M]), 2);
                n1_outside += pow(abs(psi1[k+j*M]), 2);
            }
        }
    }
    if ((n0_outside / n0_total > 0.01) or (n1_outside / n1_total > 0.01)) {
        return true;
    }
    else {
        return false;
    }
}

void exact() {
    // para
    double dx, dy, dkx, dky;
    vector<double> xarr = linspace(-L/2, L/2, M, dx);
    vector<double> yarr = linspace(-L/2, L/2, M, dy);
    dkx = 2 * M_PI / M / dx;
    dky = 2 * M_PI / M / dy;
    vector<double> kxarr = (arange(M) - M / 2) * dkx;
    vector<double> kyarr = (arange(M) - M / 2) * dky;
    xwall_left = max(xwall_left, -L/2 * 0.9);
    xwall_right = min(xwall_right, L/2 * 0.9);
    double c0 = sqrt(1 - init_s);
    double c1 = sqrt(init_s);
    // construct TU on k grid
    vector< complex<double> > TU(M * M);
    for (int k = 0; k < M; ++k) {
        double ky = kyarr[k];
        for (int j = 0; j < M; ++j) {
            double kx = kxarr[j];
            TU[k+j*M] = exp(-matrixop::IMAGIZ * dt * (kx*kx + ky*ky) / 2.0 / mass);
        }
    }
    TU = myfftshift(TU);
    // construct VU on x grid
    vector< complex<double> > V00(M*M), V01(M*M), V10(M*M), V11(M*M);
    vector< complex<double> > H00(M*M), H01(M*M), H10(M*M), H11(M*M);
    for (int k = 0; k < M; ++k) {
        double y = yarr[k];
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            vector< complex<double> > H = cal_H(x, y);
            vector<double> eva;
            vector< complex<double> > evt, evamat;
            matrixop::eigh(H, eva, evt);

            evamat.assign(4, 0.0);
            evamat[0+0*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[0]);
            evamat[1+1*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[1]);
            auto tmp = matrixop::matmatmatC(evt, evamat, evt, 2, 2);
            V00[k+j*M] = tmp[0+0*2];
            V01[k+j*M] = tmp[0+1*2];
            V10[k+j*M] = tmp[1+0*2];
            V11[k+j*M] = tmp[1+1*2];
            H00[k+j*M] = H[0+0*2];
            H01[k+j*M] = H[0+1*2];
            H10[k+j*M] = H[1+0*2];
            H11[k+j*M] = H[1+1*2];
        }
    }
    // initialized WF
    vector< complex<double> > psi0(M * M, 0.0), psi1(M * M, 0.0);
    for (int k = 0; k < M; ++k) {
        double y = yarr[k];
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            psi0[k+j*M] = c0 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
            psi1[k+j*M] = c1 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
        }
    }
    double nm = norm(psi0 | psi1);
    psi0 /= nm;
    psi1 /= nm;
    // covinience vairables
    vector<int> dim{ M, M };
    // statistics
    double KE0 = 0.0, KE1 = 0.0, PE = 0.0;
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double py0trans = 0.0, py0refl = 0.0, py1trans = 0.0, py1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    // propagate WF
    for (int istep = 0; istep < Nstep; ++istep) {
        // exp(-iVdt/2)
        auto psi_k0 = V00 * psi0 + V01 * psi1;
        auto psi_k1 = V10 * psi0 + V11 * psi1;
        // exp(-iTdt)
        psi_k0 = misc::fftn(psi_k0, dim);
        psi_k1 = misc::fftn(psi_k1, dim);
        psi_k0 *= TU;
        psi_k1 *= TU;
        // exp(-iVdt/2)
        psi_k0 = misc::ifftn(psi_k0, dim);
        psi_k1 = misc::ifftn(psi_k1, dim);
        psi0 = V00 * psi_k0 + V01 * psi_k1;
        psi1 = V10 * psi_k0 + V11 * psi_k1;
        // analysis & output
        if (istep % output_step == 0) {
            // get psi_k
            psi_k0 = myfftshift(misc::fftn(psi0, dim));
            psi_k1 = myfftshift(misc::fftn(psi1, dim));
            double nm = norm(psi_k0 | psi_k1);
            psi_k0 /= nm;
            psi_k1 /= nm;
            KE0 = KE1 = PE = 0.0;
            n0trans = n0refl = n1trans = n1refl = 0.0;
            py0trans = py0refl = py1trans = py1refl = 0.0;
            px0trans = px0refl = px1trans = px1refl = 0.0;
            for (int k = 0; k < M; ++k) {
                for (int j = 0; j < M; ++j) {
                    if (kxarr[j] >= 0.0) {
                        n0trans += pow(abs(psi_k0[k+j*M]), 2);
                        n1trans += pow(abs(psi_k1[k+j*M]), 2);
                        py0trans += pow(abs(psi_k0[k+j*M]), 2) * kyarr[k];
                        py1trans += pow(abs(psi_k1[k+j*M]), 2) * kyarr[k];
                        px0trans += pow(abs(psi_k0[k+j*M]), 2) * kxarr[j];
                        px1trans += pow(abs(psi_k1[k+j*M]), 2) * kxarr[j];
                    }
                    else {
                        n0refl += pow(abs(psi_k0[k+j*M]), 2);
                        n1refl += pow(abs(psi_k1[k+j*M]), 2);
                        py0refl += pow(abs(psi_k0[k+j*M]), 2) * kyarr[k];
                        py1refl += pow(abs(psi_k1[k+j*M]), 2) * kyarr[k];
                        px0refl += pow(abs(psi_k0[k+j*M]), 2) * kxarr[j];
                        px1refl += pow(abs(psi_k1[k+j*M]), 2) * kxarr[j];
                    }
                    KE0 += pow(abs(psi_k0[k+j*M]), 2) * (kxarr[j]*kxarr[j] + kyarr[k]*kyarr[k]) / 2.0 / mass;
                    KE1 += pow(abs(psi_k1[k+j*M]), 2) * (kxarr[j]*kxarr[j] + kyarr[k]*kyarr[k]) / 2.0 / mass;
                }
            }
            py0trans /= (n0trans + 1e-16);
            py1trans /= (n1trans + 1e-16);
            py0refl /= (n0refl + 1e-16);
            py1refl /= (n1refl + 1e-16);
            px0trans /= (n0trans + 1e-16);
            px1trans /= (n1trans + 1e-16);
            px0refl /= (n0refl + 1e-16);
            px1refl /= (n1refl + 1e-16);
            PE = real(sum(conj(psi0) * H00 * psi0 + conj(psi0) * H01 * psi1 + conj(psi1) * H10 * psi0 + conj(psi1) * H11 * psi1));
            // output
            if (istep == 0) {
                ioer::info("# EXACT 2D ");
                ioer::info("# para: ", " L = ", L, " M = ", M, " mass = ", mass, " A = ", A, " B = ", B, " W = ", W, 
                                       " xI = ", xI, " yI = ", yI, " sigmax = ", sigmax, " sigmay = ", sigmay, " kxI = ", kxI, " kyI = ", kyI, " init_s = ", init_s, " c0 = ", c0, " c1 = ", c1,
                                       " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step);
                ioer::info("# dx = ", dx, " dy = ", dy, " dkx = ", dkx, " dky = ", dky, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right);
                ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "py0trans", "px0refl", "py0refl", "px1trans", "py1trans", "px1refl", "py1refl", "Etot");
            }
            ioer::tabout('#', istep * dt, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl, KE0 + KE1 + PE);
            // check end
            if (check_end(psi0, psi1, xarr, yarr) == true) {
                ioer::info("# check_end returns true");
                break;
            }
        }
    }
    // final output
    if (output_mod == "init_s") {
        ioer::tabout(init_s, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl);
    }
    else if (output_mod == "init_px") {
        ioer::tabout(kxI, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl);
    }
    else {
    }
}

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    exact();
    ioer::info("# ", timer::toc());
    return 0;
}
