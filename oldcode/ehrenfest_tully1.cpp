#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;
using state_t = vector< complex<double> >;

// global const
const complex<double> zzero(0.0, 0.0);
const complex<double> zone(1.0, 0.0);
const complex<double> zI(0.0, 1.0);

double A = 0.01;
double B = 1.6;
double C = 0.005;
double D = 1.0;
double W = 0.0;
const double mass = 2000.0;
double init_x = -5.0;
double sigma_x = 0.5; 
double init_px = 20.0;
double sigma_px = 1.0; 
double init_y = 0.0;
double sigma_y = 0.5; 
double init_py = 0.0;
double sigma_py = 1.0; 
double init_s = 0.0;
double xwall_left = -10.0;
double xwall_right = 10.0;
int Nstep = 1000000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 2000;
int seed = 0;
string output_mod = "init_px";

vector< complex<double> > lastevt;
vector<double> eva(2);
vector< complex<double> > Fx(4), Fy(4);
vector< complex<double> > dcx(4), dcy(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("A", po::value<double>(&A), "potential para A")
        ("B", po::value<double>(&B), "potential para B")
        ("C", po::value<double>(&C), "potential para C")
        ("D", po::value<double>(&D), "potential para D")
        ("W", po::value<double>(&W), "potential W")
        ("init_x", po::value<double>(&init_x), "init x")
        ("sigma_x", po::value<double>(&sigma_x), "init sigma x")
        ("init_px", po::value<double>(&init_px), "init px")
        ("sigma_px", po::value<double>(&sigma_px), "init sigma px")
        ("init_y", po::value<double>(&init_y), "init y")
        ("sigma_y", po::value<double>(&sigma_y), "init sigma y")
        ("init_py", po::value<double>(&init_py), "init py")
        ("sigma_py", po::value<double>(&sigma_py), "init sigma py")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("xwall_left", po::value<double>(&xwall_left), "wall on x direction to check end")
        ("xwall_right", po::value<double>(&xwall_right), "wall on x direction to check end")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("seed", po::value<int>(&seed), "random seed")
        ("output_mod", po::value<string>(&output_mod), "output mode, init_s or init_px")
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

inline double cal_phi(double y) {
    return W * y;
}

inline double cal_der_phi(double y) {
    return W;
}

vector< complex<double> > cal_H(const vector<double>& r) {
    const double x = r[0];
    const double y = r[1];
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(y));

    vector< complex<double> > H(4);
    if (x >= 0.0 ) {
        H[0+0*2] = A * (1.0 - exp(-B * x));
    }
    else {
        H[0+0*2] = -A * (1.0 - exp(B * x));
    }
    H[1+1*2] = -H[0+0*2];
    H[0+1*2] = C * exp(-D * x * x) * eip;
    H[1+0*2] = conj(H[0+1*2]);
    return H;
}

vector< complex<double> > cal_nablaHx(const vector<double>& r) {
    const double x = r[0];
    const double y = r[1];
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(y));
    vector< complex<double> > nablaHx(4);

    if (x >= 0.0 ) {
        nablaHx[0+0*2] = A * B * exp(-B * x);
    }
    else {
        nablaHx[0+0*2] = A * B * exp(B * x);
    }

    nablaHx[1+1*2] = -nablaHx[0+0*2];
    nablaHx[0+1*2] = -2 * C * D * x * exp(-D * x * x) * eip;
    nablaHx[1+0*2] = conj(nablaHx[0+1*2]);

    return nablaHx;
}

vector< complex<double> > cal_nablaHy(const vector<double>& r) {
    const double x = r[0];
    const double y = r[1];
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(y));
    const double der_phi = cal_der_phi(y);
    vector< complex<double> > nablaHy(4);

    nablaHy[0+0*2] = 0.0;
    nablaHy[1+1*2] = -nablaHy[0+0*2];
    nablaHy[0+1*2] = C * exp(-D * x * x) * eip * matrixop::IMAGIZ * der_phi;
    nablaHy[1+0*2] = conj(nablaHy[0+1*2]);

    return nablaHy;
}

void cal_info_nume(const vector<double>& r)
{
    // nume
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
    Fx.assign(4, matrixop::ZEROZ);
    Fy.assign(4, matrixop::ZEROZ);
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

void cal_info(const vector<double>& r)
{
    cal_info_nume(r);
}


void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, y, vx, vy, a00, a01, a10, a11
    double x = state[0].real();
    double y = state[1].real();
    double vx = state[2].real();
    double vy = state[3].real();
    complex<double> a00 = state[4];
    complex<double> a01 = state[5];
    complex<double> a10 = state[6];
    complex<double> a11 = state[7];

    complex<double> vd01 = vx * dcx[0+1*2] + vy * dcy[0+1*2];
    complex<double> vd10 = vx * dcx[1+0*2] + vy * dcy[1+0*2];

    state_dot[0] = vx;
    state_dot[1] = vy;
    state_dot[2] = (a00 * Fx[0+0*2] + a01*Fx[1+0*2] + a10*Fx[0+1*2] + a11*Fx[1+1*2]) / mass;
    state_dot[3] = (a00 * Fy[0+0*2] + a01*Fy[1+0*2] + a10*Fy[0+1*2] + a11*Fy[1+1*2]) / mass;
    state_dot[4] = -a10 * vd01 + a01 * vd10;
    state_dot[5] = -zI * 2.0 * eva[0] * a01 + (a00 - a11) * vd01;
    state_dot[6] = -zI * 2.0 * eva[1] * a10 + (a11 - a00) * vd10;
    state_dot[7] = a10 * vd01 - a01 * vd10;
}

state_t init_state() {
    // state = x, y, vx, vy, a00, a01, a10, a11
    state_t state(8, zzero);
    state[0].real(randomer::normal( init_x, sigma_x)); 
    state[1].real(randomer::normal( init_y, sigma_y)); 
    state[2].real(randomer::normal( init_px, sigma_px) / mass); 
    state[3].real(randomer::normal( init_py, sigma_py) / mass); 

    complex<double> c0 = sqrt(1.0 - init_s);
    complex<double> c1 = sqrt(init_s);

    state[4] = c0 * conj(c0);
    state[5] = c1 * conj(c0);
    state[6] = c0 * conj(c1);
    state[7] = c1 * conj(c1);
    return state;
}

bool check_end(const state_t& state) {
    double x = state[0].real();
    double vx = state[2].real();
    return ((x > xwall_right and vx > 0.0) or (x < xwall_left and vx < 0.0));
}

void ehrenfest() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        state[itraj] = init_state();
    }

    // para
    ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, 
                " mass = ", mass, 
                " A = ", A, " B = ", B, 
                " C = ", C, " D = ", D, 
                " W = ", W,
                " init_x = ", init_x, " init_px = ", init_px, 
                " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                " init_y = ", init_y, " init_py = ", init_py, 
                " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                " init_s = ", init_s
            );

    // main loop
    double ntrans, nrefl;
    double n1trans, n1refl;
    double n0trans, n0refl;
    double vxtrans, vytrans, vxrefl, vyrefl;
    double Ep, Ek;
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            // assign last evt
            lastevt = move(lastevt_save[itraj]);
            // calc info
            cal_info(vector<double> { state[itraj][0].real(), state[itraj][1].real() } );
            // propagate
            rk4.do_step(sys, state[itraj], istep * dt, dt);
            // save last evt
            lastevt_save[itraj] = move(lastevt);
        }

        // output
        if (istep == 0) {
            ioer::tabout(
                    "#",
                    "t",
                    "n0trans", "n0refl",
                    "n1trans", "n1refl",
                    "ntrans", "nrefl",
                    "pxtrans", "pytrans",
                    "pxrefl", "pyrefl",
                    "Etot",
                    "");
        }

        if (istep % output_step == 0) {
            ntrans = nrefl = 0.0;
            n0trans = n0refl = 0.0;
            n1trans = n1refl = 0.0;
            vxtrans = vytrans = 0.0;
            vxrefl = vyrefl = 0.0;
            Ep = Ek = 0.0;
    // state = x, y, vx, vy, a00, a01, a10, a11
            for_each(state.begin(), state.end(), [&ntrans, &nrefl, &n0trans, &n0refl, &n1trans, &n1refl, &vxtrans, &vytrans, &vxrefl, &vyrefl, &Ep, &Ek](const state_t& st) {
                        if (st[2].real() > 0.0) {
                            ntrans += 1.0;
                            n0trans += st[4].real();
                            n1trans += st[7].real();

                            vxtrans += st[2].real();
                            vytrans += st[3].real();
                        }
                        else {
                            nrefl += 1.0;
                            n0refl += st[4].real();
                            n1refl += st[7].real();

                            vxrefl += st[2].real();
                            vyrefl += st[3].real();
                        }
                        Ek += 0.5 * mass * (pow(st[2].real(), 2) + pow(st[3].real(), 2));
                        matrixop::hdiag(cal_H(vector<double> { st[0].real(), st[1].real() }), eva);
                        Ep += st[4].real() * eva[0] + st[7].real() * eva[1];
                    });


            ioer::tabout(
                    "#",
                    istep * dt, 
                    n0trans / Ntraj,
                    n0refl / Ntraj,
                    n1trans / Ntraj,
                    n1refl / Ntraj,
                    ntrans / Ntraj,
                    nrefl / Ntraj,
                    ntrans > 0 ? mass * vxtrans / ntrans : 0.0, 
                    ntrans > 0 ? mass * vytrans / ntrans : 0.0, 
                    nrefl > 0 ? mass * vxrefl / nrefl : 0.0, 
                    nrefl > 0 ? mass * vyrefl / nrefl : 0.0, 
                    (Ek + Ep) / Ntraj, 
                    "");
            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                break;
            }
        }
    }
    if (output_mod == "init_px") {
        ioer::tabout_nonewline(init_px);
    }
    else if (output_mod == "init_s") {
        ioer::tabout_nonewline(init_s);
    }

    ioer::tabout(
            n0trans / Ntraj,
            n0refl / Ntraj,
            n1trans / Ntraj,
            n1refl / Ntraj,
            ntrans / Ntraj,
            nrefl / Ntraj,
            ntrans > 0 ? mass * vxtrans / ntrans : 0.0, 
            ntrans > 0 ? mass * vytrans / ntrans : 0.0, 
            nrefl > 0 ? mass * vxrefl / nrefl : 0.0, 
            nrefl > 0 ? mass * vyrefl / nrefl : 0.0, 
            (Ek + Ep) / Ntraj, 
            "");
}

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    ehrenfest();
    ioer::info("# ", timer::toc());
    return 0;
}
