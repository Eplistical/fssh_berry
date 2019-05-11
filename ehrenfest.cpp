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

double A = 0.10;
double B = 3.0;
double W = 0.3;
const double mass = 1000.0;
double init_x = -3.0;
double sigma_x = 0.5; 
double init_px = 30.0;
double sigma_px = 1.0; 
double init_y = 0.0;
double sigma_y = 0.5; 
double init_py = 0.0;
double sigma_py = 1.0; 
double init_s = 1.0;
double xwall_left = -10.0;
double xwall_right = 10.0;
int Nstep = 1000000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 5000;
int seed = 0;
bool enable_berry_force = true;
string output_mod = "init_px";

vector< complex<double> > Fx, Fy, dcx, dcy;

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("A", po::value<double>(&A), "potential para A")
        ("B", po::value<double>(&B), "potential para B")
        ("W", po::value<double>(&W), "potential W")
        ("init_x", po::value<double>(&init_x), "init x")
        ("sigma_x", po::value<double>(&sigma_x), "init sigma x")
        ("init_px", po::value<double>(&init_px), "potential para init_px")
        ("sigma_px", po::value<double>(&sigma_px), "init sigma px")
        ("init_y", po::value<double>(&init_y), "init y")
        ("sigma_y", po::value<double>(&sigma_y), "init sigma y")
        ("init_py", po::value<double>(&init_py), "potential para init_py")
        ("sigma_py", po::value<double>(&sigma_py), "init sigma py")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("xwall_left", po::value<double>(&xwall_left), "wall on x direction to check end")
        ("xwall_right", po::value<double>(&xwall_right), "wall on x direction to check end")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("seed", po::value<int>(&seed), "random seed")
        ("enable_berry_force", po::value<bool>(&enable_berry_force), "enable Berry force")
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

inline double cal_theta(double x) {
    return 0.5 * M_PI * (erf(B * x) + 1);
}

inline double cal_der_theta(double x) {
    return sqrt(M_PI) * B * exp(-B * B * x * x);
}

inline double cal_phi(double y) {
    return W * y;
}

inline double cal_der_phi(double y) {
    return W;
}

void cal_info(double x, double y)
{
    Fx.assign(4, 0.0);
    Fy.assign(4, 0.0);
    dcx.assign(4, 0.0);
    dcy.assign(4, 0.0);

    const double theta = cal_theta(x);
    const double der_theta = cal_der_theta(x);
    const double phi = cal_phi(y);
    const double der_phi = cal_der_phi(y);
    const double CC = cos(theta / 2);
    const double SS = sin(theta / 2);
    const double E1 = A, E0 = -A;

    dcx = vector< complex<double> > { 0.0, -0.5 * der_theta, 0.5 * der_theta, 0.0 };

    dcy = vector< complex<double> > { 
         matrixop::IMAGIZ * der_phi * CC * CC,
         matrixop::IMAGIZ * der_phi * SS * CC,
         matrixop::IMAGIZ * der_phi * SS * CC,
         matrixop::IMAGIZ * der_phi * SS * SS
    };

    Fx[0+0*2] = 0.0;
    Fx[0+1*2] = (E0 - E1) * dcx[0+1*2];
    Fx[1+0*2] = (E1 - E0) * dcx[1+0*2];
    Fx[1+1*2] = 0.0;

    Fy[0+0*2] = 0.0;
    Fy[0+1*2] = (E0 - E1) * dcy[0+1*2];
    Fy[1+0*2] = (E1 - E0) * dcy[1+0*2];
    Fy[1+1*2] = 0.0;
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
    // extra Berry Force
    double Fx0_berry = 0.0, Fx1_berry = 0.0, Fy0_berry = 0.0, Fy1_berry = 0.0;
    if (enable_berry_force) {
        Fx0_berry = 2 * (dcx[0+1*2] * (vx * dcx[1+0*2] + vy * dcy[1+0*2])).imag();
        Fx1_berry = 2 * (dcx[1+0*2] * (vx * dcx[0+1*2] + vy * dcy[0+1*2])).imag();
        Fy0_berry = 2 * (dcy[0+1*2] * (vx * dcx[1+0*2] + vy * dcy[1+0*2])).imag();
        Fy1_berry = 2 * (dcy[1+0*2] * (vx * dcx[0+1*2] + vy * dcy[0+1*2])).imag();
    }

    state_dot[0] = vx;
    state_dot[1] = vy;
    state_dot[2] = (a00 * (Fx[0+0*2] + Fx0_berry) + a01*Fx[1+0*2] + a10*Fx[0+1*2] + a11*(Fx[1+1*2] + Fx1_berry)) / mass;
    state_dot[3] = (a00 * (Fy[0+0*2] + Fy0_berry) + a01*Fy[1+0*2] + a10*Fy[0+1*2] + a11*(Fy[1+1*2] + Fy1_berry)) / mass;
    state_dot[4] = -a10 * vd01 + a01 * vd10;
    state_dot[5] = -zI * -2.0 * A * a01 + (a00 - a11) * vd01;
    state_dot[6] = -zI *  2.0 * A * a10 + (a11 - a00) * vd10;
    state_dot[7] = a10 * vd01 - a01 * vd10;
}

state_t init_state() {
    // state = x, y, vx, vy, a00, a01, a10, a11, s
    state_t state(8, zzero);
    state[0].real(randomer::normal(-3.0, 0.5)); 
    state[1].real(randomer::normal(-2.5, 0.5)); 
    state[2].real(randomer::normal(30.0, 1.0) / mass); 
    state[3].real(randomer::normal( 0.0, 1.0) / mass); 

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
    return ((x > 5.0 and vx > 0.0) or (x < -5.0 and vx < 0.0));
}

void erenfest() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        state[itraj] = init_state();
    }

    // para
    ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, 
                " mass = ", mass, " A = ", A, " B = ", B, " W = ", W,
                " init_x = ", init_x, " init_px = ", init_px, 
                " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                " init_y = ", init_y, " init_py = ", init_py, 
                " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                " init_s = ", init_s
            );

    // main loop
    double vx, vy, v2x, v2y, Ep;
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            // calc info
            cal_info(state[itraj][0].real(), state[itraj][1].real());
            // propagate
            rk4.do_step(sys, state[itraj], istep * dt, dt);
        }

        // output
        if (istep == 0) {
            ioer::tabout(
                    "#",
                    "t",
                    "px", "py",
                    "Ekx", "Eky",
                    "Ep", "Etot",
                    "");
        }

        if (istep % output_step == 0) {
            vx = vy = 0.0;
            v2x = v2y = 0.0;
            Ep = 0.0;

            for_each(state.begin(), state.end(), [&vx](const state_t& st) { vx += st[2].real(); });
            for_each(state.begin(), state.end(), [&vy](const state_t& st) { vy += st[3].real(); });
            for_each(state.begin(), state.end(), [&v2x](const state_t& st) { v2x += pow(st[2].real(), 2); });
            for_each(state.begin(), state.end(), [&v2y](const state_t& st) { v2y += pow(st[3].real(), 2); });
            for_each(state.begin(), state.end(), [&Ep](const state_t& st) { 
                        Ep += st[4].real() * -A + st[7].real() * A;
                    });

            ioer::tabout(
                    "#",
                    istep * dt, 
                    mass * vx / Ntraj, 
                    mass * vy / Ntraj,
                    0.5 * mass * v2x / Ntraj,
                    0.5 * mass * v2y / Ntraj,
                    Ep / Ntraj,
                    (0.5 * mass * (v2x + v2y) + Ep) / Ntraj,
                    "");
            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                break;
            }
        }
    }
    ioer::tabout(
            //init_px,
            init_s,
            mass * vx / Ntraj, 
            mass * vy / Ntraj,
            0.5 * mass * v2x / Ntraj,
            0.5 * mass * v2y / Ntraj,
            Ep / Ntraj,
            (0.5 * mass * (v2x + v2y) + Ep) / Ntraj,
            "");
}

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    erenfest();
    ioer::info("# ", timer::toc());
    return 0;
}
