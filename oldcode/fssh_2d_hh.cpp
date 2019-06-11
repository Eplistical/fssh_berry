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
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};

// BE CAREFUL! PHASE IS IMPORTANT HERE! (DC AND BERRY FORCE)

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;
using state_t = vector< complex<double> >;
const complex<double> zI(0.0, 1.0);
double A = 0.10;
double B = 3.0;
double W = 0.5;
const double mass = 1000.0;
double init_x = -3.0;
double sigma_x = 0.5; 
double init_px = 20.0;
double sigma_px = 1.0; 
double init_y = 0.0;
double sigma_y = 0.5; 
double init_py = 0.0;
double sigma_py = 1.0; 
double init_s = 1.0;
double xwall_left = -5.0;
double xwall_right = 5.0;
int Nstep = 1000000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 2000;
bool enable_hop = true;

vector< complex<double> > lastevt;
vector<double> eva(2);
vector<double> Fx(2), Fy(2);
vector< complex<double> > dcx(4), dcy(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("A", po::value<double>(&A), "potential para A")
        ("B", po::value<double>(&B), "potential para B")
        ("W", po::value<double>(&W), "potential phase para W")
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
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
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

double cal_theta(const double x) {
    return 0.5 * M_PI * (erf(B * x) + 1);
}

double cal_der_theta(const double x) {
    return sqrt(M_PI) * B * exp(-B * B * x * x);
}

double cal_phi(const double y) {
    return W * y;
}

double cal_der_phi(const double y) {
    return W;
}

vector< complex<double> > cal_H(const vector<double>& r) {
    const double theta = cal_theta(r[0]);
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
    vector< complex<double> > H {
        -cos(theta), sin(theta) * conj(eip), sin(theta) * eip, cos(theta)
    };
    return A * H;
}

vector< complex<double> > cal_nablaHx(const vector<double>& r) {
    const double theta = cal_theta(r[0]);
    const double der_theta = cal_der_theta(r[0]);
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
    vector< complex<double> > nablaHx {
        sin(theta), cos(theta) * conj(eip), cos(theta) * eip, -sin(theta)
    };
    return A * der_theta * nablaHx;
}

vector< complex<double> > cal_nablaHy(const vector<double>& r) {
    const double theta = cal_theta(r[0]);
    const double der_theta = cal_der_theta(r[0]);
    const double der_phi = cal_der_phi(r[1]);
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(r[1]));
    vector< complex<double> > nablaHy {
        0.0, -sin(theta) * conj(eip), sin(theta) * eip, 0.0
    };
    return matrixop::IMAGIZ * A * der_phi * nablaHy;
}

void init_state(state_t& state) {
    state.resize(7, matrixop::ZEROZ);
    // state = x, y, vx, vy, c0, c1, s
    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_y, sigma_y)); 
    state[2].real(randomer::normal(init_px, sigma_px) / mass); 
    state[3].real(randomer::normal(init_py, sigma_py) / mass); 
    state[4].real(sqrt(1.0 - init_s));
    state[5].real(sqrt(init_s));
    state[6].real((randomer::rand() < init_s) ? 1.0 : 0.0);

    misc::crasher::confirm(state[2].real() > 0.0, "px < 0.0 !!");
}

void cal_info_anal(const vector<double>& r)
{
    double theta = cal_theta(r[0]);
    double der_theta = cal_der_theta(r[0]);
    double der_phi= cal_der_phi(r[1]);
    double CC = cos(0.5 * theta);
    double SS = sin(0.5 * theta);

    eva = vector<double> { -A, A };

    Fx = vector<double> { 0.0, 0.0 };
    Fy = vector<double> { 0.0, 0.0 };

    // d00, d10, d01, d11
    dcx = vector< complex<double> > { 0.0, -0.5 * der_theta, 0.5 * der_theta, 0.0 };
    dcy = vector< complex<double> > { 
         matrixop::IMAGIZ * der_phi * CC * CC,
         matrixop::IMAGIZ * der_phi * SS * CC,
         matrixop::IMAGIZ * der_phi * SS * CC,
         matrixop::IMAGIZ * der_phi * SS * SS,
    };
}

void cal_info(const vector<double>& r)
{
    cal_info_anal(r);
}

void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, v, c0, c1, s
    vector<double> r { state[0].real(), state[1].real() };
    vector<double> v { state[2].real(), state[3].real() };
    vector< complex<double> > c { state[4], state[5] };
    int s = static_cast<int>(state[6].real());
    // extra Berry Force
    double Fx_berry, Fy_berry;
    Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    // state_dot
    state_dot[0] = v[0];
    state_dot[1] = v[1];
    state_dot[2] = (Fx[s] + Fx_berry) / mass;
    state_dot[3] = (Fy[s] + Fy_berry) / mass;
    state_dot[4] = -zI * c[0] * eva[0] - c[1] * (v[0] * dcx[0+1*2] + v[1] * dcy[0+1*2]) - c[0] * (v[0] * dcx[0+0*2] + v[1] * dcy[0+0*2]);
    state_dot[5] = -zI * c[1] * eva[1] - c[0] * (v[0] * dcx[1+0*2] + v[1] * dcy[1+0*2]) - c[1] * (v[0] * dcx[1+1*2] + v[1] * dcy[1+1*2]);
    state_dot[6] = matrixop::ZEROZ;
}

int hopper(state_t& state) {
    // state = x, v, c0, c1, s
    vector<double> r { state[0].real(), state[1].real() };
    vector<double> v { state[2].real(), state[3].real() };
    vector< complex<double> > c { state[4], state[5] };
    int s = static_cast<int>(state[6].real());
    // calc hop prob
    double g = -2 * dt * (c[s] * conj(c[1-s]) * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).real() / (c[s] * conj(c[s])).real();
    double dE = eva[1-s] - eva[s];
    // hop
    if (randomer::rand() < g) {
        // adjust momentum
        vector<double> n(2);

        // x-direction
        n[0] = 1.0;
        n[1] = 0.0;

        if (norm(n) > 1e-40) {
            vector<double> vn = component(v, n);
            double vn_norm(norm(vn)); 
            double tmp = vn_norm * vn_norm - 2 * dE / mass;
            if (tmp > 0.0) {
                // hop accepted
                double vn_norm_new = sqrt(tmp);
                v += (vn_norm_new - vn_norm) / vn_norm * vn;
                state[2].real(v[0]);
                state[3].real(v[1]);
                state[6].real(1.0 - s);
                return (s == 0) ? HOP_UP : HOP_DN;
            }
            else {
                return HOP_FR;
            }
        }
    }
    return HOP_RJ;
}

bool check_end(const state_t& state) {
    double x = state[0].real();
    double vx = state[2].real();
    return ((x > xwall_right and vx > 0.0) or (x < xwall_left and vx < 0.0));
}

void fssh() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        init_state(state[itraj]);
    }
    // main loop
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    vector<double> hop_count(Ntraj, 0.0);
    // statistics
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    double py0trans = 0.0, py0refl = 0.0, py1trans = 0.0, py1refl = 0.0;
    double KE = 0.0, PE = 0.0;
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info(vector<double> { state[itraj][0].real(), state[itraj][1].real() });
                // hopper
                if (enable_hop) {
                    int hopflag = hopper(state[itraj]);
                    switch (hopflag) {
                        case HOP_UP : { hopup += 1.0; hop_count[itraj] += 1.0; break; }
                        case HOP_DN : { hopdn += 1.0; hop_count[itraj] += 1.0; break; }
                        case HOP_FR : hopfr += 1.0; break;
                        case HOP_RJ : hoprj += 1.0; break;
                        default : break;
                    }
                }
                // propagate
                rk4.do_step(sys, state[itraj], istep * dt, dt);
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
            // population
            n0trans = n0refl = n1trans = n1refl = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0trans, &n0refl, &n1trans, &n1refl](const state_t& st) { 
                        int s = static_cast<int>(st[6].real());
                        double vx = st[2].real();
                        if (s == 0) {
                            (vx >= 0.0) ? n0trans += 1.0 : n0refl += 1.0;
                        }
                        else {
                            (vx >= 0.0) ? n1trans += 1.0 : n1refl += 1.0;
                        }
                    });
            // momentum
            px0trans = px0refl = px1trans = px1refl = 0.0;
            py0trans = py0refl = py1trans = py1refl = 0.0;
            for_each(state.begin(), state.end(), 
                    [&px0trans, &px0refl, &px1trans, &px1refl, &py0trans, &py0refl, &py1trans, &py1refl] (const state_t& st) {
                        int s = static_cast<int>(st[6].real());
                        double vx = st[2].real();
                        if (s == 0) {
                            if (vx >= 0.0) {
                                px0trans += mass * st[2].real();
                                py0trans += mass * st[3].real();
                            }
                            else {
                                px0refl += mass * st[2].real();
                                py0refl += mass * st[3].real();
                            }
                        }
                        else {
                            if (vx >= 0.0) {
                                px1trans += mass * st[2].real();
                                py1trans += mass * st[3].real();
                            }
                            else {
                                px1refl += mass * st[2].real();
                                py1refl += mass * st[3].real();
                            }
                        }
                    }
                    );
            // energy
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&KE, &PE] (const state_t& st) {
                        int s = static_cast<int>(st[6].real());
                        vector<double> eva;
                        matrixop::hdiag(cal_H(vector<double> { st[0].real(), st[1].real() }), eva);
                        KE += 0.5 * mass * pow(st[2].real(), 2) + 0.5 * mass * pow(st[3].real(), 2); 
                        PE += eva[s]; 
                    }
                    );
            // output
            if (istep == 0) {
                // para & header
                ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step, 
                            " mass = ", mass, " A = ", A, " B = ", B, " W = ", W,
                            " init_x = ", init_x, " init_px = ", init_px, 
                            " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                            " init_y = ", init_y, " init_py = ", init_py, 
                            " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                            " init_s = ", init_s, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right
                        );
                ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "py0trans", "px0refl", "py0refl", "px1trans", "py1trans", "px1refl", "py1refl", "Etot");
            }
            if (n0trans > 0.0) { px0trans /= n0trans; py0trans /= n0trans; }
            if (n0refl > 0.0) { px0refl /= n0refl; py0refl /= n0refl; }
            if (n1trans > 0.0) { px1trans /= n1trans; py1trans /= n1trans; }
            if (n1refl > 0.0) { px1refl /= n1refl; py1refl /= n1refl; }
            n0trans /= Ntraj;
            n0refl /= Ntraj;
            n1trans /= Ntraj;
            n1refl /= Ntraj;
            ioer::tabout('#', istep * dt, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl, (KE + PE) / Ntraj);

            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                /*
                double sum0 = 0, sum1= 0;
                for (int itraj = 0; itraj < Ntraj; ++itraj) {
                    complex<double> c0 = state[itraj][4];
                    complex<double> c1 = state[itraj][5];
                    sum0 += abs(c0) * abs(c0);
                    sum1 += abs(c1) * abs(c1);
                }
                sum0 /= Ntraj;
                sum1 /= Ntraj;
                ioer::tabout(n0trans, n1trans, sum0, sum1);
                abort();
                */
                break;
            }
        }

    }
    ioer::tabout(init_px, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl);
    // hop info
    vector<double> hop_count_summary(50, 0.0);
    for_each(hop_count.begin(), hop_count.end(), 
            [&hop_count_summary](double x) { hop_count_summary[static_cast<int>(x)] += 1.0; });
    ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
    ioer::info("# hop count: ", hop_count_summary);
}

void check_surf() {
    for (double x = xwall_left - 1; x < xwall_right + 1; x += 0.01) {
        cal_info(vector<double> {x, 0.0});
        ioer::tabout(x, eva, Fx, Fy, real(dcx), imag(dcy));
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        ioer::info("use --help for detailed info");
        return -1;
    }
    else if (string(argv[1]) == "check") {
        check_surf();
        return 0;
    }
    else {
        if (argparse(argc, argv) == false) {
            return 0;
        }
        randomer::seed(0);
        timer::tic();
        fssh();
        ioer::info("# ", timer::toc());
    }
    return 0;
}
