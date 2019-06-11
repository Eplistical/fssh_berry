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
#include "misc/MPIer.hpp"
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
bool enable_hop = true;
bool enable_nodjj = false;
string output_mod = "init_px";

vector< complex<double> > lastevt;
vector<double> eva(2);
vector<double> Fx(2), Fy(2);
vector< complex<double> > dcx(4), dcy(4);
//vector< complex<double> > Tmat(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("A", po::value<double>(&A), "potential para A")
        ("B", po::value<double>(&B), "potential para B")
        ("C", po::value<double>(&C), "potential para C")
        ("D", po::value<double>(&D), "potential para D")
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
        ("seed", po::value<int>(&seed), "random seed")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_nodjj", po::value<bool>(&enable_nodjj), "let djj = 0")
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

double cal_phi(const double y) {
    return W * y;
}

double cal_der_phi(const double y) {
    return W;
}

vector< complex<double> > cal_H(const vector<double>& r) {
    const double x = r[0];
    const double y = r[1];
    const complex<double> eip = exp(matrixop::IMAGIZ * cal_phi(y));

    vector< complex<double> > H(4);
    if (x >= 0.0) {
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

void init_state(state_t& state) {
    state.resize(7, matrixop::ZEROZ);
    // state = x, y, vx, vy, c0, c1, s
    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_y, sigma_y)); 

    state[2].real(randomer::normal(init_px, sigma_px) / mass); 
    while (state[2].real() <= 0.0) {
        state[2].real(randomer::normal(init_px, sigma_px) / mass); 
    }

    state[3].real(randomer::normal(init_py, sigma_py) / mass); 
    state[4].real(sqrt(1.0 - init_s));
    state[5].real(sqrt(init_s));
    state[6].real((randomer::rand() < init_s) ? 1.0 : 0.0);
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
    Fx.assign(2, 0.0);
    Fy.assign(2, 0.0);
    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            if (j == k) {
                Fx[j] = -dcx[j+k*2].real();
                Fy[j] = -dcy[j+k*2].real();
                dcx[j+k*2] = 0.0;
                dcy[j+k*2] = 0.0;
            }
            else {
                dcx[j+k*2] /= (eva[k] - eva[j]);
                dcy[j+k*2] /= (eva[k] - eva[j]);
            }
        }
    }
    /*
    // Tmat
    if (not lastevt.empty()) {
        // U
        matrixop::hdiag(cal_H(r), eva, evt);
        vector< complex<double> > U;
        U = matrixop::matCmat(lastevt, evt, 2);
        for (int i(0); i < 2; ++i) {
            complex<double> phase = abs(U[i+i*2]) / U[i+i*2];
            for (int j(0); j < 2; ++j) {
                evt[j+i*2] *= phase;
            }
        }
        U = matrixop::matCmat(lastevt, evt, 2);
        // Tmat
        Tmat = matrixop::logmh(U) / dt;
    }
    */
    // save evt
    lastevt = move(evt);
}

void cal_info(const vector<double>& r)
{
    cal_info_nume(r);
    if (enable_nodjj) {
        dcx[0+0*2] = 0.0;
        dcx[1+1*2] = 0.0;
        dcy[0+0*2] = 0.0;
        dcy[1+1*2] = 0.0;
    }
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


void integtrate(state_t& state, const double dt) {
    // integrate state dt forward, assuming all info needed has been evaluated
    // state = x, v, c0, c1, s
    vector<double> r { state[0].real(), state[1].real() };
    vector<double> v { state[2].real(), state[3].real() };
    vector< complex<double> > c { state[4], state[5] };
    int s = static_cast<int>(state[6].real());
    // extra Berry Force
    double Fx_berry, Fy_berry;

    // nuclear part, VV
    cal_info(r);
    Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    v[0] += 0.5 * (Fx[s] + Fx_berry) / mass * dt;
    v[1] += 0.5 * (Fy[s] + Fy_berry) / mass * dt;

    r += v * dt;

    cal_info(r);
    Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    v[0] += 0.5 * (Fx[s] + Fx_berry) / mass * dt;
    v[1] += 0.5 * (Fy[s] + Fy_berry) / mass * dt;

    state[0].real(r[0]);
    state[1].real(r[1]);
    state[2].real(v[0]);
    state[3].real(v[1]);

    // electron part, RK4
    auto rk4_func = [&v](const vector< complex<double> >& c) {
        vector< complex<double> > cdot(2);
        cdot[0] = -zI * c[0] * eva[0] - c[1] * (v[0] * dcx[0+1*2] + v[1] * dcy[0+1*2]) - c[0] * (v[0] * dcx[0+0*2] + v[1] * dcy[0+0*2]);
        cdot[1] = -zI * c[1] * eva[1] - c[0] * (v[0] * dcx[1+0*2] + v[1] * dcy[1+0*2]) - c[1] * (v[0] * dcx[1+1*2] + v[1] * dcy[1+1*2]);
        return cdot;
    };

    vector< complex<double> > k1, k2, k3, k4;
    k1 = dt * rk4_func(c);
    k2 = dt * rk4_func(c + 0.5 * k1);
    k3 = dt * rk4_func(c + 0.5 * k2);
    k4 = dt * rk4_func(c + k3);
    c += 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

    state[4] = c[0];
    state[5] = c[1];
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

        // Method #2
        const int from = s;
        const int to = 1 - s;
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

        // debug
        /*
        n[0] = 1.0;
        n[1] = 0.0;
        */

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
    // assign job
    vector<int> my_jobs = MPIer::assign_job(Ntraj);
    int my_Ntraj = my_jobs.size();
    vector<state_t> state(my_Ntraj);
    // propagation variables
    runge_kutta4<state_t> rk4;
    // initialize
    for (int itraj(0); itraj < my_Ntraj; ++itraj) {
        init_state(state[itraj]);
    }
    // statistics
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    double py0trans = 0.0, py0refl = 0.0, py1trans = 0.0, py1refl = 0.0;
    double KE = 0.0, PE = 0.0;
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    // recorders
    int Nrec = Nstep / output_step;
    vector<double> n0trans_arr(Nrec), n0refl_arr(Nrec), n1trans_arr(Nrec), n1refl_arr(Nrec);
    vector<double> px0trans_arr(Nrec), px0refl_arr(Nrec), px1trans_arr(Nrec), px1refl_arr(Nrec);
    vector<double> py0trans_arr(Nrec), py0refl_arr(Nrec), py1trans_arr(Nrec), py1refl_arr(Nrec);
    vector<double> KE_arr(Nrec), PE_arr(Nrec);
    vector<double> hop_count_summary(50, 0.0);
    // main loop
    vector< vector< complex<double> > > lastevt_save(my_Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < my_Ntraj; ++itraj) {
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
                        case HOP_FR : { hopfr += 1.0; break; }
                        case HOP_RJ : { hoprj += 1.0; break; }
                        default : break;
                    }
                }
                // propagate
                rk4.do_step(sys, state[itraj], istep * dt, dt);
                //integtrate(state[itraj], dt);
                
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
            int irec = istep / output_step;
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
            n0trans_arr[irec] = n0trans;
            n0refl_arr[irec] = n0refl;
            n1trans_arr[irec] = n1trans;
            n1refl_arr[irec] = n1refl;
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
            px0trans_arr[irec] = px0trans;
            px0refl_arr[irec] = px0refl;
            px1trans_arr[irec] = px1trans;
            px1refl_arr[irec] = px1refl;
            py0trans_arr[irec] = py0trans;
            py0refl_arr[irec] = py0refl;
            py1trans_arr[irec] = py1trans;
            py1refl_arr[irec] = py1refl;
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
            KE_arr[irec] = KE;
            PE_arr[irec] = PE;
            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                // fill the rest
                fill(n0trans_arr.begin() + irec + 1, n0trans_arr.end(), n0trans);
                fill(n0refl_arr.begin() + irec + 1, n0refl_arr.end(), n0refl);
                fill(n1trans_arr.begin() + irec + 1, n1trans_arr.end(), n1trans);
                fill(n1refl_arr.begin() + irec + 1, n1refl_arr.end(), n1refl);

                fill(px0trans_arr.begin() + irec + 1, px0trans_arr.end(), px0trans);
                fill(px0refl_arr.begin() + irec + 1, px0refl_arr.end(), px0refl);
                fill(px1trans_arr.begin() + irec + 1, px1trans_arr.end(), px1trans);
                fill(px1refl_arr.begin() + irec + 1, px1refl_arr.end(), px1refl);

                fill(py0trans_arr.begin() + irec + 1, py0trans_arr.end(), py0trans);
                fill(py0refl_arr.begin() + irec + 1, py0refl_arr.end(), py0refl);
                fill(py1trans_arr.begin() + irec + 1, py1trans_arr.end(), py1trans);
                fill(py1refl_arr.begin() + irec + 1, py1refl_arr.end(), py1refl);

                fill(KE_arr.begin() + irec + 1, KE_arr.end(), KE);
                fill(PE_arr.begin() + irec + 1, PE_arr.end(), PE);
                break;
            }
        }
    }
    // process data
    MPIer::barrier();
    for_each(hop_count.begin(), hop_count.end(),
            [&hop_count_summary](double x) { hop_count_summary[static_cast<int>(x)] += 1.0; });
    // collect data
    MPIer::barrier();
    for (int r = 1; r < MPIer::size; ++r) {
        if (MPIer::rank == r) {
            MPIer::send(0, 
                    n0trans_arr, n0refl_arr, n1trans_arr, n1refl_arr, 
                    px0trans_arr, px0refl_arr, px1trans_arr, px1refl_arr, 
                    py0trans_arr, py0refl_arr, py1trans_arr, py1refl_arr, 
                    KE_arr, PE_arr,
                    hopup, hopdn, hopfr, hoprj, hop_count_summary
                    );
        }
        else if (MPIer::master) {
            vector<double> buf;

            MPIer::recv(r, buf); n0trans_arr += buf;
            MPIer::recv(r, buf); n0refl_arr += buf;
            MPIer::recv(r, buf); n1trans_arr += buf;
            MPIer::recv(r, buf); n1refl_arr += buf;

            MPIer::recv(r, buf); px0trans_arr += buf;
            MPIer::recv(r, buf); px0refl_arr += buf;
            MPIer::recv(r, buf); px1trans_arr += buf;
            MPIer::recv(r, buf); px1refl_arr += buf;

            MPIer::recv(r, buf); py0trans_arr += buf;
            MPIer::recv(r, buf); py0refl_arr += buf;
            MPIer::recv(r, buf); py1trans_arr += buf;
            MPIer::recv(r, buf); py1refl_arr += buf;

            MPIer::recv(r, buf); KE_arr += buf;
            MPIer::recv(r, buf); PE_arr += buf;

            double dbuf;
            MPIer::recv(r, dbuf); hopup += dbuf;
            MPIer::recv(r, dbuf); hopdn += dbuf;
            MPIer::recv(r, dbuf); hopfr += dbuf;
            MPIer::recv(r, dbuf); hoprj += dbuf;

            MPIer::recv(r, buf); hop_count_summary += buf;
        }
        MPIer::barrier();
    }
    // output
    if (MPIer::master) {
        // para & header
        ioer::info("# fssh para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step, " output_mod = ", output_mod,
                " mass = ", mass, " A = ", A, " B = ", B, " W = ", W,
                " init_x = ", init_x, " init_px = ", init_px, 
                " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                " init_y = ", init_y, " init_py = ", init_py, 
                " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                " init_s = ", init_s, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right
                );
        ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "py0trans", "px0refl", "py0refl", "px1trans", "py1trans", "px1refl", "py1refl", "etot");
        for (int irec = 0; irec < Nrec; ++irec) {
            n0trans = n0trans_arr[irec];
            n0refl = n0refl_arr[irec];
            n1trans = n1trans_arr[irec];
            n1refl = n1refl_arr[irec];

            px0trans = px0trans_arr[irec];
            px0refl = px0refl_arr[irec];
            px1trans = px1trans_arr[irec];
            px1refl = px1refl_arr[irec];

            py0trans = py0trans_arr[irec];
            py0refl = py0refl_arr[irec];
            py1trans = py1trans_arr[irec];
            py1refl = py1refl_arr[irec];

            KE = KE_arr[irec];
            PE = PE_arr[irec];

            if (n0trans > 0.0) { px0trans /= n0trans; py0trans /= n0trans; }
            if (n0refl > 0.0) { px0refl /= n0refl; py0refl /= n0refl; }
            if (n1trans > 0.0) { px1trans /= n1trans; py1trans /= n1trans; }
            if (n1refl > 0.0) { px1refl /= n1refl; py1refl /= n1refl; }

            n0trans /= Ntraj;
            n0refl /= Ntraj;
            n1trans /= Ntraj;
            n1refl /= Ntraj;

            ioer::tabout('#', irec * output_step * dt, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl, (KE + PE) / Ntraj);
        }
        // final results
        ioer::tabout(init_px, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl);
        // hop info
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
    }
    MPIer::barrier();
}

void check_surf() {
    /*
    W = 5.0;
    for (double x = xwall_left - 1; x < xwall_right + 1; x += 0.01) {
        cal_info(vector<double> {x, 0.0});
        ioer::tabout(x, eva, Fx, Fy, 
                dcx[0+1*2].real(), dcx[0+1*2].imag(), abs(dcx[0+1*2]), 
                dcy[0+1*2].real(), dcy[0+1*2].imag(), abs(dcy[0+1*2])
                );
    }
    dt = 0.1;
    vector<double> r {0, 0};
    vector<double> v {0.01, 0.001};

    cal_info(r + 0.5 * v * dt);
    complex<double> vd01 = v[0] * dcx[0+1*2] + v[1] * dcy[0+1*2];
    ioer::tabout(v, dcx[0+1*2], dcy[0+1*2]);
    ioer::info("vd01 = ", vd01, " => ", abs(vd01));

    // U
    lastevt.clear();
    cal_info(r);
    vector< complex<double> > evt;

    matrixop::hdiag(cal_H(r), eva, lastevt);
    ioer::info(" ** r = ", r);
    ioer::info(" ** H = ", cal_H(r));
    ioer::info(" ** eva = ", eva);
    ioer::info(" ** evt = ", lastevt);
    matrixop::hdiag(cal_H(r + v * dt), eva, evt);
    ioer::info(" **** r + v * dt = ", r + v * dt);
    ioer::info(" **** H = ", cal_H(r + v * dt));
    ioer::info(" **** eva = ", eva);
    ioer::info(" **** evt = ", evt);

    vector< complex<double> > U;
    U = matrixop::matCmat(lastevt, evt, 2);
    for (int i(0); i < 2; ++i) {
        complex<double> phase = abs(U[i+i*2]) / U[i+i*2];
        for (int j(0); j < 2; ++j) {
            U[j+i*2] *= phase;
        }
    } 
    ioer::info("U = ", U);
    // Tmat
    //Tmat[0+1*2] = 0.5 / dt * (U[0+1*2] - U[1+0*2]);

    Tmat = matrixop::logmh(U) / dt;
    ioer::newline();
    ioer::info("logm(U)");
    ioer::info("T = ", Tmat);
    ioer::info("Tm01 = ", Tmat[0+1*2], " => ", abs(Tmat[0+1*2]));


    Tmat = 0.5 / dt * (U - matrixop::transpose(U, 2));
    ioer::newline();
    ioer::info("U-trans(U)");
    ioer::info("T = ", Tmat);
    ioer::info("Tm01 = ", Tmat[0+1*2], " => ", abs(Tmat[0+1*2]));

    Tmat = 0.5 / dt * (U - matrixop::adjoint(U, 2));
    ioer::newline();
    ioer::info("U-adj(U)");
    ioer::info("T = ", Tmat);
    ioer::info("Tm01 = ", Tmat[0+1*2], " => ", abs(Tmat[0+1*2]));
    */
}

int main(int argc, char** argv) {
    MPIer::setup();
    if (argc < 2) {
        if (MPIer::master) ioer::info("use --help for detailed info");
    }
    else if (string(argv[1]) == "check") {
        if (MPIer::master) check_surf();
    }
    else {
        if (argparse(argc, argv) == false) {
            return 0;
        }
        randomer::seed(MPIer::assign_random_seed(seed));
        if (MPIer::master) timer::tic();
        fssh();
        if (MPIer::master) ioer::info("# ", timer::toc());
    }
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
