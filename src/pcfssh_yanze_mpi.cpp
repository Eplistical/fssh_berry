#include <cstdlib>
#include <exception>
#include <iterator>
#include <set>
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
#include "boost/program_options.hpp"
#include <boost/math/quadrature/trapezoidal.hpp>

#include "2d_fssh_rescaling.hpp"
#include "2d_yanze_potential.hpp"

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};


using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::quadrature::trapezoidal;
using state_t = vector< complex<double> >;
const complex<double> zI(0.0, 1.0);

vector<double> potential_params;

double mass = 1000.0;
int Nstep = 1000000;
double dt = 0.01;
double init_s = 0.0;
double init_E = 0.3;

int output_step = 1000;
int Ntraj = 10000;
int seed = 0;
bool enable_hop = true;
bool enable_xp_output = false;
bool enable_phase_corr = false;

string init_pos = "left";
const set<string> init_pos_set { "left", "bottom" };
string rescaling_alg = "m1";
const set<string> rescaling_alg_set { "x", "m1", "m2" };
string pc_rescaling_alg = "none";
const set<string> pc_rescaling_alg_set { "none", "parallel" };
string pc_frustration_alg = "neil";
const set<string> pc_frustration_alg_set { "neil" };
string berry_force_alg = "normal";
const set<string> berry_force_alg_set { "none", "normal", "mean" };

vector< complex<double> > lastevt;
vector<double> eva(2);
vector< complex<double> > Fx(4), Fy(4);
vector< complex<double> > dcx(4), dcy(4);
vector< complex<double> > dx_dcx(4), dy_dcx(4);
vector< complex<double> > dx_dcy(4), dy_dcy(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("mass", po::value<double>(&mass), "mass")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("init_pos", po::value<string>(&init_pos), "init position, left or bottom")
        ("init_E", po::value<double>(&init_E), "init energy")
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("rescaling_alg", po::value<string>(&rescaling_alg), "rescaling algorithm")
        ("berry_force_alg", po::value<string>(&berry_force_alg), "berry force algorithm")
        ("seed", po::value<int>(&seed), "random seed")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_xp_output", po::value<bool>(&enable_xp_output), "enable detailed x, p distribution output")
        ("enable_phase_corr", po::value<bool>(&enable_phase_corr), "enable phase correction")
        ("pc_rescaling_alg", po::value<string>(&pc_rescaling_alg), "phase correction rescaling algorithm")
        ("pc_frustration_alg", po::value<string>(&pc_frustration_alg), "phase correction frustration algorithm")
        ;
    po::variables_map vm; 
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    // check args
    if (not potential_params.empty()) {
        set_potenial_params(potential_params);
    }

    misc::crasher::confirm(init_pos_set.count(init_pos), "Invalid init_pos.");
    misc::crasher::confirm(berry_force_alg_set.count(berry_force_alg), "Invalid berry_force_alg.");
    if (enable_hop) {
        misc::crasher::confirm(rescaling_alg_set.count(rescaling_alg), "Invalid rescaling_alg.");
    }

    if (enable_phase_corr) {
        misc::crasher::confirm(pc_rescaling_alg_set.count(pc_rescaling_alg), "Invalid pc_rescaling_alg.");
        if (pc_rescaling_alg != "none") {
            misc::crasher::confirm(pc_frustration_alg_set.count(pc_frustration_alg), "Invalid pc_frustration_alg.");
        }
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}


double xdist(double x) {
    if (x < param_r - 0.5 * param_R or x > 0.5 * param_R) {
        return 0.0;
    }
    else {
        const double L = param_R - param_r;
        return sqrt(2.0 / L) * sin(M_PI / L * (x + 0.5 * param_R - param_r));
    }
}


template<typename Callable> 
double wtrans(Callable func, double r, double p, double xmin = -10.0, double xmax = 10.0) {
    auto intr = [&func, &r, &p](double x) {
        return func(r-0.5*x) * func(r+0.5*x) * cos(x*p);
    };
    return trapezoidal(intr, xmin, xmax, 1e-6);
}


vector<vector<double>> rp_sample(size_t N, size_t Nstep_eql, size_t Nstep_collect, const vector<double>& rp0, const vector<double>& rpsigma) {
    // sample r&p from Wigner function
    vector<double> rpnow = rp0;
    double wnow = wtrans(xdist, rpnow[0], rpnow[1]);
    vector<double> rpnext(2);
    double wnext;
    // equilibrate
    for (size_t istep(0); istep < Nstep_eql; ++istep) {
        rpnext[0] = rpnow[0] + randomer::normal(0.0, rpsigma[0]);
        rpnext[1] = rpnow[1] + randomer::normal(0.0, rpsigma[1]);
        wnext = wtrans(xdist, rpnext[0], rpnext[1]);
        if (wnext > 0.0 and (wnext > wnow or randomer::rand() < wnext / wnow)) {
            rpnow = rpnext;
            wnow = wnext;
        }
    }
    // sampling
    vector<vector<double>> rst;
    rst.reserve(N);
    size_t Nstep_sample = Nstep_collect * N;
    for (size_t istep(0); istep < Nstep_sample; ++istep) {
        rpnext[0] = rpnow[0] + randomer::normal(0.0, rpsigma[0]);
        rpnext[1] = rpnow[1] + randomer::normal(0.0, rpsigma[1]);
        wnext = wtrans(xdist, rpnext[0], rpnext[1]);
        if (wnext > 0.0 and (wnext > wnow or randomer::rand() < wnext / wnow)) {
            rpnow = rpnext;
            wnow = wnext;
        }
        if (istep % Nstep_collect == 0) {
            rst.push_back(rpnow);
        }
    }
    return rst;
}



void init_state(vector<state_t>& state, int N) {
    state.resize(N);
    const vector<vector<double>> rps = rp_sample(N, 10000, 40, vector<double> { 0.5 * param_r, 0.0 }, vector<double> { (param_R - param_r) * 0.1, 1.0 });
    int idx = 0;
    for (auto& st : state) {
        // st => x, y, vx, vy, c0, c1, s
        st.resize(7, matrixop::ZEROZ);
        if (init_pos == "left") {
            st[0].real(-0.5 * param_R);
            st[1].real(rps[idx][0]);
            st[2].real(sqrt(init_E * 2 / mass));
            st[3].real(rps[idx][1] / mass);
            //st[3].real(0.0);
            idx += 1;
        }
        else if (init_pos == "bottom") {
            st[0].real(rps[idx][0]);
            st[1].real(-0.5 * param_R);
            st[2].real(rps[idx][1] / mass);
            //st[2].real(0.0);
            st[3].real(sqrt(init_E * 2 / mass));
            idx += 1;
        }
        st[4].real(sqrt(1.0 - init_s));
        st[5].real(sqrt(init_s));
        st[6].real((randomer::rand() < init_s) ? 1.0 : 0.0);
    }
}


/*
void init_state(vector<state_t>& state, int N) {
    state.resize(N);
    const vector<double> xs = randomer::MHsample(xdist, 2*N, 10000, 100, 0.5 * param_r, 0.1 * (param_R - param_r));
    int idx = 0;
    for (auto& st : state) {
        // st => x, y, vx, vy, c0, c1, s
        st.resize(7, matrixop::ZEROZ);
        if (init_pos == "left") {
            while (xs[idx] + 0.5*param_R > param_R - 0.0 / param_alpha or xs[idx] + 0.5*param_R < param_r + 0.0 / param_alpha) {
                // xs[idx] out of boundary
                if (idx < xs.size()) {
                    idx += 1;
                } 
                else {
                    throw out_of_range("init_state: insufficient candidates");
                }
            }
            st[0].real(-0.5 * param_R);
            st[1].real(xs[idx]);
            st[2].real(sqrt(init_E * 2 / mass));
            st[3].real(0.0);
            idx += 1;
        }
        else if (init_pos == "bottom") {
            while (xs[idx] + 0.5*param_R > param_R - 2.0 / param_alpha or xs[idx] + 0.5*param_R < param_r + 2.0 / param_alpha) {
                // xs[idx] out of boundary
                if (idx < xs.size()) {
                    idx += 1;
                } 
                else {
                    throw out_of_range("init_state: insufficient candidates");
                }
            }
            st[0].real(xs[idx]);
            st[1].real(-0.5 * param_R);
            st[2].real(0.0);
            st[3].real(sqrt(init_E * 2 / mass));
            idx += 1;
        }
        st[4].real(sqrt(1.0 - init_s));
        st[5].real(sqrt(init_s));
        st[6].real((randomer::rand() < init_s) ? 1.0 : 0.0);
    }
}
*/


void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, v, c0, c1, s
    vector<double> r { state[0].real(), state[1].real() };
    vector<double> v { state[2].real(), state[3].real() };
    vector< complex<double> > c { state[4], state[5] };
    int s = static_cast<int>(state[6].real());
    // extra Berry Force
    double Fx_berry, Fy_berry;
    if (berry_force_alg == "none") {
        Fx_berry = 0.0;
        Fy_berry = 0.0;
    }
    else if (berry_force_alg == "normal") {
        Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
        Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    }
    else if (berry_force_alg == "mean") {
        Fx_berry = 
              pow(abs(c[s]), 2) * 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag()
            + pow(abs(c[1-s]), 2) * 2 * (dcx[1-s+s*2] * (v[0] * dcx[s+(1-s)*2] + v[1] * dcy[s+(1-s)*2])).imag()
            ;
        Fy_berry = 
              pow(abs(c[s]), 2) * 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag()
            + pow(abs(c[1-s]), 2) * 2 * (dcy[1-s+s*2] * (v[0] * dcx[s+(1-s)*2] + v[1] * dcy[s+(1-s)*2])).imag()
            ;
    }
    // state_dot
    state_dot[0] = v[0];
    state_dot[1] = v[1];
    state_dot[2] = (Fx[s+s*2] + Fx_berry) / mass;
    state_dot[3] = (Fy[s+s*2] + Fy_berry) / mass;

    if (pc_rescaling_alg == "none") {
        // normal FSSH
        state_dot[4] = -zI * c[0] * eva[0] - c[1] * (v[0] * dcx[0+1*2] + v[1] * dcy[0+1*2]) - c[0] * (v[0] * dcx[0+0*2] + v[1] * dcy[0+0*2]);
        state_dot[5] = -zI * c[1] * eva[1] - c[0] * (v[0] * dcx[1+0*2] + v[1] * dcy[1+0*2]) - c[1] * (v[0] * dcx[1+1*2] + v[1] * dcy[1+1*2]);
    }
    else {
        // Phase Correction FSSH

        // rescaling direction 
        vector<double> n(2);
        if (pc_rescaling_alg == "parallel") {
            n[0] = v[0];
            n[1] = v[1];
        }
        else if (pc_rescaling_alg == "x") {
            n[0] = 1.0;
            n[1] = 0.0;
        }
        else if (pc_rescaling_alg == "red") {
            n[0] = real(dcx[s+(1-s)*2]);
            n[1] = real(dcy[s+(1-s)*2]);
            double nm = sqrt(n[0]*n[0] + n[1]*n[1]);
            n[0] /= nm;
            n[1] /= nm;
        }
        misc::crasher::confirm(norm(n) > 1e-40, "pc rescaling: norm(n) is 0.");

        if (s == 0) {
            const vector<double> v0 = v;
            // rescale to get v1
            const vector<double> vn = component(v0, n);
            const double vn_norm = norm(vn); 
            const double tmp = vn_norm * vn_norm + 2.0 / mass * (eva[0] - eva[1]);
            if (tmp > 0.0) {
                const double vn_norm_new = sqrt(tmp);
                const vector<double> v1 = v0 - vn + vn_norm_new / vn_norm * vn;
                state_dot[4] = zI * mass * (v0[0] * v0[0] + v0[1] * v0[1]) * c[0] - c[1] * (v0[0] * dcx[0+1*2] + v0[1] * dcy[0+1*2]) - c[0] * (v0[0] * dcx[0+0*2] + v0[1] * dcy[0+0*2]);
                state_dot[5] = zI * mass * (v0[0] * v1[0] + v0[1] * v1[1]) * c[1] - c[0] * (v0[0] * dcx[1+0*2] + v0[1] * dcy[1+0*2]) - c[1] * (v0[0] * dcx[1+1*2] + v0[1] * dcy[1+1*2]);
            }
            else {
                // frustrated
                if (pc_frustration_alg == "neil") {
                    const vector<double> v1 = v0 - vn;
                    state_dot[4] = zI * mass * (v0[0] * v0[0] + v0[1] * v0[1]) * c[0] - c[1] * (v0[0] * dcx[0+1*2] + v0[1] * dcy[0+1*2]) - c[0] * (v0[0] * dcx[0+0*2] + v0[1] * dcy[0+0*2]);
                    state_dot[5] = zI * mass * (v0[0] * v1[0] + v0[1] * v1[1]) * c[1] - c[0] * (v0[0] * dcx[1+0*2] + v0[1] * dcy[1+0*2]) - c[1] * (v0[0] * dcx[1+1*2] + v0[1] * dcy[1+1*2]);
                }
                else if (pc_frustration_alg == "joe") {
                    state_dot[4] = -zI * c[0] * 2.0 * eva[0] - c[1] * (v[0] * dcx[0+1*2] + v[1] * dcy[0+1*2]) - c[0] * (v[0] * dcx[0+0*2] + v[1] * dcy[0+0*2]);
                    state_dot[5] = -zI * c[1] * 2.0 * eva[1] - c[0] * (v[0] * dcx[1+0*2] + v[1] * dcy[1+0*2]) - c[1] * (v[0] * dcx[1+1*2] + v[1] * dcy[1+1*2]);
                }
            }
        }
        else if (s == 1) {
            const vector<double> v1 = v;
            // rescale to get v0
            const vector<double> vn = component(v1, n);
            const double vn_norm = norm(vn); 
            const double tmp = vn_norm * vn_norm + 2.0 / mass * (eva[1] - eva[0]);
            const double vn_norm_new = sqrt(tmp);
            const vector<double> v0 = v1 - vn + vn_norm_new / vn_norm * vn;
            state_dot[4] = zI * mass * (v1[0] * v0[0] + v1[1] * v0[1]) * c[0] - c[1] * (v1[0] * dcx[0+1*2] + v1[1] * dcy[0+1*2]) - c[0] * (v1[0] * dcx[0+0*2] + v1[1] * dcy[0+0*2]);
            state_dot[5] = zI * mass * (v1[0] * v1[0] + v1[1] * v1[1]) * c[1] - c[0] * (v1[0] * dcx[1+0*2] + v1[1] * dcy[1+0*2]) - c[1] * (v1[0] * dcx[1+1*2] + v1[1] * dcy[1+1*2]);
        }

    }
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
        // generate rescaling direction
        const int from = s;
        const int to = 1 - s;
        vector<double> n;

        n = get_rescaling_direction(rescaling_alg, r, v, c, from, to,
                dcx, dcy, dx_dcx, dy_dcx, dx_dcy, dy_dcy, Fx, Fy, eva);

        // rescale momentum
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
    const vector<double> r { state[0].real(), state[1].real() };
    const vector<double> v { state[2].real(), state[3].real() };
    const double m_theta = cal_m_theta(r);
    return (m_theta < 0.0 and v[1] < 0.0) or (m_theta > M_PI/2 and v[0] < 0.0);
}

void fssh() {
    // assign job
    vector<int> my_jobs = MPIer::assign_job(Ntraj);
    int my_Ntraj = my_jobs.size();
    // propagation variables
    runge_kutta4<state_t> rk4;
    // initialize
    vector<state_t> state;
    init_state(state, my_Ntraj);
    // statistics
    double n0left = 0.0, n0bot = 0.0, n1left = 0.0, n1bot = 0.0;
    double px0left = 0.0, px0bot = 0.0, px1left = 0.0, px1bot = 0.0;
    double py0left = 0.0, py0bot = 0.0, py1left = 0.0, py1bot = 0.0;
    double KE = 0.0, PE = 0.0;
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    // recorders
    int Nrec = Nstep / output_step;
    vector<double> n0left_arr(Nrec), n0bot_arr(Nrec), n1left_arr(Nrec), n1bot_arr(Nrec);
    vector<double> px0left_arr(Nrec), px0bot_arr(Nrec), px1left_arr(Nrec), px1bot_arr(Nrec);
    vector<double> py0left_arr(Nrec), py0bot_arr(Nrec), py1left_arr(Nrec), py1bot_arr(Nrec);
    vector<double> KE_arr(Nrec), PE_arr(Nrec);
    vector<double> hop_count_summary(50, 0.0);
    // for x,p distribution output
    vector<vector<double>> xp_arr(Nrec);
    // main loop
    vector< vector< complex<double> > > lastevt_save(my_Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < my_Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info_nume(
                        vector<double> { state[itraj][0].real(), state[itraj][1].real() },
                        Fx, Fy, dcx, dcy, eva, lastevt
                        );
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
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
            int irec = istep / output_step;
            // population
            n0left = n0bot = n1left = n1bot = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0left, &n0bot, &n1left, &n1bot](const state_t& st) { 
                        int s = static_cast<int>(st[6].real());
                        double x = st[0].real();
                        double y = st[1].real();
                        if (s == 0) {
                            (x < y) ? n0left += 1.0 : n0bot += 1.0;
                        }
                        else {
                            (x < y) ? n1left += 1.0 : n1bot += 1.0;
                        }
                    });
            n0left_arr[irec] = n0left;
            n0bot_arr[irec] = n0bot;
            n1left_arr[irec] = n1left;
            n1bot_arr[irec] = n1bot;
            // momentum
            px0left = px0bot = px1left = px1bot = 0.0;
            py0left = py0bot = py1left = py1bot = 0.0;
            for_each(state.begin(), state.end(), 
                    [&px0left, &px0bot, &px1left, &px1bot, &py0left, &py0bot, &py1left, &py1bot] (const state_t& st) {
                        int s = static_cast<int>(st[6].real());
                        double x = st[0].real();
                        double y = st[1].real();
                        if (s == 0) {
                            if (x < y) {
                                px0left += mass * st[2].real();
                                py0left += mass * st[3].real();
                            }
                            else {
                                px0bot += mass * st[2].real();
                                py0bot += mass * st[3].real();
                            }
                        }
                        else {
                            if (x < y) {
                                px1left += mass * st[2].real();
                                py1left += mass * st[3].real();
                            }
                            else {
                                px1bot += mass * st[2].real();
                                py1bot += mass * st[3].real();
                            }
                        }
                    }
                    );
            px0left_arr[irec] = px0left;
            px0bot_arr[irec] = px0bot;
            px1left_arr[irec] = px1left;
            px1bot_arr[irec] = px1bot;
            py0left_arr[irec] = py0left;
            py0bot_arr[irec] = py0bot;
            py1left_arr[irec] = py1left;
            py1bot_arr[irec] = py1bot;
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
            // x,p dist
            if (enable_xp_output) {
                for_each(state.begin(), state.end(), 
                        [irec, &xp_arr] (const state_t& st) {
                            xp_arr[irec].push_back(st[0].real()); // x
                            xp_arr[irec].push_back(st[1].real()); // y
                            xp_arr[irec].push_back(st[2].real()); // vx
                            xp_arr[irec].push_back(st[3].real()); // vy
                            xp_arr[irec].push_back(st[6].real()); // s
                        });
            }
            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                // fill the rest
                fill(n0left_arr.begin() + irec + 1, n0left_arr.end(), n0left);
                fill(n0bot_arr.begin() + irec + 1, n0bot_arr.end(), n0bot);
                fill(n1left_arr.begin() + irec + 1, n1left_arr.end(), n1left);
                fill(n1bot_arr.begin() + irec + 1, n1bot_arr.end(), n1bot);

                fill(px0left_arr.begin() + irec + 1, px0left_arr.end(), px0left);
                fill(px0bot_arr.begin() + irec + 1, px0bot_arr.end(), px0bot);
                fill(px1left_arr.begin() + irec + 1, px1left_arr.end(), px1left);
                fill(px1bot_arr.begin() + irec + 1, px1bot_arr.end(), px1bot);

                fill(py0left_arr.begin() + irec + 1, py0left_arr.end(), py0left);
                fill(py0bot_arr.begin() + irec + 1, py0bot_arr.end(), py0bot);
                fill(py1left_arr.begin() + irec + 1, py1left_arr.end(), py1left);
                fill(py1bot_arr.begin() + irec + 1, py1bot_arr.end(), py1bot);

                fill(KE_arr.begin() + irec + 1, KE_arr.end(), KE);
                fill(PE_arr.begin() + irec + 1, PE_arr.end(), PE);

                // xp dist
                if (enable_xp_output) {
                    for (int ii(irec+1); ii < Nrec; ++ii) {
                        xp_arr[ii] = xp_arr[irec];
                    }
                }
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
                    n0left_arr, n0bot_arr, n1left_arr, n1bot_arr, 
                    px0left_arr, px0bot_arr, px1left_arr, px1bot_arr, 
                    py0left_arr, py0bot_arr, py1left_arr, py1bot_arr, 
                    KE_arr, PE_arr,
                    hopup, hopdn, hopfr, hoprj, hop_count_summary
                    );

            // xp dist
            if (enable_xp_output) {
                for (int irec(0); irec < Nrec; ++irec) {
                    MPIer::send(0, xp_arr[irec]);
                }
            }
        }
        else if (MPIer::master) {
            vector<double> buf;

            MPIer::recv(r, buf); n0left_arr += buf;
            MPIer::recv(r, buf); n0bot_arr += buf;
            MPIer::recv(r, buf); n1left_arr += buf;
            MPIer::recv(r, buf); n1bot_arr += buf;

            MPIer::recv(r, buf); px0left_arr += buf;
            MPIer::recv(r, buf); px0bot_arr += buf;
            MPIer::recv(r, buf); px1left_arr += buf;
            MPIer::recv(r, buf); px1bot_arr += buf;

            MPIer::recv(r, buf); py0left_arr += buf;
            MPIer::recv(r, buf); py0bot_arr += buf;
            MPIer::recv(r, buf); py1left_arr += buf;
            MPIer::recv(r, buf); py1bot_arr += buf;

            MPIer::recv(r, buf); KE_arr += buf;
            MPIer::recv(r, buf); PE_arr += buf;

            double dbuf;
            MPIer::recv(r, dbuf); hopup += dbuf;
            MPIer::recv(r, dbuf); hopdn += dbuf;
            MPIer::recv(r, dbuf); hopfr += dbuf;
            MPIer::recv(r, dbuf); hoprj += dbuf;

            MPIer::recv(r, buf); hop_count_summary += buf;

            // xp dist
            if (enable_xp_output) {
                for (int irec(0); irec < Nrec; ++irec) {
                    MPIer::recv(r, buf); 
                    xp_arr[irec].insert(xp_arr[irec].end(), make_move_iterator(buf.begin()), make_move_iterator(buf.end()));
                }
            }
        }
        MPIer::barrier();
    }
    // output
    if (MPIer::master) {
        // para & header
        output_potential_param();
        ioer::info("# fsshx para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step,
                " mass = ", mass, 
                " init_E = ", init_E, 
                " init_pos = ", init_pos,
                " init_s = ", init_s,
                " berry_force_alg = ", berry_force_alg,
                " enable_hop = ", enable_hop, 
                " enable_xp_output = ", enable_xp_output, 
                " rescaling_alg = ", rescaling_alg,
                " enable_phase_corr = ", enable_phase_corr, 
                " pc_rescaling_alg = ", pc_rescaling_alg,
                " pc_frustration_alg = ", pc_frustration_alg,
                ""
                );
        ioer::tabout('#', "t", "n0left", "n0bot", "n1left", "n1bot", "px0left", "py0left", "px0bot", "py0bot", "px1left", "py1left", "px1bot", "py1bot", "etot");
        for (int irec = 0; irec < Nrec; ++irec) {
            n0left = n0left_arr[irec];
            n0bot = n0bot_arr[irec];
            n1left = n1left_arr[irec];
            n1bot = n1bot_arr[irec];

            px0left = px0left_arr[irec];
            px0bot = px0bot_arr[irec];
            px1left = px1left_arr[irec];
            px1bot = px1bot_arr[irec];

            py0left = py0left_arr[irec];
            py0bot = py0bot_arr[irec];
            py1left = py1left_arr[irec];
            py1bot = py1bot_arr[irec];

            KE = KE_arr[irec];
            PE = PE_arr[irec];

            if (n0left > 0.0) { px0left /= n0left; py0left /= n0left; }
            if (n0bot > 0.0) { px0bot /= n0bot; py0bot /= n0bot; }
            if (n1left > 0.0) { px1left /= n1left; py1left /= n1left; }
            if (n1bot > 0.0) { px1bot /= n1bot; py1bot /= n1bot; }

            n0left /= Ntraj;
            n0bot /= Ntraj;
            n1left /= Ntraj;
            n1bot /= Ntraj;

            ioer::tabout('#', irec * output_step * dt, n0left, n0bot, n1left, n1bot, px0left, py0left, px0bot, py0bot, px1left, py1left, px1bot, py1bot, (KE + PE) / Ntraj);
        }
        // final results
        ioer::tabout(init_E, n0left, n0bot, n1left, n1bot, px0left, py0left, px0bot, py0bot, px1left, py1left, px1bot, py1bot);
        // hop info
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
        // xp dist
        if (enable_xp_output) {
            ioer::info("## detailed distribution info for each trajectory:");
            for (int irec = 0; irec < Nrec; ++irec) {
                ioer::info("## ", xp_arr[irec]);
            }
        }
    }
    MPIer::barrier();
}

void check_surf() {
    for (double x = -5.0; x < 5.0; x += 0.05) {
        for (double y = -5.0; y < 5.0; y += 0.05) {
            cal_info_nume(vector<double> {x, y}, Fx, Fy, dcx, dcy, eva, lastevt);
            ioer::tabout(x, y, eva, 
                    dcx[0+1*2].real(), dcx[0+1*2].imag(), abs(dcx[0+1*2]), 
                    dcy[0+1*2].real(), dcy[0+1*2].imag(), abs(dcy[0+1*2])
                    );
        }
    }
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
