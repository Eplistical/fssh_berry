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
#include "boost/program_options.hpp"

#include "2d_fssh_rescaling.hpp"
#include "2d_helix_potential.hpp"

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};


using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using state_t = vector< complex<double> >;
const complex<double> zI(0.0, 1.0);

vector<double> potential_params;

const double kT = 9.5e-4;
const double mass = 2000.0;

double init_x = 0.0;
double sigma_x = 0.0; 
double init_px = 0.0;
double sigma_px = sqrt(mass * kT); 
double fric_gamma_x = 0.0;
double init_y = 5.0;
double sigma_y = 0.0; 
double init_py = 0.0;
double sigma_py = sqrt(mass * kT); 
double fric_gamma_y = 0.0;

double init_s = 0.0;
int Nstep = 200000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 10000;
int seed = 0;
bool enable_hop = true;
bool enable_bath = true;
bool enable_berry_force = true;
string rescaling_alg = "x";

vector< complex<double> > lastevt;
vector<double> eva(2);
vector< complex<double> > Fx(4), Fy(4);
vector< complex<double> > dcx(4), dcy(4);
vector< complex<double> > dx_dcx(4), dy_dcy(4);

// bath
double Fx_random = 0.0, Fy_random = 0.0;

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("init_x", po::value<double>(&init_x), "init x")
        ("sigma_x", po::value<double>(&sigma_x), "init sigma x")
        ("init_px", po::value<double>(&init_px), "potential para init_px")
        ("sigma_px", po::value<double>(&sigma_px), "init sigma px")
        ("fric_gamma_x", po::value<double>(&fric_gamma_x), "friction gamma on x direction")
        ("init_y", po::value<double>(&init_y), "init y")
        ("sigma_y", po::value<double>(&sigma_y), "init sigma y")
        ("init_py", po::value<double>(&init_py), "potential para init_py")
        ("sigma_py", po::value<double>(&sigma_py), "init sigma py")
        ("fric_gamma_y", po::value<double>(&fric_gamma_y), "friction gamma on y direction")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("rescaling_alg", po::value<string>(&rescaling_alg), "rescaling algorithm")
        ("seed", po::value<int>(&seed), "random seed")
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_bath", po::value<bool>(&enable_bath), "enable bath")
        ("enable_berry_force", po::value<bool>(&enable_berry_force), "enable Berry force")
        ;
    po::variables_map vm; 
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);    

    if (not potential_params.empty()) {
        set_potenial_params(potential_params);
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
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
}

void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, v, c0, c1, s
    vector<double> r { state[0].real(), state[1].real() };
    vector<double> v { state[2].real(), state[3].real() };
    vector< complex<double> > c { state[4], state[5] };
    int s = static_cast<int>(state[6].real());
    // nuclear friction from bath, F_random should already be initialized
    double Fx_nufric = 0.0, Fy_nufric = 0.0;
    if (enable_bath) {
        Fx_nufric = -fric_gamma_x * v[0] * mass + Fx_random;
        Fy_nufric = -fric_gamma_y * v[0] * mass + Fy_random;
    }
    // extra Berry Force
    double Fx_berry = 0.0, Fy_berry = 0.0;
    if (enable_berry_force) {
        Fx_berry = 2 * (dcx[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
        Fy_berry = 2 * (dcy[s+(1-s)*2] * (v[0] * dcx[1-s+s*2] + v[1] * dcy[1-s+s*2])).imag();
    }
    // state_dot
    state_dot[0] = v[0];
    state_dot[1] = v[1];
    state_dot[2] = (Fx[s+s*2] + Fx_nufric + Fx_berry) / mass;
    state_dot[3] = (Fy[s+s*2] + Fy_nufric + Fy_berry) / mass;
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
        // generate rescaling direction
        const int from = s;
        const int to = 1 - s;
        vector<double> n = get_rescaling_direction(rescaling_alg, r, v, c, from, to,
                                                        dcx, dcy, dx_dcx, dy_dcy, Fx, Fy, eva);
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
    return false;
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
    double n0d = 0.0, n1d = 0.0;
    double KE = 0.0, PE = 0.0;
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    vector<double> hop_count(my_Ntraj, 0.0);
    // recorders
    int Nrec = Nstep / output_step;
    vector<double> n0d_arr(Nrec), n1d_arr(Nrec);
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
                // initialize random force
                if (enable_bath) {
                    Fx_random = randomer::normal(0.0, sqrt(2.0 * fric_gamma_x * mass * kT / dt));
                    Fy_random = randomer::normal(0.0, sqrt(2.0 * fric_gamma_y * mass * kT / dt));
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
            n0d = n1d = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0d, &n1d] (const state_t& st) { 
                        vector<double> r { st[0].real(), st[1].real() };
                        vector<double> v { st[2].real(), st[3].real() };
                        vector< complex<double> > c { st[4], st[5] };
                        int s = static_cast<int>(st[6].real());

                        // calc eva & evt
                        vector<double> eva;
                        vector< complex<double> > U;
                        matrixop::hdiag(cal_H(r), eva, U);

                        // population
                        n0d += pow(abs(U[0+s*2]), 2) + 2 * (U[0+0*2] * c[0] * conj(c[1]) * conj(U[0+1*2])).real();
                        n1d += pow(abs(U[1+s*2]), 2) + 2 * (U[1+0*2] * c[0] * conj(c[1]) * conj(U[1+1*2])).real();
                    });
            n0d_arr[irec] = n0d;
            n1d_arr[irec] = n1d;
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
                fill(n0d_arr.begin() + irec + 1, n0d_arr.end(), n0d);
                fill(n1d_arr.begin() + irec + 1, n1d_arr.end(), n1d);

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
                    n0d_arr, n1d_arr,
                    KE_arr, PE_arr,
                    hopup, hopdn, hopfr, hoprj, hop_count_summary
                    );
        }
        else if (MPIer::master) {
            vector<double> buf;

            MPIer::recv(r, buf); n0d_arr += buf;
            MPIer::recv(r, buf); n1d_arr += buf;

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
        output_potential_param();
        ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step,
                " mass = ", mass, 
                " init_x = ", init_x, " init_px = ", init_px, 
                " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                " fric_gamma_x = ", fric_gamma_x,
                " init_y = ", init_y, " init_py = ", init_py, 
                " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                " fric_gamma_y = ", fric_gamma_y,
                " init_s = ", init_s, 
                " enable_hop = ", enable_hop,
                " enable_bath = ", enable_bath,
                " enable_berry_force = ", enable_berry_force, 
                " rescaling_alg = ", rescaling_alg,
                ""
                );
        ioer::tabout('#', "t", "n0d", "n1d", "KE", "PE", "Etot");
        for (int irec = 0; irec < Nrec; ++irec) {
            n0d = n0d_arr[irec] / Ntraj;
            n1d = n1d_arr[irec] / Ntraj;

            KE = KE_arr[irec] / Ntraj;
            PE = PE_arr[irec] / Ntraj;

            ioer::tabout(' ', irec * output_step * dt, n0d, n1d, KE, PE, (PE + KE));
        }
        // hop info
        ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
        ioer::info("# hop count: ", hop_count_summary);
    }
    MPIer::barrier();
}

void check_surf() {
    for (double x = -10; x < 10; x += 0.01) {
        cal_info_nume(vector<double> {x, 0.0}, Fx, Fy, dcx, dcy, eva, lastevt);
        ioer::tabout(x, eva, 
                dcx[0+1*2].real(), dcx[0+1*2].imag(), abs(dcx[0+1*2]), 
                dcy[0+1*2].real(), dcy[0+1*2].imag(), abs(dcy[0+1*2])
                );
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
