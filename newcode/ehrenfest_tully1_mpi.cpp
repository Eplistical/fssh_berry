#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/matrixop.hpp"
#include "misc/MPIer.hpp"
#include "misc/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/program_options.hpp"

#include "2d_tully1_potential.hpp"

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using state_t = vector< complex<double> >;

// global const
const complex<double> zI(0.0, 1.0);

vector<double> potential_params;

const double mass = 2000.0;
double init_x = -5.0;
double sigma_x = 0.5; 
double init_px = 20.0;
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
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("seed", po::value<int>(&seed), "random seed")
        ("output_mod", po::value<string>(&output_mod), "output mode, init_s or init_px")
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
    // state = x, y, vx, vy, a00, a01, a10, a11, s
    state_t state(8, matrixop::ZEROZ);
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
    // assign job
    vector<int> my_jobs = MPIer::assign_job(Ntraj);
    int my_Ntraj = my_jobs.size();
    vector<state_t> state(my_Ntraj);
    // propagation variables
    runge_kutta4<state_t> rk4;
    // initialize
    for (int itraj(0); itraj < my_Ntraj; ++itraj) {
        state[itraj] = init_state();
    }
    // statistics
    double n0trans, n0refl, n1trans, n1refl;
    double vxtrans, vytrans, vxrefl, vyrefl;
    double KE, PE;
    // recorders
    int Nrec = Nstep / output_step;
    vector<double> n0trans_arr(Nrec), n0refl_arr(Nrec), n1trans_arr(Nrec), n1refl_arr(Nrec);
    vector<double> vxtrans_arr(Nrec), vxrefl_arr(Nrec), vytrans_arr(Nrec), vyrefl_arr(Nrec);
    vector<double> KE_arr(Nrec), PE_arr(Nrec);
    // main loop
    vector< vector< complex<double> > > lastevt_save(my_Ntraj);
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < my_Ntraj; ++itraj) {
            // assign last evt
            lastevt = move(lastevt_save[itraj]);
            // calc info
            cal_info_nume(
                    vector<double> { state[itraj][0].real(), state[itraj][1].real() },
                    Fx, Fy, dcx, dcy, eva, lastevt
                    );
            // propagate
            rk4.do_step(sys, state[itraj], istep * dt, dt);
            // save last evt
            lastevt_save[itraj] = move(lastevt);
        }

        // data analysis
        if (istep % output_step == 0) {
            int irec = istep / output_step;
            n0trans = n0refl = n1trans = n1refl = 0.0;
            vxtrans = vxrefl = vytrans = vyrefl = 0.0;
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0trans, &n0refl, &n1trans, &n1refl, 
                     &vxtrans, &vxrefl, &vytrans, &vyrefl,
                     &KE, &PE] (const state_t& st) { 
                        if (st[2].real() > 0.0) {
                            n0trans += st[4].real();
                            n1trans += st[7].real();

                            vxtrans += st[2].real();
                            vytrans += st[3].real();
                        }
                        else {
                            n0refl += st[4].real();
                            n1refl += st[7].real();

                            vxrefl += st[2].real();
                            vyrefl += st[3].real();
                        }
                        KE += 0.5 * mass * (pow(st[2].real(), 2) + pow(st[3].real(), 2));
                        matrixop::hdiag(cal_H(vector<double> { st[0].real(), st[1].real() }), eva);
                        PE += st[4].real() * eva[0] + st[7].real() * eva[1];
                    });
            n0trans_arr[irec] = n0trans;
            n0refl_arr[irec] = n0refl;
            n1trans_arr[irec] = n1trans;
            n1refl_arr[irec] = n1refl;
            vxtrans_arr[irec] = vxtrans;
            vxrefl_arr[irec] = vxrefl;
            vytrans_arr[irec] = vytrans;
            vyrefl_arr[irec] = vyrefl;
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

                fill(vxtrans_arr.begin() + irec + 1, vxtrans_arr.end(), vxtrans);
                fill(vxrefl_arr.begin() + irec + 1, vxrefl_arr.end(), vxrefl);
                fill(vytrans_arr.begin() + irec + 1, vytrans_arr.end(), vytrans);
                fill(vyrefl_arr.begin() + irec + 1, vyrefl_arr.end(), vyrefl);

                fill(KE_arr.begin() + irec + 1, KE_arr.end(), KE);
                fill(PE_arr.begin() + irec + 1, PE_arr.end(), PE);
                break;
            }
        }
    }
    // collect data
    MPIer::barrier();
    for (int r = 1; r < MPIer::size; ++r) {
        if (MPIer::rank == r) {
            MPIer::send(0, 
                    n0trans_arr, n0refl_arr, n1trans_arr, n1refl_arr, 
                    vxtrans_arr, vxrefl_arr, vytrans_arr, vyrefl_arr, 
                    KE_arr, PE_arr
                    );
        }
        else if (MPIer::master) {
            vector<double> buf;

            MPIer::recv(r, buf); n0trans_arr += buf;
            MPIer::recv(r, buf); n0refl_arr += buf;
            MPIer::recv(r, buf); n1trans_arr += buf;
            MPIer::recv(r, buf); n1refl_arr += buf;

            MPIer::recv(r, buf); vxtrans_arr += buf;
            MPIer::recv(r, buf); vxrefl_arr += buf;
            MPIer::recv(r, buf); vytrans_arr += buf;
            MPIer::recv(r, buf); vyrefl_arr += buf;

            MPIer::recv(r, buf); KE_arr += buf;
            MPIer::recv(r, buf); PE_arr += buf;
        }
        MPIer::barrier();
    }

    // output
    if (MPIer::master) {
        // para & header
        output_potential_param();
        ioer::info("# Ehrenfest para: ", " MPIsize = ", MPIer::size, 
                    " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, 
                    " mass = ", mass, 
                    " init_x = ", init_x, " init_px = ", init_px, 
                    " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                    " init_y = ", init_y, " init_py = ", init_py, 
                    " sigma_y = ", sigma_y, " sigma_py = ", sigma_py, 
                    " init_s = ", init_s
                );
        ioer::tabout( "#", "t", "n0trans", "n0refl", "n1trans", "n1refl", "ntrans", "nrefl", "pxtrans", "pytrans", "pxrefl", "pyrefl", "Etot", "");

        for (int irec = 0; irec < Nrec; ++irec) {
            n0trans = n0trans_arr[irec];
            n0refl = n0refl_arr[irec];
            n1trans = n1trans_arr[irec];
            n1refl = n1refl_arr[irec];

            vxtrans = vxtrans_arr[irec];
            vxrefl = vxrefl_arr[irec];
            vytrans = vytrans_arr[irec];
            vyrefl = vyrefl_arr[irec];

            KE = KE_arr[irec];
            PE = PE_arr[irec];

            if (n0trans + n1trans > 0.0) {
                vxtrans /= (n0trans + n1trans);
                vytrans /= (n0trans + n1trans);
            }
            if (n0refl + n1refl > 0.0) {
                vxrefl /= (n0refl + n1refl);
                vyrefl /= (n0refl + n1refl);
            }

            n0trans /= Ntraj;
            n0refl /= Ntraj;
            n1trans /= Ntraj;
            n1refl /= Ntraj;

            KE /= Ntraj;
            PE /= Ntraj;

            ioer::tabout('#', irec * output_step * dt, 
                    n0trans, n0refl, n1trans, n1refl, 
                    n0trans + n1trans, n0refl + n1refl, 
                    vxtrans * mass, vytrans * mass, 
                    vxrefl * mass, vyrefl * mass, 
                    KE + PE);
        }
        // final results
            ioer::tabout(init_px, 
                    n0trans, n0refl, n1trans, n1refl, 
                    n0trans + n1trans, n0refl + n1refl, 
                    vxtrans * mass, vytrans * mass, 
                    vxrefl * mass, vyrefl * mass, 
                    KE + PE);
    }
    MPIer::barrier();
}

int main(int argc, char** argv) {
    MPIer::setup();
    if (argc < 2) {
        if (MPIer::master) ioer::info("use --help for detailed info");
    }
    else {
        if (argparse(argc, argv) == false) {
            return 0;
        }
        randomer::seed(MPIer::assign_random_seed(seed));
        if (MPIer::master) timer::tic();
        ehrenfest();
        if (MPIer::master) ioer::info("# ", timer::toc());
    }
    MPIer::barrier();
    MPIer::finalize();
    return 0;
}
