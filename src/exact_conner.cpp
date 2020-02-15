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
#include "boost/program_options.hpp"

#include "2d_conner_potential.hpp"

using namespace std;
namespace po = boost::program_options;

double L = 16;
int M = 256;

vector<double> potential_params;
double mass = 1000.0;

double xI = -3.0;
double yI = 0.0;
double sigmax = 1.0;
double sigmay = 1.0;
double kxI = 20.0;
double kyI = 0.0;

double init_s = 0.0;
int Nstep = 2000;
double dt = 0.1;

int output_step = 100;



inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("L", po::value<double>(&L), "grid para, the grid is [-L/2, L/2]")
        ("M", po::value<int>(&M), "grid number")
        ("mass", po::value<double>(&mass), "mass")
        ("init_x", po::value<double>(&xI), "init x")
        ("init_y", po::value<double>(&yI), "init y")
        ("sigma_x", po::value<double>(&sigmax), "init sigma x")
        ("sigma_y", po::value<double>(&sigmay), "init sigma y")
        ("init_px", po::value<double>(&kxI), "init px")
        ("init_py", po::value<double>(&kyI), "init py")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("dt", po::value<double>(&dt), "time step")
        ("output_step", po::value<int>(&output_step), "output step")
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
    return false;
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
    vector< complex<double> > evts00(M*M), evts01(M*M), evts10(M*M), evts11(M*M);
    vector< complex<double> > H00(M*M), H01(M*M), H10(M*M), H11(M*M);
    for (int k = 0; k < M; ++k) {
        double y = yarr[k];
        for (int j = 0; j < M; ++j) {
            double x = xarr[j];
            vector< complex<double> > H = cal_H(vector<double> {x, y});
            vector<double> eva;
            vector< complex<double> > evt, evamat;
            matrixop::eigh(H, eva, evt);
            
            // propagation matrix
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
    vector< complex<double> > psi0(M * M, 0.0), psi1(M * M, 0.0);
    // initialize wf on diabats
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
    double n0 = 0.0, n1 = 0.0;
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
            n0 = n1 = 0.0;
            for (int k = 0; k < M; ++k) {
                for (int j = 0; j < M; ++j) {
                    n0 += pow(abs(psi_k0[k+j*M]), 2);
                    n1 += pow(abs(psi_k1[k+j*M]), 2);
                    KE0 += pow(abs(psi_k0[k+j*M]), 2) * (kxarr[j]*kxarr[j] + kyarr[k]*kyarr[k]) / 2.0 / mass;
                    KE1 += pow(abs(psi_k1[k+j*M]), 2) * (kxarr[j]*kxarr[j] + kyarr[k]*kyarr[k]) / 2.0 / mass;
                }
            }
            PE = real(sum(conj(psi0) * H00 * psi0 + conj(psi0) * H01 * psi1 + conj(psi1) * H10 * psi0 + conj(psi1) * H11 * psi1));


            // output
            if (istep == 0) {
                output_potential_param();
                ioer::info("# EXACT 2D: ");
                ioer::info("# DIAB ");
                ioer::info("# para: ", " L = ", L, " M = ", M, " mass = ", mass, 
                                       " xI = ", xI, " yI = ", yI, " sigmax = ", sigmax, " sigmay = ", sigmay, " kxI = ", kxI, " kyI = ", kyI, " init_s = ", init_s, " c0 = ", c0, " c1 = ", c1,
                                       " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step);
                ioer::info("# dx = ", dx, " dy = ", dy, " dkx = ", dkx, " dky = ", dky);
                ioer::tabout('#', "t", "n0", "n1", "Etot");
            }
            ioer::tabout('#', istep * dt, n0, n1, KE0 + KE1 + PE);
            // check end
            if (check_end(psi0, psi1, xarr, yarr) == true) {
                ioer::info("# check_end returns true");
                break;
            }
        }
    }
    // final output
    ioer::tabout_nonewline(kxI);
    ioer::tabout(n0, n1, KE0 + KE1 + PE);
}


void test() {
    vector<double> r(3, 0.0);
    vector<double> eva;
    for (double x(-8.0); x < 8.0; x += 0.2) {
        for (double y(-8.0); y < 8.0; y += 0.2) {
            r[0] = x;
            r[1] = y;
            vector< complex<double> > H = cal_H(vector<double> {x, y});
            matrixop::eigh(H, eva);

            ioer::tabout(x, y, 
                    eva[0], eva[1], 
                    H[0+0*2].real(), H[1+1*2].real()
                    );
        }
    }
}


int main(int argc, char** argv) {
    /*
    test();
    return 0;
    */
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    exact();
    ioer::info("# ", timer::toc());
    return 0;
}
