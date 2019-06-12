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

#include "2d_flat_potential_phixplusy.hpp"

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

double xwall_left = -99999999.9;
double xwall_right = 99999999.9;

double init_s = 0.0;
int Nstep = 2000;
double dt = 0.1;

int output_step = 100;

bool enable_adiab = false;

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
        ("xwall_left", po::value<double>(&xwall_left), " left boundary x value to check end")
        ("xwall_right", po::value<double>(&xwall_right), " right boundary x value to check end")
        ("potential_params", po::value< vector<double> >(&potential_params)->multitoken(), "potential_params vector")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("dt", po::value<double>(&dt), "time step")
        ("enable_adiab", po::value<bool>(&enable_adiab), "output observables on adiabats")
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

            if (enable_adiab) {
                // adiab: phase align & record evts
                if (j == 0 and k == 0) {
                    const complex<double> phase0 = evt[0+0*2] / abs(evt[0+0*2]);
                    evt[0+0*2] *= conj(phase0);
                    evt[1+0*2] *= conj(phase0);
                    const complex<double> phase1 = evt[1+1*2] / abs(evt[1+1*2]);
                    evt[0+1*2] *= conj(phase1);
                    evt[1+1*2] *= conj(phase1);
                }
                else if (k == 0) {
                    const complex<double> phase0 = conj(evts00[k+(j-1)*M]) * evt[0+0*2] + conj(evts10[k+(j-1)*M]) * evt[1+0*2];
                    evt[0+0*2] *= conj(phase0);
                    evt[1+0*2] *= conj(phase0);
                    const complex<double> phase1 = conj(evts01[k+(j-1)*M]) * evt[0+1*2] + conj(evts11[k+(j-1)*M]) * evt[1+1*2];
                    evt[0+1*2] *= conj(phase1);
                    evt[1+1*2] *= conj(phase1);
                }
                else {
                    const complex<double> phase0 = conj(evts00[(k-1)+j*M]) * evt[0+0*2] + conj(evts10[(k-1)+j*M]) * evt[1+0*2];
                    evt[0+0*2] *= conj(phase0);
                    evt[1+0*2] *= conj(phase0);
                    const complex<double> phase1 = conj(evts01[(k-1)+j*M]) * evt[0+1*2] + conj(evts11[(k-1)+j*M]) * evt[1+1*2];
                    evt[0+1*2] *= conj(phase1);
                    evt[1+1*2] *= conj(phase1);
                }
                evts00[k+j*M] = evt[0+0*2];
                evts01[k+j*M] = evt[0+1*2];
                evts10[k+j*M] = evt[1+0*2];
                evts11[k+j*M] = evt[1+1*2];
            }
        }
    }
    // initialized WF on adiabats
    vector< complex<double> > psiad0(M * M, 0.0), psiad1(M * M, 0.0);
    vector< complex<double> > psi0(M * M, 0.0), psi1(M * M, 0.0);
    if (enable_adiab) {
        // initialize wf on adiabats
        for (int k = 0; k < M; ++k) {
            double y = yarr[k];
            for (int j = 0; j < M; ++j) {
                double x = xarr[j];
                psiad0[k+j*M] = c0 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
                psiad1[k+j*M] = c1 * exp(matrixop::IMAGIZ * (kxI * x + kyI * y)) * exp(-pow((x - xI) / sigmax, 2) - pow((y - yI) / sigmay, 2));
            }
        }
        double nm = norm(psiad0 | psiad1);
        psiad0 /= nm;
        psiad1 /= nm;
        psi0 = evts00 * psiad0 + evts01 * psiad1;
        psi1 = evts10 * psiad0 + evts11 * psiad1;
    }
    else {
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
    }

    // covinience vairables
    vector<int> dim{ M, M };
    // statistics
    double KE0 = 0.0, KE1 = 0.0, PE = 0.0;
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double py0trans = 0.0, py0refl = 0.0, py1trans = 0.0, py1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    // adiab statistics
    double nad0trans = 0.0, nad0refl = 0.0, nad1trans = 0.0, nad1refl = 0.0;
    double pady0trans = 0.0, pady0refl = 0.0, pady1trans = 0.0, pady1refl = 0.0;
    double padx0trans = 0.0, padx0refl = 0.0, padx1trans = 0.0, padx1refl = 0.0;
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

            if (enable_adiab) {
                // analysis in adiab

                // get psiad
                psiad0 = conj(evts00) * psi0 + conj(evts10) * psi1;
                psiad1 = conj(evts01) * psi0 + conj(evts11) * psi1;
                double nm = norm(psiad0 | psiad1);
                psiad0 /= nm;
                psiad1 /= nm;

                // get psiad fft
                auto psiad_k0 = myfftshift(misc::fftn(psiad0, dim));
                auto psiad_k1 = myfftshift(misc::fftn(psiad1, dim));
                nm = norm(psiad_k0 | psiad_k1);
                psiad_k0 /= nm;
                psiad_k1 /= nm;

                // analysis
                nad0trans = nad0refl = nad1trans = nad1refl = 0.0;
                pady0trans = pady0refl = pady1trans = pady1refl = 0.0;
                padx0trans = padx0refl = padx1trans = padx1refl = 0.0;
                for (int k = 0; k < M; ++k) {
                    for (int j = 0; j < M; ++j) {
                        if (kxarr[j] >= 0.0) {
                            nad0trans += pow(abs(psiad_k0[k+j*M]), 2);
                            nad1trans += pow(abs(psiad_k1[k+j*M]), 2);
                            pady0trans += pow(abs(psiad_k0[k+j*M]), 2) * kyarr[k];
                            pady1trans += pow(abs(psiad_k1[k+j*M]), 2) * kyarr[k];
                            padx0trans += pow(abs(psiad_k0[k+j*M]), 2) * kxarr[j];
                            padx1trans += pow(abs(psiad_k1[k+j*M]), 2) * kxarr[j];
                        }
                        else {
                            nad0refl += pow(abs(psiad_k0[k+j*M]), 2);
                            nad1refl += pow(abs(psiad_k1[k+j*M]), 2);
                            pady0refl += pow(abs(psiad_k0[k+j*M]), 2) * kyarr[k];
                            pady1refl += pow(abs(psiad_k1[k+j*M]), 2) * kyarr[k];
                            padx0refl += pow(abs(psiad_k0[k+j*M]), 2) * kxarr[j];
                            padx1refl += pow(abs(psiad_k1[k+j*M]), 2) * kxarr[j];
                        }
                    }
                }
                pady0trans /= (nad0trans + 1e-16);
                pady1trans /= (nad1trans + 1e-16);
                pady0refl /= (nad0refl + 1e-16);
                pady1refl /= (nad1refl + 1e-16);
                padx0trans /= (nad0trans + 1e-16);
                padx1trans /= (nad1trans + 1e-16);
                padx0refl /= (nad0refl + 1e-16);
                padx1refl /= (nad1refl + 1e-16);
            }

            // output
            if (istep == 0) {
                output_potential_param();
                ioer::info("# EXACT 2D: ");
                if (enable_adiab) {
                    ioer::info("# ADIAB ");
                }
                else {
                    ioer::info("# DIAB ");
                }
                ioer::info("# para: ", " L = ", L, " M = ", M, " mass = ", mass, 
                                       " xI = ", xI, " yI = ", yI, " sigmax = ", sigmax, " sigmay = ", sigmay, " kxI = ", kxI, " kyI = ", kyI, " init_s = ", init_s, " c0 = ", c0, " c1 = ", c1,
                                       " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step);
                ioer::info("# dx = ", dx, " dy = ", dy, " dkx = ", dkx, " dky = ", dky, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right);
                ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "py0trans", "px0refl", "py0refl", "px1trans", "py1trans", "px1refl", "py1refl", "Etot");
            }
            if (enable_adiab) {
                ioer::tabout('#', istep * dt, nad0trans, nad0refl, nad1trans, nad1refl, padx0trans, pady0trans, padx0refl, pady0refl, padx1trans, pady1trans, padx1refl, pady1refl, KE0 + KE1 + PE);
            }
            else {
                ioer::tabout('#', istep * dt, n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl, KE0 + KE1 + PE);
            }
            // check end
            if (check_end(psi0, psi1, xarr, yarr) == true) {
                ioer::info("# check_end returns true");
                break;
            }
        }
    }
    // final output
    ioer::tabout_nonewline(kxI);
    if (enable_adiab) {
        ioer::tabout(nad0trans, nad0refl, nad1trans, nad1refl, padx0trans, pady0trans, padx0refl, pady0refl, padx1trans, pady1trans, padx1refl, pady1refl, KE0 + KE1 + PE);
    }
    else {
        ioer::tabout(n0trans, n0refl, n1trans, n1refl, px0trans, py0trans, px0refl, py0refl, px1trans, py1trans, px1refl, py1refl, KE0 + KE1 + PE);
    }
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
