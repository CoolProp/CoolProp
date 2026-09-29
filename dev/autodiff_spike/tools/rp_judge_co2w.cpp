// REFPROP-only (default = Gernert) judge for CO2/water density disagreements.
// Reads: i T p rho_off rho_on phase(1ph|2ph) [z1]; prints which build matches.
// 1ph: the root that reproduces p with lower g = sum z ln f wins.
// 2ph: SS + Rachford-Rice flash, each phase on its lower-g root, seeded CO2-rich / water-rich.
#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
using SETPATH_t = void (*)(char*, long);
using SETUP_t = void (*)(int*, char*, char*, char*, int*, char*, long, long, long, long);
using FLAGS_t = void (*)(char*, int*, int*, int*, char*, long, long);
using TPRHO_t = void (*)(double*, double*, double*, int*, int*, double*, int*, char*, long);
using FG_t = void (*)(double*, double*, double*, double*, int*, char*, long);
using PR_t = void (*)(double*, double*, double*, double*);
TPRHO_t TPRHO;
FG_t FG;
PR_t PR;
double gz(double T, double D, const double* x, double* lnf) {
    double xx[20] = {x[0], x[1]}, f[20];
    int ie = 0;
    char he[256];
    FG(&T, &D, xx, f, &ie, he, 255);
    double g = 0;
    for (int i = 0; i < 2; i++) {
        lnf[i] = std::log(f[i]);
        if (x[i] > 0) g += x[i] * lnf[i];
    }
    return g;
}
// lower-g root at (T,p,x); returns D (mol/L) and ln f
bool root(double T, double pk, const double* x, double* D, double* lnf) {
    double best = 1e300;
    bool ok = false;
    for (int kph = 1; kph <= 2; kph++) {
        double xx[20] = {x[0], x[1]}, d = 0;
        int kg = 0, ie = 0;
        char he[256];
        TPRHO(&T, &pk, xx, &kph, &kg, &d, &ie, he, 255);
        if (ie > 0 || !(d > 0)) continue;
        double pc;
        PR(&T, &d, xx, &pc);
        if (std::fabs(pc / pk - 1) > 1e-8) continue;
        double lf[2];
        double g = gz(T, d, x, lf);
        if (g < best) {
            best = g;
            *D = d;
            lnf[0] = lf[0];
            lnf[1] = lf[1];
            ok = true;
        }
    }
    return ok;
}
int main(int argc, char** argv) {
    void* h = dlopen("/Users/ianbell/REFPROP10/librefprop.dylib", RTLD_NOW);
    auto SP = (SETPATH_t)dlsym(h, "SETPATHdll");
    auto SU = (SETUP_t)dlsym(h, "SETUPdll");
    auto FL = (FLAGS_t)dlsym(h, "FLAGSdll");
    TPRHO = (TPRHO_t)dlsym(h, "TPRHOdll");
    FG = (FG_t)dlsym(h, "FGCTY2dll");
    PR = (PR_t)dlsym(h, "PRESSdll");
    char path[256] = "/Users/ianbell/REFPROP10/";
    SP(path, 255);
    char h0[256] = "GERG", e0[256] = {};
    int j = 0, k = 0, ie = 0;
    FL(h0, &j, &k, &ie, e0, 255, 255);
    std::string s = "CO2.FLD|WATER.FLD";
    std::vector<char> hf(10000, 0), he(255, 0);
    memcpy(hf.data(), s.data(), s.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = 2;
    SU(&n, hf.data(), hmx, hrf, &ie, he.data(), 10000, 255, 3, 255);
    const double z[2] = {atof(argv[1]), 1 - atof(argv[1])};
    int idx;
    double T, p, roff, ron;
    char ph[8];
    int on = 0, off = 0, neither = 0, fail = 0;
    while (scanf("%d %lf %lf %lf %lf %7s", &idx, &T, &p, &roff, &ron, ph) == 6) {
        double pk = p / 1000, ref = -1;
        if (!strcmp(ph, "1ph")) {
            // candidate roots are the two builds' densities: keep those reproducing p, pick lower g
            double best = 1e300;
            for (double r : {roff, ron}) {
                double D = r / 1000, pc, xx[20] = {z[0], z[1]}, lf[2];
                PR(&T, &D, xx, &pc);
                if (std::fabs(pc / pk - 1) > 1e-5) continue;
                double g = gz(T, D, z, lf);
                if (g < best - 1e-12) {
                    best = g;
                    ref = r;
                }
            }
        } else {
            double x[2] = {1e-3, 1 - 1e-3}, y[2] = {0.995, 0.005}, K[2], DL, DV, lfL[2], lfV[2], beta = 0.5;
            bool ok = true;
            for (int it = 0; it < 1000 && ok; it++) {
                if (!root(T, pk, x, &DL, lfL) || !root(T, pk, y, &DV, lfV)) {
                    ok = false;
                    break;
                }
                double res = 0;
                // SS on K: K_i = phiL_i / phiV_i = (fL_i/x_i)/(fV_i/y_i)
                for (int i = 0; i < 2; i++) {
                    double Ki = std::exp(lfL[i] - std::log(x[i]) - lfV[i] + std::log(y[i]));
                    if (it) res = std::max(res, std::fabs(std::log(Ki / K[i])));
                    K[i] = Ki;
                }
                if (!((K[0] - 1) * (K[1] - 1) < 0)) {
                    ok = false;
                    break;
                }
                // binary: compositions from K directly
                x[0] = (1 - K[1]) / (K[0] - K[1]);
                x[1] = 1 - x[0];
                y[0] = K[0] * x[0];
                y[1] = 1 - y[0];
                if (!(x[0] > 0 && x[0] < 1 && y[0] > 0 && y[0] < 1)) {
                    ok = false;
                    break;
                }
                beta = (z[0] - x[0]) / (y[0] - x[0]);
                if (it > 5 && res < 1e-13) break;
            }
            if (ok && beta > 0 && beta < 1) ref = 1000.0 / (beta / DV + (1 - beta) / DL);
        }
        if (ref < 0) {
            fail++;
            printf("%d FAIL %s\n", idx, ph);
            continue;
        }
        double eon = std::fabs(ron / ref - 1), eoff = std::fabs(roff / ref - 1);
        if (eon < 1e-6 && eoff >= 1e-6)
            on++;
        else if (eoff < 1e-6 && eon >= 1e-6)
            off++;
        else
            neither++;
        printf("%d %s ref=%.10g on_err=%.2e off_err=%.2e\n", idx, ph, ref, eon, eoff);
    }
    fprintf(stderr, "on right %d | off right %d | neither/both %d | judge failed %d\n", on, off, neither, fail);
}
