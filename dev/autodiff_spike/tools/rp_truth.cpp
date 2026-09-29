// REFPROP-only reference answer (default model = Gernert for CO2/water) for every state of a binary.
// argv: z_CO2.  Reads "i T p[Pa]"; prints "i phase rho_ref beta tm status".
//   single-phase root: lower g = sum z ln f among the TPRHO roots that reproduce p;
//   stability: 2 near-pure + 20 random trials, both roots, SS to 1e-12 (as rp_stability.cpp);
//   split (if tm < -1e-8): SS on K, each phase on its lower-g root, seeded from the best trial and z.
#include <dlfcn.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
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
bool lnf_root(double T, double pk, const double* x, int kph, double* D, double* lnf) {
    double xx[20] = {x[0], x[1]}, d = 0;
    int kg = 0, ie = 0;
    char he[256];
    TPRHO(&T, &pk, xx, &kph, &kg, &d, &ie, he, 255);
    if (ie > 0 || !(d > 0)) return false;
    double pc;
    PR(&T, &d, xx, &pc);
    if (std::fabs(pc / pk - 1) > 1e-8) return false;
    double f[20];
    FG(&T, &d, xx, f, &ie, he, 255);
    if (ie > 0) return false;
    for (int i = 0; i < 2; i++)
        lnf[i] = std::log(f[i]);
    *D = d;
    return true;
}
// lower-g root
bool best_root(double T, double pk, const double* x, double* D, double* lnf) {
    double gb = 1e300;
    bool ok = false;
    for (int kph = 1; kph <= 2; kph++) {
        double d, lf[2];
        if (!lnf_root(T, pk, x, kph, &d, lf)) continue;
        double g = 0;
        for (int i = 0; i < 2; i++)
            if (x[i] > 0) g += x[i] * lf[i];
        if (g < gb) {
            gb = g;
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
    double T, p;
    while (scanf("%d %lf %lf", &idx, &T, &p) == 3) {
        const double pk = p / 1000;
        double Dz, lz[2];
        if (!best_root(T, pk, z, &Dz, lz)) {
            printf("%d NOREF\n", idx);
            continue;
        }
        // stability
        std::vector<std::vector<double>> trials = {{0.999, 0.001}, {0.001, 0.999}};
        std::mt19937_64 g(idx);
        std::uniform_real_distribution<double> U(0, 1);
        for (int t = 0; t < 20; t++) {
            double a = U(g);
            trials.push_back({a, 1 - a});
        }
        double tmin = 1e300, wbest[2] = {0, 0};
        for (auto w0 : trials)
            for (int kph = 1; kph <= 2; kph++) {
                double w[2] = {w0[0], w0[1]}, lw[2], D, tm = 1e300;
                for (int it = 0; it < 300; it++) {
                    if (!lnf_root(T, pk, w, kph, &D, lw)) break;
                    tm = 0;
                    for (int i = 0; i < 2; i++)
                        if (w[i] > 0) tm += w[i] * (lw[i] - lz[i]);
                    if (!std::isfinite(tm)) {
                        tm = 1e300;
                        break;
                    }
                    double W[2], S = 0, dm = 0;
                    for (int i = 0; i < 2; i++) {
                        W[i] = w[i] > 1e-300 ? w[i] * std::exp(lz[i] - lw[i]) : 0;
                        S += W[i];
                    }
                    for (int i = 0; i < 2; i++) {
                        double nw = W[i] / S;
                        dm = std::max(dm, std::fabs(nw - w[i]));
                        w[i] = nw;
                    }
                    if (dm < 1e-12) break;
                }
                if (tm < tmin) {
                    tmin = tm;
                    wbest[0] = w[0];
                    wbest[1] = w[1];
                }
            }
        if (!(tmin < -1e-8)) {
            printf("%d 1ph %.12g 0 %.3e ok\n", idx, Dz * 1000, tmin);
            continue;
        }
        // split: phase A = trial composition, phase B = z-side; x = water-richer
        double a[2] = {wbest[0], wbest[1]}, b[2] = {z[0], z[1]};
        if (std::fabs(a[0] - b[0]) < 1e-6) {
            b[0] = z[0] - 0.5 * (a[0] - z[0]);
            b[1] = 1 - b[0];
        }
        double *x = a[0] < b[0] ? a : b, *y = a[0] < b[0] ? b : a;  // x: water-rich (lower CO2)
        double K[2], DL = 0, DV = 0, lfL[2], lfV[2], beta = -1;
        bool ok = true;
        double res = 1;
        for (int it = 0; it < 2000; it++) {
            if (!best_root(T, pk, x, &DL, lfL) || !best_root(T, pk, y, &DV, lfV)) {
                ok = false;
                break;
            }
            res = 0;
            for (int i = 0; i < 2; i++) {
                double Ki = std::exp(lfL[i] - std::log(x[i]) - lfV[i] + std::log(y[i]));
                if (it) res = std::max(res, std::fabs(std::log(Ki / K[i])));
                K[i] = Ki;
            }
            if (!((K[0] - 1) * (K[1] - 1) < 0)) {
                ok = false;
                break;
            }
            x[0] = (1 - K[1]) / (K[0] - K[1]);
            x[1] = 1 - x[0];
            y[0] = K[0] * x[0];
            y[1] = 1 - y[0];
            if (!(x[0] > 0 && x[0] < 1 && y[0] > 0 && y[0] < 1)) {
                ok = false;
                break;
            }
            if (it > 5 && res < 1e-12) break;
        }
        if (ok) beta = (z[0] - x[0]) / (y[0] - x[0]);
        if (!ok || !(beta > 0 && beta < 1)) {
            printf("%d 2ph -1 -1 %.3e flashfail\n", idx, tmin);
            continue;
        }
        printf("%d 2ph %.12g %.8f %.3e %s x=%.8f y=%.8f DL=%.8g DV=%.8g\n", idx, 1000.0 / (beta / DV + (1 - beta) / DL), beta, tmin,
               res < 1e-9 ? "ok" : "loose", x[0], y[0], DL * 1000, DV * 1000);
    }
}
