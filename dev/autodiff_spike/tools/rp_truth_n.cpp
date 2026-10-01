// REFPROP-only reference answer for every state of an N-component mixture (no CoolProp code).
//
//   ./rp_truth_n 'A.FLD|B.FLD|...' 'z1,z2,...' gerg|default  < states  > truth
//
// stdin: "idx T[K] p[Pa]" per line.  stdout: "idx phase rho_ref beta tm status" where
//   phase   1ph | 2ph (NOREF when no single-phase root reproduces p)
//   rho_ref bulk molar density [mol/m3] of the reference answer (-1 when the split did not converge)
//   beta    phase fraction of the incipient-side phase (0 for 1ph)
//   tm      minimum tangent-plane distance over all trials (tm < -1e-8: a split exists)
//   status  ok | loose | flashfail
// Method:
//   single phase: TPRHO liquid and vapor roots that reproduce p; the lower g = sum z ln f wins;
//   stability:    N near-pure + 20 random trial compositions, both roots, SS to 1e-12;
//   split:        SS on K (each phase on its lower-g root) with Rachford-Rice, seeded from the
//                 lowest-tm trial (K_i = w_i / z_i).
#include <dlfcn.h>
#include <algorithm>
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
static TPRHO_t TPRHO;
static FG_t FG;
static PR_t PR;
static int N = 0;

// ln f_i at (T, p, x) on root kph (1 liquid, 2 vapor); false when that root does not reproduce p
static bool lnf_root(double T, double pk, const std::vector<double>& x, int kph, double* D, double* lnf) {
    double xx[20] = {}, d = 0;
    for (int i = 0; i < N; ++i)
        xx[i] = x[i];
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
    for (int i = 0; i < N; ++i)
        lnf[i] = x[i] > 0 ? std::log(f[i]) : -1e300;
    *D = d;
    return true;
}
// lower-g of the two roots
static bool best_root(double T, double pk, const std::vector<double>& x, double* D, double* lnf) {
    double gb = 1e300;
    bool ok = false;
    for (int kph = 1; kph <= 2; ++kph) {
        double d, lf[20];
        if (!lnf_root(T, pk, x, kph, &d, lf)) continue;
        double g = 0;
        for (int i = 0; i < N; ++i)
            if (x[i] > 0) g += x[i] * lf[i];
        if (g < gb) {
            gb = g;
            *D = d;
            std::copy(lf, lf + N, lnf);
            ok = true;
        }
    }
    return ok;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: rp_truth_n 'A.FLD|B.FLD' 'z1,z2' gerg|default < states\n");
        return 1;
    }
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
    int j = std::string(argv[3]) == "gerg" ? 1 : 0, k = 0, ie = 0;
    FL(h0, &j, &k, &ie, e0, 255, 255);
    std::string s = argv[1];
    N = 1 + static_cast<int>(std::count(s.begin(), s.end(), '|'));
    std::vector<double> z;
    {
        std::string zs = argv[2];
        std::size_t q = 0;
        while (q < zs.size()) {
            std::size_t e = zs.find(',', q);
            if (e == std::string::npos) e = zs.size();
            z.push_back(std::stod(zs.substr(q, e - q)));
            q = e + 1;
        }
        double sum = 0;
        for (double v : z)
            sum += v;
        for (double& v : z)
            v /= sum;
    }
    std::vector<char> hf(10000, 0), he(255, 0);
    std::memcpy(hf.data(), s.data(), s.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = N;
    SU(&n, hf.data(), hmx, hrf, &ie, he.data(), 10000, 255, 3, 255);

    int idx;
    double T, p;
    while (std::scanf("%d %lf %lf", &idx, &T, &p) == 3) {
        const double pk = p / 1000;
        double Dz, lz[20];
        if (!best_root(T, pk, z, &Dz, lz)) {
            std::printf("%d NOREF\n", idx);
            std::fflush(stdout);
            continue;
        }
        // stability
        std::vector<std::vector<double>> trials;
        for (int i = 0; i < N; ++i) {
            std::vector<double> w(N, 0.001 / std::max(N - 1, 1));
            w[i] = 0.999;
            trials.push_back(w);
        }
        std::mt19937_64 g(idx);
        std::gamma_distribution<double> G(1.0);
        for (int t = 0; t < 20; ++t) {
            std::vector<double> w(N);
            double S = 0;
            for (double& v : w)
                S += (v = G(g));
            for (double& v : w)
                v /= S;
            trials.push_back(w);
        }
        double tmin = 1e300;
        std::vector<double> wbest;
        for (const auto& w0 : trials)
            for (int kph = 1; kph <= 2; ++kph) {
                std::vector<double> w = w0;
                double lw[20], D, tm = 1e300;
                for (int it = 0; it < 300; ++it) {
                    if (!lnf_root(T, pk, w, kph, &D, lw)) break;
                    tm = 0;
                    for (int i = 0; i < N; ++i)
                        if (w[i] > 0) tm += w[i] * (lw[i] - lz[i]);
                    if (!std::isfinite(tm)) {
                        tm = 1e300;
                        break;
                    }
                    double W[20], S = 0, dm = 0;
                    for (int i = 0; i < N; ++i) {
                        W[i] = w[i] > 1e-300 ? w[i] * std::exp(lz[i] - lw[i]) : 0.0;
                        S += W[i];
                    }
                    for (int i = 0; i < N; ++i) {
                        const double nw = W[i] / S;
                        dm = std::max(dm, std::fabs(nw - w[i]));
                        w[i] = nw;
                    }
                    if (dm < 1e-12) break;
                }
                if (tm < tmin) {
                    tmin = tm;
                    wbest = w;
                }
            }
        if (!(tmin < -1e-8)) {
            std::printf("%d 1ph %.12g 0 %.3e ok\n", idx, Dz * 1000, tmin);
            std::fflush(stdout);
            continue;
        }
        // split: y = incipient side, x = feed side; K = y/x
        std::vector<double> K(N), x(N), y(N);
        for (int i = 0; i < N; ++i)
            K[i] = z[i] > 0 ? std::max(wbest[i], 1e-300) / z[i] : 1.0;
        double beta = 0.5, DL = 0, DV = 0, res = 1;
        bool ok = true;
        for (int it = 0; it < 5000 && ok; ++it) {
            double Kmin = 1e300, Kmax = 0;
            for (double v : K) {
                Kmin = std::min(Kmin, v);
                Kmax = std::max(Kmax, v);
            }
            if (!(Kmin < 1 && Kmax > 1)) {
                ok = false;
                break;
            }
            double lo = 1.0 / (1.0 - Kmax) + 1e-15, hi = 1.0 / (1.0 - Kmin) - 1e-15;
            lo = std::max(lo, 0.0);
            hi = std::min(hi, 1.0);
            for (int b = 0; b < 200; ++b) {
                const double m = 0.5 * (lo + hi);
                double f = 0;
                for (int i = 0; i < N; ++i)
                    f += z[i] * (K[i] - 1) / (1 + m * (K[i] - 1));
                (f > 0 ? lo : hi) = m;
            }
            beta = 0.5 * (lo + hi);
            double sx = 0, sy = 0;
            for (int i = 0; i < N; ++i) {
                x[i] = z[i] / (1 + beta * (K[i] - 1));
                y[i] = K[i] * x[i];
                sx += x[i];
                sy += y[i];
            }
            for (int i = 0; i < N; ++i) {
                x[i] /= sx;
                y[i] /= sy;
            }
            double lfx[20], lfy[20];
            if (!best_root(T, pk, x, &DL, lfx) || !best_root(T, pk, y, &DV, lfy)) {
                ok = false;
                break;
            }
            res = 0;
            for (int i = 0; i < N; ++i) {
                if (!(z[i] > 0)) continue;
                // K_i = phi_x / phi_y = (f_x/x) / (f_y/y)
                const double lK = (lfx[i] - std::log(x[i])) - (lfy[i] - std::log(y[i]));
                res = std::max(res, std::fabs(lK - std::log(K[i])));
                K[i] = std::exp(lK);
            }
            if (it > 3 && res < 1e-12) break;
        }
        if (!ok || !(beta > 0 && beta < 1)) {
            std::printf("%d 2ph -1 -1 %.3e flashfail\n", idx, tmin);
        } else {
            std::printf("%d 2ph %.12g %.10f %.3e %s\n", idx, 1000.0 / (beta / DV + (1 - beta) / DL), beta, tmin, res < 1e-9 ? "ok" : "loose");
        }
        std::fflush(stdout);
    }
}
