// Brute-force tangent-plane stability with REFPROP only (GERG mode). Reads idx\tT[K]\tp[Pa]; prints idx T p min_tm.
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
int N = 5;
// ln f_i (f in kPa) at (T,p,w) for root kph; returns false if no valid root at that pressure
bool lnf(double T, double pk, const double* w, int kph, double* out, double* Dout) {
    double x[20] = {};
    for (int i = 0; i < N; i++)
        x[i] = w[i];
    double D = 0;
    int kg = 0, ie = 0;
    char he[256];
    TPRHO(&T, &pk, x, &kph, &kg, &D, &ie, he, 255);
    if (ie > 0 || !(D > 0)) return false;
    double pc;
    PR(&T, &D, x, &pc);
    if (std::fabs(pc / pk - 1) > 1e-7) return false;
    double f[20];
    FG(&T, &D, x, f, &ie, he, 255);
    if (ie > 0) return false;
    for (int i = 0; i < N; i++)
        out[i] = std::log(f[i]);
    *Dout = D;
    return true;
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
    int j = (argc > 3 && std::string(argv[3]) == "default") ? 0 : 1, k = 0, ie = 0;
    FL(h0, &j, &k, &ie, e0, 255, 255);
    // argv: RPfile1|RPfile2|...  z1,z2,...
    std::string s = argv[1];
    N = 1;
    for (char c : s)
        if (c == '|') N++;
    std::vector<double> zv;
    {
        std::string zs = argv[2];
        size_t q = 0;
        while (q < zs.size()) {
            size_t e = zs.find(',', q);
            if (e == std::string::npos) e = zs.size();
            zv.push_back(std::stod(zs.substr(q, e - q)));
            q = e + 1;
        }
    }
    double zsum = 0;
    for (double v : zv)
        zsum += v;
    for (auto& v : zv)
        v /= zsum;
    std::vector<char> hf(10000, 0), he(255, 0);
    memcpy(hf.data(), s.data(), s.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = N;
    SU(&n, hf.data(), hmx, hrf, &ie, he.data(), 10000, 255, 3, 255);
    const double* z = zv.data();
    int idx;
    double T, p;
    while (scanf("%d %lf %lf", &idx, &T, &p) == 3) {
        double pk = p / 1000, lz[2][N], Dz[2];
        bool okz[2];
        for (int r = 0; r < 2; r++)
            okz[r] = lnf(T, pk, z, r + 1, lz[r], &Dz[r]);
        // reference = lower-G single phase: g/RT ~ sum z ln f
        int rz = -1;
        double gbest = 1e300;
        for (int r = 0; r < 2; r++)
            if (okz[r]) {
                double g = 0;
                for (int i = 0; i < N; i++)
                    g += z[i] * lz[r][i];
                if (g < gbest) {
                    gbest = g;
                    rz = r;
                }
            }
        if (rz < 0) {
            printf("%d %.10g %.10g NOREF\n", idx, T, p);
            continue;
        }
        std::vector<std::vector<double>> trials;
        for (int i = 0; i < N; i++) {
            std::vector<double> w(N, 0.02 / (N - 1));
            w[i] = 0.98;
            trials.push_back(w);
        }
        std::mt19937_64 g(idx);
        std::gamma_distribution<double> G(1.0);
        for (int t = 0; t < 20; t++) {
            std::vector<double> w(N);
            double S = 0;
            for (auto& v : w)
                S += (v = G(g));
            for (auto& v : w)
                v /= S;
            trials.push_back(w);
        }
        double tmin = 1e300;
        for (auto w0 : trials)
            for (int kph = 1; kph <= 2; kph++) {
                std::vector<double> w = w0;
                double lw[20], D, tm = 1e300;
                for (int it = 0; it < 200; it++) {
                    if (!lnf(T, pk, w.data(), kph, lw, &D)) break;
                    // tm = sum w_i (ln f_i(w) - ln f_i(z))  (dimensionless G/RT tangent-plane distance)
                    tm = 0;
                    for (int i = 0; i < N; i++)
                        if (w[i] > 0) tm += w[i] * (lw[i] - lz[rz][i]);
                    if (!std::isfinite(tm)) {
                        tm = 1e300;
                        break;
                    }
                    // SS: W_i = w_i * f_i(z)/f_i(w)
                    double W[20], S = 0, dmax = 0;
                    for (int i = 0; i < N; i++) {
                        W[i] = w[i] > 1e-300 ? w[i] * std::exp(lz[rz][i] - lw[i]) : 0.0;
                        S += W[i];
                    }
                    for (int i = 0; i < N; i++) {
                        double nw = W[i] / S;
                        dmax = std::max(dmax, std::fabs(nw - w[i]));
                        w[i] = nw;
                    }
                    if (dmax < 1e-12) break;
                }
                if (tm < tmin) tmin = tm;
            }
        printf("%d %.10g %.10g %.6e\n", idx, T, p, tmin);
        fflush(stdout);
    }
}
