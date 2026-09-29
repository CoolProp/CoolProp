#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/Helmholtz/MixtureDerivatives.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>
#include <dlfcn.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include "AbstractState.h"

namespace {
using SETPATH_t = void (*)(char*, long);
using SETUP_t = void (*)(int*, char*, char*, char*, int*, char*, long, long, long, long);
using FLAGS_t = void (*)(char*, int*, int*, int*, char*, long, long);
using TPFLSH_t = void (*)(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                          double*, double*, int*, char*, long);
SETPATH_t RP_SETPATH = nullptr;
SETUP_t RP_SETUP = nullptr;
FLAGS_t RP_FLAGS = nullptr;
TPFLSH_t RP_TPFLSH = nullptr;
bool load_rp(const std::string& dir) {
    void* h = dlopen((dir + "librefprop.dylib").c_str(), RTLD_NOW);
    if (!h) return false;
    RP_SETPATH = reinterpret_cast<SETPATH_t>(dlsym(h, "SETPATHdll"));
    RP_SETUP = reinterpret_cast<SETUP_t>(dlsym(h, "SETUPdll"));
    RP_FLAGS = reinterpret_cast<FLAGS_t>(dlsym(h, "FLAGSdll"));
    RP_TPFLSH = reinterpret_cast<TPFLSH_t>(dlsym(h, "TPFLSHdll"));
    if (!RP_SETPATH || !RP_SETUP || !RP_FLAGS || !RP_TPFLSH) return false;
    char path[256] = {};
    std::strncpy(path, dir.c_str(), 255);
    RP_SETPATH(path, 255);
    return true;
}
const std::map<std::string, std::string> RPNAME = {
  {"Methane", "METHANE"}, {"Ethane", "ETHANE"},       {"Propane", "PROPANE"},   {"HydrogenSulfide", "H2S"}, {"Nitrogen", "NITROGEN"},
  {"Oxygen", "OXYGEN"},   {"Argon", "ARGON"},         {"CarbonDioxide", "CO2"}, {"Water", "WATER"},         {"IsoButane", "ISOBUTAN"},
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"},     {"R32", "R32"},
  {"R1234yf", "R1234YF"}, {"Helium", "HELIUM"},       {"Hydrogen", "HYDROGEN"}};
bool setup_rp(const std::vector<std::string>& fluids, bool gerg) {
    char hf0[256] = {}, herr0[256] = {};
    std::strcpy(hf0, "GERG");
    int j = gerg ? 1 : 0, k = 0, ie = 0;
    RP_FLAGS(hf0, &j, &k, &ie, herr0, 255, 255);
    std::string f;
    for (std::size_t i = 0; i < fluids.size(); ++i)
        f += (i ? "|" : "") + RPNAME.at(fluids[i]) + ".FLD";
    std::vector<char> hf(10000, '\0'), herr(255, '\0');
    std::memcpy(hf.data(), f.data(), f.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = static_cast<int>(fluids.size()), ierr = 0;
    RP_SETUP(&n, hf.data(), hmx, hrf, &ierr, herr.data(), 10000, 255, 3, 255);
    return ierr <= 0;
}
std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? "&" : "") + v[i];
    return s;
}
double now_us() {
    return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
double median(std::vector<double> v) {
    if (v.empty()) return NAN;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}
double mean(const std::vector<double>& v) {
    double s = 0;
    for (double x : v)
        s += x;
    return v.empty() ? NAN : s / v.size();
}

double gRT(CoolProp::HelmholtzEOSMixtureBackend& H) {
    double g = 0;
    auto x = H.get_mole_fractions();
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i] > 0) g += x[i] * (std::log(x[i]) + CoolProp::MixtureDerivatives::ln_fugacity_coefficient(H, i, CoolProp::XN_INDEPENDENT));
    return g;
}
}  // namespace
// stdin: "k T p rho_cp" lines; k selects the mixture.  Prints TPFLSH's split evaluated in CoolProp's GERG-2008.
int main(int argc, char** argv) {
    load_rp("/Users/ianbell/REFPROP10/");
    std::vector<std::vector<std::string>> F = {{},
                                               {"Nitrogen", "Methane", "Ethane", "n-Butane", "n-Pentane"},
                                               {"CarbonDioxide", "Hydrogen"},
                                               {"CarbonDioxide", "Nitrogen", "Oxygen", "Helium"},
                                               {"Nitrogen", "Oxygen", "Argon"}};
    std::vector<std::vector<double>> Z = {{},
                                          {0.3797, 0.3225, 0.278, 0.0014, 0.0184},
                                          {0.9677, 0.0323},
                                          {0.9403 / 1.00002, 0.0582 / 1.00002, 0.00127 / 1.00002, 0.00025 / 1.00002},
                                          {0.609067, 0.370414, 0.0205193}};
    int k;
    double T, p, rhocp;
    int cur = -1;
    std::unique_ptr<CoolProp::AbstractState> AS;
    while (std::scanf("%d %lf %lf %lf", &k, &T, &p, &rhocp) == 4) {
        if (k != cur) {
            setup_rp(F[k], true);
            AS.reset(CoolProp::AbstractState::factory("GERG2008", join(F[k])));
            cur = k;
        }
        auto& H = *dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        std::vector<double> z = Z[k];
        z.resize(20, 0.0);
        double pk = p / 1000, D = 0, Dl = 0, Dv = 0, xl[20] = {}, yv[20] = {}, q = 0, e = 0, hh = 0, ss = 0, cv = 0, cp = 0, w = 0;
        int ie = 0;
        char herr[256];
        RP_TPFLSH(&T, &pk, &z[0], &D, &Dl, &Dv, xl, yv, &q, &e, &hh, &ss, &cv, &cp, &w, &ie, herr, 255);
        const std::size_t N = F[k].size();
        std::vector<double> x(xl, xl + N), y(yv, yv + N), zz(Z[k]);
        H.set_mole_fractions(zz);
        H.update_DmolarT_direct(rhocp, T);
        const double gs = gRT(H);
        H.set_mole_fractions(x);
        H.update_DmolarT_direct(Dl * 1000, T);
        const double gl = gRT(H);
        std::vector<double> lfl(N);
        for (std::size_t i = 0; i < N; ++i)
            lfl[i] = std::log(x[i]) + CoolProp::MixtureDerivatives::ln_fugacity_coefficient(H, i, CoolProp::XN_INDEPENDENT);
        const double pL = H.p();
        H.set_mole_fractions(y);
        H.update_DmolarT_direct(Dv * 1000, T);
        const double gv = gRT(H);
        double dlnf = 0;
        for (std::size_t i = 0; i < N; ++i)
            dlnf = std::max(
              dlnf, std::abs(std::log(y[i]) + CoolProp::MixtureDerivatives::ln_fugacity_coefficient(H, i, CoolProp::XN_INDEPENDENT) - lfl[i]));
        const double pV = H.p();
        const double g2 = (1 - q) * gl + q * gv;
        std::printf("k=%d T=%8.3f p=%10.4g q=%.4f rhoL=%8.1f rhoV=%8.1f (p_L/p-1=%.1e p_V/p-1=%.1e) max|dlnf|=%.1e | g_split-g_single=%+.3e | x0 L/V "
                    "%.4f/%.4f | CoolProp rho=%.1f\n",
                    k, T, p, q, Dl * 1000, Dv * 1000, pL / p - 1, pV / p - 1, dlnf, g2 - gs, x[0], y[0], rhocp);
    }
}
