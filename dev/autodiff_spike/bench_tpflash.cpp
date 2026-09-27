// Experiment 13: CoolProp PT flash (update(PT_INPUTS), TPD stability + split) vs REFPROP TPFLSH, per state.
// (derived from bench_flash.cpp; GERG-2008 only, the same model on both sides)
//
// TPFLSH answers the question a user actually asks ("what is the fluid at (T, p, x)?") with a
// stability test, a phase split when two-phase, and caloric properties at the answer.  Ours
// returns the selected homogeneous root (no phase split).  So states are split by TPFLSH's verdict:
//   single-phase : like for like -- both return the density; agreement checked
//   two-phase    : REFPROP does a flash, we do not; timed separately, not a like-for-like comparison
// GERG rows: REFPROP FLAGS("GERG", 1) (same model; agreement to 1e-7).  Reference rows: REFPROP
// default vs CoolProp HEOS (model data differ slightly; agreement to 1e-3).
// Per-state times: each call timed individually, min over 3 passes.
//
//   ./bench_flash [nstates=5000] [refprop_dir=/Users/ianbell/REFPROP10/]
#include <dlfcn.h>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <random>
#include <string>
#include "AbstractState.h"
#include "gerg_cheb_solver.hpp"

using namespace gergcheb;

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
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"}};
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
}  // namespace

int main(int argc, char** argv) {
    const int NS = argc > 1 ? std::atoi(argv[1]) : 5000;
    const std::string rpdir = argc > 2 ? argv[2] : "/Users/ianbell/REFPROP10/";
    if (!load_rp(rpdir)) {
        std::printf("REFPROP not loaded\n");
        return 1;
    }
    const std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                           "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    const std::vector<double> AMA = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};
    const std::vector<std::string> AIRW = {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"};
    const std::vector<double> XAIRW = {0.7654, 0.2053, 0.0090, 0.0003, 0.02};
    struct Mix
    {
        std::string name, backend;
        std::vector<std::string> fluids;
        std::vector<double> x;
        double Tlo, Thi;
    };
    std::vector<Mix> mixes = {
      {"C1/C2 50/50", "GERG2008", {"Methane", "Ethane"}, {0.5, 0.5}, 150, 400},
      {"C1/C2/C3 50/30/20", "GERG2008", {"Methane", "Ethane", "Propane"}, {0.5, 0.3, 0.2}, 150, 450},
      {"Amarillo (10)", "GERG2008", NG10, AMA, 150, 400},
      {"C1/H2S 50/50", "GERG2008", {"Methane", "HydrogenSulfide"}, {0.5, 0.5}, 180, 450},
      {"humid air x_w=0.02", "GERG2008", AIRW, XAIRW, 250, 450},
    };
    std::printf("%d states per mixture: T uniform in [Tlo, Thi], p log-uniform in [10 kPa, 30 MPa]; single thread.\n", NS);
    std::printf("Times: mean (median) per state in us; states classified by TPFLSH's verdict.\n\n");
    for (const auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(mx.backend, join(mx.fluids)));
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        (void)H;
        const bool gerg = mx.backend == "GERG2008";
        setup_rp(mx.fluids, gerg);
        std::vector<double> z = mx.x;
        z.resize(20, 0.0);
        std::mt19937_64 g(42);
        std::vector<double> Ts(NS), ps(NS);
        for (int i = 0; i < NS; ++i) {
            Ts[i] = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
        }
        std::vector<double> t_cp(NS, 1e300), t_rp(NS, 1e300), rho_cp(NS, NAN), rho_rp(NS, NAN), q(NS, NAN);
        std::vector<int> ierr(NS, 0), cpfail(NS, 0);
        AS->set_mole_fractions(mx.x);
        volatile double sink = 0;
        for (int pass = 0; pass < 2; ++pass)
            for (int i = 0; i < NS; ++i) {
                double t0 = now_us();
                try {
                    AS->update(CoolProp::PT_INPUTS, ps[i], Ts[i]);
                    rho_cp[i] = AS->rhomolar();
                } catch (...) {
                    cpfail[i] = 1;
                }
                t_cp[i] = std::min(t_cp[i], now_us() - t0);
                double T = Ts[i], pk = ps[i] / 1000, D = 0, Dl = 0, Dv = 0, xl[20] = {}, yv[20] = {}, qq = 0, e = 0, hh = 0, ss = 0, cv = 0, cp = 0,
                       w = 0;
                int ie = 0;
                char herr[256];
                t0 = now_us();
                RP_TPFLSH(&T, &pk, &z[0], &D, &Dl, &Dv, xl, yv, &qq, &e, &hh, &ss, &cv, &cp, &w, &ie, herr, 255);
                t_rp[i] = std::min(t_rp[i], now_us() - t0);
                rho_rp[i] = D * 1000;
                q[i] = qq;
                ierr[i] = ie;
                sink = sink + D;
            }
        for (int two = 0; two < 2; ++two) {
            std::vector<double> c, r;
            long n = 0, fail = 0, agree = 0;
            for (int i = 0; i < NS; ++i) {
                if (ierr[i] > 0) continue;
                const bool is2 = q[i] >= 0 && q[i] <= 1;
                if (is2 != (two == 1)) continue;
                ++n;
                r.push_back(t_rp[i]);
                c.push_back(t_cp[i]);
                if (cpfail[i]) {
                    ++fail;
                    continue;
                }
                agree += std::abs(rho_cp[i] - rho_rp[i]) <= 1e-6 * rho_rp[i];
            }
            std::printf(
              "%-20s %-12s | n %5ld | CoolProp PT flash %8.1f (%7.1f) us, fails %5.1f %%, same bulk density %5.1f %% | TPFLSH %8.1f (%7.1f) us\n",
              mx.name.c_str(), two ? "two-phase" : "single-phase", n, mean(c), median(c), n ? 100.0 * fail / n : 0.0, n ? 100.0 * agree / n : 0.0,
              mean(r), median(r));
        }
        std::fflush(stdout);
    }
}
