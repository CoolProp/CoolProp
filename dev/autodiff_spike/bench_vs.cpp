// Experiment 8: density-solve speed vs CoolProp and REFPROP, same states, same compositions.
//
//   ours      all roots in (0, delta_max], polished on the true equation
//   CoolProp  HelmholtzEOSMixtureBackend::solver_rho_Tp(T, p)  (one root, its own guess)
//   REFPROP   TPRHOdll(kph = 2) + TPRHOdll(kph = 1)  (vapor-like and liquid-like root, kguess = 0)
//
// Reference-EOS mixtures use CoolProp's HEOS backend (the same pure-fluid EOS and mixture
// parameters REFPROP uses, up to data-set differences); GERG-2008 rows use CoolProp's GERG2008
// backend (no REFPROP counterpart run here).  Single-threaded; min over repeats.
//
//   ./bench_vs [nstates=2000] [refprop_dir=/Users/ianbell/REFPROP10/]
#include <chrono>
#include <cstdlib>
#include <random>
#include <string>
#include <dlfcn.h>
#include <cstring>
#include <algorithm>
#include <map>
#include "AbstractState.h"
#include "Configuration.h"
#include "gerg_cheb_solver.hpp"

using namespace gergcheb;

namespace {
// REFPROP 10, loaded directly (Fortran calling convention: everything by reference, hidden string lengths)
using SETPATH_t = void (*)(char*, long);
using SETUP_t = void (*)(int*, char*, char*, char*, int*, char*, long, long, long, long);
using TPRHO_t = void (*)(double*, double*, double*, int*, int*, double*, int*, char*, long);
SETPATH_t RP_SETPATH = nullptr;
SETUP_t RP_SETUP = nullptr;
TPRHO_t RP_TPRHO = nullptr;
bool load_rp(const std::string& dir) {
    void* h = dlopen((dir + "librefprop.dylib").c_str(), RTLD_NOW);
    if (!h) return false;
    RP_SETPATH = (SETPATH_t)dlsym(h, "SETPATHdll");
    RP_SETUP = (SETUP_t)dlsym(h, "SETUPdll");
    RP_TPRHO = (TPRHO_t)dlsym(h, "TPRHOdll");
    if (!RP_SETPATH || !RP_SETUP || !RP_TPRHO) return false;
    char path[256] = {};
    std::strncpy(path, dir.c_str(), 255);
    RP_SETPATH(path, 255);
    return true;
}
const std::map<std::string, std::string> RPNAME = {
  {"Methane", "METHANE"}, {"Ethane", "ETHANE"},       {"Propane", "PROPANE"},   {"HydrogenSulfide", "H2S"}, {"Nitrogen", "NITROGEN"},
  {"Oxygen", "OXYGEN"},   {"Argon", "ARGON"},         {"CarbonDioxide", "CO2"}, {"Water", "WATER"},         {"IsoButane", "ISOBUTAN"},
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"}};
bool setup_rp(const std::vector<std::string>& fluids, std::string& msg) {
    std::string f;
    for (std::size_t i = 0; i < fluids.size(); ++i)
        f += (i ? "|" : "") + RPNAME.at(fluids[i]) + ".FLD";
    std::vector<char> hf(10000, '\0'), herr(255, '\0');
    std::memcpy(hf.data(), f.data(), f.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = static_cast<int>(fluids.size()), ierr = 0;
    RP_SETUP(&n, hf.data(), hmx, hrf, &ierr, herr.data(), 10000, 255, 3, 255);
    msg.assign(herr.data(), strnlen(herr.data(), 255));
    return ierr <= 0;
}
std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? "&" : "") + v[i];
    return s;
}
struct Mix
{
    std::string name, backend;
    std::vector<std::string> fluids;
    std::vector<double> x;
    double Tlo, Thi;
    bool refprop;
};
template <class F>
double time_us(F&& f, int n, int reps = 3) {
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < n; ++i)
            f(i);
        best = std::min(best, std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / n);
    }
    return best;
}
}  // namespace

int main(int argc, char** argv) {
    const int NS = argc > 1 ? std::atoi(argv[1]) : 2000;
    const std::string rpdir = argc > 2 ? argv[2] : "/Users/ianbell/REFPROP10/";
    if (!load_rp(rpdir)) std::printf("REFPROP not loaded from %s\n", rpdir.c_str());
    const std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                           "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    const std::vector<double> AMARILLO = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};
    const std::vector<Mix> mixes = {
      {"C1/C2 50/50", "HEOS", {"Methane", "Ethane"}, {0.5, 0.5}, 150, 400, true},
      {"C1/C2/C3 50/30/20", "HEOS", {"Methane", "Ethane", "Propane"}, {0.5, 0.3, 0.2}, 150, 450, true},
      {"Amarillo (10)", "HEOS", NG10, AMARILLO, 150, 400, true},
      {"C1/H2S 50/50", "HEOS", {"Methane", "HydrogenSulfide"}, {0.5, 0.5}, 180, 450, true},
      {"humid air, x_w=0.02",
       "HEOS",
       {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"},
       {0.7654, 0.2053, 0.0090, 0.0003, 0.02},
       250,
       450,
       true},
      {"C1/C2 50/50", "GERG2008", {"Methane", "Ethane"}, {0.5, 0.5}, 150, 400, false},
      {"Amarillo (10)", "GERG2008", NG10, AMARILLO, 150, 400, false},
      {"humid air, x_w=0.02",
       "GERG2008",
       {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"},
       {0.7654, 0.2053, 0.0090, 0.0003, 0.02},
       250,
       450,
       false},
    };
    std::printf("%d states per mixture: T uniform in [Tlo, Thi], p log-uniform in [10 kPa, 30 MPa]; single thread; times are per state\n\n", NS);
    std::printf("%-22s %-8s %6s %7s | %9s %6s | %9s %6s %8s | %10s %7s %8s %9s\n", "mixture", "model", "roots", "pieces", "ours_us", "", "CP_us",
                "fail", "CP_in", "RP_us(2x)", "fail", "RP_in", "RP_notin");
    for (const auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(mx.backend, join(mx.fluids)));
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        AS->set_mole_fractions(mx.x);
        Solver sv;
        sv.build(H, 4.0, 1.05 * H->Reducing->Tr(mx.x) / mx.Tlo, 1e-6);
        std::mt19937_64 g(42);
        std::vector<double> Ts(NS), ps(NS);
        for (int i = 0; i < NS; ++i) {
            Ts[i] = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
        }
        // ours: all roots
        std::vector<std::vector<double>> roots(NS);
        Solver::State S;
        long nroots = 0;
        for (int i = 0; i < NS; ++i) {
            sv.assemble(Ts[i], mx.x, S);
            Solver::Root r[MAXROOTS];
            const int n = sv.roots(S, ps[i], r);
            for (int k = 0; k < n; ++k)
                roots[i].push_back(r[k].rho);
            nroots += n;
        }
        volatile double sink = 0;
        const double t_ours = time_us(
          [&](int i) {
              sv.assemble(Ts[i], mx.x, S);
              Solver::Root r[MAXROOTS];
              sink = sink + sv.roots(S, ps[i], r);
          },
          NS);
        auto in_set = [&](int i, double rho, double tol) {
            for (double v : roots[i])
                if (std::abs(v - rho) <= tol * rho) return true;
            return false;
        };
        // CoolProp solver_rho_Tp
        int cp_fail = 0, cp_in = 0;
        for (int i = 0; i < NS; ++i) {
            try {
                const double rho = H->solver_rho_Tp(Ts[i], ps[i]);
                cp_in += in_set(i, rho, 1e-8);
            } catch (...) {
                ++cp_fail;
            }
        }
        const double t_cp = time_us(
          [&](int i) {
              try {
                  sink = sink + H->solver_rho_Tp(Ts[i], ps[i]);
              } catch (...) {
              }
          },
          NS);
        // REFPROP TPRHO, vapor-like and liquid-like
        double t_rp = NAN;
        int rp_fail = 0, rp_in = 0, rp_notin = 0;
        std::vector<double> rp_dev, rp_all;
        std::string rpmsg;
        if (mx.refprop && RP_TPRHO && setup_rp(mx.fluids, rpmsg)) {
            std::vector<double> z = mx.x;
            z.resize(20, 0.0);
            char herr[256];
            auto tprho = [&](double T, double p, int kph, double& D) {
                double pk = p / 1000;
                int kg = 0, ierr = 0;
                RP_TPRHO(&T, &pk, &z[0], &kph, &kg, &D, &ierr, herr, 255);
                D *= 1000;  // mol/L -> mol/m^3
                return ierr;
            };
            for (int i = 0; i < NS; ++i)
                for (int kph : {2, 1}) {
                    double D = 0;
                    if (tprho(Ts[i], ps[i], kph, D) != 0) {
                        ++rp_fail;
                        continue;
                    }
                    {
                        double best = 1e300;
                        for (double v : roots[i])
                            best = std::min(best, std::abs(v - D) / D);
                        rp_all.push_back(best);
                    }
                    if (in_set(i, D, 1e-4))  // REFPROP's own R and HMX.BNC parameters: ~1e-6..1e-3 model differences
                        ++rp_in;
                    else {
                        ++rp_notin;
                        double best = 1e300;
                        for (double v : roots[i])
                            best = std::min(best, std::abs(v - D) / D);
                        rp_dev.push_back(best);
                        if (std::getenv("RP_DEBUG") && rp_notin <= 4) {
                            std::printf("   RP kph=%d T=%.3f p=%.6g: RP rho=%.10g, ours:", kph, Ts[i], ps[i], D);
                            for (double v : roots[i])
                                std::printf(" %.10g", v);
                            std::printf("\n");
                        }
                    }
                }
            t_rp = time_us(
              [&](int i) {
                  double D1 = 0, D2 = 0;
                  tprho(Ts[i], ps[i], 2, D1);
                  tprho(Ts[i], ps[i], 1, D2);
                  sink = sink + D1 + D2;
              },
              NS);
        }
        if (!rp_all.empty()) {
            std::sort(rp_all.begin(), rp_all.end());
            std::printf("   REFPROP roots vs nearest of ours (model differences): median %.1e, 99th pct %.1e, max %.1e\n", rp_all[rp_all.size() / 2],
                        rp_all[rp_all.size() * 99 / 100], rp_all.back());
        }
        if (!rp_dev.empty()) {
            std::sort(rp_dev.begin(), rp_dev.end());
            std::printf("   REFPROP roots not within 1e-4 of ours: median rel. deviation from nearest root %.1e, 90th pct %.1e, max %.1e\n",
                        rp_dev[rp_dev.size() / 2], rp_dev[rp_dev.size() * 9 / 10], rp_dev.back());
        }
        std::printf("%-22s %-8s %6.2f %7d | %9.2f %6s | %9.2f %6d %8d | %10.2f %7d %8d %9d\n", mx.name.c_str(), mx.backend.c_str(),
                    double(nroots) / NS, sv.P, t_ours, "", t_cp, cp_fail, cp_in, t_rp, rp_fail, rp_in, rp_notin);
        std::fflush(stdout);
    }
}
