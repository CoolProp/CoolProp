// Experiment 11: root-selection policies, scored against REFPROP's GERG-2008 flash.
//
// Given ALL density roots at (T, p, x), which one is "the" density?  Candidates:
//   P0  min Gibbs among all mechanically stable roots (dp/drho > 0)
//   P1  min Gibbs among the outermost mechanically stable roots (lowest / highest density)
//   P2  P1, excluding |alphar| > 100 excursions
//   P3  spinodal branches: vapor branch = delta below the first local max of F = delta Z(delta);
//       liquid branch = delta above the last local min of F; roots in between are excluded
//       (EOS wiggles inside the dome, the deep alphar well near delta = 1 at low T); min Gibbs
//       between the branch roots.  Uses the shape of the isotherm, no thresholds.
// Arbiter: REFPROP TPFLSH with FLAGS("GERG", 1) -- same model, full phase-stability analysis.
// Where it reports a single phase (q outside [0, 1]) its density is the stable root; a policy
// is right if its pick equals it (1e-7).  Two-phase states are reported separately.
//
//   ./bench_policy [nstates=5000] [refprop_dir=/Users/ianbell/REFPROP10/]
#include <dlfcn.h>
#include <algorithm>
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
  {"n-Butane", "BUTANE"}, {"Isopentane", "IPENTANE"}, {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"},     {"n-Heptane", "HEPTANE"}};
bool setup_rp_gerg(const std::vector<std::string>& fluids) {
    char hf0[256] = {}, herr0[256] = {};
    std::strcpy(hf0, "GERG");
    int j = 1, k = 0, ie = 0;
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
// extrema of F = delta Z(delta) from the tables (roots of dG/d delta), refined on the true F'
std::vector<std::pair<double, int>> extrema(const Solver& sv, const Solver::State& S) {
    std::vector<std::pair<double, int>> out;  // (delta, +1 local max / -1 local min)
    for (int pc = 0; pc < sv.P; ++pc) {
        const VecG d = chebder(S.G[pc], NG);
        double sa = 0;
        for (double v : d)
            sa += std::abs(v);
        RootOut ro[MAXROOTS];
        const int nr = bern_roots(d, 1e-13 * sa, ro);
        for (int k = 0; k < nr; ++k) {
            double D = sv.edges[pc] + (sv.edges[pc + 1] - sv.edges[pc]) * (ro[k].u + 1) / 2;
            double G, dl, dr, sc;
            const double h = 1e-7 * D;
            sv.true_G(S, D - h, 0, G, dl, sc);
            sv.true_G(S, D + h, 0, G, dr, sc);
            if ((dl > 0) && (dr < 0)) out.push_back({D, +1});
            if ((dl < 0) && (dr > 0)) out.push_back({D, -1});
        }
    }
    std::sort(out.begin(), out.end());
    return out;
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
    struct Mix
    {
        std::string name;
        std::vector<std::string> fluids;
        std::vector<double> x;
        double Tlo, Thi;
    };
    const std::vector<Mix> mixes = {
      {"C1/C2 50/50", {"Methane", "Ethane"}, {0.5, 0.5}, 120, 400},
      {"C1/C2/C3 50/30/20", {"Methane", "Ethane", "Propane"}, {0.5, 0.3, 0.2}, 120, 450},
      {"Amarillo (10)", NG10, {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393}, 120, 400},
      {"C1/H2S 50/50", {"Methane", "HydrogenSulfide"}, {0.5, 0.5}, 150, 450},
      {"C1/nC7 70/30", {"Methane", "n-Heptane"}, {0.7, 0.3}, 200, 550},
      {"CO2/H2O 50/50", {"CarbonDioxide", "Water"}, {0.5, 0.5}, 250, 650},
    };
    const char* PN[5] = {"P0 min-G all", "P1 min-G outer", "P2 P1+|ar|<100", "P3 spinodal branches", "Solver::select"};
    std::printf("%d states per mixture (GERG-2008).  Arbiter: REFPROP TPFLSH (GERG mode), single-phase states only.\n", NS);
    std::printf("Entries: %% of single-phase multi-root states where the policy picks REFPROP's density (and count of picks that are wrong)\n\n");
    std::printf("%-20s %8s %8s | %-20s | %-20s | %-20s | %-20s | %s\n", "mixture", "1-phase", "2-phase", PN[0], PN[1], PN[2], PN[3],
                "RP fail");  // (5th column: Solver::select)
    long tot_right[5] = {}, tot_n = 0;
    for (const auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("GERG2008", join(mx.fluids)));
        std::unique_ptr<CoolProp::AbstractState> AG(CoolProp::AbstractState::factory("GERG2008", join(mx.fluids)));
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        AG->set_mole_fractions(mx.x);
        AG->specify_phase(CoolProp::iphase_gas);
        Solver sv;
        sv.build(H, 4.0, 1.05 * H->Reducing->Tr(mx.x) / mx.Tlo, 1e-6);
        setup_rp_gerg(mx.fluids);
        std::vector<double> z = mx.x;
        z.resize(20, 0.0);
        std::mt19937_64 g(11);
        long n1 = 0, n2 = 0, rpfail = 0, right[5] = {}, wrong[5] = {};
        Solver::State S;
        for (int i = 0; i < NS; ++i) {
            double T = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            const double p = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
            sv.assemble(T, mx.x, S);
            Solver::Root r[MAXROOTS];
            const int n = sv.roots(S, p, r);
            if (n < 2) continue;  // policies only matter with >1 root
            // REFPROP flash
            double pk = p / 1000, D = 0, Dl = 0, Dv = 0, xl[20] = {}, yv[20] = {}, q = 0, e = 0, hh = 0, ss = 0, cv = 0, cp = 0, w = 0;
            int ierr = 0;
            char herr[256];
            RP_TPFLSH(&T, &pk, &z[0], &D, &Dl, &Dv, xl, yv, &q, &e, &hh, &ss, &cv, &cp, &w, &ierr, herr, 255);
            if (ierr > 0) {
                ++rpfail;
                continue;
            }
            if (q >= 0 && q <= 1) {
                ++n2;
                continue;
            }
            ++n1;
            const double rho_rp = D * 1000;
            // per-root data
            std::vector<double> gk(n, NAN), ar(n, NAN);
            std::vector<bool> mech(n, false);
            for (int k = 0; k < n; ++k) {
                try {
                    AG->update(CoolProp::DmolarT_INPUTS, r[k].rho, T);
                    ar[k] = AG->alphar();
                    mech[k] = AG->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT) > 0;
                    gk[k] = std::log(r[k].rho) + ar[k] + AG->p() / (r[k].rho * AG->gas_constant() * T);
                } catch (...) {
                }
            }
            auto argmin_g = [&](const std::vector<int>& cand) {
                int b = -1;
                for (int k : cand)
                    if (b < 0 || gk[k] < gk[b]) b = k;
                return b;
            };
            int pick[5] = {-1, -1, -1, -1, -1};
            pick[4] = sv.select(S, p, r, n);
            {  // P0
                std::vector<int> c;
                for (int k = 0; k < n; ++k)
                    if (mech[k]) c.push_back(k);
                pick[0] = argmin_g(c);
            }
            for (int pol : {1, 2}) {  // P1, P2
                int lo = -1, hi = -1;
                for (int k = 0; k < n; ++k)
                    if (mech[k] && (pol == 1 || std::abs(ar[k]) <= 100)) {
                        if (lo < 0) lo = k;
                        hi = k;
                    }
                std::vector<int> c;
                if (lo >= 0) c.push_back(lo);
                if (hi >= 0 && hi != lo) c.push_back(hi);
                pick[pol] = argmin_g(c);
            }
            {  // P3
                const auto ex = extrema(sv, S);
                double dmax1 = 1e300, dminL = -1;  // first local max, last local min
                for (auto& e2 : ex)
                    if (e2.second > 0) {
                        dmax1 = e2.first;
                        break;
                    }
                for (auto it = ex.rbegin(); it != ex.rend(); ++it)
                    if (it->second < 0) {
                        dminL = it->first;
                        break;
                    }
                std::vector<int> c;
                int vap = -1, liq = -1;
                for (int k = 0; k < n; ++k) {
                    const double d = r[k].rho / S.rhor;
                    if (d < dmax1 && mech[k]) vap = k;             // on the vapor branch (take the last: at most one)
                    if (d > dminL && mech[k] && liq < 0) liq = k;  // on the liquid branch (first: at most one)
                }
                if (vap >= 0) c.push_back(vap);
                if (liq >= 0 && liq != vap) c.push_back(liq);
                pick[3] = c.empty() ? pick[2] : argmin_g(c);
            }
            for (int pol = 0; pol < 5; ++pol) {
                const bool ok = pick[pol] >= 0 && std::abs(r[pick[pol]].rho - rho_rp) <= 1e-7 * rho_rp;
                (ok ? right : wrong)[pol]++;
                if (!ok && pol == 3 && std::getenv("POLICY_DEBUG")) {
                    std::printf("  P3 wrong: %s T=%.3f p=%.6g q=%g RP rho=%.8g; roots:", mx.name.c_str(), T, p, q, rho_rp);
                    for (int k = 0; k < n; ++k)
                        std::printf(" %.6g%s", r[k].rho, k == pick[3] ? "*" : "");
                    std::printf("\n");
                }
            }
        }
        std::printf("%-20s %8ld %8ld |", mx.name.c_str(), n1, n2);
        for (int pol = 0; pol < 5; ++pol)
            std::printf(" %8.2f %% (%5ld)    |", n1 ? 100.0 * right[pol] / n1 : 0.0, wrong[pol]);
        std::printf(" %ld\n", rpfail);
        for (int pol = 0; pol < 5; ++pol)
            tot_right[pol] += right[pol];
        tot_n += n1;
    }
    std::printf("\nALL: %ld single-phase multi-root states;", tot_n);
    for (int pol = 0; pol < 5; ++pol)
        std::printf("  %s %.3f %%", PN[pol], tot_n ? 100.0 * tot_right[pol] / tot_n : 0.0);
    std::printf("\n");
}
