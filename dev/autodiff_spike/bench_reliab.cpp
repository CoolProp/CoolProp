// Experiment 9: reliability of density solving -- ours vs CoolProp solver_rho_Tp vs REFPROP TPRHO.
//
// For every state we have *all* roots (ours).  Among the mechanically stable ones (dp/drho > 0)
// the stable root is the one of lowest molar Gibbs energy (evaluated by CoolProp for the same model
// at each root).  Each solver is then scored on:
//   fail        : no result (exception / ierr != 0 / out of range)
//   not a root  : result is not one of our roots (tolerance: 1e-8 same model; REFPROP 1e-4 reference
//                 models (its own R and HMX.BNC), 1e-6 GERG-2008 (R = 8.314472 vs CODATA))
//   stable      : result is the stable root (reported for multi-root states, where it can differ)
// Rows: reference-EOS mixtures (CoolProp HEOS vs REFPROP default) and GERG-2008 mixtures
// (CoolProp GERG2008 vs REFPROP with GERG08(1)).  Every failure / non-stable pick is written to a
// CSV for the failure maps (plot_reliab.py).
//
//   ./bench_reliab [nstates=5000] [csv=reliab.csv] [refprop_dir=/Users/ianbell/REFPROP10/]
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
using TPRHO_t = void (*)(double*, double*, double*, int*, int*, double*, int*, char*, long);
using FLAGS_t = void (*)(char*, int*, int*, int*, char*, long, long);
SETPATH_t RP_SETPATH = nullptr;
SETUP_t RP_SETUP = nullptr;
TPRHO_t RP_TPRHO = nullptr;
FLAGS_t RP_FLAGS = nullptr;
bool load_rp(const std::string& dir) {
    void* h = dlopen((dir + "librefprop.dylib").c_str(), RTLD_NOW);
    if (!h) return false;
    RP_SETPATH = (SETPATH_t)dlsym(h, "SETPATHdll");
    RP_SETUP = (SETUP_t)dlsym(h, "SETUPdll");
    RP_TPRHO = (TPRHO_t)dlsym(h, "TPRHOdll");
    RP_FLAGS = (FLAGS_t)dlsym(h, "FLAGSdll");
    if (!RP_SETPATH || !RP_SETUP || !RP_TPRHO || !RP_FLAGS) return false;
    char path[256] = {};
    std::strncpy(path, dir.c_str(), 255);
    RP_SETPATH(path, 255);
    return true;
}
const std::map<std::string, std::string> RPNAME = {
  {"Methane", "METHANE"},   {"Ethane", "ETHANE"},      {"Propane", "PROPANE"},   {"HydrogenSulfide", "H2S"},
  {"Nitrogen", "NITROGEN"}, {"Oxygen", "OXYGEN"},      {"Argon", "ARGON"},       {"CarbonDioxide", "CO2"},
  {"Water", "WATER"},       {"IsoButane", "ISOBUTAN"}, {"n-Butane", "BUTANE"},   {"Isopentane", "IPENTANE"},
  {"n-Pentane", "PENTANE"}, {"n-Hexane", "HEXANE"},    {"n-Heptane", "HEPTANE"}, {"n-Decane", "DECANE"}};
bool setup_rp(const std::vector<std::string>& fluids, bool gerg) {
    {  // REFPROP 10: FLAGS("GERG", 1) switches to GERG-2008 (pure and mixture) for subsequent SETUP calls
        char hf[256] = {}, herr[256] = {};
        std::strcpy(hf, "GERG");
        int j = gerg ? 1 : 0, k = 0, ierr = 0;
        RP_FLAGS(hf, &j, &k, &ierr, herr, 255, 255);
        if (ierr != 0) std::printf("   REFPROP FLAGS(GERG,%d): ierr %d %.120s\n", j, ierr, herr);
    }
    std::string f;
    for (std::size_t i = 0; i < fluids.size(); ++i)
        f += (i ? "|" : "") + RPNAME.at(fluids[i]) + ".FLD";
    std::vector<char> hf(10000, '\0'), herr(255, '\0');
    std::memcpy(hf.data(), f.data(), f.size());
    char hmx[255] = "HMX.BNC", hrf[3] = {'D', 'E', 'F'};
    int n = static_cast<int>(fluids.size()), ierr = 0;
    RP_SETUP(&n, hf.data(), hmx, hrf, &ierr, herr.data(), 10000, 255, 3, 255);
    if (ierr > 0) std::printf("   REFPROP SETUP error %d: %.200s\n", ierr, herr.data());
    if (gerg) {  // confirm the flag took
        char hf[256] = {}, herr[256] = {};
        std::strcpy(hf, "GERG");
        int j = -999, k = 0, ie = 0;
        RP_FLAGS(hf, &j, &k, &ie, herr, 255, 255);
        std::printf("   REFPROP GERG flag after SETUP: %d\n", k);
    }
    return ierr <= 0;
}
std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? "&" : "") + v[i];
    return s;
}
// Relative roundoff uncertainty of Z from a direct double-precision term sum (as gerg_validate.cpp):
// eps sum_k |term_k| (1 + |ln|weight_k||) / |Z|, plus the non-analytic terms.
double z_rel_uncert(const Solver& sv, const std::vector<double>& x, double T, double D) {
    const double tau = sv.Red->Tr(x) / T, lt = std::log(tau);
    double Z = 1, u = 0;
    for (const auto& tm : sv.terms) {
        const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * sv.F[tm.i][tm.j];
        const double w = X * tm.kappa(tau, lt), v = w * tm.chi(D);
        Z += v;
        u += std::abs(v) * (1 + std::abs(std::log(std::abs(w) + 1e-300)));
    }
    return 2.2e-16 * u / std::max(std::abs(Z), 1e-300);
}
struct Mix
{
    std::string name, backend;
    std::vector<std::string> fluids;
    std::vector<double> x;
    double Tlo, Thi;
};
struct Score
{
    long states = 0, multi = 0, fail = 0, notroot = 0, stable_multi = 0, returned_multi = 0;
    long either_stable = 0;
};
}  // namespace

int main(int argc, char** argv) {
    const int NS = argc > 1 ? std::atoi(argv[1]) : 5000;
    const std::string csvname = argc > 2 ? argv[2] : "reliab.csv";
    const std::string rpdir = argc > 3 ? argv[3] : "/Users/ianbell/REFPROP10/";
    const bool have_rp = load_rp(rpdir);
    if (!have_rp) std::printf("REFPROP not loaded\n");
    FILE* csv = std::fopen(csvname.c_str(), "w");
    std::fprintf(csv, "mixture,model,solver,T,p,nroots,outcome\n");
    const std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                           "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    const std::vector<double> AMARILLO = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};
    const std::vector<std::string> AIRW = {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"};
    std::vector<Mix> mixes;
    for (const std::string be : {"HEOS", "GERG2008"}) {
        mixes.push_back({"C1/C2 50/50", be, {"Methane", "Ethane"}, {0.5, 0.5}, 120, 400});
        mixes.push_back({"C1/C2/C3 50/30/20", be, {"Methane", "Ethane", "Propane"}, {0.5, 0.3, 0.2}, 120, 450});
        mixes.push_back({"Amarillo (10)", be, NG10, AMARILLO, 120, 400});
        mixes.push_back({"C1/H2S 50/50", be, {"Methane", "HydrogenSulfide"}, {0.5, 0.5}, 150, 450});
        mixes.push_back({"C1/nC7 70/30", be, {"Methane", "n-Heptane"}, {0.7, 0.3}, 200, 550});
        mixes.push_back({"CO2/H2O 50/50", be, {"CarbonDioxide", "Water"}, {0.5, 0.5}, 250, 650});
        mixes.push_back({"humid air x_w=0.02", be, AIRW, {0.7654, 0.2053, 0.0090, 0.0003, 0.02}, 230, 500});
    }
    std::printf("%d states per mixture: T uniform in [Tlo, Thi], p log-uniform in [10 kPa, 30 MPa]\n", NS);
    std::printf("'stable' = returned root is the lowest-Gibbs mechanically stable root, over states with >1 root\n\n");
    std::printf("%-20s %-8s %6s | %-6s | %-24s | %-24s | %-24s | %s\n", "mixture", "model", "multi", "ours", "CoolProp solver_rho_Tp",
                "REFPROP TPRHO kph=2", "REFPROP TPRHO kph=1", "REFPROP match");
    std::printf("%-20s %-8s %6s | %-6s | %7s %7s %8s | %7s %7s %8s | %7s %7s %8s |\n", "", "", "%", "fail", "fail%", "notroot", "stable%", "fail%",
                "notroot", "stable%", "fail%", "notroot", "stable%");
    for (const auto& mx : mixes) {
        std::unique_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(mx.backend, join(mx.fluids)));
        std::unique_ptr<CoolProp::AbstractState> AG(CoolProp::AbstractState::factory(mx.backend, join(mx.fluids)));  // for Gibbs
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        AS->set_mole_fractions(mx.x);
        AG->set_mole_fractions(mx.x);
        AG->specify_phase(CoolProp::iphase_gas);  // evaluate the EOS at given (rho, T) without phase logic
        Solver sv;
        sv.build(H, 4.0, 1.05 * H->Reducing->Tr(mx.x) / mx.Tlo, 1e-6);
        const bool gerg = mx.backend == "GERG2008";
        const bool rp_ok = have_rp && setup_rp(mx.fluids, gerg);
        const double rp_tol = gerg ? 1e-6 : 1e-4;
        std::vector<double> z = mx.x;
        z.resize(20, 0.0);
        std::mt19937_64 g(7);
        Score ours, cp, rp2, rp1;
        double rp_maxdev = 0;
        long n_artifact = 0, n_interior = 0, n_excursion = 0;
        long incomplete = 0, multi_states = 0, multi_stable_found = 0, multi_two_mech = 0, rp_either = 0;
        Solver::State S;
        for (int i = 0; i < NS; ++i) {
            const double T = mx.Tlo + (mx.Thi - mx.Tlo) * std::uniform_real_distribution<double>(0, 1)(g);
            const double p = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
            sv.assemble(T, mx.x, S);
            Solver::Root r[MAXROOTS];
            const int n = sv.roots(S, p, r);
            {  // roots beyond delta_max: our set is incomplete there, so do not score this state
                double Ge, dGe, sce;
                sv.true_G(S, sv.delta_max, p * S.t_scale, Ge, dGe, sce);
                if (Ge < 0) {
                    ++incomplete;
                    continue;
                }
            }
            // stable root among mechanically stable roots
            // Stable root: among the OUTER mechanically stable roots only (lowest- and highest-density
            // roots with dp/drho > 0), the one of lower g.  Interior roots with dp/drho > 0 exist in
            // multiparameter EOS inside the two-phase dome at low T (5-root states) and can have absurd
            // g (alphar ~ -6e4 seen for C1/C2 at 132 K); they are EOS artifacts, counted separately.
            std::vector<bool> artifact(n, false);
            std::vector<double> gk(n, NAN);
            std::vector<bool> mech(n, false);
            for (int k = 0; k < n; ++k) {
                artifact[k] = z_rel_uncert(sv, mx.x, T, r[k].rho / S.rhor) > 1e-6;  // Z not resolvable in double here
                n_artifact += artifact[k];
                if (artifact[k]) continue;
                try {
                    AG->update(CoolProp::DmolarT_INPUTS, r[k].rho, T);
                    const double ar = AG->alphar();
                    if (std::abs(ar) > 100) {  // EOS excursion inside the dome (seen: alphar -1e3..-6e4 near delta = 1 at low T)
                        ++n_excursion;
                        continue;
                    }
                    mech[k] = AG->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT) > 0;
                    // g/RT up to a (T, x) constant common to all roots: ln(rho) + alphar + Z (needs no ideal-gas part)
                    gk[k] = std::log(r[k].rho) + ar + AG->p() / (r[k].rho * AG->gas_constant() * T);
                } catch (...) {
                }
            }
            int ilo = -1, ihi = -1;
            for (int k = 0; k < n; ++k)
                if (mech[k]) {
                    if (ilo < 0) ilo = k;
                    ihi = k;
                }
            for (int k = 0; k < n; ++k)
                n_interior += mech[k] && k != ilo && k != ihi;
            int istable = ilo;
            if (ilo >= 0 && ihi != ilo && gk[ihi] < gk[ilo]) istable = ihi;
            int nres = 0;
            for (int k = 0; k < n; ++k)
                nres += !artifact[k];
            const bool multi = nres > 1;
            int nmech = 0;
            for (int k = 0; k < n; ++k)
                nmech += mech[k];
            if (multi) {
                ++multi_states;
                multi_stable_found += istable >= 0;
                multi_two_mech += nmech >= 2;
            }
            int rp_hit = 0;
            static int ndbg = 0;
            const bool dbg = std::getenv("RELIAB_DEBUG") && multi && ndbg < 5 && mx.backend == "GERG2008" && mx.name == "C1/C2 50/50" && T > 185
                             && T < 215 && p > 3e6;
            if (dbg) {
                ++ndbg;
                std::printf("  DEBUG T=%.3f p=%.6g: %d roots\n", T, p, n);
                for (int k = 0; k < n; ++k) {
                    AG->update(CoolProp::DmolarT_INPUTS, r[k].rho, T);
                    std::printf("     root %d rho=%.8g delta=%.4f dpdrho=%.4g alphar=%.4g g=%.10f%s\n", k, r[k].rho, r[k].rho / S.rhor,
                                AG->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT), AG->alphar(),
                                std::log(r[k].rho) + AG->alphar() + AG->p() / (r[k].rho * AG->gas_constant() * T), k == istable ? "  <- stable" : "");
                }
                for (int kph : {2, 1}) {
                    double TT = T, pk = p / 1000, D = 0;
                    int kk = kph, kg = 0, ierr = 0;
                    char herr[256];
                    if (rp_ok) RP_TPRHO(&TT, &pk, &z[0], &kk, &kg, &D, &ierr, herr, 255);
                    std::printf("     REFPROP kph=%d -> rho=%.8g ierr=%d\n", kph, D * 1000, ierr);
                }
                try {
                    std::printf("     CoolProp -> rho=%.8g\n", H->solver_rho_Tp(T, p));
                } catch (...) {
                    std::printf("     CoolProp -> throw\n");
                }
            }
            auto which = [&](double rho, double tol) {
                int best = -1;
                double bd = 1e300;
                for (int k = 0; k < n; ++k) {
                    const double d = std::abs(r[k].rho - rho) / rho;
                    if (d < bd) {
                        bd = d;
                        best = k;
                    }
                }
                return std::make_pair(bd <= tol ? best : -1, bd);
            };
            auto score = [&](Score& s, const char* solver, bool ok, double rho, double tol) {
                ++s.states;
                s.multi += multi;
                if (!ok) {
                    ++s.fail;
                    std::fprintf(csv, "%s,%s,%s,%.6f,%.6g,%d,fail\n", mx.name.c_str(), mx.backend.c_str(), solver, T, p, n);
                    return -1.0;
                }
                const auto [k, dev] = which(rho, tol);
                if (k < 0) {
                    ++s.notroot;
                    std::fprintf(csv, "%s,%s,%s,%.6f,%.6g,%d,notroot\n", mx.name.c_str(), mx.backend.c_str(), solver, T, p, n);
                    return dev;
                }
                if (multi) {
                    ++s.returned_multi;
                    if (k == istable && std::string(solver).rfind("REFPROP", 0) == 0) rp_hit = 1;
                    if (k == istable)
                        ++s.stable_multi;
                    else
                        std::fprintf(csv, "%s,%s,%s,%.6f,%.6g,%d,nonstable\n", mx.name.c_str(), mx.backend.c_str(), solver, T, p, n);
                }
                return dev;
            };
            ++ours.states;
            ours.multi += multi;
            // CoolProp
            double rho_cp = NAN;
            bool ok_cp = true;
            try {
                rho_cp = H->solver_rho_Tp(T, p);
                ok_cp = std::isfinite(rho_cp) && rho_cp > 0;
            } catch (...) {
                ok_cp = false;
            }
            score(cp, "CoolProp", ok_cp, rho_cp, 1e-8);
            // REFPROP, vapor-like and liquid-like
            if (rp_ok)
                for (int kph : {2, 1}) {
                    double TT = T, pk = p / 1000, D = 0;
                    int kk = kph, kg = 0, ierr = 0;
                    char herr[256];
                    RP_TPRHO(&TT, &pk, &z[0], &kk, &kg, &D, &ierr, herr, 255);
                    const double dev = score(kph == 2 ? rp2 : rp1, kph == 2 ? "REFPROP_kph2" : "REFPROP_kph1", ierr == 0 && D > 0, D * 1000, rp_tol);
                    if (dev >= 0) rp_maxdev = std::max(rp_maxdev, dev);
                }
            rp_either += multi && rp_hit;
        }
        auto pct = [](long a, long b) { return b ? 100.0 * a / b : 0.0; };
        std::printf("%-20s %-8s %6.1f | %6d | %7.2f %7ld %8.1f | %7.2f %7ld %8.1f | %7.2f %7ld %8.1f | max dev %.1e\n", mx.name.c_str(),
                    mx.backend.c_str(), pct(ours.multi, ours.states), 0, pct(cp.fail, cp.states), cp.notroot, pct(cp.stable_multi, cp.returned_multi),
                    pct(rp2.fail, rp2.states), rp2.notroot, pct(rp2.stable_multi, rp2.returned_multi), pct(rp1.fail, rp1.states), rp1.notroot,
                    pct(rp1.stable_multi, rp1.returned_multi), rp_maxdev);
        if (incomplete) std::printf("   (%ld states skipped: a root lies beyond delta_max = 4)\n", incomplete);
        if (n_artifact || n_interior || n_excursion)
            std::printf(
              "   roots excluded from the stable choice: %ld numerical (Z unresolvable in double), %ld EOS excursions (|alphar| > 100), %ld interior "
              "mechanically-stable branches\n",
              n_artifact, n_excursion, n_interior);
        std::printf("   multi-root states %ld: stable root identified in %ld, >=2 mechanically stable roots in %ld; REFPROP (either kph) returned "
                    "the stable root in %.1f %%\n",
                    multi_states, multi_stable_found, multi_two_mech, multi_states ? 100.0 * rp_either / multi_states : 0.0);
        std::fflush(stdout);
    }
    std::fclose(csv);
}
