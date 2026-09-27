// Experiment 7: large-scale validation of the GERG-2008 Chebyshev all-roots density solver.
//
//   ./gerg_validate calls=1e8 threads=8 seed=1 tol=1e-6 scan_every=20000 cp_every=200000 spin=20000 log=mismatch.log
//
// Suites (component sets built once; composition, T, p random per call):
//   binaries   all 210 GERG-2008 pairs
//   multi      60 random 3..10-component subsets + the full 21-component set
//   natgas     6 predefined natural gases + air (exact or log-normally perturbed)
//   asym       12 strongly asymmetric / type-III candidate pairs (denser sampling)
//   humidair   N2/O2/Ar/CO2/H2O, water mole fraction log-uniform 1e-5..0.3
// T log-uniform in [max(60 K, 0.3 max Tc), 700 K]; p log-uniform in [100 Pa, 100 MPa].
//
// Per call:     root-count parity (G(0) = -t < 0, so #roots is odd iff G_true(delta_max) > 0);
//               residual of every polished root on the true equation; polish movement;
//               "uncertain" flags from the subdivision.
// Subsampled:   dense scan of the true equation (count + positions, mismatches adjudicated);
//               CoolProp solver_rho_Tp root must be in our set; CoolProp's own p at our root.
// spin=N:       spinodal stress test, p = p_spinodal (1 +- eps), eps = 1e-2 .. 1e-10.
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include "AbstractState.h"
#include "gerg_cheb_solver.hpp"

using namespace gergcheb;

namespace {
const std::vector<std::string> GERG21 = {"Methane",   "Nitrogen",   "CarbonDioxide",  "Ethane",    "Propane",         "n-Butane", "IsoButane",
                                         "n-Pentane", "Isopentane", "n-Hexane",       "n-Heptane", "n-Octane",        "n-Nonane", "n-Decane",
                                         "Hydrogen",  "Oxygen",     "CarbonMonoxide", "Water",     "HydrogenSulfide", "Helium",   "Argon"};

struct CompSet
{
    std::string suite, name, backend = "GERG2008";
    std::vector<std::string> fluids;
    std::vector<double> x0;  // nominal composition (natgas, humidair)
    std::unique_ptr<CoolProp::AbstractState> AS, AScp;
    Solver sv;
    double Tmin = 0, Tmax = 700;
    std::mutex cp_mtx;
};

std::string join(const std::vector<std::string>& v) {
    std::string s;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? "&" : "") + v[i];
    return s;
}

struct Stats
{
    long calls = 0, roots = 0, multi_root_calls = 0, parity_fail = 0, beyond = 0, uncertain = 0, uncertified = 0, bracket_bad = 0, polish_moved = 0;
    double max_resid = 0, max_move = 0, ns = 0;
    long scans = 0, scan_ok = 0, scan_missed_real = 0, solver_missed = 0, solver_spurious = 0;
    double scan_maxrel = 0;
    long cp_tries = 0, cp_converged = 0, cp_in_set = 0, cp_not_in_set = 0, cp_pchecks = 0, cp_illcond = 0;
    double cp_pmax = 0;
    long nodes = 0, pol_its = 0, pol_calls = 0, na_active = 0;
    void add(const Stats& o) {
        calls += o.calls;
        roots += o.roots;
        multi_root_calls += o.multi_root_calls;
        parity_fail += o.parity_fail;
        beyond += o.beyond;
        uncertain += o.uncertain;
        uncertified += o.uncertified;
        bracket_bad += o.bracket_bad;
        polish_moved += o.polish_moved;
        max_resid = std::max(max_resid, o.max_resid);
        max_move = std::max(max_move, o.max_move);
        ns += o.ns;
        scans += o.scans;
        scan_ok += o.scan_ok;
        scan_missed_real += o.scan_missed_real;
        solver_missed += o.solver_missed;
        solver_spurious += o.solver_spurious;
        scan_maxrel = std::max(scan_maxrel, o.scan_maxrel);
        cp_tries += o.cp_tries;
        cp_converged += o.cp_converged;
        cp_in_set += o.cp_in_set;
        cp_not_in_set += o.cp_not_in_set;
        cp_pchecks += o.cp_pchecks;
        cp_illcond += o.cp_illcond;
        cp_pmax = std::max(cp_pmax, o.cp_pmax);
        nodes += o.nodes;
        na_active += o.na_active;
        pol_its += o.pol_its;
        pol_calls += o.pol_calls;
    }
};

std::mutex log_mtx;
FILE* g_logfile = nullptr;
std::atomic<long> log_lines{0};
template <class... A>
void logmsg(const char* fmt, A... a) {
    if (!g_logfile || log_lines.fetch_add(1) > 5000) return;
    std::lock_guard<std::mutex> lk(log_mtx);
    std::fprintf(g_logfile, fmt, a...);
}

// ---------------------------------------------------------------- sampling
struct Rng
{
    std::mt19937_64 g;
    explicit Rng(uint64_t s) : g(s) {}
    double u() {
        return std::uniform_real_distribution<double>(0, 1)(g);
    }
    double logu(double a, double b) {
        return std::exp(std::log(a) + u() * (std::log(b) - std::log(a)));
    }
    double gamma(double a) {
        return std::gamma_distribution<double>(a, 1.0)(g);
    }
    double normal() {
        return std::normal_distribution<double>(0, 1)(g);
    }
};
void sample_x(CompSet& cs, Rng& r, std::vector<double>& x) {
    const std::size_t N = cs.fluids.size();
    x.assign(N, 0.0);
    if (cs.suite == "natgas" || cs.suite == "natgas_ref") {
        x = cs.x0;
        if (r.u() < 0.5)
            for (auto& v : x)
                v *= std::exp(0.5 * r.normal());
    } else if (cs.suite.rfind("humidair", 0) == 0 || cs.suite == "wetair_ref") {
        const double xw = cs.suite == "wetair_ref" ? r.logu(0.3, 0.99) : r.logu(1e-5, 0.3);
        for (std::size_t i = 0; i + 1 < N; ++i)
            x[i] = cs.x0[i] * (1 - xw);
        x[N - 1] = xw;
    } else if (N == 2) {
        const double v = r.u();
        double x1 = v < 0.6 ? r.u() : (v < 0.8 ? r.logu(1e-6, 0.1) : 1 - r.logu(1e-6, 0.1));
        x = {x1, 1 - x1};
    } else {
        const double a = r.u() < 0.5 ? 1.0 : 0.2;  // uniform on the simplex, or sparse
        for (auto& v : x)
            v = std::max(r.gamma(a), 1e-300);
    }
    double s = 0;
    for (double v : x)
        s += v;
    for (auto& v : x)
        v /= s;
}

// ---------------------------------------------------------------- EOS evaluation conditioning
// Roundoff bound on Z from a direct double-precision term sum: each term n tau^t delta^d exp(u)
// carries relative error ~ eps (1 + |log of its magnitude|) (exp amplifies the argument's
// rounding), so dZ ~ eps sum_k |term_k| (1 + |ln|term_k||).  Returned relative to |Z|.
double z_rel_uncert(const Solver& sv, const std::vector<double>& x, double T, double D) {
    const double tau = sv.Red->Tr(x) / T, lt = std::log(tau);
    double Z = 1, u = 0;
    for (const auto& tm : sv.terms) {
        const double X = tm.j < 0 ? x[tm.i] : x[tm.i] * x[tm.j] * sv.F[tm.i][tm.j];
        const double v = X * tm.kappa(tau, lt) * tm.chi(D);
        Z += v;
        u += std::abs(v) * (1 + std::abs(std::log(std::abs(tm.kappa(tau, lt)) + 1e-300)));
    }
    return 2.2e-16 * u / std::max(std::abs(Z), 1e-300);
}

// ---------------------------------------------------------------- dense scan of the true equation
int scan_roots(const Solver& sv, const Solver::State& S, double t, double* D, int maxn) {
    auto G = [&](double d) {
        double g, dg, sc;
        sv.true_G(S, d, t, g, dg, sc);
        return g;
    };
    int n = 0;
    double d0 = 1e-12, g0 = G(d0);
    auto step = [&](double d1) {
        const double g1 = G(d1);
        if ((g0 < 0) != (g1 < 0) && n < maxn) {
            double a = d0, b = d1, fa = g0;
            for (int it = 0; it < 64; ++it) {
                const double m = 0.5 * (a + b), fm = G(m);
                if ((fm < 0) == (fa < 0)) {
                    a = m;
                    fa = fm;
                } else
                    b = m;
            }
            D[n++] = 0.5 * (a + b);
        }
        d0 = d1;
        g0 = g1;
    };
    for (int j = 1; j <= 1000; ++j)
        step(std::pow(10.0, -10 + 8.0 * j / 1000));  // 1e-10 .. 1e-2, log-spaced
    for (int j = 1; j <= 20000; ++j)
        step(1e-2 + (sv.delta_max - 1e-2) * j / 20000);
    return n;
}

struct Opts
{
    double calls = 2e5, tol = 1e-6;
    int threads = 8, spin = 0;
    uint64_t seed = 1;
    long scan_every = 20000, cp_every = 200000;
    std::string log = "gerg_validate_mismatch.log", suites = "binaries,multi,natgas,asym,humidair,humidair_ref,natgas_ref";
};

void run_call(CompSet& cs, Rng& r, const Opts& o, Stats& st, Solver::State& S, std::vector<double>& x, long callid) {
    const Solver& sv = cs.sv;
    sample_x(cs, r, x);
    const double T = r.logu(cs.Tmin, cs.Tmax), p = r.logu(1e2, 1e8);
    Solver::Root roots[64];
    const long na0 = g_na_active;
    const long nodes0 = g_nodes, unc0 = g_uncertain, pit0 = g_pol_its, pc0 = g_pol_calls;
    const auto t0 = std::chrono::steady_clock::now();
    sv.assemble(T, x, S);
    const int n = sv.roots(S, p, roots);
    st.ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
    st.nodes += g_nodes - nodes0;
    st.na_active += g_na_active - na0;
    st.pol_its += g_pol_its - pit0;
    st.pol_calls += g_pol_calls - pc0;
    ++st.calls;
    st.roots += n;
    st.multi_root_calls += n > 1;
    const double t = p * S.t_scale;
    if (g_uncertain > unc0) {
        st.uncertain += g_uncertain - unc0;
        logmsg("UNCERTAIN %s %s T=%.17g p=%.17g x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(), T, p, x[0]);
    }
    // parity
    double Ge, dGe, sce;
    sv.true_G(S, sv.delta_max, t, Ge, dGe, sce);
    if (Ge < 0) ++st.beyond;
    if ((n % 2 == 1) != (Ge > 0)) {
        ++st.parity_fail;
        logmsg("PARITY %s %s T=%.17g p=%.17g n=%d G(dmax)=%.3e x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(), T, p, n, Ge, x[0]);
    }
    for (int k = 0; k < n; ++k) {
        const auto& rt = roots[k];
        st.max_resid = std::max(st.max_resid, rt.resid);
        st.uncertified += !rt.certified;
        st.bracket_bad += !rt.bracket_ok;
        const double mv = std::abs(rt.rho - rt.rho_cheb) / rt.rho_cheb;
        st.max_move = std::max(st.max_move, mv);
        if (mv > 1e-5 || rt.resid > 1e-10) {
            ++st.polish_moved;
            logmsg("POLISH %s %s T=%.17g p=%.17g root %d/%d move=%.2e resid=%.2e cert=%d br=%d\n", cs.suite.c_str(), cs.name.c_str(), T, p, k + 1, n,
                   mv, rt.resid, rt.certified, rt.bracket_ok);
        }
    }
    // dense scan
    if (o.scan_every > 0 && callid % o.scan_every == 0) {
        ++st.scans;
        double D[256];
        const int m = scan_roots(sv, S, t, D, 256);
        std::vector<bool> used(n, false);
        bool ok = true;
        for (int i = 0; i < m; ++i) {  // every scan root must be one of ours
            int best = -1;
            double bd = 1e300;
            for (int k = 0; k < n; ++k) {
                const double d = std::abs(roots[k].rho / S.rhor - D[i]) / D[i];
                if (!used[k] && d < bd) {
                    bd = d;
                    best = k;
                }
            }
            if (best >= 0 && bd < 1e-8) {
                used[best] = true;
                st.scan_maxrel = std::max(st.scan_maxrel, bd);
            } else {
                ++st.solver_missed;
                ok = false;
                logmsg("MISSED %s %s T=%.17g p=%.17g scan root delta=%.12g (ours %d, scan %d) x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(), T, p,
                       D[i], n, m, x[0]);
            }
        }
        for (int k = 0; k < n; ++k)
            if (!used[k]) {  // ours but not the scan's: a real sign change?
                const double d = roots[k].rho / S.rhor;
                double gl, gr, dg, sc;
                sv.true_G(S, d * (1 - 1e-7), t, gl, dg, sc);
                sv.true_G(S, d * (1 + 1e-7), t, gr, dg, sc);
                if ((gl < 0) != (gr < 0))
                    ++st.scan_missed_real;  // close pair inside one scan cell
                else {
                    ++st.solver_spurious;
                    ok = false;
                    logmsg("SPURIOUS %s %s T=%.17g p=%.17g delta=%.12g x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(), T, p, d, x[0]);
                }
            }
        st.scan_ok += ok;
    }
    // CoolProp cross-checks
    if (o.cp_every > 0 && callid % o.cp_every == 1) {
        std::lock_guard<std::mutex> lk(cs.cp_mtx);
        auto* A = cs.AScp.get();
        A->set_mole_fractions(x);
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(A);
        ++st.cp_tries;
        try {
            const double rho = H->solver_rho_Tp(T, p);
            if (std::isfinite(rho) && rho > 0 && rho / S.rhor <= sv.delta_max) {
                ++st.cp_converged;
                double bd = 1e300;
                for (int k = 0; k < n; ++k)
                    bd = std::min(bd, std::abs(roots[k].rho - rho) / rho);
                if (bd < 1e-7)
                    ++st.cp_in_set;
                else {
                    ++st.cp_not_in_set;
                    logmsg("CP_NOT_IN_SET %s %s T=%.17g p=%.17g cp_rho=%.12g nearest_rel=%.2e n=%d x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(), T,
                           p, rho, bd, n, x[0]);
                }
            }
        } catch (...) {
        }
        for (int k = 0; k < n; ++k) {  // CoolProp's own pressure at each of our roots
            try {
                if (z_rel_uncert(sv, x, T, roots[k].rho / S.rhor) > 1e-8) {  // EOS itself not resolvable here in double
                    ++st.cp_illcond;
                    continue;
                }
                A->specify_phase(CoolProp::iphase_gas);
                A->update(CoolProp::DmolarT_INPUTS, roots[k].rho, T);
                // compare as a density error: |p_CP(rho) - p| / (rho dp/drho).  A pressure difference alone is
                // amplified by rho (dp/drho) / p, which reaches ~1e7 for liquids at low pressure.
                const double dp =
                  std::abs(A->p() - p) / (roots[k].rho * std::abs(A->first_partial_deriv(CoolProp::iP, CoolProp::iDmolar, CoolProp::iT)));
                st.cp_pmax = std::max(st.cp_pmax, dp);
                if (dp > 1e-9)
                    logmsg("CP_PRESSURE %s %s T=%.17g p=%.17g rho=%.12g cp_p=%.12g drho_rel=%.2e n=%d x0=%.17g\n", cs.suite.c_str(), cs.name.c_str(),
                           T, p, roots[k].rho, A->p(), dp, n, x[0]);
                ++st.cp_pchecks;
            } catch (...) {
            }
            A->unspecify_phase();
        }
    }
}

// ---------------------------------------------------------------- spinodal stress test
struct SpinStats
{
    std::map<int, std::array<long, 7>> byexp;  // eps exponent -> {cases, agree, missed, spurious, uncertain, missed_unflagged, below_EOS_resolution}
};
void spin_case(CompSet& cs, Rng& r, SpinStats& ss, Solver::State& S, std::vector<double>& x) {
    const Solver& sv = cs.sv;
    sample_x(cs, r, x);
    const double T = r.logu(cs.Tmin, cs.Tmax);
    sv.assemble(T, x, S);
    // extrema of F = delta Z: roots of dF/d(delta) from the tables, then Newton on the true F'
    std::vector<double> ext;
    for (int pc = 0; pc < sv.P; ++pc) {
        const VecG d = chebder(S.G[pc], NG);
        double sa = 0;
        for (double v : d)
            sa += std::abs(v);
        RootOut ro[MAXROOTS];
        const int nr = bern_roots(d, 1e-12 * sa, ro);
        for (int k = 0; k < nr; ++k)
            ext.push_back(sv.edges[pc] + (sv.edges[pc + 1] - sv.edges[pc]) * (ro[k].u + 1) / 2);
    }
    for (double D : ext) {
        double G, dG, sc, Gp, dGp, Gm, dGm;
        for (int it = 0; it < 30; ++it) {
            const double h = 1e-6 * D;
            sv.true_G(S, D, 0, G, dG, sc);
            sv.true_G(S, D + h, 0, Gp, dGp, sc);
            sv.true_G(S, D - h, 0, Gm, dGm, sc);
            const double step = dG / ((dGp - dGm) / (2 * h));
            D -= step;
            if (std::abs(step) < 1e-15 * D) break;
        }
        double Fs, dF, F2p, F2m;
        sv.true_G(S, D, 0, Fs, dF, sc);
        const double h = 1e-5 * D;
        sv.true_G(S, D + h, 0, G, F2p, sc);
        sv.true_G(S, D - h, 0, G, F2m, sc);
        const double F2 = (F2p - F2m) / (2 * h);
        if (!(Fs > 0) || !(std::abs(F2) > 0) || !std::isfinite(Fs)) continue;
        const double zu = z_rel_uncert(sv, x, T, D);  // relative uncertainty of Z (= of F) at the spinodal
        for (int e = 2; e <= 10; ++e)
            for (int s : {-1, 1}) {
                const double eps = std::pow(10.0, -e), t = Fs * (1 + s * eps), p = t / S.t_scale;
                if (zu > 0.1 * eps) {  // the EOS cannot resolve this eps here: no ground truth
                    ++ss.byexp[e][6];
                    continue;
                }
                const double w = std::max(20 * std::sqrt(2 * std::abs(t - Fs) / std::abs(F2)), 1e-9 * D);
                if (D + w > sv.delta_max || D - w <= 0) continue;  // window must lie inside the solver's domain
                // reference: local scan of the true G
                int nref = 0;
                double g0;
                sv.true_G(S, D - w, t, g0, dG, sc);
                for (int j = 1; j <= 4000; ++j) {
                    double g1;
                    sv.true_G(S, D - w + 2 * w * j / 4000, t, g1, dG, sc);
                    nref += (g0 < 0) != (g1 < 0);
                    g0 = g1;
                }
                const long unc0 = g_uncertain;
                Solver::Root roots[64];
                const int n = sv.roots(S, p, roots);
                int nloc = 0;
                for (int k = 0; k < n; ++k)
                    nloc += std::abs(roots[k].rho / S.rhor - D) <= w;
                auto& a = ss.byexp[e];
                ++a[0];
                if (nloc == nref)
                    ++a[1];
                else if (nloc < nref)
                    ++a[2];
                else
                    ++a[3];
                a[4] += g_uncertain > unc0;
                a[5] += nloc < nref && !(g_uncertain > unc0);
                if (nloc != nref)
                    logmsg("SPIN %s %s T=%.17g eps=1e-%d s=%d delta_s=%.17g ours=%d ref=%d uncertain=%d x0=%.17g Fs=%.17g\n", cs.suite.c_str(),
                           cs.name.c_str(), T, e, s, D, nloc, nref, int(g_uncertain > unc0), x[0], Fs);
            }
    }
}
}  // namespace

int main(int argc, char** argv) {
    Opts o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto eq = a.find('=');
        const std::string k = a.substr(0, eq), v = eq == std::string::npos ? "" : a.substr(eq + 1);
        if (k == "calls") o.calls = std::atof(v.c_str());
        if (k == "threads") o.threads = std::atoi(v.c_str());
        if (k == "seed") o.seed = std::strtoull(v.c_str(), nullptr, 10);
        if (k == "tol") o.tol = std::atof(v.c_str());
        if (k == "scan_every") o.scan_every = std::atol(v.c_str());
        if (k == "cp_every") o.cp_every = std::atol(v.c_str());
        if (k == "spin") o.spin = std::atoi(v.c_str());
        if (k == "log") o.log = v;
        if (k == "suites") o.suites = v;
    }
    g_logfile = std::fopen(o.log.c_str(), "w");
    std::printf("degree %d, table tol %.0e, delta_max 4, %.3g calls, %d threads, seed %llu\n", NQ, o.tol, o.calls, o.threads,
                (unsigned long long)o.seed);

    // ---- component sets
    std::vector<std::unique_ptr<CompSet>> sets;
    auto add = [&](const std::string& suite, const std::string& name, std::vector<std::string> fl, std::vector<double> x0 = {},
                   const std::string& backend = "GERG2008") {
        if (("," + o.suites + ",").find("," + suite + ",") == std::string::npos) return;
        auto cs = std::make_unique<CompSet>();
        cs->suite = suite;
        cs->name = name;
        cs->fluids = std::move(fl);
        cs->x0 = std::move(x0);
        cs->backend = backend;
        sets.push_back(std::move(cs));
    };
    for (std::size_t i = 0; i < GERG21.size(); ++i)
        for (std::size_t j = i + 1; j < GERG21.size(); ++j)
            add("binaries", GERG21[i] + "/" + GERG21[j], {GERG21[i], GERG21[j]});
    {
        std::mt19937_64 g(12345);
        for (int s = 0; s < 60; ++s) {
            std::vector<std::string> pool = GERG21;
            std::shuffle(pool.begin(), pool.end(), g);
            const int k = 3 + static_cast<int>(g() % 8);
            pool.resize(k);
            add("multi", "rand" + std::to_string(s) + "(" + std::to_string(k) + ")", pool);
        }
        add("multi", "all21", GERG21);
    }
    const std::vector<std::string> NG10 = {"Methane",   "Nitrogen", "CarbonDioxide", "Ethane",    "Propane",
                                           "IsoButane", "n-Butane", "Isopentane",    "n-Pentane", "n-Hexane"};
    add("natgas", "Amarillo", NG10, {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393});
    add("natgas", "Ekofisk", {NG10.begin(), NG10.begin() + 9},
        {0.859063, 0.010068, 0.014954, 0.084919, 0.023015, 0.003486, 0.003506, 0.000509, 0.00048});
    add("natgas", "GulfCoast", NG10, {0.965222, 0.002595, 0.005956, 0.018186, 0.004596, 0.000977, 0.001007, 0.000473, 0.000324, 0.000664});
    add("natgas", "HighCO2", {NG10.begin(), NG10.begin() + 7}, {0.81212, 0.05702, 0.07585, 0.04303, 0.00895, 0.00151, 0.00152});
    add("natgas", "HighN2", {NG10.begin(), NG10.begin() + 7}, {0.81441, 0.13465, 0.00985, 0.033, 0.00605, 0.001, 0.00104});
    add("natgas", "NaturalGasSample", NG10, {0.95123, 0.00089, 0.02555, 0.01835, 0.00238, 0.0004, 0.00016, 0.00014, 0.00011, 0.00079});
    add("natgas", "Air", {"Nitrogen", "Argon", "Oxygen"}, {0.7812, 0.0092, 0.2096});
    // natural gases on the reference EOS (CO2: Span-Wagner non-analytic terms)
    add("natgas_ref", "Amarillo (HEOS)", NG10, {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393},
        "HEOS");
    add("natgas_ref", "HighCO2 (HEOS)", {NG10.begin(), NG10.begin() + 7}, {0.81212, 0.05702, 0.07585, 0.04303, 0.00895, 0.00151, 0.00152}, "HEOS");
    for (auto pr : std::vector<std::pair<std::string, std::string>>{{"Methane", "HydrogenSulfide"},
                                                                    {"Methane", "n-Heptane"},
                                                                    {"Methane", "n-Octane"},
                                                                    {"Methane", "n-Nonane"},
                                                                    {"Methane", "n-Decane"},
                                                                    {"Ethane", "Water"},
                                                                    {"CarbonDioxide", "Water"},
                                                                    {"Nitrogen", "Water"},
                                                                    {"Methane", "Water"},
                                                                    {"Hydrogen", "n-Decane"},
                                                                    {"Helium", "n-Decane"},
                                                                    {"CarbonDioxide", "n-Decane"}})
        add("asym", pr.first + "/" + pr.second, {pr.first, pr.second});
    add("humidair", "humid air", {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"}, {0.7810, 0.2095, 0.0092, 0.0003, 0.0});
    // same, with the reference EOS (IAPWS-95 water and Span-Wagner CO2 carry non-analytic terms)
    add("humidair_ref", "humid air (HEOS)", {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"}, {0.7810, 0.2095, 0.0092, 0.0003, 0.0}, "HEOS");
    // water-rich (x_w 0.3..0.99), reference EOS: exercises the IAPWS-95 non-analytic add-in (mixture tau near 1)
    add("wetair_ref", "water-rich air (HEOS)", {"Nitrogen", "Oxygen", "Argon", "CarbonDioxide", "Water"}, {0.7810, 0.2095, 0.0092, 0.0003, 0.0},
        "HEOS");

    const auto tb = std::chrono::steady_clock::now();
    std::size_t tables_bytes = 0;
    int maxP = 0;
    for (auto& cs : sets) {
        cs->AS.reset(CoolProp::AbstractState::factory(cs->backend, join(cs->fluids)));
        cs->AScp.reset(CoolProp::AbstractState::factory(cs->backend, join(cs->fluids)));
        auto* H = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(cs->AS.get());
        double Tcmax = 0;
        for (const auto& c : H->get_components())
            Tcmax = std::max(Tcmax, c.EOS().reduce.T);
        cs->Tmin = std::max(60.0, 0.3 * Tcmax);
        // bound Tr(x) over compositions (sampled) for the table error budget
        double Trmax = 0;
        Rng r(7);
        std::vector<double> x;
        for (int s = 0; s < 2000; ++s) {
            sample_x(*cs, r, x);
            Trmax = std::max(Trmax, H->Reducing->Tr(x));
        }
        for (std::size_t i = 0; i < cs->fluids.size(); ++i) {  // and the pure-component corners
            x.assign(cs->fluids.size(), 0.0);
            x[i] = 1;
            Trmax = std::max(Trmax, H->Reducing->Tr(x));
        }
        cs->sv.build(H, 4.0, 1.02 * Trmax / cs->Tmin, o.tol);
        tables_bytes += cs->sv.C.size() * sizeof(VecQ);
        maxP = std::max(maxP, cs->sv.P);
    }
    std::printf("%zu component sets built in %.1f s; tables %.1f MB total; max pieces %d\n", sets.size(),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - tb).count(), tables_bytes / 1048576.0, maxP);

    // ---- main run
    std::map<std::string, std::vector<CompSet*>> bysuite;
    for (auto& cs : sets)
        bysuite[cs->suite].push_back(cs.get());
    std::vector<std::string> suites;
    for (auto& kv : bysuite)
        suites.push_back(kv.first);
    const long total = static_cast<long>(o.calls);
    std::atomic<long> next{0};
    std::vector<std::map<std::string, Stats>> tstats(o.threads);
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<std::thread> th;
    for (int ti = 0; ti < o.threads; ++ti)
        th.emplace_back([&, ti]() {
            Rng r(o.seed * 1000003 + ti);
            Solver::State S;
            std::vector<double> x;
            for (;;) {
                const long base = next.fetch_add(1000);
                if (base >= total) break;
                for (long c = base; c < std::min(total, base + 1000); ++c) {
                    const std::string& su = suites[c % suites.size()];  // suites equally weighted
                    auto& v = bysuite[su];
                    CompSet& cs = *v[r.g() % v.size()];
                    run_call(cs, r, o, tstats[ti][su], S, x, c / static_cast<long>(suites.size()));
                }
            }
        });
    std::thread prog([&]() {
        while (next.load() < total) {
            for (int s = 0; s < 30 && next.load() < total; ++s)
                std::this_thread::sleep_for(std::chrono::seconds(1));
            if (next.load() >= total) break;
            const double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            const long d = std::min(next.load(), total);
            std::printf("  progress %.3g / %.3g calls, %.0f s, %.0f calls/s\n", double(d), double(total), el, d / el);
            std::fflush(stdout);
        }
    });
    for (auto& t : th)
        t.join();
    prog.join();
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::printf("\nwall %.0f s, %.0f calls/s\n\n", wall, total / wall);
    std::printf("%-9s %11s %8s %6s %7s %6s %7s %7s %7s %8s %8s %6s %9s %8s %8s %8s %9s %8s\n", "suite", "calls", "roots/c", ">1root", "parity",
                "uncert", "uncert", "brk_bad", "moved", "max_res", "us/call", "scans", "scan_miss", "missed", "spurious", "cp_in", "cp_notin",
                "cp_p");
    std::printf("%-9s %11s %8s %6s %7s %6s %7s %7s %7s %8s %8s %6s %9s %8s %8s %8s %9s %8s\n", "", "", "", "%", "fails", "flags", "roots", "", "", "",
                "", "", "(pairs)", "", "", "", "", "max rel");
    Stats all;
    for (auto& su : suites) {
        Stats s;
        for (auto& m : tstats)
            if (m.count(su)) s.add(m.at(su));
        all.add(s);
        std::printf("%-9s %11ld %8.3f %6.1f %7ld %6ld %7ld %7ld %7ld %8.1e %8.2f %6ld %9ld %8ld %8ld %8ld %9ld %8.1e\n", su.c_str(), s.calls,
                    double(s.roots) / s.calls, 100.0 * s.multi_root_calls / s.calls, s.parity_fail, s.uncertain, s.uncertified, s.bracket_bad,
                    s.polish_moved, s.max_resid, s.ns / s.calls / 1000, s.scans, s.scan_missed_real, s.solver_missed, s.solver_spurious, s.cp_in_set,
                    s.cp_not_in_set, s.cp_pmax);
    }
    std::printf("%-9s %11ld %8.3f %6.1f %7ld %6ld %7ld %7ld %7ld %8.1e %8.2f %6ld %9ld %8ld %8ld %8ld %9ld %8.1e\n", "ALL", all.calls,
                double(all.roots) / all.calls, 100.0 * all.multi_root_calls / all.calls, all.parity_fail, all.uncertain, all.uncertified,
                all.bracket_bad, all.polish_moved, all.max_resid, all.ns / all.calls / 1000, all.scans, all.scan_missed_real, all.solver_missed,
                all.solver_spurious, all.cp_in_set, all.cp_not_in_set, all.cp_pmax);
    for (auto& su : suites) {
        Stats s;
        for (auto& m : tstats)
            if (m.count(su)) s.add(m.at(su));
        if (s.na_active)
            std::printf("non-analytic add-in active in suite %s: %ld of %ld calls (%.2f %%)\n", su.c_str(), s.na_active, s.calls,
                        100.0 * s.na_active / s.calls);
    }
    std::printf("\nroots beyond delta_max (G(4) < 0): %ld calls; CoolProp solver_rho_Tp tried %ld, converged in range %ld; polish %.2f its/root; "
                "subdivision nodes %.1f/call\nCoolProp p(rho) at our roots: %ld compared (max density-equivalent error %.1e), %ld skipped as "
                "ill-conditioned (Z "
                "uncertain > 1e-8 in double)\n",
                all.beyond, all.cp_tries, all.cp_converged, double(all.pol_its) / std::max(1L, all.pol_calls), double(all.nodes) / all.calls,
                all.cp_pchecks, all.cp_pmax, all.cp_illcond);

    // ---- spinodal stress test
    if (o.spin > 0) {
        std::vector<SpinStats> ss(o.threads);
        std::atomic<long> nxt{0};
        std::vector<std::thread> ts;
        for (int ti = 0; ti < o.threads; ++ti)
            ts.emplace_back([&, ti]() {
                Rng r(o.seed * 7919 + ti);
                Solver::State S;
                std::vector<double> x;
                for (long c; (c = nxt.fetch_add(1)) < o.spin;) {
                    auto& v = bysuite[suites[c % suites.size()]];
                    spin_case(*v[r.g() % v.size()], r, ss[ti], S, x);
                }
            });
        for (auto& t : ts)
            t.join();
        std::map<int, std::array<long, 7>> tot;
        for (auto& s : ss)
            for (auto& kv : s.byexp)
                for (int k = 0; k < 7; ++k)
                    tot[kv.first][k] += kv.second[k];
        std::printf("\nspinodal stress test (%d (T,x) samples): p = p_spinodal (1 +- eps); root count within the local window vs a 4000-point scan\n",
                    o.spin);
        std::printf("%8s %9s %9s %8s %9s %10s %17s %14s\n", "eps", "cases", "agree", "missed", "spurious", "uncertain", "missed&unflagged",
                    "below EOS res.");
        for (auto& kv : tot)
            std::printf("%8s %9ld %9ld %8ld %9ld %10ld %17ld %14ld\n", ("1e-" + std::to_string(kv.first)).c_str(), kv.second[0], kv.second[1],
                        kv.second[2], kv.second[3], kv.second[4], kv.second[5], kv.second[6]);
    }
    if (g_logfile) std::fclose(g_logfile);
}
