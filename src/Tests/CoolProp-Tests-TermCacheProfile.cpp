// Standalone profiling harness for residual-Helmholtz term-table cache feasibility
// (optimization "E" reconnaissance for ResidualHelmholtzGeneralizedExponential::all).
//
// QUESTION:
//   For real CoolProp fluids, how many UNIQUE l_double / m_double exponents appear
//   in the residual term table vs the total term count N?  If the sharing ratio is
//   high (say 5 unique across 50 terms), caching pow(delta, l) / pow(tau, m) at the
//   top of the loop replaces N exp/log evals with just (unique + N lookups).  If the
//   ratio is low (45 unique across 50 terms), the cache buys nothing.
//
// PART 1 — Term-table structure
//   For each representative fluid: walk components[0].EOS().alphar.GenExp.elements
//   and count distinct l_double / m_double values, filtered by the gates the live
//   all() implementation actually uses (delta_li_in_u + c!=0 + l>0 for l;
//   tau_mi_in_u + |m|>0 for m).  Print a per-fluid table with the savings ratio.
//
// PART 2 — Microbench of the two evaluation paths
//   For Water, R134a, n-Propane: at a representative single-phase state, run a tight
//   loop that mimics the per-term pow/exp pattern.  Path A repeats exp(l*log_delta)
//   and exp(m*log_tau) per term.  Path B precomputes a small vector of distinct
//   values and looks them up via linear search.  std::chrono steady_clock, 1e6 reps,
//   5 trials, report mean +/- stddev for each path and the A/B ratio.
//
// This file is a PROFILING TOOL — Catch2 hidden ([.]).  It does NOT modify
// production code; it only INSPECTS public state.  Run via:
//   ./build_catch_rel/CatchTestRunner "[term_cache_profile]"

#include "CoolProp/AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/CoolPropFluid.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/fluids/Helmholtz.h"

#include <sstream>

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <cfloat>
#    include <chrono>
#    include <cmath>
#    include <cstddef>
#    include <cstdio>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

// Count distinct values in a vector under a relative tolerance.  Order-preserving:
// emits the first occurrence each time.  Used by Part 1 and Part 2 (to build the
// "small linear-search table" the cached path uses).
std::vector<double> distinct_values(const std::vector<double>& xs, double rtol = 1e-12) {
    std::vector<double> out;
    out.reserve(xs.size());
    for (double v : xs) {
        bool seen = false;
        for (double u : out) {
            const double scale = std::max(std::abs(u), 1.0);
            if (std::abs(u - v) <= rtol * scale) {
                seen = true;
                break;
            }
        }
        if (!seen) out.push_back(v);
    }
    return out;
}

// What the term table looks like from the public API.  Mirrors the gates inside
// ResidualHelmholtzGeneralizedExponential::all so the unique-count we report
// corresponds to the pow/exp calls actually executed at run time.
struct TermStats
{
    std::size_t N = 0;         // total residual term count (GenExp.elements.size())
    std::size_t N_l_used = 0;  // terms that take the exp(l*log_delta) branch
    std::size_t N_m_used = 0;  // terms that take the exp(m*log_tau) branch
    std::size_t unique_l = 0;  // distinct l_double values among the N_l_used terms
    std::size_t unique_m = 0;  // distinct m_double values among the N_m_used terms
    std::vector<double> ls;    // the l_double values per-term (for Part 2 to reuse)
    std::vector<double> ms;    // the m_double values per-term (for Part 2 to reuse)
    bool delta_li_in_u = false;
    bool tau_mi_in_u = false;
    bool ok = false;  // false => introspection failed / not a HEOS pure
};

// Walk the residual term table of fluid `name` and produce stats.  Returns ok=false
// if the fluid is unavailable, not a HEOS pure backend, or otherwise unscrutable —
// caller will report it as "skipped".
TermStats inspect_fluid(const std::string& name) {
    TermStats s;
    std::shared_ptr<CoolProp::AbstractState> AS;
    try {
        AS.reset(CoolProp::AbstractState::factory("HEOS", name));
    } catch (...) {
        return s;
    }
    auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
    if (be == nullptr) return s;
    auto& comps = be->get_components();
    if (comps.empty()) return s;
    auto& gen = comps[0].EOS().alphar.GenExp;
    s.delta_li_in_u = gen.delta_li_in_u;
    s.tau_mi_in_u = gen.tau_mi_in_u;
    s.N = gen.elements.size();

    std::vector<double> active_l, active_m;
    active_l.reserve(s.N);
    active_m.reserve(s.N);
    for (std::size_t i = 0; i < s.N; ++i) {
        const auto& el = gen.elements[i];
        // Gate for the delta^l branch — mirrors the live code in
        // src/Helmholtz.cpp around line 162 of ResidualHelmholtzGeneralizedExponential::all.
        const bool l_active = s.delta_li_in_u && std::isfinite(static_cast<double>(el.l_double)) && static_cast<double>(el.l_double) > 0.0
                              && std::abs(static_cast<double>(el.c)) > DBL_EPSILON;
        // Gate for the tau^m branch — mirrors the live code around line 177.
        const bool m_active = s.tau_mi_in_u && std::abs(static_cast<double>(el.m_double)) > 0.0;
        if (l_active) active_l.push_back(static_cast<double>(el.l_double));
        if (m_active) active_m.push_back(static_cast<double>(el.m_double));
    }
    s.N_l_used = active_l.size();
    s.N_m_used = active_m.size();
    s.ls = active_l;
    s.ms = active_m;
    s.unique_l = distinct_values(active_l).size();
    s.unique_m = distinct_values(active_m).size();
    s.ok = true;
    return s;
}

// ----- Part 2 helpers -----------------------------------------------------

// Path A: independently evaluate exp(l_i * log_delta) and exp(m_i * log_tau) for
// every term, the way the live all() does today.  Returns a sink to defeat DCE.
double path_A(const std::vector<double>& ls, const std::vector<double>& ms, double log_delta, double log_tau) {
    double acc = 0.0;
    const std::size_t Nl = ls.size(), Nm = ms.size();
    for (std::size_t i = 0; i < Nl; ++i)
        acc += std::exp(ls[i] * log_delta);
    for (std::size_t i = 0; i < Nm; ++i)
        acc += std::exp(ms[i] * log_tau);
    return acc;
}

// Path B: precompute one exp() per unique exponent, then do a small linear-search
// lookup per term to fetch the cached value.  Realistic for the typical case where
// unique <= 16 — linear search beats a hash map at that size.
double path_B(const std::vector<double>& ls, const std::vector<double>& ms, const std::vector<double>& uniq_l, const std::vector<double>& uniq_m,
              double log_delta, double log_tau) {
    // Precompute exp(u * log_x) once per unique exponent.
    double cache_l[64], cache_m[64];
    const std::size_t Ul = std::min<std::size_t>(uniq_l.size(), 64), Um = std::min<std::size_t>(uniq_m.size(), 64);
    for (std::size_t k = 0; k < Ul; ++k)
        cache_l[k] = std::exp(uniq_l[k] * log_delta);
    for (std::size_t k = 0; k < Um; ++k)
        cache_m[k] = std::exp(uniq_m[k] * log_tau);

    double acc = 0.0;
    const std::size_t Nl = ls.size(), Nm = ms.size();
    for (std::size_t i = 0; i < Nl; ++i) {
        const double li = ls[i];
        for (std::size_t k = 0; k < Ul; ++k) {
            if (uniq_l[k] == li) {  // exponents loaded from JSON are bit-exact within the table
                acc += cache_l[k];
                break;
            }
        }
    }
    for (std::size_t i = 0; i < Nm; ++i) {
        const double mi = ms[i];
        for (std::size_t k = 0; k < Um; ++k) {
            if (uniq_m[k] == mi) {
                acc += cache_m[k];
                break;
            }
        }
    }
    return acc;
}

struct TrialStats
{
    double mean_ns = 0.0;
    double stddev_ns = 0.0;
};

// Run `n_trials` of `reps` repetitions of `fn`, return mean and stddev of ns/iter.
template <typename F>
TrialStats bench(F&& fn, std::size_t reps, std::size_t n_trials) {
    std::vector<double> ns_per_iter;
    ns_per_iter.reserve(n_trials);
    volatile double sink = 0.0;
    for (std::size_t t = 0; t < n_trials; ++t) {
        const auto t0 = std::chrono::steady_clock::now();
        for (std::size_t r = 0; r < reps; ++r) {
            sink = sink + fn();
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double dt_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        ns_per_iter.push_back(dt_ns / static_cast<double>(reps));
    }
    (void)sink;
    double sum = 0.0;
    for (double x : ns_per_iter)
        sum += x;
    const double mean = sum / static_cast<double>(n_trials);
    double sq = 0.0;
    for (double x : ns_per_iter)
        sq += (x - mean) * (x - mean);
    const double stddev = (n_trials > 1) ? std::sqrt(sq / static_cast<double>(n_trials - 1)) : 0.0;
    return TrialStats{mean, stddev};
}

// Pick a representative single-phase (T, p) state for each benchmark fluid, run the
// AbstractState through update so tau/delta are well-defined, then return log(tau)
// and log(delta).  Returns false if the state can't be reached.
bool eval_state(const std::string& name, double T, double p, double& log_tau, double& log_delta) {
    std::shared_ptr<CoolProp::AbstractState> AS;
    try {
        AS.reset(CoolProp::AbstractState::factory("HEOS", name));
        AS->update(CoolProp::PT_INPUTS, p, T);
    } catch (...) {
        return false;
    }
    const double tau = AS->tau(), delta = AS->delta();
    if (!(tau > 0) || !(delta > 0)) return false;
    log_tau = std::log(tau);
    log_delta = std::log(delta);
    return true;
}

}  // anonymous namespace

TEST_CASE("Tau-mi sweep: which fluids actually populate tau_mi_in_u with non-zero m_double", "[term_cache_tau_sweep][.]") {
    // Scan the whole fluid library for fluids where N_m_used > 0.
    const std::string fluids_csv = CoolProp::get_global_param_string("FluidsList");
    std::vector<std::string> all_fluids;
    std::stringstream ss(fluids_csv);
    std::string item;
    while (std::getline(ss, item, ','))
        all_fluids.push_back(item);
    std::printf("\n=== tau_mi sweep: %zu fluids ===\n", all_fluids.size());
    std::printf("%-30s %5s %8s %10s %8s %10s\n", "fluid", "N", "N_l_used", "unique_l", "N_m_used", "unique_m");
    std::size_t with_m = 0;
    for (const auto& name : all_fluids) {
        const TermStats s = inspect_fluid(name);
        if (!s.ok) continue;
        if (s.N_m_used == 0) continue;
        ++with_m;
        std::printf("%-30s %5zu %8zu %10zu %8zu %10zu\n", name.c_str(), s.N, s.N_l_used, s.unique_l, s.N_m_used, s.unique_m);
    }
    std::printf("--- %zu / %zu fluids have N_m_used > 0 ---\n", with_m, all_fluids.size());
}

TEST_CASE("Term-cache feasibility profile for ResidualHelmholtzGeneralizedExponential::all", "[term_cache_profile][.]") {
    // -------------------- PART 1 --------------------
    const std::vector<std::string> fluids = {"Water", "CarbonDioxide", "n-Propane", "R134a", "Nitrogen", "Methane",
                                             "MM",    "Ammonia",       "Hydrogen",  "R32",   "R143a"};

    std::printf("\n=== Part 1: residual-term-table structure (per fluid) ===\n");
    std::printf("%-15s %5s %8s %10s %8s %10s %12s\n", "fluid", "N", "N_l_used", "unique_l", "N_m_used", "unique_m", "savings_pow");
    std::printf("%-15s %5s %8s %10s %8s %10s %12s\n", "-----", "-", "--------", "--------", "--------", "--------", "-----------");
    std::vector<std::string> skipped;
    for (const auto& name : fluids) {
        const TermStats s = inspect_fluid(name);
        if (!s.ok) {
            skipped.push_back(name);
            continue;
        }
        const std::size_t denom = s.N_l_used + s.N_m_used;
        const double savings = (denom > 0) ? static_cast<double>(denom - s.unique_l - s.unique_m) / static_cast<double>(denom) : 0.0;
        std::printf("%-15s %5zu %8zu %10zu %8zu %10zu %12.3f\n", name.c_str(), s.N, s.N_l_used, s.unique_l, s.N_m_used, s.unique_m, savings);
        // Sanity asserts (validation): the unique count can't exceed the participating count,
        // and savings is in [0, 1].
        REQUIRE(s.unique_l <= s.N_l_used);
        REQUIRE(s.unique_m <= s.N_m_used);
        REQUIRE(savings >= 0.0);
        REQUIRE(savings <= 1.0);
    }
    if (!skipped.empty()) {
        std::printf("(skipped: ");
        for (std::size_t i = 0; i < skipped.size(); ++i)
            std::printf("%s%s", skipped[i].c_str(), i + 1 == skipped.size() ? "" : ", ");
        std::printf(")\n");
    }

    // -------------------- PART 2 --------------------
    struct BenchFluid
    {
        std::string name;
        double T;
        double p;
    };
    // Representative single-phase states (well away from saturation / criticality).
    const std::vector<BenchFluid> benches = {
      {"Water", 400.0, 1.0e7},     // compressed liquid water, comfortably subcritical
      {"R134a", 350.0, 5.0e5},     // superheated R134a vapor
      {"n-Propane", 300.0, 2.0e5}  // gas-phase n-propane
    };
    const std::size_t reps = 1'000'000;
    const std::size_t n_trials = 5;

    std::printf("\n=== Part 2: microbench Path A (per-term exp) vs Path B (precompute+lookup) ===\n");
    std::printf("(%zu reps x %zu trials, ns per loop iteration)\n", reps, n_trials);
    std::printf("%-12s %5s %10s %10s %14s %14s %8s %14s\n", "fluid", "N", "uniq_l", "uniq_m", "A mean (ns)", "B mean (ns)", "A/B", "Delta x 25 ns");
    std::printf("%-12s %5s %10s %10s %14s %14s %8s %14s\n", "-----", "-", "------", "------", "-----------", "-----------", "---", "-------------");

    for (const auto& b : benches) {
        TermStats s = inspect_fluid(b.name);
        if (!s.ok || (s.N_l_used == 0 && s.N_m_used == 0)) {
            std::printf("%-12s [skipped — introspection failed or no l/m terms]\n", b.name.c_str());
            continue;
        }
        double log_tau = 0.0, log_delta = 0.0;
        if (!eval_state(b.name, b.T, b.p, log_tau, log_delta)) {
            std::printf("%-12s [skipped — state (%g K, %g Pa) unreachable]\n", b.name.c_str(), b.T, b.p);
            continue;
        }
        const std::vector<double> uniq_l = distinct_values(s.ls);
        const std::vector<double> uniq_m = distinct_values(s.ms);

        // Warmup once to populate caches and let the branch predictor settle.
        (void)path_A(s.ls, s.ms, log_delta, log_tau);
        (void)path_B(s.ls, s.ms, uniq_l, uniq_m, log_delta, log_tau);

        const TrialStats A = bench([&] { return path_A(s.ls, s.ms, log_delta, log_tau); }, reps, n_trials);
        const TrialStats B = bench([&] { return path_B(s.ls, s.ms, uniq_l, uniq_m, log_delta, log_tau); }, reps, n_trials);

        const double ratio = (B.mean_ns > 0) ? A.mean_ns / B.mean_ns : 0.0;
        // Rough projected solve-time saving: PT_INPUTS does on the order of ~25 EOS
        // evals per solve (residual + derivatives in the Newton iteration); each
        // residual eval calls all() once.  So saving (A - B) ns per all() call
        // translates to ~25 * (A - B) ns per flash solve.
        const double per_solve_ns = (A.mean_ns - B.mean_ns) * 25.0;

        std::printf("%-12s %5zu %10zu %10zu %8.1f +-%3.1f %8.1f +-%3.1f %8.2f %14.1f\n", b.name.c_str(), s.N, uniq_l.size(), uniq_m.size(), A.mean_ns,
                    A.stddev_ns, B.mean_ns, B.stddev_ns, ratio, per_solve_ns);
    }

    std::printf("\nNote: Path A mimics 'exp(l*log_delta) + exp(m*log_tau) per term'.\n");
    std::printf("      Path B precomputes one exp per UNIQUE exponent + linear-search lookup per term.\n");
    std::printf("      If A/B >> 1 with low unique counts, optimization (E) is worth implementing.\n");
}

#endif  // ENABLE_CATCH
