// P+X flash preconditioner — standalone unit (CoolProp-cdj prep, Phase 3a).
// Header: include/CoolProp/px_preconditioner.h.  See header for API + intent.
//
// Storage model: per-fluid MeltCurveTable is held in a static map keyed by the
// fluid's canonical name, guarded by a mutex.  This keeps the change local —
// HelmholtzEOSMixtureBackend's header does not need to grow a new member —
// and is robust to heap-slot reuse (consecutive AbstractState allocations
// can land at the same address, which would otherwise serve stale entries
// when keying by backend pointer).
// Lifetime: entries stay in the map until process exit.  The map is sized in
// proportion to the number of distinct fluid names that get queried, which
// is at most ~120 (CoolProp's pure-fluid list) for typical use.
//
// Thread safety: lazy-init takes a write lock around the build (cheap, only
// runs once per backend); the seed function takes a read lock to fetch the
// (already-built) table and otherwise operates on locals.

#include "CoolProp/px_preconditioner.h"

#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp.h"
#include "CoolPropFluid.h"
#include "DataStructures.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace CoolProp {
namespace pxprecond {

namespace {

// ---------------------------------------------------------------------------
// Per-fluid table storage.  Keyed by canonical fluid name (the first entry of
// HelmholtzEOSMixtureBackend::fluid_names() for the pure-fluid case we care
// about).  Guarded by a mutex.  The table is small (a handful of vectors);
// we never erase entries — the map grows monotonically with distinct fluids.
//
// Why fluid-name keyed (not pointer-keyed): consecutive AbstractState
// allocations can land at the same heap address, which would cause a
// pointer-keyed map to return stale entries for what is logically a different
// backend.  Fluid names are stable across backend recreation.
// ---------------------------------------------------------------------------
using table_key_t = std::string;

std::unordered_map<table_key_t, MeltCurveTable>& table_map() {
    static std::unordered_map<table_key_t, MeltCurveTable> m;
    return m;
}
std::mutex& table_mutex() {
    static std::mutex m;
    return m;
}

// Build the key (= first fluid name) for H.  Empty string on failure.
table_key_t make_key(HelmholtzEOSMixtureBackend& H) {
    try {
        const auto names = H.fluid_names();
        if (!names.empty()) return names[0];
    } catch (...) {
    }
    return {};
}

constexpr std::size_t kMeltSamples = 30;       ///< target number of (p_i, T_i, rho_i) sample triples
constexpr std::size_t kMeltMinSamples = 5;     ///< minimum successful samples required to mark has_data=true
constexpr double kSupercritIdealCutoff = 0.8;  ///< T > kSupercritIdealCutoff * Tmax  -> ideal-gas seed
constexpr double kSupercritLowPCutoff = 0.01;  ///< p < kSupercritLowPCutoff * pc    -> ideal-gas seed
constexpr double kVaporIdealCutoff = 0.5;      ///< T > kVaporIdealCutoff * Tmax     -> ideal-gas vapor seed
constexpr double kVaporLowPCutoff = 0.1;       ///< p < kVaporLowPCutoff * p_sat(T)  -> ideal-gas vapor seed
constexpr double kCompExtrapFactor = 10.0;     ///< stiff-liquid compressibility-extrap conservative factor

// Safe finite-positive check
inline bool finite_pos(double x) {
    return std::isfinite(x) && x > 0.0;
}

// Build the melt-curve table for H.  Assumes caller holds the table_mutex.
// On failure (no melt line, build error, etc.) leaves has_data=false but sets
// initialized=true so we don't retry forever.
void build_table_locked(HelmholtzEOSMixtureBackend& H, MeltCurveTable& tbl) {
    tbl.has_data = false;
    tbl.initialized = true;
    tbl.p.clear();
    tbl.T_melt.clear();
    tbl.rho_melt.clear();

    if (!H.has_melting_line()) {
        return;
    }

    // Determine melt-line p range.  The ancillary exposes pmin/pmax via the
    // iP_min / iP_max OF codes; the given arg is unused for these queries.
    double p_lo = 0, p_hi = 0;
    try {
        p_lo = H.calc_melting_line(iP_min, iT, 0.0);
        p_hi = H.calc_melting_line(iP_max, iT, 0.0);
    } catch (...) {
        return;
    }
    if (!finite_pos(p_lo) || !finite_pos(p_hi) || !(p_hi > p_lo)) {
        return;
    }

    // Defensive caps: keep the upper bound below pmax and below
    // 100 * p_critical (avoids absurd extrapolations).  Keep the lower bound
    // at p_triple if available (small numerical guard above triple_p).
    double pc = 0, pmax_eos = 0, p_triple_eos = 0;
    try {
        pc = H.p_critical();
    } catch (...) {
        pc = 0;
    }
    try {
        pmax_eos = H.pmax();
    } catch (...) {
        pmax_eos = 0;
    }
    try {
        p_triple_eos = H.p_triple();
    } catch (...) {
        p_triple_eos = 0;
    }
    if (finite_pos(p_triple_eos)) {
        p_lo = std::max(p_lo, p_triple_eos);
    }
    if (finite_pos(pmax_eos)) {
        p_hi = std::min(p_hi, pmax_eos);
    }
    if (finite_pos(pc)) {
        p_hi = std::min(p_hi, 100.0 * pc);
    }
    if (!(p_hi > p_lo)) {
        return;
    }

    // Build a fresh working state so we don't mutate H.  Use the "HEOS"
    // factory key — H is the Helmholtz mixture backend; this is a pure fluid
    // (has_melting_line() requires is_pure_or_pseudopure), so fluid_names()[0]
    // is well-defined.
    std::vector<std::string> fnames;
    try {
        fnames = H.fluid_names();
    } catch (...) {
        return;
    }
    if (fnames.empty()) {
        return;
    }
    std::shared_ptr<AbstractState> work;
    try {
        work.reset(AbstractState::factory("HEOS", fnames[0]));
    } catch (...) {
        return;
    }
    if (!work) {
        return;
    }

    // Sample log-spaced pressures.  For each, compute T_melt(p_i) from the
    // melting-line ancillary, then update(PT_INPUTS, p_i, T_melt_i) to read
    // rho.  Skip points that throw — some fluids' melting curves aren't fully
    // consistent with their EOS near the bounds.
    tbl.p.reserve(kMeltSamples);
    tbl.T_melt.reserve(kMeltSamples);
    tbl.rho_melt.reserve(kMeltSamples);
    const double log_lo = std::log(p_lo);
    const double log_hi = std::log(p_hi);
    for (std::size_t i = 0; i < kMeltSamples; ++i) {
        const double frac = (kMeltSamples == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(kMeltSamples - 1);
        const double p_i = std::exp(log_lo + frac * (log_hi - log_lo));
        double T_i = NAN, rho_i = NAN;
        try {
            T_i = H.calc_melting_line(iT, iP, p_i);
        } catch (...) {
            continue;
        }
        if (!finite_pos(T_i)) continue;
        try {
            work->update(PT_INPUTS, p_i, T_i);
            rho_i = work->rhomolar();
        } catch (...) {
            continue;
        }
        if (!finite_pos(rho_i)) continue;
        tbl.p.push_back(p_i);
        tbl.T_melt.push_back(T_i);
        tbl.rho_melt.push_back(rho_i);
    }

    if (tbl.p.size() >= kMeltMinSamples) {
        // Defend monotonicity: if any later p_i happens to be <= the previous
        // (shouldn't, but be paranoid given log-spacing+float), drop the
        // offending point.  Iterate forward; pop offenders.
        std::vector<double> p2, T2, r2;
        p2.reserve(tbl.p.size());
        T2.reserve(tbl.p.size());
        r2.reserve(tbl.p.size());
        for (std::size_t i = 0; i < tbl.p.size(); ++i) {
            if (i == 0 || tbl.p[i] > p2.back()) {
                p2.push_back(tbl.p[i]);
                T2.push_back(tbl.T_melt[i]);
                r2.push_back(tbl.rho_melt[i]);
            }
        }
        if (p2.size() >= kMeltMinSamples) {
            tbl.p = std::move(p2);
            tbl.T_melt = std::move(T2);
            tbl.rho_melt = std::move(r2);
            tbl.has_data = true;
        }
    }
}

// Snapshot of the seed-relevant invariants we need for one regime-classified
// query.  Reads them off H once so the hot path is local arithmetic.
struct SeedContext
{
    double Tc = 0, pc = 0, rho_c = 0;
    double Tmax = 0, R = 0;
    bool pure_with_superanc = false;
    std::shared_ptr<EquationOfState::SuperAncillary_t> sa;
};

SeedContext load_seed_context(HelmholtzEOSMixtureBackend& H) {
    SeedContext c;
    try {
        c.Tc = H.T_critical();
    } catch (...) {
    }
    try {
        c.pc = H.p_critical();
    } catch (...) {
    }
    try {
        c.rho_c = H.rhomolar_critical();
    } catch (...) {
    }
    try {
        c.Tmax = H.Tmax();
    } catch (...) {
    }
    try {
        c.R = H.gas_constant();
    } catch (...) {
    }
    if (H.is_pure()) {
        try {
            c.sa = H.get_superanc();
        } catch (...) {
            c.sa.reset();
        }
        c.pure_with_superanc = static_cast<bool>(c.sa);
    }
    return c;
}

}  // namespace

// ---------------------------------------------------------------------------
// MeltCurveTable::rho_at
// ---------------------------------------------------------------------------
double MeltCurveTable::rho_at(double p_query) const {
    if (!has_data || p.size() < 2 || !std::isfinite(p_query)) {
        return std::nan("");
    }
    // Clamp to table range.
    if (p_query <= p.front()) return rho_melt.front();
    if (p_query >= p.back()) return rho_melt.back();
    // Binary search for the upper bound — first p[i] > p_query.
    auto it = std::lower_bound(p.begin(), p.end(), p_query);
    if (it == p.begin()) return rho_melt.front();
    if (it == p.end()) return rho_melt.back();
    const std::size_t hi = static_cast<std::size_t>(it - p.begin());
    const std::size_t lo = hi - 1;
    const double t = (p_query - p[lo]) / (p[hi] - p[lo]);
    return rho_melt[lo] + t * (rho_melt[hi] - rho_melt[lo]);
}

// ---------------------------------------------------------------------------
// ensure_melt_curve_table / get_melt_curve_table
// ---------------------------------------------------------------------------
void ensure_melt_curve_table(HelmholtzEOSMixtureBackend& H) {
    const table_key_t key = make_key(H);
    if (key.empty()) {
        // Can't key it; nothing to do.  Subsequent get_melt_curve_table on
        // the same H will also key-empty and return a default (no-data) table.
        return;
    }
    std::lock_guard<std::mutex> lk(table_mutex());
    auto& m = table_map();
    auto it = m.find(key);
    if (it != m.end() && it->second.initialized) {
        return;
    }
    MeltCurveTable& tbl = m[key];
    if (tbl.initialized) return;  // raced; somebody else built it
    build_table_locked(H, tbl);
}

const MeltCurveTable& get_melt_curve_table(HelmholtzEOSMixtureBackend& H) {
    ensure_melt_curve_table(H);
    std::lock_guard<std::mutex> lk(table_mutex());
    auto& m = table_map();
    const table_key_t key = make_key(H);
    auto it = (!key.empty()) ? m.find(key) : m.end();
    if (it == m.end()) {
        // Return a sentinel "empty" table.  Stored in a function-local static
        // so the const reference remains valid.
        static const MeltCurveTable kEmpty{};
        return kEmpty;
    }
    return it->second;
}

// ---------------------------------------------------------------------------
// regime_classified_rho_seed
// ---------------------------------------------------------------------------
double regime_classified_rho_seed(HelmholtzEOSMixtureBackend& H, double T, double p) {
    // Ensure the melt table is built (lazy, idempotent).  We then snapshot it
    // inside the lock so the rest of the hot path can run unlocked.
    ensure_melt_curve_table(H);

    SeedContext ctx = load_seed_context(H);

    // Reasonable defaults if any constant came back bad — fall through to ideal-gas.
    auto ideal_gas = [&](double T_q, double p_q) -> double {
        if (ctx.R > 0 && T_q > 0 && std::isfinite(p_q) && p_q > 0) {
            return p_q / (ctx.R * T_q);
        }
        return std::nan("");
    };

    auto clamp01 = [](double x) {
        if (!(x >= 0.0)) return 0.0;
        if (!(x <= 1.0)) return 1.0;
        return x;
    };

    double seed = std::nan("");

    if (!finite_pos(ctx.Tc) || !finite_pos(ctx.pc) || !finite_pos(ctx.R) || !finite_pos(ctx.Tmax)) {
        // Degenerate: no critical info.  Best we can do is ideal-gas.
        seed = ideal_gas(T, p);
    } else if (T > ctx.Tc) {
        // Supercritical.
        if (T > kSupercritIdealCutoff * ctx.Tmax || p < kSupercritLowPCutoff * ctx.pc) {
            seed = ideal_gas(T, p);
        } else {
            const double alpha = clamp01((T - ctx.Tc) / ctx.Tc);
            const double rg = ideal_gas(T, p);
            const double rc = finite_pos(ctx.rho_c) ? ctx.rho_c : rg;
            if (finite_pos(rg) && finite_pos(rc)) {
                seed = (1.0 - alpha) * rc + alpha * rg;
            } else {
                seed = rg;
            }
        }
    } else {
        // Subcritical: need T_sat(p) and rho_sat,L/V(T) — requires superanc.
        if (!ctx.pure_with_superanc || !ctx.sa) {
            seed = ideal_gas(T, p);
        } else {
            EquationOfState::SuperAncillary_t& sa = *ctx.sa;
            double p_sat = NAN;
            try {
                p_sat = sa.eval_sat(T, 'P', 1);
            } catch (...) {
                p_sat = NAN;
            }
            if (!finite_pos(p_sat)) {
                seed = ideal_gas(T, p);
            } else if (p >= p_sat) {
                // Liquid branch.
                double rhoL = NAN;
                try {
                    rhoL = sa.eval_sat(T, 'D', 0);
                } catch (...) {
                    rhoL = NAN;
                }
                if (!finite_pos(rhoL)) {
                    seed = ideal_gas(T, p);
                } else {
                    // Snapshot the melt table (under lock) for the log-p interp.
                    MeltCurveTable tbl_copy;
                    {
                        std::lock_guard<std::mutex> lk(table_mutex());
                        auto& m = table_map();
                        const table_key_t key = make_key(H);
                        if (!key.empty()) {
                            auto it = m.find(key);
                            if (it != m.end()) tbl_copy = it->second;
                        }
                    }
                    if (tbl_copy.has_data && tbl_copy.p.size() >= 2) {
                        const double p_max_table = tbl_copy.p.back();
                        const double rhoM = tbl_copy.rho_at(p);
                        if (finite_pos(rhoM) && p_max_table > p_sat) {
                            const double num = std::log(std::max(p, p_sat)) - std::log(p_sat);
                            const double den = std::log(p_max_table) - std::log(p_sat);
                            const double psi = (den > 0) ? clamp01(num / den) : 0.0;
                            seed = rhoL + psi * (rhoM - rhoL);
                        } else {
                            // Melt table available but inconsistent at this query — fall through.
                            const double bulk_mod_scale = rhoL * ctx.R * T * kCompExtrapFactor;
                            seed = (bulk_mod_scale > 0) ? rhoL * (1.0 + (p - p_sat) / bulk_mod_scale) : rhoL;
                        }
                    } else {
                        // No melt table data — crude compressibility extrapolation.
                        // Conservative bulk-modulus-equivalent term (kCompExtrapFactor in
                        // the denominator damps the response for stiff liquids).
                        const double bulk_mod_scale = rhoL * ctx.R * T * kCompExtrapFactor;
                        seed = (bulk_mod_scale > 0) ? rhoL * (1.0 + (p - p_sat) / bulk_mod_scale) : rhoL;
                    }
                }
            } else {
                // Vapor branch (p < p_sat).
                double rhoV = NAN;
                try {
                    rhoV = sa.eval_sat(T, 'D', 1);
                } catch (...) {
                    rhoV = NAN;
                }
                if (!finite_pos(rhoV)) {
                    seed = ideal_gas(T, p);
                } else if (T > kVaporIdealCutoff * ctx.Tmax || p < kVaporLowPCutoff * p_sat) {
                    seed = ideal_gas(T, p);
                } else {
                    const double beta = clamp01((p_sat - p) / p_sat);
                    const double rg = ideal_gas(T, p);
                    if (finite_pos(rg)) {
                        seed = (1.0 - beta) * rhoV + beta * rg;
                    } else {
                        seed = rhoV;
                    }
                }
            }
        }
    }

    // Final safety net: NEVER return non-finite or non-positive.  Fall back
    // to ideal-gas; if that's also bad, return a tiny positive (1e-3 mol/m^3
    // — vacuum-ish but guaranteed positive).
    if (!finite_pos(seed)) {
        seed = ideal_gas(T, p);
    }
    if (!finite_pos(seed)) {
        seed = 1.0e-3;
    }
    return seed;
}

}  // namespace pxprecond
}  // namespace CoolProp
