#include "CoolProp/sbtl/SatBoundaryFactory.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/region/SuperancillaryBoundaryCurve.h"
#include "DataStructures.h"

namespace CoolProp {
namespace sbtl {

namespace {

// Try to pull a SuperAncillary handle off the source AbstractState.
// Returns null when the source isn't a HelmholtzEOSMixtureBackend
// (e.g. IF97, REFPROP) or when the fluid's HEOS model doesn't ship a
// SuperAncillary expansion (some pseudo-pure fluids / mixtures).
// Side effect: triggers caloric superancillary lazy-build so the
// downstream eval_sat('H', ...) calls don't pay that cost per-eval.
std::shared_ptr<region::SuperancillaryBoundaryCurve::SuperAncillary_t> try_get_superanc(::CoolProp::AbstractState& src, bool need_caloric) {
    auto* helmholtz = dynamic_cast<::CoolProp::HelmholtzEOSMixtureBackend*>(&src);
    if (helmholtz == nullptr) return nullptr;
    std::shared_ptr<region::SuperancillaryBoundaryCurve::SuperAncillary_t> sa;
    try {
        sa = helmholtz->get_superanc();
    } catch (...) {  // NOLINT(bugprone-empty-catch)
        // get_superanc throws for mixtures / pseudo-pures.  Treat as
        // "no SA available" and fall through to the spline path.
    }
    if (!sa) return nullptr;
    if (need_caloric) {
        try {
            helmholtz->ensure_caloric_superancillaries();
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Caloric build failed (rare: superanc rejection).  Keep
            // the SA handle for non-caloric properties but the caller
            // will null-check on the property key.
            return nullptr;
        }
    }
    return sa;
}

// Sample `f(p)` at `n_knots` log-uniform points in [p_min, p_max] and
// build a CubicSplineCurve through the resulting (p, f(p)) knots.
// This is the common skeleton every factory in this file shares;
// extracting it once keeps the per-factory code 4–5 lines.
std::unique_ptr<region::CubicSplineCurve> spline_through_log_p_samples(double p_min, double p_max, std::size_t n_knots,
                                                                       const std::function<double(double)>& f) {
    // p_min must be strictly positive for log-space sampling; otherwise
    // std::log(p_min) silently feeds -inf into the spline knots and we'd
    // cascade into hard-to-debug downstream failures.
    if (!(p_min > 0.0) || !(p_max > p_min) || n_knots < 2) {
        throw std::invalid_argument("SatBoundaryFactory: invalid p range or n_knots (need 0 < p_min < p_max and n_knots >= 2)");
    }
    std::vector<double> p_knots(n_knots);
    std::vector<double> y(n_knots);
    const double log_p_min = std::log(p_min);
    const double log_p_max = std::log(p_max);
    for (std::size_t k = 0; k < n_knots; ++k) {
        const double log_p = log_p_min + static_cast<double>(k) * (log_p_max - log_p_min) / static_cast<double>(n_knots - 1);
        p_knots[k] = std::exp(log_p);
        y[k] = f(p_knots[k]);
    }
    return region::CubicSplineCurve::build(std::move(p_knots), std::move(y));
}

}  // namespace

std::unique_ptr<region::BoundaryCurve> build_h_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/true)) {
        try {
            // SuperAncillary returns molar enthalpy (J/mol); SBTL
            // stores mass-basis (J/kg).  h_mass = h_mol / M.
            const double inv_M = 1.0 / heos.molar_mass();
            return region::SuperancillaryBoundaryCurve::build(sa, p_min, p_max, 'H', 0, inv_M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Fall through to the spline path on SA range mismatch.
        }
    }
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 0.0);
        return heos.hmass();
    });
}

std::unique_ptr<region::BoundaryCurve> build_h_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/true)) {
        try {
            const double inv_M = 1.0 / heos.molar_mass();
            return region::SuperancillaryBoundaryCurve::build(sa, p_min, p_max, 'H', 1, inv_M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    }
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 1.0);
        return heos.hmass();
    });
}

std::unique_ptr<region::BoundaryCurve> build_T_sat(::CoolProp::AbstractState& heos, double p_min, double p_max, const SatBoundaryBuildOptions& opts) {
    // No prop_key in the SuperAncillary for T directly — it's the
    // *result* of get_T_from_p, not a tabulated property.  Build a
    // small wrapper class would be possible but for now keep the
    // spline path here; T_sat is only used by IF97-source presets
    // (which don't have an SA anyway).
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 0.0);
        return heos.T();
    });
}

std::unique_ptr<region::CubicSplineCurve> build_h_isotherm_floor(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_min,
                                                                 const SatBoundaryBuildOptions& opts) {
    // For fluids with steep melting curves the requested T_min may be
    // below T_melt(p) at high p; walk T up in 0.5 K steps until HEOS
    // accepts the (T, p) state.  Matches Phase 2a's `h_lo_liq` lambda.
    const std::size_t walk_steps = opts.t_floor_walk_steps;
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&, walk_steps](double p) {
        for (std::size_t k = 0; k < walk_steps; ++k) {
            const double T_try = T_min + 0.5 * static_cast<double>(k);
            try {
                heos.update(::CoolProp::PT_INPUTS, p, T_try);
                return heos.hmass();
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // Below melting line at this p; bump T up and retry.
            }
        }
        throw std::runtime_error("build_h_isotherm_floor: floor unreachable within " + std::to_string(walk_steps / 2) + " K of T_min");
    });
}

std::unique_ptr<region::CubicSplineCurve> build_h_isotherm_ceiling(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_max,
                                                                   const SatBoundaryBuildOptions& opts) {
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PT_INPUTS, p, T_max);
        return heos.hmass();
    });
}

std::pair<double, double> subcritical_pressure_range(::CoolProp::AbstractState& heos) {
    const double p_crit = heos.p_critical();
    // Sample HEOS at triple T to get p_triple.  A tiny margin (1%)
    // helps PQ flashes converge at the boundary.
    heos.update(::CoolProp::QT_INPUTS, 0.0, heos.Ttriple() * 1.001);
    const double p_triple = heos.p();
    return {p_triple * 1.01, 0.999 * p_crit};
}

std::pair<double, double> supercritical_pressure_range(::CoolProp::AbstractState& heos) {
    const double p_crit = heos.p_critical();
    const double p_max_eos = heos.pmax();
    return {1.001 * p_crit, 0.99 * p_max_eos};
}

}  // namespace sbtl
}  // namespace CoolProp
