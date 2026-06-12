#include "CoolProp/sbtl/SatBoundaryFactory.h"

#include <algorithm>  // std::sort / std::min / std::max in find_rho_satL_extrema_T
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "boost/math/tools/toms748_solve.hpp"

#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/region/SuperancillaryBoundaryCurve.h"
#include "CoolProp/region/SuperancillaryTemperatureBoundaryCurve.h"
#include "CoolProp/DataStructures.h"

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

std::unique_ptr<region::BoundaryCurve> build_s_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/true)) {
        try {
            // SuperAncillary returns molar entropy (J/mol/K); SBTL
            // stores mass-basis (J/kg/K).  s_mass = s_mol / M.
            const double inv_M = 1.0 / heos.molar_mass();
            return region::SuperancillaryBoundaryCurve::build(sa, p_min, p_max, 'S', 0, inv_M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Fall through to the spline path on SA range mismatch.
        }
    }
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 0.0);
        return heos.smass();
    });
}

std::unique_ptr<region::BoundaryCurve> build_s_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/true)) {
        try {
            const double inv_M = 1.0 / heos.molar_mass();
            return region::SuperancillaryBoundaryCurve::build(sa, p_min, p_max, 'S', 1, inv_M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    }
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 1.0);
        return heos.smass();
    });
}

namespace {
// T-parameterized cubic-spline sampler.  Mirror of
// spline_through_log_p_samples, but linear-T spacing since the T span
// over a dome arc is at most factor ~3 and log spacing would over-
// concentrate knots at low T where rho_sat curves are nearly linear.
std::unique_ptr<region::CubicSplineCurve> spline_through_linear_T_samples(double T_min, double T_max, std::size_t n_knots,
                                                                          const std::function<double(double)>& f) {
    if (!(T_min > 0.0) || !(T_max > T_min) || n_knots < 2) {
        throw std::invalid_argument("SatBoundaryFactory: invalid T range or n_knots (need 0 < T_min < T_max and n_knots >= 2)");
    }
    std::vector<double> T_knots(n_knots);
    std::vector<double> y(n_knots);
    const double dT = (T_max - T_min) / static_cast<double>(n_knots - 1);
    for (std::size_t k = 0; k < n_knots; ++k) {
        T_knots[k] = T_min + static_cast<double>(k) * dT;
        y[k] = f(T_knots[k]);
    }
    return region::CubicSplineCurve::build(std::move(T_knots), std::move(y));
}
}  // namespace

std::unique_ptr<region::BoundaryCurve> build_rho_sat_L(::CoolProp::AbstractState& heos, double T_min, double T_max,
                                                       const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/false)) {
        try {
            // SuperAncillary returns molar density (mol/m^3); SBTL
            // stores mass-basis (kg/m^3).  rho_mass = rho_mol * M.
            const double M = heos.molar_mass();
            return region::SuperancillaryTemperatureBoundaryCurve::build(sa, T_min, T_max, 'D', 0, M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Fall through to the spline path on SA range mismatch.
        }
    }
    return spline_through_linear_T_samples(T_min, T_max, opts.n_knots, [&](double T) {
        heos.update(::CoolProp::QT_INPUTS, 0.0, T);
        return heos.rhomass();
    });
}

std::unique_ptr<region::BoundaryCurve> build_rho_sat_V(::CoolProp::AbstractState& heos, double T_min, double T_max,
                                                       const SatBoundaryBuildOptions& opts) {
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/false)) {
        try {
            const double M = heos.molar_mass();
            return region::SuperancillaryTemperatureBoundaryCurve::build(sa, T_min, T_max, 'D', 1, M);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    }
    return spline_through_linear_T_samples(T_min, T_max, opts.n_knots, [&](double T) {
        heos.update(::CoolProp::QT_INPUTS, 1.0, T);
        return heos.rhomass();
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
            } catch (...) {
                // Below melting line at this p; bump T up and retry.  Reset
                // the state first: a failed update can leave a poisoned
                // guess that makes the next (T_try, p) probe on this same
                // reused state fail spuriously — the same retry-on-melting-
                // reject hazard as build_rho_max_envelope (CoolProp-i9t8).
                heos.clear();
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

std::unique_ptr<region::CubicSplineCurve> build_s_isotherm_floor(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_min,
                                                                 const SatBoundaryBuildOptions& opts) {
    // Entropy analog of build_h_isotherm_floor — same melting-line
    // T-walk + heos.clear() reset on reject (CoolProp-i9t8), reads
    // smass() instead of hmass().
    const std::size_t walk_steps = opts.t_floor_walk_steps;
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&, walk_steps](double p) {
        for (std::size_t k = 0; k < walk_steps; ++k) {
            const double T_try = T_min + 0.5 * static_cast<double>(k);
            try {
                heos.update(::CoolProp::PT_INPUTS, p, T_try);
                return heos.smass();
            } catch (...) {
                heos.clear();
            }
        }
        throw std::runtime_error("build_s_isotherm_floor: floor unreachable within " + std::to_string(walk_steps / 2) + " K of T_min");
    });
}

std::unique_ptr<region::CubicSplineCurve> build_s_isotherm_ceiling(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_max,
                                                                   const SatBoundaryBuildOptions& opts) {
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PT_INPUTS, p, T_max);
        return heos.smass();
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

std::vector<double> find_rho_satL_extrema_T(::CoolProp::AbstractState& heos, double T_min, double T_max) {
    if (!(T_max > T_min)) {
        return {};
    }
    // SuperAncillary path — pieces are constructed to be monotonic in
    // the dependent variable, so get_x_at_extrema() returns precisely
    // the T values where drho_sat,L/dT = 0.
    if (auto sa = try_get_superanc(heos, /*need_caloric=*/false)) {
        try {
            const auto& approx = sa->get_approx1d('D', /*Q=*/0);
            std::vector<double> hits;
            for (double T : approx.get_x_at_extrema()) {
                if (T > T_min && T < T_max) {
                    hits.push_back(T);
                }
            }
            std::sort(hits.begin(), hits.end());
            return hits;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // SA may reject the property key (no rhoL expansion).
            // Fall through to the QT_INPUTS walk.
        }
    }
    // Fallback: walk QT_INPUTS at 64 T-knots, finite-difference
    // drho_sat,L/dT, bisect each sign change.  Cheap (one-shot at
    // preset construction).
    constexpr std::size_t kProbes = 64;
    std::vector<double> T_grid(kProbes);
    std::vector<double> rhoL(kProbes);
    const double dT_grid = (T_max - T_min) / static_cast<double>(kProbes - 1);
    for (std::size_t k = 0; k < kProbes; ++k) {
        T_grid[k] = T_min + static_cast<double>(k) * dT_grid;
        try {
            heos.update(::CoolProp::QT_INPUTS, 0.0, T_grid[k]);
            rhoL[k] = heos.rhomass();
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            rhoL[k] = std::nan("");
        }
    }
    // Slope sign at each interior knot (central FD).  Skip non-finite
    // entries (rare — only at very-near-critical T) by treating them
    // as "no information".
    auto slope_sign = [&](std::size_t i) -> int {
        if (i == 0 || i + 1 >= kProbes) return 0;
        const double rl = rhoL[i - 1];
        const double rr = rhoL[i + 1];
        if (!std::isfinite(rl) || !std::isfinite(rr)) return 0;
        const double dr = rr - rl;
        return (dr > 0.0) ? 1 : (dr < 0.0 ? -1 : 0);
    };
    std::vector<double> hits;
    for (std::size_t i = 1; i + 1 < kProbes; ++i) {
        const int s_lo = slope_sign(i);
        const int s_hi = slope_sign(i + 1);
        if (s_lo == 0 || s_hi == 0 || s_lo == s_hi) continue;
        // Sign change between knot i and i+1; bisect on the slope sign.
        try {
            auto residual = [&](double T) -> double {
                // Central FD of rho_sat,L(T) at T; reuse the kProbes
                // grid spacing as the step (well above the noise floor
                // of QT_INPUTS).
                const double h = 0.5 * dT_grid;
                heos.update(::CoolProp::QT_INPUTS, 0.0, T - h);
                const double rl = heos.rhomass();
                heos.update(::CoolProp::QT_INPUTS, 0.0, T + h);
                const double rr = heos.rhomass();
                return rr - rl;
            };
            boost::math::tools::eps_tolerance<double> tol(48);
            std::uintmax_t max_iter = 64;
            auto root = boost::math::tools::toms748_solve(residual, T_grid[i], T_grid[i + 1], tol, max_iter);
            hits.push_back(0.5 * (root.first + root.second));
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Bisection failed (rare; could indicate the bracketed
            // sign change is numerical noise rather than a real
            // extremum).  Skip.
        }
    }
    std::sort(hits.begin(), hits.end());
    // Fallback bracketing can return near-duplicate roots under FD
    // noise; collapse them so downstream LIQUID splitting never gets a
    // zero-width sub-region.
    hits.erase(std::unique(hits.begin(), hits.end(),
                           [](double x, double y) {
                               const double tol = 1e-8 * std::max({1.0, std::abs(x), std::abs(y)});
                               return std::abs(x - y) <= tol;
                           }),
               hits.end());
    return hits;
}

}  // namespace sbtl
}  // namespace CoolProp
