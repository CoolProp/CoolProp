#include "CoolProp/sbtl/SatBoundaryFactory.h"

#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "DataStructures.h"

namespace CoolProp {
namespace sbtl {

namespace {

// Sample `f(p)` at `n_knots` log-uniform points in [p_min, p_max] and
// build a CubicSplineCurve through the resulting (p, f(p)) knots.
// This is the common skeleton every factory in this file shares;
// extracting it once keeps the per-factory code 4–5 lines.
std::unique_ptr<region::CubicSplineCurve> spline_through_log_p_samples(double p_min, double p_max, std::size_t n_knots,
                                                                       const std::function<double(double)>& f) {
    if (!(p_max > p_min) || n_knots < 2) {
        throw std::invalid_argument("SatBoundaryFactory: invalid p range or n_knots");
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

std::unique_ptr<region::CubicSplineCurve> build_h_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                        const SatBoundaryBuildOptions& opts) {
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 0.0);
        return heos.hmass();
    });
}

std::unique_ptr<region::CubicSplineCurve> build_h_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                        const SatBoundaryBuildOptions& opts) {
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        heos.update(::CoolProp::PQ_INPUTS, p, 1.0);
        return heos.hmass();
    });
}

std::unique_ptr<region::CubicSplineCurve> build_T_sat(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                      const SatBoundaryBuildOptions& opts) {
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
    return spline_through_log_p_samples(p_min, p_max, opts.n_knots, [&](double p) {
        for (int k = 0; k < 40; ++k) {
            const double T_try = T_min + 0.5 * static_cast<double>(k);
            try {
                heos.update(::CoolProp::PT_INPUTS, p, T_try);
                return heos.hmass();
            } catch (...) {  // NOLINT(bugprone-empty-catch)
                // Below melting line at this p; bump T up and retry.
            }
        }
        throw std::runtime_error("build_h_isotherm_floor: floor unreachable within 20 K of T_min");
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

}  // namespace sbtl
}  // namespace CoolProp
