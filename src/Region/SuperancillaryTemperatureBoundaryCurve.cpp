#include "CoolProp/region/SuperancillaryTemperatureBoundaryCurve.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace CoolProp {
namespace region {

std::unique_ptr<SuperancillaryTemperatureBoundaryCurve> SuperancillaryTemperatureBoundaryCurve::build(std::shared_ptr<SuperAncillary_t> sa,
                                                                                                      double T_min, double T_max, char prop_key,
                                                                                                      short Q, double output_scale) {
    if (!sa) {
        throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve::build: null SuperAncillary handle");
    }
    if (!(T_min > 0.0) || !(T_max > T_min)) {
        throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve::build: need 0 < T_min < T_max");
    }
    if (Q != 0 && Q != 1) {
        throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve::build: Q must be 0 (liquid) or 1 (vapor)");
    }

    // Probe the curve at linear-T sample points to pre-compute the
    // (b_min, b_max) AABB.  32 probes give a tight bound for smooth
    // boundary curves.  T span is at most factor ~3, so log spacing
    // wouldn't help — linear is the right default.
    constexpr std::size_t kProbes = 32;
    double b_min = std::numeric_limits<double>::infinity();
    double b_max = -std::numeric_limits<double>::infinity();
    std::size_t finite_probe_count = 0;
    for (std::size_t k = 0; k < kProbes; ++k) {
        const double f = static_cast<double>(k) / static_cast<double>(kProbes - 1);
        const double T = T_min + f * (T_max - T_min);
        try {
            const double y = sa->eval_sat(T, prop_key, Q) * output_scale;
            if (std::isfinite(y)) {
                ++finite_probe_count;
                b_min = std::min(b_min, y);
                b_max = std::max(b_max, y);
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // SA may reject queries exotically close to the critical
            // point.  Skip and keep surrounding probes.
        }
    }
    if (finite_probe_count == 0) {
        throw std::runtime_error("SuperancillaryTemperatureBoundaryCurve::build: SuperAncillary rejected all probes in the requested T range");
    }
    if (!(b_max > b_min)) {  // CodeQL-safe form of b_min == b_max (b_min <= b_max always holds here)
        // Degenerate near-flat curve; inflate by 1 ULP per side so the
        // AABB has non-zero b-extent and Region machinery doesn't trip.
        const double eps = std::max(std::abs(b_min) * std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::min());
        b_min -= eps;
        b_max += eps;
    }
    return std::make_unique<SuperancillaryTemperatureBoundaryCurve>(std::move(sa), T_min, T_max, prop_key, Q, output_scale, b_min, b_max);
}

}  // namespace region
}  // namespace CoolProp
