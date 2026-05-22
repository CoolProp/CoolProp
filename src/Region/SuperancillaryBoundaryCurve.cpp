#include "CoolProp/region/SuperancillaryBoundaryCurve.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace CoolProp {
namespace region {

std::unique_ptr<SuperancillaryBoundaryCurve> SuperancillaryBoundaryCurve::build(std::shared_ptr<SuperAncillary_t> sa, double p_min, double p_max,
                                                                                char prop_key, short Q, double output_scale) {
    if (!sa) {
        throw std::invalid_argument("SuperancillaryBoundaryCurve::build: null SuperAncillary handle");
    }
    if (!(p_min > 0.0) || !(p_max > p_min)) {
        throw std::invalid_argument("SuperancillaryBoundaryCurve::build: need 0 < p_min < p_max");
    }
    if (Q != 0 && Q != 1) {
        throw std::invalid_argument("SuperancillaryBoundaryCurve::build: Q must be 0 (liquid) or 1 (vapor)");
    }

    // Probe the curve at log-spaced sample points to pre-compute the
    // (b_min, b_max) AABB.  Cheap (Chebyshev eval ~ns each); 32 probes
    // give a tight bound for any smooth boundary.
    constexpr std::size_t kProbes = 32;
    double b_min = std::numeric_limits<double>::infinity();
    double b_max = -std::numeric_limits<double>::infinity();
    const double log_p_min = std::log(p_min);
    const double log_p_max = std::log(p_max);
    for (std::size_t k = 0; k < kProbes; ++k) {
        const double f = static_cast<double>(k) / static_cast<double>(kProbes - 1);
        const double p = std::exp(log_p_min + f * (log_p_max - log_p_min));
        try {
            const double T_sat = sa->get_T_from_p(p);
            const double y = sa->eval_sat(T_sat, prop_key, Q) * output_scale;
            if (std::isfinite(y)) {
                b_min = std::min(b_min, y);
                b_max = std::max(b_max, y);
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // SA may reject queries near the critical point or
            // exotically close to the triple point.  Skip and keep
            // the surrounding probes.
        }
    }
    if (!std::isfinite(b_min) || !std::isfinite(b_max) || !(b_min < b_max)) {
        throw std::runtime_error("SuperancillaryBoundaryCurve::build: SuperAncillary rejected all probes in the requested p range");
    }
    return std::make_unique<SuperancillaryBoundaryCurve>(std::move(sa), p_min, p_max, prop_key, Q, output_scale, b_min, b_max);
}

}  // namespace region
}  // namespace CoolProp
