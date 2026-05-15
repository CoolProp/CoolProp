#ifndef COOLPROP_REGION_AXIS_TRANSFORM_H
#define COOLPROP_REGION_AXIS_TRANSFORM_H

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace CoolProp {
namespace region {

// Map a physical variable a in [a_lo, a_hi] to xi in [0, 1].  The map is
// either linear or applied to log(a); other monotone choices can be added
// later if needed.
enum class AxisScale : std::uint8_t
{
    LINEAR,
    LOG
};

struct AxisTransform
{
    AxisScale scale;
    double a_lo;        // physical lower bound
    double a_hi;        // physical upper bound
    double a_lo_t;      // transformed lower bound (a_lo or log(a_lo))
    double a_hi_t;      // transformed upper bound (a_hi or log(a_hi))
    double inv_span_t;  // 1 / (a_hi_t - a_lo_t); precomputed for hot path

    // Factory.  Throws on invalid inputs (a_hi <= a_lo, or LOG with non-
    // positive bound).
    static AxisTransform make(AxisScale s, double a_lo, double a_hi) {
        if (!(a_hi > a_lo)) {
            throw std::invalid_argument("AxisTransform: a_hi must exceed a_lo");
        }
        if (s == AxisScale::LOG && !(a_lo > 0.0)) {
            throw std::invalid_argument("AxisTransform: LOG requires a_lo > 0");
        }
        const double lo_t = (s == AxisScale::LOG) ? std::log(a_lo) : a_lo;
        const double hi_t = (s == AxisScale::LOG) ? std::log(a_hi) : a_hi;
        return AxisTransform{s, a_lo, a_hi, lo_t, hi_t, 1.0 / (hi_t - lo_t)};
    }

    // a -> xi in [0, 1].  Out-of-range inputs map outside [0,1] (no
    // clamp); callers can check before / after.
    [[nodiscard]] inline double forward(double a) const noexcept {
        const double a_t = (scale == AxisScale::LOG) ? std::log(a) : a;
        return (a_t - a_lo_t) * inv_span_t;
    }

    // xi in [0, 1] -> a.
    [[nodiscard]] inline double inverse(double xi) const noexcept {
        const double a_t = a_lo_t + xi * (a_hi_t - a_lo_t);
        return (scale == AxisScale::LOG) ? std::exp(a_t) : a_t;
    }

    // dxi/da at the physical point a (analytic).
    //
    //   LINEAR : dxi/da = 1 / (a_hi - a_lo) = inv_span_t
    //   LOG    : dxi/da = 1 / (a * (log a_hi - log a_lo))
    [[nodiscard]] inline double dxi_da(double a) const noexcept {
        if (scale == AxisScale::LOG) {
            return inv_span_t / a;
        }
        return inv_span_t;
    }
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_AXIS_TRANSFORM_H
