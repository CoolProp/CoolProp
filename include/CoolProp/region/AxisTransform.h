#ifndef COOLPROP_REGION_AXIS_TRANSFORM_H
#define COOLPROP_REGION_AXIS_TRANSFORM_H

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace CoolProp {
namespace region {

// Map a physical variable a in [a_lo, a_hi] to xi in [0, 1].  The map is
// linear, log(a), or critical-scaling-aware cube-root (CoolProp-4u9).
// All three are monotone in (a_lo, a_hi); other monotone choices can be
// added later if needed.
//
// POWER: forward xi = 1 - cbrt((a_hi - a) / (a_hi - a_lo)); inverse
// a = a_hi - (a_hi - a_lo) * (1 - xi)^3.  Grid lines uniform in xi
// crowd toward a_hi (intended use: NC sub-critical sub-region in p
// with a_hi = (1 - eps)·pc).  The cube-root absorbs the (pc - p)^(1/δ)
// scaling singularity at the critical point — Ising β ≈ 0.326 ≈ 1/3 —
// so the Hermite-cubic SVD interpolant stays well-conditioned all the
// way to within ~1 ppm of pc.  Diagnostic measurements on R245fa /
// Water / CO2 LIQUID + VAPOR at NT=200 NR=800 Cheb rank=20 show max
// off-grid error 1e-9 to 1e-7 for POWER vs 1e-4 to 1e-3 for LOG over
// the same band.
//
// POWER_LO: mirror of POWER that crowds toward a_lo instead of a_hi.
// Forward xi = cbrt((a - a_lo) / (a_hi - a_lo)); inverse a = a_lo +
// (a_hi - a_lo) * xi^3.  Intended use: NC super-critical sub-region
// with a_lo = (1 + eps)·pc (pc is the LOW end of the band there).
// Same accuracy story as POWER, opposite anchor.
enum class AxisScale : std::uint8_t
{
    LINEAR,
    LOG,
    POWER,
    POWER_LO
};

struct AxisTransform
{
    AxisScale scale;
    double a_lo;        // physical lower bound
    double a_hi;        // physical upper bound
    double a_lo_t;      // transformed lower bound (a_lo or log(a_lo)); unused for POWER
    double a_hi_t;      // transformed upper bound (a_hi or log(a_hi)); unused for POWER
    double inv_span_t;  // 1 / (a_hi_t - a_lo_t); precomputed for hot path.
                        // For POWER, this stores 1 / (a_hi - a_lo) — the linear
                        // span — so chain-rule derivatives have a single scale
                        // factor to multiply.

    // Factory.  Throws on invalid inputs (a_hi <= a_lo, or LOG with non-
    // positive bound).
    static AxisTransform make(AxisScale s, double a_lo, double a_hi) {
        if (!(a_hi > a_lo)) {
            throw std::invalid_argument("AxisTransform: a_hi must exceed a_lo");
        }
        if (s == AxisScale::LOG && !(a_lo > 0.0)) {
            throw std::invalid_argument("AxisTransform: LOG requires a_lo > 0");
        }
        if (s == AxisScale::POWER || s == AxisScale::POWER_LO) {
            // a_lo_t / a_hi_t intentionally unused; inv_span_t stores
            // 1 / (a_hi - a_lo) for the chain-rule derivative path.
            return AxisTransform{s, a_lo, a_hi, 0.0, 0.0, 1.0 / (a_hi - a_lo)};
        }
        const double lo_t = (s == AxisScale::LOG) ? std::log(a_lo) : a_lo;
        const double hi_t = (s == AxisScale::LOG) ? std::log(a_hi) : a_hi;
        return AxisTransform{s, a_lo, a_hi, lo_t, hi_t, 1.0 / (hi_t - lo_t)};
    }

    // a -> xi in [0, 1].  Out-of-range inputs map outside [0,1] (no
    // clamp); callers can check before / after.
    [[nodiscard]] inline double forward(double a) const noexcept {
        if (scale == AxisScale::POWER) {
            return 1.0 - std::cbrt((a_hi - a) * inv_span_t);
        }
        if (scale == AxisScale::POWER_LO) {
            return std::cbrt((a - a_lo) * inv_span_t);
        }
        const double a_t = (scale == AxisScale::LOG) ? std::log(a) : a;
        return (a_t - a_lo_t) * inv_span_t;
    }

    // xi in [0, 1] -> a.
    [[nodiscard]] inline double inverse(double xi) const noexcept {
        if (scale == AxisScale::POWER) {
            const double t = 1.0 - xi;
            return a_hi - (a_hi - a_lo) * t * t * t;
        }
        if (scale == AxisScale::POWER_LO) {
            return a_lo + (a_hi - a_lo) * xi * xi * xi;
        }
        const double a_t = a_lo_t + xi * (a_hi_t - a_lo_t);
        return (scale == AxisScale::LOG) ? std::exp(a_t) : a_t;
    }

    // dxi/da at the physical point a (analytic).
    //
    //   LINEAR   : dxi/da = 1 / (a_hi - a_lo) = inv_span_t
    //   LOG      : dxi/da = 1 / (a * (log a_hi - log a_lo))
    //   POWER    : xi = 1 - ((a_hi - a)/(a_hi - a_lo))^(1/3)
    //                dxi/da = (1/3) * inv_span_t * ((a_hi - a) * inv_span_t)^(-2/3)
    //                Diverges as a -> a_hi (anchor side — inverse has zero slope there).
    //   POWER_LO : xi = ((a - a_lo)/(a_hi - a_lo))^(1/3)
    //                dxi/da = (1/3) * inv_span_t * ((a - a_lo) * inv_span_t)^(-2/3)
    //                Diverges as a -> a_lo (anchor side).
    [[nodiscard]] inline double dxi_da(double a) const noexcept {
        if (scale == AxisScale::LOG) {
            return inv_span_t / a;
        }
        if (scale == AxisScale::POWER) {
            const double u = (a_hi - a) * inv_span_t;
            return (1.0 / 3.0) * inv_span_t / std::cbrt(u * u);
        }
        if (scale == AxisScale::POWER_LO) {
            const double u = (a - a_lo) * inv_span_t;
            return (1.0 / 3.0) * inv_span_t / std::cbrt(u * u);
        }
        return inv_span_t;
    }
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_AXIS_TRANSFORM_H
