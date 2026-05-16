#ifndef COOLPROP_REGION_CUBIC_SPLINE_CURVE_H
#define COOLPROP_REGION_CUBIC_SPLINE_CURVE_H

#include <memory>
#include <utility>
#include <vector>

#include "CoolProp/region/BoundaryCurve.h"

namespace CoolProp {
namespace region {

// Natural cubic spline interpolant through n supplied (a_i, b_i) knots,
// presented as a BoundaryCurve.
//
// Construction solves a tridiagonal system for the second derivatives
// M_i = p''(a_i) under the natural boundary condition M_0 = M_{n-1} = 0.
// Evaluation is the standard segment-cubic formula in terms of the
// (b_i, M_i) tables, which makes the analytic derivative immediate.
//
// All concrete subclasses of BoundaryCurve provide *analytic* eval_da;
// here the segment derivative is the piecewise quadratic obtained by
// differentiating the spline polynomial term-by-term — no finite
// differences at runtime.
//
// SMOOTHNESS CAVEAT.  Natural cubic spline assumes the source function
// is smooth.  When the knot sequence captures a kink or near-cusp in
// the underlying b(a), the global-C² requirement forces oscillation
// in the neighbourhood of the kink (Gibbs-like overshoot).  Phase 2a's
// e2e validation surfaced this for Propane: at the pressure where
// the LIQUID-floor lambda crosses the melting line and bumps T up,
// h_lo(p) has a slope kink, and the resulting spline shows a thin
// vertical band of elevated error on the LIQUID side at ~250-400
// kJ/kg.  Workarounds when smoothness can't be guaranteed: build a
// PiecewiseChebyshevCurve over a uniform partition (which absorbs
// the kink at one piece boundary but trades the global-C² guarantee
// for a single C⁰ point), or split the data into two splines.
class CubicSplineCurve final : public BoundaryCurve
{
   public:
    // Build from sorted knots.  a values must be strictly increasing.
    // n >= 2 required (n == 2 degenerates to the line between the two
    // points with M = 0).  Throws std::invalid_argument on malformed
    // input.
    static std::unique_ptr<CubicSplineCurve> build(std::vector<double> a, std::vector<double> b);

    // Plain-data snapshot for serialization.  Captures the natural-
    // spline state: knot positions a, knot values b, the
    // second-derivative table M, and the tight bounds (b_min, b_max)
    // computed by build() via cubic-extremum root finding.
    struct State
    {
        std::vector<double> a;
        std::vector<double> b;
        std::vector<double> M;
        double b_min;
        double b_max;
    };
    [[nodiscard]] State state() const;

    // Trusted factory: reconstructs a CubicSplineCurve from a State
    // snapshot WITHOUT re-running the tridiagonal solve.  Used by the
    // SBTL deserializer.  Validates the state shape (a/b/M all same
    // length, a strictly increasing) and rejects malformed input.
    static std::unique_ptr<CubicSplineCurve> from_state(State s);

    [[nodiscard]] double eval(double a) const noexcept override;
    [[nodiscard]] double eval_da(double a) const noexcept override;
    [[nodiscard]] std::pair<double, double> bounds() const noexcept override;
    [[nodiscard]] std::pair<double, double> a_range() const noexcept override;

   private:
    CubicSplineCurve(std::vector<double> a, std::vector<double> b, std::vector<double> M, double b_min, double b_max);

    // Locate the segment index i such that a_[i] <= a < a_[i+1], clamped
    // to [0, n-2].  Used by eval / eval_da.
    [[nodiscard]] std::size_t locate(double a) const noexcept;

    std::vector<double> a_;  // knot positions, size n
    std::vector<double> b_;  // knot values,    size n
    std::vector<double> M_;  // second derivatives at knots, size n
    double b_min_;           // tight min of the spline on [a_[0], a_[n-1]]
    double b_max_;           // tight max of the spline on [a_[0], a_[n-1]]
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_CUBIC_SPLINE_CURVE_H
