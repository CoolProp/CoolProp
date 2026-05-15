#ifndef COOLPROP_REGION_BOUNDARY_CURVE_H
#define COOLPROP_REGION_BOUNDARY_CURVE_H

#include <utility>

namespace CoolProp {
namespace region {

// 1D curve b = f(a) over a closed interval [a_lo, a_hi].  The primary
// axis variable a is the independent input; the secondary axis variable
// b is what the curve produces.
//
// Used as a building block for Region — two BoundaryCurves (lower and
// upper) plus an AxisTransform define a curved 2D region in (a, b) space
// that maps to the (xi, eta) unit square.
//
// All concrete subclasses are required to provide *analytic* eval_da
// (never a finite-difference shim).  Each subclass differentiates its
// own internal representation symbolically.
//
// Lifetime: instances are owned by Region via std::unique_ptr; the
// interface lives behind a virtual call but is only invoked inside the
// AABB-hit branch of region dispatch, so the indirection is off the
// inner SVD-eval hot path.
class BoundaryCurve
{
   public:
    BoundaryCurve() = default;
    BoundaryCurve(const BoundaryCurve&) = default;
    BoundaryCurve(BoundaryCurve&&) = default;
    BoundaryCurve& operator=(const BoundaryCurve&) = default;
    BoundaryCurve& operator=(BoundaryCurve&&) = default;
    virtual ~BoundaryCurve() = default;

    // Value of the curve at a.  Behaviour outside [a_lo, a_hi] is
    // implementation-defined (typically extrapolated by the underlying
    // basis); callers should normally not query outside the build range.
    [[nodiscard]] virtual double eval(double a) const noexcept = 0;

    // First derivative db/da at a — analytic, not finite-difference.
    [[nodiscard]] virtual double eval_da(double a) const noexcept = 0;

    // Min and max of b achieved on the build interval [a_lo, a_hi].
    // Precomputed at build; used to populate the parent Region's
    // axis-aligned bounding box for cheap region-dispatch filtering.
    [[nodiscard]] virtual std::pair<double, double> bounds() const noexcept = 0;

    // Build interval [a_lo, a_hi] on the primary axis.
    [[nodiscard]] virtual std::pair<double, double> a_range() const noexcept = 0;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_BOUNDARY_CURVE_H
