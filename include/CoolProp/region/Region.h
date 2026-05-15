#ifndef COOLPROP_REGION_REGION_H
#define COOLPROP_REGION_REGION_H

#include <memory>
#include <utility>

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/BoundaryCurve.h"

namespace CoolProp {
namespace region {

// A single region in (a, b) space.
//
// Geometry:
//   - Primary axis a is mapped to xi in [0, 1] via an AxisTransform
//     (linear or log scale).  The transform fixes the [a_lo, a_hi]
//     interval.
//   - Secondary axis b is bounded below by `b_lo(a)` and above by
//     `b_hi(a)`, two BoundaryCurves that are explicit functions of a.
//     The normalised coordinate eta = (b - b_lo(a)) / (b_hi(a) - b_lo(a))
//     lives in [0, 1] when (a, b) is inside the region.
//
// Dispatch:
//   `aabb_contains` is the cheap first-pass filter — four comparisons
//   against the precomputed bounding box; no virtual calls.  It will
//   admit false positives for curved regions (points inside the AABB
//   but outside the curve envelope), so a second-pass `curve_contains`
//   evaluates the boundary curves to make the precise decision.  This
//   matches the typical region-atlas dispatch: scan AABBs first, then
//   call into the curve check only on hits.
//
// Lifetime:
//   The Region owns its two BoundaryCurves via std::unique_ptr; copy
//   construction is deleted because BoundaryCurves are abstract.  Move
//   construction is allowed.
class Region
{
   public:
    struct BBox
    {
        double a_lo;
        double a_hi;
        double b_min;
        double b_max;
    };

    Region(AxisTransform primary, std::unique_ptr<BoundaryCurve> b_lo, std::unique_ptr<BoundaryCurve> b_hi);

    Region(const Region&) = delete;
    Region& operator=(const Region&) = delete;
    Region(Region&&) = default;
    Region& operator=(Region&&) = default;
    ~Region() = default;

    // O(1), no curve evaluation.
    [[nodiscard]] inline bool aabb_contains(double a, double b) const noexcept {
        return a >= bbox_.a_lo && a <= bbox_.a_hi && b >= bbox_.b_min && b <= bbox_.b_max;
    }

    // Evaluates both BoundaryCurves at a; only call when aabb_contains
    // already returned true.  Returns false outside the curve envelope.
    [[nodiscard]] bool curve_contains(double a, double b) const noexcept;

    // (a, b) -> (xi, eta).  Assumes (a, b) is inside the region; caller
    // is responsible for an aabb_contains / curve_contains gate.
    [[nodiscard]] std::pair<double, double> to_normalized(double a, double b) const noexcept;

    // (xi, eta) -> (a, b).
    [[nodiscard]] std::pair<double, double> from_normalized(double xi, double eta) const noexcept;

    [[nodiscard]] const AxisTransform& primary() const noexcept {
        return primary_;
    }
    [[nodiscard]] const BoundaryCurve& b_lo() const noexcept {
        return *b_lo_;
    }
    [[nodiscard]] const BoundaryCurve& b_hi() const noexcept {
        return *b_hi_;
    }
    [[nodiscard]] const BBox& bbox() const noexcept {
        return bbox_;
    }

   private:
    AxisTransform primary_;
    std::unique_ptr<BoundaryCurve> b_lo_;
    std::unique_ptr<BoundaryCurve> b_hi_;
    BBox bbox_;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_REGION_H
