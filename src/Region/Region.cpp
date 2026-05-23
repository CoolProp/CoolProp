#include "CoolProp/region/Region.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace CoolProp {
namespace region {

namespace {

Region::BBox compute_bbox(const AxisTransform& primary, const BoundaryCurve& b_lo, const BoundaryCurve& b_hi) noexcept {
    const auto lo_bounds = b_lo.bounds();
    const auto hi_bounds = b_hi.bounds();
    return {primary.a_lo, primary.a_hi, lo_bounds.first, hi_bounds.second};
}

}  // namespace

// Sentinel BBox that aabb_contains will never report a hit on (a_hi
// strictly less than a_lo).  Used for moved-from instances so they
// fail the AABB filter without anyone having to add an "is_alive"
// check to the hot path.
constexpr Region::BBox kEmptyBBox = {1.0, -1.0, 1.0, -1.0};

Region::Region(AxisTransform primary, std::unique_ptr<BoundaryCurve> b_lo, std::unique_ptr<BoundaryCurve> b_hi)
  : primary_(primary), b_lo_(std::move(b_lo)), b_hi_(std::move(b_hi)), bbox_((b_lo_ && b_hi_) ? compute_bbox(primary_, *b_lo_, *b_hi_) : kEmptyBBox) {
    if (!b_lo_ || !b_hi_) {
        throw std::invalid_argument("Region: b_lo and b_hi BoundaryCurves must be non-null");
    }
}

Region::Region(Region&& other) noexcept : primary_(other.primary_), b_lo_(std::move(other.b_lo_)), b_hi_(std::move(other.b_hi_)), bbox_(other.bbox_) {
    other.bbox_ = kEmptyBBox;
}

Region& Region::operator=(Region&& other) noexcept {
    if (this != &other) {
        primary_ = other.primary_;
        b_lo_ = std::move(other.b_lo_);
        b_hi_ = std::move(other.b_hi_);
        bbox_ = other.bbox_;
        other.bbox_ = kEmptyBBox;
    }
    return *this;
}

bool Region::curve_contains(double a, double b) const noexcept {
    // curve_contains is sign-only — it just decides whether b is
    // bracketed by the two curve values at a, not the precise value
    // of either curve.  Route through eval_fast so SA-backed curves
    // can use their precomputed linear-interp surrogate (~5 ns vs
    // ~80 ns for the full SA composition).  Precision-critical
    // callers (to_normalized / from_normalized) keep using eval().
    const double b_lo_val = b_lo_->eval_fast(a);
    const double b_hi_val = b_hi_->eval_fast(a);
    return b >= b_lo_val && b <= b_hi_val;
}

std::pair<double, double> Region::to_normalized(double a, double b) const noexcept {
    const double xi = primary_.forward(a);
    const double b_lo_val = b_lo_->eval(a);
    const double b_hi_val = b_hi_->eval(a);
    const double span = b_hi_val - b_lo_val;
    // Guard against a degenerate boundary span: if the two curves
    // converge at this a (e.g. a critical-point pinch), eta is
    // undefined.  Return NaN so callers see an obvious failure rather
    // than a silently-wrong eta = 0 that propagates into an SVD lookup
    // and produces garbage with no signal.
    const double tol = std::numeric_limits<double>::epsilon() * (1.0 + std::abs(b_lo_val) + std::abs(b_hi_val));
    const double eta = (std::abs(span) <= tol) ? std::numeric_limits<double>::quiet_NaN() : ((b - b_lo_val) / span);
    return {xi, eta};
}

std::pair<double, double> Region::from_normalized(double xi, double eta) const noexcept {
    const double a = primary_.inverse(xi);
    const double b_lo_val = b_lo_->eval(a);
    const double b_hi_val = b_hi_->eval(a);
    const double b = b_lo_val + eta * (b_hi_val - b_lo_val);
    return {a, b};
}

}  // namespace region
}  // namespace CoolProp
