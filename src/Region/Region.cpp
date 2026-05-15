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

Region::Region(AxisTransform primary, std::unique_ptr<BoundaryCurve> b_lo, std::unique_ptr<BoundaryCurve> b_hi)
  : primary_(primary),
    b_lo_(std::move(b_lo)),
    b_hi_(std::move(b_hi)),
    bbox_((b_lo_ && b_hi_) ? compute_bbox(primary_, *b_lo_, *b_hi_) : Region::BBox{0, 0, 0, 0}) {
    if (!b_lo_ || !b_hi_) {
        throw std::invalid_argument("Region: b_lo and b_hi BoundaryCurves must be non-null");
    }
}

bool Region::curve_contains(double a, double b) const noexcept {
    const double b_lo_val = b_lo_->eval(a);
    const double b_hi_val = b_hi_->eval(a);
    return b >= b_lo_val && b <= b_hi_val;
}

std::pair<double, double> Region::to_normalized(double a, double b) const noexcept {
    const double xi = primary_.forward(a);
    const double b_lo_val = b_lo_->eval(a);
    const double b_hi_val = b_hi_->eval(a);
    const double span = b_hi_val - b_lo_val;
    // Guard against a degenerate boundary span: if the two curves
    // converge at this a (e.g. a critical-point pinch), eta is
    // undefined.  Return 0 rather than inf/NaN so callers using the
    // result for grid indexing don't blow up.
    const double tol = std::numeric_limits<double>::epsilon() * (1.0 + std::abs(b_lo_val) + std::abs(b_hi_val));
    const double eta = (std::abs(span) <= tol) ? 0.0 : ((b - b_lo_val) / span);
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
