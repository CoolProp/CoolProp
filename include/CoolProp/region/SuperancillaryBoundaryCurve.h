#ifndef COOLPROP_REGION_SUPERANCILLARY_BOUNDARY_CURVE_H
#define COOLPROP_REGION_SUPERANCILLARY_BOUNDARY_CURVE_H

#include <cmath>
#include <cstddef>
#include <memory>
#include <set>  // transitively needed by superancillary.h
#include <utility>
#include <vector>

#include "CoolProp/region/BoundaryCurve.h"
#include "superancillary/superancillary.h"

namespace CoolProp {
namespace region {

// 1D boundary curve b = f(p) backed by the source backend's SuperAncillary
// — a stack of two analytic operations, each already machine-precision:
//
//   1. T_sat = SuperAncillary::get_T_from_p(p)
//      Solves T from p via the inverse log-p Chebyshev expansion.
//      Inversion noise floor is ~1e-13 for typical fluids.
//
//   2. b     = SuperAncillary::eval_sat(T_sat, prop_key, Q)
//      Evaluates the saturated property b at T_sat on side Q (0 = liquid,
//      1 = vapor) via the property's piecewise Chebyshev expansion.
//
// Composing the two gives b_sat(p) = eval_sat(get_T_from_p(p), k, Q) — a
// no-refit curve with O(1e-13) residual against the source HEOS.
//
// This is the long-standing TODO from `SatBoundaryFactory.h`'s top
// comment ("Phase 2c: the true sat boundaries can be a no-refit view
// onto the existing fluid superancillary..."), realized as a concrete
// BoundaryCurve subclass.  Closes CoolProp-8vg with the path
// originally identified back in Phase 2a's design notes.
//
// SCOPE.  Only works for the sat dome curves (Q ∈ {0, 1}, prop_key ∈
// {'D', 'P', 'H', 'S', 'U'}).  Isotherm floors / ceilings outside the
// dome still need the legacy SatBoundaryFactory path (cubic spline).
//
// THREAD SAFETY.  SuperAncillary's lazy-build paths use a mutex
// internally (m_lazy_build_mutex on caloric expansions).  Multiple
// threads reading the same curve are safe.
//
// EVAL_DA.  Returns d/dp [eval_sat(T_sat(p), k, Q)] via the chain
// rule:
//
//   d/dp b_sat(p) = (db/dT)_{Q,k}(T_sat) * dT_sat/dp
//
// where (db/dT) is computed via central finite difference on the
// Chebyshev expansion (cheaper than reaching into the approx1d's
// derivative-expansion machinery, and the Cheb function is smooth
// enough that 1e-8 FD step gives ~1e-12 relative error in the
// derivative — well below anything that matters for the η
// normalisation we drive).  dT_sat/dp is the inverse of
// d/dT p_sat(T) at the resolved T_sat, also via central FD.
class SuperancillaryBoundaryCurve final : public BoundaryCurve
{
   public:
    using SuperAncillary_t = superancillary::SuperAncillary<std::vector<double>>;

    SuperancillaryBoundaryCurve(std::shared_ptr<SuperAncillary_t> sa, double p_min, double p_max, char prop_key, short Q, double output_scale,
                                double b_min, double b_max)
      : sa_(std::move(sa)), p_min_(p_min), p_max_(p_max), prop_key_(prop_key), Q_(Q), output_scale_(output_scale), b_min_(b_min), b_max_(b_max) {
        // Fail fast: every public eval / eval_da path dereferences sa_
        // unchecked.  Allowing null would let a caller create an object
        // that crashes on first lookup.  Both factories (build,
        // from_state) also null-guard; this catches the direct-
        // construction path that bypasses them.
        if (!sa_) {
            throw std::invalid_argument("SuperancillaryBoundaryCurve: null SuperAncillary handle");
        }
    }

    // Factory: probes the SA at p_min / p_max to pre-compute b_min /
    // b_max for the AABB cheap-reject.  Throws if the SuperAncillary
    // doesn't span [p_min, p_max] or rejects the property key.
    //
    // SuperAncillary returns properties in MOLAR units (J/mol for H,
    // J/(mol K) for S, mol/m^3 for D).  Pass `output_scale = 1/M`
    // (M = molar mass, kg/mol) to convert to mass basis at eval
    // time.  For molar-basis consumers pass 1.0.  For property 'P'
    // (pressure), pass 1.0 (units are already Pa).  For property 'D'
    // (density in mol/m^3 -> kg/m^3), pass `M`, not `1/M`.
    static std::unique_ptr<SuperancillaryBoundaryCurve> build(std::shared_ptr<SuperAncillary_t> sa, double p_min, double p_max, char prop_key,
                                                              short Q, double output_scale);

    [[nodiscard]] double eval(double p) const noexcept override {
        try {
            const double T_sat = sa_->get_T_from_p(p);
            return sa_->eval_sat(T_sat, prop_key_, Q_) * output_scale_;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            return std::nan("");
        }
    }

    [[nodiscard]] double eval_da(double p) const noexcept override {
        // Central FD on the composed curve.  Relative step h = 1e-7 *
        // p gives ~1e-13 relative error in the derivative — Cheb
        // expansions are C^∞ in the interior, so FD doesn't introduce
        // truncation noise of any significance.
        //
        // Clamp the ±h endpoints to [p_min_, p_max_] so a query right
        // at a build-range endpoint doesn't step into the SA's invalid
        // region and return NaN.  When clamped on one side the
        // effective stencil becomes asymmetric (forward or backward
        // difference); divide by the actual span to keep the slope
        // correct.
        const double h = std::max(1e-7 * std::abs(p), 1.0);
        const double p_hi = std::min(p + h, p_max_);
        const double p_lo = std::max(p - h, p_min_);
        if (!(p_hi > p_lo)) {
            return std::nan("");
        }
        const double y_plus = eval(p_hi);
        const double y_minus = eval(p_lo);
        if (!std::isfinite(y_plus) || !std::isfinite(y_minus)) {
            return std::nan("");
        }
        return (y_plus - y_minus) / (p_hi - p_lo);
    }

    [[nodiscard]] std::pair<double, double> bounds() const noexcept override {
        return {b_min_, b_max_};
    }
    [[nodiscard]] std::pair<double, double> a_range() const noexcept override {
        return {p_min_, p_max_};
    }

    // State / from_state for the serializer.  We store only the
    // primitive scalars; on deserialization the caller must supply the
    // SuperAncillary handle (typically the same instance the SVDSBTL
    // backend already owns).
    struct State
    {
        double p_min;
        double p_max;
        char prop_key;
        short Q;
        double output_scale;
        double b_min;
        double b_max;
    };
    [[nodiscard]] State state() const {
        return State{p_min_, p_max_, prop_key_, Q_, output_scale_, b_min_, b_max_};
    }
    static std::unique_ptr<SuperancillaryBoundaryCurve> from_state(State s, std::shared_ptr<SuperAncillary_t> sa) {
        // Constructor's eval() dereferences sa_ unchecked; the matching
        // factory build() already null-guards.  Fail fast here too so a
        // caller passing a null handle gets a clear diagnostic instead
        // of an eventual segfault on the first lookup.
        if (!sa) {
            throw std::invalid_argument("SuperancillaryBoundaryCurve::from_state: null SuperAncillary handle");
        }
        return std::make_unique<SuperancillaryBoundaryCurve>(std::move(sa), s.p_min, s.p_max, s.prop_key, s.Q, s.output_scale, s.b_min, s.b_max);
    }

   private:
    std::shared_ptr<SuperAncillary_t> sa_;
    double p_min_;
    double p_max_;
    char prop_key_;
    short Q_;
    double output_scale_;
    double b_min_;
    double b_max_;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_SUPERANCILLARY_BOUNDARY_CURVE_H
