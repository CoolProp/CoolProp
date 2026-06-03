#ifndef COOLPROP_REGION_SUPERANCILLARY_TEMPERATURE_BOUNDARY_CURVE_H
#define COOLPROP_REGION_SUPERANCILLARY_TEMPERATURE_BOUNDARY_CURVE_H

#include <cmath>
#include <cstddef>
#include <memory>
#include <set>  // transitively needed by superancillary.h
#include <utility>
#include <vector>

#include "CoolProp/region/BoundaryCurve.h"
#include "CoolProp/superancillary/superancillary.h"

namespace CoolProp {
namespace region {

// 1D boundary curve b = f(T) backed by the source backend's SuperAncillary.
// T-parameterized sibling of SuperancillaryBoundaryCurve — that one
// composes get_T_from_p with eval_sat for p-axis presets (PH, PT); this
// one calls eval_sat(T, prop, Q) directly because the DT-indexed preset
// uses T as the primary axis and needs no inversion.
//
// SCOPE.  Sat dome curves only (Q in {0, 1}, prop_key in {'D','P','H',
// 'S','U'}).  Same scope as the pressure-parameterized sibling.
//
// EVAL_DA.  Central FD on eval(T) with relative step h = max(1e-7*T, 1)
// K — well inside the smooth interior of every property's Cheb
// expansion.  Endpoints clamped to [T_min, T_max] so a query at the
// build edge doesn't step out into the SA's invalid region.
//
// EVAL_FAST.  1024-point linear-T table (no log scaling — T span is at
// most factor ~3 from triple to critical).  ~8 KB per curve.
class SuperancillaryTemperatureBoundaryCurve final : public BoundaryCurve
{
   public:
    using SuperAncillary_t = superancillary::SuperAncillary<std::vector<double>>;

    SuperancillaryTemperatureBoundaryCurve(std::shared_ptr<SuperAncillary_t> sa, double T_min, double T_max, char prop_key, short Q,
                                           double output_scale, double b_min, double b_max)
      : sa_(std::move(sa)), T_min_(T_min), T_max_(T_max), prop_key_(prop_key), Q_(Q), output_scale_(output_scale), b_min_(b_min), b_max_(b_max) {
        // Fail fast: every public eval / eval_da / eval_fast path
        // dereferences sa_ unchecked.  Both factory build() and
        // from_state() null-guard; the direct-construction path also
        // needs the guard.
        if (!sa_) {
            throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve: null SuperAncillary handle");
        }
        if (!(T_min_ > 0.0) || !(T_max_ > T_min_)) {
            throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve: need 0 < T_min < T_max");
        }
        // Validate the deserialization-reachable invariants here, not
        // only in build(): from_state() constructs directly, so a bad
        // Q / prop_key from a corrupt blob must fail fast rather than
        // silently produce a NaN-backed curve.  Mirrors the SCOPE doc.
        if (Q_ != 0 && Q_ != 1) {
            throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve: Q must be 0 (liquid) or 1 (vapor)");
        }
        switch (prop_key_) {
            case 'D':
            case 'P':
            case 'H':
            case 'S':
            case 'U':
                break;
            default:
                throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve: unsupported prop_key (expect one of D,P,H,S,U)");
        }
        build_surrogate_();
    }

    // Factory: probes the SA across [T_min, T_max] to pre-compute b_min
    // / b_max for the AABB cheap-reject.  Throws if the SuperAncillary
    // rejects the property key or the entire requested range.
    //
    // SuperAncillary returns properties in MOLAR units (mol/m^3 for D,
    // J/mol for H, J/(mol K) for S).  For mass-basis density consumers
    // pass `output_scale = M` (kg/mol); for molar-basis pass 1.0.
    static std::unique_ptr<SuperancillaryTemperatureBoundaryCurve> build(std::shared_ptr<SuperAncillary_t> sa, double T_min, double T_max,
                                                                         char prop_key, short Q, double output_scale);

    [[nodiscard]] double eval(double T) const noexcept override {
        try {
            return sa_->eval_sat(T, prop_key_, Q_) * output_scale_;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            return std::nan("");
        }
    }

    [[nodiscard]] double eval_da(double T) const noexcept override {
        // Central FD on eval(T).  Relative step h = max(1e-7*T, 1 K)
        // gives ~1e-13 relative error in the derivative — Cheb
        // expansions are C^infinity in the interior.  Clamp the
        // stencil to [T_min_, T_max_] so an endpoint query doesn't
        // step into the SA's invalid range.  Asymmetric stencils
        // still divide by the actual span, preserving slope accuracy.
        const double h = std::max(1e-7 * std::abs(T), 1.0);
        const double T_hi = std::min(T + h, T_max_);
        const double T_lo = std::max(T - h, T_min_);
        if (!(T_hi > T_lo)) {
            return std::nan("");
        }
        const double y_plus = eval(T_hi);
        const double y_minus = eval(T_lo);
        if (!std::isfinite(y_plus) || !std::isfinite(y_minus)) {
            return std::nan("");
        }
        return (y_plus - y_minus) / (T_hi - T_lo);
    }

    // Fast 1e-6-accurate surrogate via 1024-point linear-T table.  Used
    // by Region::curve_contains for atlas dispatch.  Out-of-range clamps
    // to the nearest endpoint.
    [[nodiscard]] double eval_fast(double T) const noexcept override {
        const double idx_f = (T - T_min_) * inv_T_step_;
        if (!(idx_f > 0.0)) {
            return surrogate_table_.front();
        }
        const auto last = static_cast<double>(kSurrogatePoints - 1);
        if (idx_f >= last) {
            return surrogate_table_.back();
        }
        const auto i = static_cast<std::size_t>(idx_f);
        const double frac = idx_f - static_cast<double>(i);
        const double y0 = surrogate_table_[i];
        const double y1 = surrogate_table_[i + 1];
        return y0 + frac * (y1 - y0);
    }

    [[nodiscard]] std::pair<double, double> bounds() const noexcept override {
        return {b_min_, b_max_};
    }
    [[nodiscard]] std::pair<double, double> a_range() const noexcept override {
        return {T_min_, T_max_};
    }

    // State / from_state for the serializer.  Only primitive scalars
    // are stored; deserialization re-attaches the live SuperAncillary
    // handle (typically the same instance the SVDSBTL backend owns).
    struct State
    {
        double T_min;
        double T_max;
        char prop_key;
        short Q;
        double output_scale;
        double b_min;
        double b_max;
    };
    [[nodiscard]] State state() const {
        return State{T_min_, T_max_, prop_key_, Q_, output_scale_, b_min_, b_max_};
    }
    static std::unique_ptr<SuperancillaryTemperatureBoundaryCurve> from_state(State s, std::shared_ptr<SuperAncillary_t> sa) {
        if (!sa) {
            throw std::invalid_argument("SuperancillaryTemperatureBoundaryCurve::from_state: null SuperAncillary handle");
        }
        return std::make_unique<SuperancillaryTemperatureBoundaryCurve>(std::move(sa), s.T_min, s.T_max, s.prop_key, s.Q, s.output_scale, s.b_min,
                                                                        s.b_max);
    }

   private:
    static constexpr std::size_t kSurrogatePoints = 1024;

    // Build the 1024-point linear-T surrogate table.  NOT noexcept:
    // surrogate_table_.resize can throw std::bad_alloc.
    void build_surrogate_() {
        surrogate_table_.resize(kSurrogatePoints);
        const double T_step = (T_max_ - T_min_) / static_cast<double>(kSurrogatePoints - 1);
        inv_T_step_ = 1.0 / T_step;
        bool any_finite = false;
        for (std::size_t k = 0; k < kSurrogatePoints; ++k) {
            const double T = T_min_ + static_cast<double>(k) * T_step;
            surrogate_table_[k] = eval(T);
            any_finite = any_finite || std::isfinite(surrogate_table_[k]);
        }
        // The backfill below assumes ≥1 finite sample exists.  build()'s
        // probe guarantees that, but from_state() can construct directly
        // with a handle whose valid range doesn't overlap [T_min, T_max]
        // — every sample non-finite.  Fail fast rather than returning an
        // all-NaN table that eval_fast() would silently interpolate.
        if (!any_finite) {
            throw std::runtime_error("SuperancillaryTemperatureBoundaryCurve: no finite samples in surrogate build range");
        }
        // Backfill any non-finite samples (rare: SA can reject a probe
        // exotically close to the critical point) from the nearest
        // finite neighbour.  Storing a raw NaN would let eval_fast()
        // interpolate through it and make Region::curve_contains
        // false-miss valid points even though the exact eval() path
        // still works.  The factory's b_min/b_max probe already
        // guarantees ≥1 finite sample exists, so the backward+forward
        // scan always resolves.
        for (std::size_t k = 0; k < kSurrogatePoints; ++k) {
            if (std::isfinite(surrogate_table_[k])) continue;
            double fill = std::nan("");
            for (std::size_t j = k; j-- > 0;) {
                if (std::isfinite(surrogate_table_[j])) {
                    fill = surrogate_table_[j];
                    break;
                }
            }
            if (!std::isfinite(fill)) {
                for (std::size_t j = k + 1; j < kSurrogatePoints; ++j) {
                    if (std::isfinite(surrogate_table_[j])) {
                        fill = surrogate_table_[j];
                        break;
                    }
                }
            }
            surrogate_table_[k] = fill;  // finite by the factory's probe guarantee
        }
    }

    std::shared_ptr<SuperAncillary_t> sa_;
    double T_min_;
    double T_max_;
    char prop_key_;
    short Q_;
    double output_scale_;
    double b_min_;
    double b_max_;
    std::vector<double> surrogate_table_;
    double inv_T_step_ = 0.0;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_SUPERANCILLARY_TEMPERATURE_BOUNDARY_CURVE_H
