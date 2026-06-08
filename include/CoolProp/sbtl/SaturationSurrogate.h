#ifndef COOLPROP_SBTL_SATURATION_SURROGATE_H
#define COOLPROP_SBTL_SATURATION_SURROGATE_H

#include <cstddef>
#include <memory>

#include "CoolProp/AbstractState.h"

namespace CoolProp {
namespace region {
class CubicSplineCurve;
}  // namespace region

namespace sbtl {

// SaturationSurrogate — a small cubic-spline surrogate for the
// saturation dome, built by sampling QT_INPUTS on a source backend at
// construction.  Mirrors the eval_sat / get_T_from_p shape of
// `superancillary::SuperAncillary` so the SVDSBTL backend can use it
// as a drop-in replacement when the source backend doesn't expose a
// SuperAncillary (REFPROP today; SRK / PR / etc. in the future).
//
// **Why this exists.**  Every dome-touching query through
// SVDSBTL&<source-without-SA> falls through to source.update(PQ_INPUTS
// / QT_INPUTS) for sat endpoints — REFPROP's PQ flash is ~1-5 ms per
// call, whereas an SA-style spline eval is ~50 ns.  Closing this 4-5
// orders-of-magnitude gap is the point of the surrogate.
//
// **What is sampled.**  At each of `n_knots` temperatures spanning
// [T_triple * (1 + 1e-6), 0.9999 * T_crit], the source is queried via
// QT_INPUTS at Q=0 and Q=1, capturing P, D, H, S, U on both sides.
// The additive 1e-6 offset on T_triple keeps the first probe one ULP-
// class step inside the source's flash range — some backends report
// T_triple at the boundary of where QT_INPUTS converges, and a bare
// T_triple probe fails on those.
// The 0.9999 multiplier on T_crit keeps the highest knot well clear of
// the sqrt-singularity at the critical point — cubic splines can't
// represent the singularity, so the surrogate trades accuracy near
// T_crit for clean interpolation everywhere else.
//
// **Spacing.**  Chebyshev-Lobatto knot distribution in T over
// [T_triple, 0.9999 * T_crit] — clusters knots toward both endpoints
// to absorb the sharp density bend approaching T_crit (saturated-rho
// is the budget-binding property; log-T spacing missed 1e-6 there).
// Cubic spline interpolation between the knots.
//
// **Accuracy target.**  Relative error <= 1e-6 vs source QT_INPUTS on
// smooth-fluid probes well away from the critical point.  Accuracy
// degrades as T -> T_crit (sqrt-singularity not captured by
// polynomials); callers querying that slice should use the
// critical-patch source fallback in SVDSBTLBackend instead.
//
// **Units.**  Matches SuperAncillary's molar convention: P in Pa, D in
// mol/m^3, H/U in J/mol, S in J/(mol K).  Caller (SVDSBTLBackend)
// converts to mass basis as needed downstream.
class SaturationSurrogate
{
   public:
    // Build a surrogate by sampling `src` at `n_knots` temperatures.
    // Returns nullptr on build failure (degenerate triple/critical pair,
    // n_knots < 4, or source rejects too many QT_INPUTS probes to form a
    // valid spline).  Mutates `src` state during sampling.
    //
    // Source-backend requirements: must support QT_INPUTS, must report
    // T_triple() and T_critical() correctly.  HEOS, REFPROP, IF97 all
    // satisfy these; mixtures do not (no single sat dome).
    [[nodiscard]] static std::unique_ptr<SaturationSurrogate> build_from_source(::CoolProp::AbstractState& src, std::size_t n_knots = 96);

    SaturationSurrogate(const SaturationSurrogate&) = delete;
    SaturationSurrogate& operator=(const SaturationSurrogate&) = delete;
    SaturationSurrogate(SaturationSurrogate&&) = delete;
    SaturationSurrogate& operator=(SaturationSurrogate&&) = delete;
    ~SaturationSurrogate();

    // Inverse pressure -> temperature.  Same contract as
    // SuperAncillary::get_T_from_p — argument is pressure in Pa (NOT its
    // log), returns saturation temperature in K.  Throws std::out_of_range
    // if `p` is outside the surrogate's [p_min, p_max] build range.
    [[nodiscard]] double get_T_from_p(double p) const;

    // Evaluate a saturated property at (T, side).  Matches
    // SuperAncillary::eval_sat — `what` is in {'P','D','H','S','U'},
    // `side` is 0 (liquid) or 1 (vapor).  Throws std::invalid_argument
    // on bad keys / sides, std::out_of_range on T outside the build range.
    [[nodiscard]] double eval_sat(double T, char what, int side) const;

    // Build interval, exposed for diagnostics + tests.
    [[nodiscard]] double T_min() const noexcept;
    [[nodiscard]] double T_max() const noexcept;
    [[nodiscard]] std::size_t n_knots() const noexcept;

   private:
    SaturationSurrogate();
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SATURATION_SURROGATE_H
