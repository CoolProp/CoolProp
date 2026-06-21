#ifndef COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H
#define COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H

#include <memory>
#include <optional>
#include <vector>

#include "CoolProp/AbstractState.h"
#include "CoolProp/region/BoundaryCurve.h"
#include "CoolProp/region/CubicSplineCurve.h"

namespace CoolProp {
namespace sbtl {

// Thermodynamics-aware factories that turn HEOS state queries into
// reusable BoundaryCurve objects.  These are convenience wrappers over
// the Phase 1 generic curves — Region/Atlas know nothing about HEOS,
// but every SBTL preset needs HEOS-sampled sat curves, so the wrapper
// lives here so the Phase 1 components stay fluid-agnostic.
//
// All factories return CubicSplineCurve (C^2 continuous) sampled at
// `n_knots` points log-uniform in pressure.  C^2 is the right default
// for SBTL: smooth across knots so the η normalization doesn't
// introduce derivative kinks (which would show up as bands of
// elevated SVD error — Phase 2a verified this empirically for Water).
//
// Default knot count (64) is plenty for the smooth sat curves of
// pure fluids; raise it if a fluid's curve has unusually fast log-p
// variation.
//
// TODO (Phase 2c): the *true* sat boundaries can be a no-refit view
// onto the existing fluid superancillary (h_sat(T) chained with
// T_sat(p) — both already machine-precision Chebyshev expansions on
// the SuperAncillary class).  Doing this needs (a) a public accessor
// on HelmholtzEOSMixtureBackend (currently get_superanc() is
// protected, see HelmholtzEOSMixtureBackend.h:65), and (b) a new
// `SuperancillaryBoundaryCurve` BoundaryCurve subclass that composes
// the two expansions and carries them through eval / eval_da.
// Same factory entry points; same return type; same SVDSurface
// consumer.  Defer to Phase 2c when the backend itself will need to
// pull the superancillary out anyway.

struct SatBoundaryBuildOptions
{
    std::size_t n_knots = 64;  // CubicSplineCurve knots per curve

    // T-floor walk-up: number of 0.5 K steps to try when T_min falls
    // below T_melt(p) at high p.  Default 1000 = 500 K of slack, which
    // covers water's melting line all the way to its EOS p_max ~1 GPa
    // (melting curve reaches ~370 K vs T_triple 273 K).  Subcritical
    // fluids typically clear in a handful of steps; this is just the
    // safety net for SUPER-region builds.
    std::size_t t_floor_walk_steps = 1000;
};

// h_sat,L(p) — mass enthalpy on the saturated-liquid side of the dome.
//
// When the source backend exposes a SuperAncillary (HEOS pure
// fluids), returns a SuperancillaryBoundaryCurve composing
// eval_sat('H', 0) with get_T_from_p — machine-precision residual
// against the source HEOS, no refit (CoolProp-8vg).  Otherwise
// (REFPROP / IF97 / pseudo-pure fluids), falls back to a 64-knot
// CubicSplineCurve sampled via PQ_INPUTS with Q = 0.
std::unique_ptr<region::BoundaryCurve> build_h_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts = {});

// h_sat,V(p) — mass enthalpy on the saturated-vapor side of the dome.
// Same SuperAncillary-or-spline dispatch as build_h_sat_L; PQ_INPUTS
// with Q = 1 in the spline fallback.
std::unique_ptr<region::BoundaryCurve> build_h_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts = {});

// s_sat,L(p) — mass entropy on the saturated-liquid side of the dome.
// Entropy analog of build_h_sat_L: SuperAncillary prop_key 'S', Q = 0,
// output_scale = 1/M (J/mol/K -> J/kg/K).  Spline fallback samples via
// PQ_INPUTS with Q = 0 and reads smass().
std::unique_ptr<region::BoundaryCurve> build_s_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts = {});

// s_sat,V(p) — mass entropy on the saturated-vapor side of the dome.
// Same dispatch as build_s_sat_L; PQ_INPUTS with Q = 1 in the fallback.
std::unique_ptr<region::BoundaryCurve> build_s_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                     const SatBoundaryBuildOptions& opts = {});

// T_sat(p) — saturation temperature.  Single-valued for pure fluids
// (LIQUID and VAPOR sat curves coincide on T).  Same dispatch — the
// SuperAncillary path uses its own inverse-p Chebyshev expansion.
std::unique_ptr<region::BoundaryCurve> build_T_sat(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                   const SatBoundaryBuildOptions& opts = {});

// rho_sat,L(T) — mass density on the saturated-liquid side of the
// dome, parameterized on temperature.  Used by the DT-indexed preset
// where the secondary axis (D) is bounded below by rho_sat,V(T) and
// above by rho_sat,L(T) at fixed T.
//
// When the source backend exposes a SuperAncillary, returns a
// SuperancillaryTemperatureBoundaryCurve calling eval_sat('D', 0, T)
// directly (no inversion, machine-precision residual).  Otherwise
// falls back to a CubicSplineCurve sampled via QT_INPUTS with Q=0
// at `n_knots` linear-uniform T values across [T_min, T_max].
//
// NOTE: for fluids with a density anomaly (water at T_anom ≈ 277 K,
// heavy water at ≈ 284 K) the caller MUST split the LIQUID region so
// each sub-region's [T_min, T_max] stays inside one monotonic piece of
// rho_sat,L(T).  Querying piece boundaries on the SuperAncillary via
// `get_approx1d('D', 0).get_x_at_extrema()` is the supported way.
std::unique_ptr<region::BoundaryCurve> build_rho_sat_L(::CoolProp::AbstractState& heos, double T_min, double T_max,
                                                       const SatBoundaryBuildOptions& opts = {});

// rho_sat,V(T) — mass density on the saturated-vapor side of the
// dome, parameterized on temperature.  Same dispatch as
// build_rho_sat_L; QT_INPUTS with Q=1 in the spline fallback.  No
// anomaly handling needed — rho_sat,V is monotone in T for every
// known fluid.
std::unique_ptr<region::BoundaryCurve> build_rho_sat_V(::CoolProp::AbstractState& heos, double T_min, double T_max,
                                                       const SatBoundaryBuildOptions& opts = {});

// Isobar h-floor: h on the cold-isotherm boundary.  For fluids with
// steep melting curves (Methane / Propane / CO2 LIQUID at high p),
// T_min may fall below T_melt(p) at part of the pressure range;
// the lambda walks T up in 0.5 K increments until HEOS accepts the
// state.  Resulting floor is non-isothermal but still monotone in p.
// Matches Phase 2a's `h_lo_liq` lambda.
std::unique_ptr<region::CubicSplineCurve> build_h_isotherm_floor(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_min,
                                                                 const SatBoundaryBuildOptions& opts = {});

// Isobar h-ceiling: h on the hot-isotherm boundary (T = T_max - margin
// so we stay strictly inside the HEOS validity envelope).
std::unique_ptr<region::CubicSplineCurve> build_h_isotherm_ceiling(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_max,
                                                                   const SatBoundaryBuildOptions& opts = {});

// Isobar s-floor / s-ceiling: entropy analogs of build_h_isotherm_floor
// / build_h_isotherm_ceiling.  Same cold-isotherm melting-line T-walk
// in the floor variant; both read smass() instead of hmass().
std::unique_ptr<region::CubicSplineCurve> build_s_isotherm_floor(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_min,
                                                                 const SatBoundaryBuildOptions& opts = {});

std::unique_ptr<region::CubicSplineCurve> build_s_isotherm_ceiling(::CoolProp::AbstractState& heos, double p_min, double p_max, double T_max,
                                                                   const SatBoundaryBuildOptions& opts = {});

// Locate the interior extrema of rho_sat,L(T) inside [T_min, T_max] —
// i.e., the T values where drho_sat,L/dT = 0.  For most fluids the
// returned vector is empty (rho_sat,L is monotone decreasing in T from
// triple to critical).  Water and heavy water each have one extremum
// (the density anomaly at ~277 K and ~284 K respectively); any future
// fluid with N extrema returns N entries.
//
// Used by the DT-indexed SVDSBTL preset (CoolProp-i7j) to split the
// LIQUID region into monotonic sub-regions — required because a
// non-monotone rho_sat,L(T) creates a non-simply-connected LIQUID
// region in (D, T) space (cells straddling the anomaly contain a
// discontinuity the SVD can't represent).
//
// Dispatch:
//   1. Source has SuperAncillary: query
//      `get_approx1d('D', 0).get_x_at_extrema()` — exact, zero compute.
//   2. Else: walk QT_INPUTS on a coarse T grid (~64 points), find any
//      sign change in dρ_sat,L/dT via central difference, bisect each
//      bracket with TOMS748 to locate the extremum.
//
// Returned T values are sorted ascending and lie strictly inside
// (T_min, T_max).
std::vector<double> find_rho_satL_extrema_T(::CoolProp::AbstractState& heos, double T_min, double T_max);

// Convenience: subcritical pressure range for `fluid`.  Returns
// (p_min, p_max) = (p_triple, p_crit * 0.999), where p_triple is the
// fluid's true triple-point pressure (heos.p_triple(), with a
// QT-at-Ttriple fallback when that is not finite/positive).  The lowest
// tabulated isobar therefore sits exactly at p_triple so a PT query at
// p == p_triple() resolves (RegionAtlas containment is inclusive);
// PQ(p_triple, 0) converges, so build_T_sat can sample the bottom knot.
// Driven by the fluid's HEOS AbstractState — needs an instance already
// constructed by the caller.
//
// `p_min_override` (absolute Pa, from the backend's `pmin` option)
// replaces the p_triple floor when supplied.  It MUST be >= p_triple:
// below the triple line the liquid-vapour saturation boundary that
// bounds the subcritical regions does not exist, so a sub-triple floor
// is rejected with ValueError rather than failing obscurely deep in the
// boundary-curve sampling.
std::pair<double, double> subcritical_pressure_range(::CoolProp::AbstractState& heos, std::optional<double> p_min_override = std::nullopt);

// Convenience: supercritical pressure range for `fluid`.  Returns
// (p_min, p_max) ≈ (p_crit * 1.001, pmax_eos * 0.99) — a thin margin
// above the critical point so the SUPER region doesn't share its
// lower boundary with the subcritical regions' upper bound (the SVD
// gets ill-conditioned right at p_crit; the gap is small enough that
// caller-side dispatch falls back to direct HEOS in the gap, when /
// if that fallback ships).  Driven by the fluid's HEOS AbstractState.
std::pair<double, double> supercritical_pressure_range(::CoolProp::AbstractState& heos);

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H
