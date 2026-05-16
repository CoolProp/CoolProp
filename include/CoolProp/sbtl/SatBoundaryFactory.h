#ifndef COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H
#define COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H

#include <memory>

#include "AbstractState.h"
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
};

// h_sat,L(p) — mass enthalpy on the saturated-liquid side of the dome.
//
// HEOS calls go through PQ_INPUTS with Q = 0.  Throws if p is outside
// the fluid's saturation range or HEOS rejects the PQ flash.
std::unique_ptr<region::CubicSplineCurve> build_h_sat_L(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                        const SatBoundaryBuildOptions& opts = {});

// h_sat,V(p) — mass enthalpy on the saturated-vapor side of the dome.
// PQ_INPUTS with Q = 1.
std::unique_ptr<region::CubicSplineCurve> build_h_sat_V(::CoolProp::AbstractState& heos, double p_min, double p_max,
                                                        const SatBoundaryBuildOptions& opts = {});

// T_sat(p) — saturation temperature.  Single-valued for pure fluids
// (LIQUID and VAPOR sat curves coincide on T).  Uses PQ_INPUTS with
// Q = 0.
std::unique_ptr<region::CubicSplineCurve> build_T_sat(::CoolProp::AbstractState& heos, double p_min, double p_max,
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

// Convenience: subcritical pressure range for `fluid`.  Returns
// (p_min, p_max) ≈ (p_triple * 1.01, p_crit * 0.999) so PQ flashes
// don't fail at the exact boundary.  Driven by the fluid's HEOS
// AbstractState — needs an instance already constructed by the caller.
std::pair<double, double> subcritical_pressure_range(::CoolProp::AbstractState& heos);

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SAT_BOUNDARY_FACTORY_H
