#ifndef COOLPROP_SBTL_SVD_SURFACE_FACTORY_H
#define COOLPROP_SBTL_SVD_SURFACE_FACTORY_H

#include <cstddef>
#include <cstdint>
#include <optional>

#include "CoolProp/AbstractState.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SurfaceSpec.h"
#include "CoolProp/svd/SVDDecomposition.h"

namespace CoolProp {
namespace sbtl {

struct BuildOptions
{
    svd::SlopeSource slope_source = svd::SlopeSource::NATURAL_CUBIC_SPLINE;
    bool verbose = false;
};

// Generic factory: walk the per-region (xnorm, log_p) grid, call HEOS
// per cell via the spec's `update_state` + `read_property` callbacks,
// build one SVD per (region, property), seal into an SVDSurface.
//
// SurfaceSpec is consumed (moved-from) so the unique_ptr boundary
// curves it owns transfer cleanly into the new SVDSurface's regions.
//
// Throws on:
//   - empty spec.regions / spec.properties
//   - mismatched spec.input_pair vs the callbacks' assumed inputs
//     (only detectable if the callbacks throw at sampling time)
//   - SVDBuilder failures (rank-deficient, non-finite samples surviving
//     the row-fill, ...)
SVDSurface build_surface(::CoolProp::AbstractState& heos, SurfaceSpec spec, const BuildOptions& opts = {});

// Preset SurfaceSpec builders for the two MVP input pairs.  Each one
// constructs:
//   - 2-region atlas (LIQUID, VAPOR) over [p_triple, p_crit]
//   - CubicSpline sat boundaries via SatBoundaryFactory
//   - The standard property list for that input pair
//   - Default NT=200, NR=800, rank=20 grid sizes
//
// Caller can then mutate the SurfaceSpec before calling
// build_surface() (raise the rank for a hard fluid, swap a boundary
// curve, add an extra output property, etc.).
//
// Subsequent phases ship DU / HS / PS / DT presets the same way — no
// changes to build_surface() or downstream consumers.
namespace presets {

// Build-time knobs shared by the subcritical presets and the backend.
// Bundling them in one struct means a future knob is a new field here
// (plus where it is parsed and consumed) — NOT a new positional
// parameter rippling through every preset signature, the dispatcher,
// and every call site.
struct PresetOptions
{
    std::size_t NT = 200;    // points along the secondary (non-log) axis
    std::size_t NR = 800;    // points along the primary (log-p) axis
    std::int32_t rank = 20;  // SVD truncation rank
    // Absolute-Pa lower-pressure bound for the PT / HmassP / PSmass
    // surfaces (the backend's `pmin` option).  nullopt -> default
    // p_triple floor; must be >= p_triple (see
    // subcritical_pressure_range).  Ignored by the DmassT preset, which
    // is temperature-indexed and has no pressure floor.
    std::optional<double> p_min;
};

// HmassP_INPUTS preset.  (a, b) = (p, h).  Output properties: rho, T, s, u.
SurfaceSpec ph_subcritical(::CoolProp::AbstractState& heos, const PresetOptions& opts = {});

// PT_INPUTS preset.  (a, b) = (p, T).  Output properties: rho, h, s, u.
SurfaceSpec pt_subcritical(::CoolProp::AbstractState& heos, const PresetOptions& opts = {});

// PSmass_INPUTS preset.  (a, b) = (p, s).  Output properties: rho, T, h, u.
// Entropy analog of ph_subcritical: s is the secondary axis + query
// input; the same region geometry (LIQUID/VAPOR/NC/SUPER, IF97 R2/R3/R5
// split) applies with the build_s_* boundary curves.
SurfaceSpec ps_subcritical(::CoolProp::AbstractState& heos, const PresetOptions& opts = {});

// DmassT_INPUTS preset.  (a, b) = (T, D).  Output properties: p, h, s, u.
//
// Primary advantage over PT/PH: (D, T) is the Helmholtz EOS's native
// coordinate, so every output property is a direct evaluation — no
// inversion, no critical-region stiffness, no "cells with -760 MPa"
// failure modes like BICUBIC's invert_single_phase_y on #1301.
//
// Geometry mirrors ph_subcritical: LIQUID + VAPOR sub-critical (T <
// T_critical * 0.999), SUPER super-critical (T > T_critical * 1.001).
// Secondary axis D bounded by rho_sat,L(T) / rho_sat,V(T) and
// per-fluid HEOS-validity D-extents on the non-dome side.
//
// Anomaly handling: for fluids whose ρ_sat,L(T) curve has interior
// extrema (water at T_anom ≈ 277 K, heavy water at ≈ 284 K), the
// LIQUID region is split into N+1 sub-regions where N is the number
// of extrema — each sub-region's T span lies inside one monotonic
// piece of ρ_sat,L(T).  Auto-detected via the source backend's
// SuperAncillary `get_x_at_extrema()` on the rho_sat,L expansion;
// gracefully falls back to a single LIQUID region when no SA is
// available (REFPROP / IF97 sources, or fluids without a SA in
// HEOS).
SurfaceSpec dt_subcritical(::CoolProp::AbstractState& heos, const PresetOptions& opts = {});

}  // namespace presets

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_FACTORY_H
