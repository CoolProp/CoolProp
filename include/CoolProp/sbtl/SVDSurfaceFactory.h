#ifndef COOLPROP_SBTL_SVD_SURFACE_FACTORY_H
#define COOLPROP_SBTL_SVD_SURFACE_FACTORY_H

#include <cstddef>

#include "AbstractState.h"
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

// HmassP_INPUTS preset.  (a, b) = (p, h).  Output properties: rho, T, s, u.
SurfaceSpec ph_subcritical(::CoolProp::AbstractState& heos, std::size_t NT = 200, std::size_t NR = 800, std::int32_t rank = 20);

// PT_INPUTS preset.  (a, b) = (p, T).  Output properties: rho, h, s, u.
SurfaceSpec pt_subcritical(::CoolProp::AbstractState& heos, std::size_t NT = 200, std::size_t NR = 800, std::int32_t rank = 20);

}  // namespace presets

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_FACTORY_H
