#ifndef COOLPROP_SBTL_SURFACE_SPEC_H
#define COOLPROP_SBTL_SURFACE_SPEC_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CoolProp/AbstractState.h"
#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/BoundaryCurve.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "CoolProp/DataStructures.h"

namespace CoolProp {
namespace sbtl {

// Per-region geometry: primary axis + two boundary curves on the
// secondary axis.  Identical to what `region::Region` takes, but
// staged in a value-type aggregate so SurfaceSpec can be assembled
// incrementally before any Region exists.
struct RegionSpec
{
    // Default-state primary is filled in by builder code before any
    // build_surface() call reads it; the aggregate default is a
    // never-contains AxisTransform that fails AxisTransform::make's
    // invariants if accidentally used.
    region::AxisTransform primary{region::AxisScale::LINEAR, 0.0, 1.0, 0.0, 1.0, 1.0};
    std::unique_ptr<region::BoundaryCurve> b_lo;
    std::unique_ptr<region::BoundaryCurve> b_hi;
    // Secondary-axis (eta) normalisation scale.  LINEAR by default; set
    // LOG for wide-dynamic-range secondary axes (see region::Region).
    region::AxisScale secondary{region::AxisScale::LINEAR};

    RegionSpec() = default;
    RegionSpec(region::AxisTransform p, std::unique_ptr<region::BoundaryCurve> lo, std::unique_ptr<region::BoundaryCurve> hi,
               region::AxisScale sec = region::AxisScale::LINEAR)
      : primary(p), b_lo(std::move(lo)), b_hi(std::move(hi)), secondary(sec) {}
    RegionSpec(const RegionSpec&) = delete;
    RegionSpec& operator=(const RegionSpec&) = delete;
    RegionSpec(RegionSpec&&) = default;
    RegionSpec& operator=(RegionSpec&&) = default;
    ~RegionSpec() = default;
};

// One output property to tabulate on the surface, plus the transform
// the SVD evaluator should apply at lookup time.
//
//   IDENTITY : value = Σ_k σ_k u_k v_k                  (default)
//   EXP      : value = exp(Σ_k σ_k u_k v_k)             (log-fit)
//
// EXP is the right choice for strictly-positive properties with wide
// dynamic range (density, pressure).  IDENTITY for properties that
// can change sign or that vary by orders of magnitude only locally
// (enthalpy near reference state, entropy, internal energy).
struct PropertySpec
{
    ::CoolProp::parameters key = ::CoolProp::iDmass;
    svd::OutputTransform transform = svd::OutputTransform::IDENTITY;
};

// Complete recipe for one SVDSurface.  Captures both the geometry
// (how regions are bounded) and the thermodynamics (how to sample
// HEOS at any (a, b) and which properties to read out).
//
// The two callbacks `update_state` and `read_property` are the
// thermodynamics-aware glue that keeps the rest of the library
// input-pair-agnostic:
//
//   update_state(heos, a, b)
//     ↑ takes the surface's (a, b) inputs and puts the HEOS state
//       into the matching configuration via the right
//       AbstractState::update() input-pair call.
//
//   read_property(heos, key)
//     ↑ reads `key` off the already-updated HEOS state.
//
// PH preset uses `heos.update(HmassP_INPUTS, h, p)` then
// `heos.rhomass()` / `heos.T()` / `heos.smass()` / `heos.umass()`.
// PT preset uses `heos.update(PT_INPUTS, p, T)` then
// `heos.rhomass()` / `heos.hmass()` / `heos.smass()` / `heos.umass()`.
// DU / HS / PS presets plug in by writing their own pair of lambdas;
// no other files in this library need to change.
struct SurfaceSpec
{
    std::string fluid_name;
    // Source-backend name ("HEOS" / "REFPROP" / "IF97") — used by
    // sample_grid to spawn per-thread AbstractState instances when
    // the PARALLEL_SVDSBTL_SAMPLING config flag is set.  Empty value
    // forces serial sampling regardless of the config (no factory
    // metadata to clone from).
    std::string source_backend;
    // The unscoped CoolProp::input_pairs enum lacks a portable
    // "sentinel" value — initialise to PT_INPUTS which is always a
    // valid enumerator, then expect builders to overwrite.
    ::CoolProp::input_pairs input_pair = ::CoolProp::PT_INPUTS;
    std::vector<RegionSpec> regions;
    std::vector<PropertySpec> properties;

    // Thermodynamics-aware glue.  See SurfaceSpec docstring above.
    // update_state should throw on two-phase / HEOS failure; the
    // factory's per-cell loop catches and treats that cell as NaN.
    std::function<void(::CoolProp::AbstractState&, double a, double b)> update_state;
    std::function<double(::CoolProp::AbstractState&, ::CoolProp::parameters)> read_property;

    // SVD grid sizing.  The defaults match Phase 2a's production
    // resolution (Water max ≈ 7e-5 fractional with CubicSpline sat
    // boundaries).  Override per-fluid when needed.
    std::size_t NT = 200;
    std::size_t NR = 800;
    std::int32_t rank = 20;

    SurfaceSpec() = default;
    SurfaceSpec(const SurfaceSpec&) = delete;
    SurfaceSpec& operator=(const SurfaceSpec&) = delete;
    SurfaceSpec(SurfaceSpec&&) = default;
    SurfaceSpec& operator=(SurfaceSpec&&) = default;
    ~SurfaceSpec() = default;
};

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SURFACE_SPEC_H
