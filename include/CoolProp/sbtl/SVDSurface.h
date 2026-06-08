#ifndef COOLPROP_SBTL_SVD_SURFACE_H
#define COOLPROP_SBTL_SVD_SURFACE_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "CoolProp/region/RegionAtlas.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "CoolProp/svd/SVDEvaluator.h"
#include "CoolProp/DataStructures.h"

namespace CoolProp {
namespace sbtl {

// An SVDSurface aggregates everything needed for a one-input-pair
// property lookup on a single pure fluid:
//
//   - a RegionAtlas over the (a, b) plane (LIQUID / VAPOR / SUPER /
//     ...) that tells us which region a query falls into;
//   - one SVDDecomposition per (region, property) pair, heap-allocated
//     so addresses stay stable across vector growth (Phase 2a learned
//     this the hard way — see SVDEvaluator.h's dangling-pointer note);
//   - an SVDEvaluator per decomposition;
//   - metadata: which input pair this surface represents (PT, HP, DU,
//     ...), the fluid name, and the list of output properties.
//
// Lookups go:
//
//   user calls eval(prop, a, b)
//     ↓
//   atlas.find_region(a, b)            → region index or -1
//     ↓
//   region.to_normalized(a, b)         → (xi, eta) in [0,1]^2
//     ↓
//   SVDEvaluator::eval(xi_or_eta_x, xi_or_eta_y)
//                                       → exp(Σ S_k u_k v_k)
//     ↓
//   return ρ / h / s / ...
//
// The two-axis ordering for the SVD (which coordinate is x_grid and
// which is y_grid) is fixed by SurfaceSpec at build time — see
// SVDSurfaceFactory.h.  By the convention used in Phase 2a's e2e tool
// the SVD's x-axis is the *secondary* axis (the η normalized
// secondary coordinate, in [0, 1]) and the y-axis is the *primary*
// axis (the ξ primary coordinate, also in [0, 1]).
class SVDSurface
{
   public:
    // Construct an empty surface.  Add regions + per-region per-
    // property decompositions via add_region_property_svd().  Finalise
    // with seal() — after which no further mutations are allowed and
    // the evaluators (lazily-constructed from the heap-allocated
    // decompositions) become callable.
    //
    // The properties argument fixes the property dimension; every
    // region must register exactly one decomposition per property in
    // the same order.  Reorderable downstream via property indexing.
    SVDSurface(std::string fluid_name, ::CoolProp::input_pairs input_pair, std::vector<::CoolProp::parameters> properties);

    // Non-copyable; the heap-resident decompositions are unique_ptr
    // owned.  Movable.
    SVDSurface(const SVDSurface&) = delete;
    SVDSurface& operator=(const SVDSurface&) = delete;
    SVDSurface(SVDSurface&&) = default;
    SVDSurface& operator=(SVDSurface&&) = default;
    ~SVDSurface() = default;

    // Add a new region.  Returns the region index assigned by the
    // underlying RegionAtlas.  Once a region is added, its property
    // decompositions must be supplied via add_region_property_svd
    // before seal() is called.
    std::size_t add_region(region::Region region);

    // Register one property SVD for an already-added region.
    // `region_idx` must match the value returned by add_region();
    // `prop` must be one of the properties passed to the constructor.
    void add_region_property_svd(std::size_t region_idx, ::CoolProp::parameters prop, svd::SVDDecomposition decomp);

    // Lock the surface for evaluation.  Validates that every region
    // has registered every property.  After seal() returns, eval()
    // and resolve() are callable; add_* methods throw.
    void seal();

    // Lookup `prop` at (a, b).  Returns NaN if (a, b) is outside every
    // region of the atlas.  Throws std::invalid_argument if `prop`
    // isn't one of the properties this surface stores.
    [[nodiscard]] double eval(::CoolProp::parameters prop, double a, double b) const;

    // Two-step lookup for callers that want to query several
    // properties at the same (a, b) without paying for atlas dispatch
    // each time.  Returns (region_idx, xi, eta).  region_idx = -1 if
    // outside every region.
    struct ResolvedPoint
    {
        int region_idx;
        double svd_x;  // value to feed into SVDEvaluator's x axis
        double svd_y;  // value to feed into SVDEvaluator's y axis
    };
    [[nodiscard]] ResolvedPoint resolve(double a, double b) const noexcept;

    // Fast path after resolve(): evaluate `prop` at a pre-computed
    // (region, svd_x, svd_y).  Caller responsible for ensuring
    // region_idx >= 0.
    [[nodiscard]] double eval_with_region(::CoolProp::parameters prop, int region_idx, double svd_x, double svd_y) const;

    // Batched fast path: evaluate `n` properties at the same
    // (region, svd_x, svd_y).  All per-region per-property
    // SVDEvaluators share the region's (x_grid, y_grid), so the
    // locate() + Hermite-basis setup is done ONCE and amortized
    // across the n property evals — a measurable speedup once n is
    // ≥ 2.  out[i] receives the value for props[i].  Throws on
    // unknown property keys (i.e., props that aren't tabulated on
    // this surface); see contains_property() for the gate.
    void eval_with_region_multi(int region_idx, double svd_x, double svd_y, const ::CoolProp::parameters* props, std::size_t n, double* out) const;

    // Convenience predicates.
    [[nodiscard]] bool contains_property(::CoolProp::parameters prop) const noexcept;

    // Read-only access.
    [[nodiscard]] const std::string& fluid_name() const noexcept {
        return fluid_name_;
    }
    [[nodiscard]] ::CoolProp::input_pairs input_pair() const noexcept {
        return input_pair_;
    }
    [[nodiscard]] const std::vector<::CoolProp::parameters>& properties() const noexcept {
        return properties_;
    }
    [[nodiscard]] const region::RegionAtlas& atlas() const noexcept {
        return atlas_;
    }
    [[nodiscard]] bool sealed() const noexcept {
        return sealed_;
    }
    [[nodiscard]] std::size_t region_count() const noexcept;

    // Direct access to the heap-resident decomposition for a region /
    // property — primarily for serialization.  Throws if either
    // index is out of range.
    [[nodiscard]] const svd::SVDDecomposition& decomposition(std::size_t region_idx, ::CoolProp::parameters prop) const;

   private:
    // Lookup table: property → index into the per-region property
    // array.  Built once at construction.
    [[nodiscard]] std::size_t property_index(::CoolProp::parameters prop) const;

    std::string fluid_name_;
    ::CoolProp::input_pairs input_pair_;
    std::vector<::CoolProp::parameters> properties_;

    region::RegionAtlas atlas_;

    // Per-region storage.  Sized to atlas_.size() after seal().
    // decomps_[region_idx][property_idx] is the heap-resident SVDDecomposition.
    std::vector<std::vector<std::unique_ptr<svd::SVDDecomposition>>> decomps_;
    // Parallel evaluator vector, lazily filled in seal().
    std::vector<std::vector<std::unique_ptr<svd::SVDEvaluator>>> evaluators_;

    bool sealed_ = false;
};

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_H
