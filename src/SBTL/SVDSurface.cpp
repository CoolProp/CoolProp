#include "CoolProp/sbtl/SVDSurface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace CoolProp {
namespace sbtl {

SVDSurface::SVDSurface(std::string fluid_name, ::CoolProp::input_pairs input_pair, std::vector<::CoolProp::parameters> properties)
  : fluid_name_(std::move(fluid_name)), input_pair_(input_pair), properties_(std::move(properties)) {
    if (properties_.empty()) {
        throw std::invalid_argument("SVDSurface: properties list must be non-empty");
    }
    // Detect duplicates in the property list — silently accepting them
    // means later add_region_property_svd calls would land in the wrong
    // slot.
    auto sorted = properties_;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        throw std::invalid_argument("SVDSurface: duplicate property in properties list");
    }
}

std::size_t SVDSurface::add_region(region::Region region) {
    if (sealed_) {
        throw std::logic_error("SVDSurface: add_region called after seal()");
    }
    const std::size_t idx = atlas_.add(std::move(region));
    // Resize per-region storage to match.
    decomps_.emplace_back();
    decomps_.back().resize(properties_.size());
    evaluators_.emplace_back();
    evaluators_.back().resize(properties_.size());
    return idx;
}

void SVDSurface::add_region_property_svd(std::size_t region_idx, ::CoolProp::parameters prop, svd::SVDDecomposition decomp) {
    if (sealed_) {
        throw std::logic_error("SVDSurface: add_region_property_svd called after seal()");
    }
    if (region_idx >= decomps_.size()) {
        throw std::out_of_range("SVDSurface: region_idx out of range");
    }
    const std::size_t pidx = property_index(prop);
    if (decomps_[region_idx][pidx]) {
        throw std::logic_error("SVDSurface: decomposition for this (region, property) already registered");
    }
    decomps_[region_idx][pidx] = std::make_unique<svd::SVDDecomposition>(std::move(decomp));
}

void SVDSurface::seal() {
    if (sealed_) {
        return;
    }
    // Validate every region has every property.
    for (std::size_t r = 0; r < decomps_.size(); ++r) {
        for (std::size_t p = 0; p < properties_.size(); ++p) {
            if (!decomps_[r][p]) {
                throw std::logic_error("SVDSurface::seal: region " + std::to_string(r) + " is missing a decomposition for property index "
                                       + std::to_string(p));
            }
        }
    }
    // eval_with_region_multi reuses one SVDEvalContext across every
    // per-property evaluator in a region, which is only valid when
    // every decomposition in that region shares an identical
    // (x_grid, y_grid).  SVDSurfaceFactory builds them together so
    // this holds by construction, but a hand-built surface or one
    // loaded from disk could violate it — validate once at seal()
    // rather than trusting it forever.
    for (std::size_t r = 0; r < decomps_.size(); ++r) {
        const auto& x_ref = decomps_[r][0]->x_grid;
        const auto& y_ref = decomps_[r][0]->y_grid;
        for (std::size_t p = 1; p < properties_.size(); ++p) {
            if (decomps_[r][p]->x_grid != x_ref || decomps_[r][p]->y_grid != y_ref) {
                throw std::logic_error("SVDSurface::seal: region " + std::to_string(r) + " property index " + std::to_string(p)
                                       + " has (x_grid, y_grid) that differ from the region's first property — "
                                         "eval_with_region_multi requires identical grids across the region");
            }
        }
    }
    // Build evaluators against the heap-resident decompositions.  The
    // address of each unique_ptr's pointee is stable for the lifetime
    // of the unique_ptr (and the surface).
    for (std::size_t r = 0; r < decomps_.size(); ++r) {
        for (std::size_t p = 0; p < properties_.size(); ++p) {
            evaluators_[r][p] = std::make_unique<svd::SVDEvaluator>(*decomps_[r][p]);
        }
    }
    sealed_ = true;
}

std::size_t SVDSurface::region_count() const noexcept {
    return atlas_.size();
}

bool SVDSurface::contains_property(::CoolProp::parameters prop) const noexcept {
    return std::find(properties_.begin(), properties_.end(), prop) != properties_.end();
}

std::size_t SVDSurface::property_index(::CoolProp::parameters prop) const {
    const auto it = std::find(properties_.begin(), properties_.end(), prop);
    if (it == properties_.end()) {
        throw std::invalid_argument("SVDSurface: property is not tabulated on this surface");
    }
    return static_cast<std::size_t>(std::distance(properties_.begin(), it));
}

const svd::SVDDecomposition& SVDSurface::decomposition(std::size_t region_idx, ::CoolProp::parameters prop) const {
    if (region_idx >= decomps_.size()) {
        throw std::out_of_range("SVDSurface::decomposition: region_idx out of range");
    }
    const std::size_t pidx = property_index(prop);
    if (!decomps_[region_idx][pidx]) {
        throw std::logic_error("SVDSurface::decomposition: decomposition is null (surface not sealed?)");
    }
    return *decomps_[region_idx][pidx];
}

SVDSurface::ResolvedPoint SVDSurface::resolve(double a, double b) const noexcept {
    const int region_idx = atlas_.find_region(a, b);
    if (region_idx < 0) {
        return ResolvedPoint{-1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    }
    const auto [xi, eta] = atlas_.region(static_cast<std::size_t>(region_idx)).to_normalized(a, b);
    // SVD convention (matches Phase 2a's e2e tool): x-axis = eta
    // (secondary normalised), y-axis = xi (primary normalised).
    return ResolvedPoint{region_idx, eta, xi};
}

double SVDSurface::eval(::CoolProp::parameters prop, double a, double b) const {
    if (!sealed_) {
        throw std::logic_error("SVDSurface::eval called before seal()");
    }
    const auto resolved = resolve(a, b);
    if (resolved.region_idx < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return eval_with_region(prop, resolved.region_idx, resolved.svd_x, resolved.svd_y);
}

double SVDSurface::eval_with_region(::CoolProp::parameters prop, int region_idx, double svd_x, double svd_y) const {
    if (!sealed_) {
        throw std::logic_error("SVDSurface::eval_with_region called before seal()");
    }
    if (region_idx < 0 || static_cast<std::size_t>(region_idx) >= evaluators_.size()) {
        throw std::out_of_range("SVDSurface::eval_with_region: region_idx out of range");
    }
    const std::size_t pidx = property_index(prop);
    return evaluators_[static_cast<std::size_t>(region_idx)][pidx]->eval(svd_x, svd_y);
}

void SVDSurface::eval_with_region_multi(int region_idx, double svd_x, double svd_y, const ::CoolProp::parameters* props, std::size_t n,
                                        double* out) const {
    if (!sealed_) {
        throw std::logic_error("SVDSurface::eval_with_region_multi called before seal()");
    }
    if (region_idx < 0 || static_cast<std::size_t>(region_idx) >= evaluators_.size()) {
        throw std::out_of_range("SVDSurface::eval_with_region_multi: region_idx out of range");
    }
    if (n == 0) {
        return;
    }
    const auto& region_evals = evaluators_[static_cast<std::size_t>(region_idx)];
    // All per-property evaluators in this region share the same
    // (x_grid, y_grid) by construction (SVDSurfaceFactory builds them
    // together), so a single make_context() call covers all `n`
    // outputs.  Pick the first property's evaluator as the context
    // source — picking the same one every time is fine, but using
    // props[0] avoids surprising readers who'd expect the prop they
    // asked for to drive the context.
    const std::size_t first_pidx = property_index(props[0]);
    const svd::SVDEvalContext ctx = region_evals[first_pidx]->make_context(svd_x, svd_y);
    out[0] = region_evals[first_pidx]->eval_with_context(ctx);
    for (std::size_t k = 1; k < n; ++k) {
        const std::size_t pidx = property_index(props[k]);
        out[k] = region_evals[pidx]->eval_with_context(ctx);
    }
}

}  // namespace sbtl
}  // namespace CoolProp
