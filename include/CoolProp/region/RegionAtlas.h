#ifndef COOLPROP_REGION_REGION_ATLAS_H
#define COOLPROP_REGION_REGION_ATLAS_H

#include <cstddef>
#include <vector>

#include "CoolProp/region/Region.h"

namespace CoolProp {
namespace region {

// A registry of Regions with bounding-box-first dispatch.
//
// find_region(a, b) scans the regions in registration order:
//   1. For each region, check the cheap O(1) AABB membership.  Skip if
//      the point is outside the AABB.
//   2. On an AABB hit, call the region's curve_contains (which evaluates
//      the two BoundaryCurves at a) to confirm or reject.
//   3. Return the first region whose curve check passes, or -1 if none
//      claims the point.
//
// AABBs are stored separately as a flat SoA in a small vector so the
// first-pass scan is cache-friendly; the per-region BoundaryCurves are
// only touched when a hit needs confirming.
//
// Registration order is part of the public contract: when regions are
// expected to be disjoint but their AABBs overlap (common — the
// bounding boxes of e.g. LIQUID and VAPOR overlap in (h, p) even though
// the regions themselves are disjoint), the first match wins.  Callers
// who want to detect ambiguity can use `find_all_curve_hits` (returns a
// vector of all curve-passing region indices), which is intended for
// debug / test code, not the hot path.
class RegionAtlas
{
   public:
    RegionAtlas() = default;

    // Move only — Regions hold std::unique_ptr<BoundaryCurve>.
    RegionAtlas(const RegionAtlas&) = delete;
    RegionAtlas& operator=(const RegionAtlas&) = delete;
    RegionAtlas(RegionAtlas&&) = default;
    RegionAtlas& operator=(RegionAtlas&&) = default;
    ~RegionAtlas() = default;

    // Append a new region.  Returns its integer index.
    std::size_t add(Region region);

    // Returns the index of the first region whose AABB contains (a, b)
    // AND whose curve envelope contains (a, b); -1 if none.
    [[nodiscard]] int find_region(double a, double b) const noexcept;

    // Debug / test helper: return the indices of every region whose
    // curve envelope contains (a, b).  Useful for asserting disjointness.
    [[nodiscard]] std::vector<std::size_t> find_all_curve_hits(double a, double b) const;

    [[nodiscard]] std::size_t size() const noexcept {
        return regions_.size();
    }
    [[nodiscard]] const Region& region(std::size_t i) const noexcept {
        return regions_[i];
    }

   private:
    std::vector<Region> regions_;
    // SoA AABB cache for fast first-pass scanning.
    std::vector<double> a_lo_;
    std::vector<double> a_hi_;
    std::vector<double> b_min_;
    std::vector<double> b_max_;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_REGION_ATLAS_H
