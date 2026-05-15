#include "CoolProp/region/RegionAtlas.h"

namespace CoolProp {
namespace region {

std::size_t RegionAtlas::add(Region region) {
    const auto& bbox = region.bbox();
    a_lo_.push_back(bbox.a_lo);
    a_hi_.push_back(bbox.a_hi);
    b_min_.push_back(bbox.b_min);
    b_max_.push_back(bbox.b_max);
    regions_.push_back(std::move(region));
    return regions_.size() - 1;
}

int RegionAtlas::find_region(double a, double b) const noexcept {
    const std::size_t n = regions_.size();
    for (std::size_t i = 0; i < n; ++i) {
        // AABB filter — cache-friendly first pass.
        if (a < a_lo_[i] || a > a_hi_[i] || b < b_min_[i] || b > b_max_[i]) {
            continue;
        }
        // AABB hit — confirm against the curved envelope.
        if (regions_[i].curve_contains(a, b)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

std::vector<std::size_t> RegionAtlas::find_all_curve_hits(double a, double b) const {
    std::vector<std::size_t> hits;
    const std::size_t n = regions_.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (a < a_lo_[i] || a > a_hi_[i] || b < b_min_[i] || b > b_max_[i]) {
            continue;
        }
        if (regions_[i].curve_contains(a, b)) {
            hits.push_back(i);
        }
    }
    return hits;
}

}  // namespace region
}  // namespace CoolProp
