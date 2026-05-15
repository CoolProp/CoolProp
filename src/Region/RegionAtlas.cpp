#include "CoolProp/region/RegionAtlas.h"

namespace CoolProp {
namespace region {

std::size_t RegionAtlas::add(Region region) {
    // Transactional insertion: every container must end up with the
    // same size, or none of them grow.  A throw mid-sequence would
    // otherwise leave the SoA AABB cache and regions_ out of sync, so
    // subsequent find_region calls would mis-dispatch.
    //
    // Snapshot the BBox before the std::move below — Region's custom
    // move op invalidates the source's bbox_ (to defend against
    // dangling-curve-pointer accesses on moved-from regions), so a
    // reference taken pre-move would read sentinel values after.
    const Region::BBox bbox = region.bbox();
    const std::size_t old_size = regions_.size();
    try {
        regions_.push_back(std::move(region));
        a_lo_.push_back(bbox.a_lo);
        a_hi_.push_back(bbox.a_hi);
        b_min_.push_back(bbox.b_min);
        b_max_.push_back(bbox.b_max);
    } catch (...) {
        // Region has a deleted copy constructor (it owns std::unique_ptr
        // members), so resize() can't be used here — it would require
        // DefaultInsertable on the value type for the (potential) grow
        // path.  erase() is the right tool: shrink-only, no default
        // construction needed.
        if (regions_.size() > old_size) {
            regions_.erase(regions_.begin() + static_cast<std::ptrdiff_t>(old_size), regions_.end());
        }
        if (a_lo_.size() > old_size) {
            a_lo_.erase(a_lo_.begin() + static_cast<std::ptrdiff_t>(old_size), a_lo_.end());
        }
        if (a_hi_.size() > old_size) {
            a_hi_.erase(a_hi_.begin() + static_cast<std::ptrdiff_t>(old_size), a_hi_.end());
        }
        if (b_min_.size() > old_size) {
            b_min_.erase(b_min_.begin() + static_cast<std::ptrdiff_t>(old_size), b_min_.end());
        }
        if (b_max_.size() > old_size) {
            b_max_.erase(b_max_.begin() + static_cast<std::ptrdiff_t>(old_size), b_max_.end());
        }
        throw;
    }
    return old_size;
}

int RegionAtlas::find_region(double a, double b) const noexcept {
    const std::size_t n = regions_.size();
    for (std::size_t i = 0; i < n; ++i) {
        // AABB filter — cache-friendly first pass.  Positive form so
        // NaN inputs are correctly rejected (every comparison returns
        // false; the previous negated form `a < a_lo_[i] || ...` was
        // C++ UB on NaN — none of the `continue`s fired and locate_piece
        // cast NaN to ptrdiff_t).
        if (!(a >= a_lo_[i] && a <= a_hi_[i] && b >= b_min_[i] && b <= b_max_[i])) {
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
        if (!(a >= a_lo_[i] && a <= a_hi_[i] && b >= b_min_[i] && b <= b_max_[i])) {
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
