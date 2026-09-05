#include "CoolProp/sbtl/SVDSurfaceCache.h"

#include <algorithm>

#include "CoolProp/Configuration.h"

namespace CoolProp {
namespace sbtl {

std::size_t estimate_surface_bytes(const SVDSurface& surface) {
    std::size_t doubles = 0;
    const std::size_t n_regions = surface.region_count();
    const auto& props = surface.properties();
    for (std::size_t region_idx = 0; region_idx < n_regions; ++region_idx) {
        for (const auto prop : props) {
            const auto& decomp = surface.decomposition(region_idx, prop);
            doubles += decomp.x_grid.size() + decomp.y_grid.size() + decomp.U.size() + decomp.dU_dx.size() + decomp.V_S.size() + decomp.dV_S_dy.size()
                       + decomp.S.size();
        }
    }
    return doubles * sizeof(double);
}

SVDSurfaceCache& SVDSurfaceCache::instance() {
    static SVDSurfaceCache cache;
    return cache;
}

std::shared_ptr<const SVDSurface> SVDSurfaceCache::get(const std::string& key) {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = entries_.find(key);
    if (it == entries_.end()) {
        return nullptr;
    }
    lru_.splice(lru_.begin(), lru_, it->second.lru_it);
    return it->second.surface;
}

void SVDSurfaceCache::put(const std::string& key, std::shared_ptr<const SVDSurface> surface) {
    if (!surface) {
        return;
    }
    const std::size_t bytes = estimate_surface_bytes(*surface);

    std::lock_guard<std::mutex> lock(mutex_);

    const auto existing = entries_.find(key);
    if (existing != entries_.end()) {
        total_bytes_ -= existing->second.bytes;
        lru_.erase(existing->second.lru_it);
        entries_.erase(existing);
    }

    lru_.push_front(key);
    Entry entry;
    entry.surface = std::move(surface);
    entry.bytes = bytes;
    entry.lru_it = lru_.begin();
    total_bytes_ += bytes;
    entries_.emplace(key, std::move(entry));

    evict_to_fit_locked_();
}

void SVDSurfaceCache::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    entries_.clear();
    lru_.clear();
    total_bytes_ = 0;
}

std::size_t SVDSurfaceCache::entry_count() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return entries_.size();
}

std::size_t SVDSurfaceCache::size_bytes() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return total_bytes_;
}

void SVDSurfaceCache::evict_to_fit_locked_() {
    const std::size_t max_entries = static_cast<std::size_t>(std::max(0, ::CoolProp::get_config_int(SVDSBTL_SURFACE_CACHE_MAX_ENTRIES)));
    const std::size_t max_bytes =
      static_cast<std::size_t>(std::max(0, ::CoolProp::get_config_int(SVDSBTL_SURFACE_CACHE_MAX_SIZE_MB))) * std::size_t{1024} * std::size_t{1024};

    while (!lru_.empty() && (entries_.size() > max_entries || total_bytes_ > max_bytes)) {
        const std::string victim_key = lru_.back();
        const auto it = entries_.find(victim_key);
        if (it != entries_.end()) {
            total_bytes_ -= it->second.bytes;
            entries_.erase(it);
        }
        lru_.pop_back();
    }
}

}  // namespace sbtl
}  // namespace CoolProp
