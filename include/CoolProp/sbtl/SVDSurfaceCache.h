#ifndef COOLPROP_SBTL_SVD_SURFACE_CACHE_H
#define COOLPROP_SBTL_SVD_SURFACE_CACHE_H

#include <cstddef>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

#include "CoolProp/sbtl/SVDSurface.h"

namespace CoolProp {
namespace sbtl {

// Process-wide, thread-safe, bounded LRU cache of sealed SVDSurface
// objects shared across SVDSBTLBackend instances within one process.
//
// Without this cache, ensure_surface_() re-reads and re-deserializes
// the ~80ms-per-surface .svd.bin.z disk file on EVERY backend
// construction, even when a prior AbstractState in the same process
// already loaded the identical fluid/source/input-pair/options
// combination (the dominant case for PropsSI-style call patterns,
// which rebuild the backend per call). This cache sits in front of
// that disk load so repeat requests for the same surface are O(map
// lookup) instead of O(disk read + msgpack/zlib deserialize).
//
// Keyed by the on-disk cache path string (already encodes fluid +
// source backend + input pair + options hash via
// SVDSurfaceSerializer::default_cache_path), so cache identity is
// always exactly the disk-cache identity -- no separate key
// derivation to keep in sync with the path-construction logic.
//
// Bounded by BOTH an entry count and a total estimated byte size;
// whichever limit is hit first triggers LRU eviction. Surface size
// varies widely with the NT/NR/rank build options (a handful of MB
// to tens of MB per surface), so an entry-count cap alone can't
// bound memory. Limits are read live from the
// SVDSBTL_SURFACE_CACHE_MAX_ENTRIES / SVDSBTL_SURFACE_CACHE_MAX_SIZE_MB
// configuration keys on every put(), matching the existing
// SVDSBTL_SAMPLING_THREADS pattern (live config reads, no separate
// setter API).
//
// Surfaces are immutable once published here: add_region() /
// add_region_property_svd() / seal() are only ever called by the
// builder thread before a surface is handed to put(), so sharing one
// instance read-only across threads via shared_ptr<const SVDSurface>
// is safe without additional locking on the surface itself.
class SVDSurfaceCache
{
   public:
    static SVDSurfaceCache& instance();

    // Returns the cached surface for `key`, or nullptr on a miss. A
    // hit moves the entry to the most-recently-used position.
    std::shared_ptr<const SVDSurface> get(const std::string& key);

    // Publishes `surface` under `key`, evicting least-recently-used
    // entries first as needed to respect the current count/size
    // configuration. If `key` is already present, the existing entry
    // is replaced (last-write-wins; ensure_surface_'s own per-instance
    // idempotency check means a given backend only calls put() once
    // per key, but two backends racing to build the same missing
    // surface concurrently is a correctness-neutral, merely wasteful,
    // possibility this resolves safely rather than via UB).
    void put(const std::string& key, std::shared_ptr<const SVDSurface> surface);

    // Drops every cached surface. Test-only escape hatch so tests can
    // assert a clean cache state before exercising eviction behavior.
    void clear();

    std::size_t entry_count() const;
    std::size_t size_bytes() const;

   private:
    SVDSurfaceCache() = default;

    // Must be called with mutex_ held.
    void evict_to_fit_locked_();

    struct Entry
    {
        std::shared_ptr<const SVDSurface> surface;
        std::size_t bytes = 0;
        std::list<std::string>::iterator lru_it;
    };

    mutable std::mutex mutex_;
    std::list<std::string> lru_;  // front = most recently used, back = next to evict
    std::unordered_map<std::string, Entry> entries_;
    std::size_t total_bytes_ = 0;
};

// Estimated resident heap size of `surface` in bytes: the sum of the
// per-decomposition grid/SVD vectors (x_grid, y_grid, U, dU_dx, V_S,
// dV_S_dy, S) across every (region, property) pair. Dominates actual
// memory use at production NT/NR/rank resolutions; RegionAtlas and
// SVDEvaluator bookkeeping overhead is not counted but is small in
// comparison.
std::size_t estimate_surface_bytes(const SVDSurface& surface);

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_CACHE_H
