#ifndef COOLPROP_SBTL_SVD_SURFACE_SERIALIZER_H
#define COOLPROP_SBTL_SVD_SURFACE_SERIALIZER_H

#include <string>
#include <vector>

#include "CoolProp/sbtl/SVDSurface.h"

namespace CoolProp {
namespace sbtl {

// Persistence layer for SVDSurface.  Writes a msgpack stream wrapped
// in a zlib-compressed payload, matching the existing TabularBackends
// `.bin.z` convention.
//
// File format (compressed bytes contain a single msgpack array):
//
//   [
//     "SVDS",          // magic, string
//     1,               // revision, int
//     <fluid_name>,    // string
//     [                // surfaces, array (length n_surfaces)
//       <surface_0>,
//       ...
//     ]
//   ]
//
// Each <surface_k> is a flat positional msgpack array:
//
//   [
//     input_pair,                    // int (CoolProp::input_pairs)
//     [[prop_key], ...],             // properties (length-1 arrays;
//                                    //  per-property OutputTransform
//                                    //  lives inside each decomp_blob)
//     [region_blob, ...],            // regions
//     [decomp_blob, ...]             // per (region, property) SVDs
//   ]
//
// Region geometry + boundary curves serialize via the State POD
// snapshots exposed on ConstantCurve / CubicSplineCurve /
// PiecewiseChebyshevCurve.  Each curve is tagged with a uint8 kind
// discriminator so the deserializer knows which from_state() to
// call.
//
// Backward / forward compatibility:
//   - The magic + revision header lets the loader reject obviously
//     non-SVD files and lets future revisions add fields.
//   - A loader on rev N hitting a stream written for rev N+M with
//     extra trailing fields will silently ignore them (msgpack arrays
//     are length-prefixed).  A loader on rev N+M reading rev N data
//     will throw on the missing fields — adopters bump revision when
//     they add anything.
//
// File extension convention: `.svd.bin.z`, parallel to BicubicBackend's
// `.bin.z`.  Stored at
// `${HOME}/.CoolProp/SVDTables/<fluid>.<source>.<input_pair>.svd.bin.z`
// by the helper paths below — see default_cache_path() for the exact
// composition (source is the truth-source backend name; input_pair is
// the integer value of the enum).

class SVDSurfaceSerializer
{
   public:
    // On-wire revision.  Bumped on:
    //   rev 1: initial Phase 2b release (4 properties per surface)
    //   rev 2: speed_sound added as 5th property (rev 1 caches still
    //          deserialize as 4-property surfaces but calc_speed_sound
    //          would silently fall through to NaN; bumping the rev
    //          forces a clean rebuild so the user can't trip on stale
    //          caches missing w).
    //   rev 3: explicit source-of-truth backend in the cache filename
    //          ("<fluid>.<source>.<input_pair>.svd.bin.z").  No on-wire
    //          format change, just a path change — old caches at the
    //          rev-2 path simply won't be picked up and will be
    //          rebuilt under the new path the first time they're
    //          requested.
    static constexpr int kRevision = 3;

    // Pack one surface into a zlib-compressed msgpack blob.
    static std::vector<char> save(const SVDSurface& surface);

    // Reverse: parse a compressed blob back into a fully-sealed
    // SVDSurface.  Throws std::runtime_error on:
    //   - zlib decompression failure
    //   - magic / revision mismatch
    //   - missing or malformed fields
    //   - invalid curve / decomp dimensions surfacing through the
    //     trusted from_state() factories
    static SVDSurface load(const std::vector<char>& compressed);

    // File-level convenience.  save_to_file overwrites; load_from_file
    // throws if the file doesn't exist or is unreadable.
    static void save_to_file(const SVDSurface& surface, const std::string& path);
    static SVDSurface load_from_file(const std::string& path);

    // Returns the standard cache directory:
    //   ${HOME}/.CoolProp/SVDTables/
    // Creates the directory if it doesn't exist.  Mirrors the
    // existing BicubicBackend pattern (see TabularBackends.h:1035).
    static std::string default_cache_dir();

    // Compose default_cache_dir() with
    //   "<fluid>.<source>.<input_pair>.svd.bin.z"
    // so HEOS-built, REFPROP-built, and IF97-built tables for the
    // same fluid get distinct files, and PH and PT surfaces for the
    // same (fluid, source) get distinct files.  source_backend must
    // be non-empty and free of path-separator characters.
    static std::string default_cache_path(const std::string& fluid_name, const std::string& source_backend, ::CoolProp::input_pairs input_pair);
};

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_SERIALIZER_H
