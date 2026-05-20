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
    //   rev 4: cache filename now uses the symbolic input_pair name
    //          (e.g. "PT_INPUTS") instead of the raw enum int — closes
    //          CoolProp-b6v — and incorporates a 16-hex FNV-1a 64
    //          opthash over the canonical-JSON options blob.  No
    //          on-wire format change; just a path change.  Old rev-3
    //          int-keyed caches simply won't be picked up.
    //   rev 5: ph_properties / pt_properties extended with transport
    //          properties (η, λ) for IAPWS G13-15 Tables 12 & 13.
    //          Old rev-4 caches are missing two properties at the
    //          tail; bump to force rebuild rather than try to silently
    //          extend them in-place.
    //   rev 6: Chebyshev η-spacing on the secondary axis (cells now
    //          crowd toward the saturation curve) and IF97-source
    //          sampling Newton-refines against forward h(T, p) at
    //          build time.  Both change the *content* of the stored
    //          U-matrix and slopes at the same hashed options key, so
    //          the rev bump is the cache-invalidation mechanism — old
    //          rev-5 caches would silently serve uniform-η surfaces
    //          and undo the IF97 conformance gains.
    static constexpr int kRevision = 6;

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
    //   "<fluid>.<source>.<input_pair_name>.<opthash>.svd.bin.z"
    //
    // - input_pair_name is the symbolic enum name returned by
    //   CoolProp::get_input_pair_short_desc() (e.g. "PT_INPUTS",
    //   "HmassP_INPUTS"), so reordering the input_pairs enum in
    //   DataStructures.h can't silently misalign caches (was the
    //   raw int previously — closes CoolProp-b6v).
    // - opthash is a 16-hex string (typically FNV-1a 64 of the
    //   canonical options blob), or "no_opts" for callers that don't
    //   care about per-options caching.  Defaults to "no_opts" to
    //   keep the legacy two-arg-ish call sites compiling without
    //   touching every test.
    //
    // source_backend must be non-empty and free of path-separator
    // characters; opthash is rejected if it contains anything other
    // than [0-9a-f_].
    static std::string default_cache_path(const std::string& fluid_name, const std::string& source_backend, ::CoolProp::input_pairs input_pair,
                                          const std::string& opthash = "no_opts");
};

}  // namespace sbtl
}  // namespace CoolProp

#endif  // COOLPROP_SBTL_SVD_SURFACE_SERIALIZER_H
