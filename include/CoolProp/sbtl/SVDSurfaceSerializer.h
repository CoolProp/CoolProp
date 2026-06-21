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
    //   rev 7: SUPER region split for IF97 source backend — SUPER_R3
    //          (h < h_B23(p)) and SUPER_R2 (h > h_B23(p)) — so the
    //          IF97 R2/R3 derivative kink sits on a region edge
    //          instead of inside a cell.  Region count for IF97 goes
    //          from 3 (LIQUID/VAPOR/SUPER) to 4 or 5 (LIQUID/VAPOR/
    //          SUPER_R3/SUPER_R2 + optional SUPER_HIGH_P above 100 MPa).
    //          On-wire region list differs, so the rev bump forces a
    //          clean rebuild instead of attempting to deserialize a
    //          rev-6 cache that won't have the new boundary curves.
    //          HEOS source unchanged.
    //   rev 8: SUPER_R3 further split at the IF97 R1/R3 isotherm
    //          h_R1R3(p) = h_IF97(623.15 K, p), so R1 territory at
    //          p > pcrit gets its own SVD (SUPER_R1_super) separate
    //          from R3 proper (SUPER_R3_proper).  Region count for
    //          IF97 goes from 4-5 (post-rev-7) to 5-6.  Closes the
    //          post-rev-7 R1 conformance gap where R1 and R3 modes
    //          competed for SVD bandwidth in a single region.
    //   rev 9: IF97 sampling-side Newton replaced with TOMS748
    //          bracketed root-find (foi.9.10).  R3 cells whose Newton
    //          previously failed to converge (e.g., T_target=663.7 K,
    //          p=26.6 MPa where the T=700 fallback seed put Newton in
    //          R2 territory and it oscillated across the R2/R3 boundary)
    //          now sample correctly.  Stored T/s/ρ/w grid values for
    //          those cells change, so existing rev-8 caches must
    //          rebuild.
    //   rev 10: CoolProp-8vg.  HEOS-source presets switched their sat
    //           boundary curves (h_sat,L, h_sat,V) to
    //           region::SuperancillaryBoundaryCurve — no-refit views
    //           onto the SuperAncillary's machine-precision Chebyshev
    //           expansions.  IF97-source presets (which have no SA)
    //           keep the legacy 64-knot cubic spline path.  HEOS
    //           tables built with the new boundaries have different
    //           η normalisation at each grid sample → SVD
    //           coefficients shift → existing rev-9 caches must
    //           rebuild.  HEOS surfaces can't yet round-trip through
    //           the serializer (SuperancillaryBoundaryCurve has no
    //           CurveKind id) — they rebuild in memory each session;
    //           IF97 surfaces still cache normally.
    //   rev 11: CurveKind::SUPERANCILLARY tag added with the State
    //           POD (p_min, p_max, prop_key, Q, output_scale, b_min,
    //           b_max — 8-element array).  The SA handle itself isn't
    //           stored; load-side re-acquires it by constructing a
    //           HEOS AbstractState for the fluid in the stream and
    //           pulling its SuperAncillary.  HEOS surfaces now round-
    //           trip through disk (CoolProp-cv7); first-session cost
    //           amortises across all subsequent process invocations.
    //           Rev bump invalidates rev-10 caches, which never
    //           successfully landed any HEOS surfaces on disk anyway
    //           (save() threw on the unknown subclass).
    //   rev 12: SVDSBTL&IF97 presets gain a SUPER_R5 region covering
    //           IAPWS R7-97 Region 5 (T ∈ [1073.15, 2273.15] K, p ≤
    //           50 MPa — CoolProp-pd6).  The new region appears at
    //           the tail of the regions array (per-PH and per-PT),
    //           so the region count for IF97 caches changes from
    //           5-6 to 6-7.  Rev bump forces rebuild so the loader
    //           doesn't try to dispatch lookups for R5 cells against
    //           a rev-11 cache that has no SUPER_R5 region.  HEOS
    //           caches unchanged in content but invalidated by the
    //           rev bump too (R5 is IF97-only and HEOS presets skip
    //           the new region entirely; the cache geometry is the
    //           same as rev 11 but the rev field differs).
    //   rev 13: CoolProp-4u9.  HEOS / REFPROP source presets gain a
    //           pair of near-critical sub-regions (NC_LIQUID,
    //           NC_VAPOR) on p ∈ [0.9·pc, (1 − 1ppm)·pc] using a new
    //           AxisScale::POWER primary axis (β = 1/3, cube-root
    //           crowding toward pc).  Parent LIQUID/VAPOR p_max
    //           clipped down to 0.9·pc to hand off cleanly.  Region
    //           count for HEOS caches goes from 3 (LIQUID/VAPOR/SUPER)
    //           to 5.  POWER axis serialises via the same
    //           AxisTransform pack/unpack as LINEAR/LOG (scale enum
    //           value differs; a_lo/a_hi unchanged); existing rev-12
    //           caches don't know about the NC regions and would
    //           dispatch near-pc lookups to the parent LIQUID/VAPOR
    //           SVD where off-grid max error is ~1e-3 instead of the
    //           new ~1e-7 from the POWER NC regions, so the rev bump
    //           forces a clean rebuild.  IF97 caches unchanged in
    //           content (no NC for IF97 — G13-15 already passes) but
    //           invalidated by the rev bump too.
    //   rev 14: CoolProp-4u9 follow-up.  Add NC_SUPER region on
    //           p ∈ [(1 + 1e-10)·pc, 1.1·pc] using AxisScale::POWER_LO
    //           (mirror of POWER that crowds toward a_lo).  Closes a
    //           (T, p) coverage gap exposed by the h=350 kJ/kg R245fa
    //           sweep: cells in the thin strip just above pc with T
    //           outside the auto-cal'd patch's narrow (T) bbox fell
    //           through every region (LIQUID/VAPOR clipped to 0.9·pc,
    //           NC_LIQUID/VAPOR capped at (1−1ppm)·pc, SUPER starting
    //           at 1.001·pc) and returned NaN.  NC_SUPER claims them.
    //           HEOS region count goes from 5 to 6.  NC sub-side band
    //           also tightened from 1ppm to 1e-10 below pc — the NC
    //           POWER axis is well-behaved at a_hi=pc-ε for any ε
    //           large enough to keep the SuperAncillary sat-curve
    //           well-defined, and tightening reduces the patch-only
    //           sliver around pc.
    //   rev 16: CoolProp-wvtz.  Region gains a selectable secondary-axis
    //           (eta) scale, packed as a 9th region-blob element.  The
    //           DmassT preset's VAPOR + SUPER regions switch to
    //           AxisScale::LOG so the near-ideal-gas low-density tail
    //           (rho_hi/rho_lo ~ 1e5, p ∝ rho) gets uniform-in-decade
    //           sampling instead of collapsing below the first linear-eta
    //           grid node.  The stored grid sample positions and U/slopes
    //           for those regions change, and the on-wire region array
    //           grows from 8 to 9 elements, so rev-15 caches must rebuild.
    //           HEOS DmassT low-density pressure error drops from
    //           ~40 %–2400 % to ~1e-5.  PT / HmassP surfaces are
    //           byte-identical (their regions stay LINEAR).
    //   rev 17: CoolProp-4z79.  The DmassT preset gains three near-critical
    //           (NC) sub-regions that close the previously-uncovered
    //           ±0.1%·Tc band: POWER-axis NC_LIQUID/NC_VAPOR carry the
    //           dome-bounded sides up to ~Tc, and a POWER_LO isobar-bounded
    //           NC_SUPER spans the critical isotherm itself.  HEOS-source
    //           DmassT region count grows (3 → 6), so the on-wire region
    //           list differs and rev-16 caches must rebuild.  The whole
    //           [0.99,1.01]·Tc band is now table-served (no HEOS) to
    //           ~1e-6..2e-5 for all fluids.  PT / HmassP and the non-DmassT
    //           presets are unchanged.
    //   rev 18: CoolProp-jh6a.  The DmassT parent SUPER region's primary-T
    //           axis switches LINEAR → LOG.  For low-Tc / wide-supercritical
    //           fluids (Helium: Tc 5.2 K, Tmax 2000 K, ~380× range) LINEAR
    //           starved the near-Tc supercritical zone of grid lines, giving
    //           ~1% pressure error across the whole SUPER region; LOG drops
    //           it to ~1e-8.  No region-count change (same 6 regions), but
    //           the SUPER grid sample positions and U/slopes shift, so
    //           rev-17 caches must rebuild.  Modest-ratio fluids (Water ~3×,
    //           CO2 ~6×) are unchanged in practice (LOG ≈ LINEAR there).
    //   rev 19: CoolProp-naqt / issue #3189.  The subcritical PT / HmassP /
    //           PSmass presets lower their default p_min from
    //           ~1.01-1.03·p_triple to the true triple-point pressure
    //           (heos.p_triple()), so a PT query at p == p_triple() now
    //           resolves instead of returning NaN.  The log-p primary axis
    //           span shifts (lower bottom isobar), so the grid sample
    //           positions and SVD U/slopes change and rev-18 caches must
    //           rebuild.  (A non-default `pmin` option already produces a
    //           distinct cache via the opthash; this rev covers the
    //           all-defaults "{}" case whose opthash is unchanged.)
    static constexpr int kRevision = 19;

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
    // throws if the file doesn't exist or is unreadable.  save_to_file
    // routes through ::write_bytes_atomic (CPfilepaths.h) so concurrent
    // writers in different processes / threads never see a partial-
    // write file (CoolProp-4no.2).
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
