// SVDSBTLBackend implementation.  Architecturally inherits from
// AbstractState (mirroring IF97Backend), holds an SVDSurface per
// supported input pair, and dispatches calc_* through the active
// surface resolved at update() time.

#include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "boost/math/tools/toms748_solve.hpp"

#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Configuration.h"
#include "CoolProp/Hash.h"
#include "CoolProp/SchemaValidation.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#include "CPfilepaths.h"
#include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#include "CoolProp/schemas/SVDSBTLOptions.h"
#include "DataStructures.h"
#include "Exceptions.h"
#include "rapidjson_include.h"

namespace CoolProp {

namespace {

// The MVP set of input pairs.  Adding DU / HS / etc. in a follow-up
// is a one-line append here plus the matching preset + load helpers
// in `build_surface_for_pair_`.
constexpr std::array<CoolProp::input_pairs, 2> kSupportedPairs = {CoolProp::HmassP_INPUTS, CoolProp::PT_INPUTS};

// Resolved grid-shape knobs extracted from the validated options
// document.  Used both to size the per-input-pair surfaces and to
// stamp the canonical-options form into the cache filename.
struct ResolvedGrid
{
    std::size_t NT;
    std::size_t NR;
    std::int32_t rank;
};

constexpr ResolvedGrid kDefaultGrid = {200, 800, 20};

// Polish the patch backend's T after an HmassP_INPUTS update so the
// returned state is forward-consistent in h(T, p), not just R7-97-
// backward-consistent.  Without this, in-bbox queries inherit the
// patch source's backward-equation residual (±25 mK for IF97 R1/R3,
// ±10 mK for R2/R5) — visible as ~10-20 mK / 0.1-0.4 J/(kg·K) s
// residuals in the SVDSBTL&IF97 fail-map even inside the patch box.
// Mirrors the polish in src/SBTL/SurfacePresets.cpp's IF97 sampling
// lambda (foi.9.10).  Best-effort: on bracket-failure or out-of-range,
// leaves the state as the HmassP backward seed (still functional, just
// at the ±25 mK floor).
void polish_patch_state_(::CoolProp::AbstractState& s, double p, double h) {
    auto resid = [&s, p, h](double T) -> double {
        try {
            s.update(::CoolProp::PT_INPUTS, p, T);
            return s.hmass() - h;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            return std::nan("");
        }
    };
    const double T_seed = s.T();
    // Use backend-reported T limits rather than IF97-specific constants:
    // patch_source_ may be HEOS / REFPROP / IF97 with different envelopes,
    // and pinning to water bounds either silently clips (HEOS with T_min
    // > 273.16) or overshoots (backend with T_max < 2273.15).
    const double T_lo = std::max(T_seed - 0.5, s.Tmin());
    const double T_hi = std::min(T_seed + 0.5, s.Tmax());
    auto restore_seed = [&s, p, T_seed]() noexcept {
        try {
            s.update(::CoolProp::PT_INPUTS, p, T_seed);
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    };
    if (!(T_lo < T_hi)) {
        return;
    }
    const double r_lo = resid(T_lo);
    const double r_hi = resid(T_hi);
    if (!std::isfinite(r_lo) || !std::isfinite(r_hi) || r_lo * r_hi > 0.0) {
        // Bracket missed (rare): restore the backward seed and bail.
        restore_seed();
        return;
    }
    try {
        std::uintmax_t max_iter = 30;
        const auto bracket =
          boost::math::tools::toms748_solve(resid, T_lo, T_hi, r_lo, r_hi, boost::math::tools::eps_tolerance<double>(40), max_iter);
        s.update(::CoolProp::PT_INPUTS, p, 0.5 * (bracket.first + bracket.second));
    } catch (...) {
        // toms748_solve / final update threw: restore the backward seed
        // so the caller sees the ±25 mK floor instead of a half-mutated
        // state at the last failed PT probe.
        restore_seed();
    }
}

// Read `grid.NT/NR/rank` out of a validated options document.  Missing
// keys (and missing `grid` itself) leave the corresponding default in
// place — the validator has already rejected anything ill-formed.
ResolvedGrid resolve_grid(const rapidjson::Document& opts) {
    ResolvedGrid g = kDefaultGrid;
    if (opts.IsObject() && opts.HasMember("grid") && opts["grid"].IsObject()) {
        const auto& grid = opts["grid"];
        if (grid.HasMember("NT")) g.NT = static_cast<std::size_t>(grid["NT"].GetInt());
        if (grid.HasMember("NR")) g.NR = static_cast<std::size_t>(grid["NR"].GetInt());
        if (grid.HasMember("rank")) g.rank = grid["rank"].GetInt();
    }
    return g;
}

// Sample-and-build for a given input pair using the matching preset.
// Grid shape comes from `g` so the caller can drive NT/NR/rank from
// the options blob.  The source_backend string is threaded into the
// SurfaceSpec so sample_grid can spawn per-thread AbstractState
// instances when PARALLEL_SVDSBTL_SAMPLING is on (CoolProp-43h).
CoolProp::sbtl::SVDSurface build_surface_for_pair_(::CoolProp::AbstractState& source, const std::string& source_backend, CoolProp::input_pairs pair,
                                                   const ResolvedGrid& g) {
    namespace cp_sbtl = CoolProp::sbtl;
    cp_sbtl::SurfaceSpec spec;
    switch (pair) {
        case CoolProp::HmassP_INPUTS:
            spec = cp_sbtl::presets::ph_subcritical(source, g.NT, g.NR, g.rank);
            break;
        case CoolProp::PT_INPUTS:
            spec = cp_sbtl::presets::pt_subcritical(source, g.NT, g.NR, g.rank);
            break;
        default:
            throw ValueError("SVDSBTL backend: no preset registered for the requested input pair");
    }
    spec.source_backend = source_backend;
    return cp_sbtl::build_surface(source, std::move(spec));
}

// Parse + validate `options_json` against the SVDSBTL schema.  Returns
// the canonical-JSON form (or "{}" when options_json was empty); throws
// CoolProp::ValueError on schema violation.
std::string parse_and_canonicalise(const std::string& options_json) {
    if (options_json.empty()) {
        return "{}";
    }
    rapidjson::Document opts;
    if (opts.Parse(options_json.c_str(), options_json.size()).HasParseError()) {
        throw ValueError("SVDSBTL options: invalid JSON (offset " + std::to_string(opts.GetErrorOffset()) + ")");
    }
    CoolProp::validate_against_schema(opts, kSVDSBTLOptionsSchemaJson);
    return CoolProp::to_canonical_json(opts);
}

// FNV-1a 64 over the canonical options blob — mirrors the opthash the
// SVDSurface cache uses so the critpatch sidecar lives alongside the
// matching .svd.bin.z files and shares the cache-invalidation lifetime.
// The implementation here is intentionally a tiny duplicate; the
// canonical version lives in SVDSurfaceSerializer but isn't exposed
// for re-use, and the cost of duplicating is one short function.
std::uint64_t opthash_fnv1a_(const std::string& s) noexcept {
    constexpr std::uint64_t kFNVOffset = 0xcbf29ce484222325ULL;
    constexpr std::uint64_t kFNVPrime = 0x100000001b3ULL;
    std::uint64_t h = kFNVOffset;
    for (const unsigned char c : s) {
        h ^= static_cast<std::uint64_t>(c);
        h *= kFNVPrime;
    }
    return h;
}

// Critpatch sidecar path: lives next to the .svd.bin.z files so
// directory-level cache cleanup catches both.  The opthash binds the
// patch to the same options blob the surfaces were built against.
//
// Defense-in-depth against path traversal: rejects any
// fluid_name / source_backend containing '/', '\\', or "..".  Mirrors
// the validation in SVDSurfaceSerializer::default_cache_path; without
// it a caller passing an attacker-controlled fluid name through
// AbstractState::factory could read or write outside the cache dir,
// since `default_cache_dir()` derives the parent dir from
// `getenv("HOME")` (CodeQL tracks this as taint on fopen).
std::filesystem::path critpatch_cache_path_(const std::string& fluid_name, const std::string& source_backend, const std::string& options_canonical) {
    const auto unsafe = [](const std::string& s) {
        return s.empty() || s.find('/') != std::string::npos || s.find('\\') != std::string::npos || s.find("..") != std::string::npos;
    };
    if (unsafe(fluid_name)) {
        throw std::invalid_argument("SVDSBTLBackend: invalid fluid_name (must be a bare component name)");
    }
    if (unsafe(source_backend)) {
        throw std::invalid_argument("SVDSBTLBackend: invalid source_backend");
    }
    const std::filesystem::path dir = CoolProp::sbtl::SVDSurfaceSerializer::default_cache_dir();
    const std::uint64_t opthash = opthash_fnv1a_(options_canonical);
    // Match the SVDSurface filename pattern <Fluid>.<Source>.<...>.<opthash>.svd.bin.z
    // but with .critpatch.bin suffix (no zlib — too small to matter).
    std::ostringstream oss;
    oss << fluid_name << '.' << source_backend << ".critpatch." << std::hex << opthash << ".bin";
    return dir / oss.str();
}

}  // namespace

SVDSBTLBackend::SVDSBTLBackend(const std::string& fluid_name, const std::string& source_backend, const std::string& options_json)
  : fluid_name_(fluid_name), source_backend_(source_backend), mole_fractions_({1.0}), options_canonical_(parse_and_canonicalise(options_json)) {
    // Validate the source backend up-front.  Anything outside the
    // {HEOS, REFPROP, IF97} set is rejected — adding a new source
    // (SRK, PR, etc.) needs deliberate per-backend wiring for
    // two-phase routing and the cache-key encoding below, so failing
    // loudly is better than silently mismatching.
    if (source_backend_ != "HEOS" && source_backend_ != "REFPROP" && source_backend_ != "IF97") {
        throw ValueError(format("SVDSBTL: unsupported source backend '%s' (allowed: HEOS, REFPROP, IF97)", source_backend_.c_str()));
    }
    // IF97 is parameterized — the only fluid it knows is Water.  Catch
    // the obvious misuse here rather than at the first surface-build
    // call.
    if (source_backend_ == "IF97" && fluid_name != "Water") {
        throw ValueError(format("SVDSBTL&IF97 is Water-only; got fluid '%s'", fluid_name.c_str()));
    }
    for (const auto pair : kSupportedPairs) {
        ensure_surface_(pair);
    }
    build_critical_patch_(options_canonical_);
}

bool SVDSBTLBackend::available_in_high_level() {
    // Off by default — PropsSI rebuilds the AbstractState per call and
    // loading an SVDSurface from .svd.bin.z costs ~80 ms per
    // construction, vs ~5 us for the actual evaluation.  Throughput
    // workloads should use AbstractState directly + update() / batched
    // fast_evaluate.  Set ALLOW_SVDSBTL_IN_PROPSSI=true to opt back
    // in (e.g. for one-off interactive queries where the constructor
    // cost is irrelevant).
    return get_config_bool(ALLOW_SVDSBTL_IN_PROPSSI);
}

bool SVDSBTLBackend::CriticalPatch::contains(::CoolProp::input_pairs pair, double a, double b) const noexcept {
    if (!enabled) {
        return false;
    }
    const auto it = bbox_per_pair.find(static_cast<int>(pair));
    if (it == bbox_per_pair.end()) {
        return false;
    }
    const auto& bb = it->second;
    return (a >= bb.a_lo && a <= bb.a_hi && b >= bb.b_lo && b <= bb.b_hi);
}

std::shared_ptr<::CoolProp::AbstractState> SVDSBTLBackend::patch_source_ref_() {
    if (patch_source_) {
        return patch_source_;
    }
    // critical_patch_.source is empty when the options blob left it
    // null / default — fall back to the SVD truth source so the
    // patch IS the source by definition.
    const std::string backend = critical_patch_.source.empty() ? source_backend_ : critical_patch_.source;
    patch_source_.reset(::CoolProp::AbstractState::factory(backend, fluid_name_));
    return patch_source_;
}

void SVDSBTLBackend::build_critical_patch_(const std::string& options_canonical) {
    critical_patch_ = CriticalPatch{};  // disabled by default
    // Default behaviour when no options are supplied (mode unset) is
    // "auto", which now triggers the binary-search shrink loop in
    // auto_calibrate_critical_bbox_() (cached to a sidecar file so
    // subsequent constructions skip it).
    std::string mode = "auto";
    // Fallback multipliers: Water-sized.  Used when mode == "fixed"
    // without a bbox override, or when the auto-calibration fails
    // (source backend throws) and we can't load a cached result either.
    std::array<double, 4> bbox_mult = {0.95, 1.05, 0.75, 1.15};  // T_lo, T_hi, p_lo, p_hi
    bool user_supplied_bbox = false;
    std::string patch_source;

    if (!options_canonical.empty()) {
        rapidjson::Document opts;
        opts.Parse(options_canonical.c_str(), options_canonical.size());
        if (opts.IsObject() && opts.HasMember("critical_patch") && opts["critical_patch"].IsObject()) {
            const auto& cp = opts["critical_patch"];
            if (cp.HasMember("mode") && cp["mode"].IsString()) {
                mode = cp["mode"].GetString();
            }
            if (cp.HasMember("source") && cp["source"].IsString()) {
                patch_source = cp["source"].GetString();
            }
            if (cp.HasMember("bbox") && cp["bbox"].IsArray() && cp["bbox"].Size() == 4) {
                for (std::size_t i = 0; i < 4; ++i) {
                    bbox_mult[i] = cp["bbox"][static_cast<rapidjson::SizeType>(i)].GetDouble();
                }
                user_supplied_bbox = true;
            }
        }
    }
    if (mode == "off") {
        return;
    }
    if (mode != "auto" && mode != "fixed") {
        // schema validation should have caught this already
        return;
    }
    // Mode "auto" → try cache, then calibrate if missing.  User-supplied
    // bbox overrides auto-calibration even in "auto" mode (escape hatch
    // for users who want to pin the bbox without flipping mode to "fixed").
    if (mode == "auto" && !user_supplied_bbox) {
        const auto cached = load_critpatch_cache_(fluid_name_, source_backend_, options_canonical);
        if (cached) {
            bbox_mult = *cached;
        } else {
            // Need a patch_source for the calibration probe; ensure it's
            // set up before the calibrator pokes at SVD properties (the
            // calibration only reads from the SVD surface + source_, but
            // future calibrators may want patch_source_, so initialise
            // here to keep the calibrator simple).
            critical_patch_.source = patch_source;
            try {
                bbox_mult = auto_calibrate_critical_bbox_();
                save_critpatch_cache_(fluid_name_, source_backend_, options_canonical, bbox_mult);
            } catch (const std::exception&) {
                // Calibrator threw (source backend rejected a probe,
                // SVD surface missing a property, etc.) — fall back to
                // the Water-sized defaults rather than disabling the
                // patch entirely.  The defaults are conservative
                // (over-large) for non-water fluids, but better
                // accuracy than no patch.
            }
        }
    }

    // Set the patch source up-front so an invalid `critical_patch.source`
    // override fails at construction (factory throws on unknown backend)
    // instead of on the first in-patch query.  Also: the bbox must be
    // sampled from the SAME backend that will serve the patch — if a
    // user runs SVDSBTL&IF97 with critical_patch.source="HEOS", the
    // (p, h) envelope has to come from HEOS's h(T, p) so the in-patch
    // (h, p) classification matches what HEOS would return.
    critical_patch_.enabled = true;
    critical_patch_.source = patch_source;
    auto src = patch_source_ref_();  // factories the patch backend; throws on invalid

    const double Tc = src->T_critical();
    const double pc = src->p_critical();
    const double T_lo = bbox_mult[0] * Tc;
    const double T_hi = bbox_mult[1] * Tc;
    const double p_lo = bbox_mult[2] * pc;
    const double p_hi = bbox_mult[3] * pc;

    // PT_INPUTS — native bbox is (p, T): a = p, b = T.
    critical_patch_.bbox_per_pair[static_cast<int>(::CoolProp::PT_INPUTS)] = CriticalPatchBbox{p_lo, p_hi, T_lo, T_hi};

    // HmassP_INPUTS — bbox is (p, h).  Walk the (T, p) bbox's
    // perimeter through the PATCH backend's h(T, p) and take the
    // conservative axis-aligned envelope of the resulting h values.
    // Catches every (T, p) in the critical box; may include cells
    // slightly outside (false positive → source backend used outside
    // the strict box, correct but slightly slower; no false negatives).
    double h_lo = std::numeric_limits<double>::infinity();
    double h_hi = -std::numeric_limits<double>::infinity();
    constexpr int kPerimeterSamples = 24;
    auto try_h = [&](double T, double p) {
        try {
            src->update(::CoolProp::PT_INPUTS, p, T);
            const double h = src->hmass();
            if (std::isfinite(h)) {
                h_lo = std::min(h_lo, h);
                h_hi = std::max(h_hi, h);
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            // Source may reject (T, p) inside the dome at p just below
            // pc — those cells aren't critical-supercritical anyway.
            // Skipping them only narrows the h-bbox, which is the
            // safer direction.
        }
    };
    for (int i = 0; i <= kPerimeterSamples; ++i) {
        const double f = static_cast<double>(i) / kPerimeterSamples;
        try_h(T_lo + f * (T_hi - T_lo), p_lo);
        try_h(T_lo + f * (T_hi - T_lo), p_hi);
        try_h(T_lo, p_lo + f * (p_hi - p_lo));
        try_h(T_hi, p_lo + f * (p_hi - p_lo));
    }
    if (std::isfinite(h_lo) && std::isfinite(h_hi) && h_lo < h_hi) {
        critical_patch_.bbox_per_pair[static_cast<int>(::CoolProp::HmassP_INPUTS)] = CriticalPatchBbox{p_lo, p_hi, h_lo, h_hi};
    }
}

std::array<double, 4> SVDSBTLBackend::auto_calibrate_critical_bbox_() {
    // Binary-search-shrink the (T_lo, T_hi, p_lo, p_hi) bbox toward
    // (Tc, pc) until the strip JUST OUTSIDE the candidate patch can be
    // served by the SVD within IAPWS conformance budgets.  Each of the
    // four axes is shrunk independently; the resulting bbox is the
    // smallest axis-aligned rectangle that contains every cell where
    // the rank-r SVD's reconstruction error exceeds budget.
    //
    // Why "just outside" matters: probes INSIDE the candidate bbox
    // would be served by the patch (the source backend, by definition
    // exact), so testing them tells us nothing about SVD accuracy.
    // Probes just outside are the ones the SVD will own if we accept
    // this candidate, so those are the ones whose error budget we
    // need to respect.

    auto src = patch_source_ref_();
    const double Tc = src->T_critical();
    const double pc = src->p_critical();

    // Per-property relative-error budgets for the calibration probe.
    // INTENTIONALLY 50-100x wider than the IAPWS G13-15 perm values
    // (200 ppm rho, 100 ppm h, 200 ppm s, 1000 ppm w).  The patch
    // exists to cover the critical-singularity rank-truncation cusp
    // — where SVD reconstruction error reaches ~1%-10% in rho and
    // ~5%-50% in w.  Thin-support corner cells far from the critical
    // point also exceed the strict IAPWS budget by 10-100x (this is
    // CoolProp-3c4 accuracy ceiling territory), but they're not what
    // the patch is for — the patch can't cover scattered failures
    // with an axis-aligned bbox.  Using IAPWS-strict budgets here
    // would refuse to shrink the patch because every wide candidate
    // includes some corner outliers; relaxing to "1% relative" cleanly
    // separates the critical-cusp divergence from the corner ceiling.
    constexpr double kBudget_rho = 1.0e-2;  // 1%
    constexpr double kBudget_h = 1.0e-2;    // 1%
    constexpr double kBudget_s = 1.0e-2;    // 1%
    constexpr double kBudget_w = 5.0e-2;    // 5%

    // The PT surface is where the calibration probes go — it's the
    // canonical (T, p) → (rho, h, s, w) source the patch shape is
    // expressed in.  Bail to defaults if the surface isn't registered
    // (would indicate a non-standard preset; safer to use the
    // Water-sized fallback than to error out).
    const auto it = surfaces_.find(static_cast<int>(::CoolProp::PT_INPUTS));
    if (it == surfaces_.end() || !it->second) {
        throw ValueError("SVDSBTL auto-calibration: no PT_INPUTS surface registered");
    }
    const auto& pt_surface = *it->second;

    // Probe at (T, p): SVD vs source; return max relative error across
    // (rho, h, s, w), or NaN if either side rejects the point.
    auto probe_rel_err = [&](double T, double p) -> double {
        try {
            src->update(::CoolProp::PT_INPUTS, p, T);
            const double rho_src = src->rhomass();
            const double h_src = src->hmass();
            const double s_src = src->smass();
            double w_src = std::numeric_limits<double>::quiet_NaN();
            try {
                w_src = src->speed_sound();
            } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
                // Some sources reject speed_sound near the critical
                // point — skip and rely on (rho, h, s) alone there.
            }

            const auto resolved = pt_surface.resolve(p, T);
            if (resolved.region_idx < 0) {
                // SVD has no region containing this (p, T); not an
                // accuracy failure, just out-of-table.  Don't reject
                // the shrink based on it — the strip walker checks
                // multiple probes per axis position, and at least one
                // should land inside the SVD envelope.
                return std::numeric_limits<double>::quiet_NaN();
            }
            const double rho_svd = pt_surface.eval_with_region(::CoolProp::iDmass, resolved.region_idx, resolved.svd_x, resolved.svd_y);
            const double h_svd = pt_surface.eval_with_region(::CoolProp::iHmass, resolved.region_idx, resolved.svd_x, resolved.svd_y);
            const double s_svd = pt_surface.eval_with_region(::CoolProp::iSmass, resolved.region_idx, resolved.svd_x, resolved.svd_y);
            double w_svd = std::numeric_limits<double>::quiet_NaN();
            if (pt_surface.contains_property(::CoolProp::ispeed_sound)) {
                w_svd = pt_surface.eval_with_region(::CoolProp::ispeed_sound, resolved.region_idx, resolved.svd_x, resolved.svd_y);
            }

            double worst = 0.0;
            if (std::isfinite(rho_src) && std::isfinite(rho_svd) && std::abs(rho_src) > 0.0) {
                worst = std::max(worst, std::abs(rho_svd - rho_src) / std::abs(rho_src) / kBudget_rho);
            }
            if (std::isfinite(h_src) && std::isfinite(h_svd) && std::abs(h_src) > 0.0) {
                worst = std::max(worst, std::abs(h_svd - h_src) / std::abs(h_src) / kBudget_h);
            }
            if (std::isfinite(s_src) && std::isfinite(s_svd) && std::abs(s_src) > 0.0) {
                worst = std::max(worst, std::abs(s_svd - s_src) / std::abs(s_src) / kBudget_s);
            }
            if (std::isfinite(w_src) && std::isfinite(w_svd) && std::abs(w_src) > 0.0) {
                worst = std::max(worst, std::abs(w_svd - w_src) / std::abs(w_src) / kBudget_w);
            }
            // worst is in units of "fraction of the largest budget";
            // <= 1 means all properties pass; > 1 means at least one
            // property exceeds its budget.
            return worst;
        } catch (const std::exception&) {
            // Source rejected the probe — return NaN so the caller
            // treats this probe as a no-op (doesn't influence the
            // shrink decision either way).
            return std::numeric_limits<double>::quiet_NaN();
        }
    };

    // Strip walker: at a candidate axis multiplier, sample N probes
    // along the strip JUST OUTSIDE that axis (at axis-position ± eps,
    // varying over the other two axes' wide span).  Returns true iff
    // at least N - kStripOutlierBudget finite probes pass their budget
    // — isolated outliers (typically thin-Hermite-support cells at
    // region edges, e.g. just-below-dome cells near critical) don't
    // block the calibrator from shrinking past them.  The patch is
    // for the critical-singularity divergence, not the
    // CoolProp-3c4 corner-cell accuracy ceiling; the outlier budget
    // cleanly separates the two.
    constexpr int kStripSamples = 12;
    constexpr int kStripOutlierBudget = 2;  // tolerate up to 2/12 = 17% outliers
    constexpr double kStripEps = 0.005;     // 0.5% offset outside the candidate
    auto strip_passes = [&](double T_lo_mult, double T_hi_mult, double p_lo_mult, double p_hi_mult, int axis) {
        // axis: 0=T_lo, 1=T_hi, 2=p_lo, 3=p_hi
        const double T_lo = T_lo_mult * Tc;
        const double T_hi = T_hi_mult * Tc;
        const double p_lo = p_lo_mult * pc;
        const double p_hi = p_hi_mult * pc;
        // Outer reference for the variation axis: matches the
        // Water-default bbox so the strip walker samples the region
        // that the default patch would have covered.  Going wider
        // would sample cells well outside the default patch where
        // the SVD is the only owner anyway — those are CoolProp-3c4
        // accuracy-ceiling territory and shouldn't influence the
        // patch shape.
        constexpr double kT_outer_lo = 0.95;
        constexpr double kT_outer_hi = 1.05;
        constexpr double kp_outer_lo = 0.75;
        constexpr double kp_outer_hi = 1.15;
        int n_finite = 0;
        int n_outliers = 0;
        for (int i = 0; i < kStripSamples; ++i) {
            const double f = (i + 0.5) / kStripSamples;
            double T_probe = T_lo, p_probe = p_lo;
            if (axis == 0) {
                T_probe = T_lo * (1.0 - kStripEps);  // just below T_lo
                p_probe = (kp_outer_lo + f * (kp_outer_hi - kp_outer_lo)) * pc;
            } else if (axis == 1) {
                T_probe = T_hi * (1.0 + kStripEps);  // just above T_hi
                p_probe = (kp_outer_lo + f * (kp_outer_hi - kp_outer_lo)) * pc;
            } else if (axis == 2) {
                p_probe = p_lo * (1.0 - kStripEps);  // just below p_lo
                T_probe = (kT_outer_lo + f * (kT_outer_hi - kT_outer_lo)) * Tc;
            } else {
                p_probe = p_hi * (1.0 + kStripEps);  // just above p_hi
                T_probe = (kT_outer_lo + f * (kT_outer_hi - kT_outer_lo)) * Tc;
            }
            const double err = probe_rel_err(T_probe, p_probe);
            if (std::isfinite(err)) {
                ++n_finite;
                if (err > 1.0) {
                    ++n_outliers;
                }
            }
        }
        // Require: at least one finite probe (otherwise we're shrinking
        // past a region the source rejects entirely), AND fewer than
        // kStripOutlierBudget budget-violating probes.
        return n_finite > 0 && n_outliers <= kStripOutlierBudget;
    };

    // Binary search per axis.  Strategy: start at the Water-sized
    // defaults (which are known to give G13-15 conformant SVDSBTL&IF97
    // for water) and only SHRINK toward 1.0 (= the critical point).
    // Never widen — that way Water gets back its tested bbox, and
    // other fluids tighten if their critical regions are less
    // extended.
    //
    // Axis indices: 0=T_lo, 1=T_hi, 2=p_lo, 3=p_hi.
    // "safer" (wider patch) is the default; "tighter" is 1.0.
    std::array<double, 4> defaults = {0.95, 1.05, 0.75, 1.15};
    std::array<double, 4> tight = {1.0, 1.0, 1.0, 1.0};
    std::array<double, 4> mults = defaults;
    constexpr int kBisectSteps = 6;  // 6 = 1/64 multiplier resolution
    for (int axis = 0; axis < 4; ++axis) {
        // First confirm the default passes the strip test — if it
        // doesn't (e.g. a fluid where the default patch is itself too
        // small), bail to defaults and let the user override via
        // critical_patch.bbox.  No widening here.
        if (!strip_passes(mults[0], mults[1], mults[2], mults[3], axis)) {
            continue;
        }
        // The (safe, aggressive) pair encodes the shrink direction per
        // axis: safe = defaults[axis] (Water-sized, known to pass);
        // aggressive = 1.0 (smallest patch, may or may not pass).  Lo-
        // side axes (0=T_lo, 2=p_lo) have safe < aggressive numerically;
        // hi-side axes (1=T_hi, 3=p_hi) have safe > aggressive.  Either
        // way, the bisection moves `safe` toward `aggressive` whenever
        // the candidate passes, and pulls `aggressive` toward `safe`
        // when it doesn't.
        double safe = defaults[axis];
        double aggressive = tight[axis];
        for (int step = 0; step < kBisectSteps; ++step) {
            const double mid = 0.5 * (safe + aggressive);
            auto trial = mults;
            trial[axis] = mid;
            if (strip_passes(trial[0], trial[1], trial[2], trial[3], axis)) {
                safe = mid;
            } else {
                aggressive = mid;
            }
        }
        // safe is the latest known-passing value — commit it.
        mults[axis] = safe;
    }
    return mults;
}

std::optional<std::array<double, 4>> SVDSBTLBackend::load_critpatch_cache_(const std::string& fluid_name, const std::string& source_backend,
                                                                           const std::string& options_canonical) {
    const auto path = critpatch_cache_path_(fluid_name, source_backend, options_canonical);
    std::error_code ec;
    if (!std::filesystem::exists(path, ec) || ec) {
        return std::nullopt;
    }
    // Use FILE* (POSIX/Windows-portable; binary fixed-width 44 bytes is
    // small enough that stdio buffering overhead is irrelevant).
    std::FILE* f = std::fopen(path.string().c_str(), "rb");
    if (!f) {
        return std::nullopt;
    }
    std::array<char, 8> magic{};
    std::uint32_t version = 0;
    std::array<double, 4> mults{};
    const auto rd_n = std::fread(magic.data(), 1, 8, f);
    if (rd_n != 8 || std::string(magic.data(), 8) != "CPCRITPB") {
        (void)std::fclose(f);
        return std::nullopt;
    }
    if (std::fread(&version, sizeof(version), 1, f) != 1 || version != 1u) {
        (void)std::fclose(f);
        return std::nullopt;
    }
    if (std::fread(mults.data(), sizeof(double), 4, f) != 4) {
        (void)std::fclose(f);
        return std::nullopt;
    }
    (void)std::fclose(f);
    // Sanity-clamp: refuse obviously-wrong cached values rather than
    // letting a corrupt file silently produce a broken patch.
    if (!(mults[0] > 0.0 && mults[0] < mults[1] && mults[1] < 2.0 && mults[2] > 0.0 && mults[2] < mults[3] && mults[3] < 10.0)) {
        return std::nullopt;
    }
    return mults;
}

void SVDSBTLBackend::save_critpatch_cache_(const std::string& fluid_name, const std::string& source_backend, const std::string& options_canonical,
                                           const std::array<double, 4>& mults) {
    const auto path = critpatch_cache_path_(fluid_name, source_backend, options_canonical);
    std::error_code ec;
    std::filesystem::create_directories(path.parent_path(), ec);
    if (ec) {
        // Cache directory inaccessible — non-fatal (the calibration
        // result is still usable in-process; next construction just
        // re-runs the calibrator).
        return;
    }
    constexpr std::array<char, 8> kMagic = {'C', 'P', 'C', 'R', 'I', 'T', 'P', 'B'};
    constexpr std::uint32_t kVersion = 1;
    constexpr std::size_t kPayloadSize = kMagic.size() + sizeof(kVersion) + 4 * sizeof(double);
    std::array<char, kPayloadSize> buf{};
    std::memcpy(buf.data(), kMagic.data(), kMagic.size());
    std::memcpy(buf.data() + kMagic.size(), &kVersion, sizeof(kVersion));
    std::memcpy(buf.data() + kMagic.size() + sizeof(kVersion), mults.data(), 4 * sizeof(double));
    // Atomic write: temp + rename, with owner-only perms applied to the
    // temp file before rename so multi-process notebook builds never
    // see a torn-write critpatch sidecar (CoolProp-4no.2).  Any failure
    // is non-fatal — same semantics as the previous fopen-then-fwrite
    // path: cache miss on next load just re-runs the calibrator.
    try {
        ::write_bytes_atomic(path, buf.data(), buf.size(), /*restrict_perms=*/true);
    } catch (const std::exception&) {
        // swallowed — cache write is best-effort
    }
}

std::shared_ptr<CoolProp::AbstractState> SVDSBTLBackend::source_reference_() {
    if (!source_) {
        // The single AbstractState handle that all internal
        // truth-source calls go through: HEOS sampling for the SVD
        // tables, SuperAncillary (when source is HEOS), critical-
        // patch fallback queries, and two-phase PQ/QT routing.
        //
        // IF97's factory takes the literal "Water" (its only fluid);
        // the constructor already validated fluid_name_=="Water" for
        // this source.  REFPROP and HEOS use the same fluid string.
        source_.reset(CoolProp::AbstractState::factory(source_backend_, fluid_name_));
    }
    return source_;
}

std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> SVDSBTLBackend::superanc_() {
    if (superanc_resolved_) {
        return superanc_cached_;
    }
    superanc_resolved_ = true;
    // The SuperAncillary lives on the per-fluid CoolPropFluid record;
    // HEOSMixtureBackend wraps it via get_superanc().  Anything that's
    // not a HEOS instance simply won't have one.
    auto* helmholtz_ptr = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(source_reference_().get());
    if (helmholtz_ptr == nullptr) {
        return superanc_cached_;  // stays nullptr
    }
    try {
        superanc_cached_ = helmholtz_ptr->get_superanc();
    } catch (const std::exception&) {
        // get_superanc() throws for pseudo-pure / mixtures; SVDSBTL
        // already rejects those at construction, so this is mostly
        // defensive.  Leave superanc_cached_ as nullptr.
        return superanc_cached_;
    }
    if (superanc_cached_) {
        // Caloric expansions (h_sat, s_sat, u_sat) are lazily built --
        // trigger them eagerly here so the first two-phase calc_*
        // doesn't pay the build cost on the hot path.  Cheap (~30 ms
        // total) and runs once per backend lifetime.
        helmholtz_ptr->ensure_caloric_superancillaries();
    }
    return superanc_cached_;
}

void SVDSBTLBackend::set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
    if (mole_fractions.size() > 1) {
        throw NotImplementedError("SVDSBTL backend is pure-fluid only; mixtures are not supported");
    }
    if (mole_fractions.size() == 1) {
        mole_fractions_ = mole_fractions;
    }
}

std::vector<CoolProp::input_pairs> SVDSBTLBackend::registered_input_pairs() const {
    std::vector<CoolProp::input_pairs> out;
    out.reserve(surfaces_.size());
    for (const auto& kv : surfaces_) {
        out.push_back(static_cast<CoolProp::input_pairs>(kv.first));
    }
    return out;
}

void SVDSBTLBackend::ensure_surface_(CoolProp::input_pairs pair) {
    namespace cp_sbtl = CoolProp::sbtl;
    // Cache filename: <fluid>.<source>.<input_pair_name>.<opthash>.svd.bin.z
    //  - input_pair_name (e.g. "PT_INPUTS", "HmassP_INPUTS") instead
    //    of the raw enum int, so reordering DataStructures.h's
    //    `input_pairs` enum can no longer silently misalign caches
    //    (closes CoolProp-b6v).
    //  - opthash: 16-hex FNV-1a 64 prefix over the canonical-JSON
    //    options blob; two different option sets get different cache
    //    files automatically (canonical "{}" is treated as the all-
    //    defaults case so the no-options path still picks up its
    //    intended cache file).
    const std::string opthash = to_hex16(fnv1a_64(options_canonical_.empty() ? std::string{"{}"} : options_canonical_));
    const std::string path = cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid_name_, source_backend_, pair, opthash);
    std::unique_ptr<cp_sbtl::SVDSurface> surface;

    if (std::filesystem::exists(path)) {
        try {
            surface = std::make_unique<cp_sbtl::SVDSurface>(cp_sbtl::SVDSurfaceSerializer::load_from_file(path));
        } catch (const std::exception&) {
            // Corrupted cache — fall through to rebuild.  Logging would
            // be nice here but TabularBackend doesn't log either.
            surface.reset();
        }
    }

    if (!surface) {
        auto source = source_reference_();
        rapidjson::Document opts;
        opts.Parse(options_canonical_.empty() ? "{}" : options_canonical_.c_str());
        const ResolvedGrid grid = resolve_grid(opts);
        surface = std::make_unique<cp_sbtl::SVDSurface>(build_surface_for_pair_(*source, source_backend_, pair, grid));
        try {
            cp_sbtl::SVDSurfaceSerializer::save_to_file(*surface, path);
        } catch (const std::exception&) {
            // Cache write failure is non-fatal — the in-memory surface
            // still works for this run.
        }
    }

    surfaces_[static_cast<int>(pair)] = std::move(surface);
}

// resolve_point_: the shared kernel called by both update() and
// fast_evaluate().  Maps (input_pair, v1, v2) to a PointEvaluation
// without mutating any backend state.  The same routing logic — input-
// pair dispatch, critical-patch bbox check, atlas resolve, dome
// detection with sat-line lever-rule — lives in one place so any
// future entry point (batched PropsSI, Python update_array wrapper,
// ...) plugs in by reusing this method.
//
// Two-phase PQ_INPUTS / QT_INPUTS are folded into the DomeBlend kind
// — same sat-line endpoints, same Q, same downstream property dispatch.
SVDSBTLBackend::PointEvaluation SVDSBTLBackend::resolve_point_(CoolProp::input_pairs input_pair, double value1, double value2) {
    PointEvaluation pt;
    pt.input_pair = input_pair;
    pt.v1 = value1;
    pt.v2 = value2;

    // --- PQ / QT inputs: dome by definition, no surface needed. ---
    if (input_pair == CoolProp::PQ_INPUTS || input_pair == CoolProp::QT_INPUTS) {
        double p_sat = 0.0;
        double T_sat = 0.0;
        double Q = 0.0;
        if (input_pair == CoolProp::PQ_INPUTS) {
            p_sat = value1;
            Q = value2;
        } else {
            Q = value1;
            T_sat = value2;
        }
        if (Q < 0.0 || Q > 1.0) {
            throw ValueError("SVDSBTL backend: two-phase Q must be in [0, 1]");
        }
        auto sa = superanc_();
        if (sa) {
            if (input_pair == CoolProp::PQ_INPUTS) {
                T_sat = sa->get_T_from_p(p_sat);
            } else {
                p_sat = sa->eval_sat(T_sat, 'P', 0);
            }
            pt.hL_mol = sa->eval_sat(T_sat, 'H', 0);
            pt.hV_mol = sa->eval_sat(T_sat, 'H', 1);
        } else {
            // Source-backend PQ/QT fallback (REFPROP / IF97).
            auto src = source_reference_();
            if (input_pair == CoolProp::PQ_INPUTS) {
                src->update(CoolProp::PQ_INPUTS, p_sat, 0.0);
                T_sat = src->T();
                pt.hL_mol = src->hmolar();
                src->update(CoolProp::PQ_INPUTS, p_sat, 1.0);
                pt.hV_mol = src->hmolar();
            } else {
                src->update(CoolProp::QT_INPUTS, 0.0, T_sat);
                p_sat = src->p();
                pt.hL_mol = src->hmolar();
                src->update(CoolProp::QT_INPUTS, 1.0, T_sat);
                pt.hV_mol = src->hmolar();
            }
        }
        pt.kind = PointEvaluation::Kind::DomeBlend;
        pt.status = CoolProp::fast_evaluate_ok;
        pt.T = T_sat;
        pt.T_sat = T_sat;
        pt.p = p_sat;
        pt.Q = Q;
        const double M = calc_molar_mass();
        pt.h_mass = ((1.0 - Q) * *pt.hL_mol + Q * *pt.hV_mol) / M;
        return pt;
    }

    // --- Single-phase / dome routing through the surface atlas. ---
    const auto surface_key = (input_pair == CoolProp::HmolarP_INPUTS) ? CoolProp::HmassP_INPUTS : input_pair;
    const auto it = surfaces_.find(static_cast<int>(surface_key));
    if (it == surfaces_.end() || !it->second) {
        throw ValueError(std::string("SVDSBTL backend: no surface registered for this input pair (") + CoolProp::get_input_pair_short_desc(input_pair)
                         + ")");
    }
    const auto* surface = it->second.get();
    pt.surface = surface;

    // Extract (a, b) primary/secondary axes and the canonical h_mass
    // input shortcut where applicable.
    double a = 0.0;  // primary axis: always p for the MVP presets
    double b = 0.0;  // secondary axis: h for HmassP/HmolarP, T for PT
    switch (input_pair) {
        case CoolProp::HmassP_INPUTS:
            a = value2;
            pt.p = value2;
            pt.h_mass = value1;
            b = value1;
            break;
        case CoolProp::HmolarP_INPUTS:
            a = value2;
            pt.p = value2;
            pt.h_mass = value1 / calc_molar_mass();
            b = *pt.h_mass;
            break;
        case CoolProp::PT_INPUTS:
            a = value1;
            pt.p = value1;
            pt.T = value2;
            b = value2;
            break;
        default:
            throw ValueError(std::string("SVDSBTL backend update: unsupported input pair (") + CoolProp::get_input_pair_short_desc(input_pair) + ")");
    }

    // Critical-patch bbox routing.  Patch lookups defer their property
    // reads to evaluate_property_ via patch_source_.
    //
    // For HmassP / HmolarP inputs we follow IF97's source-backend
    // update with polish_patch_state_(): IF97's HmassP_INPUTS uses the
    // R7-97 backward equation, which is only ±25 mK forward-consistent
    // in T.  Without the polish every in-bbox query inherits that floor
    // — visible as ~10-20 mK / 0.1-0.4 J/(kg·K) s residuals in the
    // SVDSBTL&IF97 fail-map even inside the patch box.  The polish
    // TOMS748-iterates on a tight bracket around T_seed so the returned
    // state is forward-consistent (h(T_polished, p) == h_input to
    // bracket eps).  Mirrors the polish baked into ph_subcritical's
    // IF97 sampling lambda — that's the foi.9.10 sampling-side fix;
    // this is its query-side counterpart.
    //
    // Polish is gated to IF97 sources only.  HEOS's HmassP_INPUTS is
    // already iterative + forward-consistent, so the polish should be
    // a no-op — but TOMS748 in [T_seed±0.5 K] leaves a sub-ULP T
    // residual that, for low-Tc fluids (Hydrogen Tc=33.14 K — the
    // bracket is ~1.5% of Tc), the critical-region stiffness amplifies
    // into ~1e-7 ρ drift inside the patch bbox.  Visible as non-zero
    // error in 82% of bbox cells in the SVDSBTLValidation heatmap even
    // though the patch source IS the reference HEOS.  Skipping the
    // polish for HEOS also drops one extra HEOS PT_INPUTS eval per
    // in-bbox query.  See the [hydrogen] regression test in
    // CoolProp-Tests-SVDSBTLCriticalPatch.cpp.
    // Gate on the ACTIVE patch backend, not source_backend_: when the
    // user overrides critical_patch.source (e.g. SVDSBTL&HEOS with
    // patch source="IF97" for Water), the polish requirement follows
    // the backend actually serving the in-patch query.
    const std::string& patch_backend = critical_patch_.source.empty() ? source_backend_ : critical_patch_.source;
    const bool needs_polish = (patch_backend == "IF97");
    const auto patch_key = (input_pair == CoolProp::HmolarP_INPUTS) ? CoolProp::HmassP_INPUTS : input_pair;
    if (critical_patch_.contains(patch_key, a, b)) {
        try {
            if (input_pair == CoolProp::HmolarP_INPUTS) {
                patch_source_ref_()->update(CoolProp::HmassP_INPUTS, *pt.h_mass, *pt.p);
                // Same dome-aware polish gating as the HmassP branch
                // below — multi-valued h(T, p) inside the dome breaks
                // the single-phase forward-h polish bracket.
                if (needs_polish && patch_source_->phase() != CoolProp::iphase_twophase) {
                    polish_patch_state_(*patch_source_, *pt.p, *pt.h_mass);
                }
            } else if (input_pair == CoolProp::HmassP_INPUTS) {
                patch_source_ref_()->update(input_pair, value1, value2);
                // Skip the forward-h polish on two-phase points.
                // polish_patch_state_ assumes a single-phase solution
                // exists at h(T, p) = h_target on a tight T bracket
                // around T_seed.  Inside the saturation dome h(T, p)
                // is multi-valued (the function jumps from h_L to h_V
                // across T = T_sat), so the bracket [T_seed-0.5,
                // T_seed+0.5] straddles the dome and resid changes
                // sign without any single-phase root in the range.
                // TOMS748 silently converges to a nonsense T and the
                // returned state classifies as supercritical_liquid
                // with rho / Q = ±inf.  Two-phase points coming out
                // of the patch already have correct (T_sat, Q, rho)
                // from the source's PQ-blend; no polish needed.
                if (needs_polish && patch_source_->phase() != CoolProp::iphase_twophase) {
                    polish_patch_state_(*patch_source_, *pt.p, *pt.h_mass);
                }
            } else {
                // PT_INPUTS — source backend's PT path is already
                // forward-consistent; no polish needed.
                patch_source_ref_()->update(input_pair, value1, value2);
            }
            pt.kind = PointEvaluation::Kind::Patched;
            pt.status = CoolProp::fast_evaluate_ok;
            if (!pt.T) pt.T = patch_source_->T();
            if (!pt.p) pt.p = patch_source_->p();
            return pt;
        } catch (const std::exception&) {
            pt.kind = PointEvaluation::Kind::OutOfRange;
            pt.status = CoolProp::fast_evaluate_out_of_range;
            return pt;
        }
    }

    // Atlas resolve.  region_idx < 0 means the point landed inside
    // the saturation dome (for HmassP) or outside every region's bbox
    // (for PT or genuinely OOB HmassP).
    pt.resolved = surface->resolve(a, b);
    if (pt.resolved.region_idx >= 0) {
        pt.kind = PointEvaluation::Kind::SinglePhase;
        pt.status = CoolProp::fast_evaluate_ok;
        return pt;
    }

    // --- Atlas miss.  For HmassP / HmolarP try dome blend. ---
    if (input_pair == CoolProp::HmassP_INPUTS || input_pair == CoolProp::HmolarP_INPUTS) {
        const double p_val = *pt.p;
        const double h_mass = *pt.h_mass;
        const double M = calc_molar_mass();
        bool sat_resolved = false;
        double T_sat = 0.0, hL_mol = 0.0, hV_mol = 0.0;
        // Optional sat-line caloric endpoints we may also pick up from
        // the source PQ fallback (saves the redundant lookup later).
        std::optional<double> rhoL_mol_seed, rhoV_mol_seed, sL_mol_seed, sV_mol_seed;
        auto sa = superanc_();
        if (sa) {
            try {
                T_sat = sa->get_T_from_p(p_val);
                hL_mol = sa->eval_sat(T_sat, 'H', 0);
                hV_mol = sa->eval_sat(T_sat, 'H', 1);
                sat_resolved = true;
            } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
            }
        }
        if (!sat_resolved) {
            try {
                auto src = source_reference_();
                src->update(CoolProp::PQ_INPUTS, p_val, 0.0);
                T_sat = src->T();
                hL_mol = src->hmolar();
                rhoL_mol_seed = src->rhomolar();
                sL_mol_seed = src->smolar();
                src->update(CoolProp::PQ_INPUTS, p_val, 1.0);
                hV_mol = src->hmolar();
                rhoV_mol_seed = src->rhomolar();
                sV_mol_seed = src->smolar();
                sat_resolved = true;
            } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
            }
        }
        if (sat_resolved) {
            const double h_mol = h_mass * M;
            if (hL_mol <= h_mol && h_mol <= hV_mol) {
                pt.kind = PointEvaluation::Kind::DomeBlend;
                pt.status = CoolProp::fast_evaluate_ok;
                pt.T = T_sat;
                pt.T_sat = T_sat;
                pt.Q = (h_mol - hL_mol) / (hV_mol - hL_mol);
                pt.hL_mol = hL_mol;
                pt.hV_mol = hV_mol;
                pt.rhoL_mol = rhoL_mol_seed;
                pt.rhoV_mol = rhoV_mol_seed;
                pt.sL_mol = sL_mol_seed;
                pt.sV_mol = sV_mol_seed;
                return pt;
            }
        }
    }

    // --- PT-on-sat detection vs genuine OOB. ---
    if (input_pair == CoolProp::PT_INPUTS) {
        try {
            auto sa = superanc_();
            double T_sat_probe = std::numeric_limits<double>::quiet_NaN();
            if (sa) {
                T_sat_probe = sa->get_T_from_p(a);
            } else {
                auto src = source_reference_();
                src->update(CoolProp::PQ_INPUTS, a, 0.0);
                T_sat_probe = src->T();
            }
            constexpr double kSatEps = 3.3e-5;  // matches IF97Backend::set_phase
            if (std::isfinite(T_sat_probe) && std::abs(b - T_sat_probe) <= kSatEps * T_sat_probe) {
                pt.kind = PointEvaluation::Kind::TwoPhaseDisallowed;
                pt.status = CoolProp::fast_evaluate_two_phase_disallowed;
                return pt;
            }
        } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
            // No sat probe available — fall through to OOB.
        }
    }
    pt.kind = PointEvaluation::Kind::OutOfRange;
    pt.status = CoolProp::fast_evaluate_out_of_range;
    return pt;
}

void SVDSBTLBackend::ensure_dome_rho_endpoints_(PointEvaluation& pt) {
    if (pt.rhoL_mol && pt.rhoV_mol) return;
    auto sa = superanc_();
    if (sa) {
        pt.rhoL_mol = sa->eval_sat(*pt.T_sat, 'D', 0);
        pt.rhoV_mol = sa->eval_sat(*pt.T_sat, 'D', 1);
        return;
    }
    auto src = source_reference_();
    src->update(CoolProp::QT_INPUTS, 0.0, *pt.T_sat);
    pt.rhoL_mol = src->rhomolar();
    src->update(CoolProp::QT_INPUTS, 1.0, *pt.T_sat);
    pt.rhoV_mol = src->rhomolar();
}

void SVDSBTLBackend::ensure_dome_s_endpoints_(PointEvaluation& pt) {
    if (pt.sL_mol && pt.sV_mol) return;
    auto sa = superanc_();
    if (sa) {
        pt.sL_mol = sa->eval_sat(*pt.T_sat, 'S', 0);
        pt.sV_mol = sa->eval_sat(*pt.T_sat, 'S', 1);
        return;
    }
    auto src = source_reference_();
    src->update(CoolProp::QT_INPUTS, 0.0, *pt.T_sat);
    pt.sL_mol = src->smolar();
    src->update(CoolProp::QT_INPUTS, 1.0, *pt.T_sat);
    pt.sV_mol = src->smolar();
}

// evaluate_property_: the per-output kernel.  Pure dispatch over the
// PointEvaluation kind — no business logic beyond input-pair
// shortcuts, surface eval routing, dome Q-blend, and patch deref.
//
// Throws NotImplementedError for properties that are mathematically
// undefined in the two-phase region (speed_sound, cp, cv, viscosity,
// conductivity, Prandtl); fast_evaluate translates those into in-band
// NaN values, while calc_*() lets them propagate to the caller.
double SVDSBTLBackend::evaluate_property_(PointEvaluation& pt, CoolProp::parameters prop) {
    const double M = calc_molar_mass();

    // Always-available shortcuts directly from the resolved point.
    // The conditional checks matter — pt.T isn't filled for HmassP /
    // HmolarP single-phase points until we read it off the surface, so
    // we only short-circuit when the value was either supplied by the
    // input pair (PT case), set by the dome blend (T = T_sat), or
    // pulled from patch_source_->T() / ->p() in resolve_point_.
    //
    // iQ is *not* short-circuited at top level — Patched probes
    // delegate to patch_source_'s Q (which may itself be two-phase
    // near the bbox edge), and OOB / TwoPhaseDisallowed / Invalid
    // kinds return NaN through the kind switch.  Only SinglePhase
    // and DomeBlend use the pt.Q default / blend value.
    if (prop == CoolProp::iP && pt.p) {
        return *pt.p;
    }
    if (prop == CoolProp::iT && pt.T) {
        return *pt.T;
    }
    if (prop == CoolProp::iHmass && pt.h_mass) {
        return *pt.h_mass;
    }
    if (prop == CoolProp::iHmolar && pt.h_mass) {
        // pt.h_mass is the single source of truth: PQ/QT and HmassP-
        // in-dome both populate it from the lever rule at resolve
        // time, so this also handles the DomeBlend case correctly.
        return *pt.h_mass * M;
    }

    switch (pt.kind) {
        case PointEvaluation::Kind::SinglePhase: {
            // Single-phase by construction (atlas resolved to a
            // region, which excludes the dome).  Mirror the -1
            // sentinel update()'s legacy contract uses for iQ.
            if (prop == CoolProp::iQ) {
                return -1.0;
            }
            // Map molar variants to the mass-basis property the surface
            // tabulates, then scale.
            CoolProp::parameters surface_prop = prop;
            bool need_div_M = false, need_mul_M = false;
            switch (prop) {
                case CoolProp::iDmolar:
                    surface_prop = CoolProp::iDmass;
                    need_div_M = true;
                    break;
                case CoolProp::iSmolar:
                    surface_prop = CoolProp::iSmass;
                    need_mul_M = true;
                    break;
                case CoolProp::iUmolar:
                    surface_prop = CoolProp::iUmass;
                    need_mul_M = true;
                    break;
                case CoolProp::iHmolar:
                    surface_prop = CoolProp::iHmass;
                    need_mul_M = true;
                    break;
                default:
                    break;
            }
            if (!pt.surface->contains_property(surface_prop)) {
                throw ValueError(std::string("SVDSBTL backend: active surface does not expose property '")
                                 + CoolProp::get_parameter_information(surface_prop, "short") + "'");
            }
            const double v = pt.surface->eval_with_region(surface_prop, pt.resolved.region_idx, pt.resolved.svd_x, pt.resolved.svd_y);
            if (need_div_M) return v / M;
            if (need_mul_M) return v * M;
            return v;
        }
        case PointEvaluation::Kind::DomeBlend: {
            const double Q = pt.Q;
            switch (prop) {
                case CoolProp::iQ:
                    return Q;
                case CoolProp::iDmass:
                case CoolProp::iDmolar: {
                    ensure_dome_rho_endpoints_(pt);
                    const double vL = 1.0 / *pt.rhoL_mol;
                    const double vV = 1.0 / *pt.rhoV_mol;
                    const double rho_molar = 1.0 / ((1.0 - Q) * vL + Q * vV);
                    return (prop == CoolProp::iDmass) ? rho_molar * M : rho_molar;
                }
                case CoolProp::iSmass:
                case CoolProp::iSmolar: {
                    ensure_dome_s_endpoints_(pt);
                    const double s_molar = (1.0 - Q) * *pt.sL_mol + Q * *pt.sV_mol;
                    return (prop == CoolProp::iSmass) ? s_molar / M : s_molar;
                }
                case CoolProp::iUmass:
                case CoolProp::iUmolar: {
                    ensure_dome_rho_endpoints_(pt);
                    // u = h - p/rho per the thermodynamic identity.
                    const double uL_mol = *pt.hL_mol - *pt.p / *pt.rhoL_mol;
                    const double uV_mol = *pt.hV_mol - *pt.p / *pt.rhoV_mol;
                    const double u_molar = (1.0 - Q) * uL_mol + Q * uV_mol;
                    return (prop == CoolProp::iUmass) ? u_molar / M : u_molar;
                }
                // iHmass / iHmolar are served by the lever-rule short-circuit
                // above (pt.h_mass is populated in resolve_point_); no
                // duplicate arms here.
                case CoolProp::ispeed_sound:
                    throw NotImplementedError("SVDSBTL backend: speed_sound is not defined in the two-phase region");
                case CoolProp::iviscosity:
                    throw NotImplementedError("SVDSBTL backend: viscosity is not defined in the two-phase region");
                case CoolProp::iconductivity:
                    throw NotImplementedError("SVDSBTL backend: conductivity is not defined in the two-phase region");
                default:
                    throw ValueError(std::string("SVDSBTL backend: two-phase blend not implemented for property '")
                                     + CoolProp::get_parameter_information(prop, "short") + "'");
            }
        }
        case PointEvaluation::Kind::Patched: {
            // Inside the critical-patch bbox the source backend owns
            // every answer — delegate without trying to second-guess
            // it via mass-basis round-trips.  iQ in particular can be
            // genuinely two-phase near the bbox edge if the patch was
            // sized to include sat-curve excursions.
            return patch_source_->keyed_output(prop);
        }
        case PointEvaluation::Kind::OutOfRange:
        case PointEvaluation::Kind::TwoPhaseDisallowed:
        case PointEvaluation::Kind::Invalid:
        case PointEvaluation::Kind::InternalError:
            return std::nan("");
    }
    return std::nan("");
}

void SVDSBTLBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    clear();
    _phase = iphase_not_imposed;
    active_eval_ = resolve_point_(input_pair, value1, value2);

    // Copy the resolved (T, p, Q) into AbstractState's standard slots
    // so legacy callers reading _T / _p / _Q / _phase from the base
    // class still see consistent values.
    if (active_eval_.T) _T = *active_eval_.T;
    if (active_eval_.p) _p = *active_eval_.p;
    // _Q is populated for every kind that has a well-defined Q so
    // legacy callers reading state->Q() see a value consistent with
    // what fast_evaluate(iQ) would return.  SinglePhase gets the -1
    // sentinel (matches IF97Backend::update()); DomeBlend gets the
    // lever-rule fraction; Patched delegates to the source backend
    // (which can be two-phase if the bbox extends across the sat
    // curve); OOB / TwoPhaseDisallowed / Invalid leave _Q untouched
    // (clear() set it to -_HUGE) since there's no defined value.
    switch (active_eval_.kind) {
        case PointEvaluation::Kind::SinglePhase:
            _Q = -1.0;
            break;
        case PointEvaluation::Kind::DomeBlend:
            _Q = active_eval_.Q;
            _phase = iphase_twophase;
            break;
        case PointEvaluation::Kind::Patched:
            try {
                _Q = patch_source_->keyed_output(CoolProp::iQ);
            } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
                _Q = -1.0;
            }
            // Inherit the source's phase classification.  Without this
            // a two-phase patched state (which the skip-polish path now
            // correctly preserves) would surface as iphase_not_imposed
            // and state.phase() would disagree with the source.  Wrap
            // in try/catch because some sources may not report a phase
            // for every input pair (HEOS does; defensive against
            // REFPROP edge cases).
            try {
                _phase = patch_source_->phase();
            } catch (const std::exception&) {  // NOLINT(bugprone-empty-catch)
            }
            break;
        case PointEvaluation::Kind::OutOfRange:
            // Legacy back-compat: when an HmassP probe misses
            // everything (above critical pressure, etc.) we used to
            // tag iphase_twophase as a best-effort label.  Preserve
            // that to avoid surprising callers that check phase()
            // in this corner case.
            if (input_pair == CoolProp::HmassP_INPUTS || input_pair == CoolProp::HmolarP_INPUTS) {
                _phase = iphase_twophase;
            }
            break;
        case PointEvaluation::Kind::TwoPhaseDisallowed:
        case PointEvaluation::Kind::Invalid:
        case PointEvaluation::Kind::InternalError:
            break;
    }
}

// Property accessors — all route through evaluate_property_ over the
// active_eval_ kernel.  The throw cases (speed_sound / viscosity /
// conductivity in the two-phase region) propagate to the caller.
CoolPropDbl SVDSBTLBackend::calc_rhomass() {
    return evaluate_property_(active_eval_, CoolProp::iDmass);
}
CoolPropDbl SVDSBTLBackend::calc_hmass() {
    return evaluate_property_(active_eval_, CoolProp::iHmass);
}
CoolPropDbl SVDSBTLBackend::calc_smass() {
    return evaluate_property_(active_eval_, CoolProp::iSmass);
}
CoolPropDbl SVDSBTLBackend::calc_umass() {
    return evaluate_property_(active_eval_, CoolProp::iUmass);
}
CoolPropDbl SVDSBTLBackend::calc_T() {
    return evaluate_property_(active_eval_, CoolProp::iT);
}
CoolPropDbl SVDSBTLBackend::calc_speed_sound() {
    return evaluate_property_(active_eval_, CoolProp::ispeed_sound);
}
CoolPropDbl SVDSBTLBackend::calc_viscosity() {
    return evaluate_property_(active_eval_, CoolProp::iviscosity);
}
CoolPropDbl SVDSBTLBackend::calc_conductivity() {
    return evaluate_property_(active_eval_, CoolProp::iconductivity);
}
CoolPropDbl SVDSBTLBackend::calc_pressure() {
    return evaluate_property_(active_eval_, CoolProp::iP);
}

// Molar variants — defer to the mass version and scale by molar_mass.
CoolPropDbl SVDSBTLBackend::calc_rhomolar() {
    return calc_rhomass() / calc_molar_mass();
}
CoolPropDbl SVDSBTLBackend::calc_hmolar() {
    return calc_hmass() * calc_molar_mass();
}
CoolPropDbl SVDSBTLBackend::calc_smolar() {
    return calc_smass() * calc_molar_mass();
}
CoolPropDbl SVDSBTLBackend::calc_umolar() {
    return calc_umass() * calc_molar_mass();
}

// Constants — delegate to a lazily-allocated HEOS reference.
CoolPropDbl SVDSBTLBackend::calc_molar_mass() {
    if (molar_mass_cached_ < 0.0) {
        molar_mass_cached_ = source_reference_()->molar_mass();
    }
    return molar_mass_cached_;
}
CoolPropDbl SVDSBTLBackend::calc_gas_constant() {
    return source_reference_()->gas_constant();
}
CoolPropDbl SVDSBTLBackend::calc_T_critical() {
    return source_reference_()->T_critical();
}
CoolPropDbl SVDSBTLBackend::calc_p_critical() {
    return source_reference_()->p_critical();
}
CoolPropDbl SVDSBTLBackend::calc_rhomass_critical() {
    return source_reference_()->rhomass_critical();
}
CoolPropDbl SVDSBTLBackend::calc_rhomolar_critical() {
    return source_reference_()->rhomolar_critical();
}
CoolPropDbl SVDSBTLBackend::calc_Ttriple() {
    return source_reference_()->Ttriple();
}
CoolPropDbl SVDSBTLBackend::calc_p_triple() {
    return source_reference_()->p_triple();
}
CoolPropDbl SVDSBTLBackend::calc_Tmax() {
    return source_reference_()->Tmax();
}
CoolPropDbl SVDSBTLBackend::calc_Tmin() {
    return source_reference_()->Tmin();
}
CoolPropDbl SVDSBTLBackend::calc_pmax() {
    return source_reference_()->pmax();
}

// fast_evaluate: vectorized cache-bypassing batch path.  Mirrors
// update()'s routing per point but writes directly to out_buffer
// without touching any per-instance cache slot (_T, _p, _phase,
// active_*).  All allocations live outside the loop.
//
// Supported pairs:
//   HmassP_INPUTS  : val1 = h, val2 = p
//   HmolarP_INPUTS : val1 = h_molar, val2 = p  (converted to h_mass)
//   PT_INPUTS      : val1 = p,   val2 = T
//
// Per-point outcome:
//   - in-domain single-phase  → resolve + N_outputs evaluations + ok
//   - inside critical-patch   → patch_source_->update() + property reads
//   - HmassP / HmolarP in dome → Q-weighted lever-rule blend through
//                                 resolve_point_ + evaluate_property_
//   - PT exactly on sat curve  → NaN row + fast_evaluate_two_phase_disallowed
//   - out of every region      → NaN row + fast_evaluate_out_of_range
//
// In-patch points get the same patch_source_->update() routing as
// update() does — omitting the critical-patch handoff would silently
// degrade accuracy at the states the build was explicitly configured
// to cover.
void SVDSBTLBackend::fast_evaluate(CoolProp::input_pairs input_pair, const double* val1, const double* val2, std::size_t N_inputs,
                                   const CoolProp::parameters* outputs, std::size_t N_outputs, double* out_buffer, std::size_t out_buffer_size,
                                   int* status_flags, std::size_t status_flags_size, CoolProp::phases /*imposed_phase*/) {
    // Bounds + null checks first (matches IF97Backend::fast_evaluate).
    if (N_outputs != 0 && N_inputs > std::numeric_limits<std::size_t>::max() / N_outputs) {
        throw ValueError(format("fast_evaluate: N_inputs * N_outputs would overflow size_t (N_inputs=%zu, N_outputs=%zu)", N_inputs, N_outputs));
    }
    const std::size_t required_out = N_inputs * N_outputs;
    if (out_buffer_size < required_out) {
        throw ValueError(format("fast_evaluate: out_buffer_size=%zu < required %zu (N_inputs * N_outputs)", out_buffer_size, required_out));
    }
    if (status_flags_size < N_inputs) {
        throw ValueError(format("fast_evaluate: status_flags_size=%zu < required %zu (N_inputs)", status_flags_size, N_inputs));
    }
    if (N_inputs == 0) return;
    if (val1 == nullptr || val2 == nullptr || outputs == nullptr || out_buffer == nullptr || status_flags == nullptr) {
        throw ValueError("fast_evaluate: null pointer argument");
    }
    if (N_outputs == 0) {
        for (std::size_t k = 0; k < N_inputs; ++k)
            status_flags[k] = CoolProp::fast_evaluate_ok;
        return;
    }
    if (input_pair != CoolProp::HmassP_INPUTS && input_pair != CoolProp::HmolarP_INPUTS && input_pair != CoolProp::PT_INPUTS) {
        throw ValueError(format("fast_evaluate (SVDSBTL): input_pair %s not supported (use HmassP_INPUTS, HmolarP_INPUTS, or PT_INPUTS)",
                                get_input_pair_short_desc(input_pair).c_str()));
    }
    // Validate the surface exists for this input pair before the per-
    // point loop so a malformed request fails fast (resolve_point_
    // would otherwise throw on the first call).
    const auto surface_key = (input_pair == CoolProp::HmolarP_INPUTS) ? CoolProp::HmassP_INPUTS : input_pair;
    const auto sit = surfaces_.find(static_cast<int>(surface_key));
    if (sit == surfaces_.end() || !sit->second) {
        throw ValueError(std::string("fast_evaluate (SVDSBTL): no surface registered for ") + CoolProp::get_input_pair_short_desc(input_pair));
    }

    const double NaN = std::numeric_limits<double>::quiet_NaN();
    auto fill_nan_row = [&](std::size_t k) {
        for (std::size_t o = 0; o < N_outputs; ++o) {
            out_buffer[k * N_outputs + o] = NaN;
        }
    };

    // Build an output plan ONCE outside the per-point loop.
    //
    // SVDSurface::eval_with_region_multi can evaluate N properties at
    // the same (region, svd_x, svd_y) for the price of one locate() +
    // Hermite-basis-setup pair, amortizing those ~25 ns across the
    // batch.  To use it we need to know up-front which of the user's
    // outputs map to surface evals, what mass-basis surface key each
    // corresponds to (iDmolar → iDmass × M, iHmolar → iHmass × M, etc.),
    // and how to scale on output (multiply / divide by M).
    //
    // Per-output role:
    //   * Surface:        consumed from the batched surface_vals array
    //   * Shortcut:       served by evaluate_property_ (input-pair
    //                     shortcuts + iQ + kind-dependent paths)
    // We dedup the surface-prop list so requesting both iDmass and
    // iDmolar costs one eval, not two.
    //
    // Shortcut eligibility depends on which input pair is active —
    // e.g., iT is a shortcut only for PT_INPUTS (pt.T is the input);
    // for HmassP / HmolarP it has to come off the surface, in which
    // case routing it through the batch saves an extra locate().
    // Likewise iHmass / iHmolar are shortcuts only when h is supplied
    // by the input pair.
    const bool input_is_PT = (input_pair == CoolProp::PT_INPUTS);
    const bool input_has_h = (input_pair == CoolProp::HmassP_INPUTS || input_pair == CoolProp::HmolarP_INPUTS);

    // N_outputs is bounded by the caller's output schema (typically <
    // 10 for thermo-property batches), so a single heap allocation
    // per fast_evaluate call is negligible against an N_inputs-deep
    // per-point loop and avoids a hard-coded compile-time ceiling.
    std::vector<unsigned char> output_is_surface(N_outputs, 0);
    std::vector<std::size_t> output_surface_idx(N_outputs, 0);
    std::vector<double> output_M_scale(N_outputs, 1.0);
    std::vector<CoolProp::parameters> surface_props;
    surface_props.reserve(N_outputs);
    const double M = calc_molar_mass();
    for (std::size_t o = 0; o < N_outputs; ++o) {
        CoolProp::parameters mass_prop = outputs[o];
        double scale = 1.0;
        switch (outputs[o]) {
            case CoolProp::iP:
            case CoolProp::iQ:
                // Always-available shortcuts: iP comes off the input
                // pair; iQ resolves trivially (-1 for SinglePhase, Q
                // for DomeBlend) inside evaluate_property_.
                output_is_surface[o] = 0;
                continue;
            case CoolProp::iT:
                if (input_is_PT) {
                    output_is_surface[o] = 0;
                    continue;
                }
                // Falls through to surface eval (mass_prop = iT, scale = 1).
                break;
            case CoolProp::iHmass:
                if (input_has_h) {
                    // pt.h_mass is populated directly from the input.
                    output_is_surface[o] = 0;
                    continue;
                }
                // PT_INPUTS: needs a surface lookup.  mass_prop = iHmass.
                break;
            case CoolProp::iHmolar:
                if (input_has_h) {
                    // evaluate_property_ shortcut: *pt.h_mass * M.
                    output_is_surface[o] = 0;
                    continue;
                }
                // PT_INPUTS: surface returns iHmass, scale by M.
                mass_prop = CoolProp::iHmass;
                scale = M;
                break;
            case CoolProp::iDmolar:
                mass_prop = CoolProp::iDmass;
                scale = 1.0 / M;
                break;
            case CoolProp::iSmolar:
                mass_prop = CoolProp::iSmass;
                scale = M;
                break;
            case CoolProp::iUmolar:
                mass_prop = CoolProp::iUmass;
                scale = M;
                break;
            default:
                break;
        }
        // Surface-served output.  Dedup against already-listed props
        // (cheap linear scan; N_outputs is small).
        std::size_t idx = surface_props.size();
        for (std::size_t s = 0; s < surface_props.size(); ++s) {
            if (surface_props[s] == mass_prop) {
                idx = s;
                break;
            }
        }
        if (idx == surface_props.size()) {
            surface_props.push_back(mass_prop);
        }
        output_is_surface[o] = 1;
        output_surface_idx[o] = idx;
        output_M_scale[o] = scale;
    }
    const std::size_t n_surface_props = surface_props.size();

    // Per-point loop.  resolve_point_ does all the routing; for
    // SinglePhase kinds we batch the surface evals via the multi-eval
    // path; everything else (DomeBlend, Patched, shortcuts) falls
    // through evaluate_property_ per output.
    std::vector<double> surface_vals(n_surface_props, 0.0);
    for (std::size_t k = 0; k < N_inputs; ++k) {
        PointEvaluation pt;
        try {
            pt = resolve_point_(input_pair, val1[k], val2[k]);
        } catch (const std::exception&) {
            status_flags[k] = CoolProp::fast_evaluate_internal_error;
            fill_nan_row(k);
            continue;
        }
        if (pt.kind == PointEvaluation::Kind::OutOfRange || pt.kind == PointEvaluation::Kind::TwoPhaseDisallowed
            || pt.kind == PointEvaluation::Kind::Invalid || pt.kind == PointEvaluation::Kind::InternalError) {
            status_flags[k] = pt.status;
            fill_nan_row(k);
            continue;
        }

        // SinglePhase: one batched surface multi-eval covers every
        // surface-served output; shortcuts still go through
        // evaluate_property_.
        const bool can_batch = (pt.kind == PointEvaluation::Kind::SinglePhase) && n_surface_props > 0;
        if (can_batch) {
            try {
                pt.surface->eval_with_region_multi(pt.resolved.region_idx, pt.resolved.svd_x, pt.resolved.svd_y, surface_props.data(),
                                                   n_surface_props, surface_vals.data());
            } catch (const std::exception&) {
                status_flags[k] = CoolProp::fast_evaluate_internal_error;
                fill_nan_row(k);
                continue;
            }
        }

        bool row_failed = false;
        for (std::size_t o = 0; o < N_outputs; ++o) {
            double v = NaN;
            if (can_batch && output_is_surface[o]) {
                v = surface_vals[output_surface_idx[o]] * output_M_scale[o];
            } else {
                try {
                    v = evaluate_property_(pt, outputs[o]);
                } catch (const NotImplementedError&) {
                    // Two-phase undefined property (speed_sound, viscosity,
                    // conductivity, cp, cv, Prandtl).  In-band NaN signals
                    // "no equilibrium value here"; status stays ok so the
                    // caller can still read the meaningful T / D / S / U
                    // / Q outputs from the same row.
                    v = NaN;
                } catch (const std::exception&) {
                    v = NaN;
                    row_failed = true;
                }
            }
            out_buffer[k * N_outputs + o] = v;
        }
        if (row_failed) {
            status_flags[k] = CoolProp::fast_evaluate_internal_error;
            fill_nan_row(k);
        } else {
            status_flags[k] = pt.status;
        }
    }
}

}  // namespace CoolProp
