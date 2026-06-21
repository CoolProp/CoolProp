#ifndef COOLPROP_SVDSBTL_BACKEND_H
#define COOLPROP_SVDSBTL_BACKEND_H

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <set>  // transitively needed by superancillary.h
#include <string>
#include <unordered_map>
#include <vector>

#include "CoolProp/superancillary/superancillary.h"

#include "CoolProp/AbstractState.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SaturationSurrogate.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"

namespace CoolProp {

// SVDSBTL_BACKEND — fast pure-fluid property lookup backed by the
// rank-r SVD machinery from Phase 1 + 2b.  Architecturally, mirrors
// IF97Backend rather than BicubicBackend / TabularBackend: inherits
// directly from AbstractState, manages its own input_pair dispatch,
// and consumes already-built SVDSurface artifacts.
//
// One SVDSurface per supported input pair.  The MVP ships PH and PT;
// future phases register DU, HS, etc. by adding entries to the
// internal surfaces_ map — no changes to update(), calc_*, or the
// serialization layer.
//
// Lifecycle:
//
//   factory("SVDSBTL", "Water")
//     ↓
//   SVDSBTLBackend(fluid="Water")
//     ↓  (try load + build-if-missing for PH, PT)
//   surfaces_[HmassP_INPUTS] = SVDSurface{...}
//   surfaces_[PT_INPUTS]     = SVDSurface{...}
//     ↓
//   update(input_pair, v1, v2)
//     ↓  (record inputs, resolve atlas dispatch once)
//   calc_rhomass() / calc_T() / ... → eval_with_region(prop, ...)
//
// Mixtures are rejected at construction time — pure-fluid only.
// Two-phase queries return NaN (no SVDSurface covers the dome yet);
// future phases may add an in-dome lookup or fall back to HEOS.
//
// Critical-point fallback: not implemented in MVP.  The subcritical
// SVD surfaces stop short of p_crit; queries inside the critical
// micro-region currently return NaN.  A direct-HEOS fallback is a
// natural follow-up.
class SVDSBTLBackend : public AbstractState
{
   public:
    // Construct an SVDSBTL backend backed by the given truth source.
    //
    // The source argument is required (no default) — it determines:
    //   * which backend supplies the HEOS-style truth values when the
    //     SVD tables are built (or loaded from cache)
    //   * the cache key slot in <fluid>.<source>.<input_pair>.svd.bin.z
    //     so HEOS-built and REFPROP-built tables for the same fluid
    //     never collide
    //   * how two-phase queries are resolved (HEOS uses SuperAncillary;
    //     REFPROP forwards PQ/QT directly; IF97 uses psat97 + Q-blend)
    //   * the critical-patch fallback target (always the source backend,
    //     so the patch is "truth" by definition for each source)
    //
    // Supported source values: "HEOS", "REFPROP", "IF97".  Anything
    // else throws.  IF97 is restricted to Water.
    //
    // `options_json` is the raw JSON payload from the factory string's
    // `?<options>` suffix (or "" when the caller didn't supply one).
    // Validated against `kSVDSBTLOptionsSchemaJson`; unknown keys
    // throw immediately.  Defaults expanded at construction so
    // `build_options_json()` returns the canonical form.  See
    // Web/coolprop/BackendOptions.rst for the schema and grammar.
    SVDSBTLBackend(const std::string& fluid_name, const std::string& source_backend, const std::string& options_json = "");

    SVDSBTLBackend(const SVDSBTLBackend&) = delete;
    SVDSBTLBackend& operator=(const SVDSBTLBackend&) = delete;
    SVDSBTLBackend(SVDSBTLBackend&&) = delete;
    SVDSBTLBackend& operator=(SVDSBTLBackend&&) = delete;
    ~SVDSBTLBackend() override = default;

    // AbstractState identification.
    std::string backend_name() override {
        return get_backend_string(SVDSBTL_BACKEND);
    }

    // Block use through the high-level PropsSI interface by default.
    //
    // PropsSI calls AbstractState::factory() on every invocation, which
    // for SVDSBTL means loading + msgpack-deserializing the .svd.bin.z
    // cache file (~80 ms / call) before the actual SVD eval (~5 us).
    // Per-call cost is therefore ~4 orders of magnitude slower than
    // expected — a measured 77 ms / PropsSI vs 3.7 us for IF97 — which
    // turned the IF97 conformance docs script (~100k PropsSI calls)
    // into a 2+ hour job and pushed docs builds past the 3-hour
    // ceiling.  See also the BICUBIC / TTSE tabular backends which
    // are gated the same way.
    //
    // The ALLOW_SVDSBTL_IN_PROPSSI configuration key lets advanced
    // callers opt back in (e.g. for one-off interactive queries where
    // the cache load cost is fine).  For any throughput-sensitive
    // workload, use AbstractState directly + update() in a loop, or
    // fast_evaluate for a vectorized batch.
    bool available_in_high_level() override;

    // Canonical JSON of the options this instance was built with
    // (validated + schema defaults expanded).  Returns the literal "{}"
    // when no options were supplied; round-trips through factory().
    [[nodiscard]] std::string build_options_json() const override {
        return options_canonical_.empty() ? std::string{"{}"} : options_canonical_;
    }

    // Pure fluid only.  Reject mixtures via set_mole_fractions().
    bool using_mole_fractions() override {
        return true;
    }
    bool using_mass_fractions() override {
        return false;
    }
    bool using_volu_fractions() override {
        return false;
    }
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) override;
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) override {
        // SVDSBTL is pure-fluid only; for pure fluids Qmass == Qmolar
        // exactly and the only valid composition is {1.0}.  Route to
        // set_mole_fractions so we share the size/value validation
        // there instead of silently accepting any vector.
        set_mole_fractions(mass_fractions);
    }
    void set_volu_fractions(const std::vector<CoolPropDbl>& /*volu_fractions*/) override {
        throw NotImplementedError("SVDSBTL backend does not support volume fractions");
    }
    const std::vector<CoolPropDbl>& get_mole_fractions() override {
        return mole_fractions_;
    }

    // Main update dispatcher.  Hands off to resolve_point_(), stores
    // the resulting PointEvaluation as active_eval_, and copies the
    // legacy AbstractState slots (_T / _p / _Q / _phase) for back-
    // compat readers.  Subsequent calc_*() calls dispatch through
    // evaluate_property_(active_eval_, ...) — which in the SinglePhase
    // case is a single eval_with_region per output, in the DomeBlend
    // case is a lever-rule blend with lazy sat-line endpoint fills,
    // and in the Patched case forwards to patch_source_.
    void update(CoolProp::input_pairs input_pair, double value1, double value2) override;

    // Vectorized cache-bypassing batch evaluation — see
    // AbstractState::fast_evaluate for the contract.  Both update() and
    // fast_evaluate() route through the same private resolve_point_()
    // kernel; this entry point just loops over the input array, asks
    // the kernel for a PointEvaluation per probe, then dispatches each
    // requested output through evaluate_property_().  Writes directly
    // into the caller's out_buffer instead of mutating active_eval_.
    //
    // **Thread safety / interleaving:** The backend is not thread-safe
    // — two concurrent calls (any combination of update / fast_evaluate
    // / property reads) on the same instance can corrupt each other's
    // state.  Single-threaded back-to-back fast_evaluate calls and
    // single-threaded fast_evaluate interleaved with update() are
    // **only** safe when no probe in the batch falls inside the
    // critical-patch bbox.  Patched probes re-update the lazily
    // allocated patch_source_ child AbstractState in place, which
    // means a sequence like `update(state A); fast_evaluate(...batch
    // containing a Patched probe B...); state->T()` returns T from B
    // rather than A.  In practice the critical-patch bbox is small
    // (~5 % of T_c × ~30 % of p_c by default) so the corruption is
    // local to true near-critical workloads.
    //
    // HmassP / HmolarP probes inside the saturation dome are
    // Q-blended in place; PT probes exactly on the sat curve NaN-fill
    // with fast_evaluate_two_phase_disallowed; points outside every
    // region NaN-fill with fast_evaluate_out_of_range.
    void fast_evaluate(CoolProp::input_pairs input_pair, const double* val1, const double* val2, std::size_t N_inputs,
                       const CoolProp::parameters* outputs, std::size_t N_outputs, double* out_buffer, std::size_t out_buffer_size, int* status_flags,
                       std::size_t status_flags_size, CoolProp::phases imposed_phase = CoolProp::iphase_not_imposed) override;

    // Property accessors.  All route through
    // evaluate_property_(active_eval_, ...), which honours the
    // inputs the user supplied (no SVD round-trip for values the
    // user just gave us — h on HmassP, T on PT, etc.).
    CoolPropDbl calc_rhomass() override;
    CoolPropDbl calc_hmass() override;
    CoolPropDbl calc_smass() override;
    CoolPropDbl calc_umass() override;
    CoolPropDbl calc_T() override;
    CoolPropDbl calc_pressure() override;
    CoolPropDbl calc_rhomolar() override;
    CoolPropDbl calc_hmolar() override;
    CoolPropDbl calc_smolar() override;
    CoolPropDbl calc_umolar() override;

    CoolPropDbl calc_speed_sound() override;
    CoolPropDbl calc_viscosity() override;
    CoolPropDbl calc_conductivity() override;

    CoolPropDbl calc_molar_mass() override;
    CoolPropDbl calc_gas_constant() override;

    CoolPropDbl calc_T_critical() override;
    CoolPropDbl calc_p_critical() override;
    CoolPropDbl calc_rhomass_critical() override;
    CoolPropDbl calc_rhomolar_critical() override;

    CoolPropDbl calc_Ttriple() override;
    CoolPropDbl calc_p_triple() override;
    CoolPropDbl calc_Tmax() override;
    CoolPropDbl calc_Tmin() override;
    CoolPropDbl calc_pmax() override;

    phases calc_phase() override {
        return _phase;
    }

    std::vector<std::string> calc_fluid_names() override {
        return {fluid_name_};
    }
    std::string fluid_param_string(const std::string& ParamName) override {
        // HEOS exposes a range of per-fluid metadata strings via
        // fluid_param_string (name, aliases, CAS, formula, ...).  For
        // SVDSBTL we only need to make sure "name" round-trips
        // correctly; everything else routes to the HEOS reference
        // so callers don't lose access to those strings by switching
        // backends.  Keep the "name" / "aliases" fast paths inline
        // so they don't trip an unnecessary HEOS construction.
        if (ParamName == "name") {
            return fluid_name_;
        }
        return source_reference_()->fluid_param_string(ParamName);
    }

    // Test seam: exposes which surfaces were loaded vs built.  Returns
    // the set of input pairs registered in surfaces_.
    [[nodiscard]] std::vector<CoolProp::input_pairs> registered_input_pairs() const;

    // Test seam: true iff at least one sat lookup so far (any of
    // PQ/QT routing, PT-on-sat probing, HmassP dome detection, or
    // dome-endpoint fills) has gone through the SaturationSurrogate
    // cubic-spline cache rather than the source backend's QT/PQ_INPUTS
    // path or the per-fluid SuperAncillary.  Used by the CoolProp-077
    // wiring test to confirm the surrogate is actually consulted (HEOS
    // source -> false because SA is preferred; REFPROP source -> true).
    // Always false before the first sat-touching call; lazy build only
    // fires on demand.
    [[nodiscard]] bool sat_surrogate_consulted() const noexcept;

   public:
    // Self-contained per-point evaluation state — the shared currency
    // between update() (caches the result on the instance) and
    // fast_evaluate() (loops, holds it locally per probe).  Public so
    // tests + future batch callers can inspect what resolve_point_()
    // decided.
    struct PointEvaluation
    {
        enum class Kind : std::uint8_t
        {
            Invalid,             // never assigned
            SinglePhase,         // atlas resolved to a surface region
            DomeBlend,           // HmassP / HmolarP landed inside the dome OR PQ/QT input
            Patched,             // inside the critical-patch bbox; consult patch_source_
            OutOfRange,          // outside every region's bbox
            TwoPhaseDisallowed,  // PT exactly on sat curve (Q ambiguous)
            InternalError,       // per-point evaluation failure
        };

        Kind kind = Kind::Invalid;
        int status = CoolProp::fast_evaluate_internal_error;  // fast_evaluate_status code

        // The input pair this point was resolved for + the raw inputs.
        // Lets evaluate_property_() return the supplied input verbatim
        // (no surface round-trip) when the caller asks for it.  v1 /
        // v2 are always set after a valid resolve_point_ call.
        CoolProp::input_pairs input_pair = CoolProp::INPUT_PAIR_INVALID;
        double v1 = std::numeric_limits<double>::quiet_NaN();
        double v2 = std::numeric_limits<double>::quiet_NaN();

        // Universal physical coords resolved for this point.  Filled
        // for SinglePhase / DomeBlend / Patched; the OOB and
        // TwoPhaseDisallowed kinds leave them as nullopt because
        // there's no defined physical state at the input.  For
        // DomeBlend, T == T_sat and Q is in [0, 1].
        std::optional<double> T;
        std::optional<double> p;
        std::optional<double> h_mass;
        // Mass entropy input shortcut for PSmass / PSmolar inputs — the
        // entropy twin of h_mass.  s is the query input + secondary axis,
        // so evaluate_property_ returns it directly without a surface eval.
        std::optional<double> s_mass;
        // Mass density input shortcut for DmassT / DmolarT inputs.
        // Lets evaluate_property_ return iDmass directly without a
        // surface eval (D is the input, not a stored property).
        std::optional<double> rho_mass;
        double Q = -1.0;  // -1 sentinel for single-phase, [0, 1] in dome

        // SinglePhase / Patched routing.
        const CoolProp::sbtl::SVDSurface* surface = nullptr;
        CoolProp::sbtl::SVDSurface::ResolvedPoint resolved{};

        // DomeBlend endpoints (molar basis to match SuperAncillary's
        // native units; h_mass / T above is what calc_* / fast_evaluate
        // return for the row).  Sat-line h endpoints are filled at
        // resolve time (we need them to determine Q anyway); rho and
        // s endpoints fill lazily on first request via
        // ensure_dome_*_endpoints_().
        std::optional<double> T_sat;
        std::optional<double> hL_mol;
        std::optional<double> hV_mol;
        std::optional<double> rhoL_mol;
        std::optional<double> rhoV_mol;
        std::optional<double> sL_mol;
        std::optional<double> sV_mol;
    };

   private:
    // Resolve an (input_pair, v1, v2) triple into a PointEvaluation
    // without mutating any backend state.  Both update() and
    // fast_evaluate() call this — there's exactly one place where the
    // input-pair dispatch, critical-patch routing, atlas resolve, and
    // dome detection live.  Two-phase PQ_INPUTS / QT_INPUTS get the
    // same DomeBlend treatment as HmassP-inside-dome — same kernel,
    // same fields.
    PointEvaluation resolve_point_(CoolProp::input_pairs ip, double v1, double v2);

    // Evaluate one property against a resolved PointEvaluation.  Pure;
    // touches only the const data inside `pt` (the surface pointer or
    // the dome endpoints) plus the lazily-populated patch_source_ for
    // Patched rows.  Returns NaN for properties undefined in the dome
    // (speed_sound, cp, cv, viscosity, conductivity, Prandtl) when
    // kind == DomeBlend.  Throws when a request can't be served at all
    // (e.g. unsupported property on a surface that doesn't expose it).
    [[nodiscard]] double evaluate_property_(PointEvaluation& pt, CoolProp::parameters prop);

    // Lazy population of dome endpoint properties — only the ones the
    // caller actually requests are pulled from SuperAncillary / source.
    void ensure_dome_rho_endpoints_(PointEvaluation& pt);
    void ensure_dome_s_endpoints_(PointEvaluation& pt);
    void ensure_dome_h_endpoints_(PointEvaluation& pt);

    // Load <fluid>.<pair>.svd.bin.z from the default cache; if absent,
    // build via the matching preset and save it.  Inserts into surfaces_.
    void ensure_surface_(CoolProp::input_pairs pair);

    // Resolve the SuperAncillary handle for this fluid (cached after
    // first call).  Returns nullptr if no SuperAncillary ships with
    // this fluid -- two-phase queries then throw a clear error.
    std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> superanc_();

    // Resolve the SaturationSurrogate for this (fluid, source-backend).
    // Lazy: only built when superanc_() returns nullptr AND a hot-path
    // caller asks for it (sat_T_from_p_ / sat_eval_).  Cached for the
    // lifetime of the backend.  Returns nullptr if the source backend
    // can't provide enough QT_INPUTS probes to build a valid surrogate
    // — callers then fall through to the source-PQ/QT last-resort path.
    // Mirrors the superanc_() resolution pattern so adding a third
    // sat-provider source (e.g. cached blob) later is a one-line
    // addition.
    const CoolProp::sbtl::SaturationSurrogate* sat_surrogate_handle_();

    // Composite saturation lookups.  Each one tries the SuperAncillary
    // first; if unavailable, falls back to the SaturationSurrogate.
    // If BOTH are unavailable, returns NaN — callers then drop into
    // the per-site source-PQ/QT fallback path that already exists for
    // surrogate-build failures.
    //
    // Setting `record_use` to true tags the test-seam flag
    // `sat_surrogate_used_` whenever the surrogate (not the SA) is
    // consulted — drives the using_sat_surrogate() public test seam.
    // The fallback (source.update PQ/QT) does not toggle the flag.
    [[nodiscard]] double sat_T_from_p_(double p);
    [[nodiscard]] double sat_eval_(double T, char what, int side);

    // Near-dome reclassification guard for a single-phase atlas hit
    // (issue #3190).  The atlas LIQUID/VAPOR regions' dome-side boundary
    // is an interpolated saturation curve that can disagree with the
    // authoritative sat endpoints by its fit error (~1 J/kg), so a point
    // sitting on the true bubble/dew line resolves single-phase instead
    // of two-phase.  When the resolved point is within a thin eta band of
    // a region edge, this re-tests it against the same saturation
    // provider the forward PQ flash uses.  `what` is 'H' (HmassP/HmolarP)
    // or 'S' (PSmass/PSmolar); `value_mass` is the mass-basis h or s.
    // Uses ONLY the fast sat provider (SA -> surrogate): if no provider
    // is available it returns false (leaving `pt` untouched) so the
    // single-phase hot path never pays for an expensive source-PQ flash.
    // Returns true and rewrites `pt` to a DomeBlend iff the point is
    // genuinely inside the dome (yL <= y <= yV, boundary inclusive).
    [[nodiscard]] bool try_fast_dome_reclassify_(PointEvaluation& pt, char what, double p_val, double value_mass);

    std::string fluid_name_;
    std::string source_backend_;               // "HEOS" / "REFPROP" / "IF97"
    std::vector<CoolPropDbl> mole_fractions_;  // always {1.0}

    // Canonical-JSON form of the options blob (validated + defaults
    // expanded) that this instance was built with.  Used by both the
    // cache-filename hash and `build_options_json()`.  Empty for
    // backwards compatibility when no options were supplied
    // (treated as the schema's all-defaults configuration).
    std::string options_canonical_;

    // Cache of molar mass, populated lazily on the first
    // calc_molar_mass() call.  Avoids per-call shared_ptr deref +
    // virtual dispatch through source_reference_() on every
    // calc_hmolar / calc_rhomolar / calc_smolar.  Fluid is fixed at
    // construction so this never changes after first read.
    mutable CoolPropDbl molar_mass_cached_ = -1.0;

    // Reference HEOS state for constants (critical point, triple point,
    // bounds, molar mass, R).  Lazily allocated so that backend
    // construction is cheap even when only cached blobs are loaded.
    std::shared_ptr<CoolProp::AbstractState> source_reference_();
    std::shared_ptr<CoolProp::AbstractState> source_;

    // One surface per supported input pair.  Heap-allocated so the
    // ResolvedPoint pointers below stay valid across map growth.
    std::unordered_map<int /*input_pairs as int*/, std::unique_ptr<CoolProp::sbtl::SVDSurface>> surfaces_;

    // The most recent update()'s resolved point.  All calc_*() methods
    // dispatch through evaluate_property_(active_eval_, ...).
    PointEvaluation active_eval_;

    // Lazy-resolved SuperAncillary handle.  Resolved on first call to
    // superanc_(); cached for subsequent two-phase queries.  Stays
    // nullptr if the fluid ships without one.
    std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> superanc_cached_;
    bool superanc_resolved_ = false;  // distinguish "not tried yet" from "tried and got nullptr"

    // Lazy-resolved SaturationSurrogate handle.  See sat_surrogate_handle_()
    // for the build trigger.  Owned by the backend.  Built only when SA is
    // absent (REFPROP source today); HEOS source skips this entirely
    // because superanc_ provides faster + more accurate sat lookups.
    std::unique_ptr<CoolProp::sbtl::SaturationSurrogate> sat_surrogate_;
    bool sat_surrogate_resolved_ = false;

    // Test seam.  Flipped on the first hot-path call that actually
    // routes through sat_surrogate_ (as opposed to superanc_).  Set
    // only by sat_T_from_p_ / sat_eval_; never cleared during the
    // backend's lifetime.
    bool sat_surrogate_used_ = false;

    // Critical-patch HEOS-fallback.  When an active query point lands
    // inside the configured bbox (in whatever the input pair's native
    // coords are), update() forwards the call to `patch_source_`
    // instead of the SVD surface, and every calc_* below routes there
    // for the duration of the current state.  Bbox is precomputed at
    // construction per input pair so update() is a constant-time
    // axis-aligned containment check.  Disabled when
    // critical_patch.mode is "off" in the options blob.
    struct CriticalPatchBbox
    {
        double a_lo, a_hi;  // primary axis: p for both PT and HmassP
        double b_lo, b_hi;  // secondary axis: T for PT, h for HmassP
    };
    struct CriticalPatch
    {
        bool enabled = false;
        std::string source;  // "HEOS" / "REFPROP" / "IF97"
        std::unordered_map<int /*input_pairs as int*/, CriticalPatchBbox> bbox_per_pair;
        [[nodiscard]] bool contains(::CoolProp::input_pairs pair, double a, double b) const noexcept;
    };
    CriticalPatch critical_patch_{};
    std::shared_ptr<CoolProp::AbstractState> patch_source_;  // lazy
    // active_in_patch_ removed — the patch / non-patch routing now
    // lives on active_eval_.kind (Patched vs SinglePhase/etc.).

    // Configure critical_patch_ from the validated options blob.
    // Called once at construction after surfaces_ are loaded.  When
    // mode == "auto":
    //   - tries to load 4 calibrated multipliers from a sidecar cache
    //     file; if found, uses them
    //   - otherwise runs auto_calibrate_critical_bbox_() which
    //     binary-search-shrinks each of the 4 (T_lo, T_hi, p_lo, p_hi)
    //     axes until the strip just outside the candidate patch passes
    //     IAPWS conformance budgets in (rho, h, s, w) against the
    //     source backend.  The resulting multipliers are persisted to
    //     the sidecar cache for the next construction.
    // When "fixed" honours the `bbox` array verbatim; when "off"
    // leaves enabled=false and queries always go to the SVD.
    void build_critical_patch_(const std::string& options_canonical);

    // Auto-calibration shrink loop (CoolProp-dxd).  Returns the four
    // multipliers {T_lo_mult, T_hi_mult, p_lo_mult, p_hi_mult} such
    // that the SVD's reconstruction of (rho, h, s, w) at every probe
    // *just outside* the resulting (T, p) bbox satisfies the
    // calibration budgets (1% relative in rho/h/s, 5% relative in w;
    // see kBudget_* in the .cpp).  Starts from the Water-sized
    // default {0.95, 1.05, 0.75, 1.15} and binary-search-shrinks
    // each axis independently toward 1.0 in 6 steps; never widens
    // beyond the default.  Total cost dominated by ~120 PT probes
    // against the source backend (~ms each for HEOS, ~10 ms each
    // for REFPROP) → a few seconds.  Persists to the sidecar cache
    // so subsequent constructions skip the work.
    [[nodiscard]] std::array<double, 4> auto_calibrate_critical_bbox_();

    // Sidecar persistence for the 4 calibrated multipliers.  Cache
    // path: $HOME/.CoolProp/SVDTables/<Fluid>.<Source>.critpatch.<OptHash>.bin
    // (matches the SVDSurface cache directory + opthash so user-level
    // table cleanup catches both files).  Format: 8-byte magic
    // "CPCRITPB" + uint32 version + 4 doubles = 44 bytes.  Tiny
    // enough to be plain binary; no msgpack / zlib overhead.
    [[nodiscard]] static std::optional<std::array<double, 4>> load_critpatch_cache_(const std::string& fluid_name, const std::string& source_backend,
                                                                                    const std::string& options_canonical);
    static void save_critpatch_cache_(const std::string& fluid_name, const std::string& source_backend, const std::string& options_canonical,
                                      const std::array<double, 4>& mults);

    // Get the patch-side source backend, allocating on first use.
    // The patch source can differ from the SVD truth source (set via
    // `critical_patch.source` in the options blob) — e.g. SVDSBTL&IF97
    // can still route patch queries through HEOS for full IAPWS-95
    // accuracy at the critical singularity.  Defaults to the SVD
    // truth source when null.
    std::shared_ptr<CoolProp::AbstractState> patch_source_ref_();
};

}  // namespace CoolProp

#endif  // COOLPROP_SVDSBTL_BACKEND_H
