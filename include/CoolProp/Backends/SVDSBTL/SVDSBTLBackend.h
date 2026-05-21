#ifndef COOLPROP_SVDSBTL_BACKEND_H
#define COOLPROP_SVDSBTL_BACKEND_H

#include <limits>
#include <memory>
#include <optional>
#include <set>  // transitively needed by superancillary.h
#include <string>
#include <unordered_map>
#include <vector>

#include "superancillary/superancillary.h"

#include "AbstractState.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "DataStructures.h"
#include "Exceptions.h"

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

    // Main update dispatcher.  Routes by input_pair to the matching
    // surface and resolves the atlas-dispatch ResolvedPoint once so
    // subsequent calc_* calls are pure SVD evaluations.
    void update(CoolProp::input_pairs input_pair, double value1, double value2) override;

    // Vectorized cache-bypassing batch evaluation — see
    // AbstractState::fast_evaluate for the contract.  Both update() and
    // fast_evaluate() route through the same private resolve_point_()
    // kernel; this entry point just loops over the input array, asks
    // the kernel for a PointEvaluation per probe, then dispatches each
    // requested output through evaluate_property_().  Writes directly
    // into the caller's out_buffer instead of mutating active_eval_, so
    // back-to-back fast_evaluate calls can be interleaved with update()
    // on the same instance without corrupting the most-recent update()
    // — BUT the patch routing is stateful (reuses the lazily-allocated
    // patch_source_ child AbstractState across all in-patch points),
    // which makes fast_evaluate safe for single-threaded interleaving
    // with update() but NOT thread-safe against a concurrent
    // fast_evaluate/update on the same backend.  HmassP / HmolarP
    // probes inside the saturation dome are Q-blended in place; PT
    // probes exactly on the sat curve NaN-fill with
    // fast_evaluate_two_phase_disallowed; points outside every region
    // NaN-fill with fast_evaluate_out_of_range.
    void fast_evaluate(CoolProp::input_pairs input_pair, const double* val1, const double* val2, std::size_t N_inputs,
                       const CoolProp::parameters* outputs, std::size_t N_outputs, double* out_buffer, std::size_t out_buffer_size, int* status_flags,
                       std::size_t status_flags_size, CoolProp::phases imposed_phase = CoolProp::iphase_not_imposed) override;

    // Property accessors.  Each routes through lookup_(prop), which
    // honours the inputs the user supplied (no SVD round-trip for
    // values the user just gave us).
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

    // Load <fluid>.<pair>.svd.bin.z from the default cache; if absent,
    // build via the matching preset and save it.  Inserts into surfaces_.
    void ensure_surface_(CoolProp::input_pairs pair);

    // Resolve the SuperAncillary handle for this fluid (cached after
    // first call).  Returns nullptr if no SuperAncillary ships with
    // this fluid -- two-phase queries then throw a clear error.
    std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> superanc_();

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
    // mode == "auto" uses the default multipliers (PR D plugs in
    // real auto-calibration); when "fixed" honours `bbox`; when
    // "off" leaves enabled=false and queries always go to the SVD.
    void build_critical_patch_(const std::string& options_canonical);

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
