#ifndef COOLPROP_SVDSBTL_BACKEND_H
#define COOLPROP_SVDSBTL_BACKEND_H

#include <limits>
#include <memory>
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
    // AbstractState::fast_evaluate for the contract.  Mirrors update()'s
    // routing (critical-patch bbox → patch_source_, otherwise SVD
    // surface) but never touches the instance's _T / _p / _phase /
    // active_point_ slots, so back-to-back fast_evaluate batches do not
    // interfere with each other and a concurrent update() call.  No
    // heap allocations beyond what each in-patch point's
    // patch_source_->update() does internally.  Two-phase / dome-hit
    // points NaN-fill their row and set fast_evaluate_two_phase_disallowed.
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

   private:
    // Returns the prop's value at the active lookup point.  Throws
    // CoolProp::ValueError if no update() has been called or the
    // active surface doesn't expose `prop`.
    CoolPropDbl lookup_(CoolProp::parameters prop);

    // Load <fluid>.<pair>.svd.bin.z from the default cache; if absent,
    // build via the matching preset and save it.  Inserts into surfaces_.
    void ensure_surface_(CoolProp::input_pairs pair);

    // Resolve the SuperAncillary handle for this fluid (cached after
    // first call).  Returns nullptr if no SuperAncillary ships with
    // this fluid -- two-phase queries then throw a clear error.
    std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> superanc_();

    // Two-phase update for PQ_INPUTS / QT_INPUTS / dome-hit HmassP /
    // dome-hit PT.  Caches the sat-line state (T_sat, p_sat, Q, plus
    // sat-line endpoints rho_L/V, h_L/V, s_L/V, u_L/V on a molar basis)
    // so the calc_* virtuals can produce Q-weighted properties without
    // re-evaluating the superancillary.
    void update_two_phase_(double T_sat, double p_sat, double Q);

    // Per-property Q-weighted blend over the cached sat-line state.
    // Density uses 1/(Q/rho_V + (1-Q)/rho_L) (specific-volume blend);
    // h/s/u use the linear blend.
    [[nodiscard]] CoolPropDbl two_phase_property_(CoolProp::parameters prop) const;

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

    // State recorded by the most recent update() call.
    CoolProp::input_pairs active_pair_ = CoolProp::INPUT_PAIR_INVALID;
    const CoolProp::sbtl::SVDSurface* active_surface_ = nullptr;
    CoolProp::sbtl::SVDSurface::ResolvedPoint active_point_{};
    bool active_resolved_ = false;

    // Mass-basis inputs stashed by update().  Returned verbatim when
    // the user asks for the property they just supplied (no SVD
    // round-trip).  _T and _p live in AbstractState already.
    double active_hmass_ = std::numeric_limits<double>::quiet_NaN();

    // Lazy-resolved SuperAncillary handle.  Resolved on first call to
    // superanc_(); cached for subsequent two-phase queries.  Stays
    // nullptr if the fluid ships without one.
    std::shared_ptr<superancillary::SuperAncillary<std::vector<double>>> superanc_cached_;
    bool superanc_resolved_ = false;  // distinguish "not tried yet" from "tried and got nullptr"

    // Active two-phase state, populated by update_two_phase_().  All
    // four endpoint properties are stored on a MOLAR basis so the
    // Q-weighted blend math doesn't need to redo the unit conversion
    // on every calc_*.
    struct TwoPhaseState
    {
        bool active = false;
        double rhoL_mol = 0.0;
        double rhoV_mol = 0.0;
        double hL_mol = 0.0;
        double hV_mol = 0.0;
        double sL_mol = 0.0;
        double sV_mol = 0.0;
        double uL_mol = 0.0;
        double uV_mol = 0.0;
    };
    TwoPhaseState two_phase_{};

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
    bool active_in_patch_ = false;

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
