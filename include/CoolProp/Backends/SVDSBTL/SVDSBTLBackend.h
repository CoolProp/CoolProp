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
    SVDSBTLBackend(const std::string& fluid_name, const std::string& source_backend);

    SVDSBTLBackend(const SVDSBTLBackend&) = delete;
    SVDSBTLBackend& operator=(const SVDSBTLBackend&) = delete;
    SVDSBTLBackend(SVDSBTLBackend&&) = delete;
    SVDSBTLBackend& operator=(SVDSBTLBackend&&) = delete;
    ~SVDSBTLBackend() override = default;

    // AbstractState identification.
    std::string backend_name() override {
        return get_backend_string(SVDSBTL_BACKEND);
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
    std::string source_backend_;  // "HEOS" / "REFPROP" / "IF97"
    std::vector<CoolPropDbl> mole_fractions_;  // always {1.0}

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
};

}  // namespace CoolProp

#endif  // COOLPROP_SVDSBTL_BACKEND_H
