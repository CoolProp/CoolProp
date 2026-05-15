#ifndef COOLPROP_SVDSBTL_BACKEND_H
#define COOLPROP_SVDSBTL_BACKEND_H

#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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
    explicit SVDSBTLBackend(const std::string& fluid_name);

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
    void set_mass_fractions(const std::vector<CoolPropDbl>& /*mass_fractions*/) override {}
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
    std::string fluid_param_string(const std::string& /*ParamName*/) override {
        return fluid_name_;
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

    std::string fluid_name_;
    std::vector<CoolPropDbl> mole_fractions_;  // always {1.0}

    // Reference HEOS state for constants (critical point, triple point,
    // bounds, molar mass, R).  Lazily allocated so that backend
    // construction is cheap even when only cached blobs are loaded.
    std::shared_ptr<CoolProp::AbstractState> heos_reference_();
    std::shared_ptr<CoolProp::AbstractState> heos_;

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
};

}  // namespace CoolProp

#endif  // COOLPROP_SVDSBTL_BACKEND_H
