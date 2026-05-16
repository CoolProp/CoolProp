// SVDSBTLBackend implementation.  Architecturally inherits from
// AbstractState (mirroring IF97Backend), holds an SVDSurface per
// supported input pair, and dispatches calc_* through the active
// surface resolved at update() time.

#include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"

#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#include "DataStructures.h"
#include "Exceptions.h"

namespace CoolProp {

namespace {

// The MVP set of input pairs.  Adding DU / HS / etc. in a follow-up
// is a one-line append here plus the matching preset + load helpers
// in `build_surface_for_pair_`.
constexpr std::array<CoolProp::input_pairs, 2> kSupportedPairs = {CoolProp::HmassP_INPUTS, CoolProp::PT_INPUTS};

// Sample-and-build for a given input pair using the matching preset.
CoolProp::sbtl::SVDSurface build_surface_for_pair_(::CoolProp::AbstractState& source, CoolProp::input_pairs pair) {
    namespace cp_sbtl = CoolProp::sbtl;
    cp_sbtl::SurfaceSpec spec;
    switch (pair) {
        case CoolProp::HmassP_INPUTS:
            spec = cp_sbtl::presets::ph_subcritical(source);
            break;
        case CoolProp::PT_INPUTS:
            spec = cp_sbtl::presets::pt_subcritical(source);
            break;
        default:
            throw ValueError("SVDSBTL backend: no preset registered for the requested input pair");
    }
    return cp_sbtl::build_surface(source, std::move(spec));
}

}  // namespace

SVDSBTLBackend::SVDSBTLBackend(const std::string& fluid_name, const std::string& source_backend)
  : fluid_name_(fluid_name), source_backend_(source_backend), mole_fractions_({1.0}) {
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
    const std::string path = cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid_name_, source_backend_, pair);
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
        surface = std::make_unique<cp_sbtl::SVDSurface>(build_surface_for_pair_(*source, pair));
        try {
            cp_sbtl::SVDSurfaceSerializer::save_to_file(*surface, path);
        } catch (const std::exception&) {
            // Cache write failure is non-fatal — the in-memory surface
            // still works for this run.
        }
    }

    surfaces_[static_cast<int>(pair)] = std::move(surface);
}

void SVDSBTLBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    clear();
    _phase = iphase_not_imposed;
    two_phase_.active = false;

    // Two-phase inputs go to the superancillary path immediately --
    // no SVDSurface needed.  PQ_INPUTS and QT_INPUTS are degenerate
    // in the surface (the dome interior is not single-phase, so the
    // surface returns NaN) so we route around it.
    if (input_pair == CoolProp::PQ_INPUTS || input_pair == CoolProp::QT_INPUTS) {
        active_pair_ = input_pair;
        active_surface_ = nullptr;
        active_resolved_ = false;
        auto sa = superanc_();
        if (!sa) {
            throw NotImplementedError(std::string("SVDSBTL backend: two-phase inputs require a SuperAncillary for this fluid (") + fluid_name_
                                      + " has none).  Use a different backend (e.g. HEOS) for two-phase queries on this fluid.");
        }
        double p_sat = 0.0;
        double T_sat = 0.0;
        double Q = 0.0;
        if (input_pair == CoolProp::PQ_INPUTS) {
            // value1 = p, value2 = Q.
            p_sat = value1;
            Q = value2;
            T_sat = sa->get_T_from_p(p_sat);
        } else {
            // QT_INPUTS: value1 = Q, value2 = T.
            Q = value1;
            T_sat = value2;
            // p_sat from the SuperAncillary's pressure expansion.
            p_sat = sa->eval_sat(T_sat, 'P', 0);
        }
        if (Q < 0.0 || Q > 1.0) {
            throw ValueError("SVDSBTL backend: two-phase Q must be in [0, 1]");
        }
        update_two_phase_(T_sat, p_sat, Q);
        return;
    }

    const auto it = surfaces_.find(static_cast<int>(input_pair));
    if (it == surfaces_.end() || !it->second) {
        throw ValueError(std::string("SVDSBTL backend: no surface registered for this input pair (") + CoolProp::get_input_pair_short_desc(input_pair)
                         + ")");
    }

    active_pair_ = input_pair;
    active_surface_ = it->second.get();

    // Record inputs into AbstractState's standard slots so calc_T() /
    // calc_pressure() can return them directly without an SVD round-trip.
    switch (input_pair) {
        case CoolProp::HmassP_INPUTS: {
            // (a, b) = (p, h) per ph_subcritical preset; value1=h, value2=p.
            _p = value2;
            active_hmass_ = value1;
            active_point_ = active_surface_->resolve(_p, active_hmass_);
            break;
        }
        case CoolProp::HmolarP_INPUTS: {
            // Convert molar → mass and route to the PH surface.
            _p = value2;
            active_hmass_ = value1 / calc_molar_mass();
            active_point_ = active_surface_->resolve(_p, active_hmass_);
            break;
        }
        case CoolProp::PT_INPUTS: {
            // (a, b) = (p, T) per pt_subcritical preset; value1=p, value2=T.
            _p = value1;
            _T = value2;
            active_point_ = active_surface_->resolve(_p, _T);
            break;
        }
        default:
            throw ValueError(std::string("SVDSBTL backend update: unsupported input pair (") + CoolProp::get_input_pair_short_desc(input_pair) + ")");
    }
    active_resolved_ = true;

    if (active_point_.region_idx < 0) {
        // Out of every single-phase region.  For HmassP this typically
        // means the (h, p) point is inside the saturation dome; route
        // to the superancillary two-phase path and recover the user's
        // request.  Other input pairs (PT on the sat curve) currently
        // still surface as NaN -- two-phase PT is degenerate anyway
        // (no unique state at T_sat unless Q is supplied).
        if ((input_pair == CoolProp::HmassP_INPUTS || input_pair == CoolProp::HmolarP_INPUTS) && _p > 0.0) {
            auto sa = superanc_();
            if (sa) {
                const double T_sat = sa->get_T_from_p(_p);
                const double hL_mol = sa->eval_sat(T_sat, 'H', 0);
                const double hV_mol = sa->eval_sat(T_sat, 'H', 1);
                const double h_mol = active_hmass_ * calc_molar_mass();
                if (hL_mol <= h_mol && h_mol <= hV_mol) {
                    const double Q = (h_mol - hL_mol) / (hV_mol - hL_mol);
                    update_two_phase_(T_sat, _p, Q);
                    return;
                }
            }
        }
        _phase = iphase_twophase;  // best-effort tag; users typically check via the NaN return
    }
}

void SVDSBTLBackend::update_two_phase_(double T_sat, double p_sat, double Q) {
    auto sa = superanc_();
    if (!sa) {
        throw NotImplementedError("SVDSBTL backend: update_two_phase_ called without a SuperAncillary");
    }
    _p = p_sat;
    _T = T_sat;
    _Q = Q;
    _phase = iphase_twophase;
    two_phase_.rhoL_mol = sa->eval_sat(T_sat, 'D', 0);
    two_phase_.rhoV_mol = sa->eval_sat(T_sat, 'D', 1);
    two_phase_.hL_mol = sa->eval_sat(T_sat, 'H', 0);
    two_phase_.hV_mol = sa->eval_sat(T_sat, 'H', 1);
    two_phase_.sL_mol = sa->eval_sat(T_sat, 'S', 0);
    two_phase_.sV_mol = sa->eval_sat(T_sat, 'S', 1);
    // The SuperAncillary's caloric lazy-build covers H and S but not U.
    // Recover u via the thermodynamic identity u = h - p*v = h - p/rho.
    // (Equivalent in mass or molar basis; the SuperAncillary is molar.)
    two_phase_.uL_mol = two_phase_.hL_mol - p_sat / two_phase_.rhoL_mol;
    two_phase_.uV_mol = two_phase_.hV_mol - p_sat / two_phase_.rhoV_mol;
    two_phase_.active = true;
}

CoolPropDbl SVDSBTLBackend::two_phase_property_(CoolProp::parameters prop) const {
    // Specific volume blends linearly; density is its reciprocal.
    // Caloric props blend linearly directly.  Everything here is molar.
    const double Q = static_cast<double>(_Q);
    switch (prop) {
        case CoolProp::iDmolar: {
            const double vL = 1.0 / two_phase_.rhoL_mol;
            const double vV = 1.0 / two_phase_.rhoV_mol;
            return 1.0 / ((1.0 - Q) * vL + Q * vV);
        }
        case CoolProp::iHmolar:
            return (1.0 - Q) * two_phase_.hL_mol + Q * two_phase_.hV_mol;
        case CoolProp::iSmolar:
            return (1.0 - Q) * two_phase_.sL_mol + Q * two_phase_.sV_mol;
        case CoolProp::iUmolar:
            return (1.0 - Q) * two_phase_.uL_mol + Q * two_phase_.uV_mol;
        default:
            throw ValueError(std::string("SVDSBTL backend: two-phase blend not implemented for property '")
                             + CoolProp::get_parameter_information(prop, "short") + "'");
    }
}

CoolPropDbl SVDSBTLBackend::lookup_(CoolProp::parameters prop) {
    if (active_surface_ == nullptr || !active_resolved_) {
        throw ValueError("SVDSBTL backend: lookup before update()");
    }
    if (active_point_.region_idx < 0) {
        return std::nan("");
    }
    if (!active_surface_->contains_property(prop)) {
        throw ValueError(std::string("SVDSBTL backend: active surface does not expose property '")
                         + CoolProp::get_parameter_information(prop, "short") + "'");
    }
    return active_surface_->eval_with_region(prop, active_point_.region_idx, active_point_.svd_x, active_point_.svd_y);
}

// Mass-basis accessors — route through the active SVDSurface, OR
// through the cached two-phase blend when the active state is inside
// the saturation dome.
CoolPropDbl SVDSBTLBackend::calc_rhomass() {
    if (two_phase_.active) {
        return two_phase_property_(iDmolar) * calc_molar_mass();
    }
    return lookup_(iDmass);
}
CoolPropDbl SVDSBTLBackend::calc_hmass() {
    if (two_phase_.active) {
        return two_phase_property_(iHmolar) / calc_molar_mass();
    }
    // If the user supplied h via HmassP_INPUTS, return it directly; no
    // SVD round-trip needed.
    if (active_pair_ == CoolProp::HmassP_INPUTS || active_pair_ == CoolProp::HmolarP_INPUTS) {
        return active_hmass_;
    }
    return lookup_(iHmass);
}
CoolPropDbl SVDSBTLBackend::calc_smass() {
    if (two_phase_.active) {
        return two_phase_property_(iSmolar) / calc_molar_mass();
    }
    return lookup_(iSmass);
}
CoolPropDbl SVDSBTLBackend::calc_umass() {
    if (two_phase_.active) {
        return two_phase_property_(iUmolar) / calc_molar_mass();
    }
    return lookup_(iUmass);
}

CoolPropDbl SVDSBTLBackend::calc_T() {
    if (two_phase_.active || active_pair_ == CoolProp::PT_INPUTS) {
        return _T;
    }
    return lookup_(iT);
}
CoolPropDbl SVDSBTLBackend::calc_speed_sound() {
    if (two_phase_.active) {
        // Speed of sound is not uniquely defined in the two-phase
        // region (it goes to zero at the interface; equilibrium
        // thermodynamics doesn't give a meaningful bulk value).
        // Mirror what HEOS does -- throw rather than return a
        // misleading Q-weighted blend.
        throw NotImplementedError("SVDSBTL backend: speed_sound is not defined in the two-phase region");
    }
    return lookup_(ispeed_sound);
}
CoolPropDbl SVDSBTLBackend::calc_pressure() {
    return _p;
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

}  // namespace CoolProp
