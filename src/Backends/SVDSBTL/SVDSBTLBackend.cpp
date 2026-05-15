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
CoolProp::sbtl::SVDSurface build_surface_for_pair_(::CoolProp::AbstractState& heos, CoolProp::input_pairs pair) {
    namespace cp_sbtl = CoolProp::sbtl;
    cp_sbtl::SurfaceSpec spec;
    switch (pair) {
        case CoolProp::HmassP_INPUTS:
            spec = cp_sbtl::presets::ph_subcritical(heos);
            break;
        case CoolProp::PT_INPUTS:
            spec = cp_sbtl::presets::pt_subcritical(heos);
            break;
        default:
            throw ValueError("SVDSBTL backend: no preset registered for the requested input pair");
    }
    return cp_sbtl::build_surface(heos, std::move(spec));
}

}  // namespace

SVDSBTLBackend::SVDSBTLBackend(const std::string& fluid_name) : fluid_name_(fluid_name), mole_fractions_({1.0}) {
    for (const auto pair : kSupportedPairs) {
        ensure_surface_(pair);
    }
}

std::shared_ptr<CoolProp::AbstractState> SVDSBTLBackend::heos_reference_() {
    if (!heos_) {
        heos_.reset(CoolProp::AbstractState::factory("HEOS", fluid_name_));
    }
    return heos_;
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
    const std::string path = cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid_name_, pair);
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
        auto heos = heos_reference_();
        surface = std::make_unique<cp_sbtl::SVDSurface>(build_surface_for_pair_(*heos, pair));
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
        // Point outside every region (two-phase or out-of-bounds).
        _phase = iphase_twophase;  // best-effort tag; users typically check via the NaN return
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

// Mass-basis accessors — route through the active SVDSurface.
CoolPropDbl SVDSBTLBackend::calc_rhomass() {
    return lookup_(iDmass);
}
CoolPropDbl SVDSBTLBackend::calc_hmass() {
    // If the user supplied h via HmassP_INPUTS, return it directly; no
    // SVD round-trip needed.
    if (active_pair_ == CoolProp::HmassP_INPUTS || active_pair_ == CoolProp::HmolarP_INPUTS) {
        return active_hmass_;
    }
    return lookup_(iHmass);
}
CoolPropDbl SVDSBTLBackend::calc_smass() {
    return lookup_(iSmass);
}
CoolPropDbl SVDSBTLBackend::calc_umass() {
    return lookup_(iUmass);
}

CoolPropDbl SVDSBTLBackend::calc_T() {
    if (active_pair_ == CoolProp::PT_INPUTS) {
        return _T;
    }
    return lookup_(iT);
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
    return heos_reference_()->molar_mass();
}
CoolPropDbl SVDSBTLBackend::calc_gas_constant() {
    return heos_reference_()->gas_constant();
}
CoolPropDbl SVDSBTLBackend::calc_T_critical() {
    return heos_reference_()->T_critical();
}
CoolPropDbl SVDSBTLBackend::calc_p_critical() {
    return heos_reference_()->p_critical();
}
CoolPropDbl SVDSBTLBackend::calc_rhomass_critical() {
    return heos_reference_()->rhomass_critical();
}
CoolPropDbl SVDSBTLBackend::calc_rhomolar_critical() {
    return heos_reference_()->rhomolar_critical();
}
CoolPropDbl SVDSBTLBackend::calc_Ttriple() {
    return heos_reference_()->Ttriple();
}
CoolPropDbl SVDSBTLBackend::calc_p_triple() {
    return heos_reference_()->p_triple();
}
CoolPropDbl SVDSBTLBackend::calc_Tmax() {
    return heos_reference_()->Tmax();
}
CoolPropDbl SVDSBTLBackend::calc_Tmin() {
    return heos_reference_()->Tmin();
}
CoolPropDbl SVDSBTLBackend::calc_pmax() {
    return heos_reference_()->pmax();
}

}  // namespace CoolProp
