#ifndef FLUIDLIBRARYFACTORIES_H
#define FLUIDLIBRARYFACTORIES_H

// NON-INSTALLED factory helpers that build fluid-model objects from
// nlohmann::json. These live under src/ (not include/) so that the
// nlohmann::json type never appears in a public/installed header.
// Neither factory requires friendship: make_surface_tension_correlation uses
// only public members, and make_saturation_ancillary builds a plain-typed
// SaturationAncillaryFunction::Values bundle and constructs from it.
//
// Convention: every value type constructed from JSON exposes a nested public
// `struct Values` (plain/POD fields) plus a single `explicit T(const Values&)`
// ctor; these factories own all nlohmann parsing and hand a Values across, so
// no JSON type appears in any installed header. (This is distinct from the
// link-time export control via cmake/CoolPropJSONVisibility.cmake that hides
// nlohmann/valijson symbols per shared product at link time.)

#include "CoolProp/detail/json.h"
#include "CoolProp/fluids/Ancillaries.h"

namespace cpjson {

/// Build a SurfaceTensionCorrelation from its JSON description. Parses into a
/// plain-typed SurfaceTensionCorrelation::Values bundle and constructs from
/// it, so no nlohmann type crosses the installed-header boundary.
inline CoolProp::SurfaceTensionCorrelation make_surface_tension_correlation(const nlohmann::json& j) {
    CoolProp::SurfaceTensionCorrelation::Values vals;
    vals.a = cpjson::get_long_double_array(j.at("a"));
    vals.n = cpjson::get_long_double_array(j.at("n"));
    if (vals.a.size() != vals.n.size()) {
        throw CoolProp::ValueError("Surface tension 'a' and 'n' arrays must have equal length");
    }
    vals.Tc = cpjson::get_double(j, "Tc");
    vals.BibTeX = cpjson::get_string(j, "BibTeX");
    return CoolProp::SurfaceTensionCorrelation(vals);
}

/// Build a SaturationAncillaryFunction from its JSON description. Parses into a
/// plain-typed SaturationAncillaryFunction::Values bundle and constructs from
/// it, so no nlohmann type crosses the installed-header boundary.
inline CoolProp::SaturationAncillaryFunction make_saturation_ancillary(const nlohmann::json& j) {
    CoolProp::SaturationAncillaryFunction::Values vals;
    std::string type = cpjson::get_string(j, "type");
    if (!type.compare("rational_polynomial")) {
        vals.type = CoolProp::SaturationAncillaryFunction::TYPE_RATIONAL_POLYNOMIAL;
        vals.num_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("A")));
        vals.den_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("B")));
        vals.max_abs_error = cpjson::get_double(j, "max_abs_error");
        try {
            vals.Tmin = cpjson::get_double(j, "Tmin");
            vals.Tmax = cpjson::get_double(j, "Tmax");
        } catch (...) {
            vals.Tmin = _HUGE;
            vals.Tmax = _HUGE;
        }
    } else {
        if (!type.compare("rhoLnoexp"))
            vals.type = CoolProp::SaturationAncillaryFunction::TYPE_NOT_EXPONENTIAL;
        else
            vals.type = CoolProp::SaturationAncillaryFunction::TYPE_EXPONENTIAL;
        vals.n = cpjson::get_double_array(j.at("n"));
        vals.t = cpjson::get_double_array(j.at("t"));
        if (vals.n.size() != vals.t.size()) {
            throw CoolProp::ValueError("Ancillary 'n' and 't' arrays must have equal length");
        }
        vals.Tmin = cpjson::get_double(j, "Tmin");
        vals.Tmax = cpjson::get_double(j, "Tmax");
        vals.reducing_value = cpjson::get_double(j, "reducing_value");
        vals.using_tau_r = cpjson::get_bool(j, "using_tau_r");
        vals.T_r = cpjson::get_double(j, "T_r");
    }
    return CoolProp::SaturationAncillaryFunction(vals);
}

}  // namespace cpjson

#endif  // FLUIDLIBRARYFACTORIES_H
