#ifndef FLUIDLIBRARYFACTORIES_H
#define FLUIDLIBRARYFACTORIES_H

// NON-INSTALLED factory helpers that build fluid-model objects from
// nlohmann::json. These live under src/ (not include/) so that the
// nlohmann::json type never appears in a public/installed header. The
// installed Ancillaries.h friends these free functions so they can populate
// the otherwise-private SaturationAncillaryFunction fields directly.

#include "CoolProp/detail/json.h"
#include "CoolProp/fluids/Ancillaries.h"

namespace cpjson {

/// Build a SurfaceTensionCorrelation from its JSON description. Mirrors the
/// (now-removed) inline rapidjson constructor that previously lived in
/// include/CoolProp/fluids/Ancillaries.h.
inline CoolProp::SurfaceTensionCorrelation make_surface_tension_correlation(const nlohmann::json& j) {
    CoolProp::SurfaceTensionCorrelation out;
    out.a = cpjson::get_long_double_array(j.at("a"));
    out.n = cpjson::get_long_double_array(j.at("n"));
    out.Tc = cpjson::get_double(j, "Tc");
    out.BibTeX = cpjson::get_string(j, "BibTeX");
    out.N = out.n.size();
    out.s = out.n;
    return out;
}

/// Build a SaturationAncillaryFunction from its JSON description. Reproduces
/// the (now-removed) SaturationAncillaryFunction(rapidjson::Value&)
/// constructor body with the nlohmann idiom map applied. SaturationAncillaryFunction
/// friends this function so it can write the class's private members.
inline CoolProp::SaturationAncillaryFunction make_saturation_ancillary(const nlohmann::json& j) {
    CoolProp::SaturationAncillaryFunction out;
    std::string type = cpjson::get_string(j, "type");
    if (!type.compare("rational_polynomial")) {
        out.type = CoolProp::SaturationAncillaryFunction::TYPE_RATIONAL_POLYNOMIAL;
        out.num_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("A")));
        out.den_coeffs = CoolProp::vec_to_eigen(cpjson::get_double_array(j.at("B")));
        out.max_abs_error = cpjson::get_double(j, "max_abs_error");
        try {
            out.Tmin = cpjson::get_double(j, "Tmin");
            out.Tmax = cpjson::get_double(j, "Tmax");
        } catch (...) {
            out.Tmin = _HUGE;
            out.Tmax = _HUGE;
        }
    } else {
        if (!type.compare("rhoLnoexp"))
            out.type = CoolProp::SaturationAncillaryFunction::TYPE_NOT_EXPONENTIAL;
        else
            out.type = CoolProp::SaturationAncillaryFunction::TYPE_EXPONENTIAL;
        out.n = cpjson::get_double_array(j.at("n"));
        out.N = out.n.size();
        out.s = out.n;
        out.t = cpjson::get_double_array(j.at("t"));
        out.Tmin = cpjson::get_double(j, "Tmin");
        out.Tmax = cpjson::get_double(j, "Tmax");
        out.reducing_value = cpjson::get_double(j, "reducing_value");
        out.using_tau_r = cpjson::get_bool(j, "using_tau_r");
        out.T_r = cpjson::get_double(j, "T_r");
    }
    return out;
}

}  // namespace cpjson

#endif  // FLUIDLIBRARYFACTORIES_H
