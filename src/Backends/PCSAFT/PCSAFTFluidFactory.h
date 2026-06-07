#ifndef PCSAFTFLUIDFACTORY_H
#define PCSAFTFLUIDFACTORY_H

// NON-INSTALLED factory that builds a PCSAFTFluid from nlohmann::json. Lives
// under src/ so the nlohmann::json type never appears in the installed
// PCSAFTFluid.h. Mirrors the (now-removed) rapidjson constructor.

#include "CoolProp/detail/json.h"
#include "CoolProp/fluids/PCSAFTFluid.h"

namespace cpjson {

inline CoolProp::PCSAFTFluid make_pcsaft_fluid(const nlohmann::json& fluid) {
    CoolProp::PCSAFTValues params;
    params.m = cpjson::get_double(fluid, "m");
    params.sigma = cpjson::get_double(fluid, "sigma");
    params.u = cpjson::get_double(fluid, "u");

    params.uAB = (fluid.contains("uAB") && fluid.at("uAB").is_number()) ? cpjson::get_double(fluid, "uAB") : 0.;
    params.volA = (fluid.contains("volA") && fluid.at("volA").is_number()) ? cpjson::get_double(fluid, "volA") : 0.;
    params.assocScheme = fluid.contains("assocScheme") ? cpjson::get_string_array(fluid, "assocScheme") : std::vector<std::string>{};
    params.dipm = (fluid.contains("dipm") && fluid.at("dipm").is_number()) ? cpjson::get_double(fluid, "dipm") : 0.;
    params.dipnum = (fluid.contains("dipnum") && fluid.at("dipnum").is_number()) ? cpjson::get_double(fluid, "dipnum") : 0.;
    params.z = (fluid.contains("charge") && fluid.at("charge").is_number()) ? cpjson::get_double(fluid, "charge") : 0.;

    return CoolProp::PCSAFTFluid(cpjson::get_string(fluid, "name"), cpjson::get_string(fluid, "CAS"),
                                 cpjson::get_double(fluid, "molemass"), cpjson::get_string_array(fluid, "aliases"),
                                 std::move(params));
}

}  // namespace cpjson

#endif  // PCSAFTFLUIDFACTORY_H
