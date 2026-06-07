/// Out-of-line, JSON-consuming construction for SuperAncillary.
///
/// The installed header include/CoolProp/superancillary/superancillary.h keeps all
/// templated, perf-critical evaluation inline but carries NO nlohmann::json type, so
/// consumers do not need nlohmann to compile against CoolProp's public headers. The
/// string-blob constructor's definition lives here (where nlohmann is available) and is
/// explicitly instantiated for the one ArrayType used across the codebase.

#include "CoolProp/superancillary/superancillary.h"
#include "nlohmann/json.hpp"

#include <string>
#include <vector>

namespace CoolProp {
namespace superancillary {
namespace {

template <typename ArrayType>
std::vector<ChebyshevExpansion<ArrayType>> load_expansions(const nlohmann::json& j, const std::string& key) {
    std::vector<ChebyshevExpansion<ArrayType>> buf;
    for (auto& block : j.at(key)) {
        buf.emplace_back(block.at("xmin"), block.at("xmax"), block.at("coef"));
    }
    return buf;
}

template <typename ArrayType>
typename SuperAncillary<ArrayType>::LoadedData parse_loaded(const std::string& s) {
    nlohmann::json j = nlohmann::json::parse(s);
    typename SuperAncillary<ArrayType>::LoadedData d;
    d.rhoL = load_expansions<ArrayType>(j, "jexpansions_rhoL");
    d.rhoV = load_expansions<ArrayType>(j, "jexpansions_rhoV");
    d.p = load_expansions<ArrayType>(j, "jexpansions_p");
    d.Tcrit_num = j.at("meta").at("Tcrittrue / K");
    d.rhocrit_num = j.at("meta").at("rhocrittrue / mol/m^3");
    if (j.contains("check_points")) {
        for (const auto& pt : j.at("check_points")) {
            d.check_points.push_back({pt.at("T / K").get<double>(), pt.at("p(mp) / Pa").get<double>(),
                                      pt.at("rho'(mp) / mol/m^3").get<double>(), pt.at("rho''(mp) / mol/m^3").get<double>(),
                                      pt.at("p(SA)/p(mp)").get<double>(), pt.at("rho'(SA)/rho'(mp)").get<double>(),
                                      pt.at("rho''(SA)/rho''(mp)").get<double>()});
        }
    }
    return d;
}

}  // namespace

template <typename ArrayType>
SuperAncillary<ArrayType>::SuperAncillary(const std::string& s) : SuperAncillary(parse_loaded<ArrayType>(s)) {}

// The one and only ArrayType string-constructed across the codebase. Constructing
// SuperAncillary<other types> from a string is intentionally not provided.
template SuperAncillary<std::vector<double>>::SuperAncillary(const std::string&);

}  // namespace superancillary
}  // namespace CoolProp
