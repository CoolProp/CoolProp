/// Out-of-line, JSON-consuming construction for SuperAncillary.
///
/// The installed header include/CoolProp/superancillary/superancillary.h keeps all
/// templated, perf-critical evaluation inline but carries NO nlohmann::json type, so
/// consumers do not need nlohmann to compile against CoolProp's public headers. The
/// string-blob constructor's definition lives here (where nlohmann is available) and is
/// explicitly instantiated for the one ArrayType used across the codebase.

#include "CoolProp/superancillary/superancillary.h"
#include "nlohmann/json.hpp"
#include "boost/math/tools/toms748_solve.hpp"

#include <functional>
#include <string>
#include <vector>

namespace CoolProp {
namespace superancillary {

namespace detail {

// Out-of-line definition of the TOMS 748 bracketing rootfinder declared in the
// installed header.  Keeping the boost dependency here (rather than in the
// header) lets downstream consumers compile against superancillary.h with only
// Eigen on the include path, not boost (CoolProp-1tbe.14).  The callable is
// type-erased through std::function: this is the inverse (rootfinding) path,
// NOT the hot eval_sat path, and the std::function is built once per rootfind,
// so its indirection -- plus a possible one-time small-buffer heap allocation
// for the larger residual closures (they exceed libc++'s 16-byte SBO) -- is
// negligible against the tens of Chebyshev evaluations each rootfind performs.
double toms748(const std::function<double(double)>& f, double a, double b, double fa, double fb, unsigned int bits, std::size_t max_iter) {
    using namespace boost::math::tools;
    auto max_iter_ = static_cast<boost::math::uintmax_t>(max_iter);
    auto [l, r] = toms748_solve(f, a, b, fa, fb, eps_tolerance<double>(bits), max_iter_);
    return (l + r) / 2.0;
}

}  // namespace detail

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
