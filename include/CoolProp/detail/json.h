#ifndef COOLPROP_DETAIL_JSON_H
#define COOLPROP_DETAIL_JSON_H

// Internal, NON-INSTALLED wrapper around nlohmann/json.
//
// Symbol-leak strategy (see migration spec §4): nlohmann is included under
// hidden ELF/Mach-O visibility so none of its (weak/inline) symbols are
// exported from CoolProp's shared products. Combined with nlohmann's own
// versioned inline namespace (nlohmann::json_abi_v3_x_x), this prevents the
// cross-library ODR/symbol clashes that motivated the migration. Do NOT
// rename nlohmann's namespace here: Valijson's nlohmann adapter references
// nlohmann::json literally and a rename would break it.

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/tools.h"  // for CoolProp::format / CoolPropDbl

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility push(hidden)
#endif

#include <nlohmann/json.hpp>

#if defined(__GNUC__) || defined(__clang__)
#    pragma GCC visibility pop
#endif

namespace cpjson {

/// Parse a JSON-formatted string into an nlohmann::json document.
/// Throws CoolProp::ValueError (never a raw nlohmann exception) on failure.
inline nlohmann::json parse(std::string_view text) {
    try {
        return nlohmann::json::parse(text);
    } catch (const std::exception& e) {
        throw CoolProp::ValueError(std::string("Unable to parse JSON: ") + e.what());
    }
}

/// Serialize an nlohmann::json value to a pretty-printed string.
inline std::string json2string(const nlohmann::json& v) {
    return v.dump(4);
}

inline int get_integer(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_number_integer()) throw CoolProp::ValueError(format("Member [%s] is not an integer", m.c_str()));
    const std::int64_t val = it->get<std::int64_t>();
    if (val < std::numeric_limits<int>::min() || val > std::numeric_limits<int>::max())
        throw CoolProp::ValueError(format("Member [%s] is out of int range", m.c_str()));
    return static_cast<int>(val);
}

inline double get_double(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_number()) throw CoolProp::ValueError(format("Member [%s] is not a number", m.c_str()));
    return it->get<double>();
}

inline bool get_bool(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_boolean()) throw CoolProp::ValueError(format("Member [%s] is not a boolean", m.c_str()));
    return it->get<bool>();
}

inline std::string get_string(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    if (!it->is_string()) throw CoolProp::ValueError(format("Member [%s] is not a string", m.c_str()));
    return it->get<std::string>();
}

inline std::vector<double> get_double_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<double> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_number()) throw CoolProp::ValueError("input is not a number");
        out.push_back(el.get<double>());
    }
    return out;
}

inline std::vector<double> get_double_array(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    return get_double_array(*it);
}

inline std::vector<CoolPropDbl> get_long_double_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<CoolPropDbl> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_number()) throw CoolProp::ValueError("input is not a number");
        out.push_back(static_cast<CoolPropDbl>(el.get<double>()));
    }
    return out;
}

inline std::vector<CoolPropDbl> get_long_double_array(const nlohmann::json& v, const std::string& name) {
    auto it = v.find(name);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", name.c_str()));
    return get_long_double_array(*it);
}

inline std::vector<std::vector<double>> get_double_array2D(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<std::vector<double>> out;
    for (const auto& row : v) {
        if (!row.is_array()) throw CoolProp::ValueError(format("input \"%s\" is not a 2D array", json2string(v).c_str()));
        out.push_back(get_double_array(row));
    }
    return out;
}

inline std::vector<std::vector<CoolPropDbl>> get_long_double_array2D(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not a 2D array");
    std::vector<std::vector<CoolPropDbl>> out;
    for (const auto& row : v) {
        if (!row.is_array()) throw CoolProp::ValueError("input is not a 2D array");
        out.push_back(get_long_double_array(row));
    }
    return out;
}

inline std::vector<std::string> get_string_array(const nlohmann::json& v) {
    if (!v.is_array()) throw CoolProp::ValueError("input is not an array");
    std::vector<std::string> out;
    out.reserve(v.size());
    for (const auto& el : v) {
        if (!el.is_string()) throw CoolProp::ValueError("input is not a string");
        out.push_back(el.get<std::string>());
    }
    return out;
}

inline std::vector<std::string> get_string_array(const nlohmann::json& v, const std::string& m) {
    auto it = v.find(m);
    if (it == v.end()) throw CoolProp::ValueError(format("Does not have member [%s]", m.c_str()));
    return get_string_array(*it);
}

}  // namespace cpjson

#endif  // COOLPROP_DETAIL_JSON_H
