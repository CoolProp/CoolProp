#ifndef COOLPROP_DETAIL_JSON_H
#define COOLPROP_DETAIL_JSON_H

// Internal, NON-INSTALLED wrapper around nlohmann/json.
//
// Symbol-leak strategy: nlohmann and valijson symbols are hidden at LINK time
// per shared product via cmake/CoolPropJSONVisibility.cmake, which applies an
// ELF --version-script (Linux/ELF) or a Mach-O -unexported_symbols_list
// (macOS) hide-list to every shared product. This replaces the former
// compile-time GCC visibility pragma (which poisoned __assert_fail in
// non-NDEBUG loader builds). Combined with nlohmann's own versioned inline
// namespace (nlohmann::json_abi_v3_x_x), this prevents cross-library
// ODR/symbol clashes. Do NOT rename nlohmann's namespace here: Valijson's
// nlohmann adapter references nlohmann::json literally and a rename would
// break it.

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/tools.h"  // for CoolProp::format / CoolPropDbl

#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>
#include <valijson/adapters/nlohmann_json_adapter.hpp>
#include <valijson/schema.hpp>
#include <valijson/schema_parser.hpp>
#include <valijson/validator.hpp>

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

/// Decode a CBOR byte buffer into an nlohmann::json document.
/// Throws CoolProp::ValueError (never a raw nlohmann exception) on failure.
inline nlohmann::json from_cbor(const std::uint8_t* data, std::size_t size) {
    try {
        return nlohmann::json::from_cbor(data, data + size);
    } catch (const std::exception& e) {
        throw CoolProp::ValueError(std::string("Unable to decode CBOR: ") + e.what());
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
    // is_number_integer() is true for unsigned values too; a uint64 above
    // INT64_MAX would wrap to a negative int64 and slip past the range check,
    // so handle the unsigned case explicitly before the signed one.
    if (it->is_number_unsigned()) {
        const std::uint64_t uval = it->get<std::uint64_t>();
        if (uval > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
            throw CoolProp::ValueError(format("Member [%s] is out of int range", m.c_str()));
        return static_cast<int>(uval);
    }
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

/// Result of validate_schema() below.
enum schema_validation_code
{
    SCHEMA_VALIDATION_OK = 0,
    SCHEMA_INVALID_JSON,
    INPUT_INVALID_JSON,
    SCHEMA_NOT_VALIDATED
};

/// Validate a JSON-formatted input string against a JSON-formatted draft-07
/// schema string. On a validation failure, `errstr` receives a
/// human-comprehensible description. Never propagates a raw nlohmann or
/// Valijson exception to the caller.
inline schema_validation_code validate_schema(std::string_view schemaJson, std::string_view inputJson, std::string& errstr) {
    nlohmann::json schemaDoc;
    try {
        schemaDoc = nlohmann::json::parse(schemaJson);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid schema: ") + e.what();
        return SCHEMA_INVALID_JSON;
    }

    nlohmann::json inputDoc;
    try {
        inputDoc = nlohmann::json::parse(inputJson);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid input json: ") + e.what();
        return INPUT_INVALID_JSON;
    }

    valijson::Schema schema;
    valijson::SchemaParser parser;
    valijson::adapters::NlohmannJsonAdapter schemaAdapter(schemaDoc);
    try {
        parser.populateSchema(schemaAdapter, schema);
    } catch (const std::exception& e) {
        errstr = std::string("Invalid schema: ") + e.what();
        return SCHEMA_INVALID_JSON;
    }

    valijson::Validator validator;
    valijson::ValidationResults results;
    valijson::adapters::NlohmannJsonAdapter inputAdapter(inputDoc);
    try {
        if (!validator.validate(schema, inputAdapter, &results)) {
            std::string msg;
            valijson::ValidationResults::Error error;
            while (results.popError(error)) {
                for (const std::string& ctx : error.context)
                    msg += ctx;
                msg += ": " + error.description + "\n";
            }
            errstr = msg;
            return SCHEMA_NOT_VALIDATED;
        }
    } catch (const std::exception& e) {
        errstr = std::string("Schema validation internal error: ") + e.what();
        return SCHEMA_NOT_VALIDATED;
    }
    return SCHEMA_VALIDATION_OK;
}

}  // namespace cpjson

#endif  // COOLPROP_DETAIL_JSON_H
