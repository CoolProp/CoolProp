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

}  // namespace cpjson

#endif  // COOLPROP_DETAIL_JSON_H
