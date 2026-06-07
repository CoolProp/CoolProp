#ifndef COOLPROP_SCHEMA_VALIDATION_H
#define COOLPROP_SCHEMA_VALIDATION_H

#include <string>

namespace CoolProp {

// Schema-validation helpers that work on JSON strings — exposed even
// to SWIG since they take only std::string.  Internally they parse,
// validate, and (for to_canonical_json_str) re-serialise into a
// canonical form with object keys sorted at every level.
//
// Throws CoolProp::ValueError on schema violation or invalid JSON.
void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json);

std::string to_canonical_json_str(const std::string& json_str);

}  // namespace CoolProp

#endif  // COOLPROP_SCHEMA_VALIDATION_H
