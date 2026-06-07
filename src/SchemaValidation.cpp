#include "CoolProp/SchemaValidation.h"

#include <string>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/json.h"

namespace CoolProp {

void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json) {
    std::string errstr;
    // cpjson::validate_schema(schemaJson, inputJson, errstr) — Valijson-backed. Schema first.
    cpjson::schema_validation_code code = cpjson::validate_schema(schema_json, instance_json, errstr);
    if (code != cpjson::SCHEMA_VALIDATION_OK) {
        throw ValueError(std::string("schema validation failed: ") + errstr);
    }
}

std::string to_canonical_json_str(const std::string& json_str) {
    // nlohmann::json is std::map-backed, so object keys serialize sorted at every
    // level — the canonical (deterministic, key-sorted) form. Arrays keep order.
    return cpjson::parse(json_str).dump();
}

}  // namespace CoolProp
