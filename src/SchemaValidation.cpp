#include "CoolProp/SchemaValidation.h"

#include <cstddef>
#include <string>

#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/json.h"

namespace CoolProp {

void validate_json_against_schema(const std::string& instance_json, const std::string& schema_json) {
    // The instance is the options blob of a factory string: caller-supplied text that holds a handful of
    // scalars.  Refuse anything larger before it is expanded into a DOM.  nlohmann::json spends 16 bytes on
    // a two-byte array item, and its noexcept destructor allocates again to flatten the tree, so running
    // out of memory there ends in std::terminate, which no catch block at the API boundary can stop.
    constexpr std::size_t max_instance_bytes = std::size_t{64} * 1024;
    if (instance_json.size() > max_instance_bytes) {
        throw ValueError("schema validation failed: the input is " + std::to_string(instance_json.size()) + " bytes long, the limit is "
                         + std::to_string(max_instance_bytes) + " bytes");
    }
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
