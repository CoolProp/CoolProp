// Catch2 tests for CoolProp::validate_json_against_schema and
// CoolProp::to_canonical_json_str (CoolProp-bh2).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <string>

#    include "CoolProp/SchemaValidation.h"
#    include "CoolProp/Exceptions.h"

namespace {

// Common test schema modelled after the SVDSBTL options schema we're
// about to write — exercises objects, nested objects, enums, types,
// required keys, and additionalProperties:false (strict mode).
constexpr const char* kSchema = R"({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "type": "object",
  "additionalProperties": false,
  "required": ["schema"],
  "properties": {
    "schema": {"type": "integer", "const": 1},
    "grid": {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "NT": {"type": "integer", "minimum": 2},
        "NR": {"type": "integer", "minimum": 2},
        "rank": {"type": "integer", "minimum": 1}
      }
    },
    "critical_patch": {
      "type": "object",
      "additionalProperties": false,
      "properties": {
        "mode": {"type": "string", "enum": ["auto", "off", "fixed"]},
        "tolerance": {"type": "number", "exclusiveMinimum": 0},
        "bbox": {
          "type": "array",
          "minItems": 4,
          "maxItems": 4,
          "items": {"type": "number"}
        }
      }
    }
  }
})";

}  // namespace

// ---------------------------------------------------------------------------
// validate_json_against_schema
// ---------------------------------------------------------------------------

TEST_CASE("validate_against_schema: minimal valid instance passes", "[SchemaValidation]") {
    REQUIRE_NOTHROW(CoolProp::validate_json_against_schema(R"({"schema": 1})", std::string(kSchema)));
}

TEST_CASE("validate_against_schema: fully populated valid instance passes", "[SchemaValidation]") {
    REQUIRE_NOTHROW(CoolProp::validate_json_against_schema(R"({
        "schema": 1,
        "grid": {"NT": 200, "NR": 800, "rank": 20},
        "critical_patch": {"mode": "auto", "tolerance": 1e-3, "bbox": [0.95, 1.05, 0.75, 1.15]}
    })",
                                                           std::string(kSchema)));
}

TEST_CASE("validate_against_schema: unknown top-level key throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1, "frobnitz": true})", std::string(kSchema)), CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: unknown nested key throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1, "grid": {"NT": 200, "extra": 3}})", std::string(kSchema)),
                      CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: type mismatch throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1, "grid": {"NT": "two-hundred"}})", std::string(kSchema)),
                      CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: enum violation throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1, "critical_patch": {"mode": "ON"}})", std::string(kSchema)),
                      CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: required key missing throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"grid": {"NT": 200}})", std::string(kSchema)), CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: bbox wrong size throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1, "critical_patch": {"bbox": [1, 2, 3]}})", std::string(kSchema)),
                      CoolProp::ValueError);
}

TEST_CASE("validate_against_schema: invalid schema JSON throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::validate_json_against_schema(R"({"schema": 1})", std::string("not json")), CoolProp::ValueError);
}

// ---------------------------------------------------------------------------
// to_canonical_json_str
// ---------------------------------------------------------------------------

TEST_CASE("to_canonical_json: sorts object keys", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"({"b": 2, "a": 1, "c": 3})");
    REQUIRE(a == R"({"a":1,"b":2,"c":3})");
}

TEST_CASE("to_canonical_json: sorts nested object keys recursively", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"({"outer": {"z": 1, "a": 2}, "inner": 3})");
    REQUIRE(a == R"({"inner":3,"outer":{"a":2,"z":1}})");
}

TEST_CASE("to_canonical_json: preserves array order", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"([3, 1, 2])");
    REQUIRE(a == R"([3,1,2])");
}

TEST_CASE("to_canonical_json: logically-equal inputs produce identical output", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"({"grid":{"NT":200,"rank":20,"NR":800}})");
    const auto b = CoolProp::to_canonical_json_str(R"({"grid":{"rank":20,"NR":800,"NT":200}})");
    REQUIRE(a == b);
}

TEST_CASE("to_canonical_json: idempotent — re-canonicalising returns the same bytes", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"({"b":2,"a":1})");
    const auto b = CoolProp::to_canonical_json_str(a);
    REQUIRE(a == b);
}

TEST_CASE("to_canonical_json: handles all JSON scalar types", "[SchemaValidation]") {
    const auto a = CoolProp::to_canonical_json_str(R"({"i":1,"f":1.5,"s":"x","b":true,"n":null,"a":[1,2]})");
    // Sorted: a, b, f, i, n, s
    REQUIRE(a == R"({"a":[1,2],"b":true,"f":1.5,"i":1,"n":null,"s":"x"})");
}

TEST_CASE("to_canonical_json: invalid JSON throws", "[SchemaValidation]") {
    REQUIRE_THROWS_AS(CoolProp::to_canonical_json_str(std::string("{not: json}")), CoolProp::ValueError);
}

#endif  // ENABLE_CATCH
