// Catch2 tests for the nlohmann-backed cpjson helpers and Valijson schema
// validation (RapidJSON→nlohmann migration, Phase 0).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <string>
#    include <vector>

#    include "CoolProp/detail/json.h"
#    include "CoolProp/Exceptions.h"

TEST_CASE("cpjson nlohmann wrapper parses a string", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"a": 1.5})");
    REQUIRE(j.at("a").get<double>() == Catch::Approx(1.5));
}

TEST_CASE("cpjson::parse throws ValueError on malformed JSON", "[json]") {
    REQUIRE_THROWS_AS(cpjson::parse("{not json}"), CoolProp::ValueError);
}

TEST_CASE("cpjson getters extract typed members", "[json]") {
    nlohmann::json j = cpjson::parse(R"({
        "i": 7, "x": 2.5, "b": true, "s": "hi",
        "darr": [1.0, 2.0, 3.0],
        "sarr": ["a", "b"],
        "d2d": [[1.0, 2.0], [3.0]]
    })");

    REQUIRE(cpjson::get_integer(j, "i") == 7);
    REQUIRE(cpjson::get_double(j, "x") == Catch::Approx(2.5));
    REQUIRE(cpjson::get_bool(j, "b") == true);
    REQUIRE(cpjson::get_string(j, "s") == "hi");

    std::vector<double> darr = cpjson::get_double_array(j, "darr");
    REQUIRE(darr.size() == 3);
    REQUIRE(darr[2] == Catch::Approx(3.0));

    std::vector<std::string> sarr = cpjson::get_string_array(j, "sarr");
    REQUIRE(sarr.size() == 2);
    REQUIRE(sarr[1] == "b");

    std::vector<std::vector<double>> d2d = cpjson::get_double_array2D(j.at("d2d"));
    REQUIRE(d2d.size() == 2);
    REQUIRE(d2d[0][1] == Catch::Approx(2.0));
}

TEST_CASE("cpjson getters throw ValueError on missing/mistyped members", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"x": "not a number"})");
    REQUIRE_THROWS_AS(cpjson::get_double(j, "missing"), CoolProp::ValueError);
    REQUIRE_THROWS_AS(cpjson::get_double(j, "x"), CoolProp::ValueError);
}

TEST_CASE("cpjson::get_integer rejects out-of-int-range values", "[json]") {
    // Signed value above INT_MAX but within int64.
    nlohmann::json j = cpjson::parse(R"({"big": 3000000000})");
    REQUIRE_THROWS_AS(cpjson::get_integer(j, "big"), CoolProp::ValueError);
    // Unsigned value above INT64_MAX — must NOT wrap to a negative int and slip
    // through; it has to throw.
    nlohmann::json u = cpjson::parse(R"({"huge": 18446744073709551615})");
    REQUIRE_THROWS_AS(cpjson::get_integer(u, "huge"), CoolProp::ValueError);
}

namespace {
constexpr const char* kJsonSchema = R"({
  "$schema": "http://json-schema.org/draft-07/schema#",
  "type": "object",
  "required": ["name", "T"],
  "properties": {
    "name": {"type": "string"},
    "T":    {"type": "number", "exclusiveMinimum": 0}
  }
})";
}

TEST_CASE("cpjson::validate_schema accepts a conforming document", "[json]") {
    std::string err;
    cpjson::schema_validation_code code = cpjson::validate_schema(kJsonSchema, R"({"name": "water", "T": 300.0})", err);
    REQUIRE(code == cpjson::SCHEMA_VALIDATION_OK);
    REQUIRE(err.empty());
}

TEST_CASE("cpjson::validate_schema rejects a non-conforming document with a message", "[json]") {
    std::string err;
    cpjson::schema_validation_code code = cpjson::validate_schema(kJsonSchema, R"({"name": "water", "T": -5.0})", err);
    REQUIRE(code == cpjson::SCHEMA_NOT_VALIDATED);
    REQUIRE_FALSE(err.empty());  // human-comprehensible error preserved
}

TEST_CASE("cpjson::validate_schema reports malformed input JSON, never leaks nlohmann exception", "[json]") {
    std::string err;
    cpjson::schema_validation_code code = cpjson::validate_schema(kJsonSchema, R"({not valid json)", err);
    REQUIRE(code == cpjson::INPUT_INVALID_JSON);
}

TEST_CASE("cpjson::validate_schema reports malformed schema JSON", "[json]") {
    std::string err;
    cpjson::schema_validation_code code = cpjson::validate_schema("{not a schema", R"({"name":"x","T":1})", err);
    REQUIRE(code == cpjson::SCHEMA_INVALID_JSON);
    REQUIRE_FALSE(err.empty());
}

TEST_CASE("cpjson::from_cbor round-trips a document", "[json]") {
    nlohmann::json j = cpjson::parse(R"({"a": 1.5, "b": [1, 2, 3]})");
    std::vector<std::uint8_t> blob = nlohmann::json::to_cbor(j);
    nlohmann::json back = cpjson::from_cbor(blob.data(), blob.size());
    REQUIRE(back == j);
}

TEST_CASE("cpjson::from_cbor throws ValueError on garbage", "[json]") {
    std::vector<std::uint8_t> bad = {0xFF, 0xFF, 0xFF};
    REQUIRE_THROWS_AS(cpjson::from_cbor(bad.data(), bad.size()), CoolProp::ValueError);
}

#endif  // ENABLE_CATCH
