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

#endif  // ENABLE_CATCH
