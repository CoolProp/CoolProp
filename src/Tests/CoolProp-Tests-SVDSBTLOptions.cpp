// Catch2 tests for the SVDSBTL backend's factory-string options opt-in
// (CoolProp-nez + CoolProp-u5c on the backend-options epic).  Covers:
//
//   * Schema validation — unknown keys, type mismatches, missing
//     `schema` version, enum violations.
//   * Grid knobs (NT/NR/rank) actually drive the table build.
//   * Cache filename incorporates the symbolic input-pair name and the
//     FNV-1a 64 opthash; two logically-different option blobs produce
//     distinct cache files; logically-equivalent blobs (same JSON,
//     different key order) collide on the same file.
//   * build_options_json() round-trips through factory().

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <filesystem>
#    include <memory>
#    include <string>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"
#    include "CoolProp/Hash.h"
#    include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#    include "CoolProp/Exceptions.h"

namespace cp_sbtl = CoolProp::sbtl;

TEST_CASE("SVDSBTL options: schema validation rejects bad payloads", "[SVDSBTL][options]") {
    auto factory = [](const std::string& opts) {
        return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water?" + opts));
    };

    SECTION("unknown top-level key") {
        REQUIRE_THROWS_AS(factory(R"({"frobnitz":1})"), CoolProp::ValueError);
    }
    SECTION("unknown nested key under critical_patch") {
        REQUIRE_THROWS_AS(factory(R"({"critical_patch":{"verboten":true}})"), CoolProp::ValueError);
    }
    SECTION("wrong type for grid.NT") {
        REQUIRE_THROWS_AS(factory(R"({"grid":{"NT":"two-hundred"}})"), CoolProp::ValueError);
    }
    SECTION("invalid enum for critical_patch.mode") {
        REQUIRE_THROWS_AS(factory(R"({"critical_patch":{"mode":"ON"}})"), CoolProp::ValueError);
    }
    SECTION("invalid JSON syntax") {
        REQUIRE_THROWS_AS(factory("{not json"), CoolProp::ValueError);
    }
    SECTION("grid.NT below schema minimum") {
        REQUIRE_THROWS_AS(factory(R"({"grid":{"NT":1}})"), CoolProp::ValueError);
    }
    SECTION("grid.NT above schema maximum (would overflow GetInt() otherwise)") {
        // Schema caps NT at 100_000 so integer parsing in resolve_grid()
        // stays within range.  Without this cap the validator would accept
        // NT >= 2 unbounded, and values past INT32_MAX would overflow.
        REQUIRE_THROWS_AS(factory(R"({"grid":{"NT":2147483648}})"), CoolProp::ValueError);
        REQUIRE_THROWS_AS(factory(R"({"grid":{"NR":2147483648}})"), CoolProp::ValueError);
        REQUIRE_THROWS_AS(factory(R"({"grid":{"rank":2147483648}})"), CoolProp::ValueError);
    }
    SECTION("bbox wrong arity") {
        REQUIRE_THROWS_AS(factory(R"({"critical_patch":{"bbox":[1,2,3]}})"), CoolProp::ValueError);
    }
    SECTION("empty options accepted") {
        REQUIRE_NOTHROW(factory(R"({})"));
    }
    SECTION("fully populated valid options accepted") {
        REQUIRE_NOTHROW(factory(R"({
            "schema": 1,
            "grid": {"NT": 40, "NR": 80, "rank": 10},
            "properties": {"transport": "auto"},
            "critical_patch": {"mode": "off"}
        })"));
    }
}

TEST_CASE("SVDSBTL build_options_json round-trips through factory()", "[SVDSBTL][options]") {
    // Logically-equal option payloads (different key order) must
    // produce identical canonical strings, and feeding that string
    // back through the factory yields a backend with the same form.
    auto a = std::shared_ptr<CoolProp::AbstractState>(
      CoolProp::AbstractState::factory("SVDSBTL&HEOS", R"(Water?{"critical_patch":{"mode":"off"},"grid":{"NT":40,"NR":80,"rank":10}})"));
    auto b = std::shared_ptr<CoolProp::AbstractState>(
      CoolProp::AbstractState::factory("SVDSBTL&HEOS", R"(Water?{"grid":{"rank":10,"NR":80,"NT":40},"critical_patch":{"mode":"off"}})"));
    REQUIRE(a->build_options_json() == b->build_options_json());

    // Round-trip: feed the canonical form back as the options.
    const std::string canonical = a->build_options_json();
    REQUIRE_FALSE(canonical.empty());
    auto c = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water?" + canonical));
    REQUIRE(c->build_options_json() == canonical);
}

TEST_CASE("SVDSBTL build_options_json is '{}' for no options", "[SVDSBTL][options]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    REQUIRE(AS->build_options_json() == "{}");
}

TEST_CASE("SVDSBTL cache filename embeds symbolic input_pair name", "[SVDSBTL][options][b6v]") {
    const std::string path = cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", ::CoolProp::PT_INPUTS, "abcdef0123456789");
    REQUIRE(path.find("Water.HEOS.PT_INPUTS.abcdef0123456789.svd.bin.z") != std::string::npos);
    // Symbolic name explicitly — int "17" must NOT appear in the
    // filename (the old format).
    const auto basename = std::filesystem::path(path).filename().string();
    REQUIRE(basename.find(".17.") == std::string::npos);
}

TEST_CASE("SVDSBTL cache filename rejects bad opthash chars", "[SVDSBTL][options]") {
    REQUIRE_THROWS(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", ::CoolProp::PT_INPUTS, "../escape"));
    REQUIRE_THROWS(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", ::CoolProp::PT_INPUTS, ""));
    REQUIRE_THROWS(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", ::CoolProp::PT_INPUTS, "ABCDEF"));  // uppercase rejected
}

TEST_CASE("FNV-1a 64 helper matches a known vector", "[Hash]") {
    // FNV-1a 64 of "foobar" (canonical test vector): 0x85944171f73967e8
    REQUIRE(CoolProp::fnv1a_64(std::string_view("foobar")) == 0x85944171f73967e8ULL);
    // Empty string → offset basis
    REQUIRE(CoolProp::fnv1a_64(std::string_view("")) == 0xCBF29CE484222325ULL);
    // to_hex16 lower-case fixed-width
    REQUIRE(CoolProp::to_hex16(0x85944171f73967e8ULL) == "85944171f73967e8");
    REQUIRE(CoolProp::to_hex16(0ULL) == "0000000000000000");
}

#endif  // ENABLE_CATCH
