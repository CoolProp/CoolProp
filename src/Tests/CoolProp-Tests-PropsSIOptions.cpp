// Catch2 tests for the `?<options>` suffix flowing through the
// high-level CoolProp API (PropsSI, Props1SI, PhaseSI).  No backend
// opts in to options in PR A; these tests therefore verify that:
//
//   * empty / no-options strings work identically with or without the
//     `?` suffix (the suffix is parsed and dropped without affecting
//     the call);
//   * `?{...}` options on a backend that hasn't opted in throw a
//     clear NotImplementedError instead of being silently dropped.
//
// When PR B (SVDSBTL opt-in) lands, this file gains positive tests
// asserting that the SVDSBTL backend respects the supplied options.
// CoolProp-iqz.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <filesystem>
#    include <fstream>
#    include <string>

#    include "AbstractState.h"
#    include "CoolProp.h"
#    include "TestUtils.h"

namespace {

class TempJSONFile
{
   public:
    explicit TempJSONFile(const std::string& contents) {
        path_ = std::filesystem::temp_directory_path()
                / ("cp_propssi_opts_test_" + std::to_string(CoolProp::tests::test_pid()) + "_" + std::to_string(counter_++) + ".json");
        std::ofstream(path_) << contents;
    }
    ~TempJSONFile() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }
    TempJSONFile(const TempJSONFile&) = delete;
    TempJSONFile& operator=(const TempJSONFile&) = delete;
    TempJSONFile(TempJSONFile&&) = delete;
    TempJSONFile& operator=(TempJSONFile&&) = delete;

    [[nodiscard]] std::string path() const {
        return path_.string();
    }

   private:
    std::filesystem::path path_;
    static inline int counter_ = 0;
};

}  // namespace

TEST_CASE("PropsSI: '?<options>' suffix on fluid string parses + dispatches", "[FactoryOptions][PropsSI]") {
    const double T = 300.0;
    const double p = 101325.0;
    const double rho_ref = CoolProp::PropsSI("D", "T", T, "P", p, "HEOS::Water");
    REQUIRE(std::isfinite(rho_ref));

    SECTION("bare '?' is a no-op (empty options)") {
        const double rho = CoolProp::PropsSI("D", "T", T, "P", p, "HEOS::Water?");
        REQUIRE(rho == Catch::Approx(rho_ref));
    }
    SECTION("empty JSON object is a no-op") {
        const double rho = CoolProp::PropsSI("D", "T", T, "P", p, "HEOS::Water?{}");
        REQUIRE(rho == Catch::Approx(rho_ref));
    }
    // Trailing-whitespace tail is exercised by parse_factory_options
    // directly; PropsSI's string path normalises whitespace separately,
    // so don't double-test that integration here.
}

TEST_CASE("PropsSI: non-empty options on opted-out backend errors gracefully", "[FactoryOptions][PropsSI]") {
    // HEOS does not (yet) opt in to options; non-empty payload must
    // surface as a non-finite result, not a silent success at the
    // default settings.  PropsSI catches the exception and returns
    // HUGE_VAL; the errstring will carry the original message.
    const double rho = CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, R"(HEOS::Water?{"key":1})");
    REQUIRE_FALSE(std::isfinite(rho));
    const std::string err = CoolProp::get_global_param_string("errstring");
    REQUIRE_FALSE(err.empty());
}

TEST_CASE("PropsSI: '@path' indirection reads file on fluid string", "[FactoryOptions][PropsSI]") {
    TempJSONFile f("{}");
    const std::string fluid_arg = "HEOS::Water?@" + f.path();
    const double rho = CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, fluid_arg);
    REQUIRE(std::isfinite(rho));
    // Same as the no-options call.
    const double rho_ref = CoolProp::PropsSI("D", "T", 300.0, "P", 101325.0, "HEOS::Water");
    REQUIRE(rho == Catch::Approx(rho_ref));
}

TEST_CASE("Props1SI / PhaseSI: '?<options>' suffix accepted through high-level API", "[FactoryOptions][PropsSI]") {
    SECTION("Props1SI empty options") {
        const double tcrit = CoolProp::Props1SI("HEOS::Water?{}", "Tcrit");
        REQUIRE(tcrit == Catch::Approx(647.096).margin(0.5));
    }
    SECTION("PhaseSI empty options") {
        const std::string phase = CoolProp::PhaseSI("T", 300.0, "P", 101325.0, "HEOS::Water?{}");
        REQUIRE(phase == "liquid");
    }
}

TEST_CASE("AbstractState::factory: ambiguous '?<options>' on both sides throws", "[FactoryOptions]") {
    // If options are on the backend AND the fluid, that's almost
    // certainly a typo — fail loud with a clear message.
    REQUIRE_THROWS(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(R"(HEOS?{})", R"(Water?{})")));
}

TEST_CASE("AbstractState::factory: '?<options>' on the fluid side is hoisted", "[FactoryOptions]") {
    // PropsSI's extract_backend() puts the `?<options>` suffix on the
    // fluid side after the "::" split.  factory() must hoist it back
    // into the dispatch layer rather than dropping it silently.
    SECTION("empty options on fluid side accepted") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water?{}"));
        REQUIRE(AS != nullptr);
    }
    SECTION("non-empty options on fluid side reach the default overload (throws)") {
        REQUIRE_THROWS(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", R"(Water?{"x":1})")));
    }
}

TEST_CASE("AbstractState::factory: '?<options>' on the LAST mixture token is hoisted", "[FactoryOptions]") {
    // For mixtures, factory(string, string) splits on '&' before the
    // vector-overload runs — so e.g. "HEOS::R32&R125?{}" arrives as
    // fluid_names = ["R32", "R125?{}"] and the suffix lands on the
    // last token, not the first.  Make sure the parser scans the
    // whole vector and strips correctly.
    SECTION("empty options on last mixture token accepted") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32&R125?{}"));
        REQUIRE(AS != nullptr);
        REQUIRE(AS->fluid_names().size() == 2);
        REQUIRE(AS->fluid_names()[0] == "R32");
        REQUIRE(AS->fluid_names()[1] == "R125");  // trailing '?' stripped
    }
    SECTION("empty options on first mixture token accepted") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32?{}&R125"));
        REQUIRE(AS != nullptr);
        REQUIRE(AS->fluid_names()[0] == "R32");
        REQUIRE(AS->fluid_names()[1] == "R125");
    }
    SECTION("non-empty options on last mixture token reach the default overload (throws)") {
        REQUIRE_THROWS(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", R"(R32&R125?{"k":1})")));
    }
    SECTION("options on multiple mixture tokens is rejected") {
        REQUIRE_THROWS(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32?{}&R125?{}")));
    }
}

#endif  // ENABLE_CATCH
