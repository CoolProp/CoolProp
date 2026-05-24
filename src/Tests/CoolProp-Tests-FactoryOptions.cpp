// Catch2 tests for parse_factory_options() — the factory-string `?<options>`
// suffix parser.  Pure-parser tests; no backend wiring, no JSON validation
// (that lives one layer up).  See
// docs/superpowers/specs/2026-05-16-backend-options-string-design.md for
// the design.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <cstdio>
#    include <filesystem>
#    include <fstream>
#    include <string>

#    include "CoolProp/FactoryOptions.h"
#    include "Exceptions.h"
#    include "TestUtils.h"

namespace {

// RAII temp file for the @path tests.  Writes the given contents to
// a unique path under the system temp dir; removes on destruction.
class TempJSONFile
{
   public:
    explicit TempJSONFile(const std::string& contents) {
        path_ = std::filesystem::temp_directory_path()
                / ("cp_factopts_test_" + std::to_string(CoolProp::tests::test_pid()) + "_" + std::to_string(counter_++) + ".json");
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

TEST_CASE("parse_factory_options: no '?' suffix returns whole string + empty options", "[FactoryOptions]") {
    auto r = CoolProp::parse_factory_options("HEOS::Water");
    REQUIRE(r.clean_string == "HEOS::Water");
    REQUIRE(r.options_json.empty());
}

TEST_CASE("parse_factory_options: bare '?' yields empty options", "[FactoryOptions]") {
    auto r = CoolProp::parse_factory_options("HEOS::Water?");
    REQUIRE(r.clean_string == "HEOS::Water");
    REQUIRE(r.options_json.empty());
}

TEST_CASE("parse_factory_options: whitespace-only tail is treated as empty", "[FactoryOptions]") {
    auto r = CoolProp::parse_factory_options("HEOS::Water?   \t  ");
    REQUIRE(r.clean_string == "HEOS::Water");
    REQUIRE(r.options_json.empty());
}

TEST_CASE("parse_factory_options: inline JSON tail is returned verbatim", "[FactoryOptions]") {
    auto r = CoolProp::parse_factory_options(R"(SVDSBTL&HEOS::Water?{"critical_patch":"off"})");
    REQUIRE(r.clean_string == "SVDSBTL&HEOS::Water");
    REQUIRE(r.options_json == R"({"critical_patch":"off"})");
}

TEST_CASE("parse_factory_options: split is on the FIRST '?' only", "[FactoryOptions]") {
    SECTION("'?' inside a JSON string value") {
        auto r = CoolProp::parse_factory_options(R"(HEOS::Water?{"hint":"what?"})");
        REQUIRE(r.clean_string == "HEOS::Water");
        REQUIRE(r.options_json == R"({"hint":"what?"})");
    }
    SECTION("'?' is the entire string value") {
        auto r = CoolProp::parse_factory_options(R"(HEOS::Water?{"q":"?"})");
        REQUIRE(r.clean_string == "HEOS::Water");
        REQUIRE(r.options_json == R"({"q":"?"})");
    }
    SECTION("URL-style with both '?' and '&' inside the value") {
        auto r = CoolProp::parse_factory_options(R"(HEOS::Water?{"meta":{"url":"http://x.com/path?id=1&q=2"}})");
        REQUIRE(r.clean_string == "HEOS::Water");
        REQUIRE(r.options_json == R"({"meta":{"url":"http://x.com/path?id=1&q=2"}})");
    }
    SECTION("regex string with escaped '?'") {
        auto r = CoolProp::parse_factory_options(R"(HEOS::Water?{"regex":"^[A-Z]+\\?$"})");
        REQUIRE(r.clean_string == "HEOS::Water");
        REQUIRE(r.options_json == R"({"regex":"^[A-Z]+\\?$"})");
    }
    SECTION("multiple '?' chained inside one string value") {
        auto r = CoolProp::parse_factory_options(R"(HEOS::Water?{"chain":"first?second?third"})");
        REQUIRE(r.clean_string == "HEOS::Water");
        REQUIRE(r.options_json == R"({"chain":"first?second?third"})");
    }
}

TEST_CASE("parse_factory_options: '@path' reads file contents verbatim", "[FactoryOptions]") {
    SECTION("simple file path") {
        TempJSONFile f(R"({"critical_patch":"auto","grid":{"NT":200}})");
        const std::string s = "SVDSBTL&HEOS::Water?@" + f.path();
        auto r = CoolProp::parse_factory_options(s);
        REQUIRE(r.clean_string == "SVDSBTL&HEOS::Water");
        REQUIRE(r.options_json == R"({"critical_patch":"auto","grid":{"NT":200}})");
    }
    SECTION("file path with internal characters ('-', '_', '.', digits)") {
        TempJSONFile f(R"({"x":1})");
        const std::string s = "HEOS::Water?@" + f.path();
        auto r = CoolProp::parse_factory_options(s);
        REQUIRE(r.options_json == R"({"x":1})");
    }
    SECTION("missing '@path' file throws") {
        REQUIRE_THROWS(CoolProp::parse_factory_options("HEOS::Water?@/this/path/definitely/does/not/exist.json"));
    }
    SECTION("bare '@' with no path throws ValueError") {
        REQUIRE_THROWS_AS(CoolProp::parse_factory_options("HEOS::Water?@"), CoolProp::ValueError);
    }
}

TEST_CASE("parse_factory_options: clean_string preserves the full backend+fluid grammar", "[FactoryOptions]") {
    SECTION("simple backend") {
        auto r = CoolProp::parse_factory_options("HEOS::Water?{}");
        REQUIRE(r.clean_string == "HEOS::Water");
    }
    SECTION("source-backend split (& inside the cleaned half)") {
        auto r = CoolProp::parse_factory_options("SVDSBTL&HEOS::Water?{}");
        REQUIRE(r.clean_string == "SVDSBTL&HEOS::Water");
    }
    SECTION("no fluid name") {
        auto r = CoolProp::parse_factory_options("HEOS?{}");
        REQUIRE(r.clean_string == "HEOS");
        REQUIRE(r.options_json == "{}");
    }
}

// ---------------------------------------------------------------------------
// Factory entry-point integration: options are stripped from the backend
// string before dispatch, and the default generator overload rejects
// non-empty options with NotImplementedError for backends that haven't
// opted in.
// ---------------------------------------------------------------------------

#    include "AbstractState.h"

TEST_CASE("AbstractState::factory: '?<options>' is stripped before backend dispatch", "[FactoryOptions]") {
    SECTION("no options — existing behaviour unchanged") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
        REQUIRE(AS != nullptr);
        // HEOS factory returns the pure-fluid wrapper, which reports
        // "HelmholtzEOSBackend" rather than the mixture variant.
        REQUIRE(AS->backend_name() == "HelmholtzEOSBackend");
    }
    SECTION("empty options ('?{}') accepted by default overload") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS?{}", "Water"));
        REQUIRE(AS != nullptr);
    }
    SECTION("bare '?' accepted by default overload") {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS?", "Water"));
        REQUIRE(AS != nullptr);
    }
}

TEST_CASE("AbstractState::factory: non-empty options on an opted-out backend throw", "[FactoryOptions]") {
    REQUIRE_THROWS(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(R"(HEOS?{"some_key":1})", "Water")));
}

TEST_CASE("AbstractState::build_options_json: default returns empty string", "[FactoryOptions]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    REQUIRE(AS->build_options_json().empty());
}

#endif  // ENABLE_CATCH
