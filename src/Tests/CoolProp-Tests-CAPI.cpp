// Regression tests for the length checks in the exported C API (src/CoolPropLib.cpp).
//
// Each test puts a canary just past the region the caller declares, so a write
// (or a read that gets echoed back) beyond the declared length is detected even
// without a sanitizer.

#include <catch2/catch_all.hpp>

#include "CoolProp/CoolPropLib.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/numerics/numerics.h"

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace {

const double kCanary = -12345.678;

std::string errstring() {
    std::vector<char> buf(4096, '\0');
    get_global_param_string("errstring", buf.data(), static_cast<int>(buf.size()));
    return std::string(buf.data());
}

}  // namespace

TEST_CASE("Props1SImulti checks the output count against resdim1", "[CAPI]") {
    char backend[] = "HEOS";
    const double z[] = {1.0};
    std::vector<double> result(6, kCanary);

    SECTION("as many outputs as declared slots") {
        long resdim1 = 4;
        Props1SImulti("Tcrit,pcrit,molemass,Ttriple", backend, "Water", z, 1, result.data(), &resdim1);
        CHECK(resdim1 == 4);
        CHECK(result[0] == Catch::Approx(647.096));
        CHECK(result[4] == kCanary);
    }
    SECTION("one more output than declared slots is refused") {
        long resdim1 = 4;
        Props1SImulti("Tcrit,pcrit,molemass,Ttriple,Tmax", backend, "Water", z, 1, result.data(), &resdim1);
        CHECK(resdim1 == 0);
        CHECK(result[4] == kCanary);
        CHECK(result[5] == kCanary);
    }
    SECTION("negative resdim1 is refused") {
        long resdim1 = -1;
        Props1SImulti("Tcrit,pcrit", backend, "Water", z, 1, result.data(), &resdim1);
        CHECK(resdim1 == 0);
        CHECK(result[0] == kCanary);
    }
}

TEST_CASE("PropsSImulti refuses negative result dimensions", "[CAPI]") {
    char backend[] = "HEOS";
    const double z[] = {1.0};
    double T[] = {300.0};
    double p[] = {101325.0};
    std::vector<double> result(4, kCanary);
    long resdim1 = -1, resdim2 = -1;
    PropsSImulti("Dmolar,Hmolar", "T", T, 1, "P", p, 1, backend, "Water", z, 1, result.data(), &resdim1, &resdim2);
    CHECK(resdim1 == 0);
    CHECK(resdim2 == 0);
    CHECK(result[0] == kCanary);
}

TEST_CASE("HandleException refuses a negative message-buffer length", "[CAPI]") {
    // A long invalid fluid name makes the factory fail with a long error message
    const std::string fluid(200, 'A');
    std::vector<char> message(160, 'x');
    long errcode = 0;
    // Only the first 128 bytes are "owned" by the caller; the rest is canary
    long handle = AbstractState_factory("HEOS", fluid.c_str(), &errcode, message.data(), -1);
    CHECK(handle == -1);
    CHECK(errcode == 2);  // "buffer too small", nothing copied
    for (char c : message) {
        REQUIRE(c == 'x');
    }
}

TEST_CASE("AbstractState name getters respect the output-buffer length", "[CAPI]") {
    long errcode = 0;
    char message[1024] = {0};
    long handle = AbstractState_factory("HEOS", "Methane&Ethane&Propane", &errcode, message, sizeof(message));
    REQUIRE(errcode == 0);

    SECTION("fluid names: too small a buffer is refused, nothing written") {
        char names[16 + 8];
        std::memset(names, 'x', sizeof(names));
        AbstractState_fluid_names_checkedMemory(handle, names, 16, &errcode, message, sizeof(message));
        CHECK(errcode != 0);
        for (char c : names) {
            REQUIRE(c == 'x');
        }
    }
    SECTION("fluid names: adequate buffer") {
        char names[128];
        AbstractState_fluid_names_checkedMemory(handle, names, sizeof(names), &errcode, message, sizeof(message));
        CHECK(errcode == 0);
        CHECK(std::string(names) == "Methane,Ethane,n-Propane");
    }
    SECTION("fluid names: negative length is refused") {
        char names[8];
        std::memset(names, 'x', sizeof(names));
        AbstractState_fluid_names_checkedMemory(handle, names, -1, &errcode, message, sizeof(message));
        CHECK(errcode != 0);
        CHECK(names[0] == 'x');
    }
    SECTION("backend name: too small a buffer is refused") {
        char name[4 + 8];
        std::memset(name, 'x', sizeof(name));
        AbstractState_backend_name_checkedMemory(handle, name, 4, &errcode, message, sizeof(message));
        CHECK(errcode != 0);
        for (char c : name) {
            REQUIRE(c == 'x');
        }
    }
    AbstractState_free(handle, &errcode, message, sizeof(message));
}

TEST_CASE("AbstractState_set_fractions refuses a negative count", "[CAPI]") {
    long errcode = 0;
    char message[1024] = {0};
    long handle = AbstractState_factory("HEOS", "Methane&Ethane", &errcode, message, sizeof(message));
    REQUIRE(errcode == 0);
    const double z[] = {0.5, 0.5};
    AbstractState_set_fractions(handle, z, -1, &errcode, message, sizeof(message));
    CHECK(errcode != 0);
    AbstractState_free(handle, &errcode, message, sizeof(message));
}

TEST_CASE("get_parameter_information_string reads the info key only within n bytes", "[CAPI]") {
    // The declared n=16 bytes hold no terminator; the bytes after them are 'B'.
    // Before the fix the key was read with strlen and echoed into errstring.
    struct
    {
        char output[16];
        char after[32];
    } buf;
    std::memset(buf.output, 'A', sizeof(buf.output));
    std::memset(buf.after, 'B', sizeof(buf.after));
    buf.after[sizeof(buf.after) - 1] = '\0';
    CHECK(get_parameter_information_string("T", buf.output, sizeof(buf.output)) == 0);
    const std::string err = errstring();
    CAPTURE(err);
    CHECK(err.find(std::string(16, 'A')) != std::string::npos);
    CHECK(err.find("BB") == std::string::npos);

    SECTION("well-formed key still works") {
        char out[256] = "units";
        CHECK(get_parameter_information_string("T", out, sizeof(out)) == 1);
        CHECK(std::string(out) == "K");
    }
}

TEST_CASE("No exception escapes Props1 for fluid constants", "[CAPI]") {
    // convert_from_SI_to_kSI has no rule for Tcrit; this used to std::terminate
    const double v = Props1("Water", "Tcrit");
    CHECK(!ValidNumber(v));
    CHECK(!errstring().empty());
}

TEST_CASE("C_extract_backend refuses non-positive lengths", "[CAPI]") {
    char backend[8], fluid[8];
    std::memset(backend, 'x', sizeof(backend));
    std::memset(fluid, 'x', sizeof(fluid));
    CHECK(C_extract_backend("HEOS::Water", backend, -1, fluid, sizeof(fluid)) == -1);
    CHECK(backend[0] == 'x');
    CHECK(C_extract_backend("HEOS::Water", backend, sizeof(backend), fluid, 0) == -1);
    CHECK(fluid[0] == 'x');
}
