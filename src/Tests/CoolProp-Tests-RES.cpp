// Catch2 tests for the Residual Entropy Scaling (RES) transport PARAMETER TABLE.
//
// This file covers the loader only: dev/res_transport_parameters.json is compiled into the
// binary as res_transport_parameters_JSON.h and overlaid onto every CoolPropFluid as it is
// built.  Nothing consumes those parameters yet -- the RES model itself, and the tests that
// exercise it, arrive with the backend that reads them.
//
// What is worth asserting at this stage is that the table survives the trip from the published
// source tables into the fluid records intact, and that the two ways it can be wrong fail
// loudly rather than quietly:
//
//   * a fluid whose coefficients did not transfer to HEOS must carry NO HEOS block, so the
//     backend can refuse rather than return a plausible wrong number;
//   * a coefficient vector shorter than the model indexes must be rejected at load, because the
//     model reads fixed positions guarded by nothing but the `provided` flag.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/CoolPropFluid.h"
#    include "CoolProp/detail/json.h"
#    include "../Backends/Helmholtz/Fluids/FluidLibrary.h"
#    include "res_transport_parameters_JSON.h"

#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

CoolPropFluid fluid(const std::string& name) {
    return get_library().get(name);  // returns by value
}

// Dilute-gas viscosity coefficients for R134a, read straight off
// dev/RES_reference/martinek2025_viscosity/Dilute_gas_viscosity.txt so the assertion is against
// the published source table rather than against whatever the converter last happened to emit.
//
// Listed here in ASCENDING powers, as CoolProp stores and evaluates them:
//   eta0/uPa.s = n[0] + n[1]*T + n[2]*T^2 + n[3]*T^3 + n[4]*T^4
// The source table lists the same five numbers in the same order but heads the T^4 column n0,
// i.e. descending, so the converter reverses them once in dilute_ascending().  This test is what
// pins that reversal down from the other end.
const double R134A_VIS_DILUTE[5] = {-1.766998e-01, 4.267365e-02, -9.057271e-06, 4.498500e-09, -1.269154e-12};

// eta0 at 300 K from those coefficients, as an independent check on the ORDER: get the
// convention backwards and this is off by nine orders of magnitude, not by a rounding error.
const double R134A_ETA0_300K_uPas = 11.9214201626;

}  // namespace

TEST_CASE("RES parameters reach the fluid records", "[RES][transport]") {
    CoolPropFluid f = fluid("R134a");

    SECTION("viscosity") {
        const ViscosityRESData& d = f.transport.viscosity_res;
        REQUIRE(d.provided);
        REQUIRE(d.n_dilute.size() == RES_N_DILUTE);
        REQUIRE(d.n_res.size() == RES_N_RES_VISCOSITY);
        // The dilute polynomial is EOS-independent; it is the one block that must be identical
        // for every backend, so it is the sharpest check that the right row was read.
        for (std::size_t i = 0; i < RES_N_DILUTE; ++i) {
            CAPTURE(i);
            CHECK(d.n_dilute[i] == Catch::Approx(R134A_VIS_DILUTE[i]).epsilon(1e-12));
        }
        double T = 300.0;
        double eta0 = d.n_dilute[0] + T * (d.n_dilute[1] + T * (d.n_dilute[2] + T * (d.n_dilute[3] + T * d.n_dilute[4])));
        CHECK(eta0 == Catch::Approx(R134A_ETA0_300K_uPas).epsilon(1e-9));
        CHECK(d.xita > 0);
        CHECK(d.molar_mass == Catch::Approx(f.molar_mass()).epsilon(1e-12));
    }

    SECTION("conductivity") {
        const ConductivityRESData& d = f.transport.conductivity_res;
        REQUIRE(d.provided);
        REQUIRE(d.n_dilute.size() == RES_N_DILUTE);
        REQUIRE(d.n_res.size() == RES_N_RES_CONDUCTIVITY);
        CHECK(d.molar_mass == Catch::Approx(f.molar_mass()).epsilon(1e-12));
    }

    SECTION("critical enhancement is present and self-consistent") {
        const ConductivityRESData& d = f.transport.conductivity_res;
        REQUIRE(d.crit_provided);
        // crit_provided is exactly the conjunction below; an all-zero record in the source table
        // means "not fitted", and letting it through divides by t_ref and Gamma downstream.
        CHECK(d.t_ref > 0);
        CHECK(d.Gamma > 0);
        CHECK(d.phi0 > 0);
        CHECK(d.q_D > 0);
        CHECK(d.gamma_uni > 0);
    }
}

TEST_CASE("RES parameters are absent for fluids that were never fitted", "[RES][transport]") {
    // 15 of CoolProp's fluids appear in neither paper's table.  Air is one, and is a good pick
    // for the opposite reason to R134a: it is a pseudo-pure mixture, exactly the sort of entry a
    // loose alias match could wrongly bind to one of its components.  The flags must say "no
    // parameters" rather than leaving default-constructed vectors that read as usable.
    for (const char* name : {"Air", "R410A", "PropyleneGlycol"}) {
        CAPTURE(name);
        CoolPropFluid f = fluid(name);
        CHECK_FALSE(f.transport.viscosity_res.provided);
        CHECK_FALSE(f.transport.conductivity_res.provided);
    }
}

TEST_CASE("RES parameters are withheld where they do not transfer to HEOS", "[RES][transport]") {
    // These 14 fluids deliberately carry no HEOS residual coefficients: the REFPROP-fitted ones
    // do not survive the change of equation of state (>1% transport deviation at
    // residual-dominated states, or >5% in s_res itself).  Withholding them is what lets the
    // backend refuse rather than return a plausible wrong number, so the absence has to be
    // asserted -- a converter change that quietly started shipping them would otherwise be
    // invisible.  See HEOS_TRANSFER_EXCLUDE in dev/convert_RES_csv_to_json.py for the derivation.
    //
    // Ten of them still carry PR and SRK blocks.  The other four -- R1123, R1224YDZ, R1243ZF and
    // VINYLCHLORIDE -- carry no residual coefficients for ANY equation of state, because they are
    // HEOS-only fluids and the cubic gate withholds their PR/SRK blocks as unreachable.  Their
    // entries are the dilute polynomial and nothing else, so they can never produce a RES value.
    // Spelled exactly as HEOS_TRANSFER_EXCLUDE spells them, which is also how they key the JSON
    // table: every one resolves through JSONFluidLibrary's alias lookup.
    const std::vector<std::string> withheld = {"BENZENE",  "D5",      "HEPTANE", "MD3M", "MD4M", "R1123", "R1224YDZ",
                                               "R1233ZDE", "R1234YF", "R1243ZF", "R13",  "R161", "R41",   "VINYLCHLORIDE"};
    REQUIRE(withheld.size() == 14);
    for (const std::string& name : withheld) {
        CAPTURE(name);
        // FAIL rather than skip on a name CoolProp cannot resolve: these fluids are withheld from
        // the RES table, not from CoolProp, so a lookup failure is a broken test, not a skip.
        CoolPropFluid f;
        REQUIRE_NOTHROW(f = fluid(name));
        CHECK_FALSE(f.transport.viscosity_res.provided);
        CHECK_FALSE(f.transport.conductivity_res.provided);
        // The dilute polynomial is EOS-independent and IS still shipped; only the residual
        // coefficients are withheld.  Asserting this keeps the two apart.
        CHECK(f.transport.viscosity_res.n_dilute.size() == RES_N_DILUTE);
    }
}

TEST_CASE("RES loader rejects coefficient vectors the model would index past", "[RES][transport]") {
    // The shipped table is generated with fixed counts and cannot trip this, which is exactly
    // why it needs a test: a regenerated or hand-edited res_transport_parameters.json is the
    // reachable path, and without the check it is an out-of-bounds read rather than an error.
    const auto make = [](const char* n_dilute, const char* n_res) {
        return cpjson::parse(std::string(R"({"viscosity": {"R134A": {"dilute": {"n": )") + n_dilute + R"(}, "HEOS": {"n": )" + n_res
                             + R"(, "xita": 1.0, "group": 1}}},)" + R"( "conductivity": {}})");
    };

    CoolPropFluid f;
    f.name = "R134A";

    SECTION("a well-formed record loads") {
        REQUIRE_NOTHROW(JSONFluidLibrary::load_RES_transport_parameters(make("[1,2,3,4,5]", "[1,2,3]"), "HEOS", f, 0.1));
        CHECK(f.transport.viscosity_res.provided);
    }
    SECTION("too few dilute coefficients") {
        CHECK_THROWS_AS(JSONFluidLibrary::load_RES_transport_parameters(make("[1,2,3,4]", "[1,2,3]"), "HEOS", f, 0.1), ValueError);
    }
    SECTION("too few residual coefficients") {
        CHECK_THROWS_AS(JSONFluidLibrary::load_RES_transport_parameters(make("[1,2,3,4,5]", "[1,2]"), "HEOS", f, 0.1), ValueError);
    }
    SECTION("a short dilute vector is rejected even with no EOS block present") {
        // The withheld fluids have exactly this shape -- a dilute polynomial and no HEOS block --
        // and the dilute term is read off fluids whose `provided` is false, so gating the length
        // check on the residual block would leave precisely those vectors unchecked.
        const nlohmann::json short_dilute = cpjson::parse(std::string(R"({"viscosity": {"R134A": {"dilute": {"n": [1,2,3]}}}, "conductivity": {}})"));
        CHECK_THROWS_AS(JSONFluidLibrary::load_RES_transport_parameters(short_dilute, "HEOS", f, 0.1), ValueError);
    }
}

TEST_CASE("RES table is validated as a whole before any fluid is built", "[RES][transport]") {
    // A malformed record must not be discovered mid-way through add_many().  A throw there is
    // rethrown by add_one, aborts add_many's loop, and is caught by load() with nothing but a
    // stdout line -- leaving a PARTIAL fluid library, permanently, because _is_empty is already
    // false and std::call_once has committed.  Validating the whole table up front turns that
    // into a loud failure before the first fluid is added.
    SECTION("the shipped table passes") {
        REQUIRE_NOTHROW(JSONFluidLibrary::validate_RES_transport_table(cpjson::parse(std::string(res_transport_parameters_JSON))));
    }
    SECTION("a short residual vector anywhere in the table is caught") {
        const nlohmann::json bad = cpjson::parse(std::string(R"({"viscosity": {"OK": {"dilute": {"n": [1,2,3,4,5]}, "HEOS": {"n": [1,2,3]}},)"
                                                             R"( "BAD": {"dilute": {"n": [1,2,3,4,5]}, "SRK": {"n": [1,2]}}}, "conductivity": {}})"));
        CHECK_THROWS_AS(JSONFluidLibrary::validate_RES_transport_table(bad), ValueError);
    }
    SECTION("a short dilute vector anywhere in the table is caught") {
        const nlohmann::json bad = cpjson::parse(std::string(R"({"viscosity": {"BAD": {"dilute": {"n": [1,2]}}}, "conductivity": {}})"));
        CHECK_THROWS_AS(JSONFluidLibrary::validate_RES_transport_table(bad), ValueError);
    }
    SECTION("an empty table is not an error") {
        REQUIRE_NOTHROW(JSONFluidLibrary::validate_RES_transport_table(cpjson::parse(std::string(R"({})"))));
    }
}

#endif  // ENABLE_CATCH
