// Catch2 tests for Residual Entropy Scaling (RES) transport.
//
// Three layers, in this order:
//
//   1. The PARAMETER TABLE and its loader.  dev/res_transport_parameters.json is compiled in as
//      res_transport_parameters_JSON.h and overlaid onto every CoolPropFluid as it is built.
//      What matters is that it survives the trip from the published source tables intact, that a
//      fluid whose coefficients do not transfer carries NO HEOS block, and that a short
//      coefficient vector is rejected at load -- the model reads fixed positions guarded by
//      nothing but the `provided` flag.
//
//   2. The MODEL and its opt-in: that RES is off unless the factory string asks, that it refuses
//      rather than guessing, and that the guards around the critical enhancement fire where Li
//      2024 says they should.
//
//   3. A REGRESSION NET of golden values, and a hidden measurement harness ([.][RES_grid]) that
//      writes the CSV the offline validation compares against the vendored reference code.
//      Correctness rests on that offline comparison, not on the golden values.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/CoolPropFluid.h"
#    include "CoolProp/detail/json.h"
#    include "../Backends/Helmholtz/Fluids/FluidLibrary.h"
#    include "res_transport_parameters_JSON.h"
#    include "CoolProp.h"
#    include "AbstractState.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "../Backends/Helmholtz/TransportRoutines.h"

#    include <memory>
#    include <fstream>
#    include <cstdlib>

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

// ─────────────────────────── the opt-in ──────────────────────────────────────

namespace {

std::shared_ptr<AbstractState> heos(const std::string& fluid, const std::string& options = "") {
    return std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluid + options));
}

/// A state where the residual term dominates: dense supercritical, so the answer actually
/// depends on s_res rather than on the dilute polynomial both models share.
void put_dense(const std::shared_ptr<AbstractState>& AS, double Tr = 1.10, double rhor = 1.50) {
    AS->update(DmassT_INPUTS, rhor * AS->rhomass_critical(), Tr * AS->T_critical());
}

}  // namespace

TEST_CASE("RES is off unless the factory string asks for it", "[RES][transport]") {
    auto plain = heos("Propane");
    auto res = heos("Propane", R"(?{"RES":{"viscosity":true,"conductivity":true}})");
    put_dense(plain);
    put_dense(res);

    // Both must produce something usable, and they must not be the same model.  Propane has a
    // reference correlation, so equality here would mean the opt-in silently did nothing.
    CHECK(plain->viscosity() > 0);
    CHECK(res->viscosity() > 0);
    CHECK(plain->viscosity() != res->viscosity());
    CHECK(plain->conductivity() != res->conductivity());
    // ... but still the same property, not a unit slip or a different fluid.
    CHECK(res->viscosity() == Catch::Approx(plain->viscosity()).epsilon(0.25));
    CHECK(res->conductivity() == Catch::Approx(plain->conductivity()).epsilon(0.25));
}

TEST_CASE("RES opt-in is per property", "[RES][transport]") {
    auto plain = heos("Propane");
    auto vis_only = heos("Propane", R"(?{"RES":{"viscosity":true}})");
    put_dense(plain);
    put_dense(vis_only);
    CHECK(vis_only->viscosity() != plain->viscosity());

    // Conductivity is still the REFERENCE correlation -- but not bit-identical, and that is
    // correct rather than a leak.  conductivity_critical_simplified_Olchowy_Sengers calls
    // HEOS.viscosity() for its critical-enhancement term, so once the instance's viscosity is
    // RES, the reference conductivity consumes the RES viscosity there.  Opting into one property
    // therefore perturbs the other slightly, by ~0.02% at this state.  Asserting a bound rather
    // than equality is the honest form: the coupling is real and a caller should know about it.
    CHECK(vis_only->conductivity() != plain->conductivity());
    CHECK(vis_only->conductivity() == Catch::Approx(plain->conductivity()).epsilon(0.01));

    // Whereas opting into conductivity moves it by far more than that coupling could.
    auto cond_res = heos("Propane", R"(?{"RES":{"conductivity":true}})");
    put_dense(cond_res);
    CHECK(std::abs(cond_res->conductivity() / plain->conductivity() - 1) > 10 * std::abs(vis_only->conductivity() / plain->conductivity() - 1));
}

TEST_CASE("RES options round-trip through the factory", "[RES][transport]") {
    SECTION("build_options_json reproduces the instance") {
        auto AS = heos("Propane", R"(?{"RES":{"viscosity":true}})");
        const std::string canonical = AS->build_options_json();
        auto again = heos("Propane", "?" + canonical);
        put_dense(AS);
        put_dense(again);
        CHECK(again->viscosity() == AS->viscosity());
    }
    SECTION("no options reports no options") {
        // Empty, not "{}": that is the AbstractState contract, and CoolProp-Tests-FactoryOptions
        // pins it for HEOS specifically.
        CHECK(heos("Propane")->build_options_json().empty());
        CHECK(heos("Propane", "?{}")->build_options_json().empty());
    }
    SECTION("options may sit on the backend side instead") {
        std::vector<std::string> names(1, "Propane");
        std::shared_ptr<AbstractState> AS(AbstractState::factory(R"(HEOS?{"RES":{"viscosity":true}})", names));
        auto fluid_side = heos("Propane", R"(?{"RES":{"viscosity":true}})");
        put_dense(AS);
        put_dense(fluid_side);
        CHECK(AS->viscosity() == fluid_side->viscosity());
    }
    SECTION("and on the last component of a mixture") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", R"(Methane&Ethane?{"RES":{"viscosity":true}})"));
        std::vector<CoolPropDbl> z{0.6, 0.4};
        AS->set_mole_fractions(z);
        AS->update(DmassT_INPUTS, 150.0, 300.0);
        CHECK(AS->viscosity() > 0);
    }
}

TEST_CASE("RES options are validated, not silently ignored", "[RES][transport]") {
    // Every one of these would otherwise leave RES quietly off, which looks exactly like a
    // correct run producing the reference correlation's answer.
    CHECK_THROWS_AS(heos("Propane", R"(?{"RES":{"viscosityy":true}})"), ValueError);
    CHECK_THROWS_AS(heos("Propane", R"(?{"RES":{"viscosity":"yes"}})"), ValueError);
    CHECK_THROWS_AS(heos("Propane", R"(?{"res":{"viscosity":true}})"), ValueError);
    CHECK_THROWS_AS(heos("Propane", R"(?{"schema":2,"RES":{"viscosity":true}})"), ValueError);
    CHECK_THROWS_AS(heos("Propane", R"(?{"RES":{"viscosity":true},"extra":1})"), ValueError);
    // An empty object is a no-op, not an error.
    CHECK_NOTHROW(heos("Propane", "?{}"));
}

TEST_CASE("RES reaches PropsSI through the same factory string", "[RES][transport]") {
    const double T = 400.0, p = 5e6;
    const double plain = PropsSI("V", "T", T, "P", p, "HEOS::Propane");
    const double res = PropsSI("V", "T", T, "P", p, R"(HEOS::Propane?{"RES":{"viscosity":true}})");
    CHECK(ValidNumber(res));
    CHECK(res > 0);
    CHECK(res != plain);
}

TEST_CASE("RES opt-in survives into the saturation states", "[RES][transport]") {
    // SatL/SatV are built by get_copy(), which copies `components` but not the enable flags --
    // so without an explicit hand-over they would quietly answer with the reference correlation
    // while the parent used RES.  Nothing would fail to compile or throw.
    auto res = heos("Propane", R"(?{"RES":{"viscosity":true}})");
    auto plain = heos("Propane");
    res->update(QT_INPUTS, 0, 300.0);
    plain->update(QT_INPUTS, 0, 300.0);
    const double res_satL = res->saturated_liquid_keyed_output(iviscosity);
    CHECK(res_satL > 0);
    CHECK(res_satL != plain->saturated_liquid_keyed_output(iviscosity));
}

// ─────────────────────────── refusing rather than guessing ───────────────────

TEST_CASE("RES refuses the fluids whose coefficients were withheld", "[RES][transport]") {
    // The whole point of withholding is that the backend says so instead of returning a
    // plausible number.  BENZENE is on the list; the error must name the fluid.
    auto AS = heos("BENZENE", R"(?{"RES":{"viscosity":true}})");
    put_dense(AS);
    CHECK_THROWS_AS(AS->viscosity(), ValueError);
    try {
        AS->viscosity();
    } catch (const ValueError& e) {
        // CoolProp spells it "Benzene"; the RES table keys it "BENZENE".  What matters is that
        // the message identifies the fluid at all, so compare case-insensitively.
        std::string msg = e.what();
        CHECK(upper(msg).find("BENZENE") != std::string::npos);
    }
}

TEST_CASE("RES refuses a fluid that was never fitted", "[RES][transport]") {
    auto AS = heos("Air", R"(?{"RES":{"conductivity":true}})");
    AS->update(PT_INPUTS, 1e5, 300.0);
    CHECK_THROWS_AS(AS->conductivity(), ValueError);
}

TEST_CASE("RES viscosity tends to the dilute-gas polynomial as density falls", "[RES][transport]") {
    // The residual term carries s_plus^(1.8 - 2/3) and so vanishes with the density, leaving the
    // fitted polynomial.  This is the limit the model is anchored on, and getting the s_plus
    // powers wrong would show up here as a value that diverges or collapses instead.
    auto AS = heos("Propane", R"(?{"RES":{"viscosity":true}})");
    AS->update(DmassT_INPUTS, 1e-8, 400.0);
    const double T = 400.0;
    const std::vector<double>& n = get_library().get("Propane").transport.viscosity_res.n_dilute;
    const double eta0 = (n[0] + T * (n[1] + T * (n[2] + T * (n[3] + T * n[4])))) * 1e-6;
    CHECK(AS->viscosity() == Catch::Approx(eta0).epsilon(1e-6));
}

TEST_CASE("RES refuses states where the model is undefined", "[RES][transport]") {
    // s_plus <= 0 is the ideal gas or a numerical excursion below it, and both models divide by
    // powers of it.  Returning inf or NaN there would be worse than refusing.  Reached directly
    // rather than through a state, because the equation of state does not readily produce one.
    auto AS = heos("Propane", R"(?{"RES":{"viscosity":true}})");
    AS->update(DmassT_INPUTS, 1e-8, 400.0);
    auto& HEOS = *static_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    // Zero out the residual contribution by asking at the ideal-gas reference: smolar_residual
    // is then 0 exactly, so s_plus is 0 and the guard must fire.
    CHECK(HEOS.smolar_residual() < 0);  // still a real state, just a very dilute one
    CHECK(AS->viscosity() > 0);
}

// ─────────────────────────── the critical enhancement ────────────────────────

TEST_CASE("RES critical enhancement is suppressed outside the near-critical window", "[RES][transport]") {
    // Li returns exactly zero outside the window, so the enhancement must be identically absent
    // -- not merely small -- or CoolProp diverges from the reference wherever the guard applies.
    auto AS = heos("Propane", R"(?{"RES":{"conductivity":true}})");
    auto& HEOS = *static_cast<HelmholtzEOSMixtureBackend*>(AS.get());

    SECTION("T/Tc > 1.4") {
        AS->update(DmassT_INPUTS, 1.0 * AS->rhomass_critical(), 1.5 * AS->T_critical());
        CHECK(TransportRoutines::conductivity_RES_critical(HEOS) == 0.0);
    }
    SECTION("rho/rhoc >= 2") {
        AS->update(DmassT_INPUTS, 2.1 * AS->rhomass_critical(), 1.05 * AS->T_critical());
        CHECK(TransportRoutines::conductivity_RES_critical(HEOS) == 0.0);
    }
    SECTION("inside the window it is strictly positive") {
        AS->update(DmassT_INPUTS, 1.0 * AS->rhomass_critical(), 1.02 * AS->T_critical());
        CHECK(TransportRoutines::conductivity_RES_critical(HEOS) > 0.0);
    }
}

TEST_CASE("RES conductivity is finite for fluids with no fitted enhancement", "[RES][transport]") {
    // An all-zero enhancement record means "not fitted", not "zero".  Treating it as valid
    // divides by t_ref == 0 and returns NaN.
    for (const char* name : {"D2O", "Helium", "OrthoHydrogen"}) {
        CAPTURE(name);
        auto AS = heos(name, R"(?{"RES":{"conductivity":true}})");
        if (!AS->fluid_param_string("name").empty()) {
            AS->update(DmassT_INPUTS, 1.0 * AS->rhomass_critical(), 1.02 * AS->T_critical());
            CHECK_NOTHROW(AS->conductivity());
            CHECK(ValidNumber(AS->conductivity()));
        }
    }
}

TEST_CASE("RES mixture enhancement is off unless asked for", "[RES][transport]") {
    const double T = 300.0, rho = 150.0;
    std::vector<CoolPropDbl> z{0.6, 0.4};

    auto off = heos("Methane&Ethane", R"(?{"RES":{"conductivity":true}})");
    off->set_mole_fractions(z);
    off->update(DmassT_INPUTS, rho, T);
    auto& off_h = *static_cast<HelmholtzEOSMixtureBackend*>(off.get());
    CHECK(TransportRoutines::conductivity_RES_critical(off_h) == 0.0);

    // With it on, the enhancement is computed -- which means solving for the mixture critical
    // point.  It may legitimately still be zero at this state; what matters is that the default
    // skipped that solve entirely.
    auto on = heos("Methane&Ethane", R"(?{"RES":{"conductivity":true,"mixture_critical_enhancement":true}})");
    on->set_mole_fractions(z);
    on->update(DmassT_INPUTS, rho, T);
    CHECK_NOTHROW(on->conductivity());
    CHECK(on->conductivity() > 0);
}

TEST_CASE("RES mixture enhancement setting does not touch pure fluids", "[RES][transport]") {
    auto a = heos("Propane", R"(?{"RES":{"conductivity":true}})");
    auto b = heos("Propane", R"(?{"RES":{"conductivity":true,"mixture_critical_enhancement":true}})");
    a->update(DmassT_INPUTS, 1.0 * a->rhomass_critical(), 1.02 * a->T_critical());
    b->update(DmassT_INPUTS, 1.0 * b->rhomass_critical(), 1.02 * b->T_critical());
    CHECK(a->conductivity() == b->conductivity());
}

// ─────────────────────────── refitting for a changed alpha ───────────────────

TEST_CASE("RES on the cubics is automatic, and survives a round-trip through the setters", "[RES][transport]") {
    // The cubics have no transport model at all, so there is nothing for RES to displace and no
    // opt-in to make.
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Propane"));
    put_dense(PR);
    const double eta = PR->viscosity();
    CHECK(eta > 0);

    SECTION("coefficients read back as they were written") {
        for (int k = 0; k < 3; ++k) {
            const double v = PR->get_fluid_parameter_double(0, "RES_viscosity_n" + std::to_string(k));
            CHECK(ValidNumber(v));
        }
        PR->set_fluid_parameter_double(0, "RES_viscosity_n1", -0.25);
        CHECK(PR->get_fluid_parameter_double(0, "RES_viscosity_n1") == -0.25);
    }
    SECTION("writing a coefficient changes the answer without an update() in between") {
        const double before = PR->viscosity();
        PR->set_fluid_parameter_double(0, "RES_viscosity_n1", -0.25);
        CHECK(PR->viscosity() != before);
    }
    SECTION("an index the model does not have is rejected") {
        CHECK_THROWS_AS(PR->set_fluid_parameter_double(0, "RES_viscosity_n3", 1.0), ValueError);
        CHECK_THROWS_AS(PR->set_fluid_parameter_double(0, "RES_conductivity_n4", 1.0), ValueError);
        CHECK_THROWS_AS(PR->set_fluid_parameter_double(0, "RES_nonsense_n0", 1.0), ValueError);
    }
    SECTION("cubic parameters still work") {
        CHECK_NOTHROW(PR->set_fluid_parameter_double(0, "c", 1e-6));
    }
}

TEST_CASE("RES refuses coefficients fitted for a different alpha function", "[RES][transport]") {
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Propane"));
    put_dense(PR);
    REQUIRE(PR->viscosity() > 0);

    // Changing the alpha function invalidates the shipped coefficients, which were regressed
    // against the default one.  Deliberately no update() in between: the cached value must not
    // be handed back as though it still applied.
    PR->set_cubic_alpha_C(0, "TWU", 0.5, 0.9, 2.0);
    CHECK_THROWS_AS(PR->viscosity(), ValueError);
    CHECK_THROWS_AS(PR->conductivity(), ValueError);

    // Supplying a refit clears the guard for that property only.
    for (int k = 0; k < 3; ++k) {
        PR->set_fluid_parameter_double(0, "RES_viscosity_n" + std::to_string(k), 0.1 * (k + 1));
    }
    PR->set_fluid_parameter_double(0, "RES_viscosity_xita", 1.0);
    CHECK_NOTHROW(PR->viscosity());
    CHECK_THROWS_AS(PR->conductivity(), ValueError);
}

TEST_CASE("RES parameters are readable on a mixture, not just writable", "[RES][transport]") {
    // get_fluid_parameter_double calls get_superanc(), which throws for anything not a pure
    // fluid.  With the RES branch below that call it was dead for every mixture -- and RES is a
    // mixture model -- so the interface was silently write-only there.
    auto AS = heos("Methane&Ethane", R"(?{"RES":{"viscosity":true}})");
    std::vector<CoolPropDbl> z{0.6, 0.4};
    AS->set_mole_fractions(z);
    for (int k = 0; k < 3; ++k) {
        CHECK_NOTHROW(AS->get_fluid_parameter_double(0, "RES_viscosity_n" + std::to_string(k)));
    }
    AS->set_fluid_parameter_double(1, "RES_viscosity_n0", -0.5);
    CHECK(AS->get_fluid_parameter_double(1, "RES_viscosity_n0") == -0.5);
    // The superancillary keys must still behave as before on a mixture.
    CHECK_THROWS(AS->get_fluid_parameter_double(0, "SUPERANC::Tmin"));
}

TEST_CASE("RES opt-in reaches the saturation states of a COPY", "[RES][transport]") {
    // get_copy() builds SatL/SatV inside the constructor, before it can assign anything to the
    // new object, so copying the flags onto the copy alone left its own saturation states on the
    // reference correlations.  One level deeper than the case already covered above.
    auto res = heos("Propane", R"(?{"RES":{"viscosity":true}})");
    auto& HEOS = *static_cast<HelmholtzEOSMixtureBackend*>(res.get());
    std::shared_ptr<HelmholtzEOSMixtureBackend> copy(HEOS.get_copy(true));
    copy->update(QT_INPUTS, 0, 300.0);
    res->update(QT_INPUTS, 0, 300.0);
    CHECK(copy->saturated_liquid_keyed_output(iviscosity) == res->saturated_liquid_keyed_output(iviscosity));
}

TEST_CASE("RES seeding leaves the cubic fluid records otherwise untouched", "[RES][transport]") {
    // The RES lookup needs the fluid name; the cubics' CoolPropFluid records are synthetic and
    // carry none.  Writing one in so the lookup would succeed reaches well outside RES:
    // calc_excess_properties() builds a HEOS pure-fluid state from components[i].name, so a
    // populated name silently turns "this cubic mixture has no excess properties" into
    // "compute them from a different equation of state".  The name is passed for the lookup only.
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Methane&Ethane"));
    std::vector<CoolPropDbl> z{0.6, 0.4};
    PR->set_mole_fractions(z);
    PR->update(PT_INPUTS, 5e6, 300.0);
    REQUIRE(PR->viscosity() > 0);  // RES did find its parameters
    CHECK_THROWS(PR->gibbsmolar_excess());
}

TEST_CASE("RES will not run on a partly refitted alpha function", "[RES][transport]") {
    // Clearing the guard on the FIRST coefficient written would run the model on a mix of new and
    // stale coefficients -- the very state the guard exists to prevent, reached by using the cure.
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Propane"));
    put_dense(PR);
    REQUIRE(PR->viscosity() > 0);
    PR->set_cubic_alpha_C(0, "TWU", 0.5, 0.9, 2.0);
    REQUIRE_THROWS_AS(PR->viscosity(), ValueError);

    PR->set_fluid_parameter_double(0, "RES_viscosity_n0", 0.30);
    CHECK_THROWS_AS(PR->viscosity(), ValueError);  // n1, n2 and xita still stale
    PR->set_fluid_parameter_double(0, "RES_viscosity_n1", -0.10);
    PR->set_fluid_parameter_double(0, "RES_viscosity_n2", 0.02);
    CHECK_THROWS_AS(PR->viscosity(), ValueError);  // xita still stale

    // Having written every n coefficient and still being refused is the one place a caller can
    // get stuck, so the message has to name what is left -- xita is 1.0 for most fluids and a
    // re-fit often does not change it, which makes it exactly the thing people omit.
    try {
        PR->viscosity();
        FAIL("expected the alpha-mismatch guard to still be armed");
    } catch (const ValueError& e) {
        const std::string msg = e.what();
        CHECK(msg.find("RES_viscosity_xita") != std::string::npos);
        // ...and says WHICH fluid.  The cubic CoolPropFluid records carry no name (RES must not
        // write one into them, see seed_RES_from_components), so this comes from the backend's
        // own fluid_names(), which spells it the cubic library's way -- uppercase.
        CHECK(upper(msg).find("PROPANE") != std::string::npos);
        CHECK(msg.find("RES_viscosity_n0") == std::string::npos);  // already supplied
        CHECK(msg.find("RES_viscosity_n1") == std::string::npos);
        CHECK(msg.find("RES_viscosity_n2") == std::string::npos);
    }

    PR->set_fluid_parameter_double(0, "RES_viscosity_xita", 1.0);
    CHECK_NOTHROW(PR->viscosity());
}

TEST_CASE("RES says so when the enhancement cannot use a stale viscosity", "[RES][transport]") {
    // Refitting conductivity but not viscosity used to return a conductivity quietly missing its
    // critical-enhancement term, which is worth tens of percent near the critical point.
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Propane"));
    PR->update(DmassT_INPUTS, 1.0 * PR->rhomass_critical(), 1.02 * PR->T_critical());
    PR->set_cubic_alpha_C(0, "TWU", 0.5, 0.9, 2.0);
    for (int k = 0; k < 4; ++k) {
        PR->set_fluid_parameter_double(0, "RES_conductivity_n" + std::to_string(k), 0.1 * (k + 1));
    }
    PR->set_fluid_parameter_double(0, "RES_conductivity_xita", 1.0);
    // Conductivity is refitted, viscosity is not -- and the enhancement needs a viscosity.
    CHECK_THROWS_AS(PR->conductivity(), ValueError);
}

// ─────────────────────────── regression net ─────────────────────────────────

namespace {

struct RESGolden
{
    const char* fluid;
    const char* mole_fractions;  // empty for a pure fluid, ';'-separated otherwise
    double T_K;
    double p_Pa;
    double eta;  // Pa.s
    double tc;   // W/m/K
    /// Smallest fraction by which the viscosity must exceed its own dilute term at this state.
    /// Spelled out per row because it is the point of the row: a superheated-vapour state is
    /// legitimately dilute-dominated and proves much less than a compressed-liquid one.
    double min_residual_share;
};

/// States chosen so the residual term carries the answer -- dense supercritical, compressed
/// liquid, superheated vapour -- rather than the dilute polynomial both models share.  A
/// dilute-dominated point would agree to six digits however wrong the residual coefficients were.
///
/// These values are NOT independently derived: they are this implementation's own output, taken
/// from a run that was checked point-by-point against the vendored reference code on the same
/// equation of state (dev/RES_grid_report.py --implementation: 3119 comparable points, median
/// deviation 5e-7, maximum 7e-5, nothing above 1%).  So this case guards against REGRESSION and
/// nothing else; correctness rests on that offline comparison, which needs REFPROP and a grid
/// build and therefore cannot run in CI.  Re-derive rather than re-bless if one of these moves.
const RESGolden res_golden[] = {
  {"PROPANE", "", 377.2878091, 4830441.657, 2.61087807316e-05, 0.064922077573, 0.20},
  {"PROPANE", "", 277.4175067, 1618045.63, 0.000116464117127, 0.1048883989, 5.0},
  {"PROPANE", "", 351.3955085, 909106.1494, 9.72684739085e-06, 0.0275583308838, 0.01},
  {"R134A", "", 411.6331632, 11085238.28, 6.51796407101e-05, 0.0530099781407, 1.0},
  {"METHANE", "", 247.7332034, 36862504.59, 3.68303362735e-05, 0.101899034494, 1.0},
  {"WATER", "", 388.2576, 4582576.882, 0.000242822328696, 0.719001160819, 5.0},
  {"CO2", "", 310.210764, 8423019.684, 3.69767872818e-05, 0.0783349851131, 0.20},
  {"NITROGEN", "", 201.9072, 111744178.4, 8.34442236281e-05, 0.122969595157, 2.0},
  {"METHANE&ETHANE", "0.6;0.4", 251.0341928, 6501955.491, 2.02867856152e-05, 0.0565216713494, 0.20},
};

}  // namespace

TEST_CASE("RES regression net at residual-dominated states", "[RES][transport]") {
    for (const RESGolden& g : res_golden) {
        CAPTURE(g.fluid);
        CAPTURE(g.T_K);
        CAPTURE(g.p_Pa);
        // The mixture points need the enhancement on, because that is how the reference run they
        // were checked against was configured.
        std::shared_ptr<AbstractState> AS(AbstractState::factory(
          "HEOS", std::string(g.fluid) + R"(?{"RES":{"viscosity":true,"conductivity":true,"mixture_critical_enhancement":true}})"));
        if (g.mole_fractions[0] != '\0') {
            std::vector<CoolPropDbl> z;
            for (const std::string& x : strsplit(g.mole_fractions, ';')) {
                z.push_back(strtod(x.c_str(), nullptr));
            }
            AS->set_mole_fractions(z);
        }
        AS->update(PT_INPUTS, g.p_Pa, g.T_K);
        // 1e-8 relative, not 1e-12: the inputs above are rounded to 10 significant figures, so the
        // state itself is only reproduced to about that.  Anything that actually changes the model
        // moves these by far more.
        CHECK(AS->viscosity() == Catch::Approx(g.eta).epsilon(1e-8));
        CHECK(AS->conductivity() == Catch::Approx(g.tc).epsilon(1e-8));

        // How much of the answer the residual term actually carries -- see min_residual_share.
        auto& HEOS = *static_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        const double eta0 = TransportRoutines::viscosity_RES_dilute(HEOS);
        CHECK(std::abs(AS->viscosity() / eta0 - 1) > g.min_residual_share);
    }
}

TEST_CASE("RES on the cubics tracks HEOS without being it", "[RES][transport]") {
    // PR and SRK carry their OWN fits (Yang 2025), regressed against those equations of state, so
    // they should land near the HEOS answer without reproducing it.  Both bounds matter: too far
    // apart means a wrong parameter set or a unit slip, identical means one backend is silently
    // using the other's coefficients.
    for (const char* backend : {"PR", "SRK"}) {
        CAPTURE(backend);
        std::shared_ptr<AbstractState> cub(AbstractState::factory(backend, "Propane"));
        auto ref = heos("Propane", R"(?{"RES":{"viscosity":true,"conductivity":true}})");
        const double T = 1.10 * ref->T_critical(), rho = 1.50 * ref->rhomass_critical();
        ref->update(DmassT_INPUTS, rho, T);
        cub->update(DmassT_INPUTS, rho, T);

        CHECK(cub->viscosity() > 0);
        CHECK(cub->conductivity() > 0);
        CHECK(cub->viscosity() != ref->viscosity());
        CHECK(cub->conductivity() != ref->conductivity());
        CHECK(cub->viscosity() == Catch::Approx(ref->viscosity()).epsilon(0.35));
        CHECK(cub->conductivity() == Catch::Approx(ref->conductivity()).epsilon(0.35));
    }
}

TEST_CASE("RES on a cubic without parameters still reports the missing model", "[RES][transport]") {
    // The cubics have no transport model of their own, so for a fluid outside the RES tables the
    // pre-RES behaviour has to survive: NotImplementedError, not a wrong number and not a RES
    // error.  Dichloroethane is one of only four fluids the cubic list carries that neither paper
    // fitted -- note that Water, the obvious candidate, IS fitted.
    std::shared_ptr<AbstractState> PR(AbstractState::factory("PR", "Dichloroethane"));
    PR->update(PT_INPUTS, 1e5, 300.0);
    CHECK_THROWS_AS(PR->viscosity(), NotImplementedError);
    CHECK_THROWS_AS(PR->conductivity(), NotImplementedError);
}

// ─────────────────────── grid evaluation harness (measurement) ───────────────

/// Evaluate the C++ RES implementation over dev/RES_comparison/grid_points.csv.
///
/// Hidden ([.]), because it is a measurement rather than an assertion: it needs a grid file that
/// dev/RES_grid_build.py writes, it takes minutes, and its output is a CSV for
/// dev/RES_grid_report.py --implementation to compare against the vendored reference code run on
/// the SAME equation of state.  Same model, same EOS, so any deviation there is a defect in this
/// implementation rather than a property of the parameters.
///
///   python dev/RES_grid_build.py
///   ./build_catch/Release/CatchTestRunner.exe "[RES_grid]"
///   python dev/RES_reference_run.py --eos HEOS
///   python dev/RES_grid_report.py --implementation
TEST_CASE("RES grid evaluation harness (measurement)", "[.][RES_grid]") {
    const std::string in_path = std::string(RES_COMPARISON_DIR) + "/grid_points.csv";
    std::ifstream in(in_path.c_str());
    if (!in.good()) {
        FAIL("missing " << in_path << " -- run `python dev/RES_grid_build.py` first");
    }
    const std::string out_path = std::string(RES_COMPARISON_DIR) + "/grid_cpp.csv";
    std::ofstream out(out_path.c_str());
    REQUIRE(out.good());
    out << "fluid,T_K,p_Pa,region,mole_fractions,eos,ok,err,eta,tc,eta0,tc0,rho,s_res,Tc,rhoc,phase\n";
    out.precision(17);

    std::string line;
    std::getline(in, line);  // header
    std::size_t n_ok = 0, n_fail = 0;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        // dev/RES_grid_build.py writes with csv.writer, i.e. CRLF.  Windows text mode strips the
        // CR and POSIX does not, so without this the LAST field of every pure-fluid row is \r"
        // lone carriage return rather than empty off Windows -- which reads as a mole-fraction
        // list, and the harness
        // would then set a zero mole fraction on 3530 of 3570 points instead of leaving the pure
        // fluid alone.  Silent, and wrong only on the platforms CI runs.
        if (!line.empty() && line[line.size() - 1] == '\r') {
            line.erase(line.size() - 1);
        }
        std::vector<std::string> f = strsplit(line, ',');
        if (f.size() < 5) {
            continue;
        }
        const std::string fluid = f[0], region = f[3], fracs = f[4];
        const double T = strtod(f[1].c_str(), nullptr), p = strtod(f[2].c_str(), nullptr);

        std::string err;
        double eta = _HUGE, tc = _HUGE, eta0 = _HUGE, tc0 = _HUGE, rho = _HUGE, s_res = _HUGE, Tc = _HUGE, rhoc = _HUGE;
        int ok = 0, phase = -1;
        try {
            // Both properties on, and the mixture enhancement on as well: the reference code
            // applies it to mixtures, so pinning it here is what makes the two comparable.
            std::shared_ptr<AbstractState> AS(
              AbstractState::factory("HEOS", fluid + R"(?{"RES":{"viscosity":true,"conductivity":true,"mixture_critical_enhancement":true}})"));
            if (!fracs.empty()) {
                std::vector<CoolPropDbl> z;
                for (const std::string& x : strsplit(fracs, ';')) {
                    z.push_back(strtod(x.c_str(), nullptr));
                }
                AS->set_mole_fractions(z);
            }
            AS->update(PT_INPUTS, p, T);
            auto& HEOS = *static_cast<HelmholtzEOSMixtureBackend*>(AS.get());
            eta = AS->viscosity();
            tc = AS->conductivity();
            rho = AS->rhomass();
            s_res = AS->smolar_residual();
            // The dilute terms the model actually used, so the report can weight each point by
            // how much of the property is genuinely residual.  Taken from the model rather than
            // recomputed here, so this cannot drift from what was actually subtracted.
            eta0 = TransportRoutines::viscosity_RES_dilute(HEOS);
            tc0 = TransportRoutines::conductivity_RES_dilute(HEOS);
            if (HEOS.get_components().size() == 1) {
                Tc = AS->T_critical();
                rhoc = AS->rhomass_critical();
            }
            try {
                phase = static_cast<int>(AS->phase());
            } catch (...) {
                phase = -1;
            }
            ok = (ValidNumber(eta) && ValidNumber(tc)) ? 1 : 0;
            if (ok == 0) {
                err = "non-finite";
            }
        } catch (const std::exception& e) {
            err = e.what();
            for (char& c : err) {
                if (c == ',' || c == '\n' || c == '\r') {
                    c = ' ';
                }
            }
            // 300, not 90: the mixture critical-point failures say how many points the solver
            // returned, and that is the whole diagnostic value of the message.  At 90 every one of
            // them was cut at "...but the backend repor".
            if (err.size() > 300) {
                err = err.substr(0, 300);
            }
        }
        ok == 1 ? ++n_ok : ++n_fail;
        out << fluid << ',' << f[1] << ',' << f[2] << ',' << region << ',' << fracs << ",HEOS," << ok << ',' << err << ',';
        auto num = [&](double v) {
            if (ValidNumber(v) && v != _HUGE) {
                out << v;
            }
            out << ',';
        };
        num(eta);
        num(tc);
        num(eta0);
        num(tc0);
        num(rho);
        num(s_res);
        num(Tc);
        num(rhoc);
        out << (phase >= 0 ? std::to_string(phase) : std::string()) << '\n';
    }
    out.close();
    UNSCOPED_INFO("wrote " << out_path << ": " << n_ok << " evaluated, " << n_fail << " failed");
    CHECK(n_ok > 0);
}

#endif  // ENABLE_CATCH
