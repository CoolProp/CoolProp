// Focused reliability regression test for the P+{H,S,U} single-phase flash
// (HSU_P_flash_singlephase_Brent in src/Backends/Helmholtz/FlashRoutines.cpp).
//
// Runs in the default suite (tag [PXflash], NOT [.]-hidden).  Run explicitly:
//   ./CatchTestRunner "[PXflash]"
//
// Background (bd CoolProp-0nx): the solver left HEOS._p corrupted when the inner
// PT density solve threw mid-Newton; the exception-fallback guard then branched
// on the garbage pressure and re-threw instead of running the 2-D Newton
// fallback.  Two supercritical-cold (p > p_crit) Nitrogen/PS states failed.
// The fix saves p_target at entry and lets the fallback engage for p > p_crit.
//
// Input-pair argument order (verified against include/DataStructures.h):
//   HmolarP_INPUTS : update(HmolarP_INPUTS, h_J_mol, p_Pa)   (H first, P second)
//   PSmolar_INPUTS : update(PSmolar_INPUTS, p_Pa, s_J_mol_K) (P first, S second)
//   PUmolar_INPUTS : update(PUmolar_INPUTS, p_Pa, u_J_mol)   (P first, U second)

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <memory>
#    include <string>

#    include "AbstractState.h"
#    include "DataStructures.h"

// ---------------------------------------------------------------------------
// The two previously-failing Nitrogen states.  Both are supercritical-cold
// (p > p_crit ~ 3.396 MPa), the compressed-liquid corner.  Before the fix the
// P,S round-trip threw; the fix must make them solve and recover T and rho.
// ---------------------------------------------------------------------------
TEST_CASE("P+X supercritical-cold Nitrogen reliability (was failing)", "[PXflash]") {
    struct State
    {
        double T;  // K
        double p;  // Pa
    };
    const State pt = GENERATE(State{86.35, 4.754e6}, State{90.20, 7.768e6});
    CAPTURE(pt.T, pt.p);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));

    // Establish the reference state at (T, p) and read its entropy.
    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, pt.p, pt.T));
    const double s = ref->smolar();
    const double rho = ref->rhomolar();
    REQUIRE(std::isfinite(s));
    REQUIRE(std::isfinite(rho));

    // The P,S flash that used to throw must now solve and round-trip.
    REQUIRE_NOTHROW(wrk->update(CoolProp::PSmolar_INPUTS, pt.p, s));
    CHECK(wrk->T() == Catch::Approx(pt.T).epsilon(1e-5));
    CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
}

// ---------------------------------------------------------------------------
// Small, fast P+X round-trip guard across a few fluids/phases.  Sets a known
// (T, p) state, reads h/s/u, then runs each of the PH / PS / PU flashes and
// confirms T and rho recover.  Keeps CI cost low (a handful of points).
// ---------------------------------------------------------------------------
TEST_CASE("P+X single-phase round-trips (PH/PS/PU)", "[PXflash]") {
    struct Case
    {
        const char* fluid;
        double T;  // K
        double p;  // Pa
        const char* label;
    };
    const Case c = GENERATE(Case{"Water", 600.0, 1.0e5, "Water gas"},              // superheated vapor
                            Case{"Water", 320.0, 5.0e6, "Water liquid"},           // compressed liquid
                            Case{"CarbonDioxide", 320.0, 1.0e7, "CO2 supercrit"},  // p,T > critical
                            Case{"R134a", 250.0, 3.0e5, "R134a liquid"},           // subcooled liquid
                            Case{"R134a", 350.0, 2.0e5, "R134a gas"});             // superheated vapor
    CAPTURE(c.label, c.fluid, c.T, c.p);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));

    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, c.p, c.T));
    const double h = ref->hmolar();
    const double s = ref->smolar();
    const double u = ref->umolar();
    const double rho = ref->rhomolar();
    REQUIRE(std::isfinite(h));
    REQUIRE(std::isfinite(s));
    REQUIRE(std::isfinite(u));

    SECTION("PH") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::HmolarP_INPUTS, h, c.p));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
    SECTION("PS") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PSmolar_INPUTS, c.p, s));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
    SECTION("PU") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PUmolar_INPUTS, c.p, u));
        CHECK(wrk->T() == Catch::Approx(c.T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
}

#endif  // ENABLE_CATCH
