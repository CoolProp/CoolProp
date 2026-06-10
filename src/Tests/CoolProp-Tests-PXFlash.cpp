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

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

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

// ---------------------------------------------------------------------------
// High-pressure quantum-fluid P+X reliability (bd CoolProp-r1w7.1).
//
// At p >> psat_max (hundreds of MPa) with T just above Tc, the HSU_P flash's
// T-search bracket straddles Tc.  Without a melting line the search floor is
// the bare EOS triple temperature, so trial T < Tc probe the cold solid-region
// isotherm (a van der Waals loop pinned near p~0); solver_rho_Tp cannot find the
// high-density root there and threw "Inputs in Brent [...] do not bracket the
// root", aborting the flash before it could reach the valid T > Tc solution.
// This was the largest exception bucket in the HEOS consistency report and is
// perfectly correlated with a missing melting line (Hydrogen / ParaHydrogen,
// which HAVE one, do not fail; OrthoHydrogen and the Deuterium isotopes, which
// did not, failed ~684-702 times each).  The fix gives these fluids a melting
// line so the flash floors its T-search at Tmelt(p), above the loop.  Points are
// taken from the per-fluid consistency CSVs.
// ---------------------------------------------------------------------------
TEST_CASE("P+X high-pressure quantum-fluid reliability (was failing)", "[PXflash]") {
    struct Case
    {
        const char* fluid;
        double T;  // K  (> Tc, but the flash's T-search dips below Tc)
        double p;  // Pa (>> psat_max)
    };
    // Genuine fluid states (T > Tmelt(p)) at p >> psat_max that failed in the
    // consistency report.  Lower-T grid points at the same pressure are below the
    // melting line (solid) and are correctly rejected, not solved.  Adding the
    // melting line (Datchi-2000 for H2; Diatschenko-1985 for the D2 isotopes)
    // floors the HSU_P T-search at Tmelt(p) so these no longer fail.
    const Case c =
      GENERATE(Case{"OrthoHydrogen", 90.86892307692308, 404007839.5225975}, Case{"Deuterium", 79.32923076923078, 448875035.8934926},
               Case{"OrthoDeuterium", 79.32923076923078, 448875035.8934926}, Case{"ParaDeuterium", 79.32923076923078, 448875035.8934926});
    CAPTURE(c.fluid, c.T, c.p);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));

    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, c.p, c.T));
    const double h = ref->hmolar();
    const double s = ref->smolar();
    const double u = ref->umolar();
    const double rho = ref->rhomolar();
    REQUIRE(std::isfinite(rho));

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
