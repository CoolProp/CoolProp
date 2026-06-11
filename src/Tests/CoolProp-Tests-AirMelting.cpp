// Tests for the air (pseudo-pure) melting line (bd CoolProp-r1w7 Tier 3).
//
// Air lacked a melting line, so its high-pressure P+{H,S,U} consistency points
// failed ("Inputs in Brent [...] do not bracket the root"): the HSU_P flash
// floored its T-search at the EOS triple temperature and probed the cold dense
// region.  The melting line of Lemmon et al. (2000) -- the same paper as air's
// EOS, ported from REFPROP's AIR.PPF -- is added so the flash floors the search
// at Tmelt(p).  Run explicitly:  ./CatchTestRunner "[melting]"

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <memory>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

// Check values pinning the air melting line (Simon form, reproduces REFPROP's
// Lemmon-2000 curve to <0.02%).  T(p) obtained by inverting the fitted curve.
TEST_CASE("Air melting line check values (Lemmon-2000)", "[melting]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Air"));
    CHECK(AS->melting_line(CoolProp::iT, CoolProp::iP, 1.0e8) == Catch::Approx(75.920).epsilon(1e-4));
    CHECK(AS->melting_line(CoolProp::iT, CoolProp::iP, 1.0e9) == Catch::Approx(167.875).epsilon(1e-4));
}

// Regression: a genuine fluid state (T > Tmelt(p)) at p ~ 1 GPa that previously
// failed the P+{H,S,U} flash now round-trips, because the melting line floors
// the HSU_P T-search above the cold dense region.
TEST_CASE("Air high-pressure P+X reliability (was failing)", "[melting]") {
    const double T = 210.01538461538462;  // K  (> Tmelt(1.035 GPa) ~ 170 K)
    const double p = 1035411979.2369467;  // Pa
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Air"));
    auto wrk = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Air"));
    REQUIRE_NOTHROW(ref->update(CoolProp::PT_INPUTS, p, T));
    const double h = ref->hmolar(), s = ref->smolar(), u = ref->umolar(), rho = ref->rhomolar();

    SECTION("PH") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::HmolarP_INPUTS, h, p));
        CHECK(wrk->T() == Catch::Approx(T).epsilon(1e-5));
        CHECK(wrk->rhomolar() == Catch::Approx(rho).epsilon(1e-5));
    }
    SECTION("PS") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PSmolar_INPUTS, p, s));
        CHECK(wrk->T() == Catch::Approx(T).epsilon(1e-5));
    }
    SECTION("PU") {
        REQUIRE_NOTHROW(wrk->update(CoolProp::PUmolar_INPUTS, p, u));
        CHECK(wrk->T() == Catch::Approx(T).epsilon(1e-5));
    }
}

#endif  // ENABLE_CATCH
