// Tests for the neon melting line (bd CoolProp-r1w7.2).
//
// The neon melting curve was replaced (SantamariaPerez-PRB-2010 Simon fit, which
// was high by ~11-29x and whose p_min sat ~14.6 MPa above the triple pressure) by
// a polynomial-in-theta fit to the solid-liquid (crystal phase C1) melting data
// evaluated in the NIST ThermoData Engine, anchored at the triple point and
// constrained there to the Clausius-Clapeyron slope dp/dT = dS_fus/dV_fus =
// 6.28 MPa/K.  Run explicitly:  ./CatchTestRunner "[melting]"

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <memory>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

// Check values pinning the replaced neon melting line so it cannot silently
// drift.  Reference temperatures are obtained by inverting the fitted curve and
// reproduce the experimental melting temperatures to ~1 K RMS.
TEST_CASE("Neon melting line check values (TDE fit, C-C constrained)", "[melting]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Neon"));
    CHECK(AS->melting_line(CoolProp::iT, CoolProp::iP, 1.0e8) == Catch::Approx(37.083).epsilon(1e-4));
    CHECK(AS->melting_line(CoolProp::iT, CoolProp::iP, 1.0e9) == Catch::Approx(109.769).epsilon(1e-4));
}

// Regression: the prior neon curve's p_min (~14.6 MPa) sat far above the triple
// pressure, so near-triple consistency-grid points threw "unable to calculate
// melting line T(p) ... bounds are ...".  The replacement is anchored at the
// triple point (p_min == triple pressure), so these states now resolve.
TEST_CASE("Neon phase determination near the triple pressure", "[melting]") {
    const double p = 43851.4;  // Pa, just above the triple pressure
    const double T = GENERATE(24.56, 100.0);
    CAPTURE(T, p);
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Neon"));
    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, p, T));
}

#endif  // ENABLE_CATCH
