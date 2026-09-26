// Tests for the Bollengier et al. (2019) liquid-water backend.
//
// Reference values are the paper's own domain where possible; IAPWS-95 is
// used only below 100 MPa, where the two are expected to agree.  Above that
// they diverge by design -- that divergence is the reason this backend
// exists -- so agreement there would be a red flag, not a pass.

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/Exceptions.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <memory>
#    include <vector>

using namespace CoolProp;

namespace {
// Stated in Bollengier et al. (2019), NOT read back from the coefficient
// file.  A test that took its bounds from the data it is validating could
// not detect a truncated or mis-parsed coefficient set.
constexpr double kPminMPa = 0.0;
constexpr double kPmaxMPa = 2300.6;
constexpr double kTminK = 239.0;
constexpr double kTmaxK = 501.0;

std::shared_ptr<AbstractState> make() {
    return std::shared_ptr<AbstractState>(AbstractState::factory("BOLLENGIER", "Water"));
}
}  // namespace

TEST_CASE("Bollengier backend agrees with IAPWS-95 below 100 MPa", "[Bollengier][water]") {
    auto AS = make();
    // {T [K], p [MPa], rho_IAPWS95 [kg/m3]}.  Deviations are sub-3-ppm here;
    // see COO-44 for the measured table.
    const double ref[][3] = {
      {273.16, 0.101325, 999.84376},
      {300.0, 0.1, 996.55634},
      {300.0, 50.0, 1017.84630},
      {300.0, 100.0, 1037.19149},
    };
    for (const auto& r : ref) {
        INFO("T = " << r[0] << " K, p = " << r[1] << " MPa");
        AS->update(PT_INPUTS, r[1] * 1e6, r[0]);
        CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(r[2], 3e-6));
    }
}

TEST_CASE("Bollengier backend reaches the domain IAPWS-95 refuses", "[Bollengier][water]") {
    // The actual reason for this backend: CoolProp's IAPWS-95 throws
    // "below Tmelt(p)" across this whole region.
    auto AS = make();
    const double pts[][3] = {
      {250.0, 200.0, 1090.47212},
      {280.0, 1500.0, 1317.81697},
      {300.0, 2000.0, 1360.50854},
      {500.0, 2300.0, 1298.30602},
    };
    for (const auto& p : pts) {
        INFO("T = " << p[0] << " K, p = " << p[1] << " MPa");
        AS->update(PT_INPUTS, p[1] * 1e6, p[0]);
        const double rho = AS->rhomass();
        REQUIRE(std::isfinite(rho));
        CHECK_THAT(rho, Catch::Matchers::WithinRel(p[2], 1e-6));
        // The claim this backend rests on: IAPWS-95 refuses these outright.
        // Asserted through AbstractState rather than PropsSI, because the
        // C++ PropsSI catches and returns NaN with an error string instead
        // of propagating -- so CHECK_THROWS on PropsSI would pass whether
        // or not IAPWS-95 actually rejected the point.
        auto heos = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
        CHECK_THROWS(heos->update(PT_INPUTS, p[1] * 1e6, p[0]));
    }
}

TEST_CASE("Bollengier backend accepts only PT inputs", "[Bollengier][water][validation]") {
    auto AS = make();
    SECTION("positive control: PT works") {
        // Without this, an update() that threw unconditionally would satisfy
        // every other section here.
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, 10e6, 300.0));
        CHECK(std::isfinite(AS->rhomass()));
    }
    SECTION("every other input pair throws") {
        CHECK_THROWS_AS(AS->update(DmassT_INPUTS, 1000.0, 300.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(HmassP_INPUTS, 1e5, 10e6), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PSmass_INPUTS, 10e6, 100.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(HmassSmass_INPUTS, 1e5, 100.0), CoolProp::ValueError);
    }
    SECTION("saturation calls throw: this model has no vapour branch") {
        CHECK_THROWS(AS->update(PQ_INPUTS, 1e5, 0.0));
        CHECK_THROWS(AS->update(QT_INPUTS, 0.0, 300.0));
    }
}

TEST_CASE("Bollengier backend range guard throws, and is pinned open", "[Bollengier][water][validation]") {
    auto AS = make();
    SECTION("outside the published domain it throws, it does not return NaN") {
        CHECK_THROWS_AS(AS->update(PT_INPUTS, (kPmaxMPa + 1.0) * 1e6, 300.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, -1.0e6, 300.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 10e6, kTminK - 1.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 10e6, kTmaxK + 1.0), CoolProp::ValueError);
    }
    SECTION("just inside each bound it must NOT throw") {
        // Pins the guard open: an over-strict domain would fail here.
        const double eps = 1e-6;
        for (const double T : {kTminK + eps, kTmaxK - eps}) {
            for (const double P : {kPminMPa + eps, kPmaxMPa - eps}) {
                INFO("T = " << T << " K, p = " << P << " MPa");
                REQUIRE_NOTHROW(AS->update(PT_INPUTS, P * 1e6, T));
                CHECK(std::isfinite(AS->rhomass()));
            }
        }
    }
}

TEST_CASE("Bollengier backend reports the published reference state unshifted", "[Bollengier][water][reference]") {
    // This model does NOT use the IAPWS convention, and this backend does
    // not re-anchor it.  h and s are defined only up to a constant;
    // silently shifting them would make CoolProp report numbers that are
    // not the published model's, so anyone checking against the paper's own
    // tables would see an unexplained offset.  Reference state is a
    // user-level concern in CoolProp (set_reference_stateS /
    // set_reference_stateD).
    auto AS = make();

    SECTION("the model's own reference state, documented not corrected") {
        // At the triple-point saturated-liquid state IAPWS-95 has
        // h = 0.6117817 J/kg, s = 0.  This model does not, and that is
        // expected.  Pinning the values here means an accidental shift --
        // or a regenerated coefficient set with a different datum -- shows
        // up as a test failure rather than as a silent change in output.
        AS->update(PT_INPUTS, 611.657, 273.16);
        INFO("h = " << AS->hmass() << " J/kg, s = " << AS->smass() << " J/kg/K");
        CHECK_THAT(AS->hmass(), Catch::Matchers::WithinRel(71.0954117, 1e-7));
        CHECK_THAT(AS->smass(), Catch::Matchers::WithinRel(0.257680474, 1e-7));
    }

    SECTION("g = h - T s holds across the domain") {
        // Guards the property map itself: h, s and u must stay mutually
        // consistent regardless of where the energy datum sits.
        for (const double T : {250.0, 300.0, 450.0}) {
            for (const double p_MPa : {0.1, 100.0, 1500.0}) {
                AS->update(PT_INPUTS, p_MPa * 1e6, T);
                const double g_from_hs = AS->hmass() - T * AS->smass();
                const double g_from_u = AS->umass() + p_MPa * 1e6 / AS->rhomass() - T * AS->smass();
                INFO("T = " << T << " K, p = " << p_MPa << " MPa");
                CHECK_THAT(g_from_hs, Catch::Matchers::WithinRel(g_from_u, 1e-12));
            }
        }
    }

    SECTION("energies remain comparable with IAPWS-95 despite the different datum") {
        // Not a reference-state check -- it cannot be one, since the datum
        // differs by construction.  It bounds the two models' genuine
        // disagreement: cp differs by up to 1.2e-3 relative in the cold
        // liquid region, peaking near 281 K and vanishing above ~300 K,
        // which integrates to a few tens of J/kg.  Densities agree far more
        // closely, because that difference is not an integral of a cp
        // discrepancy.
        auto heos = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
        for (const double T : {280.0, 300.0, 350.0}) {
            AS->update(PT_INPUTS, 0.1e6, T);
            heos->update(PT_INPUTS, 0.1e6, T);
            INFO("T = " << T << " K: Bollengier h = " << AS->hmass() << ", IAPWS-95 h = " << heos->hmass() << " | s = " << AS->smass() << " vs "
                        << heos->smass());
            CHECK_THAT(AS->hmass(), Catch::Matchers::WithinAbs(heos->hmass(), 150.0));
            CHECK_THAT(AS->smass(), Catch::Matchers::WithinAbs(heos->smass(), 0.5));
            CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(heos->rhomass(), 3e-6));
        }
    }
}
#endif
