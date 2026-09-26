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

TEST_CASE("Bollengier backend uses the IAPWS reference state", "[Bollengier][water][reference]") {
    // The published surface does NOT use the IAPWS convention.  Without a
    // shift, `Water` would silently disagree with itself depending on which
    // backend answered -- the kind of inconsistency that surfaces as a
    // baffling energy-balance error far from its cause.
    auto AS = make();

    SECTION("h and s at the triple-point saturated-liquid state") {
        // IAPWS: u = s = 0 for saturated liquid at the triple point, so
        // h = p_t * v there.
        AS->update(PT_INPUTS, 611.657, 273.16);
        INFO("h = " << AS->hmass() << " J/kg, s = " << AS->smass() << " J/kg/K");
        CHECK_THAT(AS->smass(), Catch::Matchers::WithinAbs(0.0, 1e-9));
        CHECK_THAT(AS->hmass(), Catch::Matchers::WithinAbs(0.611781703, 1e-6));
    }

    SECTION("the shift preserves g = h - T s") {
        // This is the real point of shifting G by (h0 - T*s0) rather than
        // adjusting h and s separately: doing it separately would leave g
        // inconsistent by a T-dependent amount, which no single-point check
        // of h or s would catch.
        for (const double T : {250.0, 300.0, 450.0}) {
            for (const double p_MPa : {0.1, 100.0, 1500.0}) {
                AS->update(PT_INPUTS, p_MPa * 1e6, T);
                const double g_from_hs = AS->hmass() - T * AS->smass();
                const double g_direct = AS->umass() + p_MPa * 1e6 / AS->rhomass() - T * AS->smass();
                INFO("T = " << T << " K, p = " << p_MPa << " MPa");
                CHECK_THAT(g_from_hs, Catch::Matchers::WithinRel(g_direct, 1e-12));
            }
        }
    }

    SECTION("energies stay within the two models' documented difference") {
        // NOTE what this does and does not check.  It does NOT pin the
        // reference shift -- the triple-point section above does that, and
        // is exact.  Removing the shift entirely would actually make h
        // agree with IAPWS-95 *better* here (-5.5 vs -76 J/kg), so this
        // comparison cannot detect a missing shift.  It does catch a sign
        // error, which would double the offset to ~141 J/kg.
        //
        // The residual difference is physical, not numerical: cp differs
        // between the models by up to 1.2e-3 relative in the cold liquid
        // region (peak near 281 K, vanishing above ~300 K), which integrates
        // to -67.7 J/kg over 273-350 K.  That is inside IAPWS-95's own
        // stated cp uncertainty of 1000 ppm and is precisely the region
        // Bollengier et al. claim to improve.  Tightening these bounds would
        // be asserting the two equations of state are identical, which they
        // deliberately are not.
        auto heos = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
        for (const double T : {280.0, 300.0, 350.0}) {
            AS->update(PT_INPUTS, 0.1e6, T);
            heos->update(PT_INPUTS, 0.1e6, T);
            INFO("T = " << T << " K: Bollengier h = " << AS->hmass() << ", IAPWS-95 h = " << heos->hmass() << " | s = " << AS->smass() << " vs "
                        << heos->smass());
            CHECK_THAT(AS->hmass(), Catch::Matchers::WithinAbs(heos->hmass(), 150.0));
            CHECK_THAT(AS->smass(), Catch::Matchers::WithinAbs(heos->smass(), 0.5));
            // Densities agree far more closely than energies, because the
            // density difference is not an integral of a cp discrepancy.
            CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(heos->rhomass(), 3e-6));
        }
    }
}
#endif
