// Tests for the Bollengier et al. (2019) liquid-water backend.
//
// IAPWS-95 is used as a reference only BELOW 100 MPa, where the two models
// are expected to agree.  Above that they diverge by design -- that
// divergence is the reason this backend exists -- so agreement there would
// be a red flag, not a pass.

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
constexpr double kPmaxMPa = 2300.6;
constexpr double kTminK = 239.0;
constexpr double kTmaxK = 501.0;

std::shared_ptr<AbstractState> make() {
    return std::shared_ptr<AbstractState>(AbstractState::factory("BOLLENGIER", "Water"));
}
std::shared_ptr<AbstractState> iapws() {
    return std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Water"));
}
}  // namespace

TEST_CASE("Bollengier backend agrees with IAPWS-95 where both are well constrained", "[Bollengier][water]") {
    auto AS = make();
    auto ref = iapws();
    // T >= 300 K and p <= 50 MPa.  The models are NOT uniformly close below
    // 100 MPa: the disagreement is driven by temperature as much as
    // pressure, and grows sharply toward the cold end (see the next case).
    // Checked against IAPWS-95 live rather than against stored numbers, so
    // this cannot drift into self-reference.
    for (const double p_MPa : {0.101325, 0.1, 10.0, 50.0}) {
        for (const double T : {300.0, 350.0}) {
            AS->update(PT_INPUTS, p_MPa * 1e6, T);
            ref->update(PT_INPUTS, p_MPa * 1e6, T);
            INFO("p = " << p_MPa << " MPa, T = " << T << " K");
            CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(ref->rhomass(), 5e-6));
            // cp, cv and w had no coverage at all before; all three come from
            // second partials, where a sign or unit error is most likely and
            // least visible in density.
            CHECK_THAT(AS->cpmass(), Catch::Matchers::WithinRel(ref->cpmass(), 2e-3));
            CHECK_THAT(AS->cvmass(), Catch::Matchers::WithinRel(ref->cvmass(), 2e-3));
            CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(ref->speed_sound(), 3e-3));
            // Stability, independent of any reference: these must hold for
            // any physically sensible liquid.
            CHECK(AS->cvmass() < AS->cpmass());
            CHECK(AS->cvmass() > 0.0);
        }
    }
}

TEST_CASE("Bollengier backend departs from IAPWS-95 in the cold liquid region", "[Bollengier][water]") {
    // At 273.16 K the two models diverge rapidly with pressure -- ~7 ppm in
    // density at 10 MPa, ~46 ppm at 50 MPa -- while at 300 K the same
    // pressures agree to ~1 ppm.  That is the supercooled / cold-liquid
    // region the paper sets out to improve, and IAPWS-95's own uncertainty
    // is larger there.  Asserting tight agreement would be asserting the
    // two models are the same; this pins the divergence instead, so a
    // regression toward IAPWS-95 would be noticed.
    auto AS = make();
    auto ref = iapws();
    AS->update(PT_INPUTS, 50e6, 273.16);
    ref->update(PT_INPUTS, 50e6, 273.16);
    const double rel = AS->rhomass() / ref->rhomass() - 1.0;
    INFO("relative density difference at 273.16 K, 50 MPa = " << rel);
    CHECK(rel > 1e-5);
    CHECK(rel < 1e-4);
}

TEST_CASE("Bollengier backend diverges from IAPWS-95 above 100 MPa, as published", "[Bollengier][water]") {
    // At 100 MPa the models part company -- roughly 12 ppm in density, which
    // is the paper's claimed improvement region, not an error.  Asserting
    // tight agreement here would be asserting the two models are the same.
    auto AS = make();
    auto ref = iapws();
    AS->update(PT_INPUTS, 100e6, 300.0);
    ref->update(PT_INPUTS, 100e6, 300.0);
    const double rel = AS->rhomass() / ref->rhomass() - 1.0;
    INFO("relative density difference at 100 MPa, 300 K = " << rel);
    CHECK(std::abs(rel) > 3e-6);  // genuinely different
    CHECK(std::abs(rel) < 1e-4);  // but not wildly so
}

TEST_CASE("Bollengier backend reaches the domain IAPWS-95 refuses", "[Bollengier][water]") {
    // The actual reason for this backend.  Densities are characterisation
    // values from this implementation, not an independent source -- the
    // independent constraint is the sub-100-MPa IAPWS-95 comparison above,
    // plus the stability and consistency checks here.  Their job is to catch
    // an unintended change in the surface, e.g. a regenerated coefficient
    // set.
    auto AS = make();
    auto ref = iapws();
    const double pts[][3] = {
      {250.0, 200.0, 1090.4850002876647},
      {280.0, 1500.0, 1317.8261518857066},
      {300.0, 2000.0, 1360.5338887725438},
      {500.0, 2300.0, 1298.2923138211172},
    };
    for (const auto& q : pts) {
        INFO("T = " << q[0] << " K, p = " << q[1] << " MPa");
        AS->update(PT_INPUTS, q[1] * 1e6, q[0]);
        REQUIRE(std::isfinite(AS->rhomass()));
        CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(q[2], 1e-9));
        CHECK(AS->cvmass() < AS->cpmass());
        CHECK(std::isfinite(AS->speed_sound()));
        // IAPWS-95 refuses these outright.  Asserted through AbstractState,
        // not PropsSI: the C++ PropsSI catches and returns NaN with an error
        // string, so CHECK_THROWS on it would pass either way.
        CHECK_THROWS(ref->update(PT_INPUTS, q[1] * 1e6, q[0]));
    }
}

TEST_CASE("Bollengier backend caches nothing across updates", "[Bollengier][water][cache]") {
    // AbstractState caches speed_sound/hmolar/cpmolar and friends.  Without
    // clear() in update() they freeze at their first value, and it reaches
    // PropsSI -- successive queries return the PREVIOUS query's answer.
    // IF97Backend carries a comment describing this exact symptom.
    auto reused = make();
    reused->update(PT_INPUTS, 1e5, 300.0);
    (void)reused->speed_sound();
    (void)reused->hmolar();
    (void)reused->cpmolar();
    reused->update(PT_INPUTS, 2000e6, 250.0);

    auto fresh = make();
    fresh->update(PT_INPUTS, 2000e6, 250.0);

    INFO("reused w = " << reused->speed_sound() << ", fresh w = " << fresh->speed_sound());
    CHECK_THAT(reused->speed_sound(), Catch::Matchers::WithinRel(fresh->speed_sound(), 1e-12));
    CHECK_THAT(reused->hmolar(), Catch::Matchers::WithinRel(fresh->hmolar(), 1e-12));
    CHECK_THAT(reused->cpmolar(), Catch::Matchers::WithinRel(fresh->cpmolar(), 1e-12));
    CHECK_THAT(reused->smolar(), Catch::Matchers::WithinRel(fresh->smolar(), 1e-12));
}

TEST_CASE("Bollengier backend exposes T and p", "[Bollengier][water]") {
    // The base class's _p/_T back T(), p() and keyed_output(); if a backend
    // only stores its own copies they stay at -_HUGE and every consumer that
    // does not go through PropsSI's input short-circuit sees -inf.
    auto AS = make();
    AS->update(PT_INPUTS, 12.5e6, 321.0);
    CHECK_THAT(AS->p(), Catch::Matchers::WithinRel(12.5e6, 1e-12));
    CHECK_THAT(AS->T(), Catch::Matchers::WithinRel(321.0, 1e-12));
    CHECK_THAT(AS->keyed_output(iP), Catch::Matchers::WithinRel(12.5e6, 1e-12));
    CHECK_THAT(AS->keyed_output(iT), Catch::Matchers::WithinRel(321.0, 1e-12));
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
    SECTION("a rejected pair does not leave the previous state readable") {
        // _valid must be cleared BEFORE the input-pair check, not after it.
        AS->update(PT_INPUTS, 1e5, 300.0);
        REQUIRE(std::isfinite(AS->rhomass()));
        CHECK_THROWS(AS->update(DmassT_INPUTS, 1000.0, 300.0));
        CHECK_THROWS(AS->rhomass());
    }
}

TEST_CASE("Bollengier backend range guard throws, and is pinned open", "[Bollengier][water][validation]") {
    auto AS = make();
    SECTION("outside the published domain it throws, it does not return NaN") {
        CHECK_THROWS_AS(AS->update(PT_INPUTS, (kPmaxMPa + 1.0) * 1e6, 300.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, -1.0e6, 300.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 10e6, kTminK - 1.0), CoolProp::ValueError);
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 10e6, kTmaxK + 1.0), CoolProp::ValueError);
        // A negative pressure small enough to underflow to -0.0 under the
        // MPa conversion must still be refused.
        CHECK_THROWS_AS(AS->update(PT_INPUTS, -std::numeric_limits<double>::denorm_min(), 300.0), CoolProp::ValueError);
    }
    SECTION("the advertised bounds are the paper's, to within a ULP") {
        // Pins the bounds against being widened OR narrowed, which
        // CHECK_THROWS_AS alone cannot do -- the spline throws the same
        // exception type, so a widened backend guard would still "throw"
        // and pass.  Comparing the reported limits to the paper's literals
        // catches it.
        CHECK_THAT(AS->pmax(), Catch::Matchers::WithinRel(kPmaxMPa * 1e6, 1e-9));
        CHECK_THAT(AS->Tmin(), Catch::Matchers::WithinRel(kTminK, 1e-9));
        CHECK_THAT(AS->Tmax(), Catch::Matchers::WithinRel(kTmaxK, 1e-9));
    }
    SECTION("just inside each bound, in the stable region, it must NOT throw") {
        const double eps = 1e-6;
        for (const double T : {kTminK + eps, kTmaxK - eps}) {
            for (const double P : {eps, 1000.0}) {
                INFO("T = " << T << " K, p = " << P << " MPa");
                REQUIRE_NOTHROW(AS->update(PT_INPUTS, P * 1e6, T));
                CHECK(std::isfinite(AS->rhomass()));
            }
        }
        // Warm end of the top pressure edge, which is well inside the
        // liquid field.
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, (kPmaxMPa - eps) * 1e6, 400.0));
    }
}

TEST_CASE("Bollengier backend refuses the ice-field corner rather than returning NaN", "[Bollengier][water][nan]") {
    // The advertised domain is a RECTANGLE, but the model's real validity is
    // bounded by the melting curve.  The cold high-pressure corner lies deep
    // in the ice VI/VII field, and there the fitted surface has
    // (dv/dP)_T >= 0, which makes the speed of sound NaN and cv negative or
    // absurd.  Returning those would be precisely the silent non-finite
    // propagation TensorBSpline2D refuses to do, one layer up.
    auto AS = make();
    CHECK_THROWS_AS(AS->update(PT_INPUTS, kPmaxMPa * 1e6, kTminK), CoolProp::ValueError);
    CHECK_THROWS_WITH(AS->update(PT_INPUTS, 2290.0e6, 240.0), Catch::Matchers::ContainsSubstring("not thermodynamically stable"));
    // The same pressure at a warmer temperature is fine, so this is not a
    // pressure-bound problem in disguise.
    REQUIRE_NOTHROW(AS->update(PT_INPUTS, 2290.0e6, 400.0));
    CHECK(std::isfinite(AS->speed_sound()));
}

TEST_CASE("Bollengier backend reports the published reference state unshifted", "[Bollengier][water][reference]") {
    // This model does not use the IAPWS convention, and this backend does
    // not re-anchor it: h and s are defined only up to a constant, and
    // shifting would make CoolProp report numbers that are not the published
    // model's.  Pinning the values means an accidental shift -- or a
    // regenerated coefficient set with a different datum -- fails loudly.
    auto AS = make();
    AS->update(PT_INPUTS, 611.657, 273.16);
    INFO("h = " << AS->hmass() << " J/kg, s = " << AS->smass() << " J/kg/K");
    CHECK_THAT(AS->hmass(), Catch::Matchers::WithinRel(71.23178743304095, 1e-9));
    CHECK_THAT(AS->smass(), Catch::Matchers::WithinRel(0.25817666496269626, 1e-9));
}

TEST_CASE("Bollengier backend property map satisfies the Maxwell relations", "[Bollengier][water][consistency]") {
    // A real consistency check, by finite differences of g = h - Ts.
    //
    // NOTE: comparing (h - T*s) against (u + p*v - T*s) is NOT a consistency
    // check -- it is an algebraic identity in h and s, since u is defined as
    // h - p*v.  It passes even if h is wrong by an arbitrary constant.  The
    // derivatives below actually constrain s and v against g.
    auto AS = make();
    auto g = [&](double p_Pa, double T) {
        AS->update(PT_INPUTS, p_Pa, T);
        return AS->hmass() - T * AS->smass();
    };
    for (const double T : {260.0, 300.0, 450.0}) {
        for (const double p_MPa : {1.0, 100.0, 1000.0}) {
            const double p_Pa = p_MPa * 1e6;
            const double dT = 1e-3 * T;
            const double dp = 1e-4 * p_Pa;
            // s = -(dg/dT)_p
            const double s_fd = -(g(p_Pa, T + dT) - g(p_Pa, T - dT)) / (2 * dT);
            // v = (dg/dp)_T
            const double v_fd = (g(p_Pa + dp, T) - g(p_Pa - dp, T)) / (2 * dp);
            AS->update(PT_INPUTS, p_Pa, T);
            INFO("T = " << T << " K, p = " << p_MPa << " MPa");
            // 5e-5 is set by central-difference TRUNCATION error, not by
            // any looseness in the model: the step has to stay large enough
            // to avoid cancellation noise in g, which leaves an O(dT^2)
            // residual.  A sign error or a wrong unit factor is off by
            // orders of magnitude, so this still discriminates sharply.
            // Relative OR absolute: under compression s passes close to
            // zero (~15 J/kg/K at 1000 MPa on this datum), where a purely
            // relative tolerance measures the wrong thing.
            CHECK_THAT(s_fd, Catch::Matchers::WithinRel(AS->smass(), 5e-5) || Catch::Matchers::WithinAbs(AS->smass(), 1e-2));
            CHECK_THAT(v_fd, Catch::Matchers::WithinRel(1.0 / AS->rhomass(), 5e-5));
        }
    }
}
#endif
