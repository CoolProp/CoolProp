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
#    include <limits>
#    include <memory>
#    include <vector>

using namespace CoolProp;

namespace {
// Stated in Bollengier et al. (2019), NOT read back from the coefficient
// file.  A test that took its bounds from the data it is validating could
// not detect a truncated or mis-parsed coefficient set.
constexpr double kPmaxMPa = 2300.6;
// The paper's stated range (title and abstract), NOT the knot span --
// the knots run marginally wider, but that is a property of the fit.
constexpr double kTminK = 240.0;
constexpr double kTmaxK = 500.0;

// The box the backend excludes inside its published domain.  Duplicated
// from the backend rather than imported, so a change to either side has
// to be made deliberately on both -- and the comparisons below must stay
// IDENTICAL to the backend's (>= on p, <= on T), or the sweep's
// accepted == expected stops meaning anything.
//
// Measured on the committed coefficients, the states the fit gets wrong
// all lie at p >= 1526.4 MPa and T <= 251.6 K: thermodynamically
// inadmissible (cv <= 0 or (dv/dP)_T >= 0, with a non-finite w) over
// p in [1833.5, 2300.6] MPa and T in [240.0, 249.6] K, and admissible
// but sub-physical over the remainder.  Reported upstream to SeaFreeze.
// The box rounds that outward to round numbers.
//
// No property magnitude from inside is asserted anywhere: those values
// are not water, and pinning them would bless numbers we have told the
// authors are wrong.  What is pinned is WHERE the backend refuses.
constexpr double kExcludedPminMPa = 1500.0;
constexpr double kExcludedTmaxK = 255.0;
bool in_excluded_box(double p_MPa, double T) {
    return p_MPa >= kExcludedPminMPa && T <= kExcludedTmaxK;
}

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
            CHECK_THAT(AS->cpmass(), Catch::Matchers::WithinRel(ref->cpmass(), 2e-4));
            CHECK_THAT(AS->cvmass(), Catch::Matchers::WithinRel(ref->cvmass(), 2e-4));
            CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(ref->speed_sound(), 5e-4));
            // Stability, independent of any reference: these must hold for
            // any physically sensible liquid.
            CHECK(AS->cvmass() < AS->cpmass());
            CHECK(AS->cvmass() > 0.0);
            // Nothing else constrains molar mass: every molar accessor is
            // only ever compared against another instance of this backend,
            // so a wrong value propagates to all of them undetected.
            CHECK_THAT(AS->molar_mass(), Catch::Matchers::WithinRel(ref->molar_mass(), 1e-6));
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

TEST_CASE("Bollengier backend reproduces the authors' published tables", "[Bollengier][water]") {
    // THE primary validation for this backend.  Supplementary Material E of
    // the paper tabulates density, specific heat and sound speed on a
    // 6 x 17 (T, p) grid; these are the authors' own numbers, so unlike a
    // comparison against IAPWS-95 they remain meaningful above 100 MPa,
    // where the two models diverge by design.
    //
    // It also validates the whole chain at once: the coefficient
    // extraction, the COO-45 spline evaluator, and every MPa/Pa unit factor
    // in the property map -- none of which any other test constrains
    // end to end against an external source.
    auto AS = make();
    // p [MPa], T [K], rho [kg/m3], cp [J/kg/K], w [m/s] -- Bollengier,
    // Brown & Shaw (2019) Supplementary Material E, verbatim.  Note the
    // authors do NOT tabulate 250 K above 900 MPa -- the cold
    // high-pressure corner, where the fit is unconstrained and they
    // publish nothing.  That omission is reproduced here rather than
    // filled in, and is the clearest statement available of where the
    // representation stops being supported by data.
    const double kPublished[][5] = {
      {0.1, 250, 991.24, 4554.71, 1248.32},   {0.1, 300, 996.56, 4180.6, 1501.7},     {0.1, 350, 973.73, 4194.57, 1555.01},
      {0.1, 400, 937.41, 4256.01, 1508.94},   {0.1, 450, 889.79, 4396.46, 1396.59},   {0.1, 500, 828.9, 4686.78, 1226.63},
      {100, 250, 1048.72, 3811.19, 1440.87},  {100, 300, 1037.18, 3979.18, 1668.72},  {100, 350, 1013.57, 4025.44, 1733.91},
      {100, 400, 981.8, 4057.5, 1718.23},     {100, 450, 943.51, 4110.6, 1653.52},    {100, 500, 899.25, 4200.91, 1555.51},
      {200, 250, 1090.49, 3673.58, 1668.48},  {200, 300, 1071.06, 3883.77, 1827.92},  {200, 350, 1046.58, 3921.41, 1887.96},
      {200, 400, 1017.06, 3937.7, 1884.44},   {200, 450, 983.27, 3961.82, 1840.58},   {200, 500, 945.83, 4002.64, 1770.19},
      {300, 250, 1123.55, 3614.56, 1858.55},  {300, 300, 1100.12, 3841.68, 1974.03},  {300, 350, 1075.01, 3855.51, 2025.21},
      {300, 400, 1046.83, 3857.77, 2026.35},  {300, 450, 1015.76, 3867.18, 1993.34},  {300, 500, 982.26, 3885.6, 1937.66},
      {400, 250, 1151.38, 3546.17, 2020.41},  {400, 300, 1125.66, 3826.61, 2107.27},  {400, 350, 1100.12, 3812.09, 2149.45},
      {400, 400, 1072.85, 3801.93, 2151.95},  {400, 450, 1043.64, 3801.14, 2125.18},  {400, 500, 1012.75, 3806.9, 2078.01},
      {500, 250, 1175.61, 3467.17, 2161.04},  {500, 300, 1148.53, 3824.1, 2229.07},   {500, 350, 1122.7, 3781.99, 2263.44},
      {500, 400, 1096.11, 3760.85, 2265.67},  {500, 450, 1068.27, 3752.14, 2242.52},  {500, 500, 1039.28, 3750.48, 2200.71},
      {600, 250, 1197.21, 3392.18, 2286.71},  {600, 300, 1169.32, 3829.37, 2341.19},  {600, 350, 1143.28, 3760.18, 2368.62},
      {600, 400, 1117.24, 3729.84, 2369.7},   {600, 450, 1090.44, 3715.25, 2348.89},  {600, 500, 1062.9, 3708.25, 2310.96},
      {700, 250, 1216.79, 3332.28, 2397.11},  {700, 300, 1188.44, 3843.23, 2444.15},  {700, 350, 1162.24, 3748.07, 2465.56},
      {700, 400, 1136.64, 3706.8, 2466.32},   {700, 450, 1110.69, 3686.33, 2446.84},  {700, 500, 1084.3, 3675.89, 2412.0},
      {800, 250, 1234.76, 3291.7, 2501.71},   {800, 300, 1206.25, 3864.86, 2533.84},  {800, 350, 1179.91, 3743.18, 2553.0},
      {800, 400, 1154.64, 3689.13, 2555.86},  {800, 450, 1129.37, 3663.22, 2537.36},  {800, 500, 1103.94, 3651.53, 2504.11},
      {900, 250, 1251.34, 3275.98, 2601.76},  {900, 300, 1222.93, 3884.72, 2622.64},  {900, 350, 1196.46, 3740.83, 2638.34},
      {900, 400, 1171.46, 3675.45, 2641.68},  {900, 450, 1146.76, 3646.03, 2623.19},  {900, 500, 1122.14, 3632.73, 2590.52},
      {1000, 300, 1238.62, 3916.43, 2708.47}, {1000, 350, 1212.04, 3740.25, 2721.16}, {1000, 400, 1187.25, 3666.64, 2723.43},
      {1000, 450, 1163.04, 3634.64, 2704.45}, {1000, 500, 1139.12, 3615.72, 2672.36}, {1200, 300, 1267.64, 3997.19, 2862.02},
      {1200, 350, 1240.74, 3753.43, 2877.3},  {1200, 400, 1216.29, 3661.98, 2876.21}, {1200, 450, 1192.86, 3624.43, 2853.88},
      {1200, 500, 1170.1, 3590.14, 2823.55},  {1400, 300, 1293.8, 4018.89, 3036.22},  {1400, 350, 1266.75, 3783.45, 3007.03},
      {1400, 400, 1242.57, 3665.33, 3020.89}, {1400, 450, 1219.78, 3625.35, 2988.51}, {1400, 500, 1197.92, 3573.43, 2961.77},
      {1600, 300, 1317.73, 4059.3, 3142.26},  {1600, 350, 1290.67, 3815.78, 3152.47}, {1600, 400, 1266.61, 3680.85, 3144.39},
      {1600, 450, 1244.41, 3628.92, 3112.33}, {1600, 500, 1223.16, 3561.25, 3094.27}, {1800, 300, 1339.93, 4104.94, 3279.77},
      {1800, 350, 1312.85, 3850.81, 3267.93}, {1800, 400, 1288.9, 3704.34, 3258.19},  {1800, 450, 1267.1, 3637.35, 3234.67},
      {1800, 500, 1246.58, 3542.89, 3195.29}, {2000, 300, 1360.53, 4147.56, 3398.87}, {2000, 350, 1333.6, 3885.4, 3379.59},
      {2000, 400, 1309.75, 3731.86, 3366.34}, {2000, 450, 1288.21, 3654.84, 3342.31}, {2000, 500, 1268.32, 3536.46, 3301.84},
      {2200, 300, 1379.82, 4175.26, 3512.72}, {2200, 350, 1353.1, 3918.08, 3483.21},  {2200, 400, 1329.36, 3761.33, 3466.69},
      {2200, 450, 1308.02, 3679.22, 3443.02}, {2200, 500, 1288.6, 3539.52, 3402.27},
    };

    for (const auto& r : kPublished) {
        const double p_MPa = r[0], T = r[1];
        INFO("p = " << p_MPa << " MPa, T = " << T << " K");
        AS->update(PT_INPUTS, p_MPa * 1e6, T);
        // 2e-5 is set by the tables' own precision.  Every value is printed
        // to two decimals, so the rounding half-width is LARGEST in relative
        // terms at the smallest entry -- 6.0e-6 at 828.90, falling to 1.1e-6
        // at 4686.78.  2e-5 is ~3.3x that worst case.  Measured agreement is
        // 5.2e-6 (rho), 1.4e-6 (cp), 3.2e-6 (w), i.e. at round-off
        // throughout.
        CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(r[2], 2e-5));
        CHECK_THAT(AS->cpmass(), Catch::Matchers::WithinRel(r[3], 2e-5));
        CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(r[4], 2e-5));
    }
}

TEST_CASE("Bollengier backend reaches the domain IAPWS-95 refuses", "[Bollengier][water]") {
    // The capability this backend exists for.  Values are cross-checked
    // against the published tables above; here the point is that IAPWS-95
    // cannot answer at all.
    auto AS = make();
    auto ref = iapws();
    const double pts[][3] = {
      {250.0, 200.0, 1090.49},
      {300.0, 2000.0, 1360.53},
      {400.0, 2200.0, 1329.36},
      {500.0, 2200.0, 1288.60},
    };
    for (const auto& q : pts) {
        INFO("T = " << q[0] << " K, p = " << q[1] << " MPa");
        AS->update(PT_INPUTS, q[1] * 1e6, q[0]);
        REQUIRE(std::isfinite(AS->rhomass()));
        CHECK_THAT(AS->rhomass(), Catch::Matchers::WithinRel(q[2], 2e-5));
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
    // Warm, so the second state is outside the excluded box.  Any two
    // distinct states prove the point; this one differs from the first in
    // both p and T, which is what the cache has to notice.
    reused->update(PT_INPUTS, 2000e6, 400.0);

    auto fresh = make();
    fresh->update(PT_INPUTS, 2000e6, 400.0);

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
        CHECK_THAT(AS->pmax(), Catch::Matchers::WithinRel(kPmaxMPa * 1e6, 1e-9));
        CHECK_THAT(AS->Tmin(), Catch::Matchers::WithinRel(kTminK, 1e-9));
        CHECK_THAT(AS->Tmax(), Catch::Matchers::WithinRel(kTmaxK, 1e-9));
    }
    SECTION("the backend guards the domain itself, not the spline underneath") {
        // The fitted knots and the paper's literals differ at ULP level
        // (2300.5999999999995 vs 2300.6).  If the guard compared against
        // the literal, p = 2300.6 MPa would pass it and then throw from
        // INSIDE TensorBSpline2D.
        //
        // Asserting only CHECK_THROWS_AS(..., ValueError) cannot tell those
        // apart -- the spline throws the same type -- so this pins the
        // message SOURCE.  Verified by mutation: reverting the guard to the
        // paper literals left all other assertions in this file passing.
        CHECK_THROWS_WITH(AS->update(PT_INPUTS, kPmaxMPa * 1e6, 400.0), Catch::Matchers::ContainsSubstring("BollengierBackend"));
        // ...and the backend must accept its OWN advertised limit, which is
        // the other half: a guard that simply refused everything near the
        // top would also produce a BollengierBackend message.
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, AS->pmax(), 400.0));
        CHECK(std::isfinite(AS->rhomass()));
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, 1e5, AS->Tmin()));
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, 1e5, AS->Tmax()));
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

TEST_CASE("Bollengier backend refuses the excluded box, and nothing else", "[Bollengier][water][nan]") {
    // Contract, matching the reference implementation (SeaFreeze, by a
    // co-author): evaluate the surface wherever it is defined, and refuse
    // only where it produces nothing usable.  SeaFreeze returns NaN there;
    // CoolProp throws instead, because a silently propagating non-finite is
    // the failure this codebase refuses to ship.
    auto AS = make();

    SECTION("inside the box it throws, and says so") {
        // Interior points, far from any domain bound, so the range guard
        // cannot be what fires.  The message is matched on "excluded box"
        // rather than just the type: the post-evaluation admissibility
        // backstop throws the same ValueError, and if the box stopped
        // covering these the backstop would catch some of them and the
        // test would still pass on type alone.
        for (const auto pT : {std::make_pair(1850.0, 240.0), std::make_pair(2290.0, 240.0), std::make_pair(1700.0, 240.0),
                              std::make_pair(2185.0, 240.0), std::make_pair(1600.0, 250.0)}) {
            INFO("p = " << pT.first << " MPa, T = " << pT.second << " K");
            CHECK_THROWS_WITH(AS->update(PT_INPUTS, pT.first * 1e6, pT.second), Catch::Matchers::ContainsSubstring("excluded box"));
        }
    }
    SECTION("the box edge is exactly where it is advertised, inclusive") {
        // The corner itself is refused, and the two states a hair outside
        // it are served.  This is what pins the comparisons: flipping
        // either >= to > or <= to < serves the corner and fails here, and
        // the raw surface is perfectly healthy at all three (cv ~ 2010,
        // w ~ 3065), so nothing but the box can be deciding.
        CHECK_THROWS_WITH(AS->update(PT_INPUTS, kExcludedPminMPa * 1e6, kExcludedTmaxK), Catch::Matchers::ContainsSubstring("excluded box"));
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, kExcludedPminMPa * 1e6, kExcludedTmaxK + 0.1));
        CHECK(std::isfinite(AS->speed_sound()));
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, (kExcludedPminMPa - 0.1) * 1e6, kExcludedTmaxK));
        CHECK(std::isfinite(AS->speed_sound()));
        // ...and the box must not have swallowed the cold low-pressure
        // edge, which is ordinary supercooled liquid.
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, 1e5, kTminK));
    }
    SECTION("no ceiling on the sound speed outside the box") {
        // The box removed the states the old no-ceiling policy was written
        // around, but not the policy: w still reaches 6845 m/s at a state
        // the backend must serve, against roughly 3500 for real water.
        // The populations overlap in w -- pathological states run down to
        // 2497 m/s, below the fastest value the authors publish -- so no
        // threshold on w separates them, and adding one would start
        // refusing served states here.  Pinned by behaviour: a reinstated
        // ceiling anywhere below 6845 fails this NOTHROW.
        auto fast = make();
        REQUIRE_NOTHROW(fast->update(PT_INPUTS, 2300.0e6, 256.0));
        CHECK(fast->speed_sound() > 6000.0);
        // And the comparison that makes it unfixable: a physical state at
        // 2200 MPa / 300 K is SLOWER than that, so the fast one cannot be
        // rejected on speed without rejecting this too.
        auto warm = make();
        warm->update(PT_INPUTS, 2200.0e6, 300.0);
        CHECK(warm->cvmass() > 3000.0);
        CHECK(warm->speed_sound() < fast->speed_sound());
    }
    SECTION("the same pressures are ordinary when warm") {
        for (const double p_MPa : {1900.0, 2175.0, 2290.0}) {
            INFO("p = " << p_MPa << " MPa, T = 400 K");
            REQUIRE_NOTHROW(AS->update(PT_INPUTS, p_MPa * 1e6, 400.0));
            CHECK(std::isfinite(AS->speed_sound()));
            CHECK(AS->cvmass() > 2000.0);
        }
    }
    SECTION("every accepted state is internally consistent") {
        // Self-consistency, which is what the backend actually guarantees.
        // Deliberately NOT a physical-plausibility sweep: after the above,
        // asserting a cv floor here would contradict the documented policy.
        int accepted = 0, expected = 0;
        for (int i = 0; i <= 60; ++i) {
            for (int j = 0; j <= 40; ++j) {
                const double p_MPa = 2300.0 * i / 60.0;
                const double T = kTminK + (kTmaxK - kTminK) * j / 40.0;
                if (in_excluded_box(p_MPa, T)) {
                    continue;  // see kExcludedPminMPa
                }
                ++expected;
                try {
                    AS->update(PT_INPUTS, p_MPa * 1e6, T);
                } catch (const CoolProp::ValueError&) {
                    continue;
                }
                ++accepted;
                INFO("p = " << p_MPa << " MPa, T = " << T << " K");
                REQUIRE(std::isfinite(AS->speed_sound()));
                REQUIRE(std::isfinite(AS->cvmass()));
                CHECK(AS->cvmass() > 0.0);
                CHECK(AS->cvmass() < AS->cpmass());
                CHECK(AS->rhomass() > 700.0);
                CHECK(AS->rhomass() < 1600.0);
            }
        }
        // Binds tightly now: outside the unstable corner the guard must
        // accept EVERY state, so a single spurious refusal fails here.
        REQUIRE(accepted == expected);
    }
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
            // The v branch is not truncation-limited the way the s branch
            // is: measured residuals are 3e-11..2e-10, so 5e-5 there would
            // be ~250,000x looser than the error it is nominally set by.
            // (A 1e-7 unit-factor error would still be caught by the pinned
            // densities; this is about the tolerance meaning what it says.)
            CHECK_THAT(v_fd, Catch::Matchers::WithinRel(1.0 / AS->rhomass(), 1e-8));
        }
    }
}
#endif
