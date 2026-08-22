#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>
#    include <catch2/matchers/catch_matchers_floating_point.hpp>

#    include <algorithm>
#    include <cmath>
#    include <filesystem>
#    include <iterator>
#    include <limits>
#    include <map>
#    include <memory>
#    include <stdexcept>
#    include <string>
#    include <system_error>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/CoolProp.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/Exceptions.h"

using namespace CoolProp;

TEST_CASE("GERG backend families are registered", "[GERG]") {
    backend_families f1 = INVALID_BACKEND_FAMILY, f2 = INVALID_BACKEND_FAMILY;

    extract_backend_families("GERG2004", f1, f2);
    CHECK(f1 == GERG2004_BACKEND_FAMILY);

    extract_backend_families("GERG2008", f1, f2);
    CHECK(f1 == GERG2008_BACKEND_FAMILY);
}

#    include "../Backends/GERG/GERGBackend.h"
#    include "../Backends/GERG/GERGData.h"
#    include "../Backends/GERG/GERGReferenceValues.h"

using namespace CoolProp::GERG;

TEST_CASE("GERG component sets have the published sizes", "[GERG]") {
    CHECK(component_names(GERGModel::GERG_2004).size() == 18);
    CHECK(component_names(GERGModel::GERG_2008).size() == 21);
}

TEST_CASE("GERG pure info matches the monograph", "[GERG]") {
    // Table A3.5, GERG-2004 monograph. Tabulated in mol/dm^3 and kg/kmol;
    // our accessor returns mol/m^3 and kg/mol.
    auto methane = get_pure_info(GERGModel::GERG_2004, "methane");
    CHECK_THAT(methane.rhoc_molm3, Catch::Matchers::WithinRel(10.139342719e3, 1e-14));
    CHECK_THAT(methane.Tc_K, Catch::Matchers::WithinRel(190.564, 1e-14));
    CHECK_THAT(methane.M_kgmol, Catch::Matchers::WithinRel(16.042460e-3, 1e-14));

    // GERG-2008 changes carbon monoxide and isopentane relative to GERG-2004.
    CHECK_THAT(get_pure_info(GERGModel::GERG_2004, "carbonmonoxide").Tc_K, Catch::Matchers::WithinRel(132.800, 1e-14));
    CHECK_THAT(get_pure_info(GERGModel::GERG_2008, "carbonmonoxide").Tc_K, Catch::Matchers::WithinRel(132.860, 1e-14));

    // Added in GERG-2008 only.
    CHECK_THROWS_AS(get_pure_info(GERGModel::GERG_2004, "n-decane"), CoolProp::ValueError);
    CHECK_NOTHROW(get_pure_info(GERGModel::GERG_2008, "n-decane"));
}

TEST_CASE("GERG component resolution accepts CoolProp names and aliases", "[GERG]") {
    CHECK(resolve_component(GERGModel::GERG_2008, "Methane") == "methane");
    CHECK(resolve_component(GERGModel::GERG_2008, "METHANE") == "methane");
    CHECK(resolve_component(GERGModel::GERG_2008, "CO2") == "carbondioxide");
    CHECK(resolve_component(GERGModel::GERG_2008, "124-38-9") == "carbondioxide");
    CHECK(resolve_component(GERGModel::GERG_2008, "n-Butane") == "n-butane");

    // Known to CoolProp, outside the GERG set.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "R134a"), CoolProp::ValueError);
    // In GERG-2008 but not GERG-2004.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2004, "HydrogenSulfide"), CoolProp::ValueError);
    CHECK_NOTHROW(resolve_component(GERGModel::GERG_2008, "HydrogenSulfide"));
    // Not a fluid at all.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "NOT A FLUID"), CoolProp::ValueError);
}

TEST_CASE("GERG pure residual coefficients are internally consistent", "[GERG]") {
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            auto pc = get_pure_coeffs(model, name);
            REQUIRE(pc.n.size() > 0);
            CHECK(pc.t.size() == pc.n.size());
            CHECK(pc.d.size() == pc.n.size());
            CHECK(pc.c.size() == pc.n.size());
            CHECK(pc.l.size() == pc.n.size());
            // c_i is 1 exactly when l_i > 0
            for (std::size_t i = 0; i < pc.n.size(); ++i) {
                CHECK(pc.c[i] == ((pc.l[i] > 0) ? 1.0 : 0.0));
            }
        }
    }
}

TEST_CASE("GERG pure residual coefficients match the monograph", "[GERG]") {
    // Methane, Table A3.2 (24 terms), first and last n.
    auto ch4 = get_pure_coeffs(GERGModel::GERG_2004, "methane");
    REQUIRE(ch4.n.size() == 24);
    CHECK_THAT(ch4.n[0], Catch::Matchers::WithinRel(0.57335704239162, 1e-14));

    // The generalised 12-term set shared by most fluids.
    auto c3h8 = get_pure_coeffs(GERGModel::GERG_2004, "propane");
    REQUIRE(c3h8.n.size() == 12);
    CHECK(c3h8.t == std::vector<double>{0.250, 1.125, 1.500, 1.375, 0.250, 0.875, 0.625, 1.750, 3.625, 3.625, 14.500, 12.000});
    CHECK(c3h8.d == std::vector<double>{1, 1, 1, 2, 3, 7, 2, 5, 1, 4, 3, 4});
    CHECK(c3h8.l == std::vector<double>{0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3, 3});

    // GERG-2008 changed carbon monoxide and isopentane.
    CHECK(get_pure_coeffs(GERGModel::GERG_2004, "carbonmonoxide").n != get_pure_coeffs(GERGModel::GERG_2008, "carbonmonoxide").n);
    CHECK_THAT(get_pure_coeffs(GERGModel::GERG_2008, "carbonmonoxide").n[0], Catch::Matchers::WithinRel(0.90554, 1e-14));

    // Unchanged fluids fall through identically.
    CHECK(get_pure_coeffs(GERGModel::GERG_2004, "methane").n == get_pure_coeffs(GERGModel::GERG_2008, "methane").n);
}

namespace {

/// The reduced ideal-gas Helmholtz energy of one pure component, written out
/// exactly as teqp's GERG200XAlphaig::alphaig_pure does (GERG.hpp:395-410).
/// Note that the R*/R ratio multiplies the WHOLE bracket, and that the two
/// cosh terms are SUBTRACTED while the coefficients themselves are stored
/// positive.
double gerg_alphaig_pure(const AlphaigCoeffs& c, const PureInfo& info, double T, double rho) {
    const double x = info.Tc_K / T;
    double s = c.n0[1] + c.n0[2] * x + c.n0[3] * std::log(x);
    if (c.theta0[4] != 0) {
        s += c.n0[4] * std::log(std::abs(std::sinh(c.theta0[4] * x)));
    }
    if (c.theta0[6] != 0) {
        s += c.n0[6] * std::log(std::abs(std::sinh(c.theta0[6] * x)));
    }
    if (c.theta0[5] != 0) {
        s -= c.n0[5] * std::log(std::abs(std::cosh(c.theta0[5] * x)));
    }
    if (c.theta0[7] != 0) {
        s -= c.n0[7] * std::log(std::abs(std::cosh(c.theta0[7] * x)));
    }
    return std::log(rho / info.rhoc_molm3) + (RSTAR_GERG / R_GERG) * s;
}

}  // namespace

TEST_CASE("GERG ideal-gas coefficient tables have the padded monograph shape", "[GERG]") {
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            auto c = get_alphaig_coeffs(model, name);
            REQUIRE(c.n0.size() == 8);
            REQUIRE(c.theta0.size() == 8);
            // Padding so that the monograph's 1-based indices work.
            CHECK(c.n0[0] == 0.0);
            for (std::size_t i = 0; i < 4; ++i) {
                CHECK(c.theta0[i] == 0.0);
            }
        }
    }
    // Added in GERG-2008 only.
    CHECK_THROWS_AS(get_alphaig_coeffs(GERGModel::GERG_2004, "n-decane"), CoolProp::ValueError);
    CHECK_NOTHROW(get_alphaig_coeffs(GERGModel::GERG_2008, "n-decane"));
}

TEST_CASE("GERG ideal-gas theta values match the monograph", "[GERG]") {
    // Table A3.1: methane. n0[1] and n0[2] are recomputed, so only n0[3..7]
    // and theta0[4..7] are comparable with the published table.
    auto c = get_alphaig_coeffs(GERGModel::GERG_2004, "methane");
    CHECK_THAT(c.n0[3], Catch::Matchers::WithinRel(3.000880, 1e-14));
    CHECK_THAT(c.n0[4], Catch::Matchers::WithinRel(0.763150, 1e-14));
    CHECK_THAT(c.n0[7], Catch::Matchers::WithinRel(-4.469210000, 1e-14));
    CHECK_THAT(c.theta0[4], Catch::Matchers::WithinRel(4.306474465, 1e-14));
    CHECK_THAT(c.theta0[7], Catch::Matchers::WithinRel(5.722644361, 1e-14));

    // Sign convention: the two cosh coefficients are stored POSITIVE as
    // published, with the minus sign living in the evaluating expression.
    // Task 8 hands these straight to IdealHelmholtzGERG2004Cosh, whose all()
    // accumulates -n[i]*log(|cosh(...)|) internally.
    auto h2o = get_alphaig_coeffs(GERGModel::GERG_2004, "water");
    CHECK(h2o.n0[5] > 0.0);
    CHECK_THAT(h2o.n0[5], Catch::Matchers::WithinRel(0.987630, 1e-14));

    // GERG-2008 changed carbon monoxide and isopentane and added n-decane.
    CHECK(get_alphaig_coeffs(GERGModel::GERG_2004, "carbonmonoxide").theta0 != get_alphaig_coeffs(GERGModel::GERG_2008, "carbonmonoxide").theta0);
    CHECK_THAT(get_alphaig_coeffs(GERGModel::GERG_2008, "carbonmonoxide").theta0[4], Catch::Matchers::WithinRel(11.669802800, 1e-14));
    // Unchanged fluids fall through identically.
    CHECK(get_alphaig_coeffs(GERGModel::GERG_2004, "methane").n0 == get_alphaig_coeffs(GERGModel::GERG_2008, "methane").n0);
}

TEST_CASE("GERG ideal-gas integration constants zero h and s at the reference state", "[GERG]") {
    // teqp GERG.hpp:376-379 -- h = s = 0 for the IDEAL GAS at 298.15 K and
    // 101325 Pa, with rho0 = p0/(R*T0) (R, not R*).
    const double T0 = 298.15, p0 = 101325.0;
    const double rho0 = p0 / (R_GERG * T0);

    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            auto info = get_pure_info(model, name);
            auto c = get_alphaig_coeffs(model, name);

            // Aig10 = tau*d(alphaig)/d(tau) = -T*d(alphaig)/dT, which is
            // independent of the reducing temperature used to form tau.
            const double dT = 1e-5 * T0;
            double dalpha_dT = (gerg_alphaig_pure(c, info, T0 + dT, rho0) - gerg_alphaig_pure(c, info, T0 - dT, rho0)) / (2 * dT);
            double Aig10 = -T0 * dalpha_dT;

            double h_over_RT = 1 + Aig10;                                    // ideal gas: h/(RT) = 1 + Aig10
            double s_over_R = Aig10 - gerg_alphaig_pure(c, info, T0, rho0);  // s/R = Aig10 - Aig00

            CHECK_THAT(h_over_RT, Catch::Matchers::WithinAbs(0.0, 1e-8));
            CHECK_THAT(s_over_R, Catch::Matchers::WithinAbs(0.0, 1e-8));
        }
    }
}

TEST_CASE("GERG recomputed integration constants match an independent solve", "[GERG]") {
    // THIS IS THE TRANSCRIPTION GUARD FOR THE WHOLE IDEAL-GAS TABLE.
    //
    // {n0[1], n0[2]} is a complete fingerprint of a component's ideal-gas
    // row: n0[3..7], theta0[4..7], Tc and rhoc all feed the 2x2 solve, so a
    // single corrupted digit anywhere in the row moves these two numbers far
    // beyond the 1e-12 tolerance below.  That matters because the h = s = 0
    // test cannot catch a transcription error at all -- corrupting, say,
    // propane's theta0[6] leaves the solver and the evaluator mutually
    // consistent and h and s still vanish.  Every (model, component) pair is
    // therefore listed here: 18 for GERG-2004 and 21 for GERG-2008.
    //
    // The literals come from a completely independent implementation -- a
    // NumPy script that rebuilds teqp's 2x2 system from scratch, solves it
    // with numpy.linalg.solve (LU with partial pivoting, not Cramer's rule)
    // and verifies h and s vanish using ANALYTIC derivatives.  They were NOT
    // captured from this backend's own output, so they assert what the
    // coefficients should be rather than what the code currently does.  See
    // task-4-report.md for the script and its transcript.
    //
    // This test case is also what pins the R*/R convention, and it is the
    // ONLY thing that does.  The h = s = 0 test cannot: moving the ratio from
    // outside the bracket (GERG-2008 / teqp) to inside it (GERG-2004) simply
    // rescales n0[1] and n0[2] by exactly R*/R, leaving alpha^0 bit-identical,
    // so h and s still vanish.  What these literals catch is (a) that
    // convention swap, which changes them by 4.6e-6 relative, and (b) an
    // outright coding error that drops the ratio, which changes them by
    // 1.8e-8 to 3.4e-6 relative.  Both are enormous against a 1e-12
    // tolerance.
    struct Expect
    {
        GERGModel model;
        const char* name;
        double n1, n2;
    };
    const Expect cases[] = {
      {GERGModel::GERG_2004, "methane", 19.597508817430203, -83.95966789022567},
      {GERGModel::GERG_2004, "nitrogen", 11.083407489057498, -22.202102427578456},
      {GERGModel::GERG_2004, "carbondioxide", 11.925152757587837, -16.118762264778105},
      {GERGModel::GERG_2004, "ethane", 24.67543752664495, -77.42531376123445},
      {GERGModel::GERG_2004, "propane", 31.602908195059367, -84.46328438999367},
      {GERGModel::GERG_2004, "n-butane", 20.88414336067719, -91.63847802844411},
      {GERGModel::GERG_2004, "isobutane", 20.413726078976026, -94.46762003594714},
      {GERGModel::GERG_2004, "n-pentane", 28.587336515506692, -96.26533664852582},
      {GERGModel::GERG_2004, "isopentane", 29.158567601789613, -111.21604889349067},
      {GERGModel::GERG_2004, "n-hexane", 32.49945909536515, -103.869150116756},
      {GERGModel::GERG_2004, "n-heptane", 37.237679270562055, -105.72419451985124},
      {GERGModel::GERG_2004, "n-octane", 42.1431834637039, -106.34926315689664},
      {GERGModel::GERG_2004, "hydrogen", 13.79644339318705, -175.86448729341214},
      {GERGModel::GERG_2004, "oxygen", 10.001843585817623, -14.996095135028263},
      {GERGModel::GERG_2004, "carbonmonoxide", 10.814470255578003, -19.843695434544927},
      {GERGModel::GERG_2004, "water", 8.216535516334572, -12.002441239189446},
      {GERGModel::GERG_2004, "helium", 13.628409737317194, -143.47075960157534},
      {GERGModel::GERG_2004, "argon", 8.316631499886697, -4.946502600476912},

      {GERGModel::GERG_2008, "methane", 19.597508817430203, -83.95966789022567},
      {GERGModel::GERG_2008, "nitrogen", 11.083407489057498, -22.202102427578456},
      {GERGModel::GERG_2008, "carbondioxide", 11.925152757587837, -16.118762264778105},
      {GERGModel::GERG_2008, "ethane", 24.67543752664495, -77.42531376123445},
      {GERGModel::GERG_2008, "propane", 31.602908195059367, -84.46328438999367},
      {GERGModel::GERG_2008, "n-butane", 20.88414336067719, -91.63847802844411},
      {GERGModel::GERG_2008, "isobutane", 20.413726078976026, -94.46762003594714},
      {GERGModel::GERG_2008, "n-pentane", 28.587336515506692, -96.26533664852582},
      {GERGModel::GERG_2008, "isopentane", 29.15856192130587, -111.21604889349067},
      {GERGModel::GERG_2008, "n-hexane", 32.49945909536515, -103.869150116756},
      {GERGModel::GERG_2008, "n-heptane", 37.237679270562055, -105.72419451985124},
      {GERGModel::GERG_2008, "n-octane", 42.1431834637039, -106.34926315689664},
      {GERGModel::GERG_2008, "hydrogen", 13.79644339318705, -175.86448729341214},
      {GERGModel::GERG_2008, "oxygen", 10.001843585817623, -14.996095135028263},
      {GERGModel::GERG_2008, "carbonmonoxide", 10.813340744153283, -19.834733958634743},
      {GERGModel::GERG_2008, "water", 8.216535516334572, -12.002441239189446},
      {GERGModel::GERG_2008, "helium", 13.628409737317194, -143.47075960157534},
      {GERGModel::GERG_2008, "argon", 8.316631499886697, -4.946502600476912},
      {GERGModel::GERG_2008, "hydrogensulfide", 9.33619774177303, -16.266508993594602},
      {GERGModel::GERG_2008, "n-nonane", 46.72362520349981, -112.01770583722839},
      {GERGModel::GERG_2008, "n-decane", 50.35302335379572, -120.01206647981711},
    };
    // Guard against a component being dropped from the list above.
    CHECK(std::size(cases) == component_names(GERGModel::GERG_2004).size() + component_names(GERGModel::GERG_2008).size());

    for (const auto& e : cases) {
        CAPTURE(e.name);
        auto c = get_alphaig_coeffs(e.model, e.name);
        CHECK_THAT(c.n0[1], Catch::Matchers::WithinRel(e.n1, 1e-12));
        CHECK_THAT(c.n0[2], Catch::Matchers::WithinRel(e.n2, 1e-12));
    }

    // The published integration constants really are discarded: for the
    // heavier alkanes the monograph's own n0[1] refers to a different
    // reference state and differs by ~26 (not by roundoff).
    CHECK(std::abs(get_alphaig_coeffs(GERGModel::GERG_2004, "n-octane").n0[1] - 15.864709639) > 1.0);
}

TEST_CASE("GERG reducing parameters exist for every binary pair", "[GERG]") {
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const auto& names = component_names(model);
        for (std::size_t i = 0; i < names.size(); ++i) {
            for (std::size_t j = i + 1; j < names.size(); ++j) {
                CAPTURE(names[i], names[j]);
                CHECK_NOTHROW(get_betasgammas(model, names[i], names[j]));
                CHECK_NOTHROW(get_betasgammas(model, names[j], names[i]));
            }
        }
    }
}

TEST_CASE("GERG reducing parameters invert correctly when the pair is reversed", "[GERG]") {
    auto fwd = get_betasgammas(GERGModel::GERG_2008, "methane", "nitrogen");
    auto rev = get_betasgammas(GERGModel::GERG_2008, "nitrogen", "methane");
    CHECK_THAT(rev.betaT, Catch::Matchers::WithinRel(1.0 / fwd.betaT, 1e-14));
    CHECK_THAT(rev.betaV, Catch::Matchers::WithinRel(1.0 / fwd.betaV, 1e-14));
    // Gammas are symmetric, not reciprocal.
    CHECK_THAT(rev.gammaT, Catch::Matchers::WithinRel(fwd.gammaT, 1e-14));
    CHECK_THAT(rev.gammaV, Catch::Matchers::WithinRel(fwd.gammaV, 1e-14));
}

TEST_CASE("GERG-2008 changes some reducing parameters relative to GERG-2004", "[GERG]") {
    // Table A8 of the GERG-2008 manuscript revises a subset of pairs.
    // At least one pair must differ, and pairs GERG-2008 did not touch
    // must fall through unchanged.
    bool any_different = false;
    const auto& names = component_names(GERGModel::GERG_2004);
    for (std::size_t i = 0; i < names.size(); ++i) {
        for (std::size_t j = i + 1; j < names.size(); ++j) {
            auto a = get_betasgammas(GERGModel::GERG_2004, names[i], names[j]);
            auto b = get_betasgammas(GERGModel::GERG_2008, names[i], names[j]);
            if (a.betaT != b.betaT || a.gammaT != b.gammaT || a.betaV != b.betaV || a.gammaV != b.gammaV) {
                any_different = true;
            }
        }
    }
    CHECK(any_different);
}

TEST_CASE("GERG reducing parameter lookup rejects unknown fluids", "[GERG]") {
    CHECK_THROWS_AS(get_betasgammas(GERGModel::GERG_2004, "NOT A FLUID", "water"), CoolProp::ValueError);
}

TEST_CASE("GERG has departure functions for exactly the published pairs", "[GERG]") {
    const auto& names = component_names(GERGModel::GERG_2008);
    std::size_t with_departure = 0, with_scaled_F = 0;
    for (std::size_t i = 0; i < names.size(); ++i) {
        for (std::size_t j = i + 1; j < names.size(); ++j) {
            double F = 0;
            if (get_Fij(GERGModel::GERG_2008, names[i], names[j], F)) {
                ++with_departure;
                if (F != 1.0) ++with_scaled_F;
                CAPTURE(names[i], names[j]);
                CHECK_NOTHROW(get_departurecoeffs(GERGModel::GERG_2008, names[i], names[j]));
            } else {
                CHECK_THROWS_AS(get_departurecoeffs(GERGModel::GERG_2008, names[i], names[j]), CoolProp::ValueError);
            }
        }
    }
    // Table A3.6: 15 pairs carry a departure function; 7 of them use the
    // generalised departure function scaled by F_ij != 1.
    CHECK(with_departure == 15);
    CHECK(with_scaled_F == 7);
}

TEST_CASE("GERG F_ij is symmetric", "[GERG]") {
    double a = 0, b = 0;
    REQUIRE(get_Fij(GERGModel::GERG_2008, "methane", "isobutane", a));
    REQUIRE(get_Fij(GERGModel::GERG_2008, "isobutane", "methane", b));
    CHECK_THAT(a, Catch::Matchers::WithinRel(b, 1e-14));
    CHECK_THAT(a, Catch::Matchers::WithinRel(0.771035405688, 1e-12));
}

TEST_CASE("GERG departure coefficients are internally consistent", "[GERG]") {
    auto dc = get_departurecoeffs(GERGModel::GERG_2008, "methane", "nitrogen");
    REQUIRE(dc.n.size() > 0);
    CHECK(dc.t.size() == dc.n.size());
    CHECK(dc.d.size() == dc.n.size());
    CHECK(dc.eta.size() == dc.n.size());
    CHECK(dc.beta.size() == dc.n.size());
    CHECK(dc.gamma.size() == dc.n.size());
    CHECK(dc.epsilon.size() == dc.n.size());

    // The polynomial block comes first and is contiguous.
    std::size_t np = departure_Npower(dc);
    for (std::size_t k = 0; k < np; ++k) {
        CHECK(dc.eta[k] == 0.0);
        CHECK(dc.beta[k] == 0.0);
    }
    for (std::size_t k = np; k < dc.n.size(); ++k) {
        CHECK((dc.eta[k] != 0.0 || dc.beta[k] != 0.0));
    }
}

TEST_CASE("GERG-2004 and GERG-2008 share departure data", "[GERG]") {
    CHECK(get_departurecoeffs(GERGModel::GERG_2004, "methane", "nitrogen").n == get_departurecoeffs(GERGModel::GERG_2008, "methane", "nitrogen").n);
}

TEST_CASE("GERG departure_Npower throws on a non-contiguous polynomial block", "[GERG]") {
    // None of the 15 real departure functions violate contiguity (all 15
    // were checked by hand: see task-6-report.md), so this branch is
    // otherwise never exercised anywhere in the suite. Construct one by
    // hand: term 0 is polynomial (eta == beta == 0), term 1 is Gaussian,
    // and term 2 reverts to polynomial -- exactly the layout
    // GERG2008DepartureFunction's constructor cannot represent, since it
    // assumes every polynomial term is a contiguous prefix.
    DepartureCoeffs dc;
    dc.n = {1.0, 2.0, 3.0};
    dc.d = {1.0, 2.0, 3.0};
    dc.t = {1.0, 2.0, 3.0};
    dc.eta = {0.0, 1.0, 0.0};
    dc.beta = {0.0, 1.0, 0.0};
    dc.gamma = {0.0, 0.5, 0.0};
    dc.epsilon = {0.0, 0.5, 0.0};
    CHECK_THROWS_AS(departure_Npower(dc), CoolProp::ValueError);
}

TEST_CASE("GERG departure_Npower does not throw when every term is polynomial", "[GERG]") {
    // The adjacent off-by-one to the case above: Npower == n.size() (no
    // Gaussian block at all) is a real, valid layout -- the generalized
    // departure function and methane/hydrogen's departure function are both
    // like this (see task-6-report.md) -- and must NOT throw.
    DepartureCoeffs dc;
    dc.n = {1.0, 2.0, 3.0};
    dc.d = {1.0, 2.0, 3.0};
    dc.t = {1.0, 2.0, 3.0};
    dc.eta = {0.0, 0.0, 0.0};
    dc.beta = {0.0, 0.0, 0.0};
    dc.gamma = {0.0, 0.0, 0.0};
    dc.epsilon = {0.0, 0.0, 0.0};
    std::size_t np = 0;
    CHECK_NOTHROW(np = departure_Npower(dc));
    CHECK(np == dc.n.size());
}

// ###########################################################################
// Task 8: the ASSEMBLED pure-fluid equation of state against teqp.
//
// Everything above this line compares CoolProp's coefficient TABLES against
// the monograph.  Everything below compares CoolProp's assembled EOS against
// teqp's NUMBERS (src/Backends/GERG/GERGReferenceValues.h), which is the only
// thing that can catch a wiring error -- a swapped exponent vector, a missing
// R*/R, a negated cosh coefficient.
// ###########################################################################

namespace {

/// The reference grids below intentionally sweep a wider (T, rho) range than
/// GERGMixtureBackend::check_gerg_range_of_validity (task 10) enforces on
/// ordinary use: they include per-COMPONENT points above the model's 700 K
/// MIXTURE Tmax (propane at 739.65 K is one), and states deep in the
/// two-phase dome where the single-phase EOS legitimately returns a negative
/// or huge p. These are checks of the raw assembled EOS against teqp's
/// numbers, not claims that every grid point is a physically attainable GERG
/// state, so the sweep needs the same escape hatch CoolProp already uses
/// elsewhere for exactly this tension (Ancillaries.cpp, TabularBackends.cpp):
/// DONT_CHECK_PROPERTY_LIMITS.
class PropertyLimitsSuspender
{
   public:
    PropertyLimitsSuspender() : was_set(CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS)) {
        CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, true);
    }
    ~PropertyLimitsSuspender() {
        CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, was_set);
    }
    PropertyLimitsSuspender(const PropertyLimitsSuspender&) = delete;
    PropertyLimitsSuspender& operator=(const PropertyLimitsSuspender&) = delete;
    // Deleting the copy operations already suppresses the implicit moves, but
    // spelling them out keeps the rule of five complete and satisfies
    // cppcoreguidelines-special-member-functions.
    PropertyLimitsSuspender(PropertyLimitsSuspender&&) = delete;
    PropertyLimitsSuspender& operator=(PropertyLimitsSuspender&&) = delete;

   private:
    bool was_set;
};

/// Per-FIELD comparison tally.  GERGReferenceValues.h nulls fields
/// independently (a row can have finite p and c_v but NaN w), so a NaN must
/// skip that ONE comparison, never the whole row.  Counting both halves and
/// pinning the totals in the test below means a reference field that silently
/// turns into a NaN shows up as a failing count rather than as a quietly
/// weakened gate.
///
/// nan_skipped is broken out PER FIELD (nan_skipped_w etc), not just totalled,
/// because a cross-field total is a hazard by itself: if a NaN ever migrated
/// from `w` on one row to, say, `p` on another (a generator bug, or a field
/// getting swapped), the TOTAL nan_skipped would stay unchanged while a `p`
/// comparison silently stopped happening -- and the test below would stay
/// green throughout. Pinning nan_skipped_p == 0 (etc) separately from
/// nan_skipped_w == <the known count> is what would catch that.
struct RefTally
{
    std::size_t compared = 0;
    std::size_t rows = 0;
    std::size_t nan_skipped_alphar = 0;
    std::size_t nan_skipped_alphaig = 0;
    std::size_t nan_skipped_p = 0;
    std::size_t nan_skipped_cvmolar = 0;
    std::size_t nan_skipped_w = 0;

    /// Resolve a field label to ITS OWN counter, so `check_field` cannot be
    /// called with a label/counter pair that disagree -- task-10 review round
    /// 2: the previous signature took the counter as a separate `size_t&`
    /// parameter, and nothing stopped a call site from passing
    /// `nan_skipped_w` alongside the label `"p"`. All 10 call sites happened
    /// to be correct, but the mismatch was representable. Routing every
    /// increment through this single label-keyed lookup makes it not: there
    /// is exactly one counter reachable from a given label, full stop.
    std::size_t& nan_counter_for(const std::string& label) {
        static const std::map<std::string, std::size_t RefTally::*> table = {
          {"alphar", &RefTally::nan_skipped_alphar}, {"alphaig", &RefTally::nan_skipped_alphaig},
          {"p", &RefTally::nan_skipped_p},           {"cvmolar", &RefTally::nan_skipped_cvmolar},
          {"w", &RefTally::nan_skipped_w},
        };
        auto it = table.find(label);
        if (it == table.end()) {
            throw std::logic_error(std::string("check_field: no NaN counter registered for field [") + label + "]");
        }
        return this->*(it->second);
    }
};

/// Compare one field, or record that the reference value is NaN.
/// `computed` is only evaluated when the reference field is finite -- some of
/// the NaN'd fields (w on the mechanically unstable branch) correspond to
/// states where CoolProp's own evaluation is meaningless too.
///
/// `label` is not decoration: all five comparisons below render identically in
/// Catch2's output (`CHECK_THAT(static_cast<double>(computed()), ...)`), and
/// the surrounding CAPTURE carries only backend/name/T/rho.  Without the label
/// a failure cannot be attributed to a field, which is exactly the information
/// needed to localize a regression (an alphaig-only failure means the
/// ideal-gas constants; a p-only failure means R_u; and so on). `label` is
/// now ALSO what selects the NaN counter (via `tally.nan_counter_for`), so
/// there is only one string to keep in sync, not a string plus a separately-
/// chosen counter reference.
// `computed` is taken by const reference, not as a forwarding reference: it is
// only ever INVOKED here, never stored or moved from, so a forwarding
// reference would promise an ownership transfer that never happens
// (cppcoreguidelines-missing-std-forward).
template <typename F>
void check_field(const char* label, double reference, const F& computed, double rel_tol, RefTally& tally) {
    if (std::isnan(reference)) {
        ++tally.nan_counter_for(label);
        return;
    }
    ++tally.compared;
    INFO("field := " << label);
    CHECK_THAT(static_cast<double>(computed()), Catch::Matchers::WithinRel(reference, rel_tol));
}

void compare_pure_points(const std::string& backend, const std::vector<CoolProp::GERG::reference::PureRefPoint>& points, RefTally& tally) {
    std::shared_ptr<CoolProp::AbstractState> AS;
    std::string current_name;
    for (const auto& pt : points) {
        if (AS == nullptr || current_name != pt.name) {
            current_name = pt.name;
            AS.reset(CoolProp::AbstractState::factory(backend, std::vector<std::string>{current_name}));
            // The reference grid deliberately includes states inside the
            // two-phase dome (T = 0.7*Tc, rho = 0.5*rhoc gives p < 0 for most
            // components).  teqp evaluates the single-phase EOS there
            // analytically; imposing a phase makes CoolProp do the same
            // instead of trying to run a saturation solver the GERG backend
            // has no ancillary equations for.
            AS->specify_phase(CoolProp::iphase_gas);
        }
        CAPTURE(backend, pt.name, pt.T_K, pt.rhomolar);
        AS->update(CoolProp::DmolarT_INPUTS, pt.rhomolar, pt.T_K);
        ++tally.rows;

        // alphar and alphaig are the two halves of the EOS in isolation.
        // alphaig is the ONLY assertion here that can see a wrong ideal-gas
        // integration constant, a missing R*/R or a negated cosh coefficient:
        // p, c_v and w are all blind to those.
        check_field("alphar", pt.alphar, [&] { return AS->alphar(); }, 1e-12, tally);
        check_field("alphaig", pt.alphaig, [&] { return AS->alpha0(); }, 1e-12, tally);
        check_field("p", pt.p_Pa, [&] { return AS->p(); }, 1e-10, tally);
        check_field("cvmolar", pt.cvmolar, [&] { return AS->cvmolar(); }, 1e-10, tally);
        check_field("w", pt.w, [&] { return AS->speed_sound(); }, 1e-10, tally);
    }
}

void compare_mix_points(const std::string& backend, const std::vector<CoolProp::GERG::reference::MixRefPoint>& points, RefTally& tally,
                        std::size_t& rows_trimmed) {
    std::size_t row = 0;
    for (const auto& pt : points) {
        // Not decoration: names and z are parallel vectors, and a row whose
        // composition had been truncated would otherwise be handed to
        // set_mole_fractions and throw somewhere far from here.
        REQUIRE(pt.names.size() == pt.z.size());

        // Drop components with a mole fraction of EXACTLY zero, preserving the
        // order of the ones that remain.
        //
        // This is not a convenience: CoolProp's GERG2008ReducingFunction
        // computes f_Y_ij = x_i x_j (x_i + x_j) / (beta^2 x_i + x_j)
        // (ReducingFunctions.cpp:565-568), which is 0/0 = NaN whenever BOTH
        // mole fractions are zero -- and the AGA8 rows carry all 21 component
        // names with many exact zeros.  The whole reducing state then goes
        // NaN.  Removing a zero-fraction component is exactly equivalent
        // analytically (its x_i^2 Y_c,i term, its x_i x_j cross terms, its
        // x_i alpha^r_i and its x_i (alpha^0_i + ln x_i) are all identically
        // zero), and it preserves the AGA8 column order, so it cannot
        // reintroduce the AGA8-vs-component_names ordering bug the generator
        // warns about -- unlike re-deriving the order from a name list.
        // The count of trimmed rows is returned and pinned by the caller so
        // this cannot silently start applying to rows it should not.
        std::vector<std::string> names;
        std::vector<CoolPropDbl> z;
        for (std::size_t k = 0; k < pt.names.size(); ++k) {
            if (pt.z[k] != 0.0) {
                names.emplace_back(pt.names[k]);
                z.push_back(static_cast<CoolPropDbl>(pt.z[k]));
            }
        }
        REQUIRE(names.size() >= 2);
        if (names.size() != pt.names.size()) {
            ++rows_trimmed;
        }

        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, names));
        AS->set_mole_fractions(z);
        // Same rationale as the pure grid: the fixed (250 K, 5000 mol/m^3)
        // binary-pair state is inside the two-phase dome for many pairs, and
        // the GERG backend has no ancillary equations to run a saturation
        // solver with.  teqp evaluates the single-phase EOS analytically.
        AS->specify_phase(CoolProp::iphase_gas);
        CAPTURE(backend, row, names.size(), names.front(), names.back(), pt.T_K, pt.rhomolar);
        ++row;
        AS->update(CoolProp::DmolarT_INPUTS, pt.rhomolar, pt.T_K);
        ++tally.rows;

        // alphaig is the sharpest instrument here and the reason the column
        // was added to MixRefPoint: the mixture branch of
        // calc_alpha0_deriv_nocache is a completely different code path from
        // the pure branch, and an error in it (e.g. the spurious Tc,i/Tr
        // factor the GERG2004Sinh/Cosh terms would introduce) moves alphaig
        // while leaving alphar and p exactly right.
        check_field("alphar", pt.alphar, [&] { return AS->alphar(); }, 1e-12, tally);
        check_field("alphaig", pt.alphaig, [&] { return AS->alpha0(); }, 1e-12, tally);
        check_field("p", pt.p_Pa, [&] { return AS->p(); }, 1e-10, tally);
        check_field("cvmolar", pt.cvmolar, [&] { return AS->cvmolar(); }, 1e-10, tally);
        check_field("w", pt.w, [&] { return AS->speed_sound(); }, 1e-10, tally);
    }
}

}  // namespace

TEST_CASE("GERG pure fluids reproduce teqp", "[GERG]") {
    using namespace CoolProp::GERG::reference;
    PropertyLimitsSuspender suspend_limits;

    RefTally t2004;
    compare_pure_points("GERG2004", pure_points_2004, t2004);
    // 18 GERG-2004 components x 16 (T, rho) grid points, 5 fields each.
    // 23 of the 288 w values are NaN (mechanically unstable branch); no
    // alphar/alphaig/p/cv value is. Pinned per field, not as a 23-total, so a
    // NaN migrating from w to (say) p on a different row cannot hide behind an
    // unchanged total -- it would show up as nan_skipped_p != 0 here.
    CHECK(t2004.rows == 288);
    CHECK(t2004.nan_skipped_alphar == 0);
    CHECK(t2004.nan_skipped_alphaig == 0);
    CHECK(t2004.nan_skipped_p == 0);
    CHECK(t2004.nan_skipped_cvmolar == 0);
    CHECK(t2004.nan_skipped_w == 23);
    CHECK(t2004.compared == 288 * 5 - 23);

    RefTally t2008;
    compare_pure_points("GERG2008", pure_points_2008, t2008);
    // 21 GERG-2008 components x 16 grid points; 27 NaN w values.
    CHECK(t2008.rows == 336);
    CHECK(t2008.nan_skipped_alphar == 0);
    CHECK(t2008.nan_skipped_alphaig == 0);
    CHECK(t2008.nan_skipped_p == 0);
    CHECK(t2008.nan_skipped_cvmolar == 0);
    CHECK(t2008.nan_skipped_w == 27);
    CHECK(t2008.compared == 336 * 5 - 27);
}

TEST_CASE("GERG pure fluid uses the GERG gas constant and reducing state", "[GERG]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    // R, not R* (8.314510) and not CoolProp's own per-fluid gas constant.
    CHECK_THAT(AS->gas_constant(), Catch::Matchers::WithinRel(8.314472, 1e-14));
    CHECK_THAT(AS->T_reducing(), Catch::Matchers::WithinRel(190.564, 1e-12));
    CHECK_THAT(AS->rhomolar_reducing(), Catch::Matchers::WithinRel(10.139342719e3, 1e-12));
    CHECK_THAT(AS->molar_mass(), Catch::Matchers::WithinRel(16.042460e-3, 1e-12));

    // The GERG tables, not CoolProp's fluid library: CoolProp's own methane
    // EOS reduces at 190.564 K but 10139.128 mol/m^3, and uses R = 8.3144598.
    std::shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS", std::vector<std::string>{"Methane"}));
    CHECK(std::abs(HEOS->rhomolar_reducing() - AS->rhomolar_reducing()) > 1e-3);
}

TEST_CASE("GERG pure fluid reference state gives h = s = 0 for the ideal gas", "[GERG]") {
    // Consistency check on the whole assembled ideal-gas path, not just the
    // coefficient solve tested in Task 4: this goes through
    // calc_alpha0_deriv_nocache, so it also pins the Tc-vs-T_red convention
    // and the R*/R placement on the pure branch.
    const double T0 = 298.15, p0 = 101325.0;
    const double rho0 = p0 / (8.314472 * T0);
    for (const auto& model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const std::string backend = (model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008";
        for (const auto& name : component_names(model)) {
            CAPTURE(backend, name);
            std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{name}));
            AS->specify_phase(iphase_gas);
            AS->update(DmolarT_INPUTS, rho0, T0);
            // h_ig/(R T) = 1 + tau*dalpha0_dtau and s_ig/R = tau*dalpha0_dtau - alpha0.
            const double h_ig_over_RT = 1 + AS->tau() * AS->dalpha0_dTau();
            const double s_ig_over_R = AS->tau() * AS->dalpha0_dTau() - AS->alpha0();
            CHECK_THAT(h_ig_over_RT, Catch::Matchers::WithinAbs(0.0, 1e-8));
            CHECK_THAT(s_ig_over_R, Catch::Matchers::WithinAbs(0.0, 1e-8));
        }
    }
}

TEST_CASE("GERG pure fluid reports the GERG reducing point as its critical point", "[GERG]") {
    // Coverage for `fluid.crit = EOS.reduce` in make_gerg_fluid.  That
    // assignment is what feeds get_fluid_constant(i, iT_critical) /
    // irhomolar_critical, which is ALL the mixture branch of
    // calc_alpha0_deriv_nocache reads for T_ci / rho_ci -- but the pure branch
    // reads EOS().reduce instead, so deleting the assignment leaves every
    // other assertion in this file green while silently breaking Task 9.
    // AbstractState::T_critical()/rhomolar_critical()/p_critical() read
    // fluid.crit (HelmholtzEOSMixtureBackend::calc_*_critical, pure branch,
    // no superancillary available for a GERG fluid), so they pin it here.
    struct Case
    {
        const char* backend;
        const char* name;
        double Tc_K;
        double rhoc_molm3;
    };
    const std::vector<Case> cases = {
      {"GERG2008", "methane", 190.564, 10.139342719e3},
      {"GERG2004", "methane", 190.564, 10.139342719e3},
      {"GERG2008", "water", 647.096, 17.873716090e3},
      {"GERG2004", "helium", 5.1953, 17.399e3},
      // GERG-2008 moved the reducing point of both of these relative to
      // GERG-2004, so the pair also proves crit follows the MODEL's table
      // and not a single shared one.
      {"GERG2004", "carbonmonoxide", 132.800, 10.85e3},
      {"GERG2008", "carbonmonoxide", 132.860, 10.85e3},
      {"GERG2004", "isopentane", 460.350, 3.271018581e3},
      {"GERG2008", "isopentane", 460.350, 3.271e3},
      {"GERG2008", "n-decane", 617.7, 1.64e3},
    };
    for (const auto& c : cases) {
        CAPTURE(c.backend, c.name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory(c.backend, std::vector<std::string>{c.name}));
        CHECK_THAT(AS->T_critical(), Catch::Matchers::WithinRel(c.Tc_K, 1e-12));
        CHECK_THAT(AS->rhomolar_critical(), Catch::Matchers::WithinRel(c.rhoc_molm3, 1e-12));
        // crit must agree with reduce, since the two branches of
        // calc_alpha0_deriv_nocache read different ones.
        CHECK_THAT(AS->T_critical(), Catch::Matchers::WithinRel(AS->T_reducing(), 1e-15));
        CHECK_THAT(AS->rhomolar_critical(), Catch::Matchers::WithinRel(AS->rhomolar_reducing(), 1e-15));
        // p_critical() is inf unless reduce.p was filled in from the EOS.
        CHECK(ValidNumber(AS->p_critical()));
        CHECK(AS->p_critical() > 0);
    }

    // The pressure at the reducing point is evaluated from the GERG EOS
    // itself, so it must reproduce an independent single-phase evaluation at
    // (Tc, rhoc) -- and land near the published critical pressures.
    std::shared_ptr<AbstractState> CH4(AbstractState::factory("GERG2008", std::vector<std::string>{"methane"}));
    CH4->specify_phase(iphase_gas);
    CH4->update(DmolarT_INPUTS, CH4->rhomolar_reducing(), CH4->T_reducing());
    CHECK_THAT(CH4->p_critical(), Catch::Matchers::WithinRel(CH4->p(), 1e-12));
    CHECK_THAT(CH4->p_critical(), Catch::Matchers::WithinRel(4.5992e6, 1e-3));  // methane p_c, 4.5992 MPa
}

TEST_CASE("GERG applies the published range verbatim to every component", "[GERG]") {
    // 60-700 K is published for the MIXTURE model as a whole (Kunz & Wagner
    // 2012, section 4.1) and is applied unchanged to every component.  GERG
    // publishes no per-component lower limit, and inventing one -- capping
    // Tmin at the component's own Tc, as an earlier version of this backend
    // did -- fabricates a limit the model does not state.
    for (const auto& model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const std::string backend = (model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008";
        for (const auto& name : component_names(model)) {
            CAPTURE(backend, name);
            std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{name}));
            CHECK_THAT(AS->Tmin(), Catch::Matchers::WithinRel(60.0, 1e-14));
            CHECK_THAT(AS->Tmax(), Catch::Matchers::WithinRel(700.0, 1e-14));
        }
    }
}

TEST_CASE("GERG helium and hydrogen are supercritical-only as pure fluids", "[GERG]") {
    // The consequence of applying 60-700 K verbatim, and it is the physics
    // rather than a defect: helium critical at 5.1953 K and hydrogen at
    // 33.19 K are both SUPERCRITICAL everywhere in 60-700 K, so as pure
    // fluids they have no subcritical region here and no saturation state can
    // be requested.  A caller who asks gets the range in the error message,
    // which is a better answer than a silently empty dome.
    //
    // Both remain ordinary components in a MIXTURE, where the mixture
    // reducing temperature governs -- asserted below so this is not read as a
    // blanket restriction.
    for (const std::string& name : {std::string("helium"), std::string("hydrogen")}) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        CHECK(AS->T_critical() < AS->Tmin());
        AS->specify_phase(iphase_gas);
        CHECK_THROWS_AS(AS->update(QT_INPUTS, 0.5, 0.9 * AS->T_critical()), CoolProp::OutOfRangeError);
        CHECK_NOTHROW(AS->update(DmolarT_INPUTS, 500.0, 300.0));  // supercritical: fine
    }
    std::shared_ptr<AbstractState> mix(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "helium"}));
    mix->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    mix->specify_phase(iphase_gas);
    CHECK_NOTHROW(mix->update(DmolarT_INPUTS, 5000.0, 200.0));
}

TEST_CASE("GERG mixture range is the component-wise intersection, independent of the composition vector", "[GERG]") {
    // REPLACES an earlier test that pinned the mole-fraction-WEIGHTED range
    // as intentional.  It was not defensible: AbstractState::Tmin()/Tmax()
    // for a mixture are HelmholtzEOSMixtureBackend::calc_Tmin/calc_Tmax
    // (.cpp:1305-1319), `sum_i x_i * limits.T{min,max}` with NO normalisation
    // by sum(x) -- so a guard built on them scaled its own bounds with the
    // sum of the mole fractions and failed open.  Measured on this backend
    // before the fix, with GERG2008 methane+ethane and x = (0.5s, 0.5s):
    //
    //   s = 1   Tmin =  60.00  Tmax =  700.00   update(PT, 1e6, 900) threw
    //   s = 2   Tmin = 120.00  Tmax = 1400.00   update(PT, 1e6, 900) ACCEPTED,
    //                                           rho = 184.714 mol/m^3
    //   PropsSI("Dmolar","T",900,"P",1e6,"GERG2008::METHANE[0.6]&ETHANE[0.6]")
    //           -> "outside the GERG range of validity [72, 840] K"
    //
    // i.e. reachable straight from the string API, whose parser rejects an
    // individual fraction > 1 but not the sum.  check_gerg_range_of_validity
    // now reads the per-component EOS.limits directly and takes their
    // INTERSECTION (max of Tmin, min of Tmax); see its definition for why the
    // intersection rather than the union or a weighted average.
    //
    // Note Tmin()/Tmax() themselves are deliberately left alone -- the fix is
    // confined to src/Backends/GERG/ -- so this test asserts on what the
    // guard DOES, not on what Tmin() reports.
    {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "ethane"}));
        AS->set_mole_fractions(std::vector<CoolPropDbl>{1.0, 1.0});  // sum(x) = 2: the exploit above
        AS->specify_phase(iphase_gas);
        CHECK(AS->Tmax() > 1000);  // the weighted bound really is inflated...
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 1e6, 900.0),
                        CoolProp::OutOfRangeError);  // ...and the guard rejects 900 K anyway
    }
    {
        // The benign half of the same coupling: a gas analysis rounded to
        // sum(x) = 0.999999 used to move Tmax to 699.9993 K and throw a
        // spurious OutOfRangeError at exactly 700 K.
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "ethane"}));
        AS->set_mole_fractions(std::vector<CoolPropDbl>{0.4999995, 0.4999995});
        AS->specify_phase(iphase_gas);
        CHECK(AS->Tmax() < 700.0);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 700.0));
    }
    {
        // Same guard through the string API, where the exploit is reachable
        // without any downcast or direct set_mole_fractions call.  PropsSI
        // swallows the exception and returns _HUGE, so check both.
        const double rho = CoolProp::PropsSI("Dmolar", "T", 900, "P", 1e6, "GERG2008::METHANE[0.6]&ETHANE[0.6]");
        CHECK_FALSE(ValidNumber(rho));
        CHECK(CoolProp::get_global_param_string("errstring").find("[60, 700] K") != std::string::npos);
    }
    {
        // A 50/50 helium/methane mixture.  Every component now carries the
        // published 60 K, so the weighted average and the intersection agree
        // at 60 -- but the guard still reads the intersection, which is what
        // keeps it independent of sum(x); the methane+ethane case above is
        // where the two differ.
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"helium", "methane"}));
        AS->set_mole_fractions(std::vector<CoolPropDbl>{0.5, 0.5});
        AS->specify_phase(iphase_gas);
        CHECK_THAT(AS->Tmin(), Catch::Matchers::WithinRel(60.0, 1e-12));
        CHECK_THROWS_AS(AS->update(DmolarT_INPUTS, 5000.0, 45.0), CoolProp::OutOfRangeError);
        CHECK_THROWS_AS(AS->update(DmolarT_INPUTS, 5000.0, 30.0), CoolProp::OutOfRangeError);
        CHECK_NOTHROW(AS->update(DmolarT_INPUTS, 5000.0, 65.0));
    }
    {
        // There is no longer ANY composition whose range dips below 60 K.
        // The two lightest components together used to give [33.19, 700] K
        // via the old min(60, Tc) cap; both now carry the published 60 K, so
        // 40 K is refused here exactly as it is for methane.  Pinned because
        // the old behaviour was the last remnant of the fabricated limit.
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"helium", "hydrogen"}));
        AS->set_mole_fractions(std::vector<CoolPropDbl>{0.5, 0.5});
        AS->specify_phase(iphase_gas);
        CHECK_THROWS_AS(AS->update(DmolarT_INPUTS, 5000.0, 40.0), CoolProp::OutOfRangeError);
        CHECK_THROWS_AS(AS->update(DmolarT_INPUTS, 5000.0, 20.0), CoolProp::OutOfRangeError);
        CHECK_NOTHROW(AS->update(DmolarT_INPUTS, 5000.0, 65.0));
    }
}

TEST_CASE("GERG canonical-name fast path does not weaken model strictness", "[GERG]") {
    // resolve_component returns a canonical GERG name unchanged, which is the
    // path that makes component_names() round-trip through the factory.  That
    // early return must NOT let a GERG-2008-only component in through the
    // GERG-2004 model: the three fluids GERG-2008 adds are absent from
    // component_names(GERG_2004), so they miss the fast path and fall through
    // to the CAS route, which rejects them.
    for (const char* name : {"hydrogensulfide", "n-nonane", "n-decane"}) {
        CAPTURE(name);
        CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2004, name), CoolProp::ValueError);
        CHECK(resolve_component(GERGModel::GERG_2008, name) == std::string(name));
        // ... and the same through the factory, not just the resolver.
        CHECK_THROWS_AS(AbstractState::factory("GERG2004", std::vector<std::string>{name}), CoolProp::ValueError);
        CHECK_NOTHROW(std::shared_ptr<AbstractState>(AbstractState::factory("GERG2008", std::vector<std::string>{name})));
    }
    // A canonical-looking name that is not a GERG component at all still goes
    // down the CAS route and is still rejected.
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "r134a"), CoolProp::ValueError);
    CHECK_THROWS_AS(resolve_component(GERGModel::GERG_2008, "n-undecane"), CoolProp::ValueError);
}

TEST_CASE("GERG backend rejects an empty component list", "[GERG]") {
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{}), CoolProp::ValueError);
}

// ###########################################################################
// Mixtures (Task 9).  Two independent things are being pinned here: that the
// reducing function and the excess/departure term reproduce teqp, and that
// they were built from the GERG tables rather than from CoolProp's global
// binary-interaction-parameter library.  The second is not implied by the
// first only because the first is so broad -- it is stated separately below
// so that a future refactor deleting the override fails with a message that
// says what happened.
// ###########################################################################

TEST_CASE("GERG mixtures reproduce teqp", "[GERG]") {
    using namespace CoolProp::GERG::reference;
    PropertyLimitsSuspender suspend_limits;

    std::size_t trimmed_2004 = 0, trimmed_2008 = 0;

    RefTally t2004;
    compare_mix_points("GERG2004", mix_points_2004, t2004, trimmed_2004);
    // All C(18,2) = 153 GERG-2004 binary pairs at (250 K, 5000 mol/m^3,
    // z = [0.4, 0.6]) -- one numeric point per reducing-parameter row.  19 of
    // the 153 w values are NaN (spinodal branch); alphar/alphaig/p/cv never
    // are, which the generator asserts at generation time.
    CHECK(t2004.rows == 153);
    CHECK(t2004.nan_skipped_alphar == 0);
    CHECK(t2004.nan_skipped_alphaig == 0);
    CHECK(t2004.nan_skipped_p == 0);
    CHECK(t2004.nan_skipped_cvmolar == 0);
    CHECK(t2004.nan_skipped_w == 19);
    CHECK(t2004.compared == 153 * 5 - 19);
    // Every binary-pair row is a genuine two-component mixture with both mole
    // fractions nonzero, so none of them may be trimmed.
    CHECK(trimmed_2004 == 0);

    RefTally t2008;
    compare_mix_points("GERG2008", mix_points_2008, t2008, trimmed_2008);
    // C(21,2) = 210 binary pairs + 187 AGA8 natural-gas compositions; 40 NaN w.
    CHECK(t2008.rows == 210 + 187);
    CHECK(t2008.nan_skipped_alphar == 0);
    CHECK(t2008.nan_skipped_alphaig == 0);
    CHECK(t2008.nan_skipped_p == 0);
    CHECK(t2008.nan_skipped_cvmolar == 0);
    CHECK(t2008.nan_skipped_w == 40);
    CHECK(t2008.compared == (210 + 187) * 5 - 40);
    // All 187 AGA8 rows carry all 21 component names, most of them with a zero
    // mole fraction, and are trimmed; none of the 210 binary-pair rows is.
    CHECK(trimmed_2008 == 187);
}

TEST_CASE("GERG mixture parameters come from the GERG tables, not CoolProp's BIP library", "[GERG]") {
    // GERGMixtureBackend::set_mixture_parameters exists to keep
    // HelmholtzEOSMixtureBackend::set_mixture_parameters -- which resolves
    // every binary pair through CoolProp's GLOBAL binary-pair library by CAS
    // number -- from ever answering for a backend named GERG2004/GERG2008.
    // That library carries Bell-JCED-2016 refits that override the
    // Kunz-JCED-2012 (GERG) rows for pairs GERG also defines, so reaching it
    // would return plausible, non-GERG numbers with no error anywhere.
    //
    // This test does NOT rely on the interim "it throws because CoolPropFluid
    // ::CAS is empty" accident.  It recomputes the mixture reducing state from
    // the GERG monograph's own beta/gamma (via GERGData.h) and pins it, so
    // substituting ANY other parameter source -- throwing or not -- fails.
    struct Case
    {
        const char* backend;
        GERGModel model;
        const char* f1;
        const char* f2;
    };
    // methane/ethane and methane/n-butane are GERG pairs that CoolProp's own
    // library also defines (so the inherited path would silently succeed);
    // methane/propane exercises a pair with a fluid-specific departure
    // function; nitrogen/hydrogen has betaT != 1 in BOTH directions, which is
    // where a dropped reciprocal shows up.
    const std::vector<Case> cases = {
      {"GERG2008", GERGModel::GERG_2008, "methane", "ethane"},      {"GERG2008", GERGModel::GERG_2008, "methane", "n-butane"},
      {"GERG2008", GERGModel::GERG_2008, "methane", "propane"},     {"GERG2004", GERGModel::GERG_2004, "nitrogen", "hydrogen"},
      {"GERG2004", GERGModel::GERG_2004, "carbondioxide", "water"},
    };
    const double x1 = 0.35, x2 = 0.65;  // deliberately not 0.5/0.5: at equimolar the
                                        // beta-weighted term is symmetric in beta and
                                        // a swapped reciprocal is much harder to see.
    for (const auto& c : cases) {
        CAPTURE(c.backend, c.f1, c.f2);
        std::shared_ptr<AbstractState> AS(AbstractState::factory(c.backend, std::vector<std::string>{c.f1, c.f2}));
        AS->set_mole_fractions(std::vector<CoolPropDbl>{x1, x2});

        const BetasGammas bg = get_betasgammas(c.model, c.f1, c.f2);
        const PureInfo i1 = get_pure_info(c.model, c.f1);
        const PureInfo i2 = get_pure_info(c.model, c.f2);

        // GERG-2008 Eqs. 7.9/7.10 (Kunz & Wagner 2012), which is exactly what
        // GERG2008ReducingFunction implements.
        const double Tr = x1 * x1 * i1.Tc_K + x2 * x2 * i2.Tc_K
                          + 2 * x1 * x2 * bg.betaT * bg.gammaT * (x1 + x2) / (bg.betaT * bg.betaT * x1 + x2) * std::sqrt(i1.Tc_K * i2.Tc_K);
        const double vc12 = 1.0 / 8.0 * std::pow(std::pow(i1.rhoc_molm3, -1.0 / 3.0) + std::pow(i2.rhoc_molm3, -1.0 / 3.0), 3);
        const double vr =
          x1 * x1 / i1.rhoc_molm3 + x2 * x2 / i2.rhoc_molm3 + 2 * x1 * x2 * bg.betaV * bg.gammaV * (x1 + x2) / (bg.betaV * bg.betaV * x1 + x2) * vc12;

        CHECK_THAT(AS->T_reducing(), Catch::Matchers::WithinRel(Tr, 1e-14));
        CHECK_THAT(AS->rhomolar_reducing(), Catch::Matchers::WithinRel(1.0 / vr, 1e-14));

        // Both triangles of the beta/gamma matrices, read straight back out of
        // the reducing function.  These four-per-parameter assertions are NOT
        // redundant with the reducing-state pins above: every call site in
        // GERG2008ReducingFunction passes its index pair in ascending order
        // (f_Y_ij/dfYkidxi/dfYikdxi, ReducingFunctions.cpp:257, 306-322,
        // 550-563), so the LOWER triangle is currently read by nothing at all.
        // Mutation-verified: writing beta_T[j][i] = bg.betaT instead of
        // 1/bg.betaT leaves every other assertion in this file green.  A
        // future CoolProp change that starts reading beta[j][i] -- a
        // composition-derivative or fugacity path is the obvious candidate --
        // would silently inherit a wrong value without these.
        CHECK_THAT(AS->get_binary_interaction_double(0, 1, "betaT"), Catch::Matchers::WithinRel(bg.betaT, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(1, 0, "betaT"), Catch::Matchers::WithinRel(1.0 / bg.betaT, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(0, 1, "betaV"), Catch::Matchers::WithinRel(bg.betaV, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(1, 0, "betaV"), Catch::Matchers::WithinRel(1.0 / bg.betaV, 1e-15));
        // gammas are symmetric under exchange -- NOT reciprocated.
        CHECK_THAT(AS->get_binary_interaction_double(0, 1, "gammaT"), Catch::Matchers::WithinRel(bg.gammaT, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(1, 0, "gammaT"), Catch::Matchers::WithinRel(bg.gammaT, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(0, 1, "gammaV"), Catch::Matchers::WithinRel(bg.gammaV, 1e-15));
        CHECK_THAT(AS->get_binary_interaction_double(1, 0, "gammaV"), Catch::Matchers::WithinRel(bg.gammaV, 1e-15));
    }

    // ... and the pins above are discriminating, not tautological.  All 210
    // GERG-2008 binary pairs are present in CoolProp's global library, 194 of
    // them with the very same Kunz-JCED-2012 row GERG publishes -- but 16 with
    // a later refit (15 Gernert-Thesis-2013, 1 Tkaczuk-JPCRD-2020).  Twelve of
    // those 16 move the mixture reducing temperature by more than 0.03 K at
    // this composition; the four below move it by 1.5 K to 42 K, so the
    // inherited path could not possibly pass the pins above for them.  If any
    // of these ever stops differing, this test has to be re-pointed at a pair
    // that still does.
    //
    // The remaining four non-Kunz pairs (nitrogen/carbondioxide,
    // nitrogen/oxygen, oxygen/carbonmonoxide, carbonmonoxide/argon) carry
    // beta/gamma numerically equal to GERG's but a DIFFERENT departure
    // function, which the reducing state cannot see.  Those are caught instead
    // by the alphar column of "GERG mixtures reproduce teqp", which compares
    // every one of the 210 pairs at 1e-12.
    for (const auto& pair :
         std::vector<std::vector<std::string>>{{"water", "argon"}, {"water", "nitrogen"}, {"carbondioxide", "water"}, {"helium", "argon"}}) {
        CAPTURE(pair[0], pair[1]);
        std::shared_ptr<AbstractState> G(AbstractState::factory("GERG2008", pair));
        std::shared_ptr<AbstractState> H(AbstractState::factory("HEOS", pair));
        G->set_mole_fractions(std::vector<CoolPropDbl>{x1, x2});
        H->set_mole_fractions(std::vector<CoolPropDbl>{x1, x2});
        CHECK(std::abs(H->T_reducing() - G->T_reducing()) > 1.0);
    }
}

TEST_CASE("GERG mixture linked saturation states are GERG-typed", "[GERG]") {
    // HelmholtzEOSMixtureBackend::set_components builds SatL/SatV with an
    // EXPLICITLY qualified HelmholtzEOSMixtureBackend::get_copy(false)
    // (HelmholtzEOSMixtureBackend.cpp:142, :146), so a get_copy override alone
    // cannot change their type -- GERGMixtureBackend overrides set_components
    // instead.  Were that override removed, the base-class SatL would run the
    // INHERITED set_mixture_parameters during its own construction, i.e. the
    // fail-open this backend exists to prevent (or, with CAS still empty, a
    // throw out of the factory).  Either way this test fails.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "ethane"}));
    auto* gerg = dynamic_cast<GERGMixtureBackend*>(AS.get());
    REQUIRE(gerg != nullptr);
    CHECK(gerg->backend_name() == "GERG2008Backend");
    CHECK(gerg->get_SatL().backend_name() == "GERG2008Backend");
    CHECK(gerg->get_SatV().backend_name() == "GERG2008Backend");

    // And they evaluate identically to the parent, ideal-gas part included --
    // the linked states carry the reducing/excess objects as DATA (via
    // sync_linked_states) but the ideal-gas path is re-derived from their own
    // components, so this is the assertion that would catch a divergence
    // between parent and linked state.
    const std::vector<CoolPropDbl> z{0.35, 0.65};
    AS->set_mole_fractions(z);
    AS->specify_phase(iphase_gas);
    AS->update(DmolarT_INPUTS, 4000.0, 260.0);
    for (auto* sat : {&gerg->get_SatL(), &gerg->get_SatV()}) {
        sat->set_mole_fractions(z);
        sat->specify_phase(iphase_gas);
        sat->update(DmolarT_INPUTS, 4000.0, 260.0);
        CHECK_THAT(sat->alphar(), Catch::Matchers::WithinRel(AS->alphar(), 1e-15));
        CHECK_THAT(sat->alpha0(), Catch::Matchers::WithinRel(AS->alpha0(), 1e-15));
        CHECK_THAT(sat->p(), Catch::Matchers::WithinRel(AS->p(), 1e-15));
    }
}

TEST_CASE("GERG mixture state can be copied", "[GERG]") {
    // ExcessTerm::copy() dereferences EVERY off-diagonal departure-function
    // pointer (ExcessHEFunction.h:215-227), including for the majority of
    // pairs that carry no departure function at all -- methane/water is one,
    // F_ij is not even in the table -- so a null there is a segfault on copy,
    // not a zero.  set_mixture_parameters installs a zero-valued
    // ExponentialDepartureFunction placeholder for exactly this reason.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "water"}));
    AS->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    AS->specify_phase(iphase_gas);
    CHECK_NOTHROW(AS->update(DmolarT_INPUTS, 5000.0, 400.0));

    auto* gerg = dynamic_cast<GERGMixtureBackend*>(AS.get());
    REQUIRE(gerg != nullptr);
    std::shared_ptr<HelmholtzEOSMixtureBackend> copy;
    REQUIRE_NOTHROW(copy.reset(gerg->get_copy(true)));
    // get_copy is overridden too, so TPD_state / critical_state /
    // transient_pure_state (which call it unqualified) stay GERG-typed.
    CHECK(copy->backend_name() == "GERG2008Backend");
    copy->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    copy->specify_phase(iphase_gas);
    copy->update(DmolarT_INPUTS, 5000.0, 400.0);
    CHECK_THAT(copy->p(), Catch::Matchers::WithinRel(AS->p(), 1e-15));
    CHECK_THAT(copy->alphar(), Catch::Matchers::WithinRel(AS->alphar(), 1e-15));
    CHECK_THAT(copy->alpha0(), Catch::Matchers::WithinRel(AS->alpha0(), 1e-15));
}

TEST_CASE("GERG-2004 and GERG-2008 disagree where the models disagree", "[GERG]") {
    // Carbon monoxide's pure EOS and reducing state both changed between the
    // two models, so the two backends must not be silently aliased.
    std::shared_ptr<AbstractState> a(AbstractState::factory("GERG2004", std::vector<std::string>{"CarbonMonoxide"}));
    std::shared_ptr<AbstractState> b(AbstractState::factory("GERG2008", std::vector<std::string>{"CarbonMonoxide"}));
    a->specify_phase(iphase_gas);
    b->specify_phase(iphase_gas);
    a->update(DmolarT_INPUTS, 5000.0, 200.0);
    b->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK(a->p() != b->p());

    // Methane is unchanged between the two models.
    std::shared_ptr<AbstractState> c(AbstractState::factory("GERG2004", std::vector<std::string>{"Methane"}));
    std::shared_ptr<AbstractState> d(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    c->specify_phase(iphase_gas);
    d->specify_phase(iphase_gas);
    c->update(DmolarT_INPUTS, 5000.0, 200.0);
    d->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK_THAT(c->p(), Catch::Matchers::WithinRel(d->p(), 1e-15));

    // The same for a MIXTURE: GERG-2008 revises the isopentane/carbonmonoxide
    // reducing parameters (gammaV 1.116693501 -> 1.116694577, gammaT
    // 1.199475627 -> 1.199326059) as well as both pure EOS, so the two
    // backends must differ on this pair too -- and must NOT differ on a pair
    // GERG-2008 leaves alone.
    for (const auto& pair : std::vector<std::vector<std::string>>{{"isopentane", "carbonmonoxide"}, {"methane", "nitrogen"}}) {
        CAPTURE(pair[0], pair[1]);
        std::shared_ptr<AbstractState> m4(AbstractState::factory("GERG2004", pair));
        std::shared_ptr<AbstractState> m8(AbstractState::factory("GERG2008", pair));
        for (auto* s : {m4.get(), m8.get()}) {
            s->set_mole_fractions(std::vector<CoolPropDbl>{0.4, 0.6});
            s->specify_phase(iphase_gas);
            s->update(DmolarT_INPUTS, 3000.0, 300.0);
        }
        if (pair[0] == "methane") {
            CHECK_THAT(m4->p(), Catch::Matchers::WithinRel(m8->p(), 1e-15));
            CHECK_THAT(m4->T_reducing(), Catch::Matchers::WithinRel(m8->T_reducing(), 1e-15));
        } else {
            CHECK(m4->p() != m8->p());
            CHECK(m4->T_reducing() != m8->T_reducing());
        }
    }
}

// ###########################################################################
// Task 10: strictness.  Every guard below is verified by mutation -- see
// task-10-report.md for exactly which assertion fails, and how many, when the
// guard is removed.
// ###########################################################################

TEST_CASE("GERG rejects components outside the model", "[GERG]") {
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"R134a"}), CoolProp::ValueError);
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "R134a"}), CoolProp::ValueError);
    CHECK_THROWS_AS(AbstractState::factory("GERG2004", std::vector<std::string>{"n-Decane"}), CoolProp::ValueError);
    // AbstractState::factory returns a raw *owning* pointer, so the success
    // case has to hand it to a smart pointer or the whole backend (plus its
    // linked SatL/SatV states, reducing function and departure functions)
    // leaks.  The CHECK_THROWS_AS lines above are safe only because the
    // factory throws before it ever returns.  Same idiom as the n-decane
    // check in "GERG component resolution" above.
    CHECK_NOTHROW(std::shared_ptr<AbstractState>(AbstractState::factory("GERG2008", std::vector<std::string>{"n-Decane"})));
}

TEST_CASE("GERG rejects a name CoolProp cannot resolve at all", "[GERG]") {
    // Distinct from "known to CoolProp but outside the GERG set" above: this
    // name fails even the CAS lookup resolve_component falls back to.
    CHECK_THROWS_AS(AbstractState::factory("GERG2008", std::vector<std::string>{"NOT A FLUID"}), CoolProp::ValueError);
}

TEST_CASE("GERG refuses transport properties", "[GERG]") {
    // AbstractState::viscosity()/conductivity()/surface_tension()
    // (AbstractState.cpp:793-813) each do `if (!_cached) _cached = calc_*();
    // return _cached;` -- there is no separate cache that could satisfy the
    // call without invoking calc_viscosity() etc, so overriding the calc_*
    // hook alone is sufficient here (task-10-brief.md hazard 3).  This test
    // is what a mutation deleting one of the three overrides would fail:
    // without the override, these three calls would return CoolProp's own
    // (non-GERG) transport-correlation numbers instead of throwing.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    AS->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK_THROWS_AS(AS->viscosity(), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->conductivity(), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->surface_tension(), CoolProp::NotImplementedError);
}

TEST_CASE("GERG refuses binary interaction parameter mutation through every public route", "[GERG]") {
    // AbstractState declares FOUR public mutator routes to the BIPs:
    // index-keyed and CAS-keyed set_binary_interaction_double, and the same
    // pair for set_binary_interaction_string (AbstractState.h:950-965).
    // HelmholtzEOSMixtureBackend overrides only the two index-keyed ones
    // (HelmholtzEOSMixtureBackend.h:217,223); the two CAS-keyed overloads
    // fall through, UNOVERRIDDEN by this class too, to AbstractState's own
    // default body, which already throws NotImplementedError unconditionally.
    // All four are exercised here so that a future HelmholtzEOSMixtureBackend
    // change that overrides the CAS-keyed forms (and answers from the global
    // BIP library) would flip this test's exception type/outcome instead of
    // reopening the hole silently.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "Nitrogen"}));
    CHECK_THROWS_AS(AS->set_binary_interaction_double(0, 1, "betaT", 1.1), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_binary_interaction_double(0, 1, "gammaV", 1.1), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_binary_interaction_string(0, 1, "betaT", "1.1"), CoolProp::ValueError);
    // "betaT" is not a string set_binary_interaction_string understands even
    // in the INHERITED HelmholtzEOSMixtureBackend implementation (it only
    // recognises "function"), so that call above would throw ValueError even
    // WITHOUT this backend's override -- it does not, by itself, prove the
    // override does anything. "function" is the parameter that matters: the
    // inherited body resolves it through get_departure_function(value)
    // against CoolProp's mixture-departure-function LIBRARY
    // (MixtureParameters.cpp:513), and "Methane-Nitrogen" is a real entry in
    // that library (dev/mixtures/mixture_departure_functions.json) -- so
    // without the override this call would SUCCEED and silently swap in a
    // library departure function for the pair, exactly the kind of
    // substitution set_mixture_parameters exists to prevent.
    CHECK_THROWS_AS(AS->set_binary_interaction_string(0, 1, "function", "Methane-Nitrogen"), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->set_binary_interaction_double("74-82-8", "7727-37-9", "betaT", 1.1), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(AS->set_binary_interaction_string("74-82-8", "7727-37-9", "betaT", "1.1"), CoolProp::NotImplementedError);

    // apply_simple_mixing_rule is a fifth nominal route
    // (HelmholtzEOSMixtureBackend.cpp:383): it calls
    // set_binary_interaction_double(i, j, ...) UNQUALIFIED, so virtual
    // dispatch already sends it through the index-keyed override above with
    // no separate override needed. Covered explicitly so a refactor that
    // routes apply_simple_mixing_rule around virtual dispatch (e.g. calling
    // HelmholtzEOSMixtureBackend::set_binary_interaction_double directly)
    // would be caught here.
    CHECK_THROWS_AS(AS->apply_simple_mixing_rule(0, 1, "linear"), CoolProp::ValueError);

    // Known, documented, NOT-closed bypasses (task-10-brief.md hazard 2; full
    // enumeration in GERGBackend.h). Both are a CHOICE not to close, not an
    // impossibility -- see GERGBackend.h for the ConstantReducingFunction
    // precedent that shows the route -- and both are pinned here (not just
    // claimed in a comment) so a future change that closes or reopens them
    // shows up as a test change:
    //
    // 1. `Reducing` is a PUBLIC member of HelmholtzEOSMixtureBackend
    //    (HelmholtzEOSMixtureBackend.h:147); a caller holding a
    //    HelmholtzEOSMixtureBackend* can reach the reducing function directly
    //    and mutate it in place.
    // 2. `residual_helmholtz` is equally public, and `Excess.F[i][j]` /
    //    `Excess.DepartureFunctionMatrix[i][j]` are themselves public members
    //    with no setter at all -- this is the identical "swap in a different
    //    departure function" hazard the "function" guard above exists to
    //    prevent, reached without going through any guarded setter.
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    CHECK_NOTHROW(heos->Reducing->set_binary_interaction_double(0, 1, "betaT", 1.1));
    AS->set_mole_fractions(std::vector<CoolPropDbl>{0.4, 0.6});
    AS->specify_phase(iphase_gas);
    const double alphar_before = [&] {
        AS->update(DmolarT_INPUTS, 5000.0, 200.0);
        return AS->alphar();
    }();
    heos->residual_helmholtz->Excess.F[0][1] = 0.0;
    heos->residual_helmholtz->Excess.F[1][0] = 0.0;
    AS->update(DmolarT_INPUTS, 5000.0, 200.0);
    CHECK(AS->alphar() != alphar_before);  // silently changed the mixture's alphar with no error
}

TEST_CASE("GERG fluids carry no superancillary", "[GERG]") {
    // Load-bearing, not cosmetic: FlashRoutines::sat_superanc_path_applies
    // (FlashRoutines.cpp:558) routes pure-fluid saturation straight to the
    // Chebyshev expansion and returns THAT as the answer for any pure fluid
    // that owns one. A GERG fluid carrying CoolProp's superancillary blob
    // (e.g. by construction ever changing to build one from the CoolProp
    // fluid library instead of leaving it null) would silently return
    // Setzmann-Wagner-class saturation densities labelled GERG-2008 -- with
    // no error anywhere, since the superancillary path is a deliberate
    // shortcut, not a bug.
    // The factory result is owned FIRST and cast SECOND: casting first and
    // then adopting the cast pointer leaks the allocation (and null-derefs on
    // the next line) if the cast ever fails.
    std::shared_ptr<AbstractState> owned(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    auto* gerg = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(owned.get());
    REQUIRE(gerg != nullptr);
    CHECK(gerg->get_superanc() == nullptr);
}

TEST_CASE("GERG respects its range of validity", "[GERG]") {
    // The range guard is GERGMixtureBackend::check_gerg_range_of_validity,
    // called from the update() override -- gated behind the
    // DONT_CHECK_PROPERTY_LIMITS configuration key (default false -- checks
    // ARE performed). Asserted explicitly here (task-10-brief.md hazard 4)
    // rather than merely assumed, so this test cannot pass only because some
    // earlier test in the same binary left the global config flipped.
    REQUIRE(CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS) == false);

    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    // Original review round of this test used bare CHECK_THROWS at (1e5, 40.0)
    // / (1e5, 900.0) with no imposed phase. That PASSED even with the guard's
    // lower bound deleted entirely (`if (_T > Thi)` instead of
    // `if (_T < Tlo || _T > Thi)`): at 40 K with no phase imposed, PT_flash's
    // phase determination lands in solver_rho_Tp's LIQUID branch, which calls
    // components[0].ancillaries.rhoL.evaluate(T) on an ancillary GERG fluids
    // never populate -- the exact missing-ancillary accident documented on
    // the `update` override and in task-10-report.md section 5(b) -- so the
    // "below Tmin" case was throwing for a reason that had nothing to do with
    // check_gerg_range_of_validity, and a bare CHECK_THROWS could not tell the
    // difference. Forcing iphase_gas routes both calls through solver_rho_Tp's
    // ideal-gas-guess branch instead, which does not itself throw at either
    // bound, so only the explicit T check can make these throw -- and
    // CHECK_THROWS_AS pins the exact exception type so a differently-thrown
    // error (e.g. the ancillary accident, or the p<0 catch-all) cannot be
    // mistaken for this guard either.
    AS->specify_phase(iphase_gas);
    CHECK_THROWS_AS(AS->update(PT_INPUTS, 1e3, 59.0), CoolProp::OutOfRangeError);   // below Tmin = 60 K
    CHECK_THROWS_AS(AS->update(PT_INPUTS, 1e5, 701.0), CoolProp::OutOfRangeError);  // above Tmax = 700 K
}

TEST_CASE("GERG closes change_EOS and update_with_guesses; downcast-only bypasses stay open and pinned", "[GERG]") {
    // Task-10 review round 2: change_EOS and update_with_guesses are both
    // virtual on AbstractState itself, so both are reachable from a plain
    // AbstractState* with NO downcast -- the exact shape AbstractState::factory
    // hands back -- and both were confirmed live to accept out-of-model input
    // before these overrides existed (change_EOS installed an SRK cubic;
    // update_with_guesses accepted T = 900 K). Closed here rather than left
    // open or merely documented, per the review.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    CHECK_THROWS_AS(AS->change_EOS(0, "SRK"), CoolProp::ValueError);
    CHECK_THROWS_AS(AS->change_EOS(0, "Peng-Robinson"), CoolProp::ValueError);

    GuessesStructure guess;
    guess.rhomolar = 1e5 / (AS->gas_constant() * 900.0);  // ideal-gas guess, so the failure is the range guard, not a bad initial guess
    CHECK_THROWS_AS(AS->update_with_guesses(PT_INPUTS, 1e5, 900.0, guess), CoolProp::OutOfRangeError);

    // Known, documented, NOT-closed bypasses that DO require a downcast --
    // pinned here (not just asserted in a comment) so a future change that
    // closes or reopens them shows up as a test change, not a silent drift.
    // These are deliberately weaker exposure than the two above: reaching
    // them needs source code that already knows it holds a concrete
    // HelmholtzEOSMixtureBackend-family object, not merely an AbstractState.
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    CHECK_NOTHROW(heos->update_DmolarT_direct(5000.0, 900.0));  // T = 900 K, well above Tmax = 700 K
    CHECK_NOTHROW(heos->update_TP_guessrho(900.0, 1e5, 1e5 / (8.314472 * 900.0)));
    CHECK_NOTHROW(heos->update_TDmolarP_unchecked(900.0, 5000.0, 1e5));

    // The other documented-but-not-closed bypasses (Reducing, residual_helmholtz->Excess)
    // are pinned in "GERG refuses binary interaction parameter mutation through
    // every public route" above.
}

// ###########################################################################
// Task 11: saturation ancillaries and VLE
// ###########################################################################

namespace {

/// The saturation curve of a GERG pure fluid is only reachable where BOTH
/// limits agree it is: at or above the backend's Tmin (make_gerg_fluid's
/// min(60 K, Tc) -- so helium and hydrogen, whose reducing temperatures are
/// 5.20 K and 33.19 K, have NO reachable saturation states at all, since
/// Tmin == Tc for them) and at or above the low-temperature end of the fitted
/// ancillary range (EOS.sat_min_liquid.T, which for the heavy alkanes is far
/// above 60 K because the saturation pressure there underflows). Returning
/// false is a legitimate skip, not a silently weakened test: the count of
/// fluids that DID run is asserted below so this helper cannot quietly
/// disable the whole sweep.
bool saturation_reachable(CoolProp::AbstractState& AS, double T) {
    CoolPropDbl Tmin_satL = NAN, Tmin_satV = NAN;
    auto* heos = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(&AS);
    // Every caller passes a GERGMixtureBackend, so this cast cannot fail
    // today; REQUIRE rather than an early return so that a future caller
    // passing something else fails the test loudly instead of silently
    // reporting "saturation unreachable" and skipping the sweep.
    REQUIRE(heos != nullptr);
    heos->calc_Tmin_sat(Tmin_satL, Tmin_satV);
    const double T_lo = std::max({static_cast<double>(Tmin_satL), static_cast<double>(Tmin_satV), static_cast<double>(AS.Tmin())});
    return T > T_lo && T < AS.T_critical();
}

/// Rebuild one of a fluid's three ancillaries with its output scaled by
/// (1 + rel), leaving everything else alone.
///
/// Scaling `reducing_value` is a uniform relative shift of the ancillary's
/// OUTPUT at every temperature -- for the two exponential forms
/// (reducing_value*exp(...)) and the non-exponential one
/// (reducing_value*(1 + sum)) alike -- which is exactly the "wrong by x%
/// everywhere" perturbation a refit-gone-bad would produce, and much cleaner
/// than nudging one n_i (whose effect varies wildly with T).
CoolProp::SaturationAncillaryFunction perturbed_ancillary(const CoolProp::GERG::AncillaryCoeffs& anc, double rel) {
    CoolProp::SaturationAncillaryFunction::Values v;
    v.type = (anc.type == "rhoLnoexp") ? CoolProp::SaturationAncillaryFunction::TYPE_NOT_EXPONENTIAL
                                       : CoolProp::SaturationAncillaryFunction::TYPE_EXPONENTIAL;
    v.using_tau_r = (anc.type != "rhoLnoexp");
    v.n = anc.n;
    v.t = anc.t;
    v.Tmin = anc.Tmin;
    v.Tmax = anc.Tmax;
    v.reducing_value = anc.reducing_value * (1.0 + rel);
    v.T_r = anc.T_r;
    return CoolProp::SaturationAncillaryFunction(v);
}

}  // namespace

TEST_CASE("GERG ancillary accessors reject bad arguments", "[GERG]") {
    CHECK_THROWS_AS(get_ancillary(GERGModel::GERG_2008, "methane", "rhoM"), CoolProp::ValueError);
    CHECK_THROWS_AS(get_ancillary(GERGModel::GERG_2004, "n-decane", "rhoL"), CoolProp::ValueError);
    CHECK_THROWS_AS(get_sat_min_state(GERGModel::GERG_2004, "hydrogensulfide"), CoolProp::ValueError);
    CHECK_NOTHROW(get_ancillary(GERGModel::GERG_2008, "n-decane", "rhoL"));

    // The accessor, not the table, owns `type`, so the two can never disagree.
    CHECK(get_ancillary(GERGModel::GERG_2008, "methane", "rhoL").type == "rhoLnoexp");
    CHECK(get_ancillary(GERGModel::GERG_2008, "methane", "rhoV").type == "rhoV");
    CHECK(get_ancillary(GERGModel::GERG_2008, "methane", "pV").type == "pV");

    // Every component of both models has all three ancillaries and a
    // low-temperature end state -- no fluid may be silently missing one.
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            for (const auto* which : {"rhoL", "rhoV", "pV"}) {
                CAPTURE(which);
                AncillaryCoeffs anc;
                REQUIRE_NOTHROW(anc = get_ancillary(model, name, which));
                CHECK(anc.n.size() == anc.t.size());
                CHECK(anc.n.size() >= 4);
                CHECK(anc.reducing_value > 0);
                CHECK(anc.T_r > 0);
                CHECK(anc.Tmin > 0);
                CHECK(anc.Tmin < anc.T_r);
                // The ship gate from the task brief: 1% relative. The fitter
                // aborts at 0.5%, so this is the belt to its braces -- a
                // hand-edited table row cannot slip past both.
                CHECK(anc.max_rel_dev < 0.01);
            }
            CHECK(get_sat_min_state(model, name).p_Pa > 0);
        }
    }
}

TEST_CASE("GERG-2004 and GERG-2008 carry separately fitted ancillaries where the EOS differs", "[GERG]") {
    // GERG-2008 moved carbon monoxide's and isopentane's reducing parameters
    // (GERGData.h pure_info_2008_overrides), which moves the whole saturation
    // curve, so each model needs its own fit. A single shared table would be
    // numerically plausible and wrong; these checks are what would fail.
    for (const auto& name : {"carbonmonoxide", "isopentane"}) {
        CAPTURE(name);
        const AncillaryCoeffs a4 = get_ancillary(GERGModel::GERG_2004, name, "rhoL");
        const AncillaryCoeffs a8 = get_ancillary(GERGModel::GERG_2008, name, "rhoL");
        CHECK(a4.T_r != a8.T_r);
        CHECK(a4.n != a8.n);
        // ...and the reducing density each fit is written against must be the
        // one its own model tabulates.
        CHECK_THAT(a4.reducing_value, Catch::Matchers::WithinRel(get_pure_info(GERGModel::GERG_2004, name).rhoc_molm3, 1e-14));
        CHECK_THAT(a8.reducing_value, Catch::Matchers::WithinRel(get_pure_info(GERGModel::GERG_2008, name).rhoc_molm3, 1e-14));
    }

    // Every other component is shared, and must be byte-identical between the
    // two models -- a diverging row would mean the 2008 override table had
    // grown an entry that does not belong to it.
    for (const auto& name : component_names(GERGModel::GERG_2004)) {
        if (name == "carbonmonoxide" || name == "isopentane") continue;
        CAPTURE(name);
        CHECK(get_ancillary(GERGModel::GERG_2004, name, "pV").n == get_ancillary(GERGModel::GERG_2008, name, "pV").n);
    }
}

TEST_CASE("GERG fluids evaluate their ancillaries through CoolProp's own evaluator", "[GERG]") {
    // make_gerg_fluid builds a CoolProp SaturationAncillaryFunction from each
    // table row; GERG::evaluate_ancillary is a standalone mirror of
    // SaturationAncillaryFunction::evaluate used by the tests. If those two
    // drift apart, every ancillary-accuracy assertion below silently stops
    // describing the function the flash routines actually call. Pin them
    // together, on the fluid the backend really built.
    for (const auto& name : component_names(GERGModel::GERG_2008)) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        auto* heos = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(heos != nullptr);
        CoolProp::CoolPropFluid& fluid = heos->get_components()[0];

        const AncillaryCoeffs aL = get_ancillary(GERGModel::GERG_2008, name, "rhoL");
        const AncillaryCoeffs aV = get_ancillary(GERGModel::GERG_2008, name, "rhoV");
        const AncillaryCoeffs ap = get_ancillary(GERGModel::GERG_2008, name, "pV");
        for (double frac : {0.4, 0.6, 0.8, 0.95}) {
            const double T = aL.Tmin + frac * (aL.T_r - aL.Tmin);
            CAPTURE(T);
            CHECK_THAT(fluid.ancillaries.rhoL.evaluate(T), Catch::Matchers::WithinRel(evaluate_ancillary(aL, T), 1e-14));
            CHECK_THAT(fluid.ancillaries.rhoV.evaluate(T), Catch::Matchers::WithinRel(evaluate_ancillary(aV, T), 1e-14));
            CHECK_THAT(fluid.ancillaries.pV.evaluate(T), Catch::Matchers::WithinRel(evaluate_ancillary(ap, T), 1e-14));
            // pL and pV are the same curve: GERG has no pseudo-pure components.
            CHECK_THAT(fluid.ancillaries.pL.evaluate(T), Catch::Matchers::WithinRel(evaluate_ancillary(ap, T), 1e-14));
        }
        // Above T_r the evaluator must return NaN rather than
        // pow(negative, fractional) (GitHub #1611).
        CHECK(std::isnan(evaluate_ancillary(aL, aL.T_r * 1.001)));
    }
}

TEST_CASE("GERG carries no superancillary but does carry the saturation end states", "[GERG]") {
    // The two halves belong together: attaching ancillaries must NOT have
    // brought a superancillary along, because FlashRoutines::sat_superanc_path_applies
    // (FlashRoutines.cpp:558) returns the superancillary AS THE ANSWER, so a
    // GERG fluid carrying CoolProp's blob would report Setzmann-Wagner
    // saturation densities labelled GERG-2008.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    auto* heos = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(heos != nullptr);
    CHECK(heos->get_superanc() == nullptr);

    const SatEndState end = get_sat_min_state(GERGModel::GERG_2008, "methane");
    CoolPropDbl Tmin_satL = NAN, Tmin_satV = NAN, pmin_satL = NAN, pmin_satV = NAN;
    heos->calc_Tmin_sat(Tmin_satL, Tmin_satV);
    heos->calc_pmin_sat(pmin_satL, pmin_satV);
    CHECK_THAT(static_cast<double>(Tmin_satL), Catch::Matchers::WithinRel(end.T_K, 1e-14));
    CHECK_THAT(static_cast<double>(Tmin_satV), Catch::Matchers::WithinRel(end.T_K, 1e-14));
    CHECK_THAT(static_cast<double>(pmin_satL), Catch::Matchers::WithinRel(end.p_Pa, 1e-14));
    CHECK_THAT(static_cast<double>(pmin_satV), Catch::Matchers::WithinRel(end.p_Pa, 1e-14));

    // triple_liquid/triple_vapor get the same state. GERG publishes no triple
    // point (EOS.Ttriple stays 0); these two slots are read by
    // saturation_T_pure_Maxwell purely as the low-T end of the saturation
    // curve. Left at SimpleState's _HUGE default, every ancillary seed would
    // fail Maxwell's sanity band and every saturation call would take the
    // linear fallback -- which is why this is asserted rather than assumed.
    CHECK_THAT(AS->get_state("triple_liquid").rhomolar, Catch::Matchers::WithinRel(end.rhoL_molm3, 1e-14));
    CHECK_THAT(AS->get_state("triple_vapor").rhomolar, Catch::Matchers::WithinRel(end.rhoV_molm3, 1e-14));
    CHECK(ValidNumber(AS->get_state("triple_liquid").hmolar));

    // hs_anchor is deliberately NOT populated.  Every site that reads its
    // h/s (VLERoutines.cpp:200-294) also reads components[0].ancillaries.hL
    // or .sL, which GERG never fills, so the anchor could not be used even if
    // it were set.  Asserted rather than left implicit: a future change that
    // starts populating it should have to justify itself here.
    const CoolProp::SimpleState& anchor = AS->get_state("hs_anchor");
    CHECK_FALSE(ValidNumber(anchor.T));
    CHECK_FALSE(ValidNumber(anchor.hmolar));

    // reduce/critical h and s ARE filled, because get_state() exposes them.
    CHECK(ValidNumber(AS->get_state("reducing").smolar));
    CHECK(ValidNumber(AS->get_state("critical").hmolar));
}

TEST_CASE("GERG constructs the fluids whose anchor temperature was out of range", "[GERG]") {
    // REPLACES a test that pinned hs_anchor's clamping.  Dropping the anchor
    // removed the clamp and the hazard together, so what is left worth
    // asserting is that the three fluids the anchor used to break on still
    // construct and evaluate.
    //
    //   water    1.1*647.096 K = 711.8 K, ABOVE the 700 K Tmax
    //   hydrogen 1.1*33.19 K   =  36.5 K, BELOW the 60 K Tmin
    //   helium   1.1*5.1953 K  =   5.7 K, BELOW the 60 K Tmin
    //
    // Every one of them made update_states() throw OutOfRangeError during
    // construction, so the fluid could not be built at all -- water via the
    // upper bound, the other two via the lower bound once the fabricated
    // per-component Tmin cap was removed.
    for (const std::string& name : {std::string("Water"), std::string("Hydrogen"), std::string("Helium")}) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS;
        REQUIRE_NOTHROW(AS.reset(AbstractState::factory("GERG2008", std::vector<std::string>{name})));
        CHECK_FALSE(ValidNumber(AS->get_state("hs_anchor").T));
        CHECK(ValidNumber(AS->get_state("reducing").hmolar));
        AS->specify_phase(iphase_gas);
        REQUIRE_NOTHROW(AS->update(DmolarT_INPUTS, 500.0, 400.0));
        CHECK(ValidNumber(AS->hmolar()));
    }
}

TEST_CASE("GERG pure fluid saturation converges", "[GERG]") {
    // The ancillaries seed the iteration; the converged answer is pure GERG.
    // The consistency check is what distinguishes the two: evaluating the EOS
    // independently at (T, rhoL) must reproduce the saturation pressure the
    // flash returned. An implementation that handed the ancillary back as the
    // answer would miss by the fit's own deviation, ~1e-3 relative, which is
    // five orders of magnitude outside the tolerance below.
    for (const auto& name : {"Methane", "Ethane", "Propane", "Nitrogen", "CarbonDioxide"}) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        const double T = 0.8 * AS->T_critical();
        REQUIRE(saturation_reachable(*AS, T));
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.0, T));
        const double p_sat = AS->p();
        const double rhoL = AS->rhomolar();
        CHECK(rhoL > AS->rhomolar_critical());

        std::shared_ptr<AbstractState> chk(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        chk->specify_phase(iphase_liquid);
        chk->update(DmolarT_INPUTS, rhoL, T);
        CHECK_THAT(chk->p(), Catch::Matchers::WithinRel(p_sat, 1e-8));

        // Q = 1 must land on the same pressure from the other side of the dome.
        std::shared_ptr<AbstractState> vap(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        vap->update(QT_INPUTS, 1.0, T);
        CHECK_THAT(vap->p(), Catch::Matchers::WithinRel(p_sat, 1e-8));
        CHECK(vap->rhomolar() < AS->rhomolar_critical());
        std::shared_ptr<AbstractState> chkv(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        chkv->specify_phase(iphase_gas);
        chkv->update(DmolarT_INPUTS, vap->rhomolar(), T);
        CHECK_THAT(chkv->p(), Catch::Matchers::WithinRel(p_sat, 1e-8));
    }
}

TEST_CASE("GERG saturation converges for every component whose dome is in range", "[GERG]") {
    // Sweeps both models. Helium and hydrogen are expected to be skipped
    // entirely: make_gerg_fluid clamps Tmin to min(60 K, Tc), which for those
    // two IS Tc, so they have no reachable subcritical state at all. The
    // ran/skipped counts are asserted so that a future change which quietly
    // makes MOST fluids unreachable fails here instead of passing vacuously.
    int ran = 0, skipped = 0;
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const std::string backend = (model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008";
        for (const auto& name : component_names(model)) {
            CAPTURE(backend);
            CAPTURE(name);
            std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{name}));
            const double T = 0.8 * AS->T_critical();
            if (!saturation_reachable(*AS, T)) {
                ++skipped;
                continue;
            }
            ++ran;
            REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.0, T));
            const double p_sat = AS->p(), rhoL = AS->rhomolar();
            CHECK(p_sat > 0);
            CHECK(rhoL > AS->rhomolar_critical());

            std::shared_ptr<AbstractState> chk(AbstractState::factory(backend, std::vector<std::string>{name}));
            chk->specify_phase(iphase_liquid);
            chk->update(DmolarT_INPUTS, rhoL, T);
            CHECK_THAT(chk->p(), Catch::Matchers::WithinRel(p_sat, 1e-7));
        }
    }
    CHECK(ran == 35);     // 16 of 18 GERG-2004 + 19 of 21 GERG-2008 components
    CHECK(skipped == 4);  // helium and hydrogen, in both models
}

TEST_CASE("GERG ancillaries are close to the converged saturation state", "[GERG]") {
    // The ancillary is only a guess, but a bad guess means slow or failed
    // convergence -- and this is also the test that fails if the generated
    // coefficient table is edited (the Task 11 mutation experiment). The
    // tolerance is each fit's OWN recorded max_rel_dev with headroom, not a
    // blanket number, so a fluid whose fit degrades cannot hide behind a
    // loose global bound.
    int checks = 0;
    for (auto model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const std::string backend = (model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008";
        for (const auto& name : component_names(model)) {
            CAPTURE(backend);
            CAPTURE(name);
            std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{name}));
            const AncillaryCoeffs aL = get_ancillary(model, name, "rhoL");
            const AncillaryCoeffs aV = get_ancillary(model, name, "rhoV");
            const AncillaryCoeffs ap = get_ancillary(model, name, "pV");
            const double tol_L = 2.0 * aL.max_rel_dev;
            const double tol_V = 2.0 * aV.max_rel_dev;
            const double tol_p = 2.0 * ap.max_rel_dev;
            for (double frac : {0.35, 0.5, 0.65, 0.8, 0.95}) {
                const double T = aL.Tmin + frac * (AS->T_critical() - aL.Tmin);
                if (!saturation_reachable(*AS, T)) continue;
                CAPTURE(T);
                REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.0, T));
                const double rhoL = AS->rhomolar(), p_sat = AS->p();
                AS->update(QT_INPUTS, 1.0, T);
                const double rhoV = AS->rhomolar();
                CHECK_THAT(evaluate_ancillary(aL, T), Catch::Matchers::WithinRel(rhoL, tol_L));
                CHECK_THAT(evaluate_ancillary(aV, T), Catch::Matchers::WithinRel(rhoV, tol_V));
                CHECK_THAT(evaluate_ancillary(ap, T), Catch::Matchers::WithinRel(p_sat, tol_p));
                ++checks;
            }
        }
    }
    // Pins the sweep's own coverage: if `saturation_reachable` starts skipping
    // everything, this fails rather than the suite passing with zero work.
    CHECK(checks >= 150);
}

TEST_CASE("GERG saturation end state agrees with the traced VLE point", "[GERG]") {
    // sat_min_liquid/sat_min_vapor were traced with teqp, not guessed. Solving
    // GERG's own VLE at that temperature must reproduce them -- which checks
    // the traced values, the ancillary fit at the very bottom of its range
    // (its worst-conditioned end), and the units of the stored p all at once.
    // The `continue` below MUST be accounted for. Three of these six fluids --
    // methane (traced end state 57.169 K), nitrogen (37.858 K) and argon
    // (45.237 K) -- sit BELOW the enforced Tmin of 60 K, so they skip today
    // and assert nothing. The three that do run are ethane (91.60 K),
    // n-butane (127.73 K) and water (239.96 K).
    //
    // Without this bookkeeping, raising the enforced lower temperature limit
    // would silently shrink the sweep while the TEST_CASE still reported
    // green -- and this is the ONLY GERG saturation test with an external
    // teqp-traced reference, so it is the worst one to lose quietly. Above
    // 91.60 K ethane drops out; above 239.96 K all six do and the case would
    // verify nothing at all.
    //
    // Names, not just a count: an exact count alone still passes if one
    // reachable fluid is swapped for one unreachable one, which is precisely
    // the "changes in either direction" this pin is supposed to catch.
    std::vector<std::string> ran_names;
    for (const auto& name : {"Methane", "Nitrogen", "Ethane", "n-Butane", "Water", "Argon"}) {
        CAPTURE(name);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        const SatEndState end = get_sat_min_state(GERGModel::GERG_2008, resolve_component(GERGModel::GERG_2008, name));
        // A hair above the stored T: calc_Tmin_sat is an inclusive lower bound
        // and QT_flash subtracts only 1e-13 from it.
        const double T = end.T_K * (1 + 1e-9);
        if (!saturation_reachable(*AS, T)) continue;
        ran_names.emplace_back(name);
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.0, T));
        CHECK_THAT(AS->rhomolar(), Catch::Matchers::WithinRel(end.rhoL_molm3, 1e-6));
        // 1e-5 on p, not 1e-6: at the bottom of the fitted range the
        // saturation pressure is sub-Pascal for the heavier alkanes
        // (n-butane's end state is 0.1515 Pa), and Maxwell's convergence
        // criterion is on the density residual, so p carries a bit more
        // slack than rhoL/rhoV do. Still three orders tighter than the
        // ancillary's own accuracy, so it cannot be satisfied by a seed.
        CHECK_THAT(AS->p(), Catch::Matchers::WithinRel(end.p_Pa, 1e-5));
        AS->update(QT_INPUTS, 1.0, T);
        CHECK_THAT(AS->rhomolar(), Catch::Matchers::WithinRel(end.rhoV_molm3, 1e-6));
    }
    // Exact set, not a count: see the note at the top of this case.
    REQUIRE(ran_names == std::vector<std::string>{"Ethane", "n-Butane", "Water"});
}

TEST_CASE("GERG saturation is independent of the ancillary seed", "[GERG]") {
    // THE PROPERTY THAT PROVES ANCILLARIES ARE SEEDS AND NOT ANSWERS.
    //
    // Perturb all three ancillaries by 1% -- roughly three times the worst
    // fit deviation in the shipped table, i.e. comfortably outside the band
    // the fit claims -- and the converged saturation state must not move at
    // all beyond solver tolerance. If it moves by anything resembling the
    // perturbation, the flash is returning the seed (or stopping one step
    // short of convergence), and every ancillary-accuracy assertion above is
    // really measuring the ancillary rather than GERG.
    for (const auto& name : {"Methane", "Propane", "CarbonDioxide", "n-Octane"}) {
        CAPTURE(name);
        const std::string gerg = resolve_component(GERGModel::GERG_2008, name);

        std::shared_ptr<AbstractState> ref(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        const double T = 0.75 * ref->T_critical();
        REQUIRE(saturation_reachable(*ref, T));
        ref->update(QT_INPUTS, 0.0, T);
        const double p0 = ref->p(), rhoL0 = ref->rhomolar();
        ref->update(QT_INPUTS, 1.0, T);
        const double rhoV0 = ref->rhomolar();

        for (double rel : {-0.01, +0.01}) {
            CAPTURE(rel);
            std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
            auto* heos = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
            REQUIRE(heos != nullptr);
            CoolProp::CoolPropFluid& fluid = heos->get_components()[0];
            fluid.ancillaries.rhoL = perturbed_ancillary(get_ancillary(GERGModel::GERG_2008, gerg, "rhoL"), rel);
            fluid.ancillaries.rhoV = perturbed_ancillary(get_ancillary(GERGModel::GERG_2008, gerg, "rhoV"), rel);
            fluid.ancillaries.pV = perturbed_ancillary(get_ancillary(GERGModel::GERG_2008, gerg, "pV"), rel);
            fluid.ancillaries.pL = fluid.ancillaries.pV;
            // Confirm the perturbation actually took effect -- otherwise this
            // test would pass by doing nothing at all.
            REQUIRE_THAT(fluid.ancillaries.rhoL.evaluate(T),
                         Catch::Matchers::WithinRel(evaluate_ancillary(get_ancillary(GERGModel::GERG_2008, gerg, "rhoL"), T) * (1 + rel), 1e-12));

            AS->update(QT_INPUTS, 0.0, T);
            CHECK_THAT(AS->p(), Catch::Matchers::WithinRel(p0, 1e-9));
            CHECK_THAT(AS->rhomolar(), Catch::Matchers::WithinRel(rhoL0, 1e-9));
            AS->update(QT_INPUTS, 1.0, T);
            CHECK_THAT(AS->rhomolar(), Catch::Matchers::WithinRel(rhoV0, 1e-9));
        }
    }
}

TEST_CASE("GERG saturation keeps GERG-typed linked states", "[GERG]") {
    // saturation_T_pure_Maxwell does all its work through HEOS.SatL/SatV. If
    // either were base-typed, calc_gas_constant would return CODATA R
    // (8.31446261815324) instead of R_GERG (8.314472) -- a 1.1e-6 relative
    // error that would land in every saturation pressure with nothing to
    // indicate it.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    auto* gerg = dynamic_cast<GERGMixtureBackend*>(AS.get());
    REQUIRE(gerg != nullptr);
    AS->update(QT_INPUTS, 0.0, 0.8 * AS->T_critical());
    for (auto* sat : {&gerg->get_SatL(), &gerg->get_SatV()}) {
        CHECK(sat->backend_name() == "GERG2008Backend");
        CHECK_THAT(sat->gas_constant(), Catch::Matchers::WithinRel(8.314472, 1e-15));
        CHECK(dynamic_cast<GERGMixtureBackend*>(sat) != nullptr);
    }
    // Not CODATA: the whole point of the check above.
    CHECK(std::abs(8.314472 - CoolProp::get_config_double(R_U_CODATA)) > 1e-6);
}

TEST_CASE("GERG handles a full 21-component list with mostly-zero mole fractions", "[GERG]") {
    // THE BACKEND'S HEADLINE USE CASE, and until now the one thing it could
    // not do. A natural-gas analysis is normally handed over as the full
    // GERG-2008 component list with a mole fraction for each, most of them
    // exactly zero. Every such composition used to return a silent NaN --
    // GERG2008ReducingFunction::f_Y_ij is x_i*x_j*(x_i+x_j)/(beta^2 x_i + x_j),
    // which is 0/0 as soon as TWO mole fractions are exactly zero -- so
    // T_reducing, rhomolar_reducing and hence every property came back NaN
    // with no error raised. Pre-existing for every multi-fluid backend
    // (GitHub #1677 / bd CoolProp-8psx); Task 9 routed around it by trimming
    // compositions to their nonzero components, which left this case untested.
    //
    // The fix is the removable-singularity guard in ReducingFunctions.cpp. The
    // assertion that matters is not merely "finite" but "identical to the
    // trimmed composition": a guard returning some other constant would still
    // be finite and still be wrong.
    const std::vector<std::string>& all = component_names(GERGModel::GERG_2008);
    REQUIRE(all.size() == 21);

    // A realistic pipeline-gas analysis over the full component list: 8
    // nonzero, 13 exactly zero. The zeros are what this test is about.
    std::map<std::string, double> analysis = {{"methane", 0.90644}, {"nitrogen", 0.03134},  {"carbondioxide", 0.00466}, {"ethane", 0.04528},
                                              {"propane", 0.00828}, {"isobutane", 0.00156}, {"n-butane", 0.00104},      {"isopentane", 0.00140}};
    std::vector<CoolPropDbl> z_full;
    std::vector<std::string> names_trim;
    std::vector<CoolPropDbl> z_trim;
    double total = 0;
    for (const auto& name : all) {
        auto it = analysis.find(name);
        const double zi = (it == analysis.end()) ? 0.0 : it->second;
        z_full.push_back(zi);
        total += zi;
        if (zi != 0.0) {
            names_trim.push_back(name);
            z_trim.push_back(zi);
        }
    }
    REQUIRE_THAT(total, Catch::Matchers::WithinRel(1.0, 1e-12));
    REQUIRE(names_trim.size() == 8);
    // 13 exactly-zero entries -> C(13,2) = 78 both-zero pairs reach f_Y_ij.
    REQUIRE(z_full.size() - names_trim.size() == 13);

    std::shared_ptr<AbstractState> full(AbstractState::factory("GERG2008", all));
    full->set_mole_fractions(z_full);
    std::shared_ptr<AbstractState> trim(AbstractState::factory("GERG2008", names_trim));
    trim->set_mole_fractions(z_trim);

    // The reducing state is where the NaN was born.
    CHECK(ValidNumber(full->T_reducing()));
    CHECK(ValidNumber(full->rhomolar_reducing()));
    CHECK_THAT(full->T_reducing(), Catch::Matchers::WithinRel(trim->T_reducing(), 1e-14));
    CHECK_THAT(full->rhomolar_reducing(), Catch::Matchers::WithinRel(trim->rhomolar_reducing(), 1e-14));

    // ...and every property follows. Note alphar is NOT trivially equal:
    // the residual sum runs over all 21 components and all 210 pairs in the
    // full case versus 8 and 28 in the trimmed one, so agreeing to 1e-14 is a
    // real statement that the zero-mole-fraction terms contribute exactly zero.
    for (auto* AS : {full.get(), trim.get()}) {
        AS->specify_phase(iphase_gas);
        AS->update(DmolarT_INPUTS, 4000.0, 300.0);
    }
    CHECK(ValidNumber(full->p()));
    CHECK_THAT(full->alphar(), Catch::Matchers::WithinRel(trim->alphar(), 1e-14));
    CHECK_THAT(full->alpha0(), Catch::Matchers::WithinRel(trim->alpha0(), 1e-14));
    CHECK_THAT(full->p(), Catch::Matchers::WithinRel(trim->p(), 1e-13));
    CHECK_THAT(full->molar_mass(), Catch::Matchers::WithinRel(trim->molar_mass(), 1e-14));
    CHECK_THAT(full->cvmolar(), Catch::Matchers::WithinRel(trim->cvmolar(), 1e-12));

    // A single zero mole fraction was never affected (f_Y_ij's denominator is
    // nonzero then), and must stay exact -- the guard fires on the both-zero
    // corner only, so infinite-dilution behaviour is untouched.
    std::shared_ptr<AbstractState> one_zero(AbstractState::factory("GERG2008", std::vector<std::string>{"methane", "ethane"}));
    one_zero->set_mole_fractions(std::vector<CoolPropDbl>{1.0, 0.0});
    std::shared_ptr<AbstractState> pure(AbstractState::factory("GERG2008", std::vector<std::string>{"methane"}));
    for (auto* AS : {one_zero.get(), pure.get()}) {
        AS->specify_phase(iphase_gas);
        AS->update(DmolarT_INPUTS, 4000.0, 300.0);
    }
    CHECK_THAT(one_zero->p(), Catch::Matchers::WithinRel(pure->p(), 1e-13));
}

TEST_CASE("GERG saturation round-trips through pressure and quality inputs", "[GERG]") {
    // PQ_flash and DQ_flash reach the saturation solver by different routes
    // than QT_flash (via the pV ancillary's invert() and a Brent sweep
    // respectively), so both are exercised.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    const double T = 0.8 * AS->T_critical();
    AS->update(QT_INPUTS, 0.0, T);
    const double p_sat = AS->p(), rhoL = AS->rhomolar();

    std::shared_ptr<AbstractState> pq(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    REQUIRE_NOTHROW(pq->update(PQ_INPUTS, p_sat, 0.0));
    CHECK_THAT(pq->T(), Catch::Matchers::WithinRel(T, 1e-8));
    CHECK_THAT(pq->rhomolar(), Catch::Matchers::WithinRel(rhoL, 1e-7));

    std::shared_ptr<AbstractState> dq(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    REQUIRE_NOTHROW(dq->update(DmolarQ_INPUTS, rhoL, 0.0));
    CHECK_THAT(dq->T(), Catch::Matchers::WithinRel(T, 1e-6));
}

TEST_CASE("GERG refuses set_reference_stateS instead of silently ignoring it", "[GERG]") {
    // set_reference_stateS dispatches on the backend prefix and, before the
    // GERG arm existed, had no final else: "GERG2008::Methane" matched
    // neither the REFPROP nor the HEOS branch, so the call returned having
    // done NOTHING -- not even validating the reference-state string.  That
    // is the worst possible outcome for this particular function, because
    // GERG's reference state (h = s = 0 for the IDEAL GAS at 298.15 K /
    // 101325 Pa) differs from every other CoolProp backend's by a large
    // offset, so a user reconciling the two would get silence and keep
    // computing in the old convention.
    for (const std::string& backend : {std::string("GERG2004"), std::string("GERG2008")}) {
        const std::string fluid = backend + "::Methane";
        CHECK_THROWS_AS(CoolProp::set_reference_stateS(fluid, "NBP"), CoolProp::NotImplementedError);
        // Including for a string that is not a valid reference state at all:
        // the throw must not depend on the argument being recognisable.
        CHECK_THROWS_AS(CoolProp::set_reference_stateS(fluid, "NOT_A_REFERENCE_STATE"), CoolProp::NotImplementedError);
    }
    // Both accepted spellings.  DataStructures.cpp registers a family name
    // ("GERG2008") AND a backend name ("GERG2008Backend") for each backend and
    // AbstractState::factory takes either, so a guard written as a literal
    // string comparison against the family name lets "GERG2008Backend::"
    // through -- which is exactly what the first version of this guard did.
    CHECK_THROWS_AS(CoolProp::set_reference_stateS("GERG2008Backend::Methane", "NBP"), CoolProp::NotImplementedError);
    CHECK_THROWS_AS(CoolProp::set_reference_stateS("GERG2004Backend::Methane", "NBP"), CoolProp::NotImplementedError);
    // The "?<options>" factory suffix too.  extract_backend splits on "::"
    // BEFORE the suffix is stripped, so the raw token reaching the guard still
    // carries it -- "GERG2008?::Methane" is a working factory string that
    // resolved to no family at all and slipped past the first version of it.
    CHECK_THROWS_AS(CoolProp::set_reference_stateS("GERG2008?::Methane", "NBP"), CoolProp::NotImplementedError);
    // ...and a composed tabular string, whose SECOND half is the GERG family.
    CHECK_THROWS_AS(CoolProp::set_reference_stateS("BICUBIC&GERG2008::Methane", "NBP"), CoolProp::NotImplementedError);
    // Non-GERG backends must NOT start throwing: a guard that over-matches
    // would break set_reference_stateS for everyone else.  Note these two are
    // currently silent no-ops for a SEPARATE pre-existing reason (they match
    // no arm of the dispatch chain at all -- bd CoolProp-mh1q); the assertion
    // here is only that the GERG guard does not claim them.  If mh1q is fixed
    // so unrecognised prefixes throw, update these two lines rather than the
    // GERG guard.
    CHECK_NOTHROW(CoolProp::set_reference_stateS("SRK::Propane", "NBP"));
    CHECK_NOTHROW(CoolProp::set_reference_stateS("BICUBIC&HEOS::Methane", "NBP"));

    // ...and the throw is scoped to GERG: HEOS is untouched.  RESET is the
    // restore idiom, chosen because it zeroes the offset rather than setting
    // one.  It is NOT a complete no-op -- it flips the fluid entry's
    // EnthalpyEntropyOffset.enabled flag and rewrites reduce.hmolar by about
    // 1.6e-4 relative -- but h and s at any state are bit-identical
    // afterwards, which is what the rest of the suite depends on.
    CHECK_NOTHROW(CoolProp::set_reference_stateS("HEOS::Methane", "RESET"));
    CHECK_THROWS_AS(CoolProp::set_reference_stateS("HEOS::Methane", "NOT_A_REFERENCE_STATE"), CoolProp::ValueError);
}

TEST_CASE("GERG refuses set_reference_stateD as well as the S variant", "[GERG]") {
    // The D variant is the sibling of the function above and ends in the same
    // set_fluid_enthalpy_entropy_offset call, so it carries the same hazard --
    // but it was left unguarded in the first round.  Its pre-guard behaviour
    // was arguably worse than the S variant's, because it failed in two
    // DIFFERENT unhelpful ways depending on spelling:
    // set_reference_stateD("GERG2008::Methane", ...) threw
    // "key [GERG2008::Methane] was not found in string_to_index_map in
    // JSONFluidLibrary" -- an internal lookup error naming neither GERG nor
    // reference states -- while set_reference_stateD("Methane", ...) succeeded
    // and adjusted HEOS with no effect on any GERG state.
    for (const char* fluid : {"GERG2004::Methane", "GERG2008::Methane", "GERG2008Backend::Methane", "GERG2004Backend::Methane", "GERG2008?::Methane",
                              "BICUBIC&GERG2008::Methane"}) {
        CAPTURE(fluid);
        CHECK_THROWS_AS(CoolProp::set_reference_stateD(fluid, 300.0, 1000.0, 0.0, 0.0), CoolProp::NotImplementedError);
    }
    // Scoped to GERG: an ordinary HEOS call still works, and is restored with
    // RESET exactly as in the S-variant case above.  Note the BARE name: unlike
    // set_reference_stateS, set_reference_stateD passes FluidName through to
    // HelmholtzEOSMixtureBackend WITHOUT splitting off the backend prefix, so
    // even "HEOS::Methane" fails here with "key [HEOS::Methane] was not found
    // in string_to_index_map in JSONFluidLibrary".  That is pre-existing and
    // deliberately not touched -- it is the non-GERG half of bd CoolProp-mh1q
    // -- but it is why the GERG guard has to do its own extract_backend rather
    // than inheriting a split the function never performs.
    CHECK_THROWS(CoolProp::set_reference_stateD("HEOS::Methane", 300.0, 1000.0, 0.0, 0.0));
    CHECK_NOTHROW(CoolProp::set_reference_stateD("Methane", 300.0, 1000.0, 0.0, 0.0));
    CHECK_NOTHROW(CoolProp::set_reference_stateS("HEOS::Methane", "RESET"));
}

TEST_CASE("GERG rejects a non-finite pure-info row", "[GERG]") {
    // get_acentric_factor has an explicit ValidNumber drift-guard on the
    // reasoning that a hand-edit made after the generator ran cannot be
    // policed by the generator.  pure_info_2004()/pure_info_2008_overrides()
    // are an equally hand-transcribed table and had no such guard.
    //
    // The specific reason it matters here: make_gerg_fluid's next use of
    // Tc_K is `EOS.limits.Tmin = std::min(60.0, info.Tc_K)`, and std::min
    // ABSORBS a NaN -- it returns 60.0, a perfectly ordinary-looking limit on
    // a fluid whose reducing temperature is NaN.  Asserted directly on the
    // helper, since the shipped tables are (and must remain) all finite, so
    // there is no fluid name that reaches it.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_NOTHROW(require_finite_pure_info("methane", PureInfo{10.139342719, 190.564, 16.04246}));
    CHECK_THROWS_AS(require_finite_pure_info("methane", PureInfo{nan, 190.564, 16.04246}), CoolProp::ValueError);
    CHECK_THROWS_AS(require_finite_pure_info("methane", PureInfo{10.139342719, nan, 16.04246}), CoolProp::ValueError);
    CHECK_THROWS_AS(require_finite_pure_info("methane", PureInfo{10.139342719, 190.564, nan}), CoolProp::ValueError);
    // ...and infinities, not just NaN.
    const double inf = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(require_finite_pure_info("methane", PureInfo{10.139342719, inf, 16.04246}), CoolProp::ValueError);
    // The std::min-absorbs-NaN mechanism this guard is protecting, pinned so
    // the comment above cannot rot into a false claim.
    CHECK(std::min(60.0, nan) == 60.0);
    // Every shipped row passes, both models.
    for (const auto& model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const auto& name : component_names(model)) {
            CAPTURE(name);
            CHECK_NOTHROW(get_pure_info(model, name));
        }
    }
}

TEST_CASE("GERG derives its acentric factor from its own equation", "[GERG]") {
    // This case used to assert the opposite: acentric_factor() threw
    // NotImplementedError, because GERG publishes no acentric factor and
    // make_gerg_fluid left EquationOfState::acentric at the _HUGE sentinel.
    // Throwing was the right answer for a sentinel; it is the wrong answer now
    // that the backend HAS the quantity.  omega is computed from GERG's own
    // equation by Pitzer's defining relation (GERGAcentric.h,
    // dev/gerg/compute_acentric.py) rather than borrowed from CoolProp's
    // per-fluid library, which belongs to a different equation of state -- so
    // nothing is borrowed and there is nothing left to refuse.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    const double omega = AS->acentric_factor();
    CHECK(ValidNumber(omega));
    CHECK_THAT(omega, Catch::Matchers::WithinRel(get_acentric_factor(GERGModel::GERG_2008, "methane"), 1e-15));
    // The keyed/trivial-output route reaches the same value, not a sentinel.
    CHECK_THAT(AS->trivial_keyed_output(iacentric_factor), Catch::Matchers::WithinRel(omega, 1e-15));
    // ...and so does PropsSI, which used to answer _HUGE (== +inf) here
    // because it swallows the exception.
    CHECK_THAT(CoolProp::PropsSI("acentric", "T", 0.0, "P", 0.0, "GERG2008::Methane"), Catch::Matchers::WithinRel(omega, 1e-15));

    // The two models disagree for the two components whose reducing
    // parameters GERG-2008 moved, and the backend must not answer one from
    // the other's table.  Carbon monoxide and isopentane are the whole reason
    // the table has 23 rows rather than 21.
    for (const char* name : {"CarbonMonoxide", "Isopentane"}) {
        std::shared_ptr<AbstractState> a4(AbstractState::factory("GERG2004", std::vector<std::string>{name}));
        std::shared_ptr<AbstractState> a8(AbstractState::factory("GERG2008", std::vector<std::string>{name}));
        INFO("component := " << name);
        CHECK(a4->acentric_factor() != a8->acentric_factor());
    }

    // A MIXTURE still has no acentric factor -- but that is CoolProp's own
    // rule for every backend (HelmholtzEOSMixtureBackend::calc_acentric_factor
    // throws ValueError for a mixture), not a GERG restriction, so the GERG
    // override was removed rather than narrowed.  Asserting the same shape on
    // HEOS is what documents that this is shared behaviour.
    for (const std::string& backend : {std::string("GERG2008"), std::string("HEOS")}) {
        std::shared_ptr<AbstractState> mix(AbstractState::factory(backend, std::vector<std::string>{"Methane", "Ethane"}));
        mix->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
        INFO("backend := " << backend);
        CHECK_THROWS_AS(mix->acentric_factor(), CoolProp::ValueError);
    }

    // The accessor rejects a component the model does not have, rather than
    // falling through to the GERG-2004 base table.  n-nonane is GERG-2008-only.
    CHECK_THROWS_AS(get_acentric_factor(GERGModel::GERG_2004, "n-nonane"), CoolProp::ValueError);
    CHECK_NOTHROW(get_acentric_factor(GERGModel::GERG_2008, "n-nonane"));
}

TEST_CASE("GERG acentric factors satisfy their own definition", "[GERG]") {
    // GERGAcentric.h is a GENERATED table: dev/gerg/compute_acentric.py
    // evaluates omega = -1 - log10(p_sat(0.7 T_c)/p_c) against teqp's GERG
    // implementation and writes the scalars out.  Nothing about a committed
    // scalar is self-checking, so this case re-derives every row from THIS
    // backend -- CoolProp's own saturation solver, on CoolProp's own
    // transcription of the coefficient tables -- and checks the two agree.
    //
    // That makes it a genuine cross-implementation check rather than a
    // restatement: teqp's path is a spinodal scan plus geometric bisection on
    // p with equal chemical potential, then a damped 2-D Newton;
    // CoolProp's is SaturationSolvers::saturation_T_pure_Maxwell seeded from
    // the fitted ancillaries.  Nothing but the underlying equation is shared.
    //
    // It is also the tripwire for a stale table: edit a pure-fluid n/t/d/l row
    // or a Table A3.5 reducing parameter without rerunning compute_acentric.py
    // and this case fails, naming the fluid.
    //
    // Tolerance: the worst observed disagreement across all 23 is 8.6e-12
    // absolute in omega.  1e-9 leaves three orders of magnitude of headroom
    // for solver-tolerance noise while still catching any real change --
    // a one-digit typo in the table moves omega by ~1e-3 at least.
    const double tol = 1e-9;
    std::size_t checked = 0;
    for (const GERGModel model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        const std::string backend = (model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008";
        for (const std::string& gerg_name : component_names(model)) {
            std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{gerg_name}));
            // p_c is the pressure at GERG's TABULATED reducing state, which is
            // exactly what p_critical() reports -- reading it back off the
            // backend here (rather than recomputing it) is what pins the
            // consistency the Wilson K-factor seed depends on: it reads p_c,
            // T_c and omega off the same fluid.
            const double Tc = AS->T_critical(), pc = AS->p_critical();
            const double T = 0.7 * Tc;

            // Helium (T_c = 5.1953 K) and hydrogen (T_c = 33.19 K) have
            // EOS.limits.Tmin == T_c, so 0.7*T_c is below the range
            // check_gerg_range_of_validity lets a caller reach.  That limit is
            // the published operating envelope of the MIXTURE model, not a
            // statement that the pure equation stops working, and omega is a
            // definitional constant of the equation rather than a state a
            // caller asked for -- so it is computed there deliberately, using
            // CoolProp's own documented opt-out rather than by weakening the
            // guard.  See dev/gerg/compute_acentric.py and Web/coolprop/GERG.rst.
            const bool below_enforced_range = T < AS->Tmin();
            // BIDIRECTIONAL, on purpose.  Asserting only "if we needed the
            // opt-out then it was helium or hydrogen" pins nothing in the
            // other direction: if a future change gave those two a reachable
            // 0.7*T_c, the opt-out would quietly stop being exercised and
            // nothing would say so.  Equality between the two sets is what
            // keeps the documented exception exactly two fluids wide.
            const bool expected_below = (gerg_name == "helium" || gerg_name == "hydrogen");
            CHECK(below_enforced_range == expected_below);
            const bool was_set = CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS);
            if (below_enforced_range) {
                CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, true);
            }
            double psat = std::numeric_limits<double>::quiet_NaN();
            try {
                std::shared_ptr<AbstractState> S(AbstractState::factory(backend, std::vector<std::string>{gerg_name}));
                S->update(QT_INPUTS, 0.0, T);
                psat = S->p();
            } catch (...) {  // NOLINT(bugprone-empty-catch) -- psat stays NaN and the CHECKs below fail, naming the fluid
            }
            // Restore before asserting: a failing CHECK must not leave the
            // global config flipped for every later test case.  Restoring the
            // SAVED value rather than hard-coding false is what keeps this
            // correct if the suite is ever run with the flag already on.
            if (below_enforced_range) {
                CoolProp::set_config_bool(DONT_CHECK_PROPERTY_LIMITS, was_set);
            }
            INFO(backend << " :: " << gerg_name);
            // ValidNumber first: a NaN psat makes log10 NaN, and
            // WithinAbs(NaN, tol) is false -- but so is every other matcher,
            // so without this the failure would not say WHY.
            REQUIRE(ValidNumber(psat));
            REQUIRE(psat > 0);
            const double omega_rederived = -1.0 - std::log10(psat / pc);
            CHECK_THAT(omega_rederived, Catch::Matchers::WithinAbs(get_acentric_factor(model, gerg_name), tol));
            ++checked;
        }
    }
    // 18 GERG-2004 components + 21 GERG-2008 components.  The table itself has
    // 23 distinct rows (GERG-2008 overrides carbon monoxide and isopentane);
    // this loop exercises every (model, component) pair that resolves to one.
    CHECK(checked == 39);
}

TEST_CASE("GERG acentric factors are physically sensible", "[GERG]") {
    // A COARSE sanity net, not a validation target.  These published values
    // belong to each fluid's REFERENCE equation of state -- exactly the
    // provenance the backend refuses to carry -- so they live here, in a test,
    // and nothing in src/Backends/GERG reads them.  They cannot agree exactly:
    // GERG's reducing point is a fitted quantity rather than the true critical
    // point, and the shortened technical form is not the reference equation.
    //
    // What this catches is a gross error -- a sign flip, a fluid's row copied
    // from its neighbour, a factor of ten -- which is the class of mistake a
    // generated table is actually exposed to.  The largest observed
    // disagreement across all 23 rows is 3.1e-3 (n-octane); the 0.01 band is
    // ~3x that, tight enough that no two fluids in the table are within it of
    // each other except by physics.
    //
    // Source: Poling, Prausnitz & O'Connell, "The Properties of Gases and
    // Liquids", 5th ed., Appendix A.
    const std::map<std::string, double> published = {
      {"methane", 0.011},         {"nitrogen", 0.037},  {"carbondioxide", 0.224},  {"ethane", 0.099},   {"propane", 0.152},   {"n-butane", 0.200},
      {"isobutane", 0.186},       {"n-pentane", 0.251}, {"isopentane", 0.227},     {"n-hexane", 0.299}, {"n-heptane", 0.349}, {"n-octane", 0.398},
      {"hydrogen", -0.219},       {"oxygen", 0.022},    {"carbonmonoxide", 0.050}, {"water", 0.344},    {"helium", -0.385},   {"argon", -0.002},
      {"hydrogensulfide", 0.100}, {"n-nonane", 0.443},  {"n-decane", 0.489}};
    CHECK(published.size() == 21);
    for (const GERGModel model : {GERGModel::GERG_2004, GERGModel::GERG_2008}) {
        for (const std::string& gerg_name : component_names(model)) {
            auto it = published.find(gerg_name);
            REQUIRE(it != published.end());
            INFO(((model == GERGModel::GERG_2004) ? "GERG2004" : "GERG2008") << " :: " << gerg_name);
            CHECK_THAT(get_acentric_factor(model, gerg_name), Catch::Matchers::WithinAbs(it->second, 0.01));
        }
    }
}

TEST_CASE("GERG mixture VLE, phase envelopes and gas-phase DmolarP work", "[GERG]") {
    // Everything in this case failed before make_gerg_fluid was given a real
    // acentric factor: CoolProp's Wilson K-factor seed and
    // FlashRoutines::T_DP_PengRobinson read EquationOfState::acentric
    // DIRECTLY, so the _HUGE sentinel produced NaN initial guesses that
    // surfaced several solver frames later as e.g. "solver_rho_Tp was unable
    // to find a solution for T=150, p=0, with guess value nan".
    //
    // The numbers below are NOT reference values -- they are this backend's
    // own converged answers, recorded so that a regression in the seeding
    // shows up as a failure here rather than as a silently different flash.
    // What is asserted about them independently is thermodynamic consistency:
    // the bubble and dew states at the same T bracket the pressure, and a
    // round trip through PQ returns the temperature it started from.
    auto mix = [](const char* backend) {
        std::shared_ptr<AbstractState> p(AbstractState::factory(backend, std::vector<std::string>{"Methane", "Ethane"}));
        p->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
        return p;
    };

    // QT, the exact case named in the bug report.
    std::shared_ptr<AbstractState> bub = mix("GERG2008");
    REQUIRE_NOTHROW(bub->update(QT_INPUTS, 0.0, 150.0));
    std::shared_ptr<AbstractState> dew = mix("GERG2008");
    REQUIRE_NOTHROW(dew->update(QT_INPUTS, 1.0, 150.0));
    CHECK(ValidNumber(bub->p()));
    CHECK(ValidNumber(dew->p()));
    // Bubble pressure above dew pressure at fixed T and composition.
    CHECK(bub->p() > dew->p());
    CHECK_THAT(bub->p(), Catch::Matchers::WithinRel(924079.0, 1e-4));
    CHECK_THAT(dew->p(), Catch::Matchers::WithinRel(93875.4, 1e-4));

    // PQ, and a round trip: the bubble temperature at the bubble pressure is
    // the temperature we started from.
    std::shared_ptr<AbstractState> pq = mix("GERG2008");
    REQUIRE_NOTHROW(pq->update(PQ_INPUTS, bub->p(), 0.0));
    CHECK_THAT(pq->T(), Catch::Matchers::WithinRel(150.0, 1e-6));

    // GERG-2004 reaches the same paths.
    std::shared_ptr<AbstractState> bub04 = mix("GERG2004");
    REQUIRE_NOTHROW(bub04->update(QT_INPUTS, 0.0, 150.0));
    CHECK(ValidNumber(bub04->p()));

    // build_phase_envelope, which the Wilson seed also gates.
    std::shared_ptr<AbstractState> env = mix("GERG2008");
    REQUIRE_NOTHROW(env->build_phase_envelope(""));
    const PhaseEnvelopeData& ped = env->get_phase_envelope_data();
    CHECK(ped.T.size() > 100);
    CHECK(ped.T.size() == ped.p.size());
    for (std::size_t i = 0; i < ped.T.size(); ++i) {
        INFO("envelope point " << i);
        REQUIRE(ValidNumber(ped.T[i]));
        REQUIRE(ValidNumber(ped.p[i]));
    }

    // DmolarP on a gas-phase PURE state: this went through
    // T_DP_PengRobinson, which reads the acentric factor directly.
    std::shared_ptr<AbstractState> ref(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    REQUIRE_NOTHROW(ref->update(PT_INPUTS, 1e6, 300.0));
    const double rho = ref->rhomolar();
    std::shared_ptr<AbstractState> dp(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    REQUIRE_NOTHROW(dp->update(DmolarP_INPUTS, rho, 1e6));
    CHECK_THAT(dp->T(), Catch::Matchers::WithinRel(300.0, 1e-8));

    // STILL UNAVAILABLE, and NOT a GERG limitation: CoolProp implements
    // neither DQ nor DP for ANY mixture backend.  Asserting the identical
    // failure on HEOS is what identifies it as shared code rather than
    // something the GERG backend broke -- and if CoolProp ever implements
    // them, BOTH halves start failing together and this comment is the
    // pointer to what to update.
    //
    // The exception TYPE is asserted, not just "something threw": these are
    // the NotImplementedError raised by FlashRoutines::DQ_flash /
    // FlashRoutines::DP_flash (FlashRoutines.cpp:500,568).  A bare
    // CHECK_THROWS would keep passing if the call started failing for an
    // entirely different reason -- a broken mixture, a bad composition -- and
    // the pin would silently stop meaning what it says.
    for (const char* backend : {"GERG2008", "HEOS"}) {
        INFO("backend := " << backend);
        CHECK_THROWS_AS(mix(backend)->update(DmolarQ_INPUTS, 20000.0, 0.0), CoolProp::NotImplementedError);
        std::shared_ptr<AbstractState> g = mix(backend);
        REQUIRE_NOTHROW(g->update(PT_INPUTS, 1e6, 300.0));
        const double rho_mix = g->rhomolar();
        CHECK_THROWS_AS(mix(backend)->update(DmolarP_INPUTS, rho_mix, 1e6), CoolProp::NotImplementedError);
    }

    // DQ on a PURE GERG fluid does work -- the mixture gap above is about
    // mixtures, not about the GERG backend.
    std::shared_ptr<AbstractState> pure_dq(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    CHECK_NOTHROW(pure_dq->update(DmolarQ_INPUTS, 20000.0, 0.0));
}

TEST_CASE("GERG pins which zero-mole-fraction properties work and which do not", "[GERG]") {
    // Where the both-zero guard stops, and WHY it stops there.  f_Y_ij is a
    // cubic over a linear denominator, i.e. homogeneous of degree 2, so its
    // k-th composition derivative is homogeneous of degree 2 - k.  Over
    // non-negative mole fractions the denominator is bounded away from zero on
    // the unit directions, so each derivative is bounded there, and
    // |g(t*u)| <= M * t^(2-k):
    //
    //   k = 0, 1  ->  vanishes with t uniformly in direction u.  The limit is
    //                 0 however the corner is approached, so 0 is the unique
    //                 CONTINUOUS extension and guarding is exact, not a fudge.
    //   k >= 2    ->  degree 0 or less.  Constant along each ray but different
    //                 between rays, so NO limit exists and NO value is correct.
    //
    // Hence this split, which is a property of the function rather than a
    // choice: the value and the first derivatives are repaired (so fugacity
    // works), and the second derivatives are deliberately left NaN (so
    // dlnphi_dxj, and therefore phase envelopes, still fail).  Web/coolprop/
    // GERG.rst states the same split; this test is what stops the code and the
    // documentation drifting apart in either direction.
    //
    // The second-derivative assertion below must NOT be "fixed" by returning
    // some value: doing so would replace an honest NaN with a silently
    // direction-dependent number.  See both_mole_fractions_zero in
    // ReducingFunctions.cpp.
    const std::vector<std::string> names = {"Methane", "Nitrogen", "Ethane", "Propane"};
    const std::vector<CoolPropDbl> z_full = {0.9, 0.1, 0.0, 0.0};  // two exact zeros
    const std::vector<std::string> names_trim = {"Methane", "Nitrogen"};
    const std::vector<CoolPropDbl> z_trim = {0.9, 0.1};

    // Both backends, because this is a shared-code behaviour and not a GERG
    // peculiarity: asserting it on HEOS too is what documents that the branch
    // neither introduced nor GERG-scoped it.
    for (const std::string& backend : {std::string("GERG2008"), std::string("HEOS")}) {
        std::shared_ptr<AbstractState> full(AbstractState::factory(backend, names));
        full->set_mole_fractions(z_full);
        std::shared_ptr<AbstractState> trim(AbstractState::factory(backend, names_trim));
        trim->set_mole_fractions(z_trim);
        for (AbstractState* AS : {full.get(), trim.get()}) {
            AS->update(PT_INPUTS, 1e6, 300.0);
        }
        INFO("backend := " << backend);
        // WORKS: everything that flows from the reducing state.
        CHECK(ValidNumber(full->rhomolar()));
        CHECK_THAT(full->rhomolar(), Catch::Matchers::WithinRel(trim->rhomolar(), 1e-12));
        CHECK_THAT(full->alphar(), Catch::Matchers::WithinRel(trim->alphar(), 1e-12));
        CHECK_THAT(full->cvmolar(), Catch::Matchers::WithinRel(trim->cvmolar(), 1e-11));
        // WORKS NOW: first composition derivatives, and therefore fugacity.
        // Checked against the trimmed composition rather than merely for
        // finiteness -- a finite but wrong number would pass ValidNumber.
        CHECK(ValidNumber(full->fugacity_coefficient(0)));
        CHECK_THAT(full->fugacity_coefficient(0), Catch::Matchers::WithinRel(trim->fugacity_coefficient(0), 1e-12));
        CHECK_THAT(full->fugacity_coefficient(1), Catch::Matchers::WithinRel(trim->fugacity_coefficient(1), 1e-12));
        // The trimmed composition, with no exact zeros, is fine -- which is
        // what identifies the zeros as the cause rather than the mixture.
        CHECK(ValidNumber(trim->fugacity_coefficient(0)));

        // STILL DOES NOT WORK, and cannot: the SECOND composition derivative.
        // Reached directly rather than through a phase envelope so the pin is
        // unambiguous -- an envelope can fail for unrelated reasons, which
        // would make this assertion pass for the wrong cause.
        auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(full.get());
        REQUIRE(heos != nullptr);
        CHECK(ValidNumber(heos->Reducing->dTrdxi__constxj(z_full, 0, XN_DEPENDENT)));
        CHECK(!ValidNumber(heos->Reducing->d2Trdxidxj(z_full, 0, 1, XN_DEPENDENT)));
    }
}

TEST_CASE("GERG pins which pure-fluid input pairs work and which do not", "[GERG]") {
    // Web/coolprop/GERG.rst has a section "Which input pairs actually work"
    // asserting that HmolarP/PSmolar/PUmolar do NOT work on a PURE GERG fluid.
    // Nothing mechanically links documentation to code, so this case is that
    // link: it fails when the behaviour changes and points at the doc section
    // to update.
    //
    // THE ORIGINAL DIAGNOSIS WAS WRONG AND WAS CORRECTED HERE, which is worth
    // recording because it is the reason this case did NOT start failing when
    // the acentric factor was fixed. It used to say the cause was the _HUGE
    // acentric sentinel collapsing solver_rho_Tp_SRK's guess to the ideal-gas
    // density. The sentinel is gone -- make_gerg_fluid now writes a
    // GERG-derived omega (GERGAcentric.h) and the SRK seed is good: at the
    // failing point (T = 60 K, p = 1 MPa) it now guesses 30499 mol/m^3 against
    // a true 30581 mol/m^3, and update(PT_INPUTS, 1e6, 60) succeeds. These
    // three flashes still fail, for two OTHER reasons, neither of which is the
    // acentric factor and neither of which is in that fix's scope:
    //
    //  1. EOS.ptriple is still the _HUGE sentinel (GERG publishes no triple
    //     point). FlashRoutines' gas-phase branch tests `p < p_triple()` to
    //     decide whether a saturation temperature exists; with p_triple()
    //     infinite that test is always true, so the T bracket's lower end
    //     falls back to Tmin() = 60 K instead of the saturation temperature
    //     (149.1 K at 1 MPa). At 60 K / 1 MPa methane is a compressed liquid
    //     while the solver has been told the phase is gas, and the density
    //     solve diverges to a negative root.
    //  2. Behind that sits a second, independent block: the bracket's UPPER
    //     end is 1.5*Tmax() = 1050 K, and GERGMixtureBackend::update()'s own
    //     range guard rejects any T above 700 K -- including the internal
    //     update() the bracketing solver makes. Verified by experiment:
    //     giving EOS.ptriple a finite value moves the failure from (1) to (2).
    //
    // Both are tracked as bd CoolProp-kr4o.  If either is
    // fixed, the CHECK_THROWS below start failing -- that is the signal to
    // update GERG.rst's "Which input pairs actually work", not to relax them.
    std::shared_ptr<AbstractState> ref(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"}));
    const double p = 1e6, T = 300.0;
    REQUIRE_NOTHROW(ref->update(PT_INPUTS, p, T));
    const double h = ref->hmolar(), s = ref->smolar(), u = ref->umolar();

    // PT is the reference point and must keep working -- otherwise the three
    // CHECK_THROWS below could "pass" simply because the fluid is broken.
    CHECK(ValidNumber(h));

    auto fresh = [] { return std::shared_ptr<AbstractState>(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane"})); };
    CHECK_THROWS(fresh()->update(HmolarP_INPUTS, h, p));
    CHECK_THROWS(fresh()->update(PSmolar_INPUTS, p, s));
    CHECK_THROWS(fresh()->update(PUmolar_INPUTS, p, u));

    // Diagnosis (1) above, asserted rather than asserted-in-a-comment: the
    // acentric factor is a real number now, and the state the density solve
    // fails at is reachable directly through PT. If this ever stops holding,
    // the comment above is stale.
    CHECK(ValidNumber(ref->acentric_factor()));
    CHECK(!ValidNumber(ref->p_triple()));
    std::shared_ptr<AbstractState> cold = fresh();
    REQUIRE_NOTHROW(cold->update(PT_INPUTS, 1e6, 60.0));
    CHECK(cold->rhomolar() > 30000.0);

    // DmolarP used to fail here too, by the separate T_DP_PengRobinson route.
    // The acentric fix DID unblock that one -- pinned so it cannot silently
    // regress back into the list above.
    std::shared_ptr<AbstractState> dp = fresh();
    REQUIRE_NOTHROW(dp->update(DmolarP_INPUTS, ref->rhomolar(), p));
    CHECK_THAT(dp->T(), Catch::Matchers::WithinRel(T, 1e-8));

    // Through PropsSI the same failure degrades to _HUGE plus an errstring
    // rather than an exception, which is the user-visible shape and the reason
    // the docs call it out. Pinned so that a change in THAT behaviour is also
    // caught.
    const double h_hi = CoolProp::PropsSI("Hmolar", "P", p, "T", T, "GERG2008::Methane");
    CHECK(ValidNumber(h_hi));
    CHECK_FALSE(ValidNumber(CoolProp::PropsSI("T", "P", p, "Hmolar", h_hi, "GERG2008::Methane")));

    // The same flash on a GERG MIXTURE succeeds -- this is specific to pure
    // fluids, and the doc says so. If this ever starts throwing, the doc's
    // "pure fluids only" framing is wrong.
    std::shared_ptr<AbstractState> mix(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "Ethane"}));
    mix->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    REQUIRE_NOTHROW(mix->update(PT_INPUTS, p, T));
    const double h_mix = mix->hmolar();
    std::shared_ptr<AbstractState> mix2(AbstractState::factory("GERG2008", std::vector<std::string>{"Methane", "Ethane"}));
    mix2->set_mole_fractions(std::vector<CoolPropDbl>{0.9, 0.1});
    CHECK_NOTHROW(mix2->update(HmolarP_INPUTS, h_mix, p));
    CHECK_THAT(mix2->T(), Catch::Matchers::WithinRel(T, 1e-6));
}

TEST_CASE("GERG tabular wrappers build a table and then reject every lookup", "[GERG]") {
    // Web/coolprop/GERG.rst has a "Tabular backends wrapping GERG" section
    // stating that BICUBIC&GERG2008 / TTSE&GERG2008 are not supported, that
    // they do NOT fail loudly, and that they persist a cache directory before
    // rejecting every lookup.  That was the only claim on the page with no
    // code or test behind it; this case is the link, in the same spirit as
    // "GERG pins which pure-fluid input pairs work and which do not" above.
    //
    // ALTERNATIVE_TABLES_DIRECTORY is redirected first, and restored after, so
    // running the suite does not leave a stale GERG table cache in the
    // developer's ~/.CoolProp/Tables -- which is exactly what happens without
    // the redirect, and exactly what the documentation tells readers to go
    // delete.  Note the TRAILING SLASH: TabularBackend::path_to_tables
    // concatenates this string with the backend name directly, so omitting it
    // silently produces a sibling "<dir>GERG2008Backend(...)" path instead of
    // a subdirectory.
    const std::string prev = CoolProp::get_config_string(ALTERNATIVE_TABLES_DIRECTORY);
    const std::string dir = (std::filesystem::temp_directory_path() / "CoolProp-GERG-tabular-test").string() + "/";
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
    CoolProp::set_config_string(ALTERNATIVE_TABLES_DIRECTORY, dir);

    for (const char* backend : {"BICUBIC&GERG2008", "TTSE&GERG2008"}) {
        CAPTURE(backend);
        // The factory itself succeeds -- this is the "does not fail loudly" half.
        std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, std::vector<std::string>{"Methane"}));
        REQUIRE(AS != nullptr);
        // ...and then a lookup at a perfectly ordinary state well inside the
        // GERG range of validity is rejected.
        CHECK_THROWS_AS(AS->update(PT_INPUTS, 1e6, 300.0), CoolProp::ValueError);
    }
    // The documented side effect: a cache directory really is written.
    CHECK(std::filesystem::exists(dir));

    CoolProp::set_config_string(ALTERNATIVE_TABLES_DIRECTORY, prev);
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("GERG mixture properties are invariant under component reordering", "[GERG]") {
    // Ian's review check.  The XN_DEPENDENT composition-derivative formulation
    // singles out the LAST component: dYrdxi__constxj has a dedicated (i, N-1)
    // term and a loop over (k, N-1) pairs, both of which inlined the 0/0
    // expression rather than calling the guarded helpers.  So permuting the
    // component list is the sharpest available probe of those guards -- it moves
    // which component is "N-1" and therefore which pairs take the special path,
    // while the physics must not move at all.
    //
    // Non-unity beta is what gives the test teeth.  For beta_ij == 1 the
    // reducing function is symmetric in i and j term-by-term, so an ordering
    // bug can hide.  These four components are chosen for strongly non-unity
    // betas -- CO2/n-butane is betaV = 1.174761, betaT = 1.018171, and
    // get_betasgammas RECIPROCATES the betas (not the gammas) when the pair is
    // looked up in reversed order, so a reversed lookup that forgot to
    // reciprocate would shift the reducing state rather than cancel out.
    //
    // Both compositions matter and for different reasons:
    //   - no zeros  : baseline invariance, which held before the guards too.
    //   - two zeros : the case the guards address.  Before them this was
    //                 genuinely order-DEPENDENT -- whether a zero landed in the
    //                 last slot decided whether the derivative came back NaN,
    //                 which is why the documentation used to say the behaviour
    //                 depended on component order.  It must not any more.
    const std::vector<std::string> base_names = {"CarbonDioxide", "n-Butane", "Methane", "Nitrogen"};

    struct Composition
    {
        const char* label;
        std::vector<CoolPropDbl> z;
    };
    const std::vector<Composition> compositions = {
      {"no zeros", {0.10, 0.05, 0.80, 0.05}},
      {"two exact zeros", {0.10, 0.00, 0.90, 0.00}},
      // Three zeros -> C(3,2) = 3 both-zero pairs, the densest exercise of the
      // guards available in a quaternary.  The permutation loop below already
      // walks every arrangement, so this does not need its own "zero in the
      // last slot" variant: some permutation of each composition puts a zero
      // there, which is the arrangement that used to poison the whole first
      // derivative (the (k, N-1) loop runs over every component).
      //
      // Kept methane-dominant deliberately.  A CO2/n-butane-rich composition
      // is two-phase at this state and fugacity_coefficient then throws
      // "not well-defined in the two-phase region", which would make this
      // case fail for a reason that has nothing to do with reordering.
      {"three exact zeros", {0.00, 0.00, 1.00, 0.00}},
    };

    for (const std::string& backend : {std::string("GERG2008"), std::string("HEOS")}) {
        for (const Composition& comp : compositions) {
            INFO("backend := " << backend << ", composition := " << comp.label);

            std::shared_ptr<AbstractState> ref(AbstractState::factory(backend, base_names));
            ref->set_mole_fractions(comp.z);
            ref->update(PT_INPUTS, 2e6, 320.0);

            // Every one of the 4! orderings, not a hand-picked few: the special
            // (i, N-1) path only fires for whichever component is last, so a
            // sample could miss the arrangement that breaks.
            std::vector<std::size_t> perm = {0, 1, 2, 3};
            std::size_t n_perms = 0;
            do {
                std::vector<std::string> names;
                std::vector<CoolPropDbl> z;
                for (std::size_t idx : perm) {
                    names.push_back(base_names[idx]);
                    z.push_back(comp.z[idx]);
                }
                INFO("permutation := " << names[0] << "," << names[1] << "," << names[2] << "," << names[3]);

                std::shared_ptr<AbstractState> AS(AbstractState::factory(backend, names));
                AS->set_mole_fractions(z);
                AS->update(PT_INPUTS, 2e6, 320.0);

                // Reducing state first -- if this moves, everything else does.
                CHECK_THAT(AS->T_reducing(), Catch::Matchers::WithinRel(ref->T_reducing(), 1e-13));
                CHECK_THAT(AS->rhomolar_reducing(), Catch::Matchers::WithinRel(ref->rhomolar_reducing(), 1e-13));
                CHECK_THAT(AS->rhomolar(), Catch::Matchers::WithinRel(ref->rhomolar(), 1e-12));
                CHECK_THAT(AS->alphar(), Catch::Matchers::WithinRel(ref->alphar(), 1e-12));
                CHECK_THAT(AS->hmolar(), Catch::Matchers::WithinRel(ref->hmolar(), 1e-11));
                CHECK_THAT(AS->smolar(), Catch::Matchers::WithinRel(ref->smolar(), 1e-11));
                CHECK_THAT(AS->cvmolar(), Catch::Matchers::WithinRel(ref->cvmolar(), 1e-11));
                CHECK_THAT(AS->speed_sound(), Catch::Matchers::WithinRel(ref->speed_sound(), 1e-11));
                CHECK_THAT(AS->molar_mass(), Catch::Matchers::WithinRel(ref->molar_mass(), 1e-13));

                // Fugacity coefficients are the FIRST composition derivatives,
                // i.e. exactly what the newly guarded inlined terms feed.
                // Compared per component THROUGH the permutation, so a
                // consistent-but-shuffled answer cannot pass.
                for (std::size_t k = 0; k < perm.size(); ++k) {
                    if (comp.z[perm[k]] == 0.0) {
                        continue;  // fugacity of an absent component is not meaningful
                    }
                    INFO("component := " << names[k]);
                    CHECK(ValidNumber(AS->fugacity_coefficient(k)));
                    CHECK_THAT(AS->fugacity_coefficient(k), Catch::Matchers::WithinRel(ref->fugacity_coefficient(perm[k]), 1e-11));
                }
                ++n_perms;
            } while (std::next_permutation(perm.begin(), perm.end()));
            CHECK(n_perms == 24);
        }
    }
}

#endif /* ENABLE_CATCH */
