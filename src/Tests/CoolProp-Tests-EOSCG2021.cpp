// Validation tests for the EOS-CG-2021 binary mixture models added for the CCS/flue-gas
// components (Neumann et al., Int. J. Thermophys. 44, 178 (2023), doi:10.1007/s10765-023-03263-6).
//
// The two binary-specific departure functions that come with published computer-implementation
// test values are checked at full precision against their primary sources:
//   - CO2 + Ar : Loevseth et al., Fluid Phase Equilib. 466, 48 (2018), Table 15 (single phase).
//   - CO2 + CO : Souza et al., Appl. Energy 251, 113398 (2019), in-text verification point.
// The CO2 + SO2 reducing model (no departure function) has no test table, so a few bubble-point
// pressures at 263 K are checked against the paper's Fig. 12 (looser, data-scatter tolerance).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

TEST_CASE("EOS-CG-2021 CO2-Ar reproduces Loevseth (2018) Table 15", "[mixtures][EOS-CG][CCS]") {
    // Single-phase test values (T [K], rho [mol/m3], z_CO2) -> p [MPa], cv [J/mol/K], w [m/s].
    struct Pt
    {
        double T, rho, z, p_MPa, cv, w;
    };
    const std::vector<Pt> pts = {
      {273.15, 915, 0.25, 2.00274, 16.7461, 286.969},   {273.15, 4000, 0.50, 6.98643, 24.6084, 255.285},
      {273.15, 18400, 0.70, 18.0108, 33.3998, 415.969}, {323.15, 1175, 0.50, 2.99922, 22.0310, 293.881},
      {323.15, 40, 0.75, 0.107170, 25.6280, 287.368},   {323.15, 24000, 0.95, 89.4303, 39.5474, 922.065},
    };
    for (const auto& pt : pts) {
        CAPTURE(pt.T, pt.rho, pt.z);
        // CoolProp has no DmolarT flash for mixtures; T and rho fix the state directly.
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&Argon"));
        auto* be = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(be != nullptr);
        be->set_mole_fractions({pt.z, 1 - pt.z});
        be->specify_phase(iphase_gas);  // single-phase test states; label only unblocks speed_sound
        be->update_DmolarT_direct(pt.rho, pt.T);
        CHECK(be->p() / 1e6 == Catch::Approx(pt.p_MPa).epsilon(3e-5));
        CHECK(be->cvmolar() == Catch::Approx(pt.cv).epsilon(3e-5));
        CHECK(be->speed_sound() == Catch::Approx(pt.w).epsilon(3e-5));
    }
}

TEST_CASE("EOS-CG-2021 CO2-CO reproduces Souza (2019) test value", "[mixtures][EOS-CG][CCS]") {
    // x = 0.5, T = 250 K, rho = 20000 mol/m3 -> p = 32.803 MPa, w = 597.04 m/s, a^r = -1.02820177.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&CarbonMonoxide"));
    auto* be = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    REQUIRE(be != nullptr);
    be->set_mole_fractions({0.5, 0.5});
    be->specify_phase(iphase_gas);  // single-phase test state; label only unblocks speed_sound
    be->update_DmolarT_direct(20000.0, 250.0);
    CHECK(be->p() / 1e6 == Catch::Approx(32.803).epsilon(3e-5));
    CHECK(be->speed_sound() == Catch::Approx(597.04).epsilon(3e-5));
    CHECK(be->alphar() == Catch::Approx(-1.02820177).epsilon(1e-5));
}

TEST_CASE("EOS-CG-2021 CO2-SO2 bubble/dew match Fig. 12 (263 K)", "[mixtures][EOS-CG][CCS]") {
    // Bubble- and dew-point pressures digitized from Fig. 12 at 263 K.  Reference values carry
    // both figure-reading and experimental scatter, so a 5% relative tolerance is used.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&SulfurDioxide"));
    SECTION("bubble (liquid composition -> bubble pressure)") {
        struct B
        {
            double xCO2, p_MPa;
        };
        for (const auto& b : std::vector<B>{{0.40, 1.29}, {0.60, 1.70}, {0.80, 2.10}}) {
            CAPTURE(b.xCO2);
            AS->set_mole_fractions({b.xCO2, 1 - b.xCO2});
            AS->update(QT_INPUTS, 0.0, 263.15);
            CHECK(AS->p() / 1e6 == Catch::Approx(b.p_MPa).epsilon(0.05));
        }
    }
    SECTION("dew (vapor composition -> dew pressure)") {
        // Lower branch of Fig. 12; the vapour is CO2-rich at low pressure.
        struct D
        {
            double yCO2, p_MPa;
        };
        for (const auto& d : std::vector<D>{{0.60, 0.24}, {0.80, 0.52}}) {
            CAPTURE(d.yCO2);
            AS->set_mole_fractions({d.yCO2, 1 - d.yCO2});
            AS->update(QT_INPUTS, 1.0, 263.15);  // dew point
            CHECK(AS->p() / 1e6 == Catch::Approx(d.p_MPa).epsilon(0.08));
        }
    }
}

TEST_CASE("EOS-CG-2021 CO2-Ar saturated liquid/vapor reproduce Loevseth (2018) Table 16", "[mixtures][EOS-CG][CCS]") {
    // Table 16: saturated vapour and liquid at the SAME composition (not in equilibrium with each
    // other).  Each is a single-phase state pinned by (T, rho, z), so the pressure is checked with
    // update_DmolarT_direct rather than a QT flash -- the mixture bubble/dew flash does not yet
    // converge for CO2-Ar at these high liquid pressures on this backend (see the DHSU_T/HSU flash
    // work in #3148, not yet on master).  z = the phase composition; rho = that phase's density.
    struct Sat
    {
        double T, z, rho, p_MPa;
    };
    const std::vector<Sat> pts = {
      {223.15, 0.60, 716.791, 1.20861},   // saturated vapour
      {223.15, 0.60, 23240.3, 14.2684},   // saturated liquid
      {273.15, 0.70, 4038.89, 6.13145},   // saturated vapour
      {273.15, 0.70, 14337.0, 11.6669},   // saturated liquid
    };
    for (const auto& pt : pts) {
        CAPTURE(pt.T, pt.z, pt.rho);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&Argon"));
        auto* be = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(be != nullptr);
        be->set_mole_fractions({pt.z, 1 - pt.z});
        be->update_DmolarT_direct(pt.rho, pt.T);
        CHECK(be->p() / 1e6 == Catch::Approx(pt.p_MPa).epsilon(3e-5));
    }
}

TEST_CASE("EOS-CG-2021 SO2-HCl bubble matches Fig. 31 (290 K)", "[mixtures][EOS-CG][CCS]") {
    // Bubble line of Fig. 31 (x-axis is x_HCl); read at mid composition, figure-scatter tolerance.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "SulfurDioxide&HydrogenChloride"));
    struct B
    {
        double xHCl, p_MPa;
    };
    for (const auto& b : std::vector<B>{{0.40, 1.49}, {0.60, 2.20}}) {
        CAPTURE(b.xHCl);
        AS->set_mole_fractions({1 - b.xHCl, b.xHCl});  // component order: SO2, HCl
        AS->update(QT_INPUTS, 0.0, 290.0);
        CHECK(AS->p() / 1e6 == Catch::Approx(b.p_MPa).epsilon(0.08));
    }
}

TEST_CASE("EOS-CG-2021 Cl2-HCl bubble matches Fig. 63 (293 K)", "[mixtures][EOS-CG][CCS]") {
    // Bubble line of Fig. 63 (x-axis is x_HCl).  Pure endpoints: Cl2 ~0.68 MPa, HCl ~4.2 MPa;
    // mid-composition bubble read against the "This work" curve, figure-scatter tolerance.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Chlorine&HydrogenChloride"));
    struct B
    {
        double xHCl, p_MPa;
    };
    for (const auto& b : std::vector<B>{{0.20, 1.60}, {0.40, 2.35}, {0.60, 2.98}, {0.80, 3.57}}) {
        CAPTURE(b.xHCl);
        AS->set_mole_fractions({1 - b.xHCl, b.xHCl});  // component order: Cl2, HCl
        AS->update(QT_INPUTS, 0.0, 293.15);
        CHECK(AS->p() / 1e6 == Catch::Approx(b.p_MPa).epsilon(0.06));
    }
}

TEST_CASE("EOS-CG-2021 SO2-N2 pT-rho reproduces Koester-Vrabec Table S5", "[mixtures][EOS-CG][CCS]") {
    // Single-phase simulated pT-rho points (EOS-CG-2021 Supplementary Material, Table S5, Koester &
    // Vrabec).  State fixed by (T, rho, x); pressure checked with update_DmolarT_direct.  These
    // points were part of the SO2-N2 fit and are reproduced to <1% in the gas/supercritical region
    // sampled here (the SO2-rich dense-liquid points deviate more and are not used).
    struct Pt
    {
        double x, rho_moldm3, T, p_MPa;
    };
    const std::vector<Pt> pts = {
      {0.05, 0.531292, 450.0, 2.0}, {0.05, 0.598917, 400.0, 2.0}, {0.25, 0.619653, 300.0, 1.5},
      {0.25, 0.480516, 500.0, 2.0}, {0.50, 2.173751, 380.0, 6.0}, {0.50, 0.486666, 500.0, 2.0},
      {0.75, 0.496122, 500.0, 2.0}, {0.75, 0.561091, 450.0, 2.0},
    };
    for (const auto& pt : pts) {
        CAPTURE(pt.T, pt.rho_moldm3, pt.x);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "SulfurDioxide&Nitrogen"));
        auto* be = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
        REQUIRE(be != nullptr);
        be->set_mole_fractions({pt.x, 1 - pt.x});
        be->update_DmolarT_direct(pt.rho_moldm3 * 1000.0, pt.T);  // mol/dm^3 -> mol/m^3
        CHECK(be->p() / 1e6 == Catch::Approx(pt.p_MPa).epsilon(0.015));
    }
}

TEST_CASE("EOS-CG-2021 SO2-HCl and Cl2-HCl reproduce Koester-Vrabec VLE (Table S5)", "[mixtures][EOS-CG][CCS]") {
    // Simulated bubble points (EOS-CG-2021 Supplementary Material, Table S5).  At (T, x_liquid) the
    // QT flash gives the bubble pressure and incipient-vapour composition; both are checked against
    // the simulation, which the models represent within molecular-simulation scatter (~5% on p, a
    // few mole-% on y).  Mid-composition points at 273.15 K are used (cleanest simulation region).
    SECTION("SO2 + HCl") {
        struct V
        {
            double x_SO2, p_MPa, y_SO2;
        };
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "SulfurDioxide&HydrogenChloride"));
        for (const auto& v : std::vector<V>{{0.7974, 0.55468, 0.2709}, {0.6691, 0.86487, 0.1749},
                                            {0.4869, 1.27608, 0.0755}, {0.2952, 1.73630, 0.0350}}) {
            CAPTURE(v.x_SO2);
            AS->set_mole_fractions({v.x_SO2, 1 - v.x_SO2});  // component order: SO2, HCl
            AS->update(QT_INPUTS, 0.0, 273.15);
            CHECK(AS->p() / 1e6 == Catch::Approx(v.p_MPa).epsilon(0.05));
            CHECK(AS->mole_fractions_vapor()[0] == Catch::Approx(v.y_SO2).margin(0.04));
        }
    }
    SECTION("Cl2 + HCl") {
        struct V
        {
            double x_Cl2, p_MPa, y_Cl2;
        };
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Chlorine&HydrogenChloride"));
        for (const auto& v : std::vector<V>{{0.8688, 0.80, 0.4454}, {0.8010, 0.99, 0.3464},
                                            {0.7274, 1.20, 0.2816}, {0.5999, 1.44, 0.1802}}) {
            CAPTURE(v.x_Cl2);
            AS->set_mole_fractions({v.x_Cl2, 1 - v.x_Cl2});  // component order: Cl2, HCl
            AS->update(QT_INPUTS, 0.0, 273.15);
            CHECK(AS->p() / 1e6 == Catch::Approx(v.p_MPa).epsilon(0.05));
            CHECK(AS->mole_fractions_vapor()[0] == Catch::Approx(v.y_Cl2).margin(0.04));
        }
    }
}

#endif
