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

TEST_CASE("EOS-CG-2021 CO2-SO2 bubble pressures match Fig. 12 (263 K)", "[mixtures][EOS-CG][CCS]") {
    // Bubble-point pressures read from Fig. 12 at 263 K (digitized; experimental-scatter tolerance).
    struct B
    {
        double xCO2, p_MPa;
    };
    const std::vector<B> bs = {{0.40, 1.28}, {0.60, 1.70}, {0.80, 2.10}};
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&SulfurDioxide"));
    for (const auto& b : bs) {
        CAPTURE(b.xCO2);
        AS->set_mole_fractions({b.xCO2, 1 - b.xCO2});
        AS->update(QT_INPUTS, 0.0, 263.15);  // bubble point
        CHECK(AS->p() / 1e6 == Catch::Approx(b.p_MPa).epsilon(0.04));
    }
}

#endif
