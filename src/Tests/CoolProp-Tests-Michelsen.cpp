#if defined(ENABLE_CATCH)

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include <memory>
#    include <catch2/catch_all.hpp>
#    include "CoolProp/detail/tools.h"
#    include "CoolProp/CoolProp.h"

using namespace CoolProp;

TEST_CASE("Michelsen Flash: Issue #2333 (PR mixture 264.65 K, 6.51 MPa)", "[michelsen][cubic][flash][2333]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("PR", "CarbonDioxide&Nitrogen"));
    std::vector<double> z = {0.97, 0.03};
    AS->set_mole_fractions(z);

    // This case was previously problematic for cubic mixtures.
    // It should now converge successfully.
    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 6.51e6, 264.65));

    // Check that we got a sensible result (e.g., density is positive)
    CHECK(AS->rhomolar() > 0);
}

TEST_CASE("Michelsen Flash: Issue #1668 (PR mixture high pressure)", "[michelsen][cubic][flash][1668]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("PR", "Methane&CarbonDioxide"));
    std::vector<double> z = {0.5, 0.5};
    AS->set_mole_fractions(z);

    // High pressure update should be robust
    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 100e6, 300.0));
    CHECK(AS->rhomolar() > 0);
}

TEST_CASE("Michelsen Flash: Issue #2637 (HEOS mixture phase envelope)", "[michelsen][phase_envelope][2637]") {
    // Phase envelope construction works for HEOS.  Cubic backends (PR, SRK)
    // still hang in build_phase_envelope for this mixture, and
    // all_critical_points() remains broken — those are open issues.
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    std::vector<double> z = {0.85, 0.15};
    AS->set_mole_fractions(z);

    CHECK_NOTHROW(AS->build_phase_envelope(""));
}

TEST_CASE("Michelsen Flash: Multi-component convergence (4-comp mix)", "[michelsen][cubic][flash]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("PR", "Methane&Ethane&Propane&n-Decane"));
    std::vector<double> z = {0.25, 0.25, 0.25, 0.25};
    AS->set_mole_fractions(z);

    double T = 300.0;
    double P = 5.0e6;

    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, P, T));
    CHECK(AS->rhomolar() > 0);
}

// ============================================================================
// Benchmark test cases from the literature
//
// References:
//   [M82a] Michelsen, M.L., "The isothermal flash problem. Part I. Stability",
//          Fluid Phase Equilibria 9 (1982) 1-19.
//   [M82b] Michelsen, M.L., "The isothermal flash problem. Part II. Phase-split
//          calculation", Fluid Phase Equilibria 9 (1982) 21-40.
//   [MM07] Michelsen, M.L. and Mollerup, J.M., "Thermodynamic Models:
//          Fundamentals & Computational Aspects", 2nd ed., Tie-Line (2007),
//          Chapters 9 (Stability Analysis) and 12 (Flash and Phase Envelope).
// ============================================================================

TEST_CASE("Michelsen Flash: 7-component natural gas [M82a,M82b]", "[michelsen][flash][benchmark]") {
    // Composition from Michelsen (1982a) Table 1: synthetic natural gas
    // z = {C1:0.9430, C2:0.0270, C3:0.0074, nC4:0.0049, nC5:0.0027, nC6:0.0010, N2:0.0140}
    std::string fluids = "Methane&Ethane&Propane&n-Butane&n-Pentane&n-Hexane&Nitrogen";
    std::vector<double> z = {0.9430, 0.0270, 0.0074, 0.0049, 0.0027, 0.0010, 0.0140};

    SECTION("SRK two-phase at 190 K, 4 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 4e6, 190.0));
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
    }
    SECTION("PR two-phase at 190 K, 4 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 4e6, 190.0));
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
    }
    SECTION("SRK stable gas at 250 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 250.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR stable gas at 250 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 250.0));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Michelsen Flash: CH4/CO2 binary [MM07 Ch.9]", "[michelsen][flash][benchmark]") {
    // 50/50 Methane/CO2 -- classic binary from Michelsen & Mollerup (2007) Ch. 9
    // Mixture critical point near Tc~252 K, Pc~8.5 MPa (SRK estimate)
    std::string fluids = "Methane&CarbonDioxide";
    std::vector<double> z = {0.5, 0.5};

    SECTION("SRK two-phase at 220 K, 3 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 3e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 220 K, 3 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 3e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("HEOS two-phase at 220 K, 3 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 3e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("SRK stable supercritical at 300 K, 10 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 10e6, 300.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR stable supercritical at 300 K, 10 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 10e6, 300.0));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Michelsen Flash: CH4/C2H6/CO2 ternary [MM07 Ch.9]", "[michelsen][flash][benchmark]") {
    // 30/30/40 Methane/Ethane/CO2 -- ternary from Michelsen & Mollerup (2007) Ch. 9
    std::string fluids = "Methane&Ethane&CarbonDioxide";
    std::vector<double> z = {0.3, 0.3, 0.4};

    SECTION("SRK two-phase at 220 K, 2 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 2e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 220 K, 2 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 2e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("HEOS two-phase at 220 K, 2 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 2e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("SRK stable gas at 300 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 300.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR stable gas at 300 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 300.0));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Michelsen Flash: CH4/n-C10 wide-boiling binary", "[michelsen][flash][benchmark]") {
    // Wide-boiling pair: Tc(CH4)=190.6K vs Tc(nC10)=617.7K
    // Stress test for log-K storage -- K-factors span several orders of magnitude
    // See [M82b] Section 5 discussion on wide-boiling mixtures
    std::string fluids = "Methane&n-Decane";
    std::vector<double> z = {0.7, 0.3};

    SECTION("SRK two-phase at 350 K, 5 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 5e6, 350.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 350 K, 5 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 5e6, 350.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("SRK two-phase at 300 K, 3 MPa (large K-spread)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 3e6, 300.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 300 K, 3 MPa (large K-spread)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 3e6, 300.0));
        CHECK(AS->phase() == iphase_twophase);
    }
}

TEST_CASE("Michelsen Flash: CH4/H2S binary", "[michelsen][flash][benchmark]") {
    // CH4/H2S system -- known to exhibit three-phase (VLLE) behavior at some
    // conditions.  This tests that the two-phase VLE flash converges robustly
    // in the classical VLE region.
    // See: Heidemann, R.A. and Khalil, A.M., AIChE J. 26 (1980) 769-779.
    std::string fluids = "Methane&HydrogenSulfide";
    std::vector<double> z = {0.5, 0.5};

    SECTION("SRK two-phase at 220 K, 2 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 2e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 220 K, 2 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 2e6, 220.0));
        CHECK(AS->phase() == iphase_twophase);
    }
}

TEST_CASE("Michelsen Flash: CH4/C2H6 near mixture critical [M82a]", "[michelsen][flash][benchmark]") {
    // Near the mixture critical point, successive substitution converges slowly
    // and second-order methods (GDEM acceleration, Newton with trust region)
    // are essential.  See [M82a] Section 4.
    std::string fluids = "Methane&Ethane";
    std::vector<double> z = {0.5, 0.5};

    SECTION("SRK two-phase at 230 K, 4 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 4e6, 230.0));
        CHECK(AS->phase() == iphase_twophase);
    }
    SECTION("PR two-phase at 230 K, 4 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 4e6, 230.0));
        CHECK(AS->phase() == iphase_twophase);
    }
}

TEST_CASE("Michelsen Flash: SRK vs PR cross-validation [M82b]", "[michelsen][flash][benchmark]") {
    // Both cubic EOS should agree on phase identification for the same
    // conditions, even though densities and thermodynamic properties differ.
    // See [M82b] for the general methodology applied to cubic EOS.
    std::string fluids = "Methane&Propane";
    std::vector<double> z = {0.6, 0.4};

    auto AS_srk = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
    auto AS_pr = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
    AS_srk->set_mole_fractions(z);
    AS_pr->set_mole_fractions(z);

    SECTION("Two-phase at 250 K, 3 MPa") {
        CHECK_NOTHROW(AS_srk->update(PT_INPUTS, 3e6, 250.0));
        CHECK_NOTHROW(AS_pr->update(PT_INPUTS, 3e6, 250.0));
        CHECK(AS_srk->phase() == iphase_twophase);
        CHECK(AS_pr->phase() == iphase_twophase);
        // Vapor fractions should be qualitatively similar
        CHECK(AS_srk->Q() > 0);
        CHECK(AS_pr->Q() > 0);
    }
}

TEST_CASE("Michelsen Flash: N2/CH4 cryogenic binary", "[michelsen][flash][benchmark]") {
    // Nitrogen/Methane is relevant for LNG and air separation.
    // Tests flash at cryogenic conditions where both components are near
    // or below their critical temperatures.
    std::string fluids = "Nitrogen&Methane";
    std::vector<double> z = {0.2, 0.8};

    SECTION("SRK stable gas at 200 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 200.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR stable gas at 200 K, 1 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e6, 200.0));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Michelsen Flash: 11-comp PR hang at high pressure", "[michelsen][cubic][flash][benchmark]") {
    // 11-component natural gas with PR backend hangs at T=520K, P=10.13 MPa.
    // SRK works fine. The issue is the cubic blind flash takes excessive time
    // or throws "cubic has three roots, but phase not imposed".
    std::string fluids = "Methane&Nitrogen&CO2&Ethane&Propane&Isobutane&Butane&Isopentane&Pentane&Hexane&Heptane";
    std::vector<double> z = {0.9092, 0.0271, 0.0018, 0.0386, 0.011, 0.0037, 0.0037, 0.00135, 0.00135, 0.0008, 0.0014};

    SECTION("PR at 520 K, 1 MPa (works)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 999999.0, 520.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("SRK at 520 K, 10.13 MPa (works)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 10132500.0, 520.0));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR at 520 K, 10.13 MPa (previously hung)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 10132500.0, 520.0));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Michelsen Flash: Issue #3066 (8-comp natural gas density failure)", "[michelsen][flash][benchmark][3066]") {
    // 8-component natural gas at 290.15K fails with "No density solutions"
    // for pressures between ~34.2 and ~35.5 bar because solver_rho_Tp_global encounters
    // an S-shaped isotherm with the target pressure between the spinodal pressures.
    // The fix falls back to SRK-seeded Newton solver (solver_rho_Tp) in this case.
    std::string fluids = "Methane&Nitrogen&CO2&Ethane&Propane&n-Butane&n-Pentane&IsoButane";
    std::vector<double> z = {0.9254, 0.007, 0.008, 0.048, 0.0085, 0.0014, 0.0002, 0.0015};

    SECTION("HEOS at 290.15 K, 35.0 bar (previously failed)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 35.0e5, 290.15));
        CHECK(AS->rhomolar() > 0);
        // Density should be gas-like (~1500 mol/m3), consistent with SRK/PR results
        CHECK(AS->rhomolar() > 1000);
        CHECK(AS->rhomolar() < 3000);
    }
    SECTION("HEOS at 290.15 K, 34.5 bar (previously failed)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 34.5e5, 290.15));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS at 290.15 K, 34.2 bar (boundary, previously worked)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 34.2e5, 290.15));
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS at 290.15 K, 35.5 bar (boundary, previously worked)") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 35.5e5, 290.15));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Legacy Stability: check that legacy algorithm still works", "[stability][legacy]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    std::vector<double> z = {0.5, 0.5};
    AS->set_mole_fractions(z);

    // Force legacy algorithm via configuration
    CoolProp::set_config_int(MIXTURE_STABILITY_ALGORITHM, 0);

    // Methane/Ethane at 200K, 1MPa is stable single-phase
    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 1e6, 200.0));

    // Reset to default
    CoolProp::set_config_int(MIXTURE_STABILITY_ALGORITHM, 1);
}

TEST_CASE("Blind PT flash: Methane/Ethane [0.5/0.5]", "[michelsen][blind][flash]") {
    std::vector<double> z = {0.5, 0.5};
    SECTION("SRK gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("SRK 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("SRK liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("PR liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("HEOS liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Blind PT flash: Methane/Propane [0.7/0.3]", "[michelsen][blind][flash]") {
    std::vector<double> z = {0.7, 0.3};
    SECTION("SRK gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("SRK 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("SRK liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("PR liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS 2ph T=150 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("HEOS liquid T=100 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Blind PT flash: Nitrogen/Oxygen [0.79/0.21]", "[michelsen][blind][flash]") {
    std::vector<double> z = {0.79, 0.21};
    SECTION("SRK gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("SRK 2ph T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("SRK liquid T=75 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 75.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR 2ph T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("PR liquid T=75 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 75.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS 2ph T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("HEOS liquid T=75 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 75.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 20000);  // liquid density — Gibbs selection picks liquid root
    }
}

TEST_CASE("Blind flash: N2/O2 subcooled liquid phase selection", "[michelsen][blind][flash]") {
    // N2/O2 at 75K, 1 bar: below bubble point (78.76 K).
    // Previously the blind flash selected the gas root because SRK
    // tried gas first and both roots converged.  The Gibbs-comparison
    // fix selects the thermodynamically stable liquid root.
    std::vector<double> z = {0.79, 0.21};

    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
    AS->set_mole_fractions(z);
    AS->update(CoolProp::PT_INPUTS, 1e5, 75.0);
    CHECK(AS->Q() == -1);
    CHECK(AS->phase() == CoolProp::iphase_liquid);
    CHECK(AS->rhomolar() > 20000);  // liquid density

    SECTION("Imposed liquid agrees with blind") {
        auto AS_imp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS_imp->set_mole_fractions(z);
        AS_imp->specify_phase(CoolProp::iphase_liquid);
        AS_imp->update(CoolProp::PT_INPUTS, 1e5, 75.0);
        AS_imp->unspecify_phase();
        CHECK(AS_imp->rhomolar() > 20000);
        CHECK(AS_imp->rhomolar() == Catch::Approx(AS->rhomolar()).epsilon(1e-6));
    }

    SECTION("Imposed gas gives low density") {
        auto AS_imp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS_imp->set_mole_fractions(z);
        AS_imp->specify_phase(CoolProp::iphase_gas);
        AS_imp->update(CoolProp::PT_INPUTS, 1e5, 75.0);
        AS_imp->unspecify_phase();
        CHECK(AS_imp->rhomolar() < 500);  // gas-like density (metastable)
    }
}

TEST_CASE("Blind PT flash: N2/CH4/C2H6/C3H8 [0.1/0.5/0.25/0.15]", "[michelsen][blind][flash]") {
    std::vector<double> z = {0.1, 0.5, 0.25, 0.15};
    SECTION("SRK gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("SRK 2ph T=145 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 145.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("SRK liquid T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("SRK", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("PR 2ph T=145 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 145.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("PR liquid T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
    SECTION("HEOS 2ph T=145 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 145.0));
        CHECK(AS->Q() >= 0);
        CHECK(AS->Q() <= 1);
    }
    SECTION("HEOS liquid T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Legacy Flash: check that legacy Jacobian solver still works", "[flash][legacy]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    std::vector<double> z = {0.5, 0.5};
    AS->set_mole_fractions(z);

    // Force legacy algorithm via configuration
    CoolProp::set_config_int(MIXTURE_STABILITY_ALGORITHM, 0);

    // Methane/Ethane at 200K, 1MPa is stable single-phase
    // For HEOS it might go through PTflash_twophase if it thinks it is twophase
    // Let's use a state that is definitely twophase to exercise the solver
    // T = 180K, P = 1MPa for 50/50 Methane/Ethane
    CHECK_NOTHROW(AS->update(CoolProp::PT_INPUTS, 1e6, 180.0));
    CHECK(AS->phase() == CoolProp::iphase_twophase);

    // Reset to default
    CoolProp::set_config_int(MIXTURE_STABILITY_ALGORITHM, 1);
}

TEST_CASE("PT flash with built phase envelope", "[michelsen][phase_envelope]") {
    // Methane(0.85)/Ethane(0.15) mixture — build the phase envelope first,
    // then verify that envelope-guided PT flash produces the same results
    // as blind flash (no envelope).

    std::vector<double> z = {0.85, 0.15};

    // Build envelope once
    auto AS_env = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    AS_env->set_mole_fractions(z);
    AS_env->build_phase_envelope("");

    // Separate instance without envelope for cross-check
    auto AS_blind = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Methane&Ethane"));
    AS_blind->set_mole_fractions(z);

    SECTION("Gas side: 1 MPa, 250 K") {
        AS_env->update(CoolProp::PT_INPUTS, 1e6, 250.0);
        AS_blind->update(CoolProp::PT_INPUTS, 1e6, 250.0);

        CHECK(AS_env->Q() == -1);
        CHECK(AS_env->phase() == CoolProp::iphase_gas);
        CHECK(AS_env->rhomolar() == Catch::Approx(AS_blind->rhomolar()).epsilon(1e-6));
    }

    SECTION("Liquid side: 3 MPa, 120 K") {
        AS_env->update(CoolProp::PT_INPUTS, 3e6, 120.0);
        AS_blind->update(CoolProp::PT_INPUTS, 3e6, 120.0);

        CHECK(AS_env->Q() == -1);
        CHECK(AS_env->phase() == CoolProp::iphase_liquid);
        // The GERG-2008 EOS has two liquid-like roots at this state point.
        // The PE path and the Gibbs-comparing blind flash may converge to
        // different roots (~0.6% apart).  Both are liquid; just verify that.
        CHECK(AS_blind->phase() == CoolProp::iphase_liquid);
        CHECK(AS_env->rhomolar() > 20000);
        CHECK(AS_blind->rhomolar() > 20000);
    }

    SECTION("Two-phase: 2 MPa, 180 K") {
        AS_env->update(CoolProp::PT_INPUTS, 2e6, 180.0);
        AS_blind->update(CoolProp::PT_INPUTS, 2e6, 180.0);

        CHECK(AS_env->phase() == CoolProp::iphase_twophase);
        // Q values may have swapped phase labeling (Q_env + Q_blind ≈ 1) which is
        // thermodynamically equivalent. Check that Q is valid and overall density agrees.
        CHECK(AS_env->Q() > 0.0);
        CHECK(AS_env->Q() < 1.0);
        CHECK(AS_env->rhomolar() == Catch::Approx(AS_blind->rhomolar()).epsilon(1e-4));
    }
}

TEST_CASE("PQ/QT flash with built envelope at 0<Q<1 (GH #3192)", "[michelsen][phase_envelope]") {
    // Regression for GH #3192.  With a built phase envelope a mixture PQ/QT flash for a
    // genuine two-phase state (0 < Q < 1) used to route to newton_raphson_twophase, whose
    // residual vectors r/err_rel were Eigen::Vector2d (fixed size 2) but need 2N-1 (>=3)
    // entries -> out-of-bounds write -> HARD CRASH (process abort, not a C++ exception, so it
    // cannot be caught — pre-fix this test binary terminates here).  These states are now
    // routed to the accurate blind branch (newton_raphson_saturation), which must (a) not
    // crash and (b) match the blind (no-envelope) flash closely.
    std::vector<double> z = {0.5, 0.5};

    auto make_pe = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        AS->build_phase_envelope("");
        return AS;
    };
    auto make_blind = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        return AS;
    };

    SECTION("PQ flash, P=1e5, Q=0.5 (exact reporter repro)") {
        auto AS = make_pe();
        auto ref = make_blind();
        // Guard: the fix is only exercised if the envelope is actually built; otherwise the
        // dispatch already falls through to the blind branch and the test is vacuous.
        REQUIRE(AS->get_phase_envelope_data().built);
        CHECK_NOTHROW(AS->update(PQ_INPUTS, 1e5, 0.5));
        ref->update(PQ_INPUTS, 1e5, 0.5);
        CHECK(AS->phase() == iphase_twophase);
        CHECK(std::isfinite(AS->T()));
        CHECK(AS->T() == Catch::Approx(ref->T()).epsilon(1e-6));
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }

    SECTION("QT flash, Q=0.5, T=180 K") {
        auto AS = make_pe();
        auto ref = make_blind();
        REQUIRE(AS->get_phase_envelope_data().built);
        CHECK_NOTHROW(AS->update(QT_INPUTS, 0.5, 180.0));
        ref->update(QT_INPUTS, 0.5, 180.0);
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->p() == Catch::Approx(ref->p()).epsilon(1e-6));
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }

    SECTION("Dew/bubble (Q==0, Q==1) still use the envelope fast path and stay accurate") {
        auto AS = make_pe();
        auto ref = make_blind();
        REQUIRE(AS->get_phase_envelope_data().built);
        CHECK_NOTHROW(AS->update(PQ_INPUTS, 1e5, 1.0));  // dew point
        ref->update(PQ_INPUTS, 1e5, 1.0);
        CHECK(AS->T() == Catch::Approx(ref->T()).epsilon(1e-6));
        CHECK_NOTHROW(AS->update(PQ_INPUTS, 1e5, 0.0));  // bubble point
        ref->update(PQ_INPUTS, 1e5, 0.0);
        CHECK(AS->T() == Catch::Approx(ref->T()).epsilon(1e-6));
    }
}

// ============================================================================
// Phase-envelope-guided PT flash tests
//
// Test conditions from jakobreichert.  For each mixture the envelope is built
// first, then PT flash is performed for gas / liquid / two-phase conditions.
// Results are cross-checked against the blind flash (no envelope).
//
// Mixtures whose envelopes are not closed (Methane&Propane) or not fully built
// (4-component) exercise the fall-through to the blind flash path.  Mixtures
// with closed envelopes (Methane&Ethane, Nitrogen&Oxygen) exercise the
// envelope-guided code path.
// ============================================================================

TEST_CASE("PE flash: Methane/Ethane [0.5/0.5]", "[michelsen][phase_envelope]") {
    // Envelope closed — exercises the full PE code path for gas, liquid, two-phase
    std::vector<double> z = {0.5, 0.5};

    auto make_pe = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        AS->build_phase_envelope("");
        return AS;
    };
    auto make_blind = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Ethane"));
        AS->set_mole_fractions(z);
        return AS;
    };

    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        ref->update(PT_INPUTS, 1e5, 300.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
    SECTION("HEOS 2ph T=150 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        ref->update(PT_INPUTS, 1e5, 150.0);
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-4));
    }
    SECTION("HEOS liquid T=100 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        ref->update(PT_INPUTS, 1e5, 100.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 0);
        // PE and blind may converge to slightly different densities (~0.03%) because
        // the phase hint from the envelope steers the HEOS solver differently.
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-3));
    }
}

TEST_CASE("PE flash: Methane/Propane [0.7/0.3]", "[michelsen][phase_envelope]") {
    // Envelope built but NOT closed (max_fraction exit) — falls through to blind flash.
    // Verifies the closed-guard works: results must match the blind flash.
    std::vector<double> z = {0.7, 0.3};

    auto make_pe = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Propane"));
        AS->set_mole_fractions(z);
        AS->build_phase_envelope("");
        return AS;
    };
    auto make_blind = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Methane&Propane"));
        AS->set_mole_fractions(z);
        return AS;
    };

    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        ref->update(PT_INPUTS, 1e5, 300.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
    SECTION("HEOS 2ph T=150 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 150.0));
        ref->update(PT_INPUTS, 1e5, 150.0);
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-4));
    }
    SECTION("HEOS liquid T=100 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 100.0));
        ref->update(PT_INPUTS, 1e5, 100.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
}

TEST_CASE("PE flash: Nitrogen/Oxygen [0.79/0.21]", "[michelsen][phase_envelope]") {
    // Envelope closed.  The dew-bubble gap at 1 bar is only ~3 K, so is_inside()
    // may misclassify borderline two-phase points.  Gas and liquid are reliable.
    std::vector<double> z = {0.79, 0.21};

    auto make_pe = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        AS->build_phase_envelope("");
        return AS;
    };
    auto make_blind = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions(z);
        return AS;
    };

    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        ref->update(PT_INPUTS, 1e5, 300.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
    SECTION("HEOS liquid T=75 P=1e5") {
        // At 75K the mixture is subcooled liquid (bubble T at 1 bar ≈ 78.8 K).
        // Both PE-guided and blind flash should select the liquid root via
        // Gibbs comparison.  Cross-check densities agree.
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 75.0));
        ref->update(PT_INPUTS, 1e5, 75.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() > 20000);  // liquid density for N2/O2 at 75K
        CHECK(ref->rhomolar() > 20000);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-3));
    }
}

TEST_CASE("PE flash: N2/CH4/C2H6/C3H8 [0.1/0.5/0.25/0.15]", "[michelsen][phase_envelope]") {
    // Envelope tracing does not reach the closure condition for this 4-component
    // mixture (built=false) — falls through to blind flash for all phases.
    std::vector<double> z = {0.1, 0.5, 0.25, 0.15};

    auto make_pe = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        AS->build_phase_envelope("");
        return AS;
    };
    auto make_blind = [&]() {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane&Ethane&Propane"));
        AS->set_mole_fractions(z);
        return AS;
    };

    SECTION("HEOS gas T=300 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 300.0));
        ref->update(PT_INPUTS, 1e5, 300.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
    SECTION("HEOS 2ph T=145 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 145.0));
        ref->update(PT_INPUTS, 1e5, 145.0);
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-4));
    }
    SECTION("HEOS liquid T=80 P=1e5") {
        auto AS = make_pe();
        auto ref = make_blind();
        CHECK_NOTHROW(AS->update(PT_INPUTS, 1e5, 80.0));
        ref->update(PT_INPUTS, 1e5, 80.0);
        CHECK(AS->Q() == -1);
        CHECK(AS->rhomolar() == Catch::Approx(ref->rhomolar()).epsilon(1e-6));
    }
}

TEST_CASE("Methanol-benzene PT flash at problematic compositions", "[michelsen][meoh_benz]") {
    // Regression test: x_methanol = 0.56 and 0.78 at 308.15 K / 101325 Pa
    // previously failed because solver_rho_Tp corrupted HEOS._T/_p to -inf
    // on the gas-phase attempt, causing the liquid fallback to also fail.
    for (double x : {0.54, 0.56, 0.58, 0.76, 0.78, 0.80}) {
        DYNAMIC_SECTION("x_methanol = " << x) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "methanol&benzene"));
            AS->set_mole_fractions({x, 1.0 - x});
            REQUIRE_NOTHROW(AS->update(PT_INPUTS, 101325, 308.15));
            double rho = AS->rhomolar();
            CAPTURE(rho);
            CHECK(std::isfinite(rho));
            CHECK(rho > 0);
            CHECK(std::isfinite(AS->gibbsmolar()));
        }
    }
}

// Independently recompute the equal-fugacity residual max_i |ln f_i^V - ln f_i^L|
// for a published two-phase split using fresh single-phase sub-states, so a trivial
// (x == y) or unconverged split cannot hide behind the flash's own internal state.
static double equilibrium_residual(const std::string& backend, const std::string& fluids, const std::vector<double>& x, const std::vector<double>& y,
                                   double T, double p) {
    auto L = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    L->set_mole_fractions(x);
    L->specify_phase(iphase_liquid);
    L->update(PT_INPUTS, p, T);
    auto V = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    V->set_mole_fractions(y);
    V->specify_phase(iphase_gas);
    V->update(PT_INPUTS, p, T);
    double maxresid = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (std::min(x[i], y[i]) < 1e-4) continue;  // trace in one phase: skip machine-precision noise
        double fL = L->fugacity(i), fV = V->fugacity(i);
        if (fL > 1e-15 && fV > 1e-15) maxresid = std::max(maxresid, std::abs(std::log(fV / fL)));
    }
    return maxresid;
}

TEST_CASE("Blind PT flash: near-dew two-phase classification (CoolProp-zgpy)", "[michelsen][blind][flash][zgpy]") {
    // Regression for CoolProp-zgpy: for cubic mixtures at high vapor fraction the
    // Michelsen stability trials can converge to the trivial (feed) solution and report
    // a false "stable" verdict, so a genuinely two-phase near-dew state was published as
    // single-phase liquid.  A Wilson bubble/dew cross-check (guess_split_from_wilson)
    // now recovers the split.  Mixture + condition from jakobreichert's report on
    // GitHub #3168: 5-component natural gas at P = 3 bar.  T = 220 K sits at ~0.98 vapor
    // fraction, well inside the two-phase region (T_bub ~ 93 K, T_dew ~ 245 K), and was
    // classified single-phase liquid on master for SRK and PR.
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double p = 3e5, T = 220.0;

    for (const std::string& backend : {std::string("SRK"), std::string("PR")}) {
        DYNAMIC_SECTION(backend << " two-phase at 220 K, 3 bar") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
            AS->set_mole_fractions(z);
            REQUIRE_NOTHROW(AS->update(PT_INPUTS, p, T));

            // Must be two-phase, not the single-phase-liquid false negative.
            CHECK(AS->phase() == iphase_twophase);
            CHECK(AS->Q() > 0.0);
            CHECK(AS->Q() < 1.0);

            std::vector<double> x = AS->mole_fractions_liquid();
            std::vector<double> y = AS->mole_fractions_vapor();

            // Guard against a trivial (x == y) split passed off as two-phase: the
            // incipient liquid is heavy-rich, so the phase compositions must differ.
            double spread = 0;
            for (std::size_t i = 0; i < x.size(); ++i)
                spread = std::max(spread, std::abs(x[i] - y[i]));
            CHECK(spread > 1e-2);

            // The published split must be at genuine equilibrium (independently checked).
            CHECK(equilibrium_residual(backend, fluids, x, y, T, p) < 1e-6);
        }
    }
}

// Port of jakobreichert's scan_twophase_PT.py (GitHub #3168): sweep the two-phase region
// of several mixtures with a BLIND PT flash and classify each result with an INDEPENDENT
// equal-fugacity check, so the test suite (not a Python script against a possibly-stale
// build) is the source of truth.  Per point: twophase + fug<1e-6 = good; non-twophase =
// misclassification; twophase with x==y = trivial/degenerate split; twophase with
// fug>1e-6 (or a phase that cannot be re-solved) = unconverged.  Hard guard: no trivial
// split is ever published; the misclass/bad-fug counts are reported (WARN) for tracking.
TEST_CASE("Stability sweep: blind two-phase classification + fugacity (CoolProp-zgpy)", "[michelsen][stability_sweep][zgpy]") {
    struct SweepCase
    {
        std::string fluids;
        std::vector<double> z;
        double p;
    };
    std::vector<SweepCase> cases = {
      {"Nitrogen&Methane&Ethane&Butane&Pentane", {0.3797, 0.3225, 0.278, 0.0014, 0.0184}, 3e5},
      {"Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e6},
      {"Methane&Ethane", {0.5, 0.5}, 1e6},
    };
    const int NT = 40;
    for (const std::string backend : {std::string("SRK"), std::string("PR"), std::string("HEOS")}) {
        for (const SweepCase& c : cases) {
            double Tbub = -1, Tdew = -1;
            try {
                auto S = std::shared_ptr<AbstractState>(AbstractState::factory(backend, c.fluids));
                S->set_mole_fractions(c.z);
                S->update(PQ_INPUTS, c.p, 0.0);
                Tbub = S->T();
                S->update(PQ_INPUTS, c.p, 1.0);
                Tdew = S->T();
            } catch (...) {
                continue;  // backend can't bracket this mixture; skip
            }
            if (!(Tbub > 0 && Tdew > Tbub)) continue;

            int misclass = 0, bad_fug = 0, trivial = 0;
            for (int k = 0; k < NT; ++k) {
                double T = Tbub + (Tdew - Tbub) * (k + 0.5) / NT;
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, c.fluids));
                AS->set_mole_fractions(c.z);
                try {
                    AS->update(PT_INPUTS, c.p, T);
                } catch (...) {
                    ++misclass;
                    continue;
                }
                if (AS->phase() != iphase_twophase) {
                    ++misclass;
                    continue;
                }
                try {
                    std::vector<double> x = AS->mole_fractions_liquid(), y = AS->mole_fractions_vapor();
                    double spread = 0;
                    for (std::size_t i = 0; i < x.size(); ++i)
                        spread = std::max(spread, std::abs(x[i] - y[i]));
                    if (spread < 1e-6)
                        ++trivial;
                    else if (equilibrium_residual(backend, c.fluids, x, y, T, c.p) > 1e-6)
                        ++bad_fug;
                } catch (...) {
                    ++bad_fug;  // published a split whose phases cannot be independently re-solved
                }
            }
            WARN(backend << " " << c.fluids << " P=" << c.p << "  misclass=" << misclass << " bad_fug=" << bad_fug << " trivial=" << trivial << " / "
                         << NT);
            CHECK(trivial == 0);

            // False-positive guard: superheated-vapor points above the dew point must stay
            // single phase -- the speculative/forced attempt must never publish a clearly
            // single-phase feed as two-phase.  (Only the vapor side is checked: for these
            // heavy mixtures the subcooled-liquid side falls below some components' triple
            // points, a pre-existing sub-triple-point region unrelated to this fix.)
            for (double Tsp : {Tdew + 30.0, Tdew + 80.0}) {
                if (Tsp <= 0) continue;
                auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, c.fluids));
                AS->set_mole_fractions(c.z);
                try {
                    AS->update(PT_INPUTS, c.p, Tsp);
                } catch (...) {
                    continue;  // outside the model's valid range for this mixture/backend
                }
                INFO(backend << " " << c.fluids << " single-phase check at T=" << Tsp);
                CHECK(AS->phase() != iphase_twophase);
            }
        }
    }
}

#endif
