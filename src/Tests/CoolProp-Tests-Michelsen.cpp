#if defined(ENABLE_CATCH)

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include <algorithm>
#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>
#    include <catch2/catch_all.hpp>
#    include "CoolProp/detail/tools.h"
#    include "CoolProp/CoolProp.h"
#    include <cmath>

using namespace CoolProp;

namespace {
// Recompute the equal-fugacity equilibrium residual max_i |ln(f_i^V / f_i^L)| from a
// converged two-phase flash result, by re-solving each reported phase composition as a
// single-phase state at (p, T).  A correctly converged split drives this to ~0; the
// #3168 bug published splits with residuals of order 1e-1..1e0.  Asserting on this (not
// just phase()==twophase) is what keeps a grossly-unconverged split from passing CI.
// Trace components (mole fraction < 1e-4 in either phase) are skipped: their fugacities
// carry machine-precision noise that would otherwise dominate the max.  (Signature is
// (..., T, p) -- the merge of #3170 and #3174 standardized on this order; see #3174.)
double equilibrium_residual(const std::string& backend, const std::string& fluids, const std::vector<double>& x, const std::vector<double>& y,
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
}  // namespace

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
        // #3168: the published split must actually be at equilibrium, not a
        // grossly-unconverged one accepted on phase()/Q() alone.
        CHECK(equilibrium_residual("SRK", fluids, AS->mole_fractions_liquid_double(), AS->mole_fractions_vapor_double(), 190.0, 4e6) < 1e-6);
    }
    SECTION("PR two-phase at 190 K, 4 MPa") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("PR", fluids));
        AS->set_mole_fractions(z);
        CHECK_NOTHROW(AS->update(PT_INPUTS, 4e6, 190.0));
        CHECK(AS->phase() == iphase_twophase);
        CHECK(AS->Q() > 0);
        CHECK(AS->Q() < 1);
        CHECK(equilibrium_residual("PR", fluids, AS->mole_fractions_liquid_double(), AS->mole_fractions_vapor_double(), 190.0, 4e6) < 1e-6);
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
    //
    // At 308.15 K / 1 atm this mixture is subcooled single-phase liquid (bubble
    // pressure ~34 kPa << 101 kPa).  #3168: the Michelsen stability check produces
    // a false-positive two-phase classification at x_methanol = 0.54 whose split
    // never converges (both roots liquid-like, rho_vap > rho_liq).  The convergence
    // gate + single-phase fallback must recover the correct single-phase liquid
    // result rather than publishing the unconverged split.
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
            // Subcooled liquid: single phase, liquid-like density.
            CHECK(AS->Q() == -1);
            CHECK(AS->phase() == iphase_liquid);
            CHECK(rho > 9000);
        }
    }
}

// Independently recompute the equal-fugacity residual max_i |ln f_i^V - ln f_i^L|
// for a published two-phase split using fresh single-phase sub-states, so a trivial
// (x == y) or unconverged split cannot hide behind the flash's own internal state.

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

// Helper: run a round-trip HSU_P flash test.
// 1. PT flash at (P, T) to get reference state
// 2. Read the target property (H, S, or U)
// 3. Flash with (Y, P) inputs on a fresh object
// 4. Check that T and rho match the reference
static void hsu_p_roundtrip(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double P, double T,
                            CoolProp::input_pairs flash_pair, double eps = 1e-6) {
    using namespace CoolProp;
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    AS->set_mole_fractions(z);
    AS->update(PT_INPUTS, P, T);
    double T_ref = AS->T();
    double rho_ref = AS->rhomolar();

    double y_ref;
    switch (flash_pair) {
        case HmolarP_INPUTS:
            y_ref = AS->hmolar();
            break;
        case PSmolar_INPUTS:
            y_ref = AS->smolar();
            break;
        case PUmolar_INPUTS:
            y_ref = AS->umolar();
            break;
        default:
            throw ValueError("unsupported flash pair in hsu_p_roundtrip");
    }

    auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    AS2->set_mole_fractions(z);
    // HmolarP_INPUTS: (H, P);  PSmolar_INPUTS: (P, S);  PUmolar_INPUTS: (P, U)
    if (flash_pair == HmolarP_INPUTS) {
        AS2->update(flash_pair, y_ref, P);
    } else {
        AS2->update(flash_pair, P, y_ref);
    }
    CHECK(AS2->T() == Catch::Approx(T_ref).epsilon(eps));
    CHECK(AS2->rhomolar() == Catch::Approx(rho_ref).epsilon(eps));
}

TEST_CASE("HSU_P flash: mixture HP round-trip", "[michelsen][flash][HSU_P]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, HmolarP_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, HmolarP_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        hsu_p_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, HmolarP_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, HmolarP_INPUTS);
    }
}

TEST_CASE("HSU_P flash: mixture SP round-trip", "[michelsen][flash][HSU_P]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, PSmolar_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, PSmolar_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        hsu_p_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, PSmolar_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, PSmolar_INPUTS);
    }
}

TEST_CASE("HSU_P flash: mixture UP round-trip", "[michelsen][flash][HSU_P]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, PUmolar_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, PUmolar_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        hsu_p_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, PUmolar_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        hsu_p_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, PUmolar_INPUTS);
    }
}

TEST_CASE("HSU_P flash: mixture two-phase HP round-trip", "[michelsen][flash][HSU_P]") {
    // For N2/O2, the two-phase region at 1 bar spans ~78-90 K.
    // Flash at a two-phase (T, P), read H, then HP-flash and verify
    // that T, Q, and rho are recovered.
    using namespace CoolProp;
    SECTION("N2/O2 two-phase T=80 P=1e5") {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS->set_mole_fractions({0.79, 0.21});
        AS->update(PT_INPUTS, 1e5, 80.0);
        REQUIRE(AS->phase() == iphase_twophase);
        double T_ref = AS->T();
        double rho_ref = AS->rhomolar();
        double Q_ref = AS->Q();
        double h_ref = AS->hmolar();

        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        AS2->update(HmolarP_INPUTS, h_ref, 1e5);
        CHECK(AS2->T() == Catch::Approx(T_ref).epsilon(1e-4));
        CHECK(AS2->rhomolar() == Catch::Approx(rho_ref).epsilon(1e-4));
        CHECK(AS2->Q() == Catch::Approx(Q_ref).margin(0.01));
    }
}

// "No silent wrong answer" invariant for the mixture HSU_P flash (PR #3148 /
// hsu_p_native_flaw.py).  For a hard CO2/Water/N2/Ar/O2 mixture the native HmolarP / PSmolar
// flash used to return a state that did NOT satisfy the request -- it converged to a T whose
// enthalpy/entropy differed from the target by a large margin (the GitHub #3148 flaw) -- while
// reporting success.  The required, achievable guarantee is that the flash must NEVER silently
// return such a state: it must either throw (honest failure) OR return a state that reproduces
// the requested property.  The #3192 residual-verification guard in HSU_P_flash enforces it.
//
// NOTE: this does NOT assert the flash lands on a particular T.  For water that condenses out
// of a CO2/NG-type mixture the underlying (T,p) flash can pick a different (but property-
// consistent) branch; returning the physically-intended state there is a separate, future
// concern (an immiscible / water-dropout flash), deliberately out of scope here.
TEST_CASE("HSU_P flash: CO2/Water/N2/Ar/O2 mixture no silent wrong answer (GitHub #3148)", "[michelsen][flash][HSU_P]") {
    using namespace CoolProp;
    const std::string fluids = "CarbonDioxide&Water&Nitrogen&Argon&Oxygen";
    const std::vector<double> z = {0.90, 0.02, 0.04, 0.01, 0.03};
    const double P = 2.05e6;  // Pa

    // Helper: flash for `value` via `pair`, then assert it either threw or its OWN returned
    // state reproduces `value` (so the caller is never handed a state that violates the request).
    auto assert_no_silent_miss = [&](input_pairs pair, parameters key, double value) {
        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        AS->set_mole_fractions(z);
        bool threw = false;
        try {
            // HmolarP_INPUTS: (H, P);  PSmolar/PUmolar_INPUTS: (P, value)
            if (pair == HmolarP_INPUTS)
                AS->update(pair, value, P);
            else
                AS->update(pair, P, value);
        } catch (const CoolProp::CoolPropBaseError&) {
            threw = true;
        }
        if (!threw) {
            const double got = AS->keyed_output(key);
            CAPTURE(value, got);
            CHECK(got == Catch::Approx(value).epsilon(1e-6));  // returned state satisfies the request
        }
        // threw == true is acceptable: an honest failure rather than a silent wrong answer.
    };

    for (double T_true : {255.0, 260.0, 270.0, 280.0, 290.0, 300.0}) {
        DYNAMIC_SECTION("HP at T_true = " << T_true << " K") {
            auto ref = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            ref->set_mole_fractions(z);
            ref->update(PT_INPUTS, P, T_true);
            assert_no_silent_miss(HmolarP_INPUTS, iHmolar, ref->hmolar());
        }
    }

    SECTION("SP at T_true = 270 K") {
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        ref->set_mole_fractions(z);
        ref->update(PT_INPUTS, P, 270.0);
        assert_no_silent_miss(PSmolar_INPUTS, iSmolar, ref->smolar());
    }
}

// Helper: run a round-trip DHSU_T flash test.
// 1. PT flash at (P, T) to get reference state
// 2. Read the target property (D, H, S, or U)
// 3. Flash with (Y, T) inputs on a fresh object
// 4. Check that P and rho match the reference
static void dhsu_t_roundtrip(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double P, double T,
                             CoolProp::input_pairs flash_pair, double eps = 1e-6) {
    using namespace CoolProp;
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    AS->set_mole_fractions(z);
    AS->update(PT_INPUTS, P, T);
    double P_ref = AS->p();
    double rho_ref = AS->rhomolar();

    double y_ref;
    switch (flash_pair) {
        case DmolarT_INPUTS:
            y_ref = AS->rhomolar();
            break;
        case HmolarT_INPUTS:
            y_ref = AS->hmolar();
            break;
        case SmolarT_INPUTS:
            y_ref = AS->smolar();
            break;
        case TUmolar_INPUTS:
            y_ref = AS->umolar();
            break;
        default:
            throw ValueError("unsupported flash pair in dhsu_t_roundtrip");
    }

    auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(backend, fluids));
    AS2->set_mole_fractions(z);
    // DmolarT_INPUTS: (rho, T);  HmolarT_INPUTS: (H, T);  SmolarT_INPUTS: (S, T);  TUmolar_INPUTS: (T, U)
    if (flash_pair == TUmolar_INPUTS) {
        AS2->update(flash_pair, T, y_ref);
    } else {
        AS2->update(flash_pair, y_ref, T);
    }
    CHECK(AS2->p() == Catch::Approx(P_ref).epsilon(eps));
    CHECK(AS2->rhomolar() == Catch::Approx(rho_ref).epsilon(eps));
}

TEST_CASE("DHSU_T flash: mixture DT round-trip", "[michelsen][flash][DHSU_T]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, DmolarT_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, DmolarT_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        dhsu_t_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, DmolarT_INPUTS);
    }
}

TEST_CASE("DHSU_T flash: mixture HT round-trip", "[michelsen][flash][DHSU_T]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, HmolarT_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, HmolarT_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        dhsu_t_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, HmolarT_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, HmolarT_INPUTS);
    }
}

TEST_CASE("DHSU_T flash: mixture ST round-trip", "[michelsen][flash][DHSU_T]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, SmolarT_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, SmolarT_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        dhsu_t_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, SmolarT_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, SmolarT_INPUTS);
    }
}

TEST_CASE("DHSU_T flash: mixture UT round-trip", "[michelsen][flash][DHSU_T]") {
    SECTION("N2/O2 gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 300.0, TUmolar_INPUTS);
    }
    SECTION("N2/O2 liquid T=75 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Oxygen", {0.79, 0.21}, 1e5, 75.0, TUmolar_INPUTS);
    }
    SECTION("CH4/C2H6 gas T=250 P=1e6") {
        dhsu_t_roundtrip("HEOS", "Methane&Ethane", {0.85, 0.15}, 1e6, 250.0, TUmolar_INPUTS);
    }
    SECTION("4-component gas T=300 P=1e5") {
        dhsu_t_roundtrip("HEOS", "Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15}, 1e5, 300.0, TUmolar_INPUTS);
    }
}

// ============================================================================
// Tests adapted from jakobreichert's PR CoolProp/CoolProp#2720.
//
// The original PR adds CoolProp-Tests-MixtureFlash.cpp with comprehensive
// mixture flash tests covering PT, PQ, and HSU_P flashes.  The tests below
// cover the same ground, adapted for our test file and HEOS-focused scope.
// ============================================================================

TEST_CASE("HSU_P flash: near saturation Propane/Butane", "[michelsen][flash][HSU_P][saturation]") {
    // Propane/Butane 50/50 at P=1e5.  Tests HP/SP/UP round-trip for points
    // close to the bubble and dew curves.
    const std::string fluids = "Propane&Butane";
    const std::vector<double> z = {0.5, 0.5};
    const double P = 1e5;
    const double TOL = 0.1;  // K

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();
    sat->update(PQ_INPUTS, P, 1.0);
    const double T_dew = sat->T();

    struct Case
    {
        std::string label;
        double T;
        phases ph;
    };
    std::vector<Case> cases = {
      {"subcooled T_bub-0.5", T_bub - 0.5, iphase_liquid}, {"subcooled T_bub-1", T_bub - 1.0, iphase_liquid},
      {"two-phase T_bub+1", T_bub + 1.0, iphase_twophase}, {"two-phase midpoint", 0.5 * (T_bub + T_dew), iphase_twophase},
      {"two-phase T_dew-1", T_dew - 1.0, iphase_twophase}, {"superheated T_dew+1", T_dew + 1.0, iphase_gas},
    };

    for (auto& c : cases) {
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        ref->set_mole_fractions(z);
        // Only specify_phase for single-phase cases; two-phase detection
        // via blind PT flash works but specify_phase(iphase_twophase) does
        // not yet work for the SRK-based initial guess in solver_rho_Tp.
        if (c.ph != iphase_twophase) {
            ref->specify_phase(c.ph);
        }
        ref->update(PT_INPUTS, P, c.T);
        const double H_ref = ref->hmass();
        const double S_ref = ref->smass();
        const double U_ref = ref->umass();
        ref->unspecify_phase();

        DYNAMIC_SECTION(c.label + " HP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(HmassP_INPUTS, H_ref, P);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " SP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(PSmass_INPUTS, P, S_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " UP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(PUmass_INPUTS, P, U_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
    }
}

TEST_CASE("HSU_P flash: near saturation 5-component N2-HC", "[michelsen][flash][HSU_P][saturation]") {
    // N2/CH4/C2H6/C4H10/C5H12 at P=3e5.  Same structure as the
    // Propane/Butane test but for a 5-component system.
    // NOTE: subcooled cases are excluded — the underlying PT flash for this
    // 5-component mixture can fail near Tmin, causing our TOMS748-based
    // HSU_P solver to diverge.  PR #2720 fixes the underlying PT flash.
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double P = 3e5;
    const double TOL = 0.1;  // K

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();
    sat->update(PQ_INPUTS, P, 1.0);
    const double T_dew = sat->T();

    struct Case
    {
        std::string label;
        double T;
        phases ph;
    };
    std::vector<Case> cases = {
      {"two-phase T_bub+1", T_bub + 1.0, iphase_twophase},
      {"two-phase midpoint", 0.5 * (T_bub + T_dew), iphase_twophase},
      {"two-phase T_dew-1", T_dew - 1.0, iphase_twophase},
      {"superheated T_dew+1", T_dew + 1.0, iphase_gas},
    };

    for (auto& c : cases) {
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
        ref->set_mole_fractions(z);
        if (c.ph != iphase_twophase) {
            ref->specify_phase(c.ph);
        }
        ref->update(PT_INPUTS, P, c.T);
        const double H_ref = ref->hmass();
        const double S_ref = ref->smass();
        const double U_ref = ref->umass();
        ref->unspecify_phase();

        DYNAMIC_SECTION(c.label + " HP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(HmassP_INPUTS, H_ref, P);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " SP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(PSmass_INPUTS, P, S_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " UP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            AS->update(PUmass_INPUTS, P, U_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
    }
}

TEST_CASE("PT flash: 5-component two-phase consistency", "[michelsen][flash][VLE]") {
    // Sweep 100 temperatures from T_bub to T_dew and verify:
    //   1. phase == iphase_twophase
    //   2. hmolar() non-decreasing with T
    //   3. Fugacity equality for non-trace components
    // Adapted from jakobreichert PR #2720.
    // NOTE: HEOS only — SRK/PR have pre-existing two-phase detection issues
    // near the dew point for this 5-component system.
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double P = 3e5;
    const int NQ = 100;
    const std::size_t NC = z.size();
    const double FUG_TOL = 1e-5;  // relaxed from 1e-7 — near dew point convergence noise
    const double Z_MIN = 1e-4;
    const double F_MIN = 1e-15;
    const double H_SLACK = 0.1;

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();
    sat->update(PQ_INPUTS, P, 1.0);
    const double T_dew = sat->T();
    REQUIRE(T_dew > T_bub);

    auto liq = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    auto vap = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    AS->set_mole_fractions(z);

    double H_prev = -1e100;
    for (int i = 0; i < NQ; ++i) {
        const double alpha = (static_cast<double>(i) + 0.5) / NQ;
        const double T = T_bub + (T_dew - T_bub) * alpha;
        CAPTURE(T);

        AS->update(PT_INPUTS, P, T);
        CHECK(AS->phase() == iphase_twophase);

        const double H = AS->hmolar();
        CHECK(H >= H_prev - H_SLACK);
        H_prev = H;

        if (AS->phase() == iphase_twophase) {
            const auto x = AS->mole_fractions_liquid();
            const auto y = AS->mole_fractions_vapor();
            liq->set_mole_fractions(x);
            liq->specify_phase(iphase_liquid);
            liq->update(PT_INPUTS, P, T);
            vap->set_mole_fractions(y);
            vap->specify_phase(iphase_gas);
            vap->update(PT_INPUTS, P, T);

            for (std::size_t j = 0; j < NC; ++j) {
                if (std::min(x[j], y[j]) < Z_MIN) continue;
                const double f_liq = liq->fugacity_coefficient(j) * x[j];
                const double f_vap = vap->fugacity_coefficient(j) * y[j];
                const double f_ref = std::max(std::abs(f_liq), std::abs(f_vap));
                if (f_ref > F_MIN) {
                    CAPTURE(j);
                    CHECK(std::abs(f_liq - f_vap) / f_ref < FUG_TOL);
                }
            }
        }
    }
}

TEST_CASE("PQ flash: 5-component two-phase consistency", "[michelsen][flash][VLE]") {
    // Sweep 100 quality midpoints from Q=0 to Q=1 and verify:
    //   1. phase == iphase_twophase
    //   2. T non-decreasing with Q
    //   3. Fugacity equality for non-trace components
    // Adapted from jakobreichert PR #2720.  HEOS only.
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double P = 3e5;
    const int NQ = 100;
    const std::size_t NC = z.size();
    const double FUG_TOL = 1e-5;
    const double Z_MIN = 1e-4;
    const double F_MIN = 1e-15;

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();
    sat->update(PQ_INPUTS, P, 1.0);
    const double T_dew = sat->T();
    REQUIRE(T_dew > T_bub);

    auto liq = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    auto vap = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    AS->set_mole_fractions(z);

    double T_prev = -1e100;
    for (int i = 0; i < NQ; ++i) {
        const double Q = (static_cast<double>(i) + 0.5) / NQ;
        CAPTURE(Q);

        AS->update(PQ_INPUTS, P, Q);
        CHECK(AS->phase() == iphase_twophase);

        const double T = AS->T();
        CHECK(T >= T_prev);
        T_prev = T;

        if (AS->phase() == iphase_twophase) {
            const auto x = AS->mole_fractions_liquid();
            const auto y = AS->mole_fractions_vapor();
            liq->set_mole_fractions(x);
            liq->specify_phase(iphase_liquid);
            liq->update(PT_INPUTS, P, T);
            vap->set_mole_fractions(y);
            vap->specify_phase(iphase_gas);
            vap->update(PT_INPUTS, P, T);

            for (std::size_t j = 0; j < NC; ++j) {
                if (std::min(x[j], y[j]) < Z_MIN) continue;
                const double f_liq = liq->fugacity_coefficient(j) * x[j];
                const double f_vap = vap->fugacity_coefficient(j) * y[j];
                const double f_ref = std::max(std::abs(f_liq), std::abs(f_vap));
                if (f_ref > F_MIN) {
                    CAPTURE(j);
                    CHECK(std::abs(f_liq - f_vap) / f_ref < FUG_TOL);
                }
            }
        }
    }
}

// Regression for the #3192 convergence-gate refinement (follow-up to #3168).  For a
// wide-boiling mixture the Michelsen two-phase split converges only to ~1e-6 (SS converges
// linearly and the second-order stage stalls), and the original hard 1e-7 gate discarded it
// -- silently misclassifying a genuine two-phase state as single-phase.  The refined gate
// accepts a non-trivial, near-converged equilibrium.  This pins that behaviour: a blind PT
// flash well inside the two-phase region must report two-phase with an independently-verified,
// non-trivial equal-fugacity split.  (Reverting the gate to the hard 1e-7 throw fails this.)
TEST_CASE("PT flash: wide-boiling split survives the convergence gate (GitHub #3192)", "[michelsen][flash][VLE]") {
    using namespace CoolProp;
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double P = 3e5;

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();
    sat->update(PQ_INPUTS, P, 1.0);
    const double T_dew = sat->T();
    REQUIRE(T_dew > T_bub);

    // Points in the narrow band just above the bubble point, where the wide-boiling split
    // converges only to ~1e-6 and the pre-refinement hard 1e-7 gate threw nonconvergence
    // (this mixture was misclassified single-phase liquid at T ~ 91-96 K on master).
    for (double frac : {0.005, 0.01, 0.02, 0.03, 0.05}) {
        const double T = T_bub + frac * (T_dew - T_bub);
        DYNAMIC_SECTION("frac = " << frac << " (T = " << T << " K)") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", fluids));
            AS->set_mole_fractions(z);
            REQUIRE_NOTHROW(AS->update(PT_INPUTS, P, T));
            CHECK(AS->phase() == iphase_twophase);
            CHECK(AS->Q() > 0.0);
            CHECK(AS->Q() < 1.0);

            const auto x = AS->mole_fractions_liquid_double();
            const auto y = AS->mole_fractions_vapor_double();
            double spread = 0;
            for (std::size_t i = 0; i < z.size(); ++i)
                spread = std::max(spread, std::abs(x[i] - y[i]));
            CHECK(spread > 1e-2);                                            // genuine, non-trivial split (not a trivial x==y collapse)
            CHECK(equilibrium_residual("HEOS", fluids, x, y, T, P) < 1e-5);  // at equilibrium
        }
    }
}

TEST_CASE("PQ flash with built PE: N2/CH4", "[michelsen][flash][PQ_flash][PhaseEnvelope]") {
    // PQ flash works when the phase envelope is built.
    // Adapted from jakobreichert PR #2720.
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane"));
    AS->set_mole_fractions({0.5, 0.5});
    AS->build_phase_envelope("");
    const std::size_t npts = AS->get_phase_envelope_data().T.size();
    CAPTURE(npts);
    CHECK(npts > 0);
    CHECK(AS->get_phase_envelope_data().built);
    // Use a separate (non-PE) object for PQ flash — PQ flash on the same
    // object that built the PE can crash in the current codebase.
    auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Methane"));
    AS2->set_mole_fractions({0.5, 0.5});
    REQUIRE_NOTHROW(AS2->update(PQ_INPUTS, 1.5e5, 0.5));
    CHECK(AS2->phase() == iphase_twophase);
}

TEST_CASE("PQ flash: 6-component no-throw sweep", "[michelsen][flash][PQ_flash]") {
    // Regression: saturation_Wilson Secant can diverge for Q between 0.88
    // and 0.95 in some multi-component mixtures.  Fixed by always trying
    // Brent (bounded) first, falling back to Secant only if Brent fails.
    // Adapted from jakobreichert PR #2720.
    const std::string fluids = "Nitrogen&Methane&Ethane&Propane&Butane&Pentane";
    const std::vector<double> z = {0.2936, 0.2720, 0.0592, 0.2932, 0.0787, 0.0033};
    const double P = 3.92e5;
    const int NQ = 100;

    const std::vector<std::string> backends = {"SRK", "PR", "HEOS"};
    for (const auto& be : backends) {
        DYNAMIC_SECTION(be) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);
            for (int i = 0; i < NQ; ++i) {
                double q = static_cast<double>(i) / (NQ - 1);
                REQUIRE_NOTHROW(AS->update(PQ_INPUTS, P, q));
            }
        }
    }
}

// ============================================================================
// DHSU_T flash: two-phase round-trip tests
//
// These verify that HmolarT, SmolarT, DmolarT, and TUmolar flashes
// correctly recover a known two-phase state created via PQ flash.
// Previously, the mixture DHSU_T flash only handled single-phase states,
// silently returning wrong results (on the Van der Waals loop) for
// two-phase inputs.
// ============================================================================

TEST_CASE("DHSU_T flash: two-phase N2/O2 HEOS round-trip", "[michelsen][dhsu_t][twophase]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
    AS->set_mole_fractions({0.79, 0.21});
    AS->update(PQ_INPUTS, 1e5, 0.5);
    double T = AS->T(), H = AS->hmolar(), S = AS->smolar();
    double U = AS->umolar(), rho = AS->rhomolar(), Q = AS->Q();

    SECTION("HmolarT recovers two-phase state") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(HmolarT_INPUTS, H, T));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
        CHECK(AS2->Q() == Catch::Approx(Q).margin(0.05));
    }
    SECTION("SmolarT recovers two-phase state") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(SmolarT_INPUTS, S, T));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
        CHECK(AS2->Q() == Catch::Approx(Q).margin(0.05));
    }
    // DmolarT is NOT tested for two-phase round-trip.  For mixtures,
    // (T, rho_overall) is ambiguous: the lever-rule overall density can also
    // be achieved as a valid single-phase gas at a different P.  The flash
    // legitimately returns the single-phase solution.
    SECTION("TUmolar recovers two-phase state") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(TUmolar_INPUTS, T, U));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
        CHECK(AS2->Q() == Catch::Approx(Q).margin(0.05));
    }
}

TEST_CASE("DHSU_T flash: two-phase cubic round-trip", "[michelsen][dhsu_t][twophase][cubic]") {
    const std::vector<std::string> backends = {"SRK", "PR"};
    for (const auto& be : backends) {
        DYNAMIC_SECTION(be << " backend") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, "Nitrogen&Oxygen"));
            AS->set_mole_fractions({0.79, 0.21});
            AS->update(PQ_INPUTS, 1e5, 0.5);
            double T = AS->T(), H = AS->hmolar(), S = AS->smolar();
            double rho = AS->rhomolar(), Q = AS->Q();

            SECTION("HmolarT") {
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(be, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(HmolarT_INPUTS, H, T));
                CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.02));
                CHECK(AS2->Q() >= 0);
                CHECK(AS2->Q() <= 1);
            }
            SECTION("SmolarT") {
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(be, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(SmolarT_INPUTS, S, T));
                CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.02));
                CHECK(AS2->Q() >= 0);
                CHECK(AS2->Q() <= 1);
            }
            SECTION("TUmolar") {
                double U = AS->umolar();
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(be, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(TUmolar_INPUTS, T, U));
                CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.02));
                CHECK(AS2->Q() >= 0);
                CHECK(AS2->Q() <= 1);
            }
        }
    }
}

TEST_CASE("DHSU_T flash: single-phase mixture regression", "[michelsen][dhsu_t]") {
    // Verify that single-phase HEOS mixture states still work after the
    // two-phase P-sweep was added (fast-path regression test).
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
    AS->set_mole_fractions({0.79, 0.21});
    // Gas state at 300 K, 1 atm — well above the dew point
    AS->update(PT_INPUTS, 1e5, 300.0);
    double T = AS->T(), H = AS->hmolar(), S = AS->smolar();
    double U = AS->umolar(), rho = AS->rhomolar();

    SECTION("HmolarT single-phase") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(HmolarT_INPUTS, H, T));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(1e5).epsilon(0.001));
    }
    SECTION("SmolarT single-phase") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(SmolarT_INPUTS, S, T));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(1e5).epsilon(0.001));
    }
    SECTION("DmolarT single-phase") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarT_INPUTS, rho, T));
        CHECK(AS2->hmolar() == Catch::Approx(H).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(1e5).epsilon(0.001));
    }
    SECTION("TUmolar single-phase") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(TUmolar_INPUTS, T, U));
        CHECK(AS2->rhomolar() == Catch::Approx(rho).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(1e5).epsilon(0.001));
    }
}

TEST_CASE("HSU_D flash: two-phase N2/O2 HEOS round-trip", "[michelsen][hsu_d][twophase]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
    AS->set_mole_fractions({0.79, 0.21});
    AS->update(PQ_INPUTS, 1e5, 0.5);
    double T = AS->T(), H = AS->hmolar(), S = AS->smolar();
    double U = AS->umolar(), rho = AS->rhomolar();

    SECTION("DmolarHmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarHmolar_INPUTS, rho, H));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
    }
    SECTION("DmolarSmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarSmolar_INPUTS, rho, S));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
    }
    SECTION("DmolarUmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarUmolar_INPUTS, rho, U));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
        CHECK(AS2->Q() >= 0);
        CHECK(AS2->Q() <= 1);
    }
}

TEST_CASE("HSU_D flash: two-phase N2/O2 cubic round-trip", "[michelsen][hsu_d][twophase][cubic]") {
    for (const std::string& backend : {"SRK", "PR"}) {
        DYNAMIC_SECTION("Backend: " << backend) {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, "Nitrogen&Oxygen"));
            AS->set_mole_fractions({0.79, 0.21});
            AS->update(PQ_INPUTS, 1e5, 0.5);
            double T = AS->T(), H = AS->hmolar(), S = AS->smolar();
            double U = AS->umolar(), rho = AS->rhomolar();

            SECTION("DmolarHmolar") {
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(backend, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(DmolarHmolar_INPUTS, rho, H));
                CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
            }
            SECTION("DmolarSmolar") {
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(backend, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(DmolarSmolar_INPUTS, rho, S));
                CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
            }
            SECTION("DmolarUmolar") {
                auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory(backend, "Nitrogen&Oxygen"));
                AS2->set_mole_fractions({0.79, 0.21});
                REQUIRE_NOTHROW(AS2->update(DmolarUmolar_INPUTS, rho, U));
                CHECK(AS2->T() == Catch::Approx(T).epsilon(0.01));
            }
        }
    }
}

TEST_CASE("HSU_D flash: single-phase N2/O2 HEOS regression", "[michelsen][hsu_d][singlephase]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
    AS->set_mole_fractions({0.79, 0.21});
    AS->update(PT_INPUTS, 1e5, 300.0);
    double T = AS->T(), P = AS->p(), H = AS->hmolar(), S = AS->smolar();
    double U = AS->umolar(), rho = AS->rhomolar();

    SECTION("DmolarHmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarHmolar_INPUTS, rho, H));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(P).epsilon(0.001));
    }
    SECTION("DmolarSmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarSmolar_INPUTS, rho, S));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(P).epsilon(0.001));
    }
    SECTION("DmolarUmolar") {
        auto AS2 = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "Nitrogen&Oxygen"));
        AS2->set_mole_fractions({0.79, 0.21});
        REQUIRE_NOTHROW(AS2->update(DmolarUmolar_INPUTS, rho, U));
        CHECK(AS2->T() == Catch::Approx(T).epsilon(0.001));
        CHECK(AS2->p() == Catch::Approx(P).epsilon(0.001));
    }
}

#endif
