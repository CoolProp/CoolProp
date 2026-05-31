#if defined(ENABLE_CATCH)

#    include "AbstractState.h"
#    include "DataStructures.h"
#    include "../Backends/Cubics/CubicBackend.h"
#    include <memory>
#    include <catch2/catch_all.hpp>
#    include "CoolPropTools.h"
#    include "CoolProp.h"

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

TEST_CASE("Michelsen Flash: Issue #2637 (PR mixture phase envelope)", "[michelsen][cubic][phase_envelope][2637]") {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("PR", "Methane&Ethane"));
    std::vector<double> z = {0.85, 0.15};
    AS->set_mole_fractions(z);

    // Phase envelope and critical points should work correctly
    CHECK_NOTHROW(AS->build_phase_envelope(""));
    std::vector<CoolProp::CriticalState> crit_points;
    CHECK_NOTHROW(crit_points = AS->all_critical_points());
    CHECK(crit_points.size() > 0);
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

#endif
