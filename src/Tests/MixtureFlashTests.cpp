/**
 * Comprehensive tests for mixture PT flash and HSU_P flash routines.
 *
 * Tests all combinations of:
 *   - Backends:  SRK, PR, HEOS
 *   - Phases:    gas, liquid, two-phase
 *   - Modes:     PE (phase envelope), imposed, blind
 *
 * PT flash: verifies Q round-trip for two-phase and property self-
 * consistency for single-phase, plus cross-mode agreement.
 *
 * HSU_P flash: round-trips PT → HP (and SP, UP for one exemplary case)
 * and checks that the recovered temperature matches the reference.
 */

#if defined(ENABLE_CATCH)

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <catch2/catch_all.hpp>

#include "AbstractState.h"
#include "DataStructures.h"
#include "CoolProp.h"

using namespace CoolProp;

// ===================================================================
// Mixture / condition table shared by all tests
// ===================================================================

struct PhaseCondition {
    std::string label;   // "gas", "liquid", "2ph"
    double P;            // Pa
    double T;            // K
    phases imposed;      // iphase_gas, iphase_liquid, iphase_twophase
};

struct MixtureDef {
    std::string fluids;
    std::vector<double> z;
    std::vector<PhaseCondition> conditions;
};

static const std::vector<MixtureDef> mixtures = {
    {"Methane&Ethane", {0.5, 0.5},
     {{"gas", 1e5, 300, iphase_gas},
      {"2ph", 1e5, 150, iphase_twophase},
      {"liquid", 1e5, 100, iphase_liquid}}},

    {"Methane&Propane", {0.7, 0.3},
     {{"gas", 1e5, 300, iphase_gas},
      {"2ph", 1e5, 150, iphase_twophase},
      {"liquid", 1e5, 100, iphase_liquid}}},

    {"Nitrogen&Oxygen", {0.79, 0.21},
     {{"gas", 1e5, 300, iphase_gas},
      {"2ph", 1e5, 80, iphase_twophase},
      {"liquid", 1e5, 75, iphase_liquid}}},

    {"Nitrogen&Methane&Ethane&Propane", {0.1, 0.5, 0.25, 0.15},
     {{"gas", 1e5, 300, iphase_gas},
      {"2ph", 1e5, 145, iphase_twophase},
      {"liquid", 1e5, 80, iphase_liquid}}},
};

static const std::vector<std::string> test_backends = {"SRK", "PR", "HEOS"};

enum class Mode { PE, Imposed, Blind };

static const std::vector<std::pair<std::string, Mode>> modes = {
    {"PE", Mode::PE},
    {"imposed", Mode::Imposed},
    {"blind", Mode::Blind},
};

// ===================================================================
// Helpers
// ===================================================================

/// Create an AbstractState with the given mode applied.
static std::shared_ptr<AbstractState>
make_state(const std::string& backend, const MixtureDef& mix,
           Mode mode, phases imposed_phase)
{
    auto AS = std::shared_ptr<AbstractState>(
        AbstractState::factory(backend, mix.fluids));
    AS->set_mole_fractions(mix.z);
    switch (mode) {
        case Mode::PE:
            AS->build_phase_envelope("dummy");
            break;
        case Mode::Imposed:
            AS->specify_phase(imposed_phase);
            break;
        case Mode::Blind:
            break;
    }
    return AS;
}

// ===================================================================
//  PT FLASH TESTS
// ===================================================================

TEST_CASE("PT flash mixture per-mode", "[mixture][PT_flash]") {

    for (auto& mix : mixtures) {
        for (auto& be : test_backends) {
            for (auto& cond : mix.conditions) {
                for (auto& [mode_name, mode] : modes) {

                    std::string tag = be + " " + mode_name + " "
                                    + mix.fluids + " " + cond.label;
                    DYNAMIC_SECTION(tag) {
                        auto AS = make_state(be, mix, mode, cond.imposed);
                        AS->update(PT_INPUTS, cond.P, cond.T);

                        double Q   = AS->Q();
                        double H   = AS->hmolar();
                        double rho = AS->rhomolar();

                        if (cond.label == "2ph") {
                            // Q must be in [0, 1]
                            REQUIRE(Q >= 0.0);
                            REQUIRE(Q <= 1.0);

                            // Q round-trip via PQ → T → PT → Q
                            auto AS_pq = std::shared_ptr<AbstractState>(
                                AbstractState::factory(be, mix.fluids));
                            AS_pq->set_mole_fractions(mix.z);
                            AS_pq->update(PQ_INPUTS, cond.P, Q);
                            double T_from_pq = AS_pq->T();
                            CHECK(std::abs(T_from_pq - cond.T) < 1.0);

                            // Repeat PT flash recovers same Q
                            auto AS2 = make_state(be, mix, mode, cond.imposed);
                            AS2->update(PT_INPUTS, cond.P, cond.T);
                            CHECK(std::abs(AS2->Q() - Q) < 0.01);
                        } else {
                            // Single-phase: Q must be -1
                            CHECK(Q == -1.0);

                            // Repeat PT flash is self-consistent
                            auto AS2 = make_state(be, mix, mode, cond.imposed);
                            AS2->update(PT_INPUTS, cond.P, cond.T);
                            double rho2 = AS2->rhomolar();
                            double H2   = AS2->hmolar();
                            double rel_rho = std::abs(rho2 - rho)
                                           / std::max(std::abs(rho), 1e-10);
                            CHECK(rel_rho < 0.01);
                            CHECK(std::abs(H2 - H) < 1.0);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("PT flash mixture cross-mode consistency", "[mixture][PT_flash][cross-mode]") {

    for (auto& mix : mixtures) {
        for (auto& be : test_backends) {
            for (auto& cond : mix.conditions) {

                std::string tag = be + " x-mode " + mix.fluids
                                + " " + cond.label;
                DYNAMIC_SECTION(tag) {
                    // Evaluate all three modes
                    double Qs[3], rhos[3];
                    int idx = 0;
                    for (auto& [mode_name, mode] : modes) {
                        auto AS = make_state(be, mix, mode, cond.imposed);
                        AS->update(PT_INPUTS, cond.P, cond.T);
                        Qs[idx]   = AS->Q();
                        rhos[idx] = AS->rhomolar();
                        idx++;
                    }

                    if (cond.label == "2ph") {
                        double Q_min = *std::min_element(Qs, Qs + 3);
                        double Q_max = *std::max_element(Qs, Qs + 3);
                        CHECK(Q_max - Q_min < 0.01);
                    } else {
                        double rho_min = *std::min_element(rhos, rhos + 3);
                        double rho_max = *std::max_element(rhos, rhos + 3);
                        double spread  = (rho_max - rho_min)
                                       / std::max(rho_max, 1e-10);
                        CHECK(spread < 0.01);
                    }
                }
            }
        }
    }
}

// ===================================================================
//  HSU_P FLASH TESTS
// ===================================================================

TEST_CASE("HSU_P flash mixture HP round-trip", "[mixture][HSU_P_flash]") {

    const double TOL = 1.0;  // K

    for (auto& mix : mixtures) {
        for (auto& be : test_backends) {
            for (auto& cond : mix.conditions) {
                for (auto& [mode_name, mode] : modes) {

                    std::string tag = be + " " + mode_name + " "
                                    + mix.fluids + " HP " + cond.label;
                    DYNAMIC_SECTION(tag) {
                        // 1. Reference state via PT
                        auto AS_ref = make_state(be, mix, mode, cond.imposed);
                        AS_ref->update(PT_INPUTS, cond.P, cond.T);
                        double H_ref = AS_ref->hmolar();

                        // 2. HP flash round-trip
                        auto AS = make_state(be, mix, mode, cond.imposed);
                        AS->update(HmolarP_INPUTS, H_ref, cond.P);
                        double T_got = AS->T();
                        CHECK(std::abs(T_got - cond.T) < TOL);
                    }
                }
            }
        }
    }
}

TEST_CASE("HSU_P flash mixture SP and UP round-trip",
          "[mixture][HSU_P_flash]") {

    // Exemplary case: SRK, blind, Methane&Ethane, two-phase
    const std::string be = "SRK";
    const MixtureDef& mix = mixtures[0];  // Methane&Ethane
    const PhaseCondition& cond = mix.conditions[1];  // 2ph
    const double TOL = 1.0;

    auto AS_ref = make_state(be, mix, Mode::Blind, cond.imposed);
    AS_ref->update(PT_INPUTS, cond.P, cond.T);

    SECTION("SP round-trip") {
        double S_ref = AS_ref->smolar();
        auto AS = make_state(be, mix, Mode::Blind, cond.imposed);
        AS->update(PSmolar_INPUTS, cond.P, S_ref);
        CHECK(std::abs(AS->T() - cond.T) < TOL);
    }

    SECTION("UP round-trip") {
        double U_ref = AS_ref->umolar();
        auto AS = make_state(be, mix, Mode::Blind, cond.imposed);
        AS->update(PUmolar_INPUTS, cond.P, U_ref);
        CHECK(std::abs(AS->T() - cond.T) < TOL);
    }
}

#endif // ENABLE_CATCH
