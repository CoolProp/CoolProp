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

#    include <cmath>
#    include <string>
#    include <vector>
#    include <memory>
#    include <catch2/catch_all.hpp>

#    include "AbstractState.h"
#    include "DataStructures.h"
#    include "CoolProp.h"

using namespace CoolProp;

// ===================================================================
// Mixture / condition table shared by all tests
// ===================================================================

struct PhaseCondition
{
    std::string label;  // "gas", "liquid", "2ph"
    double P;           // Pa
    double T;           // K
    phases imposed;     // iphase_gas, iphase_liquid, iphase_twophase
};

struct MixtureDef
{
    std::string fluids;
    std::vector<double> z;
    std::vector<PhaseCondition> conditions;
};

static const std::vector<MixtureDef> mixtures = {
  {"Methane&Ethane", {0.5, 0.5}, {{"gas", 1e5, 300, iphase_gas}, {"2ph", 1e5, 150, iphase_twophase}, {"liquid", 1e5, 100, iphase_liquid}}},

  {"Methane&Propane", {0.7, 0.3}, {{"gas", 1e5, 300, iphase_gas}, {"2ph", 1e5, 150, iphase_twophase}, {"liquid", 1e5, 100, iphase_liquid}}},

  {"Nitrogen&Oxygen", {0.79, 0.21}, {{"gas", 1e5, 300, iphase_gas}, {"2ph", 1e5, 80, iphase_twophase}, {"liquid", 1e5, 75, iphase_liquid}}},

  {"Nitrogen&Methane&Ethane&Propane",
   {0.1, 0.5, 0.25, 0.15},
   {{"gas", 1e5, 300, iphase_gas}, {"2ph", 1e5, 145, iphase_twophase}, {"liquid", 1e5, 80, iphase_liquid}}},
};

static const std::vector<std::string> test_backends = {"SRK", "PR", "HEOS"};

enum class Mode
{
    PE,
    Imposed,
    Blind
};

static const std::vector<std::pair<std::string, Mode>> modes = {
  {"PE", Mode::PE},
  {"imposed", Mode::Imposed},
  {"blind", Mode::Blind},
};

// ===================================================================
// Helpers
// ===================================================================

/// Create an AbstractState with the given mode applied.
static std::shared_ptr<AbstractState> make_state(const std::string& backend, const MixtureDef& mix, Mode mode, phases imposed_phase) {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(backend, mix.fluids));
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

                    std::string tag = be + " " + mode_name + " " + mix.fluids + " " + cond.label;
                    DYNAMIC_SECTION(tag) {
                        auto AS = make_state(be, mix, mode, cond.imposed);
                        AS->update(PT_INPUTS, cond.P, cond.T);

                        double Q = AS->Q();
                        double H = AS->hmolar();
                        double rho = AS->rhomolar();

                        if (cond.label == "2ph") {
                            // Q must be in [0, 1]
                            REQUIRE(Q >= 0.0);
                            REQUIRE(Q <= 1.0);

                            // Q round-trip via PQ → T → PT → Q
                            auto AS_pq = std::shared_ptr<AbstractState>(AbstractState::factory(be, mix.fluids));
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
                            double H2 = AS2->hmolar();
                            double rel_rho = std::abs(rho2 - rho) / std::max(std::abs(rho), 1e-10);
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

                std::string tag = be + " x-mode " + mix.fluids + " " + cond.label;
                DYNAMIC_SECTION(tag) {
                    // Evaluate all three modes
                    double Qs[3], rhos[3];
                    int idx = 0;
                    for (auto& [mode_name, mode] : modes) {
                        auto AS = make_state(be, mix, mode, cond.imposed);
                        AS->update(PT_INPUTS, cond.P, cond.T);
                        Qs[idx] = AS->Q();
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
                        double spread = (rho_max - rho_min) / std::max(rho_max, 1e-10);
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

                    std::string tag = be + " " + mode_name + " " + mix.fluids + " HP " + cond.label;
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

TEST_CASE("HSU_P flash mixture SP and UP round-trip", "[mixture][HSU_P_flash]") {

    // Exemplary case: SRK, blind, Methane&Ethane, two-phase
    const std::string be = "SRK";
    const MixtureDef& mix = mixtures[0];             // Methane&Ethane
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

// ===================================================================
//  NEAR-SATURATION AND CONSISTENCY TESTS
// ===================================================================

TEST_CASE("HSU_P flash close to saturation for Propane&Butane", "[mixture][HSU_P_flash][saturation]") {
    // Propane&Butane 50/50 at P=1e5 Pa.  The stability analysis returns
    // gas from ~0.7 K below T_bubble, so the PQ-routed liquid bracket
    // with phase imposition is essential for correctness.  Tests cover
    // HP, SP, and UP for all five regions including T_bub-0.5 K which
    // falls inside the false-gas zone.
    const std::string be = "HEOS";
    const std::string fluids = "Propane&Butane";
    const std::vector<double> z = {0.5, 0.5};
    const double P = 1e5;
    const double TOL = 1e-6;  // K

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
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
    const std::vector<Case> cases = {
      {"subcooled T_bub-0.5", T_bub - 0.5, iphase_liquid}, {"subcooled T_bub-1", T_bub - 1.0, iphase_liquid},
      {"two-phase T_bub+1", T_bub + 1.0, iphase_twophase}, {"two-phase midpoint", 0.5 * (T_bub + T_dew), iphase_twophase},
      {"two-phase T_dew-1", T_dew - 1.0, iphase_twophase}, {"superheated T_dew+1", T_dew + 1.0, iphase_gas},
    };

    for (auto& c : cases) {
        // Reference properties from a phase-imposed PT flash
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
        ref->set_mole_fractions(z);
        ref->specify_phase(c.ph);
        ref->update(PT_INPUTS, P, c.T);
        const double H_ref = ref->hmass();
        const double S_ref = ref->smass();
        const double U_ref = ref->umass();
        ref->unspecify_phase();

        DYNAMIC_SECTION(c.label + " HP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);
            AS->update(HmassP_INPUTS, H_ref, P);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " SP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);
            AS->update(PSmass_INPUTS, P, S_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
        DYNAMIC_SECTION(c.label + " UP") {
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);
            AS->update(PUmass_INPUTS, P, U_ref);
            CHECK(std::abs(AS->T() - c.T) < TOL);
        }
    }
}

TEST_CASE("HSU_P flash catches EOS inconsistency in 5-component mixture", "[mixture][HSU_P_flash][consistency]") {
    // Nitrogen&Methane&Ethane&Propane&Butane at P=1e5 Pa exhibits two
    // known EOS limitations (possibly due to VLLE) detected by HSU_P_flash_mixtures
    // Throw instead of returning unphysical results:
    //   * subcooled (T_bub-1): mismatch of PQ and phase imposed PT at bubble point indicates bad liquid region
    //     -> reclassified to two-phase but cannot be bracketed -> ValueError.
    //   * two-phase1 (T_bub+1): PT solution is not confirmed by PQ re-verification
    //     -> ValueError.
    const std::string be = "HEOS";
    const std::string fluids = "Nitrogen&Methane&Ethane&Propane&Butane";
    const std::vector<double> z = {0.1, 0.2, 0.3, 0.2, 0.2};
    const double P = 1e5;

    auto sat = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
    sat->set_mole_fractions(z);
    sat->update(PQ_INPUTS, P, 0.0);
    const double T_bub = sat->T();

    SECTION("subcooled detects bubble-point inconsistency") {
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
        ref->set_mole_fractions(z);
        ref->specify_phase(iphase_liquid);
        ref->update(PT_INPUTS, P, T_bub - 1.0);
        const double H_ref = ref->hmass();
        ref->unspecify_phase();

        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
        AS->set_mole_fractions(z);
        REQUIRE_THROWS_WITH(AS->update(HmassP_INPUTS, H_ref, P), Catch::Matchers::ContainsSubstring("bubble-point inconsistency detected"));
    }

    SECTION("two-phase1 catches PQ verification failure") {
        auto ref = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
        ref->set_mole_fractions(z);
        ref->specify_phase(iphase_twophase);
        ref->update(PT_INPUTS, P, T_bub + 1.0);
        const double H_ref = ref->hmass();
        ref->unspecify_phase();

        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
        AS->set_mole_fractions(z);
        REQUIRE_THROWS_WITH(AS->update(HmassP_INPUTS, H_ref, P), Catch::Matchers::ContainsSubstring("not confirmed by PQ flash"));
    }
}

#endif  // ENABLE_CATCH
