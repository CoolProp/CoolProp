/**
 * Tests for mixture flash routines: PT, PQ, and HSU_P.
 *
 * Tests all combinations of:
 *   - Backends:  SRK, PR, HEOS
 *   - Phases:    gas, liquid, two-phase
 *   - Modes:     PE (phase envelope), imposed, blind
 *
 * PT flash:   Q round-trip for two-phase; property self-consistency and
 *             cross-mode agreement for all phases.
 * PQ flash:   no-throw sweep; VLE consistency (fugacity equality + T
 *             monotonicity); phase-envelope integration.
 * HSU_P flash: PT → HP/SP/UP round-trips near and away from saturation.
 */

#if defined(ENABLE_CATCH)

#    include <cmath>
#    include <string>
#    include <vector>
#    include <memory>
#    include <algorithm>
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
//  PT FLASH - PER-MODE AND CROSS-MODE AGREEMENT
// ===================================================================

TEST_CASE("PT flash mixture per-mode", "[mixture][PT_flash]") {
    // Two-phase: Q in [0,1], Q round-trip via PQ, repeat-flash Q agreement.
    // Single-phase: Q==-1, repeat-flash density and enthalpy agreement.
    const double T_TOL = 0.01;    // K     - PQ round-trip temperature
    const double Q_TOL = 1e-4;    // -     - repeat-flash vapor quality
    const double H_TOL = 0.1;     // J/mol - repeat-flash enthalpy
    const double RHO_TOL = 1e-4;  // -     - relative density agreement

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
                            CHECK(std::abs(AS_pq->T() - cond.T) < T_TOL);

                            // Repeat PT flash recovers same Q
                            auto AS2 = make_state(be, mix, mode, cond.imposed);
                            AS2->update(PT_INPUTS, cond.P, cond.T);
                            CHECK(std::abs(AS2->Q() - Q) < Q_TOL);
                        } else {
                            // Single-phase: Q must be -1
                            CHECK(Q == -1.0);

                            // Repeat PT flash is self-consistent
                            auto AS2 = make_state(be, mix, mode, cond.imposed);
                            AS2->update(PT_INPUTS, cond.P, cond.T);
                            double rel_rho = std::abs(AS2->rhomolar() - rho) / std::max(std::abs(rho), 1e-10);
                            CHECK(rel_rho < RHO_TOL);
                            CHECK(std::abs(AS2->hmolar() - H) < H_TOL);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("PT flash mixture cross-mode consistency", "[mixture][PT_flash][cross-mode]") {
    // All three modes (PE, imposed, blind) must agree on Q (two-phase) or
    // density (single-phase) to within tolerance.
    const double Q_TOL = 1e-4;    // - vapor quality spread across modes
    const double RHO_TOL = 1e-3;  // - relative density spread across modes

    for (auto& mix : mixtures) {
        for (auto& be : test_backends) {
            for (auto& cond : mix.conditions) {

                std::string tag = be + " x-mode " + mix.fluids + " " + cond.label;
                DYNAMIC_SECTION(tag) {
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
                        CHECK(Q_max - Q_min < Q_TOL);
                    } else {
                        double rho_min = *std::min_element(rhos, rhos + 3);
                        double rho_max = *std::max_element(rhos, rhos + 3);
                        CHECK((rho_max - rho_min) / std::max(rho_max, 1e-10) < RHO_TOL);
                    }
                }
            }
        }
    }
}

// ===================================================================
//  HSU_P FLASH - ROUND-TRIP AND NEAR-SATURATION
// ===================================================================

TEST_CASE("HSU_P flash mixture HP round-trip", "[mixture][HSU_P_flash]") {

    const double TOL = 0.1;  // K

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
    const double TOL = 0.1;

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

TEST_CASE("HSU_P flash close to saturation for Propane&Butane", "[mixture][HSU_P_flash][saturation]") {
    // Propane&Butane 50/50 at P=1e5 Pa.
    const std::string be = "HEOS";
    const std::string fluids = "Propane&Butane";
    const std::vector<double> z = {0.5, 0.5};
    const double P = 1e5;
    const double TOL = 0.1;  // K

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

TEST_CASE("HSU_P flash close to saturation for Nitrogen&Methane&Ethane&Butane&Pentane", "[mixture][HSU_P_flash][saturation]") {
    // 5-component mixture: Nitrogen/Methane/Ethane/Butane/Pentane at P=3e5 Pa.
    const std::string be = "HEOS";
    const std::string fluids = "Nitrogen&Methane&Ethane&Butane&Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    const double P = 3e5;
    const double TOL = 0.1;  // K

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

// ===================================================================
//  VLE CONSISTENCY - PT AND PQ SWEEP
// ===================================================================

// Shared mixture definition for the VLE consistency tests
namespace {
constexpr double VLE_P = 3e5;         // Pa
constexpr int VLE_NPTS = 100;         // interior scan points
constexpr double VLE_Z_MIN = 1e-4;    // skip trace-phase components
constexpr double VLE_FUG_TOL = 1e-7;  // |f_L - f_V| / max(f_L, f_V)
constexpr double VLE_F_MIN = 1e-15;   // Pa - fugacity reference floor
constexpr double VLE_H_SLACK = 0.1;   // J/mol - allowed H regression
static const std::string VLE_FLUIDS = "Nitrogen&Methane&Ethane&Butane&Pentane";
static const std::vector<double> VLE_Z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
}  // namespace

TEST_CASE("PT flash 5-component N2-HC two-phase consistency", "[mixture][PT_flash]") {
    // Sweeps 100 temperatures uniformly from T_bub to T_dew.  Checks:
    //   1. phase == iphase_twophase at every point
    //   2. hmolar() non-decreasing with T (slack: VLE_H_SLACK J/mol)
    //   3. Fugacity equality for non-trace components (< VLE_FUG_TOL)
    const std::string fluids = VLE_FLUIDS;
    const std::vector<double> z = VLE_Z;
    const double P = VLE_P;
    const int NQ = VLE_NPTS;
    const std::size_t NC = z.size();

    for (const auto& be : test_backends) {
        DYNAMIC_SECTION(be) {
            // Bubble and dew temperatures for this backend
            auto sat = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            sat->set_mole_fractions(z);
            sat->update(PQ_INPUTS, P, 0.0);
            const double T_bub = sat->T();
            sat->update(PQ_INPUTS, P, 1.0);
            const double T_dew = sat->T();
            REQUIRE(T_dew > T_bub);

            // Reused states for per-component fugacity sub-calculations
            auto liq = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            auto vap = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);

            double H_prev = -1e100;
            for (int i = 0; i < NQ; ++i) {
                // Interior points - skip exact bubble/dew boundaries
                const double alpha = (static_cast<double>(i) + 0.5) / NQ;
                const double T = T_bub + (T_dew - T_bub) * alpha;
                CAPTURE(T);

                AS->update(PT_INPUTS, P, T);

                // 1. Phase
                CHECK(AS->phase() == iphase_twophase);

                // 2. Enthalpy monotonicity
                const double H = AS->hmolar();
                CHECK(H >= H_prev - VLE_H_SLACK);
                H_prev = H;

                // 3. Fugacity equality (only meaningful while phase is two-phase)
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
                        if (std::min(x[j], y[j]) < VLE_Z_MIN) continue;
                        const double f_liq = liq->fugacity_coefficient(j) * x[j];
                        const double f_vap = vap->fugacity_coefficient(j) * y[j];
                        const double f_ref = std::max(std::abs(f_liq), std::abs(f_vap));
                        if (f_ref > VLE_F_MIN) {
                            CAPTURE(j);
                            CHECK(std::abs(f_liq - f_vap) / f_ref < VLE_FUG_TOL);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("PQ flash 5-component N2-HC two-phase consistency", "[mixture][PQ_flash]") {
    // Sweeps 100 quality midpoints from Q=0 to Q=1.  Checks:
    //   1. phase == iphase_twophase at every point
    //   2. T non-decreasing with Q
    //   3. Fugacity equality for non-trace components (< VLE_FUG_TOL)
    const std::string fluids = VLE_FLUIDS;
    const std::vector<double> z = VLE_Z;
    const double P = VLE_P;
    const int NQ = VLE_NPTS;
    const std::size_t NC = z.size();

    for (const auto& be : test_backends) {
        DYNAMIC_SECTION(be) {
            // Bubble and dew temperatures bracket
            auto sat = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            sat->set_mole_fractions(z);
            sat->update(PQ_INPUTS, P, 0.0);
            const double T_bub = sat->T();
            sat->update(PQ_INPUTS, P, 1.0);
            const double T_dew = sat->T();
            REQUIRE(T_dew > T_bub);

            // Reused states for per-component fugacity sub-calculations
            auto liq = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            auto vap = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            auto AS = std::shared_ptr<AbstractState>(AbstractState::factory(be, fluids));
            AS->set_mole_fractions(z);

            double T_prev = -1e100;
            for (int i = 0; i < NQ; ++i) {
                // Interior midpoints - skip exact Q=0 and Q=1 boundaries
                const double Q = (static_cast<double>(i) + 0.5) / NQ;
                CAPTURE(Q);

                AS->update(PQ_INPUTS, P, Q);

                // 1. Phase
                CHECK(AS->phase() == iphase_twophase);

                // 2. Temperature monotonicity
                const double T = AS->T();
                CHECK(T >= T_prev);
                T_prev = T;

                // 3. Fugacity equality (only meaningful while phase is two-phase)
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
                        if (std::min(x[j], y[j]) < VLE_Z_MIN) continue;
                        const double f_liq = liq->fugacity_coefficient(j) * x[j];
                        const double f_vap = vap->fugacity_coefficient(j) * y[j];
                        const double f_ref = std::max(std::abs(f_liq), std::abs(f_vap));
                        if (f_ref > VLE_F_MIN) {
                            CAPTURE(j);
                            CHECK(std::abs(f_liq - f_vap) / f_ref < VLE_FUG_TOL);
                        }
                    }
                }
            }
        }
    }
}

// ===================================================================
//  PQ FLASH - PHASE ENVELOPE AND REGRESSION
// ===================================================================

TEST_CASE("PQ flash with built phase envelope - N2+CH4", "[mixture][PQ_flash][PhaseEnvelope]") {
    // Test that PQ flash works when the phase envelope is built.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen&Methane"));
    AS->set_mole_fractions({0.5, 0.5});
    AS->build_phase_envelope("");
    const std::size_t npts = AS->get_phase_envelope_data().T.size();
    // Ensure PE built
    CAPTURE(npts);
    CHECK(npts > 0);
    // Check built flag is actually true (no early exit)
    CHECK(AS->get_phase_envelope_data().built);
    // Calculate point inside PE
    REQUIRE_NOTHROW(AS->update(CoolProp::PQ_INPUTS, 1.5e5, 0.5));
    CHECK(AS->phase() == CoolProp::iphase_twophase);
}

TEST_CASE("PQ flash 6-component N2-HC mixture does not throw", "[mixture][PQ_flash]") {
    // Regression test: saturation_Wilson Secant diverges with preconditioning guess for q between 0.88 and 0.95.
    // New logic tries Brent first.
    const std::string fluids = "Nitrogen&Methane&Ethane&Propane&Butane&Pentane";
    const std::vector<double> z = {0.2936, 0.2720, 0.0592, 0.2932, 0.0787, 0.0033};
    const double P = 3.92e5;  // Pa
    const int NQ = 100;

    for (const auto& be : test_backends) {
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

#endif  // ENABLE_CATCH
