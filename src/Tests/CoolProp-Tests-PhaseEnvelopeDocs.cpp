// Every phase-envelope example that appears in the documentation must actually build.
//
// The docs are executed at build time (the .plot:: and .ipython:: directives run the code), so a
// mixture whose envelope silently produces nothing turns into a blank or truncated figure rather
// than an error.  These cases are transcribed from:
//
//   Web/fluid_properties/methane-ethane.py   six methane/ethane isopleths, level "dummy"
//   Web/coolprop/Cubics.rst                  SRK with two kij values, and the SRK-transformed
//                                            multi-fluid model, at a 1e4 Pa starting pressure
//   Web/coolprop/REFPROP.rst                 R32/R125 through CoolProp's own routines
//   Web/coolprop/GERG.rst                    GERG2008 methane/ethane, documented as 222 points
//
// If a documented example is changed, changed here too; if one starts failing, the figure in the
// docs is already wrong.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/DataStructures.h"
#    include "Backends/Helmholtz/PhaseEnvelopeTracers.h"

#    include "CoolProp/detail/tools.h"
#    include <cmath>
#    include <iostream>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

/// Restores the configured starting pressure when a section ends, pass or fail.
struct StartingPressureGuard
{
    explicit StartingPressureGuard(double p) : saved(get_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA)) {
        set_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA, p);
    }
    ~StartingPressureGuard() {
        set_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA, saved);
    }
    StartingPressureGuard(const StartingPressureGuard&) = delete;
    StartingPressureGuard& operator=(const StartingPressureGuard&) = delete;
    StartingPressureGuard(StartingPressureGuard&&) = delete;
    StartingPressureGuard& operator=(StartingPressureGuard&&) = delete;

   private:
    double saved;
};

/// A documented envelope has to be plottable: enough points to draw, all finite, and a pressure
/// range wide enough that a log-scale figure shows a curve rather than a dot.
void check_plottable(const PhaseEnvelopeData& env, const std::string& what, std::size_t min_points = 30) {
    std::cout << format("[docs] %-46s n=%4d closed=%d stop=%s\n", what.c_str(), static_cast<int>(env.T.size()), static_cast<int>(env.closed),
                        PhaseEnvelopeTracers::last_stop_reason().c_str());
    REQUIRE(env.built);
    REQUIRE(env.T.size() >= min_points);
    REQUIRE(env.p.size() == env.T.size());
    double pmin = env.p[0], pmax = env.p[0];
    for (std::size_t k = 0; k < env.T.size(); ++k) {
        CAPTURE(k, env.T[k], env.p[k]);
        REQUIRE(std::isfinite(env.T[k]));
        REQUIRE(std::isfinite(env.p[k]));
        REQUIRE(env.T[k] > 0);
        REQUIRE(env.p[k] > 0);
        for (const auto& xj : env.x) {
            REQUIRE(std::isfinite(xj[k]));
        }
        pmin = std::min(pmin, env.p[k]);
        pmax = std::max(pmax, env.p[k]);
    }
    CAPTURE(pmin, pmax, env.T.size());
    CHECK(pmax / pmin > 100);  // spans more than two decades, i.e. it is a curve
}

}  // namespace

TEST_CASE("Docs: methane/ethane isopleths from methane-ethane.py", "[phase_envelope][docs]") {
    // Web/fluid_properties/methane-ethane.py builds one envelope per composition and swallows
    // ValueError.  Before the start-pressure retry and the partial-envelope change, a failure
    // here was invisible: the script caught the error and simply plotted one fewer curve.
    for (double x0 : {0.02, 0.2, 0.4, 0.6, 0.8, 0.98}) {
        SECTION("x(methane) = " + std::to_string(x0)) {
            std::shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS", "Methane&Ethane"));
            HEOS->set_mole_fractions({x0, 1 - x0});
            REQUIRE_NOTHROW(HEOS->build_phase_envelope("dummy"));
            const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
            check_plottable(env, "methane-ethane.py x=" + std::to_string(x0));
            // A methane/ethane isopleth is a closed loop through the critical point for every one
            // of these compositions; anything else means the figure is missing a branch.
            CHECK(env.closed);
        }
    }
}

TEST_CASE("Docs: cubic and SRK-transformed envelopes from Cubics.rst", "[phase_envelope][docs]") {
    // The doc raises the starting pressure first: "behavior at very low pressure is problematic".
    StartingPressureGuard guard(1e4);

    SECTION("SRK with kij") {
        for (double kij : {0.0, 0.1}) {
            SECTION("kij = " + std::to_string(kij)) {
                std::shared_ptr<AbstractState> SRK(AbstractState::factory("SRK", "Methane&Ethane"));
                SRK->set_mole_fractions({0.5, 0.5});
                SRK->set_binary_interaction_double(0, 1, "kij", kij);
                REQUIRE_NOTHROW(SRK->build_phase_envelope(""));
                check_plottable(SRK->get_phase_envelope_data(), "Cubics.rst SRK kij=" + std::to_string(kij));
            }
        }
    }
    SECTION("SRK transformations inside the multi-fluid model") {
        std::shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS", "Methane-SRK&Ethane-SRK"));
        HEOS->set_mole_fractions({0.5, 0.5});
        REQUIRE_NOTHROW(HEOS->build_phase_envelope("none"));
        check_plottable(HEOS->get_phase_envelope_data(), "Cubics.rst Methane-SRK&Ethane-SRK");
    }
}

TEST_CASE("Docs: R32/R125 envelope from REFPROP.rst", "[phase_envelope][docs]") {
    // Only the CoolProp half; the REFPROP half of that example exercises REFPROP's own routine.
    std::shared_ptr<AbstractState> HEOS(AbstractState::factory("HEOS", "R32&R125"));
    HEOS->set_mole_fractions({0.5, 0.5});
    REQUIRE_NOTHROW(HEOS->build_phase_envelope(""));
    const PhaseEnvelopeData& env = HEOS->get_phase_envelope_data();
    check_plottable(env, "REFPROP.rst R32&R125");
    CHECK(env.closed);
}

TEST_CASE("Docs: GERG2008 methane/ethane envelope from GERG.rst", "[phase_envelope][docs]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", "Methane&Ethane"));
    AS->set_mole_fractions({0.9, 0.1});
    REQUIRE_NOTHROW(AS->build_phase_envelope(""));
    const PhaseEnvelopeData& env = AS->get_phase_envelope_data();
    check_plottable(env, "GERG.rst GERG2008 Methane&Ethane", 100);  // the doc quotes 222 points
    CHECK(env.closed);
}

#endif  // ENABLE_CATCH
