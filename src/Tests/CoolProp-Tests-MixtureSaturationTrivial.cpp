// Regression tests for two silent wrong answers from the mixture flash routines.
//
// Runs in the default suite (tag [mixsat], NOT [.]-hidden).  Run explicitly:
//   ./CatchTestRunner "[mixsat]"
//
// Both were found while building the isoline tracer for the Python plotting
// package (GH #3344, discussion #3269), where they surfaced as excursions and
// gaps in mixture property plots.  Neither raised an error at the time: each
// returned a plausible-looking state that was simply not the one asked for.
//
//  GH #3346 / bd CoolProp-mojf -- QT_flash and PQ_flash could converge on the
//      trivial solution near the critical point of a mixture, returning a state
//      whose saturated liquid and vapour are the same root.  Successive
//      substitution and the Newton-Raphson saturation solver both satisfy their
//      residual when every K-value reaches unity, and near the critical point
//      that is where they land.
//
//  bd CoolProp-1gth -- DHSU_T_flash excluded density from its own
//      "did this converge to what was asked for" check, on the reasoning that
//      D+T is a direct evaluation.  True of the fast path, but the P-sweep
//      fallback solves rho(P) = value through PT flashes and can settle on a
//      different density entirely.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>
#    include <catch2/catch_approx.hpp>

#    include "AbstractState.h"
#    include "CoolProp.h"

#    include <cmath>
#    include <memory>
#    include <string>

using namespace CoolProp;

TEST_CASE("QT_flash does not return the trivial solution for a mixture", "[mixsat]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "R513A.mix"));

    SECTION("a genuine two-phase point is unaffected") {
        // 364.05 K is comfortably below the critical region: the two phases are
        // separated by a factor of about 2.6 in density.
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.5, 364.0501));
        const double rho_liq = AS->saturated_liquid_keyed_output(iDmolar);
        const double rho_vap = AS->saturated_vapor_keyed_output(iDmolar);
        CHECK(rho_liq > rho_vap * 1.5);
        CHECK(AS->smass() == Catch::Approx(1581.676).epsilon(1e-4));
    }

    SECTION("a collapsed pair is reported rather than returned") {
        // At these temperatures the solver used to return rho_liq == rho_vap to
        // seven digits, which put the Q=0.5 line 93 J/kg/K out of place.
        for (double T : {367.4330, 368.0}) {
            CAPTURE(T);
            bool threw = false;
            try {
                AS->update(QT_INPUTS, 0.5, T);
            } catch (const std::exception&) {
                threw = true;
            }
            if (!threw) {
                // Converging here is fine, provided the two phases are distinct.
                const double rho_liq = AS->saturated_liquid_keyed_output(iDmolar);
                const double rho_vap = AS->saturated_vapor_keyed_output(iDmolar);
                CHECK(std::abs(rho_liq - rho_vap) > 1e-4 * std::max(rho_liq, rho_vap));
            }
        }
    }
}

TEST_CASE("DmassT_INPUTS returns the density it was given for a mixture", "[mixsat]") {
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "R513A.mix"));
    const double T = 288.853439;

    // 1198.85 and 1200.0 used to come back as 617.6 and 1256.6 kg/m3, the first
    // of them with h = -1.02e7 J/kg, while their neighbours round tripped exactly.
    for (double rho : {1190.0, 1195.0, 1198.85, 1200.0, 1210.0}) {
        CAPTURE(rho);
        bool threw = false;
        try {
            AS->update(DmassT_INPUTS, rho, T);
        } catch (const std::exception&) {
            threw = true;  // refusing is acceptable; answering wrongly is not
        }
        if (!threw) {
            CHECK(AS->rhomass() == Catch::Approx(rho).epsilon(1e-9));
            // and the enthalpy has to be a compressed-liquid value, not -1e7
            CHECK(AS->hmass() > 0.0);
            CHECK(AS->hmass() < 1e6);
        }
    }
}

TEST_CASE("the trivial-solution guard does not disturb phase envelope tracing", "[mixsat]") {
    // The guard sits at the QT/PQ flash boundary rather than inside the
    // saturation solvers precisely because PhaseEnvelopeRoutines drives those
    // same solvers up to and through the critical point, where the two phases
    // collapsing together is the answer it is looking for.  If that reasoning
    // were wrong, envelope construction would be the first thing to break.
    for (std::string fluid :
         {"R513A.mix", "R410A.mix", "R404A.mix", "R407C.mix", "Air.mix", "R454B.mix", "R448A.mix", "R441A.mix", "R504.mix", "R507A.mix"}) {
        CAPTURE(fluid);
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", fluid));
        REQUIRE_NOTHROW(AS->build_phase_envelope("none"));
        CHECK(AS->get_phase_envelope_data().T.size() > 10);
    }
}

TEST_CASE("imposing the phase still gives the correct compressed liquid", "[mixsat]") {
    // The EOS always had the right answer; the defect was in the free path's
    // phase determination.  This pins the reference the fix is measured against.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "R513A.mix"));
    AS->specify_phase(iphase_liquid);
    AS->update(DmassT_INPUTS, 1198.85, 288.853439);
    CHECK(AS->rhomass() == Catch::Approx(1198.85).epsilon(1e-9));
    CHECK(AS->p() == Catch::Approx(6.4474e6).epsilon(1e-3));
    CHECK(AS->hmass() == Catch::Approx(224455.0).epsilon(1e-3));
}

#endif
