// Blind mixture PQ/QT flash at interior quality (GH #3585).
//
// The blind (no phase envelope) mixture QT and PQ flashes solve interior Q with
// newton_raphson_twophase, which enforces z = (1-Q) x + Q y.  When that solve fails they fall
// back to newton_raphson_saturation, a bubble/dew solver with no mass-balance condition, whose
// split can be a genuine VLE state of a DIFFERENT overall composition.
//
// Two defects combined in #3585:
//   * the two-phase solve failed far more often than it should: its iso-fugacity residual went
//     through ln Z = log(1 + delta dalphar/ddelta), which cancels in a low-pressure liquid
//     (Z ~ 2e-6 for R436A at 140 K, 37 Pa) and left ~1e-8 of noise, above the solver's 1e-9
//     convergence gate, so it ran out of iterations after Newton had already converged;
//   * the fallback state was returned whether or not it balanced.
// The first is fixed by evaluating the residual as ln f^V - ln f^L straight from (T, rho, x);
// the second by refusing a fallback state that does not satisfy the mass balance.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/CoolProp.h"
#    include "CoolProp/DataStructures.h"

#    include <algorithm>
#    include <cmath>
#    include <limits>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

/// max_i |(1-Q) x_i + Q y_i - z_i| from the flash's published output, NaN-propagating.
double mass_balance_error(AbstractState& AS) {
    const std::vector<double> z = AS.get_mole_fractions();
    const std::vector<double> x = AS.mole_fractions_liquid();
    const std::vector<double> y = AS.mole_fractions_vapor();
    REQUIRE(x.size() == z.size());
    REQUIRE(y.size() == z.size());
    const double Q = AS.Q();
    double err = 0;
    for (std::size_t i = 0; i < z.size(); ++i) {
        const double e = std::abs((1 - Q) * x[i] + Q * y[i] - z[i]);
        if (!std::isfinite(e)) return std::numeric_limits<double>::infinity();  // a NaN must fail every check
        if (e > err) err = e;
    }
    return err;
}

/// Reading the warning slot clears it.
std::string take_warning() {
    return get_global_param_string("warnstring");
}

}  // namespace

// The issue's reproducer: two neighbouring floats of T.  The first used to exhaust the two-phase
// solver's iterations and return the bubble/dew fallback (mass balance off by 1e-2, p off by
// 2.8 %); the second converged.  Both must now converge, to the same state.
TEST_CASE("Mixture QT flash: R436A at neighbouring floats of T (#3585)", "[flash][QT_flash][twophase][massbalance]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "R436A.mix"));
    std::vector<double> p;
    for (double T : {140.00109271678056, 140.00109271678062}) {
        CAPTURE(T);
        take_warning();
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.4, T));
        // An empty warning slot proves the two-phase solver answered, not the fallback.
        CHECK(take_warning().empty());
        CHECK(mass_balance_error(*AS) < 1e-12);
        p.push_back(AS->p());
    }
    CAPTURE(p[0], p[1]);
    CHECK(p[0] == Catch::Approx(p[1]).epsilon(1e-10));
}

// The issue's sweep: R436A at Q = 0.4 from 100 to 200 K.  156 of 400 states fell back on
// 8.1.0 and every fallback was wrong.  The cold end of this range is the low-pressure liquid
// (tens of Pa) where ln Z cancels.
TEST_CASE("Mixture QT flash: R436A cold-end sweep converges without fallback (#3585)", "[flash][QT_flash][twophase][massbalance]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "R436A.mix"));
    const int n = 101;
    for (int k = 0; k < n; ++k) {
        const double T = 100.0 + 100.0 * k / (n - 1);
        CAPTURE(T);
        take_warning();
        REQUIRE_NOTHROW(AS->update(QT_INPUTS, 0.4, T));
        CHECK(take_warning().empty());
        CHECK(mass_balance_error(*AS) < 1e-12);
    }
}

// The PQ flash uses the same two-phase solver with p imposed, and the same low-pressure
// liquid.  R436A at Q = 0.4 from 10 Pa to 1 MPa.
TEST_CASE("Mixture PQ flash: R436A low-pressure sweep converges without fallback (#3585)", "[flash][PQ_flash][twophase][massbalance]") {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", "R436A.mix"));
    const int n = 51;
    for (int k = 0; k < n; ++k) {
        const double p = 10.0 * std::pow(1e5, static_cast<double>(k) / (n - 1));
        CAPTURE(p);
        take_warning();
        REQUIRE_NOTHROW(AS->update(PQ_INPUTS, p, 0.4));
        CHECK(take_warning().empty());
        CHECK(mass_balance_error(*AS) < 1e-12);
    }
}

// Whatever the blind QT flash returns at interior Q must balance: either the two-phase solver
// answered, or the fallback happened to land on the right mixture, or the flash threw.  The
// blends are the 32 the issue lists plus the five whose remaining fallbacks (near the top of the
// envelope, and HighCO2 below CO2's triple point) returned states off by up to 0.15 before the
// fallback was gated.  Temperatures span each blend's phase envelope, built on a separate
// object so the flash under test stays on the blind path.
TEST_CASE("Mixture QT flash never returns an unbalanced interior-Q state (#3585)", "[flash][QT_flash][massbalance]") {
    const std::vector<std::string> blends = {
      "Ekofisk", "HighCO2", "R411A", "R411B", "R415A", "R415B", "R418A", "R430A", "R431A", "R432A",     "R433A",  "R433B", "R433C",
      "R436A",   "R436B",   "R436C", "R440A", "R441A", "R443A", "R454A", "R454B", "R454C", "R457A",     "R457B",  "R457C", "R465A",
      "R500",    "R501",    "R509A", "R510A", "R511A", "R512A", "R421B", "R422A", "R434A", "GulfCoast", "HighN2", "R507A", "NaturalGasSample",
    };
    int n_returned = 0;
    for (const auto& blend : blends) {
        const std::string name = blend + ".mix";
        CAPTURE(name);
        auto env = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", name));
        env->build_phase_envelope("");
        const auto& PE = env->get_phase_envelope_data();
        // Only the temperature span is used; `built` is false for some blends (Ekofisk) whose
        // envelope trace still covers the dome.
        REQUIRE(PE.T.size() > 10);
        const double Tlo = *std::min_element(PE.T.begin(), PE.T.end());
        const double Thi = *std::max_element(PE.T.begin(), PE.T.end());

        auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", name));
        const int n = 39;
        for (int k = 0; k < n; ++k) {
            const double T = Tlo + (Thi - Tlo) * k / (n - 1);
            for (double Q : {0.1, 0.5, 0.9}) {
                CAPTURE(T, Q);
                try {
                    AS->update(QT_INPUTS, Q, T);
                } catch (const std::exception&) {
                    continue;  // refusing is allowed; returning the wrong mixture is not
                }
                ++n_returned;
                CHECK(mass_balance_error(*AS) <= 1e-10);
            }
        }
    }
    take_warning();
    // Guard against the test going vacuous by refusing too much: 4255 of the 4563 flashes
    // (39 blends x 117) return; the rest throw, nearly all of them independently of #3585.
    CAPTURE(n_returned);
    CHECK(n_returned > 4150);
}

#endif
