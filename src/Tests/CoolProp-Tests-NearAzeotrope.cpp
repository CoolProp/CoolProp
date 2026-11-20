// PT flash inside the two-phase dome of NEAR-AZEOTROPIC refrigerant blends (COO-193).
//
// For a near-azeotrope the incipient phase is close to the feed and its tangent-plane distance is shallow.
// The stability test's successive-substitution proximity heuristic ("trial converging to the trivial
// solution") also fires on such a trial.  When that heuristic alone was taken as proof of stability
// (#3427), the PT flash reported single-phase for most dome states of R410A and R507A.  Each state here is
// put on the dome with a QT flash, and the PT flash at that same (T, p) must find the split with the same
// quality.  The PT flash only reaches the stability test because no phase envelope is built here (a built,
// closed envelope bypasses it in PT_flash_mixtures).

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"

#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

struct Blend
{
    const char* name;
    const char* fluids;
    std::vector<double> mass_fractions;
    double T_min, T_max;  // [K]
};

}  // namespace

TEST_CASE("Near-azeotrope: PT flash finds the split at every QT dome state", "[flash][mixture]") {
    // T ranges stay below the blends' critical regions.  R507A below 235 K is left out: some of those
    // states were not recovered before #3427 either (COO-195).
    const std::vector<Blend> blends = {
      {"R410A", "R32&R125", {0.5, 0.5}, 215.0, 305.0},
      {"R454B", "R32&R1234yf", {0.689, 0.311}, 215.0, 305.0},
      {"R404A", "R125&R143a&R134a", {0.44, 0.52, 0.04}, 215.0, 305.0},
      {"R507A", "R125&R143a", {0.5, 0.5}, 235.0, 305.0},
    };
    for (const auto& b : blends) {
        DYNAMIC_SECTION(b.name) {
            std::unique_ptr<AbstractState> qt(AbstractState::factory("HEOS", b.fluids));
            std::unique_ptr<AbstractState> pt(AbstractState::factory("HEOS", b.fluids));
            qt->set_mass_fractions(b.mass_fractions);
            pt->set_mass_fractions(b.mass_fractions);
            int n = 0;
            for (int k = 0; b.T_min + 10.0 * k <= b.T_max + 1e-9; ++k) {
                const double T = b.T_min + 10.0 * k;
                for (double Q : {0.05, 0.25, 0.5, 0.75, 0.95}) {
                    CAPTURE(T, Q);
                    qt->update(QT_INPUTS, Q, T);
                    REQUIRE(qt->phase() == iphase_twophase);
                    const double p = qt->p();
                    CAPTURE(p);
                    pt->update(PT_INPUTS, p, T);
                    CHECK(pt->phase() == iphase_twophase);
                    if (pt->phase() == iphase_twophase) {
                        CHECK(std::abs(pt->Q() - Q) < 1e-3);
                    }
                    ++n;
                }
            }
            CHECK(n >= 40);  // the grid was not empty
        }
    }
}

#endif
