// PT/QT consistency for WIDE-BOILING mixtures whose incipient phase is small and nearly pure
// (GitHub: flash-trial-compositions).
//
// The mixture PT flash seeds its stability test / phase split from ideal Wilson K-factors and two
// trial phases (z*K vapour-like, z/K liquid-like), refined by successive substitution.  A QT flash
// instead imposes the quality and lands directly on the phase boundary.  They must agree: if a QT
// flash produces a genuine two-phase state at (T, P), a PT flash at that same (T, P) must also find
// two phases.  Each test proves the split with a QT flash, then asserts the PT flash agrees.
//
// Two physical archetypes (the two extra trial compositions a general fix should add):
//   * Water/CO2 -> incipient near-pure WATER LIQUID condensing out of a CO2-rich gas.
//        CURRENTLY FAILS on SRK / PR / HEOS: water's Wilson K is so extreme that the z/K liquid
//        trial collapses back to the trivial root before reaching the near-pure-water composition,
//        so the PT flash reports single-phase.  This is the primary bug the near-pure-LIQUID trial
//        (+ re-seed of the split K from the detected incipient phase) is meant to fix.
//   * H2/CO2  -> incipient near-pure H2 VAPOUR evaporating from a CO2-rich liquid.
//        CURRENTLY PASSES: CoolProp's SS refinement already reaches this split for representative
//        dilute-H2 feeds.  Kept as a GUARD so the general near-pure-trial fix does not regress the
//        light-gas-vapour case.  (H2/CO2 only misses in a degenerate ~100 ppm-H2 corner where the
//        "vapour" is essentially pure-CO2 saturation -- not a representative wide-boiling split.)
//
// The tests assert the SAME correct behaviour everywhere (PT must match QT); water/CO2 is red until
// the trial compositions are generalised, H2/CO2 is green now and must stay green.

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

struct SplitProbe
{
    bool qt_two_phase;   // did the reference QT flash land on a genuine two-phase state?
    double P;            // boundary pressure from the QT flash [Pa]
    bool pt_two_phase;   // did the PT flash at (T, P) also find two phases?
    double pt_Q;         // vapour quality reported by the PT flash (< 0 or > 1 => single-phase)
    double xL_light;     // incipient-liquid mole fraction of the light (first) component
    double yV_light;     // incipient-vapour mole fraction of the light (first) component
};

// Reference the split with a QT flash (imposes Q -> boundary), then test PT at that same (T, P).
SplitProbe probe_split(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double Q) {
    SplitProbe r{};
    std::shared_ptr<AbstractState> ref(AbstractState::factory(backend, fluids));
    ref->set_mole_fractions(z);
    ref->update(QT_INPUTS, Q, T);
    r.P = ref->p();
    r.qt_two_phase = (ref->phase() == iphase_twophase);
    r.xL_light = ref->mole_fractions_liquid()[0];
    r.yV_light = ref->mole_fractions_vapor()[0];

    std::shared_ptr<AbstractState> pt(AbstractState::factory(backend, fluids));
    pt->set_mole_fractions(z);
    try {
        pt->update(PT_INPUTS, r.P, T);
        r.pt_two_phase = (pt->phase() == iphase_twophase);
        r.pt_Q = r.pt_two_phase ? pt->Q() : -1.0;
    } catch (...) {
        r.pt_two_phase = false;  // a throw from the flash is also a "missed split" for our purposes
        r.pt_Q = -1.0;
    }
    return r;
}

void check_pt_matches_qt(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double Q) {
    SplitProbe r = probe_split(backend, fluids, z, T, Q);
    CAPTURE(backend, fluids, T, Q, r.P, r.pt_Q, r.xL_light, r.yV_light);
    REQUIRE(r.qt_two_phase);   // sanity: the QT reference really is two-phase here
    CHECK(r.pt_two_phase);     // PT flash must agree; near-pure incipient phase must not be missed
    // The PT flash must land on the SAME physical state as the QT reference, not merely on "a" split:
    // its vapor fraction must match the imposed QT quality.  This guards the material-balance fix --
    // a published split whose (x, y, beta) did not reconstruct the feed would report the wrong Q here.
    if (r.pt_two_phase) {
        CHECK(std::abs(r.pt_Q - Q) <= 0.05);
    }
}

// False-positive guard: BELOW the water dew pressure a CO2-rich gas is genuinely single-phase.  The
// near-pure recovery must NOT publish a spurious two-phase split there -- its verify step must reject
// a pressure-inconsistent SS seed (an earlier revision accepted a "liquid" root from the negative-
// pressure spinodal region and published a bogus Q=0.5 split).  These (T, P) states are confirmed
// single-phase (water partial pressure below its saturation pressure) and each reproduced the bogus
// split before the phase-pressure-consistency check was added.
void check_single_phase_at(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double P) {
    std::shared_ptr<AbstractState> pt(AbstractState::factory(backend, fluids));
    pt->set_mole_fractions(z);
    pt->update(PT_INPUTS, P, T);
    CAPTURE(backend, fluids, T, P, pt->Q(), pt->phase());
    CHECK(pt->phase() != iphase_twophase);
}

}  // namespace

// --- Water/CO2: incipient near-pure WATER LIQUID just inside the water dew point (CURRENTLY FAILS) ---
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (SRK)", "[flash][mixture]") {
    check_pt_matches_qt("SRK", "CarbonDioxide&Water", {0.995, 0.005}, 298.0, 0.9999);
}
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (PR)", "[flash][mixture]") {
    check_pt_matches_qt("PR", "CarbonDioxide&Water", {0.995, 0.005}, 298.0, 0.9999);
}
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (HEOS)", "[flash][mixture]") {
    check_pt_matches_qt("HEOS", "CarbonDioxide&Water", {0.98, 0.02}, 298.0, 0.9999);
}

// --- H2/CO2: incipient near-pure H2 VAPOUR just inside the bubble point (GUARD: currently passes) ---
TEST_CASE("Wide-boiling split: PT flash finds near-pure H2 vapour in H2/CO2 (SRK)", "[flash][mixture]") {
    check_pt_matches_qt("SRK", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0001);
}
TEST_CASE("Wide-boiling split: PT flash finds near-pure H2 vapour in H2/CO2 (PR)", "[flash][mixture]") {
    check_pt_matches_qt("PR", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0001);
}

// --- False-positive guard: Water/CO2 below the water dew must stay single-phase (no bogus split) ---
TEST_CASE("Wide-boiling split: Water/CO2 below the water dew stays single-phase (PR)", "[flash][mixture]") {
    const std::vector<double> z = {0.995, 0.005};
    check_single_phase_at("PR", "CarbonDioxide&Water", z, 295.0, 4.069e5);
    check_single_phase_at("PR", "CarbonDioxide&Water", z, 300.0, 2.326e5);
    check_single_phase_at("PR", "CarbonDioxide&Water", z, 305.0, 1.530e5);
    check_single_phase_at("PR", "CarbonDioxide&Water", z, 330.0, 6.612e4);
}
TEST_CASE("Wide-boiling split: Water/CO2 below the water dew stays single-phase (SRK)", "[flash][mixture]") {
    const std::vector<double> z = {0.995, 0.005};
    check_single_phase_at("SRK", "CarbonDioxide&Water", z, 300.0, 2.326e5);
    check_single_phase_at("SRK", "CarbonDioxide&Water", z, 305.0, 1.530e5);
    check_single_phase_at("SRK", "CarbonDioxide&Water", z, 330.0, 6.612e4);
}

#endif
