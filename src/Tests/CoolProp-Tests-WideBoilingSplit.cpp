// PT/QT consistency for WIDE-BOILING mixtures whose incipient phase is small and nearly pure
// (GitHub: flash-trial-compositions).
//
// The mixture PT flash seeds its stability test / phase split from ideal Wilson K-factors and two
// trial phases (z*K vapour-like, z/K liquid-like), refined by successive substitution.  A QT flash
// instead imposes the quality and lands directly on the phase boundary.  They must agree: if a QT
// flash produces a genuine two-phase state at (T, P), a PT flash at that same (T, P) must also find
// two phases.  Each test proves the split with a QT flash, then asserts the PT flash agrees.
//
// Two physical archetypes:
//   * Water/CO2 -> incipient near-pure WATER LIQUID condensing out of a CO2-rich gas.  Before the
//        near-pure recovery (guess_split_from_wilson, Strategy 2) this failed on SRK / PR / HEOS:
//        water's Wilson K is so extreme that the z/K liquid trial collapses back to the trivial root
//        before reaching the near-pure-water composition, so the PT flash reported single-phase.
//   * H2/CO2  -> incipient near-pure H2 VAPOUR evaporating from a CO2-rich liquid.  CoolProp's SS
//        refinement already reached this split for representative dilute-H2 feeds; kept as a GUARD so
//        the near-pure recovery does not regress the light-gas-vapour case.  (H2/CO2 only misses in a
//        degenerate ~100 ppm-H2 corner where the "vapour" is essentially pure-CO2 saturation -- not a
//        representative wide-boiling split.)
//
// The tests assert the SAME correct behaviour everywhere: PT must match QT.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/DataStructures.h"
#    include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include "../Backends/Helmholtz/VLERoutines.h"

#    include <algorithm>
#    include <cmath>
#    include <memory>
#    include <string>
#    include <vector>

using namespace CoolProp;

namespace {

struct SplitProbe
{
    bool qt_two_phase;         // did the reference QT flash land on a genuine two-phase state?
    double P;                  // boundary pressure from the QT flash [Pa]
    bool pt_two_phase;         // did the PT flash at (T, P) also find two phases?
    double pt_Q;               // vapour quality reported by the PT flash (< 0 or > 1 => single-phase)
    std::vector<double> qt_x;  // QT reference liquid composition (all components)
    std::vector<double> qt_y;  // QT reference vapour composition (all components)
    std::vector<double> pt_x;  // PT-flash liquid composition (filled only when two-phase)
    std::vector<double> pt_y;  // PT-flash vapour composition (filled only when two-phase)
    double pt_mass_balance;    // max_i |(1-Q) x_i + Q y_i - z_i| of the PT split; NaN if any term is NaN
};

// Reference the split with a QT flash (imposes Q -> boundary), then test PT at that same (T, P).
SplitProbe probe_split(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double Q) {
    SplitProbe r{};
    std::shared_ptr<AbstractState> ref(AbstractState::factory(backend, fluids));
    ref->set_mole_fractions(z);
    ref->update(QT_INPUTS, Q, T);
    r.P = ref->p();
    r.qt_two_phase = (ref->phase() == iphase_twophase);
    for (auto v : ref->mole_fractions_liquid())
        r.qt_x.push_back(static_cast<double>(v));
    for (auto v : ref->mole_fractions_vapor())
        r.qt_y.push_back(static_cast<double>(v));

    std::shared_ptr<AbstractState> pt(AbstractState::factory(backend, fluids));
    pt->set_mole_fractions(z);
    try {
        pt->update(PT_INPUTS, r.P, T);
        r.pt_two_phase = (pt->phase() == iphase_twophase);
        r.pt_Q = r.pt_two_phase ? pt->Q() : -1.0;
        if (r.pt_two_phase) {
            for (auto v : pt->mole_fractions_liquid())
                r.pt_x.push_back(static_cast<double>(v));
            for (auto v : pt->mole_fractions_vapor())
                r.pt_y.push_back(static_cast<double>(v));
            // NaN-propagating max: std::max(0, NaN) would silently return 0 and pass the check below.
            r.pt_mass_balance = (r.pt_x.size() == z.size() && r.pt_y.size() == z.size()) ? 0.0 : NAN;
            for (std::size_t i = 0; i < z.size() && r.pt_x.size() == z.size() && r.pt_y.size() == z.size(); ++i) {
                const double err = std::abs((1.0 - r.pt_Q) * r.pt_x[i] + r.pt_Q * r.pt_y[i] - z[i]);
                if (std::isnan(err) || std::isnan(r.pt_mass_balance)) {
                    r.pt_mass_balance = NAN;  // once NaN, stays NaN: a later finite term must not overwrite it
                } else if (err > r.pt_mass_balance) {
                    r.pt_mass_balance = err;
                }
            }
        }
    } catch (...) {
        r.pt_two_phase = false;  // a throw from the flash is also a "missed split" for our purposes
        r.pt_Q = -1.0;
    }
    return r;
}

// Every component of a PT-flash phase must match the QT reference phase RELATIVELY, so the trace
// (minority) component -- n-decane in a near-pure H2 vapour, dissolved H2 in a CO2-rich liquid -- is
// held to the same 1e-3 relative accuracy as the major one; an absolute or major-component check would
// let a wrong incipient phase through.  The tiny margin only covers exact-zero references.
void check_phase_matches(const char* phase, const std::vector<double>& pt, const std::vector<double>& qt) {
    REQUIRE(pt.size() == qt.size());
    for (std::size_t i = 0; i < qt.size(); ++i) {
        CAPTURE(phase, i, pt[i], qt[i]);
        CHECK(std::isfinite(pt[i]));
        CHECK(pt[i] == Catch::Approx(qt[i]).epsilon(1e-3).margin(1e-12));
    }
}

void check_pt_matches_qt(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double Q) {
    SplitProbe r = probe_split(backend, fluids, z, T, Q);
    CAPTURE(backend, fluids, z, T, Q, r.P, r.pt_Q, r.qt_x, r.qt_y, r.pt_x, r.pt_y, r.pt_mass_balance);
    REQUIRE(r.qt_two_phase);  // sanity: the QT reference really is two-phase here
    CHECK(r.pt_two_phase);    // PT flash must agree; near-pure incipient phase must not be missed
    // The PT flash must land on the SAME physical state as the QT reference, not merely on "a" split.
    // For a binary at fixed (T, P) the phase compositions are fixed (phase rule), so they must match
    // the QT reference's; and the published (x, y, Q) must reconstruct the feed -- this guards the
    // material-balance fix directly.  The vapour fraction itself is only checked loosely: the mixture
    // QT flash's reported Q is not always consistent with its own x, y at extreme Q (e.g. H2/n-decane,
    // QT Q = 1e-4 where its x, y imply 2.1e-4), so a tight Q comparison would test the reference.
    if (r.pt_two_phase) {
        check_phase_matches("liquid", r.pt_x, r.qt_x);
        check_phase_matches("vapour", r.pt_y, r.qt_y);
        CHECK(r.pt_mass_balance <= 1e-6);  // false for NaN
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

// Same guard, with the state placed relative to the backend's OWN phase boundary so one case holds on
// every equation of state: the boundary pressure comes from a QT flash at Q_boundary (1 = dew, 0 =
// bubble), scaled by `factor` -- below a dew pressure (factor < 1) the light-rich gas is single-phase
// vapour, above a bubble pressure (factor > 1) the heavy-rich feed is compressed liquid.
void check_single_phase_off_boundary(const std::string& backend, const std::string& fluids, const std::vector<double>& z, double T, double Q_boundary,
                                     double factor) {
    std::shared_ptr<AbstractState> ref(AbstractState::factory(backend, fluids));
    ref->set_mole_fractions(z);
    ref->update(QT_INPUTS, Q_boundary, T);
    // The mixture QT flash always reports two-phase, so check what can actually go wrong: a usable
    // boundary pressure.  (A throwing QT flash fails the case as an error, not silently.)
    REQUIRE(std::isfinite(ref->p()));
    REQUIRE(ref->p() > 0);
    check_single_phase_at(backend, fluids, z, T, factor * ref->p());
}

struct SplitCase
{
    const char* backend;
    const char* fluids;
    std::vector<double> z;
    double T;
    double Q;
};

struct OffBoundaryCase
{
    const char* backend;
    const char* fluids;
    std::vector<double> z;
    double T;
    double Q_boundary;
    double factor;
};

}  // namespace

// --- Water/CO2: incipient near-pure WATER LIQUID just inside the water dew point ---
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (SRK)", "[flash][mixture]") {
    check_pt_matches_qt("SRK", "CarbonDioxide&Water", {0.995, 0.005}, 298.0, 0.9999);
}
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (PR)", "[flash][mixture]") {
    check_pt_matches_qt("PR", "CarbonDioxide&Water", {0.995, 0.005}, 298.0, 0.9999);
}
TEST_CASE("Wide-boiling split: PT flash finds near-pure water liquid in Water/CO2 (HEOS)", "[flash][mixture]") {
    check_pt_matches_qt("HEOS", "CarbonDioxide&Water", {0.98, 0.02}, 298.0, 0.9999);
}

// --- H2/CO2: incipient near-pure H2 VAPOUR just inside the bubble point (guard) ---
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

// --- Broader wide-boiling coverage.  Each case was chosen from a sweep that logged which path the PT
// flash took, so the groups below exercise what their names say rather than passing incidentally. ---

// Near-pure HEAVY liquid out of a light-rich gas, reached through the near-pure recovery (Strategy 2,
// heavy seed first: the feed is vapour-like on the Wilson Rachford-Rice side).  Light gases other than
// CO2, and a heavy other than water.  On HEOS the n-decane cases are caught by the stability test
// instead and are kept as guards.
TEST_CASE("Wide-boiling split: near-pure heavy liquid out of a light gas", "[flash][mixture]") {
    const std::vector<SplitCase> cases = {
      {"HEOS", "Nitrogen&Water", {0.99, 0.01}, 320.0, 0.9999},        {"SRK", "Nitrogen&Water", {0.99, 0.01}, 320.0, 0.9999},
      {"PR", "Nitrogen&Water", {0.99, 0.01}, 320.0, 0.9999},          {"HEOS", "Methane&Water", {0.99, 0.01}, 320.0, 0.9999},
      {"SRK", "Methane&Water", {0.99, 0.01}, 320.0, 0.9999},          {"PR", "Methane&Water", {0.99, 0.01}, 320.0, 0.9999},
      {"HEOS", "Methane&n-Decane", {0.99, 0.01}, 300.0, 0.9999},      {"SRK", "Methane&n-Decane", {0.99, 0.01}, 300.0, 0.9999},
      {"PR", "Methane&n-Decane", {0.99, 0.01}, 300.0, 0.9999},        {"HEOS", "CarbonDioxide&n-Decane", {0.99, 0.01}, 300.0, 0.9999},
      {"SRK", "CarbonDioxide&n-Decane", {0.99, 0.01}, 300.0, 0.9999}, {"PR", "CarbonDioxide&n-Decane", {0.99, 0.01}, 300.0, 0.9999},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.backend << " " << c.fluids << " z0=" << c.z[0] << " T=" << c.T << " Q=" << c.Q) {
            check_pt_matches_qt(c.backend, c.fluids, c.z, c.T, c.Q);
        }
    }
}

// Near-pure LIGHT vapour out of a heavy-rich liquid.  HEOS H2/n-decane is the case that reaches the
// near-pure recovery with the light (vapour) seed first, and needs its ideal-gas density fallback: the
// global density search throws for near-pure H2 at 300 K (no van der Waals loop at ~9 Tc).  The other
// rows are caught by the stability test and guard the bubble side against regressions.
TEST_CASE("Wide-boiling split: near-pure light vapour out of a heavy liquid", "[flash][mixture]") {
    const std::vector<SplitCase> cases = {
      {"HEOS", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0001},
      {"HEOS", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.01},
      {"SRK", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0001},
      {"PR", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0001},
      {"HEOS", "Methane&n-Decane", {0.3, 0.7}, 350.0, 0.0001},
      {"SRK", "Methane&n-Decane", {0.3, 0.7}, 350.0, 0.0001},
      {"PR", "Methane&n-Decane", {0.3, 0.7}, 350.0, 0.0001},
      {"HEOS", "Nitrogen&n-Decane", {0.05, 0.95}, 300.0, 0.0001},
      {"HEOS", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0001},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.backend << " " << c.fluids << " z0=" << c.z[0] << " T=" << c.T << " Q=" << c.Q) {
            check_pt_matches_qt(c.backend, c.fluids, c.z, c.T, c.Q);
        }
    }
}

// Narrow-boiling controls: Kmax/Kmin is below the near-pure recovery's 1e3 gate, so the ordinary
// stability-test / Strategy 1 path must keep handling these unchanged.
TEST_CASE("Wide-boiling split: narrow-boiling controls are unaffected", "[flash][mixture]") {
    const std::vector<SplitCase> cases = {
      {"HEOS", "Methane&Ethane", {0.5, 0.5}, 200.0, 0.5},
      {"SRK", "Methane&Ethane", {0.5, 0.5}, 200.0, 0.5},
      {"PR", "Methane&Ethane", {0.5, 0.5}, 200.0, 0.5},
      {"HEOS", "CarbonDioxide&Nitrogen", {0.5, 0.5}, 230.0, 0.5},
      {"SRK", "CarbonDioxide&Nitrogen", {0.5, 0.5}, 230.0, 0.9999},
      {"PR", "CarbonDioxide&Nitrogen", {0.5, 0.5}, 230.0, 0.9999},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.backend << " " << c.fluids << " z0=" << c.z[0] << " T=" << c.T << " Q=" << c.Q) {
            check_pt_matches_qt(c.backend, c.fluids, c.z, c.T, c.Q);
        }
    }
}

// False-positive guards on both sides of the envelope.  The near-pure recovery now tries a second seed
// on the side the ideal estimate calls wrong, so a genuinely single-phase state must still come back
// single-phase: a light-rich gas below its water dew pressure, and a heavy-rich liquid compressed above
// its bubble pressure (where the second seed is a near-pure heavy "liquid" against a liquid feed).
TEST_CASE("Wide-boiling split: single-phase states off the envelope stay single-phase", "[flash][mixture]") {
    const std::vector<OffBoundaryCase> cases = {
      // below the dew pressure: single-phase vapour
      {"HEOS", "CarbonDioxide&Water", {0.98, 0.02}, 298.0, 1.0, 0.8},
      {"HEOS", "Nitrogen&Water", {0.99, 0.01}, 320.0, 1.0, 0.8},
      {"SRK", "Nitrogen&Water", {0.99, 0.01}, 320.0, 1.0, 0.8},
      {"PR", "Nitrogen&Water", {0.99, 0.01}, 320.0, 1.0, 0.8},
      {"HEOS", "Methane&Water", {0.99, 0.01}, 320.0, 1.0, 0.8},
      // above the bubble pressure: compressed liquid
      {"HEOS", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0, 1.2},
      {"SRK", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0, 1.2},
      {"PR", "CarbonDioxide&Hydrogen", {0.999, 0.001}, 250.0, 0.0, 1.2},
      {"HEOS", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0, 1.2},
      {"SRK", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0, 1.2},
      {"PR", "Hydrogen&n-Decane", {0.02, 0.98}, 300.0, 0.0, 1.2},
      {"HEOS", "Methane&n-Decane", {0.3, 0.7}, 350.0, 0.0, 1.2},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.backend << " " << c.fluids << " z0=" << c.z[0] << " T=" << c.T << " x" << c.factor << " of Q=" << c.Q_boundary
                                  << " pressure") {
            check_single_phase_off_boundary(c.backend, c.fluids, c.z, c.T, c.Q_boundary, c.factor);
        }
    }
}

// --- Sub-backend phase hygiene and state independence -------------------------------------------------
// The SatL / SatV sub-backends are created with a liquid / gas phase imposed, and that imposed phase acts as a
// branch selector for every density solve on them.  A flash that leaves them with no (or the wrong) imposed
// phase -- e.g. a bare specify_phase()/unspecify_phase() pair, or one interrupted by a throw -- silently
// changes every LATER flash on the same backend.  These tests pin both halves of that contract.

namespace {

// Natural-gas-like 5-component mixture whose benchmark sequence exposed the leak (GERG-2008).
const char* const N2MIX_FLUIDS = "Nitrogen&Methane&Ethane&n-Butane&n-Pentane";
const std::vector<double> N2MIX_Z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
const char* const AMARILLO_FLUIDS = "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane";
const std::vector<double> AMARILLO_Z = {0.906724, 0.031284, 0.004676, 0.045279, 0.00828, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393};

struct TP
{
    double T, p;
};
// Exact states from the 2000-state benchmark sequence (seed 42).  N2MIX_THROWER is the flash whose
// feed-density fallback threw and, before the fix, left SatL with a GAS phase imposed for good; the three
// N2MIX_SPLITS are two-phase states that were then published single-phase.
const TP N2MIX_THROWER = {193.6161018699774, 3264240.6120160879};
const TP N2MIX_SPLITS[] = {
  {105.53449217556991, 12240237.163369412}, {112.56308294554636, 12202073.496515781}, {115.91656425615635, 5090636.8002782827}};
// Amarillo at ~180 K, 5-10 MPa: a near-pure n-hexane liquid trial reports instability (possibly a genuine
// methane/n-hexane liquid-liquid split) that the split solver cannot follow ("lost a phase density solve").
const TP AMARILLO_LLE[] = {{179.71102122230218, 9845253.2507431675},
                           {186.17175133436808, 5170463.8953728043},
                           {182.40775807259934, 6712188.8902760586},
                           {181.83245805823532, 6415294.9663369628}};

HelmholtzEOSMixtureBackend& as_heos(AbstractState& AS) {
    auto* H = dynamic_cast<HelmholtzEOSMixtureBackend*>(&AS);
    REQUIRE(H != nullptr);
    return *H;
}

void check_sub_backend_phases(HelmholtzEOSMixtureBackend& H) {
    REQUIRE(H.SatL);
    REQUIRE(H.SatV);
    CHECK(H.SatL->imposed_phase() == iphase_liquid);
    CHECK(H.SatV->imposed_phase() == iphase_gas);
}

void flash_ignoring_errors(AbstractState& AS, double T, double p) {
    try {
        AS.update(PT_INPUTS, p, T);
    } catch (const CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
        // a throwing flash is allowed here; what is under test is the state it leaves behind
    }
}

}  // namespace

TEST_CASE("Wide-boiling split: PT flash leaves the SatL / SatV imposed phases unchanged", "[flash][mixture]") {
    SECTION("HEOS CO2/water, near-pure water liquid out of the gas") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&Water"));
        AS->set_mole_fractions({0.98, 0.02});
        check_sub_backend_phases(as_heos(*AS));
        for (const double p : {1.0e5, 1.6e5, 5.0e5, 6.0e6}) {
            flash_ignoring_errors(*AS, 298.0, p);
            CAPTURE(p);
            check_sub_backend_phases(as_heos(*AS));
        }
    }
    SECTION("HEOS H2/n-decane, near-pure H2 vapour out of the liquid") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Hydrogen&n-Decane"));
        AS->set_mole_fractions({0.02, 0.98});
        for (const double p : {1.0e6, 3.5e6, 7.3e6, 9.0e6}) {
            flash_ignoring_errors(*AS, 300.0, p);
            CAPTURE(p);
            check_sub_backend_phases(as_heos(*AS));
        }
    }
    SECTION("GERG-2008 N2/C1/C2/nC4/nC5, including a flash that throws inside the stability test") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", N2MIX_FLUIDS));
        AS->set_mole_fractions(N2MIX_Z);
        // On the reference platform (Linux, GCC) this flash throws inside the stability test's feed-density
        // fallback, the path that used to leak; elsewhere it may not throw, in which case this section only
        // checks the normal path (the ScopedImposedPhase test above covers the unwinding contract directly).
        flash_ignoring_errors(*AS, N2MIX_THROWER.T, N2MIX_THROWER.p);
        check_sub_backend_phases(as_heos(*AS));
        for (const auto& s : N2MIX_SPLITS) {
            flash_ignoring_errors(*AS, s.T, s.p);
            CAPTURE(s.T, s.p);
            check_sub_backend_phases(as_heos(*AS));
        }
    }
    SECTION("GERG-2008 Amarillo at the liquid-liquid-like states") {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", AMARILLO_FLUIDS));
        AS->set_mole_fractions(AMARILLO_Z);
        for (const auto& s : AMARILLO_LLE) {
            flash_ignoring_errors(*AS, s.T, s.p);
            CAPTURE(s.T, s.p);
            check_sub_backend_phases(as_heos(*AS));
        }
    }
}

TEST_CASE("Wide-boiling split: ScopedImposedPhase restores the imposed phase, also when unwinding", "[flash][mixture]") {
    // The exception-safety contract the stability test and the flash rely on, tested directly (the flash-level
    // tests below depend on a particular state actually throwing, which can vary across platforms).
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "Methane&Ethane"));
    auto& H = as_heos(*AS);
    REQUIRE(H.SatL);
    HelmholtzEOSMixtureBackend& L = *H.SatL;
    REQUIRE(L.imposed_phase() == iphase_liquid);
    SECTION("normal exit") {
        {
            const ScopedImposedPhase gas(L, iphase_gas);
            CHECK(L.imposed_phase() == iphase_gas);
        }
        CHECK(L.imposed_phase() == iphase_liquid);
    }
    SECTION("exit by exception") {
        try {
            const ScopedImposedPhase gas(L, iphase_gas);
            throw CoolProp::ValueError("simulated density-solve failure");
        } catch (const CoolProp::CoolPropBaseError&) {  // NOLINT(bugprone-empty-catch)
        }
        CHECK(L.imposed_phase() == iphase_liquid);
    }
    SECTION("lifting the imposed phase, nested scopes, inactive scope") {
        {
            const ScopedImposedPhase none(L, iphase_not_imposed);
            CHECK(L.imposed_phase() == iphase_not_imposed);
            {
                const ScopedImposedPhase gas(L, iphase_gas);
                CHECK(L.imposed_phase() == iphase_gas);
            }
            CHECK(L.imposed_phase() == iphase_not_imposed);
            const ScopedImposedPhase inactive(L, iphase_gas, false);
            CHECK(L.imposed_phase() == iphase_not_imposed);
        }
        CHECK(L.imposed_phase() == iphase_liquid);
    }
}

TEST_CASE("Wide-boiling split: PT flash result does not depend on the previous flash", "[flash][mixture]") {
    // Fresh backend per state vs one backend that first ran the flash that used to leave SatL with a gas
    // phase imposed: the three two-phase states must come out identically.
    std::shared_ptr<AbstractState> seq(AbstractState::factory("GERG2008", N2MIX_FLUIDS));
    seq->set_mole_fractions(N2MIX_Z);
    flash_ignoring_errors(*seq, N2MIX_THROWER.T, N2MIX_THROWER.p);
    for (const auto& s : N2MIX_SPLITS) {
        std::shared_ptr<AbstractState> fresh(AbstractState::factory("GERG2008", N2MIX_FLUIDS));
        fresh->set_mole_fractions(N2MIX_Z);
        CAPTURE(s.T, s.p);
        REQUIRE_NOTHROW(fresh->update(PT_INPUTS, s.p, s.T));
        REQUIRE_NOTHROW(seq->update(PT_INPUTS, s.p, s.T));
        CHECK(fresh->phase() == iphase_twophase);
        CHECK(seq->phase() == fresh->phase());
        CHECK(seq->Q() == Catch::Approx(fresh->Q()).epsilon(1e-9));
        CHECK(seq->rhomolar() == Catch::Approx(fresh->rhomolar()).epsilon(1e-9));
    }
}

TEST_CASE("Wide-boiling split: an instability the split solver cannot follow does not make the flash throw", "[flash][mixture]") {
    // Before the near-pure stability trials these states were answered single-phase; a near-pure trial now
    // reports an instability there, and the flash must fall back to an answer rather than throw.  Whether the
    // state is truly single-phase (or a methane/n-hexane liquid-liquid split) is NOT asserted here.
    for (const auto& s : AMARILLO_LLE) {
        std::shared_ptr<AbstractState> AS(AbstractState::factory("GERG2008", AMARILLO_FLUIDS));
        AS->set_mole_fractions(AMARILLO_Z);
        CAPTURE(s.T, s.p);
        REQUIRE_NOTHROW(AS->update(PT_INPUTS, s.p, s.T));
        CHECK(std::isfinite(AS->rhomolar()));
        CHECK(AS->rhomolar() > 0);
    }
}

TEST_CASE("Wide-boiling split: the stability test's density solve never returns an unstable-branch root", "[flash][mixture]") {
    // GitHub #3448: for a water-rich liquid, solver_rho_Tp_global (on a backend with a phase imposed, as on
    // SatL / SatV) returns the mechanically unstable middle root (~47.8 kmol/m3, dp/drho < 0) instead of the
    // liquid root (~55.3 kmol/m3).  The guarded solve used by the stability test must return the latter, and
    // must leave the backend's imposed phase as it found it.
    std::shared_ptr<AbstractState> AS(AbstractState::factory("HEOS", "CarbonDioxide&Water"));
    auto& H = as_heos(*AS);
    H.set_mole_fractions({1e-6, 1 - 1e-6});
    const double T = 298.0, p = 1.5964e5;
    double rho_liquid = -1;
    {
        const ScopedImposedPhase liquid(H, iphase_liquid);
        rho_liquid = H.solver_rho_Tp(T, p);
    }
    REQUIRE(rho_liquid > 5.0e4);  // sanity: the phase-specified solve finds liquid water
    for (const phases imposed : {iphase_liquid, iphase_gas}) {
        CAPTURE(imposed);
        const ScopedImposedPhase scope(H, imposed);
        bool replaced = false;
        const double rho = SaturationSolvers::solve_rho_Tp_global_stable(H, T, p, &replaced);
        H.update_DmolarT_direct(rho, T);
        // While #3448 is open the plain global solver returns the unstable root here, so the guard replaces
        // it.  Not asserted, so this test keeps passing once #3448 is fixed; the INFO is reported only when a
        // check below fails.
        INFO("guard replaced the global root: " << replaced);
        CHECK(H.first_partial_deriv(iP, iDmolar, iT) > 0);
        CHECK(rho == Catch::Approx(rho_liquid).epsilon(1e-6));
        CHECK(H.imposed_phase() == imposed);
    }
}

#endif
