// Characterization tests for D + {H, S, U} flash inputs (HSU_D_flash).
//
// These tests pin the superancillary "happy path" ON and verify the
// density-plus-caloric flashes reproduce their inputs across the regions
// that were historically fragile; they are the acceptance harness for the
// hybrid superancillary/fallback rewrite tracked in bd CoolProp-j3n.  (The
// legacy ancillary fallback is known-broken on exactly these points -- that
// is the rewrite's motivation -- so it is not asserted here; the opt-in
// [HSU_D_bench] case exercises it under COOLPROP_DISABLE_SUPERANC_HSU_D.)
//
// Every case is mined from a real (mostly closed) GitHub issue and tagged
// with its number so a failure points straight at the report:
//
//   [2486]  Water  DmolarSmolar round-trip along the saturation boundary
//   [2157]  R1233ZD(E)/R32  D,U two-phase quality-band failures
//   [1698]  R134a  D,U two-phase "Brent reached max steps"
//   [1054]  R245fa D,S two-phase point that fails while neighbours pass
//   [2154]  Hydrogen  D,U near the critical density -> "p is not valid"
//   [2173]  H2/He  D,U near-critical-density sweep (PR #2173 fix)
//   [1965]  Nitrogen  DmassUmass -> "matrix is singular"
//   [2022]  Air  Hmass,Smass (related HS pair, not D+; characterised too)
//   [2685]  REFPROP::MM  DmassUmass -> smass == inf
//   [2426]  REFPROP::Hydrogen  DmassUmass -> ungraceful failure / wrong h
//
// The checks are deliberately self-consistency based: after a D+X flash the
// resulting state must reproduce the input density AND the input caloric
// property (and, where a saturation/PT reference exists, the temperature).
// Mismatches are recorded with FAIL_CHECK rather than aborting, so a single
// run characterises the entire failure surface instead of stopping at the
// first bad point.
//
// Built into CatchTestRunner when COOLPROP_CATCH_MODULE=ON.  Run the whole
// suite with `CatchTestRunner "[HSU_D]"`.

#include "AbstractState.h"
#include "CoolProp.h"
#include "CPstrings.h"
#include "DataStructures.h"
#include "Configuration.h"
#include "../Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "superancillary/superancillary.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <cstdlib>
#    include <limits>
#    include <memory>
#    include <string>
#    include <vector>

namespace {

// RAII toggle of the ENABLE_SUPERANCILLARIES config flag, restored on scope
// exit.  These tests pin it ON (SuperancGuard guard(true)) so the
// superancillary "happy path" is exercised hermetically regardless of the
// ambient config.
struct SuperancGuard
{
    bool saved;
    explicit SuperancGuard(bool on) : saved(CoolProp::get_config_bool(ENABLE_SUPERANCILLARIES)) {
        CoolProp::set_config_bool(ENABLE_SUPERANCILLARIES, on);
    }
    ~SuperancGuard() {
        CoolProp::set_config_bool(ENABLE_SUPERANCILLARIES, saved);
    }
    SuperancGuard(const SuperancGuard&) = delete;
    SuperancGuard& operator=(const SuperancGuard&) = delete;
    SuperancGuard(SuperancGuard&&) = delete;
    SuperancGuard& operator=(SuperancGuard&&) = delete;
};

std::vector<double> linspace(double a, double b, std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) {
        v[i] = (n == 1) ? a : a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
    }
    return v;
}

// The mass-based caloric accessor matching a D+X input pair.
double other_mass(CoolProp::AbstractState& AS, CoolProp::input_pairs pair) {
    switch (pair) {
        case CoolProp::DmassHmass_INPUTS:
            return AS.hmass();
        case CoolProp::DmassSmass_INPUTS:
            return AS.smass();
        case CoolProp::DmassUmass_INPUTS:
            return AS.umass();
        default:
            return _HUGE;
    }
}

const char* pair_name(CoolProp::input_pairs pair) {
    switch (pair) {
        case CoolProp::DmassHmass_INPUTS:
            return "DmassHmass";
        case CoolProp::DmassSmass_INPUTS:
            return "DmassSmass";
        case CoolProp::DmassUmass_INPUTS:
            return "DmassUmass";
        default:
            return "?";
    }
}

// Re-flash a fresh state via (D, other) and verify it reproduces the inputs.
// T_expect > 0 adds a temperature-recovery check.  Records failures (throw or
// mismatch) without aborting so a sweep maps the whole surface.
void check_DX(CoolProp::AbstractState& rt, CoolProp::input_pairs pair, double d_mass, double other_in, double T_expect) {
    INFO("pair=" << pair_name(pair) << " d_mass=" << d_mass << " other_in=" << other_in << " T_expect=" << T_expect);
    try {
        rt.update(pair, d_mass, other_in);
    } catch (const std::exception& e) {
        FAIL_CHECK("update threw: " << e.what());
        return;
    } catch (...) {
        FAIL_CHECK("update threw non-std exception");
        return;
    }
    const double d_out = rt.rhomass();
    const double other_out = other_mass(rt, pair);
    CHECK(std::isfinite(d_out));
    CHECK(std::isfinite(other_out));
    CHECK(d_out == Catch::Approx(d_mass).epsilon(1e-5));
    // Caloric properties can pass through zero (reference-state dependent),
    // so guard the relative compare with a margin.
    CHECK(other_out == Catch::Approx(other_in).epsilon(1e-5).margin(1e-3 * std::abs(other_in) + 1e-6));
    if (T_expect > 0) {
        // Input recovery (d_out, other_out above) is the PRIMARY correctness
        // check -- it is exactly what the flash solves for.  T is a derived,
        // secondary guard: loose enough to tolerate the critical-region
        // ill-conditioning (dT/d(rho,X) blows up near T_crit, so a converged
        // (rho, X) root can sit ~1e-3 off in T) while still catching a flash
        // that lands on the wrong phase/branch (those miss T by whole percent).
        CHECK(rt.T() == Catch::Approx(T_expect).epsilon(1e-3));
    }
}

}  // namespace

// ---------------------------------------------------------------------------
// #2486 : round-trip on the saturation boundary (Q = 0 and Q = 1).
// This is the hardest region for HSU_D: the input density sits exactly where
// the two-phase and single-phase branches meet.  Establish the reference with
// QT_INPUTS, then re-flash through each of D+{H,S,U}.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: saturation-boundary round-trip (#2486)", "[HSU_D][HSU_D_satbound][2486]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Water", "n-Propane", "R134a", "CarbonDioxide", "Nitrogen");
    const CoolProp::input_pairs pair = GENERATE(CoolProp::DmassHmass_INPUTS, CoolProp::DmassSmass_INPUTS, CoolProp::DmassUmass_INPUTS);
    SuperancGuard guard(true);
    CAPTURE(fluid, pair_name(pair));

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tt = CoolProp::Props1SI(fluid, "Ttriple");
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit");

    for (double T : linspace(Tt + 0.5, Tc - 0.5, 60)) {
        for (double Q : {0.0, 1.0}) {
            CAPTURE(T, Q);
            ref->update(CoolProp::QT_INPUTS, Q, T);
            const double d = ref->rhomass();
            double other = 0;
            switch (pair) {
                case CoolProp::DmassHmass_INPUTS:
                    other = ref->hmass();
                    break;
                case CoolProp::DmassSmass_INPUTS:
                    other = ref->smass();
                    break;
                default:
                    other = ref->umass();
                    break;
            }
            check_DX(*rt, pair, d, other, T);
        }
    }
}

// ---------------------------------------------------------------------------
// #2157 / #1698 / #1054 : two-phase interior round-trips.
// Reports describe failures (Brent max-steps, invalid secant residual) for
// D,U and D,S inside the dome, worst at low quality.  Sweep T and Q.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: two-phase interior round-trip (#2157,#1698,#1054)", "[HSU_D][HSU_D_twophase][2157][1698][1054]") {
    // R1233ZD(E)=#2157, R32=#2157, R134a=#1698, R245fa=#1054.
    const std::string fluid = GENERATE(as<std::string>{}, "R1233zd(E)", "R32", "R134a", "R245fa");
    const CoolProp::input_pairs pair = GENERATE(CoolProp::DmassUmass_INPUTS, CoolProp::DmassSmass_INPUTS);
    SuperancGuard guard(true);
    CAPTURE(fluid, pair_name(pair));

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tt = CoolProp::Props1SI(fluid, "Ttriple");
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit");

    for (double T : linspace(Tt + 1.0, Tc - 1.0, 40)) {
        for (double Q : {0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9}) {
            CAPTURE(T, Q);
            ref->update(CoolProp::QT_INPUTS, Q, T);
            const double d = ref->rhomass();
            const double other = (pair == CoolProp::DmassUmass_INPUTS) ? ref->umass() : ref->smass();
            check_DX(*rt, pair, d, other, T);
        }
    }
}

// ---------------------------------------------------------------------------
// #2154 / #2173 : near the critical density.  saturation_D_pure destabilises
// when rho ~ rho_crit; PR #2173 patched the HEOS path.  This region is the
// hardest part of the surface, so it gets three dedicated sweeps:
//   (A) a dense (T, rho) box around the critical point, all three D+X pairs,
//   (B) the critical-isochore crossing that #2154/#1965 actually exercise,
//   (C) the near-critical two-phase dome edge (T within ~0.01 K of T_crit),
//       where HSU_D_flash_twophase's hard 0.01 K bracket limit bites (#2173).
// The grids are concentrated near rho_crit / T_crit, where the failures live.
// ---------------------------------------------------------------------------

// Density and temperature multipliers clustered around the critical point.
static const std::vector<double> kRhoFactors = {0.50,  0.75,  0.90, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0,
                                                1.001, 1.005, 1.01, 1.02, 1.05, 1.10, 1.25,  1.50,  2.00};
static const std::vector<double> kTFactors = {0.90, 0.95, 0.98, 0.99, 0.995, 1.0, 1.005, 1.01, 1.02, 1.05, 1.10};

// (A) Dense (T, rho) box around the critical point, all three D+X pairs.
TEST_CASE("HSU_D: near-critical (T,rho) box, all pairs (#2154,#2173)", "[HSU_D][HSU_D_critdens][2154][2173]") {
    const std::string fluid =
      GENERATE(as<std::string>{}, "Water", "CarbonDioxide", "n-Propane", "R134a", "Nitrogen", "Hydrogen", "Helium", "Methane");
    const CoolProp::input_pairs pair = GENERATE(CoolProp::DmassHmass_INPUTS, CoolProp::DmassSmass_INPUTS, CoolProp::DmassUmass_INPUTS);
    SuperancGuard guard(true);
    CAPTURE(fluid, pair_name(pair));

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit");
    const double rhoc = CoolProp::Props1SI(fluid, "rhomass_critical");

    for (double tf : kTFactors) {
        for (double rf : kRhoFactors) {
            // The exact critical point (T_c, rho_c) is a singular point of the
            // EOS where the density<->caloric map degenerates; flashing there is
            // ill-posed, so skip it (all near-critical neighbours are kept).
            if (tf == 1.0 && rf == 1.0) continue;
            const double T = tf * Tc, rho = rf * rhoc;
            CAPTURE(tf, rf, T, rho);
            // Reference state from (rho, T): always single-phase well-posed here.
            ref->update(CoolProp::DmassT_INPUTS, rho, T);
            check_DX(*rt, pair, rho, other_mass(*ref, pair), T);
        }
    }
}

// (B) Crossing the critical isochore: hold rho ~ rho_crit and sweep T (and
// hence internal energy) from subcritical to supercritical, exactly the
// tank-fill/defuel scenario in #2154 and #1965.
TEST_CASE("HSU_D: critical-isochore crossing (#2154,#1965)", "[HSU_D][HSU_D_critdens][2154][1965]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Hydrogen", "Helium", "Nitrogen", "CarbonDioxide", "Water", "Methane");
    const double rf = GENERATE(0.97, 0.99, 1.0, 1.01, 1.03);
    SuperancGuard guard(true);
    CAPTURE(fluid, rf);

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit");
    const double rhoc = CoolProp::Props1SI(fluid, "rhomass_critical");
    const double rho = rf * rhoc;

    for (double T : linspace(0.92 * Tc, 1.12 * Tc, 40)) {
        CAPTURE(T);
        ref->update(CoolProp::DmassT_INPUTS, rho, T);
        check_DX(*rt, CoolProp::DmassUmass_INPUTS, rho, ref->umass(), T);
    }
}

// (C) Two-phase dome edge as T -> T_crit from below.  Densities on both
// branches collapse toward rho_crit and HSU_D_flash_twophase's 0.01 K bracket
// limit on T can prevent convergence (#2173 noted this as a residual failure).
TEST_CASE("HSU_D: near-critical two-phase dome edge (#2173)", "[HSU_D][HSU_D_critdens][2173]") {
    const std::string fluid = GENERATE(as<std::string>{}, "Hydrogen", "Helium", "CarbonDioxide", "Water", "n-Propane");
    const CoolProp::input_pairs pair = GENERATE(CoolProp::DmassHmass_INPUTS, CoolProp::DmassSmass_INPUTS, CoolProp::DmassUmass_INPUTS);
    SuperancGuard guard(true);
    CAPTURE(fluid, pair_name(pair));

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
    const double Tc = CoolProp::Props1SI(fluid, "Tcrit");
    const double Tt = CoolProp::Props1SI(fluid, "Ttriple");

    // Approach Tc from below; the tightest point sits well within 0.01 K of Tc.
    // Skip offsets that would fall below the triple point (Helium's Tc is only
    // ~5.2 K, so the coarse 10 K offset is invalid there).
    for (double dT : {10.0, 1.0, 0.1, 0.01, 0.001}) {
        const double T = Tc - dT;
        if (T <= Tt + 0.1) continue;
        for (double Q : {0.1, 0.3, 0.5, 0.7, 0.9}) {
            CAPTURE(dT, T, Q);
            ref->update(CoolProp::QT_INPUTS, Q, T);
            check_DX(*rt, pair, ref->rhomass(), other_mass(*ref, pair), T);
        }
    }
}

// Exact point from #2154: PropsSI("T","D",31.46258141,"U",2391760.261,"Hydrogen").
TEST_CASE("HSU_D: #2154 reported hydrogen D,U point", "[HSU_D][2154]") {
    SuperancGuard guard(true);
    double T = 0;
    REQUIRE_NOTHROW(T = CoolProp::PropsSI("T", "D", 31.46258141, "U", 2391760.261, "Hydrogen"));
    CAPTURE(T);
    CHECK(std::isfinite(T));
    // Issue states the answer should be roughly -40..140 degC.
    CHECK(T > 233.0);
    CHECK(T < 414.0);
}

// ---------------------------------------------------------------------------
// #1965 : nitrogen DmassUmass -> "Zero occurred in row 1, the matrix is
// singular" (open).  The reported point sits near rho_crit on the liquid side.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: #1965 nitrogen DmassUmass singular matrix", "[HSU_D][1965]") {
    SuperancGuard guard(true);
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Nitrogen"));
    check_DX(*rt, CoolProp::DmassUmass_INPUTS, 313.305125, 154834.285193, /*T_expect=*/-1);
}

// ---------------------------------------------------------------------------
// #1054 : the specific R245fa D,S point that fails while neighbours pass.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: #1054 R245fa D,S isolated failure point", "[HSU_D][1054]") {
    SuperancGuard guard(true);
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R245fa"));
    // Neighbours that the issue reports as working - they must keep working.
    for (double d : {17.04, 17.05, 17.048396706125065}) {
        for (double s : {1.6174e3, 1.6175e3, 1617.463750828963}) {
            CAPTURE(d, s);
            check_DX(*rt, CoolProp::DmassSmass_INPUTS, d, s, /*T_expect=*/-1);
        }
    }
}

// ---------------------------------------------------------------------------
// #2022 : HEOS::Air with Hmass,Smass.  NOT a D+X pair, but the same family of
// flash-robustness failures; characterised here so the rewrite can confirm it
// did not regress (or, ideally, fixed) the reported points.  Air is a
// pseudo-pure mixture, so the superancillary happy path never runs on it --
// these points are governed entirely by the (unchanged) HmassSmass flash.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D-adjacent: #2022 HEOS::Air Hmass,Smass", "[HSU_D][2022]") {
    SuperancGuard guard(true);
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Air"));

    // Round-trip check for one (h, s) input; returns true iff the flash
    // reproduces both inputs at a physical (finite, positive) pressure.
    auto roundtrips = [&](double h, double s) -> bool {
        try {
            rt->update(CoolProp::HmassSmass_INPUTS, h, s);
            return std::isfinite(rt->p()) && rt->p() > 0 && rt->hmass() == Catch::Approx(h).epsilon(1e-4)
                   && rt->smass() == Catch::Approx(s).epsilon(1e-4);
        } catch (const std::exception&) {
            return false;
        }
    };

    // Points that round-trip today: assert they keep working (no-regression guard).
    const std::vector<std::pair<double, double>> hs_ok = {
      {429667.3064, 2735.612747}, {430867.4193, 2732.980272}, {427500.0, 2740.0}, {435000.0, 2724.0}};
    for (auto& [h, s] : hs_ok) {
        CAPTURE(h, s);
        CHECK(roundtrips(h, s));
    }

    // Known-broken: {h=435000, s=2654} is a physically valid Air state
    // (~T=320 K, p~8.3 MPa, just lower entropy / higher pressure than the
    // {435000, 2724} neighbour above) but the HmassSmass flash diverges,
    // driving p toward ~1e11 Pa (returns hmass~4.5e7 or throws "Brent ... do
    // not bracket the root").  Pre-existing on master and out of scope for
    // this D+{H,S,U} PR; characterised non-fatally so a future fix flips it
    // green.  Tracked in bd CoolProp-l34.
    {
        const double h = 435000.0, s = 2654.0;
        CAPTURE(h, s);
        if (!roundtrips(h, s)) {
            WARN("known-broken Air HmassSmass point (h=" << h << ", s=" << s << ") still diverges; see bd CoolProp-l34 / gh #2022");
        }
    }
}

// ---------------------------------------------------------------------------
// REFPROP-backed reports.  REFPROP is available locally (see CLAUDE.md), so
// these run rather than skip on Ian's machine; CI skips cleanly.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: #2685 REFPROP::MM DmassUmass returns finite entropy", "[HSU_D][2685][REFPROP]") {
    CoolProp::Skip_if_No_REFPROP();
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "MM"));
    // Reported point returned smass == inf while T/p/h were fine.
    check_DX(*rt, CoolProp::DmassUmass_INPUTS, 126.64425, 413304.53, /*T_expect=*/-1);
}

TEST_CASE("HSU_D: #2426 REFPROP::Hydrogen DmassUmass round-trip", "[HSU_D][2426][REFPROP]") {
    CoolProp::Skip_if_No_REFPROP();
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Hydrogen"));
    // Reported failure (ungraceful) plus wrong hmass (~2.94e4 vs expected ~3.89e6).
    check_DX(*rt, CoolProp::DmassUmass_INPUTS, 20.0, 2.5e6, /*T_expect=*/-1);
}

// ---------------------------------------------------------------------------
// Water's saturated-liquid density has a maximum near T ~ 277 K (the famous
// density anomaly): rho_satL(T) rises from the triple point to a maximum at
// T_anom, then falls toward rho_crit.  For densities just below that maximum
// the liquid saturation branch is non-monotonic, so get_all_intersections('D',
// rho) returns MULTIPLE liquid-side roots -- the topology the happy-path
// interval classification has to get right (it cannot assume "two-phase below
// the single root").  Round-trip two-phase states whose overall density sits in
// the anomaly band, plus genuinely compressed-liquid states on either side.
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: water two-phase around the saturated-liquid density maximum", "[HSU_D][HSU_D_wateranom]") {
    const CoolProp::input_pairs pair = GENERATE(CoolProp::DmassHmass_INPUTS, CoolProp::DmassSmass_INPUTS, CoolProp::DmassUmass_INPUTS);
    SuperancGuard guard(true);
    CAPTURE(pair_name(pair));

    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Locate the saturated-liquid density maximum (T_anom ~ 277 K for water).
    double T_anom = 277.0, rhoL_max = 0.0;
    for (double T : linspace(274.0, 295.0, 600)) {
        ref->update(CoolProp::QT_INPUTS, 0.0, T);
        if (ref->rhomass() > rhoL_max) {
            rhoL_max = ref->rhomass();
            T_anom = T;
        }
    }
    CAPTURE(T_anom, rhoL_max);

    // (1) Two-phase round-trips spanning the anomaly temperature, including very
    // liquid-rich states whose mixture density lands in the non-monotonic band.
    for (double T : linspace(274.0, 300.0, 40)) {
        for (double Q : {0.0, 0.001, 0.01, 0.1, 0.5, 0.9, 1.0}) {
            CAPTURE(T, Q);
            ref->update(CoolProp::QT_INPUTS, Q, T);
            check_DX(*rt, pair, ref->rhomass(), other_mass(*ref, pair), T);
        }
    }

    // (2) Genuinely single-phase compressed liquid (high pressure, so densities
    // above rho_satL,max) across the anomaly temperature band.
    for (double p : {20e6, 50e6, 100e6}) {
        for (double T : linspace(274.0, 320.0, 16)) {
            CAPTURE(p, T);
            ref->update(CoolProp::PT_INPUTS, p, T);
            check_DX(*rt, pair, ref->rhomass(), other_mass(*ref, pair), T);
        }
    }
}

// ---------------------------------------------------------------------------
// HSU_D_TWOPHASE_EOS_POLISH config toggle: both modes must reproduce the inputs;
// polish-off (raw superancillary) should agree with polish-on (EOS-exact) to the
// superancillary precision (~1e-8).
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: two-phase EOS-polish config toggle", "[HSU_D][HSU_D_polishcfg]") {
    SuperancGuard sa(true);
    const bool saved = CoolProp::get_config_bool(HSU_D_TWOPHASE_EOS_POLISH);
    struct Restore
    {
        bool v;
        ~Restore() {
            CoolProp::set_config_bool(HSU_D_TWOPHASE_EOS_POLISH, v);
        }
        Restore(const Restore&) = delete;
        Restore& operator=(const Restore&) = delete;
        Restore(Restore&&) = delete;
        Restore& operator=(Restore&&) = delete;
    } restore{saved};

    for (const auto& fluid : {std::string("Water"), std::string("CarbonDioxide"), std::string("n-Propane")}) {
        auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        const double Tt = CoolProp::Props1SI(fluid, "Ttriple"), Tc = CoolProp::Props1SI(fluid, "Tcrit");
        for (double T : linspace(Tt + 1.0, Tc - 1.0, 12)) {
            for (double Q : {0.05, 0.5, 0.95}) {
                CAPTURE(fluid, T, Q);
                ref->update(CoolProp::QT_INPUTS, Q, T);
                const double d = ref->rhomass(), s = ref->smass();
                CoolProp::set_config_bool(HSU_D_TWOPHASE_EOS_POLISH, true);
                rt->update(CoolProp::DmassSmass_INPUTS, d, s);
                const double T_on = rt->T();
                CoolProp::set_config_bool(HSU_D_TWOPHASE_EOS_POLISH, false);
                rt->update(CoolProp::DmassSmass_INPUTS, d, s);
                const double T_off = rt->T();
                // Both reproduce the saturation temperature; they agree to ~SA precision.
                CHECK(T_on == Catch::Approx(T).epsilon(1e-6));
                CHECK(T_off == Catch::Approx(T).epsilon(1e-5));
                CHECK(std::abs(T_on - T_off) / T < 1e-6);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Comprehensive single-phase sweep: for every pure fluid that has a
// superancillary, walk a dense (log p, linear T) grid -- the same coordinate
// layout as the consistency plots -- establish each state with PT_INPUTS
// (always single phase), and round-trip it through D+{H,S,U} with the
// superancillary "happy path" enabled.  This blankets the liquid, vapor and
// supercritical regions of the whole fluid set.
//
// Heavy (hundreds of fluids x grid x 3 pairs), so it is hidden from the default
// run via the [.] tag.  Invoke explicitly:
//     CatchTestRunner "[HSU_D_ptsweep]"
//     CatchTestRunner "[HSU_D_ptsweep]" -c Water      # single fluid section
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: dense PT sweep over all superancillary fluids", "[HSU_D_ptsweep][.]") {
    SuperancGuard guard(true);

    static const std::vector<std::string> fluids = strsplit(CoolProp::get_global_param_string("fluids_list"), ',');
    const std::string fluid = GENERATE(from_range(fluids));
    CAPTURE(fluid);

    std::shared_ptr<CoolProp::AbstractState> ref, rt;
    try {
        ref.reset(CoolProp::AbstractState::factory("HEOS", fluid));
        rt.reset(CoolProp::AbstractState::factory("HEOS", fluid));
    } catch (...) {
        return;  // un-constructible fluid -> nothing to test
    }
    auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(ref.get());
    if (be == nullptr || !be->is_pure()) return;  // mixtures/pseudo-pure: no superancillary
    std::shared_ptr<CoolProp::EquationOfState::SuperAncillary_t> sa;
    try {
        sa = be->get_superanc();
    } catch (...) {
        sa.reset();  // treat a failed lookup as "no superancillary"
    }
    if (!sa) return;  // pure fluid without a superancillary -> happy path n/a, skip

    const double Tmin = ref->Tmin(), Tmax = ref->Tmax(), pmax = ref->pmax();
    double plo;
    try {
        plo = std::max(ref->p_triple(), pmax * 1e-8);
    } catch (...) {
        plo = pmax * 1e-8;
    }
    if (!(Tmax > Tmin) || !(pmax > plo)) return;

    // Grid density is tunable via env vars (no recompile needed); defaults keep
    // the opt-in run reasonable.  e.g. HSU_PTSWEEP_NT=200 HSU_PTSWEEP_NP=200 for
    // a high-confidence sweep.
    auto env_or = [](const char* name, std::size_t def) -> std::size_t {
        const char* v = std::getenv(name);
        if (v == nullptr) return def;
        try {
            const long n = std::stol(v);
            if (n > 1) return static_cast<std::size_t>(n);
        } catch (...) {
            return def;  // malformed env value -> default
        }
        return def;
    };
    const std::size_t NT = env_or("HSU_PTSWEEP_NT", 40), NP = env_or("HSU_PTSWEEP_NP", 30);
    const CoolProp::input_pairs pairs[] = {CoolProp::DmassHmass_INPUTS, CoolProp::DmassSmass_INPUTS, CoolProp::DmassUmass_INPUTS};

    for (double T : linspace(Tmin + 0.1, Tmax, NT)) {
        for (std::size_t j = 0; j < NP; ++j) {
            // Logarithmic spacing in pressure (linear in T), as in the
            // consistency plots.
            const double p = plo * std::pow(pmax / plo, static_cast<double>(j) / static_cast<double>(NP - 1));
            try {
                ref->update(CoolProp::PT_INPUTS, p, T);  // single-phase reference
            } catch (...) {
                continue;  // (p, T) outside the EOS domain (e.g. solid) -> skip
            }
            const double d = ref->rhomass();
            if (!std::isfinite(d) || d <= 0) continue;
            for (auto pr : pairs) {
                CAPTURE(T, p);
                // T_expect = -1: the temperature check is intentionally skipped
                // here.  At fixed density h/s/u(rho, T) is NOT monotonic for some
                // fluids (e.g. compressed cryogenic hydrogen has a minimum of
                // h along the isochore), so (rho, X) can map to two valid
                // single-phase temperatures.  The flash's contract is only to
                // return a state reproducing the inputs (rho and X), which is
                // exactly what the recovery checks in check_DX assert.
                check_DX(*rt, pr, d, other_mass(*ref, pr), /*T_expect=*/-1);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Regression: the matrix-balancing helper used by the superancillary
// Chebyshev root finder (CoolProp::superancillary::detail::balance_matrix) must
// always terminate, even when a companion-matrix row/column norm is non-finite.
//
// A near-zero leading Chebyshev coefficient makes the companion matrix's last
// column (-coeffs.head(N)/(2*coeffs(N))) overflow to +/-inf. With an inf row
// norm r, the scaling loop condition (c < r/beta) is perpetually true and
// c *= beta never catches up, so balance_matrix spun forever -- the root cause
// of the devdocs consistency-plot hang (Propyne D+H flash on x86 Linux, exposed
// when #2995 routed D+{H,S,U} through get_all_intersections -> get_x_for_y ->
// balance_matrix). The hang was FP-rounding dependent, hence x86-only.
//
// This drives balance_matrix directly with an inf entry so the guard is
// exercised hermetically on every platform. Without the fix the test hangs.
// ---------------------------------------------------------------------------
TEST_CASE("balance_matrix terminates on a non-finite companion entry", "[HSU_D][balance_matrix]") {
    Eigen::MatrixXd A(3, 3), Aprime, D;
    A.setZero();
    A(1, 0) = 1.0;
    A(2, 1) = 1.0;
    // Last column as produced by a near-zero leading coefficient: overflows to inf.
    A(0, 2) = -1.0 / (2.0 * std::numeric_limits<double>::min() * std::numeric_limits<double>::min());
    REQUIRE(std::isinf(A(0, 2)));

    // Must return (terminate) rather than spin in the scaling loops.
    CoolProp::superancillary::detail::balance_matrix(A, Aprime, D);

    // The inf row/column is left at unit scaling; finite rows stay finite.
    CHECK(std::isfinite(D(1, 1)));
    CHECK(D(1, 1) == 1.0);
}

// ---------------------------------------------------------------------------
// Speed comparison: the superancillary "happy path" vs the legacy ancillary
// path on representative D+{H,S,U} flashes.  Same binary, toggled per process
// via the kill-switch, so the timings are directly comparable:
//     CatchTestRunner "[HSU_D_bench]"                                  (happy)
//     COOLPROP_DISABLE_SUPERANC_HSU_D=1 CatchTestRunner "[HSU_D_bench]" (legacy)
// Hidden by default ([.]).
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: happy-vs-legacy speed", "[HSU_D_bench][.]") {
    auto rt = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto co2 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
    auto setup = [](CoolProp::AbstractState& AS, CoolProp::input_pairs ref_pair, double a, double b, CoolProp::input_pairs out) {
        AS.update(ref_pair, a, b);
        return std::pair<double, double>{AS.rhomass(), other_mass(AS, out)};
    };

    // Two-phase (mid-dome), Water, D+S
    {
        auto [d, s] = setup(*rt, CoolProp::QT_INPUTS, 0.4, 450.0, CoolProp::DmassSmass_INPUTS);
        BENCHMARK("Water two-phase D+S") {
            rt->update(CoolProp::DmassSmass_INPUTS, d, s);
            return rt->T();
        };
    }
    // Two-phase (mid-dome), Water, D+U
    {
        auto [d, u] = setup(*rt, CoolProp::QT_INPUTS, 0.4, 450.0, CoolProp::DmassUmass_INPUTS);
        BENCHMARK("Water two-phase D+U") {
            rt->update(CoolProp::DmassUmass_INPUTS, d, u);
            return rt->T();
        };
    }
    // Single-phase compressed liquid, Water, D+H
    {
        auto [d, h] = setup(*rt, CoolProp::PT_INPUTS, 50e6, 320.0, CoolProp::DmassHmass_INPUTS);
        BENCHMARK("Water liquid D+H") {
            rt->update(CoolProp::DmassHmass_INPUTS, d, h);
            return rt->T();
        };
    }
    // Superheated vapor, Water, D+U
    {
        auto [d, u] = setup(*rt, CoolProp::PT_INPUTS, 1e5, 600.0, CoolProp::DmassUmass_INPUTS);
        BENCHMARK("Water vapor D+U") {
            rt->update(CoolProp::DmassUmass_INPUTS, d, u);
            return rt->T();
        };
    }
    // Supercritical, CO2, D+H
    {
        auto [d, h] = setup(*co2, CoolProp::PT_INPUTS, 12e6, 330.0, CoolProp::DmassHmass_INPUTS);
        BENCHMARK("CO2 supercritical D+H") {
            co2->update(CoolProp::DmassHmass_INPUTS, d, h);
            return co2->T();
        };
    }
    // Near-critical two-phase, CO2, D+S
    {
        auto [d, s] = setup(*co2, CoolProp::QT_INPUTS, 0.5, 300.0, CoolProp::DmassSmass_INPUTS);
        BENCHMARK("CO2 near-crit two-phase D+S") {
            co2->update(CoolProp::DmassSmass_INPUTS, d, s);
            return co2->T();
        };
    }
    // Near-critical DENSITY, supercritical, CO2, D+U.  rho ~ rho_crit just above
    // T_crit is where legacy saturation_D_pure is least stable (omega-halving
    // retry loop), so the happy path should pull ahead here.
    {
        const double rhoc = CoolProp::Props1SI("CarbonDioxide", "rhomass_critical");
        const double Tc = CoolProp::Props1SI("CarbonDioxide", "Tcrit");
        co2->update(CoolProp::DmassT_INPUTS, rhoc, 1.003 * Tc);
        const double d = co2->rhomass(), u = co2->umass();
        BENCHMARK("CO2 near-crit density D+U") {
            co2->update(CoolProp::DmassUmass_INPUTS, d, u);
            return co2->T();
        };
    }
}

// ---------------------------------------------------------------------------
// Component profile of the two-phase residual call(): how is the time split
// between the superancillary Chebyshev evals and the EOS property evaluations?
//     CatchTestRunner "[HSU_D_prof]"
// ---------------------------------------------------------------------------
TEST_CASE("HSU_D: two-phase call() component profile", "[HSU_D_prof][.]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto* be = dynamic_cast<CoolProp::HelmholtzEOSMixtureBackend*>(AS.get());
    if (be == nullptr || !be->is_pure()) return;  // not a pure HEOS backend -> no superancillary
    std::shared_ptr<CoolProp::EquationOfState::SuperAncillary_t> sa;
    try {
        sa = be->get_superanc();
    } catch (...) {
        sa.reset();  // treat a failed lookup as "no superancillary"
    }
    if (!sa) return;  // Water lacks a superancillary -> nothing to profile
    const double T = 450.0;
    const double rhoL = sa->eval_sat(T, 'D', 0), rhoV = sa->eval_sat(T, 'D', 1), p = sa->eval_sat(T, 'P', 1);

    BENCHMARK("3x eval_sat (Chebyshev only)") {
        return sa->eval_sat(T, 'D', 0) + sa->eval_sat(T, 'D', 1) + sa->eval_sat(T, 'P', 1);
    };
    BENCHMARK("1x update_TDmolarP_unchecked (no prop)") {
        be->SatL->update_TDmolarP_unchecked(T, rhoL, p);
        return be->SatL->rhomolar();
    };
    BENCHMARK("1x keyed_output(Hmolar) after set") {
        be->SatL->update_TDmolarP_unchecked(T, rhoL, p);
        return be->SatL->keyed_output(CoolProp::iHmolar);
    };
    BENCHMARK("FULL 2phase residual body (3 eval_sat + 2 set + 2 keyed_output)") {
        const double rl = sa->eval_sat(T, 'D', 0), rv = sa->eval_sat(T, 'D', 1), pp = sa->eval_sat(T, 'P', 1);
        be->SatL->update_TDmolarP_unchecked(T, rl, pp);
        be->SatV->update_TDmolarP_unchecked(T, rv, pp);
        const double yL = be->SatL->keyed_output(CoolProp::iHmolar);
        const double yV = be->SatV->keyed_output(CoolProp::iHmolar);
        return yL + yV + rl + rv;
    };
    BENCHMARK("baseline: update_DmolarT_direct + hmolar (1 EOS eval)") {
        be->update_DmolarT_direct(900.0, T);
        return be->keyed_output(CoolProp::iHmolar);
    };
    // Setup + commit costs that bracket the rootfind in the real flash.
    const double rho_2ph = 0.5 * (rhoL + rhoV);  // a two-phase density at T=450
    BENCHMARK("get_all_intersections('D', rho)") {
        return static_cast<double>(sa->get_all_intersections('D', rho_2ph, 48, 100, 1e-13).size());
    };
    BENCHMARK("commit: update(QT_INPUTS, Q, T)") {
        AS->update(CoolProp::QT_INPUTS, 0.5, T);
        return AS->T();
    };
}

#endif  // ENABLE_CATCH
