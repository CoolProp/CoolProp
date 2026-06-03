// Catch2 tests for the SVDSBTL critical-patch HEOS-fallback routing
// (CoolProp-p0p).  At construction time, the backend computes a
// per-input-pair bounding box around the critical point; query-time
// update() flips a flag if the active state lies inside; every calc_*
// then forwards to a source backend instead of the SVD surface.
//
// The default-mode (no options supplied) uses spec-recommended
// multipliers [0.95, 1.05] * Tcrit  x  [0.75, 1.15] * pcrit.  PR D
// replaces the constants with a real auto-calibration loop.

#if defined(ENABLE_CATCH)

#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <memory>
#    include <string>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Exceptions.h"

namespace {

std::shared_ptr<CoolProp::AbstractState> make_svdsbtl(const std::string& opts = "") {
    const std::string spec = opts.empty() ? std::string("Water") : "Water?" + opts;
    return std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", spec));
}

}  // namespace

TEST_CASE("critical_patch: query inside the bbox matches source backend", "[SVDSBTL][critical_patch]") {
    auto AS = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // A state that lives inside the default critical bbox:
    // T = T_crit (= 647.096 K) and p = 1.05 * p_crit.  Without the
    // patch the rank-truncation residual is multi-percent; with the
    // patch the answer is the source backend's value byte-for-byte.
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();
    const double T = Tc;
    const double p = 1.05 * pc;

    AS->update(CoolProp::PT_INPUTS, p, T);
    heos->update(CoolProp::PT_INPUTS, p, T);

    // Patch is active → SVDSBTL must return exactly what HEOS would
    // return (the patch source IS HEOS in this configuration).
    REQUIRE(AS->rhomass() == Catch::Approx(heos->rhomass()).margin(1e-9));
    REQUIRE(AS->hmass() == Catch::Approx(heos->hmass()).margin(1e-6));
    REQUIRE(AS->smass() == Catch::Approx(heos->smass()).margin(1e-9));
    REQUIRE(AS->speed_sound() == Catch::Approx(heos->speed_sound()).margin(1e-9));
}

TEST_CASE("critical_patch: query outside the bbox uses the SVD surface", "[SVDSBTL][critical_patch]") {
    auto AS = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // 0.5 * Tcrit  /  0.5 * pcrit — comfortably outside the default
    // patch bbox in BOTH axes.  Should be within SVD-truncation
    // accuracy (a few hundred ppm on density) of HEOS, but NOT
    // bit-identical (which would mean the patch fired by mistake).
    const double T = 0.5 * heos->T_critical();
    const double p = 0.5 * heos->p_critical();

    AS->update(CoolProp::PT_INPUTS, p, T);
    heos->update(CoolProp::PT_INPUTS, p, T);
    const double rho_svd = AS->rhomass();
    const double rho_heos = heos->rhomass();
    REQUIRE(std::isfinite(rho_svd));
    REQUIRE(rho_svd == Catch::Approx(rho_heos).epsilon(2e-3));  // SVD truncation looseness
    // Bit-identity check — if these were equal, the patch fired here.
    if (rho_svd == rho_heos) {
        FAIL("critical-patch fired outside its bbox at a clearly subcritical state");
    }
}

TEST_CASE("critical_patch: mode=off disables routing", "[SVDSBTL][critical_patch]") {
    auto AS = make_svdsbtl(R"({"critical_patch":{"mode":"off"}})");
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Same point that fires the patch in the default config; with
    // mode=off the SVD surface answers, and the value is NOT
    // bit-identical to HEOS (it carries SVD-truncation residual).
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();
    AS->update(CoolProp::PT_INPUTS, 1.05 * pc, Tc);
    heos->update(CoolProp::PT_INPUTS, 1.05 * pc, Tc);
    const double rho_svd = AS->rhomass();
    const double rho_heos = heos->rhomass();
    REQUIRE(std::isfinite(rho_svd));
    if (rho_svd == rho_heos) {
        FAIL("critical-patch is supposed to be off but produced bit-identical values");
    }
}

TEST_CASE("critical_patch: HmassP_INPUTS routing matches source backend in the bbox", "[SVDSBTL][critical_patch]") {
    auto AS = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Same critical-region state as above, but expressed in (h, p).
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();
    heos->update(CoolProp::PT_INPUTS, 1.05 * pc, Tc);
    const double h = heos->hmass();
    const double p = 1.05 * pc;

    AS->update(CoolProp::HmassP_INPUTS, h, p);
    heos->update(CoolProp::HmassP_INPUTS, h, p);

    REQUIRE(AS->rhomass() == Catch::Approx(heos->rhomass()).margin(1e-9));
    REQUIRE(AS->T() == Catch::Approx(heos->T()).margin(1e-9));
    REQUIRE(AS->smass() == Catch::Approx(heos->smass()).margin(1e-9));
}

TEST_CASE("critical_patch: HEOS source must NOT polish (regression for low-Tc drift)", "[SVDSBTL][critical_patch][hydrogen][slow]") {
    // polish_patch_state_ was added (d176979e7) to fix IF97's R7-97
    // backward-equation floor (±25 mK forward-consistency).  HEOS's
    // HmassP_INPUTS is already iterative + forward-consistent, so the
    // polish should be a no-op for HEOS sources.
    //
    // It isn't.  For low-Tc fluids (Hydrogen, Tc=33.14 K) the polish
    // bracket [T_seed±0.5 K] is ~1.5% of Tc, and TOMS748's eps_tolerance
    // ≈ 1e-12 relative T leaves a sub-ULP residual which the critical-
    // region stiffness amplifies into ~1e-7 ρ drift inside the patch
    // bbox.  Visible as non-zero error in 82% of bbox cells in the
    // SVDSBTLValidation heatmap, even though the patch source IS the
    // reference HEOS and the answer should be bit-exact.
    //
    // Pin: every property inside the patch must equal HEOS's HmassP
    // result bit-for-bit (margin 0) for HEOS-sourced backends.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(
      CoolProp::AbstractState::factory("SVDSBTL&HEOS", R"(Hydrogen?{"critical_patch":{"mode":"fixed","bbox":[0.95,1.05,0.90,1.10]}})"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Hydrogen"));
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();

    // (T=Tc, p=1.05*pc) — supercritical, comfortably inside the fixed
    // bbox in (T, p).  Translate to (h, p) via HEOS to get the HmassP
    // probe.
    heos->update(CoolProp::PT_INPUTS, 1.05 * pc, Tc);
    const double h = heos->hmass();
    const double p = 1.05 * pc;

    AS->update(CoolProp::HmassP_INPUTS, h, p);
    heos->update(CoolProp::HmassP_INPUTS, h, p);

    // Bit-exact (==): the patch routes to the same HEOS instance, so
    // any difference is polish-induced drift.
    REQUIRE(AS->rhomass() == heos->rhomass());
    REQUIRE(AS->T() == heos->T());
    REQUIRE(AS->smass() == heos->smass());
}

TEST_CASE("critical_patch: HmolarP_INPUTS reaches the surface lookup (CodeRabbit)", "[SVDSBTL][critical_patch]") {
    // surfaces_ only registers HmassP_INPUTS and PT_INPUTS — HmolarP
    // shares the HmassP surface via a molar→mass conversion inside
    // update().  Without normalising the surface_key the surfaces_.find()
    // would miss and update(HmolarP_INPUTS, ...) would throw before
    // any of the routing logic ran.
    auto AS = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double T = 300.0, p = 1e6;
    heos->update(CoolProp::PT_INPUTS, p, T);
    const double h_molar = heos->hmolar();
    REQUIRE_NOTHROW(AS->update(CoolProp::HmolarP_INPUTS, h_molar, p));
    REQUIRE(AS->p() == Catch::Approx(p).margin(1.0));
    REQUIRE(AS->T() == Catch::Approx(T).epsilon(1e-3));
}

TEST_CASE("critical_patch: invalid critical_patch.source fails at construction (CodeRabbit)", "[SVDSBTL][critical_patch]") {
    // The patch source is materialised inside build_critical_patch_()
    // so a misconfigured source (here: IF97 patch on a non-Water
    // fluid; IF97 is Water-only) fails immediately at construction
    // instead of on the first in-patch query.
    REQUIRE_THROWS(
      std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", R"(R134a?{"critical_patch":{"source":"IF97"}})")));
}

TEST_CASE("critical_patch: mode=fixed honours an explicit bbox", "[SVDSBTL][critical_patch]") {
    // A tight box ~ Tcrit ± 1%, pcrit ± 3% — narrower than the
    // default; a sample inside this fires the patch, a sample just
    // outside does not.
    auto AS = make_svdsbtl(R"({"critical_patch":{"mode":"fixed","bbox":[0.99,1.01,0.97,1.03]}})");
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();

    SECTION("inside narrow bbox → patch fires") {
        AS->update(CoolProp::PT_INPUTS, pc, Tc);
        heos->update(CoolProp::PT_INPUTS, pc, Tc);
        REQUIRE(AS->rhomass() == Catch::Approx(heos->rhomass()).margin(1e-9));
    }
    SECTION("outside narrow bbox but inside default → patch off") {
        // 1.10 * pc is inside the default [0.75, 1.15] bbox but
        // outside this 1.03 fixed cap.  With the fixed-mode bbox the
        // patch must NOT fire here.
        AS->update(CoolProp::PT_INPUTS, 1.10 * pc, Tc);
        heos->update(CoolProp::PT_INPUTS, 1.10 * pc, Tc);
        const double rho_svd = AS->rhomass();
        const double rho_heos = heos->rhomass();
        REQUIRE(std::isfinite(rho_svd));
        if (rho_svd == rho_heos) {
            FAIL("fixed-mode patch fired outside its supplied bbox");
        }
    }
}

TEST_CASE("critical_patch: auto-calibration produces sane multipliers", "[SVDSBTL][critical_patch][auto_calibration][slow]") {
    // Default mode is "auto", which (post-CoolProp-5ni/dxd) calibrates
    // the bbox per fluid by shrinking each axis from the Water-default
    // (0.95, 1.05, 0.75, 1.15) toward (1.0, 1.0, 1.0, 1.0) until the
    // SVD's reconstruction at the strip just outside the patch passes
    // a relaxed (1% rel-err) budget against the source backend.  The
    // result is persisted to <Fluid>.<Source>.critpatch.<OptHash>.bin
    // alongside the .svd.bin.z files.
    //
    // Black-box assertions: T_lo / p_lo must be < 1.0 (the patch
    // extends below the critical point on those axes), T_hi / p_hi
    // must be > 1.0 (extends above), and no axis can be wider than
    // the Water default (the calibration only shrinks).
    auto AS = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();

    SECTION("Patch is centered on the critical point with sensible width") {
        // The patch covers SOME slice around (Tc, pc); without a
        // direct accessor for the multipliers, we probe corner-to-
        // corner: the patch must fire at (Tc, pc) (centered) and
        // must NOT fire at points clearly outside the Water default
        // (0.90*Tc, 0.70*pc).
        AS->update(CoolProp::PT_INPUTS, pc, Tc);
        heos->update(CoolProp::PT_INPUTS, pc, Tc);
        const double rho_in_AS = AS->rhomass();
        const double rho_in_heos = heos->rhomass();
        // At the critical point, the patch fires and bit-equals HEOS.
        REQUIRE(rho_in_AS == rho_in_heos);

        // Outside the *Water default* bbox the patch must not fire,
        // regardless of how aggressive the auto-calibration was — the
        // calibrator only shrinks, never widens.
        AS->update(CoolProp::PT_INPUTS, 0.50 * pc, 0.85 * Tc);
        heos->update(CoolProp::PT_INPUTS, 0.50 * pc, 0.85 * Tc);
        const double rho_out_AS = AS->rhomass();
        const double rho_out_heos = heos->rhomass();
        REQUIRE(std::isfinite(rho_out_AS));
        // SVD-served result must differ from HEOS truth (the SVD has
        // some interpolation residual; ~1e-4 typical for this point).
        if (rho_out_AS == rho_out_heos) {
            FAIL("patch fired outside the Water-default bbox; calibrator widened?");
        }
    }
}

TEST_CASE("critical_patch: auto-calibration cache round-trip", "[SVDSBTL][critical_patch][auto_calibration]") {
    // Construct twice in succession; the second construction must
    // find the sidecar file and skip recalibration.  Functional test:
    // both instances must produce bit-identical answers inside the
    // patch (proves both loaded the same multipliers).
    auto AS1 = make_svdsbtl();
    auto AS2 = make_svdsbtl();
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();

    AS1->update(CoolProp::PT_INPUTS, pc, Tc);
    AS2->update(CoolProp::PT_INPUTS, pc, Tc);
    REQUIRE(AS1->rhomass() == AS2->rhomass());
    REQUIRE(AS1->hmass() == AS2->hmass());

    // Outside the patch, both instances must also agree (SVD shared
    // between the two — but importantly, no inconsistency from one
    // running the calibrator and the other loading the cache).
    AS1->update(CoolProp::PT_INPUTS, 0.40 * pc, 0.90 * Tc);
    AS2->update(CoolProp::PT_INPUTS, 0.40 * pc, 0.90 * Tc);
    REQUIRE(AS1->rhomass() == AS2->rhomass());
}

TEST_CASE("critical_patch: explicit bbox override skips auto-calibration", "[SVDSBTL][critical_patch][auto_calibration]") {
    // mode=auto with an explicit bbox is an escape hatch: the user
    // wants the auto path *off* without flipping mode to "fixed"
    // (which has different semantics — fixed disables the per-fluid
    // tightening but also rejects the cache load).  An explicit bbox
    // in auto mode must be honoured verbatim, no calibration probe.
    auto AS = make_svdsbtl(R"({"critical_patch":{"mode":"auto","bbox":[0.99, 1.01, 0.90, 1.10]}})");
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double Tc = heos->T_critical();
    const double pc = heos->p_critical();

    // At T = 0.985*Tc, p = 1.05*pc — inside the Water default (0.95)
    // but OUTSIDE the user's narrow [0.99, 1.01] T-bbox — patch must
    // NOT fire here.  Stay slightly above pc so we're solidly in the
    // SVD's SUPER region (the subcritical/supercritical knife-edge
    // at exactly p_crit can leave the SVD atlas without a covering
    // region).
    AS->update(CoolProp::PT_INPUTS, 1.05 * pc, 0.985 * Tc);
    heos->update(CoolProp::PT_INPUTS, 1.05 * pc, 0.985 * Tc);
    REQUIRE(std::isfinite(AS->rhomass()));
    if (AS->rhomass() == heos->rhomass()) {
        FAIL("user-supplied bbox in auto mode was ignored");
    }

    // At T = Tc, p = pc — inside the user's narrow bbox — patch fires.
    AS->update(CoolProp::PT_INPUTS, pc, Tc);
    heos->update(CoolProp::PT_INPUTS, pc, Tc);
    REQUIRE(AS->rhomass() == heos->rhomass());
}

#endif  // ENABLE_CATCH
