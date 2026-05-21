// Catch2 tests for the SVDSBTL backend (Phase 2c).  Exercise the
// public CoolProp::AbstractState::factory path, the high-level
// PropsSI surface, and accuracy vs HEOS at probe points.

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

using Catch::Approx;

#    include <algorithm>
#    include <filesystem>
#    include <memory>
#    include <random>
#    include <string>
#    include <vector>

#    include "AbstractState.h"
#    include "Configuration.h"
#    include "CoolProp.h"
#    include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"
#    include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#    include "DataStructures.h"

namespace {

// Returns true if a cache file matching the SVDSBTL on-disk pattern
// exists for both HmassP_INPUTS and PT_INPUTS for `fluid`.
//
// We can't just call default_cache_path() and stat — the backend
// embeds an FNV-1a 64 opthash over the canonical options blob in the
// filename, and that hash isn't exposed from SVDSBTLBackend.cpp (it
// lives in an anonymous-namespace helper).  Glob the cache dir for
//
//   <fluid>.<source>.<input_pair_name>.*.svd.bin.z
//
// instead, which is robust to whatever opthash the build path happens
// to compute for the no-options default.
bool cache_present_for(const std::string& fluid, const std::string& source_backend = "HEOS") {
    namespace cp_sbtl = CoolProp::sbtl;
    namespace fs = std::filesystem;
    const fs::path dir = cp_sbtl::SVDSurfaceSerializer::default_cache_dir();
    std::error_code ec;
    if (!fs::exists(dir, ec)) return false;
    const auto match = [&](::CoolProp::input_pairs pair) {
        const std::string prefix = fluid + "." + source_backend + "." + CoolProp::get_input_pair_short_desc(pair) + ".";
        // Match the full SVDSBTL suffix, not just .z — the cache dir is
        // shared with BicubicBackend's .bin.z files and could accumulate
        // other .z artifacts.  Checking the literal contract avoids
        // false positives across unrelated extensions.
        static const std::string kSuffix = ".svd.bin.z";
        for (const auto& entry : fs::directory_iterator(dir, ec)) {
            const std::string name = entry.path().filename().string();
            const bool has_prefix = name.rfind(prefix, 0) == 0 && name.size() > prefix.size();
            const bool has_suffix = name.size() >= kSuffix.size() && name.compare(name.size() - kSuffix.size(), kSuffix.size(), kSuffix) == 0;
            if (has_prefix && has_suffix) {
                return true;
            }
        }
        return false;
    };
    return match(::CoolProp::HmassP_INPUTS) && match(::CoolProp::PT_INPUTS);
}

}  // namespace

TEST_CASE("SVDSBTL backend constructs via AbstractState::factory", "[SVDSBTL][factory][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    REQUIRE(AS != nullptr);
    REQUIRE(AS->backend_name() == "SVDSBTLBackend");
    REQUIRE(AS->fluid_names().size() == 1);
    REQUIRE(AS->fluid_names()[0] == "Water");
    // Cache-write failures are non-fatal by design (the backend
    // explicitly catches and swallows them in ensure_surface_), so
    // a missing cache file isn't an error -- the in-memory surface
    // still works.  Log the absence as INFO and only require a
    // non-null backend.
    if (!cache_present_for("Water")) {
        WARN("SVDSBTL Water cache not present after construction -- in-memory only (likely read-only home dir)");
    }
}

TEST_CASE("SVDSBTL backend reports critical / triple constants from HEOS reference", "[SVDSBTL][constants][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    REQUIRE(AS->T_critical() == Approx(heos->T_critical()).epsilon(1e-12));
    REQUIRE(AS->p_critical() == Approx(heos->p_critical()).epsilon(1e-12));
    REQUIRE(AS->Ttriple() == Approx(heos->Ttriple()).epsilon(1e-12));
    REQUIRE(AS->p_triple() == Approx(heos->p_triple()).epsilon(1e-12));
    REQUIRE(AS->molar_mass() == Approx(heos->molar_mass()).epsilon(1e-12));
    REQUIRE(AS->Tmax() == Approx(heos->Tmax()).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend PT lookup matches HEOS within tolerance", "[SVDSBTL][pt][water][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Compressed-liquid spot check matching Phase 2b's water PT preset assertion.
    const double p = 1.0e6;
    const double T = 350.0;
    AS->update(CoolProp::PT_INPUTS, p, T);
    heos->update(CoolProp::PT_INPUTS, p, T);
    INFO("rho HEOS=" << heos->rhomass() << "  SVDSBTL=" << AS->rhomass());
    INFO("h   HEOS=" << heos->hmass() << "    SVDSBTL=" << AS->hmass());
    REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
    REQUIRE(AS->hmass() == Approx(heos->hmass()).epsilon(5e-3));
    REQUIRE(AS->smass() == Approx(heos->smass()).epsilon(5e-3));
    REQUIRE(AS->umass() == Approx(heos->umass()).epsilon(5e-3));

    // Stored T / p round-trip the user's inputs verbatim.
    REQUIRE(AS->T() == Approx(T).epsilon(1e-12));
    REQUIRE(AS->p() == Approx(p).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend PH lookup matches HEOS within tolerance", "[SVDSBTL][ph][water][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Build a (p, T) probe via HEOS, then ask SVDSBTL for the same point via (h, p).
    const double p = 5.0e5;
    const double T = 400.0;
    heos->update(CoolProp::PT_INPUTS, p, T);
    const double h = heos->hmass();

    AS->update(CoolProp::HmassP_INPUTS, h, p);
    INFO("at (p, T)=(" << p << ", " << T << ")  h=" << h);
    INFO("rho HEOS=" << heos->rhomass() << "  SVDSBTL=" << AS->rhomass());
    INFO("T   HEOS=" << T << "  SVDSBTL=" << AS->T());
    REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
    REQUIRE(AS->T() == Approx(T).epsilon(5e-3));
    REQUIRE(AS->smass() == Approx(heos->smass()).epsilon(5e-3));

    // The user-supplied h flows through directly.
    REQUIRE(AS->hmass() == Approx(h).epsilon(1e-12));
    REQUIRE(AS->p() == Approx(p).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend rejects mixtures", "[SVDSBTL][reject][mixture]") {
    using Catch::Matchers::ContainsSubstring;
    REQUIRE_THROWS_WITH(
      std::unique_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", std::vector<std::string>{"Water", "Ethanol"})),
      ContainsSubstring("pure-fluid"));
}

TEST_CASE("SVDSBTL backend is rejected from PropsSI by default", "[SVDSBTL][propsi][gate]") {
    // Make sure the gate is closed (it is by default; this just guards
    // against an earlier test forgetting to restore the flag).
    CoolProp::set_config_bool(ALLOW_SVDSBTL_IN_PROPSSI, false);
    // PropsSI catches internally and returns _HUGE on failure, with the
    // exception message left in get_global_param_string("errstring").
    const double v = CoolProp::PropsSI("D", "T", 350.0, "P", 1.0e6, "SVDSBTL&HEOS::Water");
    const std::string err = CoolProp::get_global_param_string("errstring");
    INFO("PropsSI returned " << v << "; errstring=" << err);
    REQUIRE(!ValidNumber(v));
    REQUIRE(err.find("high-level") != std::string::npos);
}

TEST_CASE("SVDSBTL backend usable via PropsSI when ALLOW_SVDSBTL_IN_PROPSSI is set", "[SVDSBTL][propsi][slow]") {
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    heos->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_truth = heos->rhomass();

    // Opt in for this case, then restore the default so subsequent
    // tests see the gate as closed.
    CoolProp::set_config_bool(ALLOW_SVDSBTL_IN_PROPSSI, true);
    const double rho_via_propsi = CoolProp::PropsSI("D", "T", 350.0, "P", 1.0e6, "SVDSBTL&HEOS::Water");
    CoolProp::set_config_bool(ALLOW_SVDSBTL_IN_PROPSSI, false);
    INFO("PropsSI = " << rho_via_propsi << "  HEOS = " << rho_truth);
    REQUIRE(rho_via_propsi == Approx(rho_truth).epsilon(5e-3));
}

TEST_CASE("SVDSBTL backend cache reload produces the same result", "[SVDSBTL][cache][slow]") {
    // First instance populates the cache (or no-ops if already present).
    auto first = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    if (!cache_present_for("Water")) {
        WARN("SVDSBTL cache absent after first construction (read-only home dir?); reload test will exercise rebuild path instead");
    }
    first->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_first = first->rhomass();

    // Second instance hits the cache.
    auto second = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    second->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_second = second->rhomass();

    REQUIRE(rho_first == Approx(rho_second).epsilon(1e-15));
}

TEST_CASE("SVDSBTL backend PQ_INPUTS two-phase blend matches HEOS", "[SVDSBTL][twophase][pq][water][slow]") {
    // SVDSBTL and HEOS both end up consulting the same SuperAncillary
    // for sat-line properties, so the agreement is bit-identical on
    // T and rho (the SA is the *source* HEOS uses internally for PQ
    // flashes on fluids that ship one).  h is computed via a Q-weighted
    // blend on both sides; tiny floating-point roundoff differences
    // (~1e-15 relative) are the only delta.  margin(...) lets us
    // express the tolerance in absolute terms where relative
    // comparisons break down at zero-crossing references.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Sweep a few (p, Q) combinations across the subcritical range.
    for (const double p : {1.0e5, 1.0e6, 5.0e6, 1.5e7}) {
        for (const double Q : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            AS->update(CoolProp::PQ_INPUTS, p, Q);
            heos->update(CoolProp::PQ_INPUTS, p, Q);
            INFO("p=" << p << "  Q=" << Q << "  T_AS=" << AS->T() << "  T_HEOS=" << heos->T());
            REQUIRE(AS->T() == Approx(heos->T()).margin(1e-9));
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).margin(1e-9));
            // h has a near-zero reference state at the triple point,
            // so use an absolute margin (in J/kg).
            REQUIRE(AS->hmass() == Approx(heos->hmass()).margin(1e-6));
            REQUIRE(AS->smass() == Approx(heos->smass()).margin(1e-9));
            // p / Q flow through verbatim.
            REQUIRE(AS->p() == Approx(p).margin(0.0));
            REQUIRE(AS->Q() == Approx(Q).margin(0.0));
        }
    }
}

TEST_CASE("SVDSBTL backend QT_INPUTS two-phase blend matches HEOS", "[SVDSBTL][twophase][qt][water][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    for (const double T : {280.0, 350.0, 450.0, 600.0}) {
        for (const double Q : {0.0, 0.5, 1.0}) {
            AS->update(CoolProp::QT_INPUTS, Q, T);
            heos->update(CoolProp::QT_INPUTS, Q, T);
            INFO("T=" << T << "  Q=" << Q);
            REQUIRE(AS->p() == Approx(heos->p()).margin(1e-3));
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).margin(1e-9));
            REQUIRE(AS->hmass() == Approx(heos->hmass()).margin(1e-6));
            REQUIRE(AS->T() == Approx(T).margin(0.0));
            REQUIRE(AS->Q() == Approx(Q).margin(0.0));
        }
    }
}

TEST_CASE("SVDSBTL backend HmassP_INPUTS dome-hit routes to two-phase blend", "[SVDSBTL][twophase][hp_dome][water][slow]") {
    // Pick a sat state, take h half-way between hL and hV at p_sat, and
    // confirm SVDSBTL recovers the same Q as HEOS via the in-dome
    // HmassP path.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    for (const double p : {2.0e5, 1.0e6, 5.0e6}) {
        heos->update(CoolProp::PQ_INPUTS, p, 0.0);
        const double hL = heos->hmass();
        heos->update(CoolProp::PQ_INPUTS, p, 1.0);
        const double hV = heos->hmass();
        const double h_mid = 0.5 * (hL + hV);
        AS->update(CoolProp::HmassP_INPUTS, h_mid, p);
        heos->update(CoolProp::HmassP_INPUTS, h_mid, p);
        INFO("p=" << p << "  h_mid=" << h_mid);
        REQUIRE(AS->phase() == CoolProp::iphase_twophase);
        REQUIRE(AS->Q() == Approx(heos->Q()).margin(1e-9));
        REQUIRE(AS->T() == Approx(heos->T()).margin(1e-9));
        REQUIRE(AS->rhomass() == Approx(heos->rhomass()).margin(1e-9));
        REQUIRE(AS->hmass() == Approx(h_mid).margin(1e-9));  // input round-trip
    }
}

TEST_CASE("SVDSBTL backend speed of sound matches HEOS in single phase", "[SVDSBTL][speed_sound][water][slow]") {
    // Speed of sound is the 5th property on both PH and PT surfaces.
    // For typical single-phase Water states, SVDSBTL should agree with
    // HEOS to ~SVD fitting tolerance (median ~1e-7 relative for ρ-like
    // quantities; w roughly tracks that).
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    struct Pt
    {
        double p;
        double T;
        const char* label;
    };
    for (const auto& pt : std::vector<Pt>{
           {1.0e5, 320.0, "subcooled liquid, near-ambient"},
           {1.0e6, 400.0, "subcooled liquid, mid-p"},
           {5.0e6, 350.0, "compressed liquid"},
           {1.0e5, 500.0, "superheated vapor, low-p"},
           {1.0e6, 700.0, "superheated vapor, mid-p"},
           {3.0e7, 800.0, "supercritical-ish (just past p_crit)"},
         }) {
        AS->update(CoolProp::PT_INPUTS, pt.p, pt.T);
        heos->update(CoolProp::PT_INPUTS, pt.p, pt.T);
        INFO(pt.label << "  p=" << pt.p << " T=" << pt.T << "  w_AS=" << AS->speed_sound() << "  w_HEOS=" << heos->speed_sound());
        REQUIRE(AS->speed_sound() == Approx(heos->speed_sound()).epsilon(5e-3));  // 0.5% — well under IF97's 0.1% budget at most points
    }
}

TEST_CASE("SVDSBTL backend speed of sound throws in two-phase", "[SVDSBTL][speed_sound][twophase]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    AS->update(CoolProp::PQ_INPUTS, 1.0e6, 0.5);
    REQUIRE_THROWS(AS->speed_sound());
}

TEST_CASE("SVDSBTL backend works across the Phase 2a fluid set", "[SVDSBTL][multi_fluid][slow]") {
    // For each Phase 2a fluid, exercise the full SVDSBTL pipeline:
    // PT single-phase, HmassP round-trip, PQ two-phase (via the
    // SuperAncillary cascade), and speed_sound.  Looser tolerances
    // than the Water-only suite because some fluids have less
    // smooth h_sat curves and the cache build runs at production
    // resolution (NT=200, NR=800, rank=20) -- a 0.5% bar still
    // confirms each backend wires up correctly.
    for (const auto* fluid : {"R134a", "Ammonia", "Methane", "Propane", "CarbonDioxide", "D6"}) {
        SECTION(fluid) {
            auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
            auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
            REQUIRE(AS->backend_name() == "SVDSBTLBackend");

            // Single-phase PT probe in compressed liquid.  0.85 * T_crit
            // keeps us comfortably above the melting line for CO2 (which
            // has T_melt(p_crit/2) > T_triple by ~1 K) while still well
            // into the liquid phase across the rest of the fluid set.
            const double p_pt = 0.5 * heos->p_critical();
            const double T_pt = 0.85 * heos->T_critical();
            AS->update(CoolProp::PT_INPUTS, p_pt, T_pt);
            heos->update(CoolProp::PT_INPUTS, p_pt, T_pt);
            INFO(fluid << " PT  rho_AS=" << AS->rhomass() << "  rho_HEOS=" << heos->rhomass());
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
            REQUIRE(AS->hmass() == Approx(heos->hmass()).epsilon(5e-3));

            // HmassP round-trip from that same (p, T) state.
            const double h_pt = heos->hmass();
            AS->update(CoolProp::HmassP_INPUTS, h_pt, p_pt);
            INFO(fluid << " PH  T_AS=" << AS->T() << "  T_HEOS=" << T_pt);
            REQUIRE(AS->T() == Approx(T_pt).epsilon(5e-3));
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
            REQUIRE(AS->hmass() == Approx(h_pt).epsilon(1e-12));  // input passes through

            // Speed of sound at the same point.
            AS->update(CoolProp::PT_INPUTS, p_pt, T_pt);
            INFO(fluid << " w  AS=" << AS->speed_sound() << "  HEOS=" << heos->speed_sound());
            REQUIRE(AS->speed_sound() == Approx(heos->speed_sound()).epsilon(1e-2));  // 1% bar > IF97 budget but absorbs any fluid-specific quirks

            // Two-phase PQ at half-critical pressure.  All Phase 2a
            // fluids ship a SuperAncillary, so the two-phase path
            // must succeed -- match HEOS to machine precision.
            const double p_pq = 0.3 * heos->p_critical();
            AS->update(CoolProp::PQ_INPUTS, p_pq, 0.5);
            heos->update(CoolProp::PQ_INPUTS, p_pq, 0.5);
            INFO(fluid << " PQ  rho_AS=" << AS->rhomass() << "  rho_HEOS=" << heos->rhomass());
            REQUIRE(AS->T() == Approx(heos->T()).margin(1e-9));
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).margin(1e-9));
            REQUIRE(AS->Q() == Approx(0.5).margin(0.0));
        }
    }
}

TEST_CASE("SVDSBTL backend Q out of [0, 1] is rejected", "[SVDSBTL][twophase][reject]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    REQUIRE_THROWS(AS->update(CoolProp::PQ_INPUTS, 1.0e6, -0.1));
    REQUIRE_THROWS(AS->update(CoolProp::PQ_INPUTS, 1.0e6, 1.1));
    REQUIRE_THROWS(AS->update(CoolProp::QT_INPUTS, -0.1, 350.0));
    REQUIRE_THROWS(AS->update(CoolProp::QT_INPUTS, 1.1, 350.0));
}

// -----------------------------------------------------------------------------
// Non-HEOS source-of-truth backends
// -----------------------------------------------------------------------------

TEST_CASE("SVDSBTL backend without explicit source throws", "[SVDSBTL][source][reject]") {
    // factory("SVDSBTL", ...) without &<source> must throw.  The error
    // text mentions SVDSBTL&HEOS so the user has an obvious fix.
    REQUIRE_THROWS_WITH(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water")),
                        Catch::Matchers::ContainsSubstring("explicit source"));
}

TEST_CASE("SVDSBTL backend rejects unsupported source names", "[SVDSBTL][source][reject]") {
    REQUIRE_THROWS_WITH(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&SRK", "Water")),
                        Catch::Matchers::ContainsSubstring("unsupported source"));
}

TEST_CASE("SVDSBTL&IF97 rejects non-Water fluids", "[SVDSBTL][source][if97][reject]") {
    REQUIRE_THROWS_WITH(std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "CarbonDioxide")),
                        Catch::Matchers::ContainsSubstring("Water-only"));
}

// Try-and-skip helper: REFPROP and IF97 may not be available on every
// CI machine.  Construct once at TEST_CASE entry and SKIP if the
// source-backend itself can't be built.
namespace {
bool source_backend_available(const std::string& source, const std::string& fluid) {
    try {
        const std::string spec = std::string("SVDSBTL&") + source;
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(spec, fluid));
        return true;
    } catch (const std::exception&) {
        return false;
    }
}
}  // namespace

TEST_CASE("SVDSBTL&REFPROP single-phase PT lookup matches REFPROP truth", "[SVDSBTL][source][refprop][pt][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    const double T = 0.85 * refprop->T_critical();
    const double p = 0.50 * refprop->p_critical();
    AS->update(CoolProp::PT_INPUTS, p, T);
    refprop->update(CoolProp::PT_INPUTS, p, T);
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-3));
    REQUIRE(AS->hmass() == Approx(refprop->hmass()).epsilon(1e-3));
}

TEST_CASE("SVDSBTL&REFPROP two-phase PQ matches REFPROP truth", "[SVDSBTL][source][refprop][twophase][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP two-phase test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    const double p_sat = 0.4 * refprop->p_critical();
    AS->update(CoolProp::PQ_INPUTS, p_sat, 0.5);
    refprop->update(CoolProp::PQ_INPUTS, p_sat, 0.5);
    // Two-phase: SVDSBTL routes through the source backend's QT_INPUTS
    // for sat endpoints (REFPROP has no SuperAncillary), then does the
    // standard Q-weighted blend.  Should agree to within REFPROP's own
    // numerical precision on the sat curves.
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-8));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-6));
    REQUIRE(AS->Q() == Approx(0.5));
}

TEST_CASE("SVDSBTL&IF97 single-phase PT lookup matches IF97 truth", "[SVDSBTL][source][if97][pt][slow]") {
    if (!source_backend_available("IF97", "Water")) {
        SKIP("IF97 backend not available; skipping SVDSBTL&IF97 test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    auto if97 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("IF97", "Water"));
    // 500 K, 3 MPa — IAPWS Table 5 region (compressed liquid)
    const double T = 500.0;
    const double p = 3.0e6;
    AS->update(CoolProp::PT_INPUTS, p, T);
    if97->update(CoolProp::PT_INPUTS, p, T);
    REQUIRE(AS->rhomass() == Approx(if97->rhomass()).epsilon(1e-3));
    REQUIRE(AS->hmass() == Approx(if97->hmass()).epsilon(1e-3));
}

TEST_CASE("SVDSBTL&HEOS and SVDSBTL&REFPROP produce distinct cache files", "[SVDSBTL][source][cache]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping cache-key disambiguation test");
    }
    namespace cp_sbtl = CoolProp::sbtl;
    const std::string heos_path = cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", CoolProp::HmassP_INPUTS);
    const std::string refprop_path = cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "REFPROP", CoolProp::HmassP_INPUTS);
    REQUIRE(heos_path != refprop_path);
    REQUIRE(heos_path.find("Water.HEOS.") != std::string::npos);
    REQUIRE(refprop_path.find("Water.REFPROP.") != std::string::npos);
}

TEST_CASE("SVDSBTL fast_evaluate matches per-call update() bit-for-bit", "[SVDSBTL][fast_evaluate][water][slow]") {
    // Build a small (T, p) probe set via IF97, pull h(T, p), then evaluate
    // each probe via both update() + property reads AND via fast_evaluate;
    // assert the two paths return identical numbers (no rounding tolerance).
    auto if97 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("IF97", "Water"));
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    std::vector<double> h_v, p_v;
    for (double T = 320.0; T <= 800.0; T += 60.0) {
        for (double lp = std::log(1.0e5); lp <= std::log(3.0e7); lp += 0.5 * std::log(10.0)) {
            const double p = std::exp(lp);
            try {
                if97->update(CoolProp::PT_INPUTS, p, T);
                const double h = if97->hmass();
                if (std::isfinite(h)) {
                    h_v.push_back(h);
                    p_v.push_back(p);
                }
            } catch (...) {  // NOLINT(bugprone-empty-catch)
            }
        }
    }
    REQUIRE(h_v.size() > 10);
    const std::size_t N = h_v.size();
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iT, CoolProp::iDmass, CoolProp::iSmass, CoolProp::ispeed_sound};
    std::vector<double> out_fe(N * outputs.size(), std::nan(""));
    std::vector<int> status(N, -1);
    svd->fast_evaluate(CoolProp::HmassP_INPUTS, h_v.data(), p_v.data(), N, outputs.data(), outputs.size(), out_fe.data(), out_fe.size(),
                       status.data(), status.size());
    for (std::size_t k = 0; k < N; ++k) {
        REQUIRE(status[k] == CoolProp::fast_evaluate_ok);
        svd->update(CoolProp::HmassP_INPUTS, h_v[k], p_v[k]);
        const double T_up = svd->T();
        const double rho_up = svd->rhomass();
        const double s_up = svd->smass();
        const double w_up = svd->speed_sound();
        // Bit-identical: both paths call the same SVDSurface::eval_with_region
        // on the same ResolvedPoint, so there's no floating-point freedom
        // between them.  No epsilon tolerance.
        REQUIRE(out_fe[k * outputs.size() + 0] == T_up);
        REQUIRE(out_fe[k * outputs.size() + 1] == rho_up);
        REQUIRE(out_fe[k * outputs.size() + 2] == s_up);
        REQUIRE(out_fe[k * outputs.size() + 3] == w_up);
    }
}

TEST_CASE("SVDSBTL fast_evaluate Q-blends inside the saturation dome", "[SVDSBTL][fast_evaluate][twophase][water][slow]") {
    // Probe inside the dome at three subcritical pressures, vary Q.
    // fast_evaluate must (a) report fast_evaluate_ok, (b) recover Q,
    // and (c) blend mass-basis properties linearly between the
    // sat-line endpoints.
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    const std::vector<double> ps = {1.0e5, 1.0e6, 5.0e6};
    const std::vector<double> qs = {0.05, 0.25, 0.5, 0.75, 0.95};
    std::vector<double> h_v, p_v, q_v;
    for (double p : ps) {
        ref->update(CoolProp::PQ_INPUTS, p, 0.0);
        const double hL = ref->hmass();
        ref->update(CoolProp::PQ_INPUTS, p, 1.0);
        const double hV = ref->hmass();
        for (double q : qs) {
            h_v.push_back(hL + q * (hV - hL));
            p_v.push_back(p);
            q_v.push_back(q);
        }
    }
    const std::size_t N = h_v.size();
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iT, CoolProp::iDmass, CoolProp::iSmass, CoolProp::iUmass, CoolProp::iQ};
    std::vector<double> out(N * outputs.size(), std::nan(""));
    std::vector<int> status(N, -1);
    svd->fast_evaluate(CoolProp::HmassP_INPUTS, h_v.data(), p_v.data(), N, outputs.data(), outputs.size(), out.data(), out.size(), status.data(),
                       status.size());
    for (std::size_t k = 0; k < N; ++k) {
        INFO("k=" << k << "  p=" << p_v[k] << "  h=" << h_v[k]);
        REQUIRE(status[k] == CoolProp::fast_evaluate_ok);
        // Q recovered to within a few ULP of the lever-rule construction;
        // any residual is the HEOS sat-line round-trip noise (~1e-5).
        REQUIRE(out[k * outputs.size() + 4] == Approx(q_v[k]).epsilon(1e-4));
        // Q-blend matches update_two_phase_() bit-identically — both paths
        // use the same SuperAncillary endpoints.
        svd->update(CoolProp::HmassP_INPUTS, h_v[k], p_v[k]);
        REQUIRE(out[k * outputs.size() + 0] == Approx(svd->T()).epsilon(1e-12));
        REQUIRE(out[k * outputs.size() + 1] == Approx(svd->rhomass()).epsilon(1e-12));
        REQUIRE(out[k * outputs.size() + 2] == Approx(svd->smass()).epsilon(1e-12));
        REQUIRE(out[k * outputs.size() + 3] == Approx(svd->umass()).epsilon(1e-12));
    }
}

TEST_CASE("SVDSBTL fast_evaluate returns Q = -1 on single-phase rows", "[SVDSBTL][fast_evaluate][water]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    // Three solid single-phase points well inside the surface envelope:
    // compressed liquid, sub-saturated gas, supercritical.  All should
    // report iQ = -1.0 (matching update()'s sentinel for single-phase).
    const std::vector<double> h_v = {5.0e5, 3.0e6, 2.0e6};
    const std::vector<double> p_v = {1.0e7, 1.0e5, 3.0e7};
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iT, CoolProp::iQ};
    std::vector<double> out(h_v.size() * outputs.size(), std::nan(""));
    std::vector<int> status(h_v.size(), -1);
    svd->fast_evaluate(CoolProp::HmassP_INPUTS, h_v.data(), p_v.data(), h_v.size(), outputs.data(), outputs.size(), out.data(), out.size(),
                       status.data(), status.size());
    for (std::size_t k = 0; k < h_v.size(); ++k) {
        INFO("k=" << k << " p=" << p_v[k] << " h=" << h_v[k]);
        REQUIRE(status[k] == CoolProp::fast_evaluate_ok);
        REQUIRE(out[k * outputs.size() + 1] == -1.0);
    }
}

TEST_CASE("SVDSBTL fast_evaluate routes through patch_source_ inside the critical-patch bbox", "[SVDSBTL][fast_evaluate][water][slow]") {
    // The default auto-calibrated patch bbox for Water sits at
    // [0.95, 1.05] * T_c x [0.75, 1.15] * p_c.  Pick a state well
    // inside that box and assert fast_evaluate's per-output values
    // match what patch_source_->update + property reads return on
    // the same input.  This exercises the Patched-kind branch in
    // resolve_point_ + evaluate_property_ that nothing else in the
    // suite hits.
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto ref = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    const double T_c = ref->T_critical();
    const double p_c = ref->p_critical();
    // 99% T_c, 95% p_c — comfortably inside [0.95, 1.05] x [0.75, 1.15].
    const double T_probe = 0.99 * T_c;
    const double p_probe = 0.95 * p_c;
    // Get h(p, T) from HEOS so we can route the probe through both
    // PT_INPUTS (direct) and HmassP_INPUTS (which exercises the same
    // bbox via the (p, h) projection).
    ref->update(CoolProp::PT_INPUTS, p_probe, T_probe);
    const double h_probe = ref->hmass();
    for (const auto& routing : std::vector<std::pair<CoolProp::input_pairs, std::pair<double, double>>>{
           {CoolProp::PT_INPUTS, {p_probe, T_probe}},
           {CoolProp::HmassP_INPUTS, {h_probe, p_probe}},
         }) {
        const double v1 = routing.second.first;
        const double v2 = routing.second.second;
        const std::vector<CoolProp::parameters> outputs = {CoolProp::iT, CoolProp::iDmass, CoolProp::iSmass, CoolProp::iHmass, CoolProp::iQ};
        std::vector<double> out(outputs.size(), std::nan(""));
        std::vector<int> status(1, -1);
        svd->fast_evaluate(routing.first, &v1, &v2, 1, outputs.data(), outputs.size(), out.data(), out.size(), status.data(), status.size());
        REQUIRE(status[0] == CoolProp::fast_evaluate_ok);
        // Bit-identical to update() + property reads — both code paths
        // route through the same patch_source_ AbstractState.
        svd->update(routing.first, v1, v2);
        INFO("input_pair=" << routing.first << " v1=" << v1 << " v2=" << v2);
        REQUIRE(out[0] == svd->T());
        REQUIRE(out[1] == svd->rhomass());
        REQUIRE(out[2] == svd->smass());
        REQUIRE(out[3] == svd->hmass());
        // iQ on a single-phase supercritical probe should mirror the
        // patch_source_'s answer (typically iphase_supercritical, Q
        // returned as a sentinel-shaped finite value from HEOS).
        REQUIRE(out[4] == svd->Q());
    }
}

TEST_CASE("SVDSBTL fast_evaluate flags out-of-range PT misses cleanly", "[SVDSBTL][fast_evaluate][water]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    // PT below the IF97 triple-point floor — genuinely outside every
    // region of the atlas.  Must report fast_evaluate_out_of_range
    // (not two-phase-disallowed, which is reserved for atlas misses
    // that the sat-curve probe can identify; the SVDSBTL PT atlas
    // registers LIQUID and VAPOR as adjacent at T_sat so a sat-curve
    // hit resolves as liquid rather than missing).
    const std::vector<double> v1 = {1.0e6};
    const std::vector<double> v2 = {200.0};  // K, below T_triple
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iDmass};
    std::vector<double> out(1, std::nan(""));
    std::vector<int> status(1, -1);
    svd->fast_evaluate(CoolProp::PT_INPUTS, v1.data(), v2.data(), 1, outputs.data(), outputs.size(), out.data(), out.size(), status.data(),
                       status.size());
    REQUIRE(status[0] == CoolProp::fast_evaluate_out_of_range);
    REQUIRE(!std::isfinite(out[0]));
}

TEST_CASE("SVDSBTL fast_evaluate rejects unsupported inputs cleanly", "[SVDSBTL][fast_evaluate][reject]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    const std::vector<double> v1 = {1.0};
    const std::vector<double> v2 = {1.0};
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iT};
    std::vector<double> out(1, std::nan(""));
    std::vector<int> status(1, -1);
    // PQ_INPUTS isn't a fast-path input — caller should get a clear error,
    // not a NaN row.
    REQUIRE_THROWS_AS(svd->fast_evaluate(CoolProp::PQ_INPUTS, v1.data(), v2.data(), 1, outputs.data(), outputs.size(), out.data(), out.size(),
                                         status.data(), status.size()),
                      CoolProp::ValueError);
}

#endif  // ENABLE_CATCH
