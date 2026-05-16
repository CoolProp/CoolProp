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
#    include "CoolProp.h"
#    include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"
#    include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#    include "DataStructures.h"

namespace {

// Returns true if both files at the SVDSBTL default cache path exist
// for `fluid`.  Used as a "did the constructor populate the cache?"
// gate without re-running the build path itself.
bool cache_present_for(const std::string& fluid, const std::string& source_backend = "HEOS") {
    namespace cp_sbtl = CoolProp::sbtl;
    return std::filesystem::exists(cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid, source_backend, CoolProp::HmassP_INPUTS))
           && std::filesystem::exists(cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid, source_backend, CoolProp::PT_INPUTS));
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

TEST_CASE("SVDSBTL backend supports the high-level PropsSI surface", "[SVDSBTL][propsi][slow]") {
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    heos->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_truth = heos->rhomass();

    const double rho_via_propsi = CoolProp::PropsSI("D", "T", 350.0, "P", 1.0e6, "SVDSBTL&HEOS::Water");
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

#endif  // ENABLE_CATCH
