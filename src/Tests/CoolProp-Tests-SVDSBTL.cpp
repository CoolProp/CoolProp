// Catch2 tests for the SVDSBTL backend (Phase 2c).  Exercise the
// public CoolProp::AbstractState::factory path, the high-level
// PropsSI surface, and accuracy vs HEOS at probe points.

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

using Catch::Approx;

#    include <algorithm>
#    include <atomic>
#    include <chrono>
#    include <cmath>
#    include <cstdint>
#    include <filesystem>
#    include <fstream>
#    include <memory>
#    include <random>
#    include <string>
#    include <thread>
#    include <vector>

#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/CoolProp.h"
#    include "CoolProp/Backends/SVDSBTL/SVDSBTLBackend.h"
#    include "CoolProp/sbtl/SatBoundaryFactory.h"
#    include "CoolProp/sbtl/SVDSurfaceFactory.h"
#    include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#    include "CoolProp/detail/atomic_write.h"
#    include "CoolProp/DataStructures.h"
#    include "TestUtils.h"

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

TEST_CASE("SVDSBTL backend resolves PT down to the triple-point pressure (issue #3189)", "[SVDSBTL][pt][low_pressure][issue_3189][slow]") {
    // Regression for issue #3189 / CoolProp-naqt.  The subcritical PT
    // surface used to floor its lowest isobar at ~1.01-1.03·p_triple
    // (subcritical_pressure_range returned p_triple*1.01, and the
    // p_triple it used was itself a QT sample at Ttriple*1.001, ~1-2%
    // high).  A low-pressure gas query at p == p_triple() with T well
    // above the saturation temperature therefore landed below the grid
    // and returned NaN.  The floor is now the true triple-point pressure
    // (heos.p_triple()), so p == p_triple() resolves and matches HEOS.
    const std::vector<std::string> fluids = {"Water", "Nitrogen", "Methane", "Hydrogen"};
    for (const auto& fluid : fluids) {
        SECTION(fluid) {
            auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
            auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));

            const double p_tri = AS->p_triple();
            const double T = 400.0;  // gas: T >> T_sat(p_triple) for all four fluids

            // Exactly at the triple pressure — the boundary that used to fail.
            AS->update(CoolProp::PT_INPUTS, p_tri, T);
            heos->update(CoolProp::PT_INPUTS, p_tri, T);
            INFO(fluid << " at p=p_triple=" << p_tri << " T=" << T);
            INFO("rho HEOS=" << heos->rhomass() << "  SVDSBTL=" << AS->rhomass());
            INFO("h   HEOS=" << heos->hmass() << "    SVDSBTL=" << AS->hmass());
            REQUIRE(std::isfinite(AS->hmass()));
            REQUIRE(std::isfinite(AS->rhomass()));
            REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
            REQUIRE(AS->hmass() == Approx(heos->hmass()).epsilon(5e-3));

            // A point just inside the floor (1% above p_triple) also resolves.
            const double p_hi = 1.01 * p_tri;
            AS->update(CoolProp::PT_INPUTS, p_hi, T);
            heos->update(CoolProp::PT_INPUTS, p_hi, T);
            REQUIRE(AS->hmass() == Approx(heos->hmass()).epsilon(5e-3));

            // Strictly below the triple line remains out of the tabulated
            // domain by design (p_triple is the documented floor): a
            // property read on an out-of-range point is NaN.
            AS->update(CoolProp::PT_INPUTS, 0.9 * p_tri, T);
            REQUIRE_FALSE(std::isfinite(AS->hmass()));
        }
    }
}

TEST_CASE("SVDSBTL pmin option overrides the subcritical lower pressure bound", "[SVDSBTL][pmin][options][slow]") {
    // The `pmin` option (absolute Pa) raises the PT/PH/PS table floor.
    // Use Water (p_triple ~= 611 Pa) with a small grid to keep the fresh
    // surface build cheap — this test exercises geometry, not accuracy.
    const std::string opts = R"({"pmin": 2000.0, "grid": {"NT": 40, "NR": 80, "rank": 8}})";
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water?" + opts));

    const double T = 400.0;
    // Below the raised floor (but above p_triple) -> out of range -> NaN.
    AS->update(CoolProp::PT_INPUTS, 1000.0, T);
    REQUIRE_FALSE(std::isfinite(AS->hmass()));
    // Above the floor -> resolves.
    AS->update(CoolProp::PT_INPUTS, 5000.0, T);
    REQUIRE(std::isfinite(AS->hmass()));
}

TEST_CASE("SVDSBTL pmin below the triple pressure is rejected", "[SVDSBTL][pmin][options][reject]") {
    // A sub-triple pmin has no liquid-vapour saturation boundary to bound
    // the subcritical regions, so construction (which eagerly builds the
    // PT surface) throws rather than failing obscurely deep in sampling.
    const std::string opts = R"({"pmin": 100.0})";  // < Water p_triple (~611 Pa)
    REQUIRE_THROWS_AS(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water?" + opts), CoolProp::ValueError);
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
            // h / u for a PQ dome point now flow through the LAZY DomeBlend
            // arms (the resolve no longer pre-fills enthalpy, CoolProp-qp0n).
            REQUIRE(AS->umass() == Approx(heos->umass()).margin(1e-6));
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
            // h / u via the LAZY DomeBlend arms (CoolProp-qp0n).
            REQUIRE(AS->umass() == Approx(heos->umass()).margin(1e-6));
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

TEST_CASE("SVDSBTL backend HmassP_INPUTS on the bubble/dew line stays two-phase (issue #3190)",
          "[SVDSBTL][twophase][hp_dome][regression][issue_3190][water][slow]") {
    // Round-trip h(p, Q) -> (p, h) on the saturation boundary itself.
    // The atlas single-phase regions' dome-side boundary is an
    // interpolated sat curve that overshoots the true bubble/dew line by
    // its fit error, so a point sitting EXACTLY on the boundary used to
    // resolve single-phase (Q=-1, iphase_not_imposed).  The near-dome
    // eta guard must reclassify it as two-phase with Q == Q_in.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));

    const double p_tri = AS->p_triple();
    const double p_crit = AS->p_critical();
    for (const double Q_in : {0.0, 1.0}) {
        for (int i = 1; i <= 12; ++i) {
            const double p = p_tri + (p_crit - p_tri) * (static_cast<double>(i) / 13.0);
            AS->update(CoolProp::PQ_INPUTS, p, Q_in);
            const double h = AS->hmass();
            AS->update(CoolProp::HmassP_INPUTS, h, p);
            INFO("Q_in=" << Q_in << "  p=" << p << "  h=" << h << "  -> phase=" << AS->phase() << "  Q=" << AS->Q());
            REQUIRE(AS->phase() == CoolProp::iphase_twophase);
            REQUIRE(AS->Q() == Approx(Q_in).margin(1e-6));
        }
    }
}

TEST_CASE("SVDSBTL backend PSmass_INPUTS on the bubble/dew line stays two-phase (issue #3190)",
          "[SVDSBTL][twophase][ps_dome][regression][issue_3190][water][slow]") {
    // Entropy twin of the HmassP bubble/dew-line round-trip above.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));

    const double p_tri = AS->p_triple();
    const double p_crit = AS->p_critical();
    for (const double Q_in : {0.0, 1.0}) {
        for (int i = 1; i <= 12; ++i) {
            const double p = p_tri + (p_crit - p_tri) * (static_cast<double>(i) / 13.0);
            AS->update(CoolProp::PQ_INPUTS, p, Q_in);
            const double s = AS->smass();
            AS->update(CoolProp::PSmass_INPUTS, p, s);
            INFO("Q_in=" << Q_in << "  p=" << p << "  s=" << s << "  -> phase=" << AS->phase() << "  Q=" << AS->Q());
            REQUIRE(AS->phase() == CoolProp::iphase_twophase);
            REQUIRE(AS->Q() == Approx(Q_in).margin(1e-6));
        }
    }
}

TEST_CASE("SVDSBTL backend bubble/dew-line round-trip across fluids (issue #3190, atlas-miss path)",
          "[SVDSBTL][twophase][regression][issue_3190][multi_fluid][slow]") {
    // Multi-fluid guard for the atlas-MISS branch of the #3190 fix.  The
    // Water-only cases above exercise the atlas-hit near-dome guard; these
    // non-water fluids surfaced low-p dew-line points that go through the
    // atlas-miss dome blend, where a Q=0/1 round-trip lands a ULP outside
    // [yL, yV] (the y_mass = yV/M then y_mol = y_mass*M re-scaling is not
    // bit-exact) and a strict containment check dropped them to OutOfRange
    // (reported Q=-inf).  Sweep p in [p_triple, p_critical] and require the
    // re-flash stays two-phase with Q == Q_in for both HmassP and PSmass.
    for (const char* fluid : {"n-Propane", "R245fa"}) {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
        const double p_tri = AS->p_triple();
        const double p_crit = AS->p_critical();
        // Sweep stops short of p_crit (max i/N = 0.975): as p -> p_crit the
        // dome (yV - yL) collapses and the lever-rule Q noise grows, which is
        // the critical patch's domain, not this dome-blend round-trip's.
        const int N = 40;
        for (const double Q_in : {0.0, 1.0}) {
            for (int i = 1; i < N; ++i) {
                const double p = p_tri + (p_crit - p_tri) * (static_cast<double>(i) / static_cast<double>(N));
                AS->update(CoolProp::PQ_INPUTS, p, Q_in);
                const double h = AS->hmass();
                AS->update(CoolProp::HmassP_INPUTS, h, p);
                INFO(fluid << " HmassP Q_in=" << Q_in << " p/pc=" << (p / p_crit) << " -> phase=" << AS->phase() << " Q=" << AS->Q());
                REQUIRE(AS->phase() == CoolProp::iphase_twophase);
                REQUIRE(AS->Q() == Approx(Q_in).margin(1e-6));

                AS->update(CoolProp::PQ_INPUTS, p, Q_in);
                const double s = AS->smass();
                AS->update(CoolProp::PSmass_INPUTS, p, s);
                INFO(fluid << " PSmass Q_in=" << Q_in << " p/pc=" << (p / p_crit) << " -> phase=" << AS->phase() << " Q=" << AS->Q());
                REQUIRE(AS->phase() == CoolProp::iphase_twophase);
                REQUIRE(AS->Q() == Approx(Q_in).margin(1e-6));
            }
        }
    }
}

TEST_CASE("SVDSBTL backend PS lookup matches HEOS within tolerance", "[SVDSBTL][ps][water][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // Build a (p, T) probe via HEOS, then ask SVDSBTL for the same point
    // via (p, s) — entropy twin of the PH lookup test.
    const double p = 5.0e5;
    const double T = 400.0;
    heos->update(CoolProp::PT_INPUTS, p, T);
    const double s = heos->smass();

    AS->update(CoolProp::PSmass_INPUTS, p, s);
    INFO("at (p, T)=(" << p << ", " << T << ")  s=" << s);
    INFO("rho HEOS=" << heos->rhomass() << "  SVDSBTL=" << AS->rhomass());
    INFO("T   HEOS=" << T << "  SVDSBTL=" << AS->T());
    REQUIRE(AS->rhomass() == Approx(heos->rhomass()).epsilon(5e-3));
    REQUIRE(AS->T() == Approx(T).epsilon(5e-3));
    REQUIRE(AS->hmass() == Approx(heos->hmass()).epsilon(5e-3));

    // The user-supplied s flows through directly; p round-trips.
    REQUIRE(AS->smass() == Approx(s).epsilon(1e-12));
    REQUIRE(AS->p() == Approx(p).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend PSmolar lookup matches HEOS within tolerance", "[SVDSBTL][ps][psmolar][water][slow]") {
    // Molar-basis mirror of the PSmass lookup — exercises the
    // s_molar -> s_mass conversion in resolve_point_ and the iSmolar /
    // iHmolar / iDmolar output scaling.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    const double p = 5.0e5;
    const double T = 400.0;
    heos->update(CoolProp::PT_INPUTS, p, T);
    const double s_molar = heos->smolar();

    AS->update(CoolProp::PSmolar_INPUTS, p, s_molar);
    INFO("at (p, T)=(" << p << ", " << T << ")  s_molar=" << s_molar);
    REQUIRE(AS->rhomolar() == Approx(heos->rhomolar()).epsilon(5e-3));
    REQUIRE(AS->T() == Approx(T).epsilon(5e-3));
    REQUIRE(AS->hmolar() == Approx(heos->hmolar()).epsilon(5e-3));

    // The user-supplied molar entropy flows through directly; p round-trips.
    REQUIRE(AS->smolar() == Approx(s_molar).epsilon(1e-12));
    REQUIRE(AS->p() == Approx(p).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend PSmass_INPUTS dome-hit routes to two-phase blend", "[SVDSBTL][twophase][ps_dome][water][slow]") {
    // Entropy twin of the HmassP dome test: take s half-way between sL
    // and sV at p_sat and confirm SVDSBTL recovers the same Q as HEOS.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    for (const double p : {2.0e5, 1.0e6, 5.0e6}) {
        heos->update(CoolProp::PQ_INPUTS, p, 0.0);
        const double sL = heos->smass();
        heos->update(CoolProp::PQ_INPUTS, p, 1.0);
        const double sV = heos->smass();
        const double s_mid = 0.5 * (sL + sV);
        AS->update(CoolProp::PSmass_INPUTS, p, s_mid);
        heos->update(CoolProp::PSmass_INPUTS, p, s_mid);
        INFO("p=" << p << "  s_mid=" << s_mid);
        REQUIRE(AS->phase() == CoolProp::iphase_twophase);
        REQUIRE(AS->Q() == Approx(heos->Q()).margin(1e-9));
        REQUIRE(AS->T() == Approx(heos->T()).margin(1e-9));
        REQUIRE(AS->rhomass() == Approx(heos->rhomass()).margin(1e-9));
        REQUIRE(AS->smass() == Approx(s_mid).margin(1e-9));  // input round-trip
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

TEST_CASE("SVDSBTL backend tracks HEOS for R32 vapor at low superheat (issue #2247)", "[SVDSBTL][regression][issue_2247][r32][low_superheat][slow]") {
    // Regression for GH #2247: BICUBIC / TTSE smear density and
    // viscosity for R32 vapor at low superheat because a rectangular
    // interpolation stencil straddles the saturation dome.  SVDSBTL
    // bounds its VAPOR region by the actual saturated-vapor curve, so a
    // point ~1 K above T_sat resolves to the VAPOR surface and is never
    // blended across the dome.  This test reproduces the issue's sweep:
    // for T_sat in [-30, +30] degC, evaluate the vapor 1 K above
    // saturation via PT_INPUTS and require SVDSBTL to track HEOS.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "R32"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "R32"));
    REQUIRE(AS->backend_name() == "SVDSBTLBackend");

    constexpr double kSuperheat = 1.0;  // K above saturation (the issue's "SH")
    double max_rel_rho = 0.0;
    double max_rel_mu = 0.0;
    // 61-point sweep matches the issue's np.linspace(243.15, 303.15, 61).
    for (int i = 0; i < 61; ++i) {
        const double T_sat = 243.15 + i * (303.15 - 243.15) / 60.0;
        // Saturation pressure from HEOS (the issue's SSP), then the vapor
        // state at T_sat + 1 K (the issue's rho at P = SSP, T = SST + SH).
        heos->update(CoolProp::QT_INPUTS, 1.0, T_sat);
        const double p_sat = heos->p();
        const double T = T_sat + kSuperheat;

        heos->update(CoolProp::PT_INPUTS, p_sat, T);
        const double rho_truth = heos->rhomass();
        const double mu_truth = heos->viscosity();

        AS->update(CoolProp::PT_INPUTS, p_sat, T);

        INFO("T_sat=" << T_sat << " K  p_sat=" << p_sat << " Pa  T=" << T << " K");
        INFO("rho HEOS=" << rho_truth << "  SVDSBTL=" << AS->rhomass());
        INFO("mu  HEOS=" << mu_truth << "  SVDSBTL=" << AS->viscosity());

        // The point is genuinely superheated vapor, not two-phase: the
        // dome-straddling bug would mis-route it here.
        REQUIRE(AS->phase() != CoolProp::iphase_twophase);

        // 0.5 % on density — far inside the (often double-digit) BICUBIC
        // deviations the issue reports, yet comfortably above SVDSBTL's
        // fit residual for smooth single-phase vapor.
        REQUIRE(AS->rhomass() == Approx(rho_truth).epsilon(5e-3));
        // Viscosity is the second property #2247 calls out.  Transport
        // correlations are less smooth than density near the sat curve,
        // so allow 1 %.
        REQUIRE(AS->viscosity() == Approx(mu_truth).epsilon(1e-2));

        max_rel_rho = std::max(max_rel_rho, std::abs(AS->rhomass() - rho_truth) / rho_truth);
        max_rel_mu = std::max(max_rel_mu, std::abs(AS->viscosity() - mu_truth) / mu_truth);
    }
    INFO("max rel-err over sweep: rho=" << max_rel_rho << "  mu=" << max_rel_mu);
    REQUIRE(max_rel_rho < 5e-3);
    REQUIRE(max_rel_mu < 1e-2);
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

TEST_CASE("SVDSBTL&IF97 single-phase PS lookup matches IF97 truth", "[SVDSBTL][source][if97][ps][slow]") {
    if (!source_backend_available("IF97", "Water")) {
        SKIP("IF97 backend not available; skipping SVDSBTL&IF97 PS test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    auto if97 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("IF97", "Water"));
    // 500 K, 3 MPa — IAPWS Table 5 region (compressed liquid).  Build s
    // from IF97, then query SVDSBTL via (p, s).
    const double T = 500.0;
    const double p = 3.0e6;
    if97->update(CoolProp::PT_INPUTS, p, T);
    const double s = if97->smass();
    AS->update(CoolProp::PSmass_INPUTS, p, s);
    REQUIRE(AS->rhomass() == Approx(if97->rhomass()).epsilon(1e-3));
    REQUIRE(AS->T() == Approx(T).epsilon(1e-3));
    REQUIRE(AS->hmass() == Approx(if97->hmass()).epsilon(1e-3));

    // Supercritical point across the R2/R3 kink: 700 K, 30 MPa (R3).
    const double T2 = 700.0;
    const double p2 = 3.0e7;
    if97->update(CoolProp::PT_INPUTS, p2, T2);
    const double s2 = if97->smass();
    AS->update(CoolProp::PSmass_INPUTS, p2, s2);
    REQUIRE(AS->rhomass() == Approx(if97->rhomass()).epsilon(1e-3));
    REQUIRE(AS->T() == Approx(T2).epsilon(1e-3));
}

TEST_CASE("SVDSBTL&REFPROP PSmass single-phase matches REFPROP truth", "[SVDSBTL][source][refprop][psmass][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP PSmass test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    // Subcritical liquid: T=400 K, p=10 MPa.  Compute s from REFPROP,
    // then feed (p, s) to SVDSBTL.
    refprop->update(CoolProp::PT_INPUTS, 1.0e7, 400.0);
    const double s_liq = refprop->smass();
    AS->update(CoolProp::PSmass_INPUTS, 1.0e7, s_liq);
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-3));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-3));
    // Superheated vapor: T=600 K, p=0.5 MPa.
    refprop->update(CoolProp::PT_INPUTS, 5.0e5, 600.0);
    const double s_vap = refprop->smass();
    AS->update(CoolProp::PSmass_INPUTS, 5.0e5, s_vap);
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-3));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-3));
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

// HmassP single-phase via SVDSBTL&REFPROP — the dominant call pattern,
// exercises a different routing path from PT (uses log-p as primary axis,
// h as secondary, atlas-resolves the LIQUID/VAPOR/SUPER region directly
// from the inputs).  REFPROP-truth comparison so we catch any HEOS-only
// assumption that slipped past the PT test.
TEST_CASE("SVDSBTL&REFPROP HmassP single-phase matches REFPROP truth", "[SVDSBTL][source][refprop][hmassp][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP HmassP test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    // Subcritical liquid (well inside LIQUID region): T=400 K, p=10 MPa.
    // Compute h from REFPROP, then feed (h, p) to SVDSBTL.  Both should
    // recover the same density / temperature.
    refprop->update(CoolProp::PT_INPUTS, 1.0e7, 400.0);
    const double h_liq = refprop->hmass();
    AS->update(CoolProp::HmassP_INPUTS, h_liq, 1.0e7);
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-3));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-3));
    // Superheated vapor (well inside VAPOR region): T=600 K, p=0.5 MPa.
    refprop->update(CoolProp::PT_INPUTS, 5.0e5, 600.0);
    const double h_vap = refprop->hmass();
    AS->update(CoolProp::HmassP_INPUTS, h_vap, 5.0e5);
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-3));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-3));
}

// HmassP dome-blend via SVDSBTL&REFPROP — the harder path.  Inputs (h, p)
// land inside the saturation dome; the source PQ fallback runs (REFPROP
// has no SuperAncillary), computes sat endpoints via PQ_INPUTS, and the
// lever-rule blend should match REFPROP's own two-phase Q determination.
TEST_CASE("SVDSBTL&REFPROP HmassP dome-blend matches REFPROP truth", "[SVDSBTL][source][refprop][hmassp][twophase][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP dome-blend test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    // Pick a sat pressure away from the critical point so the lever rule
    // isn't near-degenerate.  Q=0.3 puts us comfortably inside the dome
    // (h_L + 0.3*(h_V - h_L)).
    const double p_sat = 0.3 * refprop->p_critical();
    refprop->update(CoolProp::PQ_INPUTS, p_sat, 0.3);
    const double h_dome = refprop->hmass();
    AS->update(CoolProp::HmassP_INPUTS, h_dome, p_sat);
    REQUIRE(AS->T() == Approx(refprop->T()).epsilon(1e-6));
    REQUIRE(AS->Q() == Approx(0.3).epsilon(1e-6));
    REQUIRE(AS->rhomass() == Approx(refprop->rhomass()).epsilon(1e-4));
}

// fast_evaluate via SVDSBTL&REFPROP — exercises the batched path with the
// REFPROP-backed source.  The fast_evaluate kernel is source-agnostic
// (uses the SVDSurface directly + source for dome endpoints) so the
// principal risk is the REFPROP-via-source_reference_ path for dome
// queries; both single-phase and dome probes are included.
TEST_CASE("SVDSBTL&REFPROP fast_evaluate works for mixed single-phase + dome batch", "[SVDSBTL][source][refprop][fast_evaluate][slow]") {
    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping SVDSBTL&REFPROP fast_evaluate test");
    }
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto refprop = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water"));
    // Build a (h, p) probe set spanning LIQUID (subcooled), VAPOR
    // (superheated), and DOME (two-phase) so each Kind of the kernel is
    // exercised.
    std::vector<double> h_v, p_v;
    for (double T_K : {350.0, 600.0}) {
        for (double p_Pa : {1.0e5, 5.0e6}) {
            refprop->update(CoolProp::PT_INPUTS, p_Pa, T_K);
            h_v.push_back(refprop->hmass());
            p_v.push_back(p_Pa);
        }
    }
    // One dome probe: p = 0.3*p_crit, Q = 0.4.
    refprop->update(CoolProp::PQ_INPUTS, 0.3 * refprop->p_critical(), 0.4);
    h_v.push_back(refprop->hmass());
    p_v.push_back(0.3 * refprop->p_critical());
    const std::size_t N = h_v.size();
    const std::vector<CoolProp::parameters> outputs = {CoolProp::iT, CoolProp::iDmass, CoolProp::iQ};
    std::vector<double> out_buffer(N * outputs.size(), std::nan(""));
    std::vector<int> status(N, -1);
    AS->fast_evaluate(CoolProp::HmassP_INPUTS, h_v.data(), p_v.data(), N, outputs.data(), outputs.size(), out_buffer.data(), out_buffer.size(),
                      status.data(), status.size());
    // Compare against update() + per-property reads point-by-point.
    // Bit-exact would be ideal but fails on the dome-blend probe
    // because REFPROP's PQ flash on the two sat endpoints isn't
    // order-deterministic at the last bits — fast_evaluate's
    // resolve_point_ does both endpoints in one call, while
    // update() does them in a different order across the two
    // round-trips here.  4-ULP-class differences observed in
    // practice; use a generous-but-tight relative tolerance.
    for (std::size_t k = 0; k < N; ++k) {
        INFO("k=" << k << " p=" << p_v[k] << " h=" << h_v[k]);
        // Check status BEFORE comparing — a non-OK status leaves
        // out_buffer holding NaN / garbage, and the comparison
        // failure would mask the real issue.
        REQUIRE(status[k] == CoolProp::fast_evaluate_ok);
        AS->update(CoolProp::HmassP_INPUTS, h_v[k], p_v[k]);
        REQUIRE(out_buffer[k * 3 + 0] == Approx(AS->T()).epsilon(1e-10));        // iT
        REQUIRE(out_buffer[k * 3 + 1] == Approx(AS->rhomass()).epsilon(1e-10));  // iDmass
        // iQ: -1 sentinel for single-phase, [0,1] in dome (per
        // SVDSBTL's PointEvaluation convention; see SVDSBTLBackend.h
        // PointEvaluation::Q docs).  The dome probe inherits the same
        // last-bit divergence as T / rho via REFPROP's not-quite-
        // order-deterministic PQ flash on the sat endpoints (see the
        // top-of-loop comment) — Q = (h_mol - h_L_mol) / (h_V_mol -
        // h_L_mol) is a function of those endpoints, so the bit-exact
        // compare that worked for HEOS fails for REFPROP at ~4-ULP.
        // The -1.0 single-phase sentinel is bit-exact under Approx
        // with margin 0 (Approx default), so this margin only relaxes
        // the dome case.
        REQUIRE(out_buffer[k * 3 + 2] == Approx(AS->Q()).margin(1e-10));
    }
}

TEST_CASE("SVDSBTL fast_evaluate matches per-call update() bit-for-bit", "[SVDSBTL][fast_evaluate][water][slow]") {
    // Build a small (T, p) probe set via IF97, pull h(T, p), then evaluate
    // each probe via both update() + property reads AND via fast_evaluate;
    // assert the two paths return identical numbers (no rounding tolerance).
    auto if97 = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("IF97", "Water"));
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&IF97", "Water"));
    std::vector<double> h_v, p_v;
    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter,cert-flp30-c)
    for (double T = 320.0; T <= 800.0; T += 60.0) {
        // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter,cert-flp30-c)
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

// -----------------------------------------------------------------------------
// ALTERNATIVE_SVDTABLES_DIRECTORY config key (CoolProp-fhp)
// -----------------------------------------------------------------------------

// Fast resolver-only test: no table build.  Verifies the contract that
// default_cache_dir() consults the config key, creates the directory,
// normalizes the trailing separator, and falls back to ~/.CoolProp/SVDTables
// when the override is empty.  Both the .svd.bin.z surfaces and the
// .critpatch.bin sidecars route through default_cache_dir(), so this single
// check covers both file kinds by construction.
TEST_CASE("SVDSBTL ALTERNATIVE_SVDTABLES_DIRECTORY routes default_cache_dir", "[SVDSBTL][cache][fhp]") {
    namespace fs = std::filesystem;
    namespace cp_sbtl = CoolProp::sbtl;
    const std::string saved = CoolProp::get_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY);
    const fs::path tmpdir = fs::temp_directory_path() / ("coolprop_svdtables_fhp_resolver_" + std::to_string(CoolProp::tests::test_pid()));
    std::error_code ec;
    fs::remove_all(tmpdir, ec);

    // Override -> resolved path starts with override and ends with a separator.
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, tmpdir.string());
    const std::string resolved = cp_sbtl::SVDSurfaceSerializer::default_cache_dir();
    REQUIRE(resolved.rfind(tmpdir.string(), 0) == 0);
    REQUIRE_FALSE(resolved.empty());
    REQUIRE((resolved.back() == '/' || resolved.back() == '\\'));
    REQUIRE(fs::is_directory(tmpdir));

    // Trailing separator already present in override -> no double slash.
    const std::string with_slash = tmpdir.string() + "/";
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, with_slash);
    const std::string resolved2 = cp_sbtl::SVDSurfaceSerializer::default_cache_dir();
    REQUIRE(resolved2.find("//") == std::string::npos);

    // Empty -> default $HOME/.CoolProp/SVDTables path.
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, "");
    const std::string defaulted = cp_sbtl::SVDSurfaceSerializer::default_cache_dir();
    REQUIRE(defaulted.find("/.CoolProp/SVDTables") != std::string::npos);

    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, saved);
    fs::remove_all(tmpdir, ec);
}

// End-to-end test: with the override set to a fresh tmpdir, building an
// SVDSBTL AbstractState must land its cache files inside that tmpdir (not in
// $HOME), and a second construction with the same options must hit-and-load
// from there (cache file's mtime must NOT advance).  Uses the smallest schema-
// valid grid so the build is a few seconds; tagged [slow] so the default
// preflight tag sweep doesn't pay the cost on every iteration.
TEST_CASE("SVDSBTL ALTERNATIVE_SVDTABLES_DIRECTORY end-to-end build + reload", "[SVDSBTL][cache][fhp][slow]") {
    namespace fs = std::filesystem;
    const std::string saved = CoolProp::get_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY);
    const fs::path tmpdir = fs::temp_directory_path() / ("coolprop_svdtables_fhp_e2e_" + std::to_string(CoolProp::tests::test_pid()));
    std::error_code ec;
    fs::remove_all(tmpdir, ec);
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, tmpdir.string());

    // Smallest valid grid + critical_patch off so we don't pay the
    // bbox-calibration loop.  The options blob is the cache key (via
    // its FNV-1a 64 opthash) so first and second construction must use
    // an identical blob to share a cache file.
    const std::string opts = R"({"grid":{"NT":40,"NR":80,"rank":10},"critical_patch":{"mode":"off"}})";
    const std::string spec = "SVDSBTL&HEOS";
    const std::string spec_with_opts = std::string("Water?") + opts;

    auto cache_files_in_tmpdir = [&]() {
        std::vector<fs::path> out;
        if (!fs::is_directory(tmpdir, ec)) return out;
        for (const auto& entry : fs::directory_iterator(tmpdir, ec)) {
            const std::string name = entry.path().filename().string();
            if (name.find(".svd.bin.z") != std::string::npos) out.push_back(entry.path());
        }
        return out;
    };

    // First construction -> build + write cache files into tmpdir.
    {
        auto first = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(spec, spec_with_opts));
        first->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    }
    REQUIRE(fs::is_directory(tmpdir));
    auto files = cache_files_in_tmpdir();
    INFO("cache files in tmpdir: " << files.size());
    REQUIRE_FALSE(files.empty());
    // Stash mtimes as plain int64 nanoseconds: Catch2's stringifier on
    // macOS can't print std::filesystem::file_time_type (the underlying
    // duration uses __int128) so a failing REQUIRE on the raw type
    // produces an ambiguous operator<< compile error.
    auto to_ns = [](const fs::file_time_type& t) {
        return static_cast<std::int64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t.time_since_epoch()).count());
    };
    std::vector<std::int64_t> mtimes_before_ns;
    mtimes_before_ns.reserve(files.size());
    for (const auto& f : files)
        mtimes_before_ns.push_back(to_ns(fs::last_write_time(f)));

    // Second construction with the same options -> must read from cache,
    // not rebuild.  A rebuild would rewrite the .svd.bin.z files and
    // bump their mtime.
    {
        auto second = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(spec, spec_with_opts));
        second->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    }
    const auto files_after = cache_files_in_tmpdir();
    REQUIRE(files_after.size() == files.size());
    for (std::size_t i = 0; i < files.size(); ++i) {
        INFO("file: " << files[i].filename().string());
        REQUIRE(to_ns(fs::last_write_time(files[i])) == mtimes_before_ns[i]);
    }

    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, saved);
    fs::remove_all(tmpdir, ec);
}

// CoolProp-i7j: a DmassT (DT-indexed) SVDSBTL surface must round-trip
// through the on-disk cache.  The DT preset's rho_sat,L / rho_sat,V dome
// boundaries are region::SuperancillaryTemperatureBoundaryCurve on HEOS
// sources; before the serializer learned CurveKind::SUPERANCILLARY_T,
// pack_curve threw "unknown BoundaryCurve subclass" on them.  save_to_file
// swallows that throw, so the DT surface silently never cached and the
// dense SVD rebuilt every session.  This test catches that regression:
// pre-fix the first construction writes NO DmassT cache file (REQUIRE_FALSE
// fails), or any cache that lands fails to load and rebuilds (mtime bumps).
TEST_CASE("SVDSBTL DmassT surface caches + reloads from disk (CoolProp-i7j)", "[SBTL][SVDSBTL][dt][cache][slow]") {
    namespace fs = std::filesystem;
    const std::string saved = CoolProp::get_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY);
    const fs::path tmpdir = fs::temp_directory_path() / ("coolprop_svdtables_dt_e2e_" + std::to_string(CoolProp::tests::test_pid()));
    std::error_code ec;
    fs::remove_all(tmpdir, ec);
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, tmpdir.string());

    // Smallest valid grid + critical_patch off, same as the PT e2e test.
    const std::string opts = R"({"grid":{"NT":40,"NR":80,"rank":10},"critical_patch":{"mode":"off"}})";
    const std::string spec = "SVDSBTL&HEOS";
    const std::string spec_with_opts = std::string("Water?") + opts;

    // A valid single-phase compressed-liquid (D, T): take D from HEOS at
    // (5 MPa, 350 K) so the probe lands inside the LIQUID region.
    double D_query = 0.0;
    {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
        heos->update(CoolProp::PT_INPUTS, 5.0e6, 350.0);
        D_query = heos->rhomass();
    }

    auto dt_cache_files = [&]() {
        std::vector<fs::path> out;
        if (!fs::is_directory(tmpdir, ec)) return out;
        for (const auto& entry : fs::directory_iterator(tmpdir, ec)) {
            const std::string name = entry.path().filename().string();
            if (name.find("DmassT_INPUTS") != std::string::npos && name.find(".svd.bin.z") != std::string::npos) out.push_back(entry.path());
        }
        return out;
    };

    // First construction -> build + write the DmassT cache file.
    {
        auto first = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(spec, spec_with_opts));
        first->update(CoolProp::DmassT_INPUTS, D_query, 350.0);
    }
    auto files = dt_cache_files();
    INFO("DmassT cache files in tmpdir: " << files.size());
    REQUIRE_FALSE(files.empty());  // pre-fix: pack_curve threw, save swallowed -> no DT cache

    auto to_ns = [](const fs::file_time_type& t) {
        return static_cast<std::int64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(t.time_since_epoch()).count());
    };
    std::vector<std::int64_t> mtimes_before_ns;
    mtimes_before_ns.reserve(files.size());
    for (const auto& f : files)
        mtimes_before_ns.push_back(to_ns(fs::last_write_time(f)));

    // Second construction with identical options -> must reload from cache.
    // A failed load (e.g. unknown curve kind) would fall through to rebuild
    // and bump the mtime.
    {
        auto second = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(spec, spec_with_opts));
        second->update(CoolProp::DmassT_INPUTS, D_query, 350.0);
    }
    const auto files_after = dt_cache_files();
    REQUIRE(files_after.size() == files.size());
    for (std::size_t i = 0; i < files.size(); ++i) {
        INFO("file: " << files[i].filename().string());
        REQUIRE(to_ns(fs::last_write_time(files[i])) == mtimes_before_ns[i]);
    }

    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, saved);
    fs::remove_all(tmpdir, ec);
}

// CoolProp-4no.2: parallel docs-build invocations (joblib loky workers)
// race on writes to the same cache file, producing torn-write artifacts.
// ::write_bytes_atomic does write-temp + rename so readers never see a
// partial-write file on the visible path.  This test spawns many threads
// contending on the same target and asserts the post-write file matches
// exactly one writer's payload (last writer wins) — never a torn /
// interleaved mix.
TEST_CASE("write_bytes_atomic is race-safe across threads", "[SVDSBTL][cache][race][4no.2]") {
    namespace fs = std::filesystem;
    const fs::path tmpdir = fs::temp_directory_path() / ("coolprop_svdtables_atomic_race_" + std::to_string(CoolProp::tests::test_pid()));
    std::error_code ec;
    fs::remove_all(tmpdir, ec);
    fs::create_directories(tmpdir);
    const fs::path target = tmpdir / "race.bin";

    constexpr int kThreads = 16;
    // Payload large enough that an ofstream::write of one payload would
    // be observable mid-stream by another writer if the writes were not
    // serialized via rename.  Several pages worth of bytes.
    constexpr std::size_t kPayloadSize = static_cast<std::size_t>(1) << 14U;  // 16 KiB
    std::vector<std::vector<char>> payloads(kThreads);
    for (int i = 0; i < kThreads; ++i) {
        payloads[i].assign(kPayloadSize, static_cast<char>(i + 1));  // distinct fill bytes per thread
    }

    std::atomic<bool> go{false};
    std::vector<std::thread> threads;
    threads.reserve(kThreads);
    for (int i = 0; i < kThreads; ++i) {
        threads.emplace_back([&, i]() {
            while (!go.load(std::memory_order_acquire)) {
                std::this_thread::yield();
            }
            ::write_bytes_atomic(target, payloads[i].data(), payloads[i].size(), /*restrict_perms=*/false);
        });
    }
    go.store(true, std::memory_order_release);
    for (auto& t : threads) {
        t.join();
    }

    // Read back: file must exist and match EXACTLY one writer's payload.
    REQUIRE(fs::exists(target));
    std::ifstream in(target, std::ios::binary);
    REQUIRE(in.good());
    std::vector<char> result((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    REQUIRE(result.size() == kPayloadSize);

    bool matched_one = false;
    int matched_idx = -1;
    for (int i = 0; i < kThreads; ++i) {
        if (result == payloads[i]) {
            matched_one = true;
            matched_idx = i;
            break;
        }
    }
    // cppcheck-suppress shiftNegative // Catch2 INFO macro: << is ostream insertion, not bit-shift
    INFO("Result must match exactly one writer's payload; instead matched index " << matched_idx);
    REQUIRE(matched_one);

    // No leftover .tmp.* sibling files (write_bytes_atomic must rename or
    // unlink on its way out — never leak temps on success).
    int leftover_temps = 0;
    for (const auto& entry : fs::directory_iterator(tmpdir)) {
        if (entry.path().filename().string().find(".tmp.") != std::string::npos) {
            ++leftover_temps;
        }
    }
    REQUIRE(leftover_temps == 0);

    fs::remove_all(tmpdir, ec);
}

// -- SaturationSurrogate (CoolProp-077) -----------------------------------

#    include "CoolProp/sbtl/SaturationSurrogate.h"

// Standalone surrogate accuracy: build from HEOS Water (always available,
// has a SuperAncillary as ground truth elsewhere; here we treat HEOS's own
// QT_INPUTS as truth and check the spline interpolant against it).
//
// Bead acceptance: rel err <= 1e-6 on smooth-fluid probes away from T_crit.
// Probe set deliberately picks T in [0.55, 0.95] * T_crit to stay well
// clear of the sqrt-singularity near T_crit and the iteration-sensitive
// strip near T_triple.
TEST_CASE("SaturationSurrogate eval_sat matches source QT_INPUTS on HEOS water", "[SBTL][SVDSBTL][sat_surrogate][slow]") {
    auto src = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto surrogate = CoolProp::sbtl::SaturationSurrogate::build_from_source(*src);
    REQUIRE(surrogate != nullptr);
    REQUIRE(surrogate->n_knots() >= 4);

    const double T_crit = src->T_critical();
    std::mt19937 rng(2026'05'23);
    std::uniform_real_distribution<double> tfrac(0.55, 0.95);

    for (int trial = 0; trial < 32; ++trial) {
        const double T = tfrac(rng) * T_crit;
        for (int side : {0, 1}) {
            src->update(CoolProp::QT_INPUTS, static_cast<double>(side), T);
            const double p_ref = src->p();
            const double h_ref = src->hmolar();
            const double s_ref = src->smolar();
            const double rho_ref = src->rhomolar();
            const double u_ref = src->umolar();

            REQUIRE(surrogate->eval_sat(T, 'P', side) == Approx(p_ref).epsilon(1e-6));
            REQUIRE(surrogate->eval_sat(T, 'H', side) == Approx(h_ref).epsilon(1e-6));
            REQUIRE(surrogate->eval_sat(T, 'S', side) == Approx(s_ref).epsilon(1e-6));
            REQUIRE(surrogate->eval_sat(T, 'D', side) == Approx(rho_ref).epsilon(1e-6));
            REQUIRE(surrogate->eval_sat(T, 'U', side) == Approx(u_ref).epsilon(1e-6));
        }
    }
}

// Near-critical guard: at T = 0.99 * T_crit (well above the eval_sat
// accuracy test's 0.95 upper bound, but inside the surrogate's build
// range of 0.9999 * T_crit), the surrogate must still return finite,
// physically-sensible values.  The bead's accuracy contract doesn't
// hold here — the sqrt-singularity at T_crit pushes cubic-spline error
// well past 1e-6 — but the surrogate must not return NaN / inf or
// overshoot below zero, since callers blend rho/h/s/u under the lever
// rule and a non-finite endpoint corrupts the whole dome row.
TEST_CASE("SaturationSurrogate stays finite and positive near T_crit", "[SBTL][SVDSBTL][sat_surrogate][slow]") {
    auto src = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto surrogate = CoolProp::sbtl::SaturationSurrogate::build_from_source(*src);
    REQUIRE(surrogate != nullptr);
    const double T = 0.99 * src->T_critical();
    for (int side : {0, 1}) {
        for (char k : {'P', 'D', 'H', 'S', 'U'}) {
            const double y = surrogate->eval_sat(T, k, side);
            REQUIRE(std::isfinite(y));
            // P/D/H positive at the dome interior for water (H_L is
            // positive in the IAPWS reference state above triple liquid;
            // P and D are positive by physics).  U/S are signed but
            // finite-only here.
            if (k == 'P' || k == 'D' || k == 'H') {
                REQUIRE(y > 0.0);
            }
        }
    }
}

// Inverse: p -> T_sat.  The surrogate stores T(log p) so the inversion is
// just one spline evaluation.  Same well-behaved probe interval as above.
TEST_CASE("SaturationSurrogate get_T_from_p inverts source PQ_INPUTS on HEOS water", "[SBTL][SVDSBTL][sat_surrogate][slow]") {
    auto src = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    auto surrogate = CoolProp::sbtl::SaturationSurrogate::build_from_source(*src);
    REQUIRE(surrogate != nullptr);

    const double p_crit = src->p_critical();
    std::mt19937 rng(2026'05'24);
    std::uniform_real_distribution<double> log_pfrac(std::log(1e-3), std::log(0.95));

    for (int trial = 0; trial < 16; ++trial) {
        const double p = std::exp(log_pfrac(rng)) * p_crit;
        src->update(CoolProp::PQ_INPUTS, p, 0.0);
        const double T_ref = src->T();
        REQUIRE(surrogate->get_T_from_p(p) == Approx(T_ref).epsilon(1e-6));
    }
}

// Backend wiring: when SVDSBTL is built on a source without a SuperAncillary
// (REFPROP today), the backend must lazily build a SaturationSurrogate.
// HEOS source must NOT trigger the surrogate (the SA path is preferred).
TEST_CASE("SVDSBTLBackend uses surrogate when source has no SuperAncillary", "[SBTL][SVDSBTL][sat_surrogate][source]") {
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto* heos_be = dynamic_cast<CoolProp::SVDSBTLBackend*>(heos.get());
    REQUIRE(heos_be != nullptr);

    // Hot path that flushes through the sat resolver — pick a clean dome
    // probe so the backend definitely went through the sat-provider code.
    heos->update(CoolProp::PQ_INPUTS, 0.4 * heos->p_critical(), 0.5);
    REQUIRE_FALSE(heos_be->sat_surrogate_consulted());

    if (!source_backend_available("REFPROP", "Water")) {
        SKIP("REFPROP not available; skipping surrogate-wired-in assertion for REFPROP source");
    }
    auto rp = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&REFPROP", "Water"));
    auto* rp_be = dynamic_cast<CoolProp::SVDSBTLBackend*>(rp.get());
    REQUIRE(rp_be != nullptr);
    rp->update(CoolProp::PQ_INPUTS, 0.4 * rp->p_critical(), 0.5);
    REQUIRE(rp_be->sat_surrogate_consulted());
}

// CoolProp-rhb5: the opt-in `prebuild` option materializes every
// supported input-pair surface at construction instead of lazy-loading
// the secondary pairs (DmassT / PSmass) on first query.
TEST_CASE("SVDSBTL prebuild option eagerly builds all surfaces", "[SBTL][SVDSBTL][prebuild][slow]") {
    auto contains = [](const std::vector<CoolProp::input_pairs>& v, CoolProp::input_pairs p) { return std::find(v.begin(), v.end(), p) != v.end(); };
    // Small grid keeps the four dense SVD builds fast for a unit test.
    const char* small_grid = R"({"grid":{"NT":40,"NR":80,"rank":10}})";
    const char* small_grid_prebuild = R"({"prebuild":true,"grid":{"NT":40,"NR":80,"rank":10}})";

    // Default (lazy): only the eager pairs PT + HmassP are registered at
    // construction; DmassT / PSmass are absent until first queried.
    {
        auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", std::string("Water?") + small_grid));
        auto* be = dynamic_cast<CoolProp::SVDSBTLBackend*>(AS.get());
        REQUIRE(be != nullptr);
        const auto pairs = be->registered_input_pairs();
        REQUIRE(contains(pairs, CoolProp::PT_INPUTS));
        REQUIRE(contains(pairs, CoolProp::HmassP_INPUTS));
        REQUIRE_FALSE(contains(pairs, CoolProp::DmassT_INPUTS));
        REQUIRE_FALSE(contains(pairs, CoolProp::PSmass_INPUTS));
    }
    // prebuild=true: all four supported pairs are registered eagerly.
    {
        auto AS =
          std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", std::string("Water?") + small_grid_prebuild));
        auto* be = dynamic_cast<CoolProp::SVDSBTLBackend*>(AS.get());
        REQUIRE(be != nullptr);
        const auto pairs = be->registered_input_pairs();
        REQUIRE(contains(pairs, CoolProp::PT_INPUTS));
        REQUIRE(contains(pairs, CoolProp::HmassP_INPUTS));
        REQUIRE(contains(pairs, CoolProp::DmassT_INPUTS));
        REQUIRE(contains(pairs, CoolProp::PSmass_INPUTS));
    }
}

TEST_CASE("SVDSBTL&IF97 prebuild skips the unbuildable DmassT surface", "[SBTL][SVDSBTL][prebuild][if97][slow]") {
    if (!source_backend_available("IF97", "Water")) {
        SKIP("IF97 backend not available; skipping SVDSBTL&IF97 prebuild test");
    }
    auto contains = [](const std::vector<CoolProp::input_pairs>& v, CoolProp::input_pairs p) { return std::find(v.begin(), v.end(), p) != v.end(); };
    // DmassT can't be sampled on a (D,T) grid from IF97 (all-NaN matrix),
    // so prebuild builds PT + HmassP + PSmass but leaves DmassT out rather
    // than throwing at construction.
    auto AS = std::shared_ptr<CoolProp::AbstractState>(
      CoolProp::AbstractState::factory("SVDSBTL&IF97", std::string(R"(Water?{"prebuild":true,"grid":{"NT":40,"NR":80,"rank":10}})")));
    auto* be = dynamic_cast<CoolProp::SVDSBTLBackend*>(AS.get());
    REQUIRE(be != nullptr);
    const auto pairs = be->registered_input_pairs();
    REQUIRE(contains(pairs, CoolProp::PT_INPUTS));
    REQUIRE(contains(pairs, CoolProp::HmassP_INPUTS));
    REQUIRE(contains(pairs, CoolProp::PSmass_INPUTS));
    REQUIRE_FALSE(contains(pairs, CoolProp::DmassT_INPUTS));
}

// prebuild must NOT change the cache key: a {"prebuild":true} build has
// to share its serialized surfaces with a plain no-prebuild instance (the
// opthash is stripped of the build-eagerness flag).  Otherwise the docs
// pre-warm would cache under a different opthash than the panels load.
TEST_CASE("SVDSBTL prebuild shares the surface cache with a plain instance", "[SBTL][SVDSBTL][prebuild][cache][slow]") {
    namespace fs = std::filesystem;
    const std::string saved = CoolProp::get_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY);
    const fs::path tmpdir = fs::temp_directory_path() / ("coolprop_svdtables_prebuild_" + std::to_string(CoolProp::tests::test_pid()));
    std::error_code ec;
    fs::remove_all(tmpdir, ec);
    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, tmpdir.string());

    auto count_ps_files = [&]() {
        int n = 0;
        for (const auto& e : fs::directory_iterator(tmpdir, ec)) {
            const auto name = e.path().filename().string();
            if (name.find("PSmass_INPUTS") != std::string::npos && name.find(".svd.bin.z") != std::string::npos) ++n;
        }
        return n;
    };

    // critical_patch off so we don't pay the bbox-calibration loop; small
    // grid so the four builds are a few seconds.
    const char* prebuild_opts = R"({"prebuild":true,"critical_patch":{"mode":"off"},"grid":{"NT":40,"NR":80,"rank":10}})";
    const char* plain_opts = R"({"critical_patch":{"mode":"off"},"grid":{"NT":40,"NR":80,"rank":10}})";

    // Prebuild instance: writes one PSmass surface to the tmpdir.
    { auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", std::string("Water?") + prebuild_opts)); }
    REQUIRE(count_ps_files() == 1);

    // Plain instance (same grid, NO prebuild): a PS query must LOAD the
    // already-cached surface — the opthash strip means it looks for the
    // same file, so no second PSmass cache file appears.
    auto plain = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", std::string("Water?") + plain_opts));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    heos->update(CoolProp::PT_INPUTS, 5.0e5, 400.0);
    plain->update(CoolProp::PSmass_INPUTS, 5.0e5, heos->smass());
    REQUIRE(plain->rhomass() == Approx(heos->rhomass()).epsilon(1e-1));  // coarse grid; just confirms it resolved
    REQUIRE(count_ps_files() == 1);                                      // shared, not a second opthash

    CoolProp::set_config_string(ALTERNATIVE_SVDTABLES_DIRECTORY, saved);
    fs::remove_all(tmpdir, ec);
}

// CoolProp-4u9: the tiny p-strip immediately around pc (sub: p ∈
// [(1−1e-10)·pc, pc]; super: p ∈ [pc, (1+1e-10)·pc]) sits outside
// the NC sub-regions (which stop at (1−1e-10)·pc on the sub-side
// and start at (1+1e-10)·pc on the super-side) and outside SUPER
// (which starts at 1.1·pc when NC is enabled).  Its natural home
// is the critical patch — cells in the strip route to the source
// backend directly and should return HEOS-exact values.  This test
// guards against (a) patch HmassP envelope misses (e.g. perimeter
// walk skipping cells at T > T_max for fluids with T_max ≈ Tc + ε)
// and (b) runtime gaps where neither NC, SUPER, nor patch claims a
// cell.
TEST_CASE("SVDSBTL pc-strip routes through critical patch", "[SVDSBTL][CoolProp-4u9][pc_strip][slow]") {
    // Include R245fa (T_max/Tc = 1.030 — the narrowest envelope in
    // the supported set), Water (T_max/Tc = 3.51 — comfortable),
    // and CO2 (T_max/Tc = 6.58 — comfortable).  R245fa is the one
    // that exposed the T_hi_mult clamp bug; Water/CO2 are regression
    // anchors.
    for (const auto* fluid : {"R245fa", "Water", "CarbonDioxide"}) {
        SECTION(fluid) {
            auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
            auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
            const double pc = heos->p_critical();
            const double Tc = heos->T_critical();

            // Probe both sides of pc at p-offsets that step through
            // the strip.  Use T very close to Tc on each side so the
            // cell lands inside the patch's auto-cal'd (T, p) bbox.
            // (The auto-cal bisects T_lo aggressively; for R245fa
            // T_lo ends at 0.9992·Tc.  A T 0.5 K below Tc would fall
            // BELOW patch.T_lo and route to NC at its ξ=1 boundary
            // edge — 1e-6 error, not the machine-precision the strip
            // is supposed to give.  Testing physically-meaningful
            // critical-region cells means T near Tc anyway.)
            const double T_sub = (1.0 - 1.0e-4) * Tc;
            const double T_sup = (1.0 + 1.0e-4) * Tc;
            for (const double frac : {1e-6, 1e-7, 1e-8}) {
                for (const auto& side : {std::make_pair(pc * (1.0 - frac), T_sub), std::make_pair(pc * (1.0 + frac), T_sup)}) {
                    const double p = side.first;
                    const double T = side.second;
                    // Forward (p, T) -> h via HEOS truth.
                    heos->update(CoolProp::PT_INPUTS, p, T);
                    const double h = heos->hmass();
                    const double rho_truth = heos->rhomass();

                    // 1) PT lookup: SVD must route to patch and
                    //    return HEOS-exact rho.
                    svd->update(CoolProp::PT_INPUTS, p, T);
                    INFO(fluid << " pc-strip PT: p=" << p << " T=" << T);
                    REQUIRE_FALSE(std::isnan(svd->rhomass()));
                    REQUIRE(svd->rhomass() == Approx(rho_truth).epsilon(1e-12));

                    // 2) HmassP lookup: same cell via (h, p) must
                    //    also route to patch and return what HEOS
                    //    returns for the SAME (h, p) input.  This
                    //    is the path that fails for R245fa without
                    //    the T_hi_mult clamp.  Compare against
                    //    heos->HmassP, not heos->PT, because the
                    //    polish gate (#2966) means HEOS-source
                    //    patches no longer iterate T to match
                    //    rho_truth_PT — they just call HEOS's
                    //    HmassP_INPUTS, which is iterative+forward-
                    //    consistent but not bit-exact with PT in
                    //    the ill-conditioned critical region
                    //    (~1e-7 residual for R245fa).
                    heos->update(CoolProp::HmassP_INPUTS, h, p);
                    const double rho_truth_hp = heos->rhomass();
                    svd->update(CoolProp::HmassP_INPUTS, h, p);
                    INFO(fluid << " pc-strip HmassP: p=" << p << " h=" << h << " (T_truth=" << T << ")");
                    REQUIRE_FALSE(std::isnan(svd->rhomass()));
                    REQUIRE(svd->rhomass() == Approx(rho_truth_hp).epsilon(1e-12));
                }
            }
        }
    }
}

// -- DT-indexed surface (CoolProp-i7j) ------------------------------------

// Regression gate for CoolProp-4z79: the near-critical band [0.99,1.01]*Tc
// used to be a structural GAP in the DmassT atlas (the dome-bounded
// LIQUID/VAPOR regions can't reach Tc — rho_sat,L/V pinch to rho_c).  The
// NC sub-regions (POWER/POWER_LO T axis + an isobar-bounded region that
// spans Tc) now carry the *table* through the whole band, including the
// critical isotherm, at all densities — no HEOS fallback.  This test
// asserts BOTH: (1) no gap — every single-phase (T,rho) in the band is
// served (svd->update does not throw); (2) accuracy — p(rho,T) matches
// HEOS to <1e-4 (the observed worst is ~2e-5; the IAPWS budget is 2e-4).
// Helium is included on purpose: its critical-region PT solver fails
// within ~1 nanoK of Tc, which earlier aborted the surface build (the
// 1e-6 NC handoff margin fixes it).
TEST_CASE("SVDSBTL DT near-critical band is gap-free and accurate (CoolProp-4z79)", "[SBTL][SVDSBTL][dt][near_critical][regression][slow]") {
    for (const std::string& fluid : {std::string("Water"), std::string("CarbonDioxide"), std::string("Helium")}) {
        auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        const double Tc = heos->T_critical();
        const double rho_c = heos->rhomass_critical();
        int single_phase_points = 0;
        for (double tf : {0.99, 0.999, 0.9999, 1.0, 1.0001, 1.001, 1.01}) {
            const double T = tf * Tc;
            for (double rf : {0.1, 0.3, 0.6, 0.9, 1.0, 1.1, 1.5, 2.0}) {
                const double rho = rf * rho_c;
                double p_ref = 0.0;
                try {
                    heos->update(CoolProp::DmassT_INPUTS, rho, T);
                    if (heos->Q() > 0.0 && heos->Q() < 1.0) continue;  // dome: lever-rule routed
                    p_ref = heos->p();
                } catch (...) {
                    continue;  // outside HEOS validity
                }
                INFO(fluid << " T=" << T << " (" << tf << "Tc)  rho=" << rho << " (" << rf << "rho_c)  p_ref=" << p_ref);
                // (1) no gap: the surrogate MUST serve this single-phase point.
                REQUIRE_NOTHROW(svd->update(CoolProp::DmassT_INPUTS, rho, T));
                const double p_svd = svd->p();
                // (2) accuracy: table holds <1e-4 across the band, incl. the cusp.
                CHECK(std::isfinite(p_svd));
                CHECK(p_svd == Approx(p_ref).epsilon(1e-4));
                CHECK(svd->rhomass() == rho);  // density input echoed back
                ++single_phase_points;
            }
        }
        INFO(fluid << " single-phase points checked: " << single_phase_points);
        REQUIRE(single_phase_points > 20);  // guard against the sweep silently skipping everything
    }
}

// CoolProp-4z79 dome-routing guard: the isobar-bounded NC_SUPER region
// dips below Tc, so its AABB spans the two-phase dome.  A two-phase DT
// query in that sub-Tc sliver must still be routed to the lever rule
// (p = P_sat(T), Q from the lever rule), NOT read as single-phase off the
// surface's interpolated-across-the-dome cells.  This is the exact hole
// the resolve-order dome PRE-CHECK closes; the [near_critical] sweep above
// skips dome points and can't see it, so it is gated explicitly here.
TEST_CASE("SVDSBTL DT two-phase routing holds inside the NC sub-Tc sliver (CoolProp-4z79)",
          "[SBTL][SVDSBTL][dt][near_critical][dome][regression][slow]") {
    for (const std::string& fluid : {std::string("Water"), std::string("CarbonDioxide")}) {
        auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        const double Tc = heos->T_critical();
        const double rho_c = heos->rhomass_critical();
        // Just below Tc but inside the NC_SUPER span (Tc*(1-1e-6) .. Tc),
        // at rho_c — squarely inside the (vanishing) dome there.
        for (double tf : {1.0 - 5.0e-7, 1.0 - 1.0e-7}) {
            const double T = tf * Tc;
            heos->update(CoolProp::DmassT_INPUTS, rho_c, T);
            if (!(heos->Q() > 0.0 && heos->Q() < 1.0)) continue;  // require it actually be two-phase
            const double p_sat = heos->p();
            const double Q_ref = heos->Q();
            INFO(fluid << " T=" << T << " (" << tf << "Tc) rho=rho_c=" << rho_c << " expect two-phase p_sat=" << p_sat);
            REQUIRE_NOTHROW(svd->update(CoolProp::DmassT_INPUTS, rho_c, T));
            CHECK(svd->Q() > 0.0);
            CHECK(svd->Q() < 1.0);
            CHECK(svd->p() == Approx(p_sat).epsilon(1e-3));  // P_sat(T), not an interpolated surface value
            CHECK(svd->Q() == Approx(Q_ref).epsilon(1e-2));
        }
    }
}

// Regression gate for CoolProp-jh6a: the DmassT parent SUPER region uses a
// LOG primary-T axis so that wide-supercritical fluids stay accurate.
// Helium is the witness: Tc=5.2 K, Tmax=2000 K (~380x), where a LINEAR
// axis gave ~1% (p90 ~4%) pressure error across the whole SUPER region.
// LOG holds ~1e-7.  This test sweeps the supercritical region well away
// from the critical point and the dome and asserts <1e-3 -- a tolerance
// that PASSES under LOG (~1e-7) and FAILS under the old LINEAR axis (~1e-2).
TEST_CASE("SVDSBTL DT supercritical accuracy holds for wide-T-range fluids (CoolProp-jh6a)", "[SBTL][SVDSBTL][dt][supercritical][regression][slow]") {
    for (const std::string& fluid : {std::string("Helium"), std::string("Hydrogen")}) {
        auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", fluid));
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluid));
        const double Tc = heos->T_critical();
        const double rho_c = heos->rhomass_critical();
        int checked = 0;
        for (double tf : {1.05, 1.15, 1.3}) {
            const double T = tf * Tc;
            for (double rf : {0.3, 0.6, 1.0, 1.5, 2.0}) {
                const double rho = rf * rho_c;
                double p_ref = 0.0;
                try {
                    heos->update(CoolProp::DmassT_INPUTS, rho, T);
                    p_ref = heos->p();
                } catch (...) {
                    continue;  // outside HEOS validity (e.g. above the melting line)
                }
                if (!(p_ref > 0)) continue;
                INFO(fluid << " T=" << T << " (" << tf << "Tc) rho=" << rho << " (" << rf << "rho_c) p_ref=" << p_ref);
                REQUIRE_NOTHROW(svd->update(CoolProp::DmassT_INPUTS, rho, T));
                CHECK(svd->p() == Approx(p_ref).epsilon(1e-3));
                ++checked;
            }
        }
        INFO(fluid << " supercritical points checked: " << checked);
        REQUIRE(checked > 8);
    }
}

// Tight tolerance comes in the multi-fluid sweep below — this one just
// proves the wiring is end-to-end.
TEST_CASE("SVDSBTL DT-indexed surface builds and evaluates for CO2", "[SBTL][SVDSBTL][dt][smoke][slow]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "CarbonDioxide"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
    // Supercritical CO2 at moderate density — the regime that broke
    // BICUBIC #1301 with rectangular-band errors.
    // Probes inside the atlas coverage envelope for CO2.  SUPER is the
    // primary fix target (yufang67 #1301: low-D supercritical returned
    // -760 MPa under BICUBIC); LIQUID + VAPOR probes pick D above /
    // below the sat dome respectively to land in the single-phase
    // sub-regions the preset emits.
    struct Probe
    {
        double D, T;
        const char* region;
    };
    for (const auto& p : std::vector<Probe>{
           {10.0, 400.0, "SUPER mid-D"},
           {500.0, 350.0, "SUPER high-D"},
           {1100.0, 250.0, "LIQUID (D > rho_sat,L)"},
           {15.0, 250.0, "VAPOR (D < rho_sat,V)"},
         }) {
        heos->update(CoolProp::DmassT_INPUTS, p.D, p.T);
        svd->update(CoolProp::DmassT_INPUTS, p.D, p.T);
        INFO(p.region << " D=" << p.D << " T=" << p.T << " svd_p=" << svd->p() << " heos_p=" << heos->p() << " svd_h=" << svd->hmass()
                      << " heos_h=" << heos->hmass());
        CHECK(std::isfinite(svd->p()));
        CHECK(svd->p() > 0.0);
        CHECK(svd->p() == Approx(heos->p()).epsilon(1e-2));
        CHECK(svd->hmass() == Approx(heos->hmass()).epsilon(1e-2));
        CHECK(svd->rhomass() == p.D);
        CHECK(svd->T() == p.T);
    }
}

// Regression gate for CoolProp-wvtz: the low-density (near-ideal-gas)
// edge of the DT preset's VAPOR/SUPER regions.  The secondary (density)
// axis is normalised η = (ρ − ρ_lo)/(ρ_hi − ρ_lo); with a LINEAR scale
// and ρ_hi/ρ_lo ~ 1e5 the entire ideal-gas tail (where p ∝ ρ sweeps
// several decades) collapses below the first sampled grid node, so the
// EXP-SVD extrapolates and p(ρ, T) was off by 40 %–2400 % there
// (SVDSBTLValidation DT panels showed a high-error band at low ρ).
// A LOG secondary axis (this fix) gives the tail real resolution and
// makes log p ≈ log ρ + const linear in η, so the EXP-SVD is near-exact.
// Probe densities are derived from a low p_target via HEOS PT so they
// land in genuine single-phase low-ρ territory, mirroring the validation
// methodology.
TEST_CASE("SVDSBTL DT low-density pressure matches HEOS (CoolProp-wvtz)", "[SBTL][SVDSBTL][dt][low_density][regression][slow]") {
    struct Case
    {
        const char* fluid;
        double T_factor;  // probe temperature as a multiple of Tc
        const char* region;
    };
    // Reliably-broken pre-fix cases (rev-15 linear-η): Water's wide SUPER/
    // VAPOR regions show 20×–2000× error in the first cell; Argon's SUPER
    // shows ~17 %.  Each probe sits just above the region's low-density
    // floor ρ_lo — exactly where the validation heatmap's left edge lives.
    for (const Case& c : {Case{"Water", 1.05, "SUPER"}, Case{"Water", 0.85, "VAPOR"}, Case{"Argon", 1.05, "SUPER"}}) {
        auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", c.fluid));
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", c.fluid));
        const double T = c.T_factor * heos->T_critical();

        // Replicate the preset's low-density floor: p_min_eos ≈ p_sat at
        // the triple point.  ρ_lo = ρ(T, p_floor) is the region's bottom
        // boundary; probe a hair above it (1.3·ρ_lo) so the point is
        // unambiguously inside the atlas yet inside the first η cell, the
        // near-ideal-gas tail where the linear-η surface failed.
        heos->update(CoolProp::QT_INPUTS, 0.0, heos->Ttriple() * 1.001);
        const double p_floor = heos->p() * 1.02;
        heos->update(CoolProp::PT_INPUTS, p_floor, T);
        const double rho = 1.3 * heos->rhomass();

        heos->update(CoolProp::DmassT_INPUTS, rho, T);
        const double p_ref = heos->p();
        svd->update(CoolProp::DmassT_INPUTS, rho, T);
        const double p_svd = svd->p();
        INFO(c.fluid << " " << c.region << " T=" << T << " ρ=" << rho << " p_ref=" << p_ref << " p_svd=" << p_svd
                     << " rel_err=" << std::abs(p_svd - p_ref) / p_ref);
        CHECK(std::isfinite(p_svd));
        CHECK(p_svd > 0.0);
        CHECK(p_svd == Approx(p_ref).epsilon(5e-3));
        CHECK(svd->rhomass() == rho);  // density input echoed back
    }
}

// Regression gate for issue #1301 (CoolProp-i7j): a random (T, ρ)
// sweep over yufang67's exact bounds (T ∈ [220, 500] K, ρ ∈ [0.01,
// 1200] kg/m³) spanning subcritical liquid / vapor / two-phase dome
// AND supercritical.  Asserts the headline contrast: SVDSBTL's
// DT-indexed P(ρ, T) is clean across the envelope it covers, while
// BICUBIC reproduces the documented failures (large near-saturation
// error + negative-pressure cells).  The browsable figure is rendered
// separately by Web/coolprop/_gen/gen_DT_validation_fig1301.py (pure
// Python; same sweep).
// Two-phase dome routing: a DT probe with rho_sat,V(T) < ρ < rho_sat,L(T)
// at subcritical T lands in the dome.  Pressure must be P_sat(T) and the
// caloric props must come from the lever rule — in particular hmass() /
// hmolar() must NOT throw (regression: the dome block originally left
// pt.h_mass unset and the DomeBlend switch has no enthalpy arm, so the
// enthalpy query hit `default` and threw ValueError).
TEST_CASE("SVDSBTL DT two-phase dome returns lever-rule props (CoolProp-i7j)", "[SBTL][SVDSBTL][dt][dome]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "CarbonDioxide"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
    // T=250 K subcritical; ρ=500 sits between ρ_sat,V (~18) and
    // ρ_sat,L (~1046), squarely two-phase.
    const double D = 500.0;
    const double T = 250.0;
    heos->update(CoolProp::DmassT_INPUTS, D, T);
    svd->update(CoolProp::DmassT_INPUTS, D, T);
    INFO("dome probe D=" << D << " T=" << T << " svd_p=" << svd->p() << " heos_p=" << heos->p());
    REQUIRE(svd->Q() > 0.0);
    REQUIRE(svd->Q() < 1.0);
    REQUIRE(svd->p() == Approx(heos->p()).epsilon(1e-3));  // P_sat(T)
    // hmass() / hmolar() / umass() / umolar() must return the lever-rule
    // blend, not throw.  For a DmassT dome point these now flow through the
    // LAZY DomeBlend enthalpy arms (ensure_dome_h_endpoints_, CoolProp-qp0n)
    // — the resolve no longer eagerly fills h — so this both checks the
    // values and exercises the lazy path for h and u.
    REQUIRE_NOTHROW(svd->hmass());
    REQUIRE(svd->hmass() == Approx(heos->hmass()).epsilon(1e-3));
    REQUIRE(svd->hmolar() == Approx(heos->hmolar()).epsilon(1e-3));
    REQUIRE_NOTHROW(svd->umass());
    REQUIRE(svd->umass() == Approx(heos->umass()).epsilon(1e-3));
    REQUIRE(svd->umolar() == Approx(heos->umolar()).epsilon(1e-3));
    REQUIRE(svd->rhomass() == D);  // density input echoed back
}

// fast_evaluate on a DT surface requesting density back: density is the
// INPUT, not a tabulated property, so the batched output plan must
// short-circuit iDmass / iDmolar to the input echo.  Regression: it
// originally routed them to a surface lookup, property_index(iDmass)
// threw, and the whole output row was NaN'd with internal_error.
TEST_CASE("SVDSBTL DT fast_evaluate echoes density + matches HEOS (CoolProp-i7j)", "[SBTL][SVDSBTL][dt][fast_evaluate]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "CarbonDioxide"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
    // Supercritical single-phase batch.
    const std::vector<double> Dvals = {10.0, 100.0, 500.0, 800.0};  // DmassT: val1 = D
    const double T = 350.0;
    std::vector<double> v2(Dvals.size(), T);  // val2 = T
    std::vector<CoolProp::parameters> outs = {CoolProp::iP, CoolProp::iDmass, CoolProp::iHmass};
    std::vector<double> out(Dvals.size() * outs.size(), std::nan(""));
    std::vector<int> status(Dvals.size(), CoolProp::fast_evaluate_internal_error);
    svd->fast_evaluate(CoolProp::DmassT_INPUTS, Dvals.data(), v2.data(), Dvals.size(), outs.data(), outs.size(), out.data(), out.size(),
                       status.data(), status.size());
    for (std::size_t k = 0; k < Dvals.size(); ++k) {
        INFO("fast_evaluate k=" << k << " D=" << Dvals[k] << " status=" << static_cast<int>(status[k]));
        REQUIRE(status[k] == CoolProp::fast_evaluate_ok);
        const double p_fe = out[k * outs.size() + 0];
        const double D_fe = out[k * outs.size() + 1];
        const double h_fe = out[k * outs.size() + 2];
        heos->update(CoolProp::DmassT_INPUTS, Dvals[k], T);
        REQUIRE(D_fe == Dvals[k]);  // density echoed, not NaN
        REQUIRE(p_fe == Approx(heos->p()).epsilon(1e-2));
        REQUIRE(h_fe == Approx(heos->hmass()).epsilon(1e-2));
    }
}

TEST_CASE("CoolProp-i7j: SVDSBTL DT clean vs BICUBIC fails on #1301 sweep", "[SBTL][SVDSBTL][dt][fig1301][slow]") {
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "CarbonDioxide"));
    auto bic = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("BICUBIC&HEOS", "CarbonDioxide"));

    constexpr double T_lo = 220.0;
    constexpr double T_hi = 500.0;
    constexpr double rho_lo = 1.0e-2;
    constexpr double rho_hi = 1200.0;
    constexpr int N_SAMPLES = 20000;  // ample to characterize the contrast in CI
    std::mt19937 rng(1301);
    std::uniform_real_distribution<double> T_dist(T_lo, T_hi);
    std::uniform_real_distribution<double> rho_dist(rho_lo, rho_hi);

    int n_heos_valid = 0;   // cells where HEOS gives a valid positive reference p
    int n_bicubic_bad = 0;  // BICUBIC cells with |%-err| > 1% or sign-flip, inside HEOS validity
    int n_svd_bad = 0;      // SVDSBTL cells with |%-err| > 1% OR a finite non-positive p
    int n_svd_covered = 0;  // cells SVDSBTL resolved to a FINITE p (positive or not)
    auto pressure = [](CoolProp::AbstractState& s, double rho, double T) -> double {
        try {
            s.update(CoolProp::DmassT_INPUTS, rho, T);
            return s.p();
        } catch (...) {  // NOLINT(bugprone-empty-catch)
            return std::nan("");
        }
    };
    for (int k = 0; k < N_SAMPLES; ++k) {
        const double T = T_dist(rng);
        const double rho = rho_dist(rng);
        const double p_heos = pressure(*heos, rho, T);
        if (!std::isfinite(p_heos) || p_heos == 0.0) continue;
        ++n_heos_valid;
        const double p_svd = pressure(*svd, rho, T);
        const double p_bic = pressure(*bic, rho, T);
        // %-error as yufang67 defined it: |P_bk - P_EOS|/P_EOS*100.
        if (std::isfinite(p_bic) && (std::abs(p_bic - p_heos) / std::abs(p_heos) * 100.0 > 1.0 || p_bic <= 0.0)) {
            ++n_bicubic_bad;
        }
        if (std::isfinite(p_svd)) {
            ++n_svd_covered;
            // A finite non-positive pressure is itself a failure (the
            // #1301 sign-flip mode), not just a >1% deviation — count
            // both so a regression that returns negative P can't hide
            // by being dropped from the denominator.
            if (p_svd <= 0.0 || std::abs(p_svd - p_heos) / std::abs(p_heos) * 100.0 > 1.0) ++n_svd_bad;
        }
    }
    UNSCOPED_INFO("n_heos_valid=" << n_heos_valid << " n_svd_covered=" << n_svd_covered << " n_svd_bad=" << n_svd_bad
                                  << " n_bicubic_bad=" << n_bicubic_bad);
    // SVDSBTL must cover most of the HEOS-valid envelope (the only
    // uncovered band is the deep low-ρ ideal-gas corner below the
    // p_triple sampling floor — a documented follow-up), and within
    // what it covers the failure fraction (>1% error OR non-positive
    // p) must be tiny.  The handful that remain are at the extreme
    // low-ρ atlas edge where the η-normalisation sits at ξ≈0 — NOT the
    // near-saturation region that motivated #1301.
    CHECK(n_heos_valid > 0);
    CHECK(static_cast<double>(n_svd_covered) > 0.90 * static_cast<double>(n_heos_valid));
    CHECK(static_cast<double>(n_svd_bad) < 0.001 * static_cast<double>(n_svd_covered));
    CHECK(n_bicubic_bad > 0);  // #1301 failure pattern must reproduce
}

// Anomaly detection: water has a single rho_sat,L extremum at
// T ≈ 277 K (4 °C density anomaly).  CO2 / methane have none.
// Verifies find_rho_satL_extrema_T does its job before the preset
// uses it to split the LIQUID region.
TEST_CASE("find_rho_satL_extrema_T returns water anomaly, none for other fluids", "[SBTL][SVDSBTL][dt][anomaly]") {
    SECTION("water has anomaly near 277 K") {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
        auto extrema = CoolProp::sbtl::find_rho_satL_extrema_T(*heos, heos->Ttriple() * 1.001, heos->T_critical() * 0.999);
        REQUIRE(extrema.size() == 1);
        REQUIRE(extrema[0] == Approx(277.13).margin(2.0));
    }
    SECTION("CO2 has no anomaly") {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
        auto extrema = CoolProp::sbtl::find_rho_satL_extrema_T(*heos, heos->Ttriple() * 1.001, heos->T_critical() * 0.999);
        REQUIRE(extrema.empty());
    }
    SECTION("methane has no anomaly") {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Methane"));
        auto extrema = CoolProp::sbtl::find_rho_satL_extrema_T(*heos, heos->Ttriple() * 1.001, heos->T_critical() * 0.999);
        REQUIRE(extrema.empty());
    }
}

// LIQUID-split geometry: the dt_subcritical preset must emit one extra
// LIQUID sub-region per rho_sat,L(T) extremum, plus the three fixed NC
// near-critical regions (NC_LIQUID + NC_VAPOR + NC_SUPER, CoolProp-4z79).
// Water (one anomaly) gets LIQUID_LO + LIQUID_HI + VAPOR + SUPER + 3 NC = 7
// regions; CO2 (no anomaly) gets LIQUID + VAPOR + SUPER + 3 NC = 6.  This
// is the structural guard that the anomaly split actually happens —
// without it, a single LIQUID region would straddle the density-maximum
// hairpin and the SVD couldn't represent the resulting discontinuity.
TEST_CASE("dt_subcritical splits LIQUID at the water density anomaly", "[SBTL][SVDSBTL][dt][anomaly]") {
    SECTION("water emits 7 regions (LIQUID split LO + HI, + 3 NC)") {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
        auto spec = CoolProp::sbtl::presets::dt_subcritical(*heos);
        REQUIRE(spec.regions.size() == 7);
        REQUIRE(spec.input_pair == CoolProp::DmassT_INPUTS);
    }
    SECTION("CO2 emits 6 regions (single LIQUID, + 3 NC)") {
        auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "CarbonDioxide"));
        auto spec = CoolProp::sbtl::presets::dt_subcritical(*heos);
        REQUIRE(spec.regions.size() == 6);
    }
}

// End-to-end: water DT lookups across the 4 °C density anomaly.  Probe
// compressed-liquid states at D just above the saturation dome, sweeping
// T from below the anomaly (LIQUID_LO sub-region) through the maximum to
// above it (LIQUID_HI).  D ≈ 1000-1050 kg/m³ stays just above
// rho_sat,L(T) (≈ 999.9 at the maximum) across the whole sweep, so every
// probe is single-phase liquid sitting right on top of the anomaly's
// hairpin — exactly the cells a single un-split LIQUID region would
// botch.  All must match HEOS density-input pressure / enthalpy.
TEST_CASE("SVDSBTL DT water across the density anomaly matches HEOS", "[SBTL][SVDSBTL][dt][anomaly][water][slow]") {
    auto svd = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL&HEOS", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));

    // T sweep brackets the anomaly at ~277.13 K: two points below, the
    // maximum, two above.  D values are compressed-liquid, comfortably
    // above rho_sat,L (which peaks at ~999.97 kg/m³) yet inside the
    // HEOS validity envelope.
    for (const double T : {274.0, 276.0, 277.13, 279.0, 283.0}) {
        for (const double D : {1000.5, 1020.0, 1050.0}) {
            heos->update(CoolProp::DmassT_INPUTS, D, T);
            const double p_ref = heos->p();
            const double h_ref = heos->hmass();
            svd->update(CoolProp::DmassT_INPUTS, D, T);
            INFO("water anomaly probe D=" << D << " T=" << T << " svd_p=" << svd->p() << " heos_p=" << p_ref << " svd_h=" << svd->hmass()
                                          << " heos_h=" << h_ref);
            REQUIRE_FALSE(std::isnan(svd->p()));
            REQUIRE(svd->p() > 0.0);
            // 0.1% on p across the anomaly — the split keeps each
            // sub-region's rho_sat,L(T) boundary monotone, so the
            // η-normalisation stays well-conditioned right at the
            // density maximum.
            REQUIRE(svd->p() == Approx(p_ref).epsilon(1e-3));
            REQUIRE(svd->hmass() == Approx(h_ref).epsilon(1e-3));
            REQUIRE(svd->rhomass() == D);
            REQUIRE(svd->T() == T);
        }
    }
}

#endif  // ENABLE_CATCH
