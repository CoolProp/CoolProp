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
bool cache_present_for(const std::string& fluid) {
    namespace cp_sbtl = CoolProp::sbtl;
    return std::filesystem::exists(cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid, CoolProp::HmassP_INPUTS))
           && std::filesystem::exists(cp_sbtl::SVDSurfaceSerializer::default_cache_path(fluid, CoolProp::PT_INPUTS));
}

}  // namespace

TEST_CASE("SVDSBTL backend constructs via AbstractState::factory", "[SVDSBTL][factory][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
    REQUIRE(AS != nullptr);
    REQUIRE(AS->backend_name() == "SVDSBTLBackend");
    REQUIRE(AS->fluid_names().size() == 1);
    REQUIRE(AS->fluid_names()[0] == "Water");
    REQUIRE(cache_present_for("Water"));
}

TEST_CASE("SVDSBTL backend reports critical / triple constants from HEOS reference", "[SVDSBTL][constants][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    REQUIRE(AS->T_critical() == Approx(heos->T_critical()).epsilon(1e-12));
    REQUIRE(AS->p_critical() == Approx(heos->p_critical()).epsilon(1e-12));
    REQUIRE(AS->Ttriple() == Approx(heos->Ttriple()).epsilon(1e-12));
    REQUIRE(AS->p_triple() == Approx(heos->p_triple()).epsilon(1e-12));
    REQUIRE(AS->molar_mass() == Approx(heos->molar_mass()).epsilon(1e-12));
    REQUIRE(AS->Tmax() == Approx(heos->Tmax()).epsilon(1e-12));
}

TEST_CASE("SVDSBTL backend PT lookup matches HEOS within tolerance", "[SVDSBTL][pt][water][slow]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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
      std::unique_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", std::vector<std::string>{"Water", "Ethanol"})),
      ContainsSubstring("pure-fluid"));
}

TEST_CASE("SVDSBTL backend supports the high-level PropsSI surface", "[SVDSBTL][propsi][slow]") {
    auto heos = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water"));
    heos->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_truth = heos->rhomass();

    const double rho_via_propsi = CoolProp::PropsSI("D", "T", 350.0, "P", 1.0e6, "SVDSBTL::Water");
    INFO("PropsSI = " << rho_via_propsi << "  HEOS = " << rho_truth);
    REQUIRE(rho_via_propsi == Approx(rho_truth).epsilon(5e-3));
}

TEST_CASE("SVDSBTL backend cache reload produces the same result", "[SVDSBTL][cache][slow]") {
    // First instance populates the cache (or no-ops if already present).
    auto first = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
    REQUIRE(cache_present_for("Water"));
    first->update(CoolProp::PT_INPUTS, 1.0e6, 350.0);
    const double rho_first = first->rhomass();

    // Second instance hits the cache.
    auto second = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
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

TEST_CASE("SVDSBTL backend Q out of [0, 1] is rejected", "[SVDSBTL][twophase][reject]") {
    auto AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("SVDSBTL", "Water"));
    REQUIRE_THROWS(AS->update(CoolProp::PQ_INPUTS, 1.0e6, -0.1));
    REQUIRE_THROWS(AS->update(CoolProp::PQ_INPUTS, 1.0e6, 1.1));
    REQUIRE_THROWS(AS->update(CoolProp::QT_INPUTS, -0.1, 350.0));
    REQUIRE_THROWS(AS->update(CoolProp::QT_INPUTS, 1.1, 350.0));
}

#endif  // ENABLE_CATCH
