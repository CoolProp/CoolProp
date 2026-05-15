// Catch2 tests for the SBTL adapter layer (Phase 2b): SatBoundaryFactory,
// SurfaceSpec, SVDSurfaceFactory, SVDSurface, presets.
//
// The serializer round-trip tests are added in a separate file once
// the SVDSurfaceSerializer lands.  For now this exercises the
// build-and-eval end-to-end at a SMALL grid so the suite stays fast.

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/ConstantCurve.h"
#include "CoolProp/sbtl/SVDSurface.h"
#include "CoolProp/sbtl/SVDSurfaceFactory.h"
#include "CoolProp/sbtl/SatBoundaryFactory.h"
#include "CoolProp/sbtl/SurfaceSpec.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <memory>
#    include <random>
#    include <vector>

#    include "AbstractState.h"

namespace cp_sbtl = CoolProp::sbtl;
namespace cp_region = CoolProp::region;

TEST_CASE("SVDSurface construction + seal", "[SBTL][SVDSurface]") {
    cp_sbtl::SVDSurface s("Water", ::CoolProp::PT_INPUTS, {::CoolProp::iDmass, ::CoolProp::iHmass});
    REQUIRE(s.fluid_name() == "Water");
    REQUIRE(s.input_pair() == ::CoolProp::PT_INPUTS);
    REQUIRE(s.properties().size() == 2);
    REQUIRE(s.contains_property(::CoolProp::iDmass));
    REQUIRE_FALSE(s.contains_property(::CoolProp::iT));
    REQUIRE_FALSE(s.sealed());

    // seal() with no regions registered is a no-op success (empty atlas).
    s.seal();
    REQUIRE(s.sealed());
    REQUIRE(std::isnan(s.eval(::CoolProp::iDmass, 1e5, 300.0)));  // no regions → NaN
}

TEST_CASE("SVDSurface rejects duplicate properties", "[SBTL][SVDSurface]") {
    REQUIRE_THROWS_AS(cp_sbtl::SVDSurface("X", ::CoolProp::PT_INPUTS, {::CoolProp::iDmass, ::CoolProp::iDmass}), std::invalid_argument);
}

TEST_CASE("SVDSurface rejects empty property list", "[SBTL][SVDSurface]") {
    REQUIRE_THROWS_AS(cp_sbtl::SVDSurface("X", ::CoolProp::PT_INPUTS, {}), std::invalid_argument);
}

TEST_CASE("SatBoundaryFactory builds water sat curves matching HEOS", "[SBTL][SatBoundaryFactory][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    const auto [p_min, p_max] = cp_sbtl::subcritical_pressure_range(*heos);

    auto h_sat_L = cp_sbtl::build_h_sat_L(*heos, p_min, p_max);
    auto h_sat_V = cp_sbtl::build_h_sat_V(*heos, p_min, p_max);

    // Sample at 10 random log-p points, compare to direct HEOS.
    std::mt19937 rng(31);
    std::uniform_real_distribution<double> u(std::log(p_min * 1.05), std::log(p_max * 0.99));
    for (int k = 0; k < 10; ++k) {
        const double p = std::exp(u(rng));
        heos->update(::CoolProp::PQ_INPUTS, p, 0.0);
        const double hL_truth = heos->hmass();
        heos->update(::CoolProp::PQ_INPUTS, p, 1.0);
        const double hV_truth = heos->hmass();
        // CubicSpline through 64 HEOS knots should match HEOS within
        // ~1e-4 relative on the smooth water sat curve.
        REQUIRE(h_sat_L->eval(p) == Catch::Approx(hL_truth).epsilon(1e-3));
        REQUIRE(h_sat_V->eval(p) == Catch::Approx(hV_truth).epsilon(1e-3));
    }
}

TEST_CASE("SVDSurface PH preset builds + evals against HEOS", "[SBTL][SVDSurface][preset_ph][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    // Tiny grid for unit-test runtime — accuracy is checked
    // qualitatively (must reproduce HEOS to a few percent on a
    // handful of single-phase probes).  Phase 2a's e2e tool covers
    // production-resolution accuracy.
    auto spec = cp_sbtl::presets::ph_subcritical(*heos, /*NT=*/40, /*NR=*/80, /*rank=*/10);
    REQUIRE(spec.input_pair == ::CoolProp::HmassP_INPUTS);
    REQUIRE(spec.regions.size() == 2);
    REQUIRE(spec.properties.size() == 4);

    auto surface = cp_sbtl::build_surface(*heos, std::move(spec));
    REQUIRE(surface.sealed());
    REQUIRE(surface.region_count() == 2);

    // Probe single-phase Water states in (T, p) → look up via (h, p).
    std::mt19937 rng(57);
    const auto [p_min, p_max] = cp_sbtl::subcritical_pressure_range(*heos);
    std::uniform_real_distribution<double> uT(290.0, heos->Tmax() - 1.0);
    std::uniform_real_distribution<double> u_log_p(std::log(p_min * 1.5), std::log(p_max * 0.9));

    int kept = 0;
    double max_rel = 0;
    for (int trial = 0; trial < 200 && kept < 30; ++trial) {
        const double T = uT(rng);
        const double p = std::exp(u_log_p(rng));
        try {
            heos->update(::CoolProp::PT_INPUTS, p, T);
            if (heos->Q() > 0.0 && heos->Q() < 1.0) {
                continue;
            }
            const double h = heos->hmass();
            const double rho_truth = heos->rhomass();
            const double rho_pred = surface.eval(::CoolProp::iDmass, p, h);
            if (std::isnan(rho_pred)) {
                continue;
            }
            const double rel = std::abs(rho_pred - rho_truth) / rho_truth;
            max_rel = std::max(max_rel, rel);
            ++kept;
        } catch (...) {
            // skip
        }
    }
    INFO("kept=" << kept << " max_rel=" << max_rel);
    REQUIRE(kept >= 10);
    // 40x80 rank-10 SVD over the whole subcritical range is coarse —
    // anything under 5% on Water is a sanity-passing build.
    REQUIRE(max_rel < 5e-2);
}

#endif  // ENABLE_CATCH
