// Catch2 tests for the SBTL adapter layer (Phase 2b): SatBoundaryFactory,
// SurfaceSpec, SVDSurfaceFactory, SVDSurface, presets.
//
// The serializer round-trip tests are added in a separate file once
// the SVDSurfaceSerializer lands.  For now this exercises the
// build-and-eval end-to-end at a SMALL grid so the suite stays fast.

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <memory>
#    include <random>
#    include <vector>

#    include "AbstractState.h"
#    include "CoolProp/region/AxisTransform.h"
#    include "CoolProp/region/ConstantCurve.h"
#    include "CoolProp/sbtl/SVDSurface.h"
#    include "CoolProp/sbtl/SVDSurfaceFactory.h"
#    include "CoolProp/sbtl/SVDSurfaceSerializer.h"
#    include "CoolProp/sbtl/SatBoundaryFactory.h"
#    include "CoolProp/sbtl/SurfaceSpec.h"
#    include "miniz.h"

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

TEST_CASE("SVDSurface save/load round-trip is bit-identical", "[SBTL][serializer][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    auto spec = cp_sbtl::presets::ph_subcritical(*heos, /*NT=*/40, /*NR=*/80, /*rank=*/10);
    auto surface_before = cp_sbtl::build_surface(*heos, std::move(spec));

    // Save to a byte buffer and load back.
    const auto blob = cp_sbtl::SVDSurfaceSerializer::save(surface_before);
    auto surface_after = cp_sbtl::SVDSurfaceSerializer::load(blob);

    // Round-trip checks: metadata identity.
    REQUIRE(surface_after.fluid_name() == surface_before.fluid_name());
    REQUIRE(surface_after.input_pair() == surface_before.input_pair());
    REQUIRE(surface_after.region_count() == surface_before.region_count());
    REQUIRE(surface_after.properties() == surface_before.properties());

    // Same query at multiple probes must give bit-identical results.
    // (Floating-point bit-identical because Region/SVDEvaluator math
    // is deterministic and we packed the exact ULP-precise doubles.)
    std::mt19937 rng(101);
    const auto [p_min, p_max] = cp_sbtl::subcritical_pressure_range(*heos);
    std::uniform_real_distribution<double> u_log_p(std::log(p_min * 1.5), std::log(p_max * 0.9));
    std::uniform_real_distribution<double> uT(290.0, heos->Tmax() - 1.0);

    int compared = 0;
    for (int trial = 0; trial < 200 && compared < 30; ++trial) {
        const double T = uT(rng);
        const double p = std::exp(u_log_p(rng));
        try {
            heos->update(::CoolProp::PT_INPUTS, p, T);
            if (heos->Q() > 0.0 && heos->Q() < 1.0) {
                continue;
            }
            const double h = heos->hmass();
            const double rho_before = surface_before.eval(::CoolProp::iDmass, p, h);
            const double rho_after = surface_after.eval(::CoolProp::iDmass, p, h);
            if (std::isnan(rho_before) || std::isnan(rho_after)) {
                continue;
            }
            // Same bytes → same evaluator → bit-identical output.
            REQUIRE(rho_before == rho_after);
            ++compared;
        } catch (...) {  // NOLINT(bugprone-empty-catch)
        }
    }
    REQUIRE(compared >= 10);
}

TEST_CASE("SVDSurfaceSerializer rejects corrupt input", "[SBTL][serializer]") {
    // Empty / truncated / wrong-magic bytes should all throw, not segfault.
    const std::vector<char> empty;
    REQUIRE_THROWS_AS(cp_sbtl::SVDSurfaceSerializer::load(empty), std::runtime_error);

    const std::vector<char> garbage(64, '\0');
    REQUIRE_THROWS_AS(cp_sbtl::SVDSurfaceSerializer::load(garbage), std::runtime_error);

    // Valid zlib-compressed payload that decompresses to "hello world" —
    // not a valid msgpack stream, so msgpack::unpack throws.  The
    // serializer wraps that in std::runtime_error.
    std::vector<char> bogus_payload;
    {
        const char hello[] = "hello world";
        bogus_payload.resize(64);
        mz_ulong out_size = static_cast<mz_ulong>(bogus_payload.size());
        compress((unsigned char*)bogus_payload.data(), &out_size, (const unsigned char*)hello, sizeof(hello) - 1);
        bogus_payload.resize(out_size);
    }
    REQUIRE_THROWS(cp_sbtl::SVDSurfaceSerializer::load(bogus_payload));
}

// Multi-fluid coverage: build and round-trip PH + PT surfaces for a
// handful of fluids at a tiny grid.  Tiny grid (NT=30, NR=50, rank=8)
// keeps total per-fluid runtime ≈ 1–2 s; the test set runs in
// ~10 s with all 5 fluids.  Loose accuracy tolerance — production
// resolution is what Phase 2a's SVDSBTL_E2E covers.
TEST_CASE("SVDSurface PH preset across multi-fluid set", "[SBTL][SVDSurface][preset_ph][multi_fluid][slow]") {
    for (const auto* fluid : {"R134a", "Ammonia", "Methane", "Propane", "CarbonDioxide"}) {
        SECTION(fluid) {
            auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", fluid));
            auto spec = cp_sbtl::presets::ph_subcritical(*heos, /*NT=*/30, /*NR=*/50, /*rank=*/8);
            cp_sbtl::SVDSurface surface = cp_sbtl::build_surface(*heos, std::move(spec));
            REQUIRE(surface.sealed());
            REQUIRE(surface.region_count() == 2);

            // Round-trip through serializer.  HEOS canonicalises aliases
            // (e.g. "Propane" → "n-Propane"), so compare against what the
            // surface actually stored, not the user-facing input string.
            const auto blob = cp_sbtl::SVDSurfaceSerializer::save(surface);
            auto reloaded = cp_sbtl::SVDSurfaceSerializer::load(blob);
            REQUIRE(reloaded.fluid_name() == surface.fluid_name());
            REQUIRE(reloaded.region_count() == 2);

            // Spot-check eval matches HEOS at a handful of single-phase probes.
            std::mt19937 rng(91);
            const auto [p_min, p_max] = cp_sbtl::subcritical_pressure_range(*heos);
            std::uniform_real_distribution<double> uT(std::max(heos->Ttriple() + 5.0, heos->Tmin() + 5.0), heos->Tmax() - 5.0);
            std::uniform_real_distribution<double> u_log_p(std::log(p_min * 2.0), std::log(p_max * 0.8));
            int kept = 0;
            double max_rel = 0;
            for (int trial = 0; trial < 200 && kept < 10; ++trial) {
                const double T = uT(rng);
                const double p = std::exp(u_log_p(rng));
                try {
                    heos->update(::CoolProp::PT_INPUTS, p, T);
                    if (heos->Q() > 0.0 && heos->Q() < 1.0) {
                        continue;
                    }
                    const double h = heos->hmass();
                    const double rho_truth = heos->rhomass();
                    const double rho_pred = reloaded.eval(::CoolProp::iDmass, p, h);
                    if (std::isnan(rho_pred)) {
                        continue;
                    }
                    const double rel = std::abs(rho_pred - rho_truth) / rho_truth;
                    max_rel = std::max(max_rel, rel);
                    ++kept;
                } catch (...) {  // NOLINT(bugprone-empty-catch)
                }
            }
            INFO(fluid << " kept=" << kept << " max_rel=" << max_rel);
            REQUIRE(kept >= 3);
            // Tiny grid → loose tolerance; just confirm no order-of-magnitude regressions.
            REQUIRE(max_rel < 0.20);
        }
    }
}

TEST_CASE("SVDSurface PT preset water round-trip", "[SBTL][SVDSurface][preset_pt][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    auto spec = cp_sbtl::presets::pt_subcritical(*heos, /*NT=*/40, /*NR=*/80, /*rank=*/10);
    REQUIRE(spec.input_pair == ::CoolProp::PT_INPUTS);
    REQUIRE(spec.regions.size() == 2);

    cp_sbtl::SVDSurface surface = cp_sbtl::build_surface(*heos, std::move(spec));
    REQUIRE(surface.sealed());

    // Round-trip serializer.
    const auto blob = cp_sbtl::SVDSurfaceSerializer::save(surface);
    auto reloaded = cp_sbtl::SVDSurfaceSerializer::load(blob);
    REQUIRE(reloaded.input_pair() == ::CoolProp::PT_INPUTS);

    // Spot-check on water at a known compressed-liquid state.
    heos->update(::CoolProp::PT_INPUTS, 1e6, 350.0);
    const double rho_truth = heos->rhomass();
    const double rho_pred = reloaded.eval(::CoolProp::iDmass, 1e6, 350.0);
    INFO("rho_truth=" << rho_truth << " rho_pred=" << rho_pred);
    REQUIRE_FALSE(std::isnan(rho_pred));
    REQUIRE(std::abs(rho_pred - rho_truth) / rho_truth < 0.05);
}

TEST_CASE("SVDSurfaceSerializer file save/load round-trip", "[SBTL][serializer][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    auto spec = cp_sbtl::presets::ph_subcritical(*heos, /*NT=*/40, /*NR=*/80, /*rank=*/10);
    auto surface = cp_sbtl::build_surface(*heos, std::move(spec));

    // Save to a tmp file and load back via the file API.
    const std::string path = "/tmp/svd_test_water_ph.svd.bin.z";
    cp_sbtl::SVDSurfaceSerializer::save_to_file(surface, path);
    auto loaded = cp_sbtl::SVDSurfaceSerializer::load_from_file(path);
    REQUIRE(loaded.fluid_name() == "Water");
    REQUIRE(loaded.region_count() == 2);
    // Spot-check a single eval matches.
    heos->update(::CoolProp::PT_INPUTS, 1e5, 350.0);
    const double h = heos->hmass();
    REQUIRE(loaded.eval(::CoolProp::iDmass, 1e5, h) == surface.eval(::CoolProp::iDmass, 1e5, h));
}

#endif  // ENABLE_CATCH
