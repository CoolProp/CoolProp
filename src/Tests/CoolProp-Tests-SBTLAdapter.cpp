// Catch2 tests for the SBTL adapter layer (Phase 2b): SatBoundaryFactory,
// SurfaceSpec, SVDSurfaceFactory, SVDSurface, presets, and the
// SVDSurfaceSerializer round-trip + corruption-rejection coverage.
// Build-and-eval tests use a SMALL grid so the suite stays fast;
// production-resolution accuracy is covered by Phase 2a's SVDSBTL_E2E.

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <filesystem>
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
#    include "TestUtils.h"
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

TEST_CASE("SuperancillaryBoundaryCurve::eval_fast tracks eval within surrogate budget", "[SBTL][SatBoundaryFactory][eval_fast][slow]") {
    // CoolProp-akt: eval_fast is a 1024-point log-spaced linear-interp
    // surrogate for the SA-backed boundary curve, used by Region::
    // curve_contains where sign-only accuracy is sufficient.  Atlas-
    // dispatch correctness depends on eval_fast giving the same
    // bracketing decision as eval at every probe.  Check tightness
    // empirically over Water's whole subcritical p range.
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    const auto [p_min, p_max] = cp_sbtl::subcritical_pressure_range(*heos);

    auto h_sat_L = cp_sbtl::build_h_sat_L(*heos, p_min, p_max);
    auto h_sat_V = cp_sbtl::build_h_sat_V(*heos, p_min, p_max);

    std::mt19937 rng(89);
    std::uniform_real_distribution<double> u(std::log(p_min * 1.05), std::log(p_max * 0.99));
    double max_rel_L = 0.0;
    double max_rel_V = 0.0;
    for (int k = 0; k < 1024; ++k) {
        const double p = std::exp(u(rng));
        const double yL = h_sat_L->eval(p);
        const double yV = h_sat_V->eval(p);
        const double yL_fast = h_sat_L->eval_fast(p);
        const double yV_fast = h_sat_V->eval_fast(p);
        // Fast-path finiteness: eval_fast() interpolates from a
        // surrogate table populated by eval(); if any cell in the
        // table is NaN/Inf the interp can propagate that.  Assert
        // explicitly so std::max doesn't silently absorb a non-finite
        // value and let the 1e-3 budget pass with garbage in
        // max_rel_L/max_rel_V.
        if (std::isfinite(yL) && std::abs(yL) > 0.0) {
            REQUIRE(std::isfinite(yL_fast));
            const double relL = std::abs(yL_fast - yL) / std::abs(yL);
            REQUIRE(std::isfinite(relL));
            max_rel_L = std::max(max_rel_L, relL);
        }
        if (std::isfinite(yV) && std::abs(yV) > 0.0) {
            REQUIRE(std::isfinite(yV_fast));
            const double relV = std::abs(yV_fast - yV) / std::abs(yV);
            REQUIRE(std::isfinite(relV));
            max_rel_V = std::max(max_rel_V, relV);
        }
    }
    INFO("max rel-err: h_sat_L=" << max_rel_L << " h_sat_V=" << max_rel_V);
    // The saturation enthalpies h_sat,L/V have a sqrt-singularity at
    // p_crit (latent heat -> 0 like sqrt(p_crit - p)), so linear
    // interp on a log-uniform grid bottoms out at ~5e-4 near the
    // singularity even at 1024 knots.  Sign-only callers (atlas
    // curve_contains) tolerate this: a probe within 5e-4 relative of
    // the true sat curve is essentially ON the curve, and either
    // routing (single-phase via the matched region's SVD vs dome blend
    // via the lever rule) returns a numerically consistent answer at
    // the boundary.  The 1e-3 budget is the contract the atlas relies
    // on — outside the near-critical strip the empirical error is
    // 1e-7 ish.
    REQUIRE(max_rel_L < 1e-3);
    REQUIRE(max_rel_V < 1e-3);
}

TEST_CASE("SVDSurface PH preset builds + evals against HEOS", "[SBTL][SVDSurface][preset_ph][slow]") {
    auto heos = std::shared_ptr<::CoolProp::AbstractState>(::CoolProp::AbstractState::factory("HEOS", "Water"));
    // Tiny grid for unit-test runtime — accuracy is checked
    // qualitatively (must reproduce HEOS to a few percent on a
    // handful of single-phase probes).  Phase 2a's e2e tool covers
    // production-resolution accuracy.
    auto spec = cp_sbtl::presets::ph_subcritical(*heos, /*NT=*/40, /*NR=*/80, /*rank=*/10);
    REQUIRE(spec.input_pair == ::CoolProp::HmassP_INPUTS);
    // 6 regions: LIQUID, VAPOR, NC_LIQUID, NC_VAPOR, NC_SUPER, SUPER.
    // The NC_* regions are the near-critical sub-regions added in
    // CoolProp-4u9 (POWER / POWER_LO β=1/3 primary axis crowding
    // toward pc).
    REQUIRE(spec.regions.size() == 6);
    // 5 thermodynamic properties (ρ, T, s, u, w) plus 2 optional
    // transport properties (η, λ) when the source backend exposes them
    // — HEOS Water does, so we expect 7 here.  Fluids without transport
    // correlations get only the 5 thermodynamic properties.
    REQUIRE(spec.properties.size() == 7);
    auto spec_has = [&spec](::CoolProp::parameters key) {
        return std::any_of(spec.properties.begin(), spec.properties.end(), [key](const cp_sbtl::PropertySpec& ps) { return ps.key == key; });
    };
    REQUIRE(spec_has(::CoolProp::iDmass));
    REQUIRE(spec_has(::CoolProp::iT));
    REQUIRE(spec_has(::CoolProp::iSmass));
    REQUIRE(spec_has(::CoolProp::iUmass));
    REQUIRE(spec_has(::CoolProp::ispeed_sound));
    REQUIRE(spec_has(::CoolProp::iviscosity));
    REQUIRE(spec_has(::CoolProp::iconductivity));

    auto surface = cp_sbtl::build_surface(*heos, std::move(spec));
    REQUIRE(surface.sealed());
    REQUIRE(surface.region_count() == 6);  // CoolProp-4u9: NC_LIQUID + NC_VAPOR + NC_SUPER

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
        } catch (...) {  // NOLINT(bugprone-empty-catch)
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

TEST_CASE("SVDSurfaceSerializer::default_cache_path rejects path-traversal fluid names", "[SBTL][serializer]") {
    // Defense-in-depth: anything that could escape the cache directory
    // must throw rather than silently writing outside ~/.CoolProp/SVDTables.
    using Catch::Matchers::ContainsSubstring;
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("../etc/passwd", "HEOS", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid fluid_name"));
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("foo/bar", "HEOS", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid fluid_name"));
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("foo\\bar", "HEOS", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid fluid_name"));
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("..", "HEOS", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid fluid_name"));
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("", "HEOS", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid fluid_name"));
    // Same checks apply to the source_backend slot.
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "foo/bar", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid source_backend"));
    REQUIRE_THROWS_WITH(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "", ::CoolProp::HmassP_INPUTS),
                        ContainsSubstring("invalid source_backend"));
    // Legitimate names pass.
    REQUIRE_NOTHROW(cp_sbtl::SVDSurfaceSerializer::default_cache_path("Water", "HEOS", ::CoolProp::HmassP_INPUTS));
    REQUIRE_NOTHROW(cp_sbtl::SVDSurfaceSerializer::default_cache_path("R134a", "REFPROP", ::CoolProp::PT_INPUTS));
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
        const int rc = compress((unsigned char*)bogus_payload.data(), &out_size, (const unsigned char*)hello, sizeof(hello) - 1);
        // miniz returns MZ_OK == Z_OK == 0.  Failing here means the
        // test setup itself is broken (not what we're trying to
        // validate), so REQUIRE rather than silently resizing.
        REQUIRE(rc == 0);
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
            REQUIRE(surface.region_count() == 6);  // CoolProp-4u9: NC_LIQUID + NC_VAPOR + NC_SUPER

            // Round-trip through serializer.  HEOS canonicalises aliases
            // (e.g. "Propane" → "n-Propane"), so compare against what the
            // surface actually stored, not the user-facing input string.
            const auto blob = cp_sbtl::SVDSurfaceSerializer::save(surface);
            auto reloaded = cp_sbtl::SVDSurfaceSerializer::load(blob);
            REQUIRE(reloaded.fluid_name() == surface.fluid_name());
            REQUIRE(reloaded.region_count() == 6);

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
    REQUIRE(spec.regions.size() == 6);  // CoolProp-4u9: NC_LIQUID + NC_VAPOR + NC_SUPER

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

    // Save to a tmp file and load back via the file API.  Use
    // std::filesystem::temp_directory_path() rather than hard-coding
    // /tmp -- the latter doesn't exist on Windows CI.  Suffix with the
    // process PID (see CoolProp-8ft) so two concurrent CatchTestRunner
    // instances don't collide on the same file.
    const auto tmp_path = std::filesystem::temp_directory_path()
                          / ("svd_test_water_ph_" + std::to_string(CoolProp::tests::test_pid()) + ".svd.bin.z");
    const std::string path = tmp_path.string();
    cp_sbtl::SVDSurfaceSerializer::save_to_file(surface, path);
    auto loaded = cp_sbtl::SVDSurfaceSerializer::load_from_file(path);
    REQUIRE(loaded.fluid_name() == "Water");
    REQUIRE(loaded.region_count() == 6);  // CoolProp-4u9: NC_LIQUID + NC_VAPOR + NC_SUPER
    // Spot-check a single eval matches.
    heos->update(::CoolProp::PT_INPUTS, 1e5, 350.0);
    const double h = heos->hmass();
    REQUIRE(loaded.eval(::CoolProp::iDmass, 1e5, h) == surface.eval(::CoolProp::iDmass, 1e5, h));
}

#endif  // ENABLE_CATCH
