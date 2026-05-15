// Correctness and benchmark tests for the generic SVD + Region
// components (include/CoolProp/svd, include/CoolProp/region).
//
// Built into the existing CatchTestRunner target when
// COOLPROP_CATCH_MODULE=ON.  Correctness tests run by default; the
// benchmark cases are tagged [!benchmark] so they only execute when
// invoked explicitly (e.g. `CatchTestRunner "[!benchmark]"`).

#include "CoolProp/region/AxisTransform.h"
#include "CoolProp/region/BoundaryCurve.h"
#include "CoolProp/region/ConstantCurve.h"
#include "CoolProp/region/CubicSplineCurve.h"
#include "CoolProp/region/PiecewiseChebyshevCurve.h"
#include "CoolProp/region/Region.h"
#include "CoolProp/region/RegionAtlas.h"
#include "CoolProp/svd/Hermite1D.h"
#include "CoolProp/svd/SVDBuilder.h"
#include "CoolProp/svd/SVDDecomposition.h"
#include "CoolProp/svd/SVDEvaluator.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>

#    include <cmath>
#    include <random>
#    include <vector>

namespace cp_svd = CoolProp::svd;
namespace cp_region = CoolProp::region;

namespace {
// Local portability constant — M_PI is non-standard on MSVC without
// _USE_MATH_DEFINES.
constexpr double kPi = 3.141592653589793238462643383279502884;
}  // namespace

// ============================================================
// AxisTransform
// ============================================================

TEST_CASE("AxisTransform linear round-trip", "[SVDComponents][Region][axis_transform]") {
    const auto t = cp_region::AxisTransform::make(cp_region::AxisScale::LINEAR, 1.5, 9.5);
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u(1.5, 9.5);
    for (int k = 0; k < 100; ++k) {
        const double a = u(rng);
        const double xi = t.forward(a);
        REQUIRE(xi >= -1e-15);
        REQUIRE(xi <= 1.0 + 1e-15);
        const double a_back = t.inverse(xi);
        REQUIRE(std::abs(a_back - a) / a < 1e-14);
    }
}

TEST_CASE("AxisTransform log round-trip", "[SVDComponents][Region][axis_transform]") {
    const auto t = cp_region::AxisTransform::make(cp_region::AxisScale::LOG, 1e-3, 1e3);
    std::mt19937 rng(2);
    std::uniform_real_distribution<double> u(std::log(1e-3), std::log(1e3));
    for (int k = 0; k < 100; ++k) {
        const double a = std::exp(u(rng));
        const double xi = t.forward(a);
        REQUIRE(xi >= -1e-15);
        REQUIRE(xi <= 1.0 + 1e-15);
        const double a_back = t.inverse(xi);
        REQUIRE(std::abs(a_back - a) / a < 1e-14);
    }
}

TEST_CASE("AxisTransform analytic Jacobian", "[SVDComponents][Region][axis_transform]") {
    const auto t_lin = cp_region::AxisTransform::make(cp_region::AxisScale::LINEAR, 2.0, 10.0);
    REQUIRE(t_lin.dxi_da(5.0) == Catch::Approx(1.0 / 8.0).margin(1e-15));

    const auto t_log = cp_region::AxisTransform::make(cp_region::AxisScale::LOG, 1.0, 100.0);
    // Analytic: dxi/da = 1 / (a * log(100)).  Cross-check by central FD.
    const double a = 10.0;
    const double fd = (t_log.forward(a * (1 + 1e-7)) - t_log.forward(a * (1 - 1e-7))) / (2.0 * a * 1e-7);
    REQUIRE(t_log.dxi_da(a) == Catch::Approx(fd).epsilon(1e-6));
}

// ============================================================
// BoundaryCurve: ConstantCurve
// ============================================================

TEST_CASE("ConstantCurve eval / eval_da / bounds", "[SVDComponents][Region][boundary_curve_constant]") {
    cp_region::ConstantCurve c(0.0, 5.0, 3.14);
    REQUIRE(c.eval(2.7) == 3.14);
    REQUIRE(c.eval_da(2.7) == 0.0);
    REQUIRE(c.bounds().first == 3.14);
    REQUIRE(c.bounds().second == 3.14);
    REQUIRE(c.a_range().first == 0.0);
    REQUIRE(c.a_range().second == 5.0);
}

// ============================================================
// BoundaryCurve: PiecewiseChebyshevCurve
// ============================================================

TEST_CASE("PiecewiseChebyshevCurve fits a smooth function (LINEAR axis)", "[SVDComponents][Region][boundary_curve_cheb]") {
    const auto f = [](double a) { return std::sin(a) + 0.5 * a; };
    auto curve =
      cp_region::PiecewiseChebyshevCurve::build(0.5, 5.5, /*n_pieces=*/4, /*degree=*/12, cp_region::PiecewiseChebyshevCurve::ParamScale::LINEAR, f);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(0.6, 5.4);
    for (int k = 0; k < 100; ++k) {
        const double a = u(rng);
        REQUIRE(curve->eval(a) == Catch::Approx(f(a)).epsilon(1e-10));
        // Analytic eval_da cross-checked against central FD.
        const double h = 1e-6;
        const double fd = (curve->eval(a + h) - curve->eval(a - h)) / (2.0 * h);
        REQUIRE(curve->eval_da(a) == Catch::Approx(fd).epsilon(1e-6));
    }
}

TEST_CASE("PiecewiseChebyshevCurve fits a smooth function (LOG axis)", "[SVDComponents][Region][boundary_curve_cheb]") {
    // Use a function that is smooth in log space.
    const auto f = [](double a) { return std::log(a) * std::log(a); };
    auto curve =
      cp_region::PiecewiseChebyshevCurve::build(1.0, 1e4, /*n_pieces=*/3, /*degree=*/10, cp_region::PiecewiseChebyshevCurve::ParamScale::LOG, f);

    std::mt19937 rng(9);
    std::uniform_real_distribution<double> u(std::log(1.05), std::log(1e4 * 0.99));
    for (int k = 0; k < 100; ++k) {
        const double a = std::exp(u(rng));
        REQUIRE(curve->eval(a) == Catch::Approx(f(a)).epsilon(1e-8));
        const double h = a * 1e-6;
        const double fd = (curve->eval(a + h) - curve->eval(a - h)) / (2.0 * h);
        REQUIRE(curve->eval_da(a) == Catch::Approx(fd).epsilon(1e-5));
    }
}

// ============================================================
// BoundaryCurve: CubicSplineCurve
// ============================================================

TEST_CASE("CubicSplineCurve through known knots", "[SVDComponents][Region][boundary_curve_spline]") {
    // Knots from y = x^2 (a strict-convex test).
    std::vector<double> x = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
    std::vector<double> y;
    for (double xi : x) {
        y.push_back(xi * xi);
    }
    auto curve = cp_region::CubicSplineCurve::build(x, y);

    // Spline interpolates each knot exactly (to roundoff).
    for (std::size_t i = 0; i < x.size(); ++i) {
        REQUIRE(curve->eval(x[i]) == Catch::Approx(y[i]).margin(1e-13));
    }
    // Analytic eval_da matches central FD.
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(0.05, 2.95);
    for (int k = 0; k < 50; ++k) {
        const double a = u(rng);
        const double h = 1e-6;
        const double fd = (curve->eval(a + h) - curve->eval(a - h)) / (2.0 * h);
        REQUIRE(curve->eval_da(a) == Catch::Approx(fd).epsilon(1e-5));
    }
    // Tight bounds enclose the knot values.
    const auto bnds = curve->bounds();
    REQUIRE(bnds.first <= 0.0);
    REQUIRE(bnds.second >= 9.0);
}

// ============================================================
// Region + RegionAtlas: hand-built thermodynamics-free atlas
// ============================================================

namespace {

// Build a deliberately overlapping-AABB but disjoint-curve atlas.
//
//   Primary axis a in [0, 10], linear.
//   Region L:  b_lo = 0,                  b_hi = 0.5 * a       (lower triangle)
//   Region V:  b_lo = 0.5 * a,            b_hi = 10            (upper region above the line)
//
// The two regions are disjoint by construction, but their AABBs both
// span (a in [0,10], b in [0,10]), so the AABB filter admits both.
// curve_contains is what disambiguates.
cp_region::RegionAtlas build_test_atlas() {
    cp_region::RegionAtlas atlas;
    // LIQUID-like.
    {
        auto axis = cp_region::AxisTransform::make(cp_region::AxisScale::LINEAR, 0.0, 10.0);
        auto lo = std::make_unique<cp_region::ConstantCurve>(0.0, 10.0, 0.0);
        auto hi = cp_region::PiecewiseChebyshevCurve::build(0.0, 10.0, 1, 4, cp_region::PiecewiseChebyshevCurve::ParamScale::LINEAR,
                                                            [](double a) { return 0.5 * a; });
        atlas.add(cp_region::Region(axis, std::move(lo), std::move(hi)));
    }
    // VAPOR-like.
    {
        auto axis = cp_region::AxisTransform::make(cp_region::AxisScale::LINEAR, 0.0, 10.0);
        auto lo = cp_region::PiecewiseChebyshevCurve::build(0.0, 10.0, 1, 4, cp_region::PiecewiseChebyshevCurve::ParamScale::LINEAR,
                                                            [](double a) { return 0.5 * a; });
        auto hi = std::make_unique<cp_region::ConstantCurve>(0.0, 10.0, 10.0);
        atlas.add(cp_region::Region(axis, std::move(lo), std::move(hi)));
    }
    return atlas;
}

}  // namespace

TEST_CASE("RegionAtlas dispatch picks the right region", "[SVDComponents][Region][atlas_dispatch]") {
    auto atlas = build_test_atlas();
    REQUIRE(atlas.size() == 2);

    // Clearly inside region 0 (LIQUID-like).
    REQUIRE(atlas.find_region(2.0, 0.3) == 0);
    // Clearly inside region 1 (VAPOR-like).
    REQUIRE(atlas.find_region(2.0, 5.0) == 1);
    // Outside both (a > 10).
    REQUIRE(atlas.find_region(11.0, 5.0) == -1);
    // Outside both (b < 0).
    REQUIRE(atlas.find_region(2.0, -0.1) == -1);
}

TEST_CASE("RegionAtlas AABB filter rejects strictly more than curve filter", "[SVDComponents][Region][aabb_first]") {
    auto atlas = build_test_atlas();
    // For each region we count how many random (a, b) probes pass the
    // AABB filter (would *trigger* a curve check) vs how many pass the
    // curve filter.  Since the AABBs of both regions completely overlap
    // but the curve envelopes are disjoint, AABB hits = 2 * curve hits
    // for points that land inside the union of the two regions.
    std::mt19937 rng(13);
    std::uniform_real_distribution<double> ua(0.0, 10.0);
    std::uniform_real_distribution<double> ub(0.0, 10.0);
    int aabb_hits = 0;
    int curve_hits = 0;
    for (int k = 0; k < 10000; ++k) {
        const double a = ua(rng);
        const double b = ub(rng);
        for (std::size_t i = 0; i < atlas.size(); ++i) {
            if (atlas.region(i).aabb_contains(a, b)) {
                ++aabb_hits;
                if (atlas.region(i).curve_contains(a, b)) {
                    ++curve_hits;
                }
            }
        }
    }
    INFO("aabb_hits=" << aabb_hits << " curve_hits=" << curve_hits);
    REQUIRE(aabb_hits > curve_hits);
    // Sanity: every point that passes a curve passed an AABB.
    REQUIRE(aabb_hits >= curve_hits);
}

TEST_CASE("Region (a, b) <-> (xi, eta) round-trip", "[SVDComponents][Region][round_trip]") {
    auto atlas = build_test_atlas();
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int round_tripped = 0;
    for (int k = 0; k < 200; ++k) {
        const double xi = u(rng);
        const double eta = u(rng);
        for (std::size_t i = 0; i < atlas.size(); ++i) {
            const auto ab = atlas.region(i).from_normalized(xi, eta);
            const auto xe = atlas.region(i).to_normalized(ab.first, ab.second);
            REQUIRE(xe.first == Catch::Approx(xi).epsilon(1e-12));
            REQUIRE(xe.second == Catch::Approx(eta).epsilon(1e-12));
            ++round_tripped;
        }
    }
    REQUIRE(round_tripped > 0);
}

// ============================================================
// SVD end-to-end + slope source comparison
// ============================================================

namespace {

// Test function: f(x, y) = sin(pi x) * (1 + 0.3 sin(pi y)) on [0, 1]^2.
//
// Two reasons for this choice:
//   1. The second derivative of sin(pi z) vanishes at z = 0 and z = 1,
//      so the natural-cubic-spline boundary condition (M = 0 at the
//      endpoints) is *correct* for this function.  That makes spline
//      slopes strictly better than central-FD slopes here, which is
//      what we want to assert.
//   2. The function is exactly separable, so its SVD is dominated by a
//      single mode plus a small correction — rank 4 is overkill for
//      tiny error.
cp_svd::SVDDecomposition build_test_decomposition(cp_svd::SlopeSource slope_source, std::int32_t rank) {
    const int N = 64;
    std::vector<double> xg(N), yg(N), M(N * N);
    for (int i = 0; i < N; ++i) {
        xg[i] = static_cast<double>(i) / (N - 1);
    }
    for (int j = 0; j < N; ++j) {
        yg[j] = static_cast<double>(j) / (N - 1);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            M[i * N + j] = std::sin(kPi * xg[i]) * (1.0 + 0.3 * std::sin(kPi * yg[j]));
        }
    }
    cp_svd::SVDBuildOptions opts;
    opts.rank = rank;
    opts.out_transform = cp_svd::OutputTransform::IDENTITY;
    opts.slope_source = slope_source;
    return cp_svd::build_svd(xg, yg, M, opts);
}

double test_function_truth(double x, double y) {
    return std::sin(kPi * x) * (1.0 + 0.3 * std::sin(kPi * y));
}

}  // namespace

TEST_CASE("SVD round-trip on a smooth analytic 2D function", "[SVDComponents][SVD][roundtrip]") {
    auto d = build_test_decomposition(cp_svd::SlopeSource::NATURAL_CUBIC_SPLINE, 4);
    cp_svd::SVDEvaluator ev(d);
    std::mt19937 rng(19);
    std::uniform_real_distribution<double> u(0.05, 0.95);
    double max_rel = 0.0;
    for (int k = 0; k < 2000; ++k) {
        const double x = u(rng);
        const double y = u(rng);
        const double truth = test_function_truth(x, y);
        const double pred = ev.eval(x, y);
        const double rel = std::abs(pred - truth) / std::max(std::abs(truth), 1e-12);
        if (rel > max_rel) {
            max_rel = rel;
        }
    }
    INFO("max_rel = " << max_rel);
    REQUIRE(max_rel < 1e-4);
}

TEST_CASE("SVD natural-cubic-spline slopes outperform Hermite-FD", "[SVDComponents][SVD][slope_source]") {
    // Compare the two slope sources on a function whose second
    // derivative vanishes at the build-domain boundaries.  Spline's
    // natural BC is then exact, FD's local-stencil approach is not, so
    // spline should produce strictly smaller error.
    auto d_fd = build_test_decomposition(cp_svd::SlopeSource::HERMITE_FD, 4);
    auto d_spline = build_test_decomposition(cp_svd::SlopeSource::NATURAL_CUBIC_SPLINE, 4);
    cp_svd::SVDEvaluator ev_fd(d_fd);
    cp_svd::SVDEvaluator ev_spline(d_spline);
    std::mt19937 rng(21);
    std::uniform_real_distribution<double> u(0.05, 0.95);
    double max_fd = 0.0, max_spline = 0.0;
    for (int k = 0; k < 2000; ++k) {
        const double x = u(rng);
        const double y = u(rng);
        const double truth = test_function_truth(x, y);
        const double rel_fd = std::abs(ev_fd.eval(x, y) - truth) / std::max(std::abs(truth), 1e-12);
        const double rel_spline = std::abs(ev_spline.eval(x, y) - truth) / std::max(std::abs(truth), 1e-12);
        if (rel_fd > max_fd) {
            max_fd = rel_fd;
        }
        if (rel_spline > max_spline) {
            max_spline = rel_spline;
        }
    }
    INFO("max_fd=" << max_fd << " max_spline=" << max_spline);
    REQUIRE(max_spline <= max_fd);
}

TEST_CASE("SVD EXP transform recovers values from log-space build", "[SVDComponents][SVD][roundtrip]") {
    // Build the SVD of log(f) and ask the evaluator to apply exp on the
    // way out — the round-trip should still hit the original f values.
    const int N = 32;
    std::vector<double> xg(N), yg(N), M(N * N);
    for (int i = 0; i < N; ++i) {
        xg[i] = 0.1 + 4.0 * i / (N - 1);
    }
    for (int j = 0; j < N; ++j) {
        yg[j] = 0.1 + 4.0 * j / (N - 1);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // f is strictly positive; use log here.
            M[i * N + j] = std::log(xg[i] * xg[i] + yg[j] * yg[j] + 1.0);
        }
    }
    cp_svd::SVDBuildOptions opts;
    opts.rank = 6;
    opts.out_transform = cp_svd::OutputTransform::EXP;
    auto d = cp_svd::build_svd(xg, yg, M, opts);
    cp_svd::SVDEvaluator ev(d);
    std::mt19937 rng(23);
    std::uniform_real_distribution<double> u(0.2, 3.9);
    double max_rel = 0.0;
    for (int k = 0; k < 1000; ++k) {
        const double x = u(rng);
        const double y = u(rng);
        const double truth_log = std::log(x * x + y * y + 1.0);
        const double truth = std::exp(truth_log);
        const double pred = ev.eval(x, y);
        const double rel = std::abs(pred - truth) / truth;
        if (rel > max_rel) {
            max_rel = rel;
        }
    }
    REQUIRE(max_rel < 1e-2);
}

// ============================================================
// Benchmarks (opt-in via [!benchmark])
// ============================================================

TEST_CASE("SVDEvaluator::eval timing", "[SVDComponents][!benchmark]") {
    auto d = build_test_decomposition(cp_svd::SlopeSource::NATURAL_CUBIC_SPLINE, 20);
    cp_svd::SVDEvaluator ev(d);
    std::mt19937 rng(101);
    std::uniform_real_distribution<double> u(0.05, 0.95);
    std::vector<std::pair<double, double>> probes;
    probes.reserve(1000);
    for (int i = 0; i < 1000; ++i) {
        probes.emplace_back(u(rng), u(rng));
    }

    BENCHMARK("eval rank-20") {
        double sink = 0.0;
        for (const auto& p : probes) {
            sink += ev.eval(p.first, p.second);
        }
        return sink;
    };
}

TEST_CASE("RegionAtlas::find_region timing", "[SVDComponents][!benchmark]") {
    // Use the same hand-built 2-region atlas as the correctness tests.
    //
    //   Region 0 (LIQUID-like):  b_lo = 0,         b_hi = 0.5*a  on a in [0, 10]
    //                            AABB = [0, 10] x [0, 5]
    //   Region 1 (VAPOR-like):   b_lo = 0.5*a,     b_hi = 10     on a in [0, 10]
    //                            AABB = [0, 10] x [0, 10]
    //
    // Three benchmark scenarios that isolate the real per-call work:
    //
    //   1. atlas_miss     — point outside both AABBs: scans 2 AABBs,
    //                       no curve_contains call, returns -1.
    //   2. first_region   — point inside region 0's curve envelope:
    //                       1 AABB hit + 1 curve_contains, returns 0.
    //   3. last_region    — point inside region 1's curve envelope
    //                       AND inside region 0's AABB: 2 AABB hits +
    //                       2 curve_contains, returns 1.
    //
    // The cost of a curve_contains call is dominated by two
    // PiecewiseChebyshevCurve Clenshaw evaluations (one per side).
    auto atlas = build_test_atlas();
    constexpr int N_PROBES = 1000;
    std::mt19937 rng(103);
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    std::vector<std::pair<double, double>> miss_probes;
    std::vector<std::pair<double, double>> first_probes;
    std::vector<std::pair<double, double>> last_probes;
    miss_probes.reserve(N_PROBES);
    first_probes.reserve(N_PROBES);
    last_probes.reserve(N_PROBES);
    while (static_cast<int>(miss_probes.size()) < N_PROBES) {
        // Outside both AABBs: a in [11, 20], b anything.
        miss_probes.emplace_back(11.0 + 9.0 * u01(rng), 10.0 * u01(rng));
    }
    while (static_cast<int>(first_probes.size()) < N_PROBES) {
        // Inside region 0: 0 <= b <= 0.5 a.  Sample a uniformly, then
        // b uniformly in [0, 0.5 a] so probes are guaranteed in region 0.
        const double a = 10.0 * u01(rng);
        const double b = 0.5 * a * u01(rng);
        first_probes.emplace_back(a, b);
    }
    while (static_cast<int>(last_probes.size()) < N_PROBES) {
        // Inside region 1 but ALSO inside region 0's AABB ([0,10] x [0,5]):
        //   pick a in [0, 10], b in (0.5 a, 5].  If 0.5 a >= 5 the interval
        //   is empty; resample.
        const double a = 10.0 * u01(rng);
        const double b_lo = 0.5 * a;
        if (b_lo >= 5.0) {
            continue;
        }
        const double b = b_lo + (5.0 - b_lo) * u01(rng);
        last_probes.emplace_back(a, b);
    }

    BENCHMARK("find_region — atlas-miss (AABB filter only)") {
        int sink = 0;
        for (const auto& p : miss_probes) {
            sink += atlas.find_region(p.first, p.second);
        }
        return sink;
    };

    BENCHMARK("find_region — first-region hit (1 AABB + 1 curve)") {
        int sink = 0;
        for (const auto& p : first_probes) {
            sink += atlas.find_region(p.first, p.second);
        }
        return sink;
    };

    BENCHMARK("find_region — last-region hit (2 AABBs + 2 curves)") {
        int sink = 0;
        for (const auto& p : last_probes) {
            sink += atlas.find_region(p.first, p.second);
        }
        return sink;
    };
}

#endif  // ENABLE_CATCH
