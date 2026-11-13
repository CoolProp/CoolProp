// Tests for the certified Chebyshev-series root isolation in CoolProp/numerics/ChebyshevBernstein.h.
//
// The property every test leans on: every root of the series in [-1, 1] lies in [ua, ub] of some reported
// entry, and a certified entry holds exactly one simple root.  A root must never disappear silently -- an
// interval the subdivision cannot resolve (depth cap, node budget, coefficients within the tolerance) is
// reported as uncertified instead.

#if defined(ENABLE_CATCH)

#    include "CoolProp/numerics/ChebyshevBernstein.h"

#    include <catch2/catch_all.hpp>

#    include <algorithm>
#    include <cmath>
#    include <cstdint>
#    include <limits>
#    include <random>
#    include <vector>

namespace CB = CoolProp::ChebyshevBernstein;

namespace {

/// Chebyshev coefficients of the polynomial scale * prod (u - r_i), padded to degree N (exact up to roundoff:
/// interpolation at N + 1 Chebyshev points of a polynomial of degree <= N)
template <int N>
CB::Coeffs<N> from_roots(const std::vector<double>& roots, double scale = 1.0) {
    const double PI = 3.14159265358979323846;
    CB::Coeffs<N> fv{}, c{};
    for (int j = 0; j <= N; ++j) {
        const double u = std::cos(PI * (j + 0.5) / (N + 1));
        double v = scale;
        for (double r : roots)
            v *= (u - r);
        fv[j] = v;
    }
    for (int k = 0; k <= N; ++k) {
        double s = 0;
        for (int j = 0; j <= N; ++j)
            s += fv[j] * std::cos(PI * k * (j + 0.5) / (N + 1));
        c[k] = s * (k == 0 ? 1.0 : 2.0) / (N + 1);
    }
    return c;
}

/// Bernstein polynomial value at u in [-1, 1] by de Casteljau
template <int N>
double bernstein_value(const CB::Coeffs<N>& b, double u) {
    const double t = 0.5 * (u + 1);
    CB::Coeffs<N> w = b;
    for (int r = 1; r <= N; ++r)
        for (int i = 0; i <= N - r; ++i)
            w[i] = (1 - t) * w[i] + t * w[i + 1];
    return w[0];
}

bool covered(const std::vector<CB::Root>& out, double r) {
    return std::any_of(out.begin(), out.end(), [r](const CB::Root& e) { return e.ua <= r && r <= e.ub; });
}

}  // namespace

TEST_CASE("Chebyshev-to-Bernstein matrices are exact for low degree", "[cheb_bernstein]") {
    // degree 2 on [-1, 1], t = (u + 1)/2: T0 = 1 -> (1, 1, 1); T1 = 2t - 1 -> (-1, 0, 1);
    // T2 = 8t^2 - 8t + 1 -> (f(0), a0 + a1/2, f(1)) = (1, -3, 1)
    const auto& M = CB::cheb2bern_matrix<2>();
    const double expect[3][3] = {{1, -1, 1}, {1, 0, -3}, {1, 1, 1}};  // M[j][k]
    for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
            CHECK(M[j][k] == expect[j][k]);
}

TEST_CASE("Bernstein conversion reproduces the Chebyshev series", "[cheb_bernstein]") {
    std::mt19937_64 g(7);
    std::uniform_real_distribution<double> U(-1, 1);
    auto check = [&](auto tag) {
        constexpr int N = decltype(tag)::value;
        for (int trial = 0; trial < 50; ++trial) {
            CB::Coeffs<N> c{};
            double sa = 0;
            for (auto& v : c)
                sa += std::abs(v = U(g));
            double err = 0;
            const auto b = CB::to_bernstein<N>(c, err);
            CHECK(err < 1e-10 * std::max(1.0, sa));  // the bound is small ...
            const auto& M = CB::cheb2bern_matrix<N>();
            for (int j = 0; j <= N; ++j) {  // ... and it bounds the actual error (reference: double-double M c)
                double hi = 0, lo = 0;
                for (int k = 0; k <= N; ++k) {
                    const double p = M[j][k] * c[k], pe = std::fma(M[j][k], c[k], -p);   // exact product = p + pe
                    const double t = hi + p, z = t - hi, se = (hi - (t - z)) + (p - z);  // TwoSum
                    hi = t;
                    lo += se + pe;
                }
                CHECK(std::abs(b[j] - (hi + lo)) <= err);
            }
            for (int k = 0; k <= 20; ++k) {
                const double u = -1 + 2.0 * k / 20;
                CHECK(std::abs(bernstein_value<N>(b, u) - CB::clenshaw<N>(c, u)) <= 1e-12 * (1 + sa));
            }
        }
    };
    check(std::integral_constant<int, 1>{});
    check(std::integral_constant<int, 5>{});
    check(std::integral_constant<int, 12>{});
    check(std::integral_constant<int, 17>{});
}

TEST_CASE("Chebyshev derivative and Clenshaw derivative agree", "[cheb_bernstein]") {
    const auto c = from_roots<17>({-0.9, -0.3, 0.2, 0.7, 0.95});
    const auto d = CB::derivative<17>(c);
    for (int k = 0; k <= 40; ++k) {
        const double u = -1 + 2.0 * k / 40;
        double f, df;
        CB::clenshaw_fd<17>(c, u, f, df);
        CHECK(f == Catch::Approx(CB::clenshaw<17>(c, u)).margin(1e-14));
        CHECK(df == Catch::Approx(CB::clenshaw<17>(d, u)).margin(1e-12));
        const double h = 1e-6;
        CHECK(df == Catch::Approx((CB::clenshaw<17>(c, u + h) - CB::clenshaw<17>(c, u - h)) / (2 * h)).margin(1e-6));
    }
}

TEST_CASE("Simple roots are certified and accurate", "[cheb_bernstein]") {
    const std::vector<double> roots = {-0.7, 0.1, 0.55};
    std::vector<CB::Root> out;
    CB::real_roots<17>(from_roots<17>(roots), 0.0, out);
    REQUIRE(out.size() == roots.size());
    for (std::size_t i = 0; i < roots.size(); ++i) {
        CHECK(out[i].certified);
        CHECK(out[i].sign_change);
        CHECK(out[i].u == Catch::Approx(roots[i]).margin(1e-12));
        CHECK((out[i].ua <= roots[i] && roots[i] <= out[i].ub));
    }
}

TEST_CASE("No roots: nothing is reported", "[cheb_bernstein]") {
    std::vector<CB::Root> out;
    CB::Coeffs<17> c{};  // u^2 + 0.5 = 0.5 T2 + 1
    c[0] = 1.0;
    c[2] = 0.5;
    CB::real_roots<17>(c, 0.0, out);
    CHECK(out.empty());
}

TEST_CASE("Close pairs are separated or reported, never dropped", "[cheb_bernstein]") {
    for (double gap : {1e-3, 1e-6, 1e-9, 1e-12, 1e-15}) {
        CAPTURE(gap);
        const std::vector<double> roots = {0.3, 0.3 + gap, -0.4};
        std::vector<CB::Root> out;
        CB::real_roots<17>(from_roots<17>(roots), 0.0, out);
        for (double r : roots)
            CHECK(covered(out, r));
        // Between the pair the polynomial peaks at ~0.7 (gap/2)^2: ~2e-13 for gap = 1e-6, but ~2e-19 for
        // gap = 1e-9 -- below the roundoff of the coefficients themselves, so no double-precision method can
        // certify that pair.  Resolvable pairs must come out as two certified roots; the others only covered.
        if (gap >= 1e-6) {
            const auto n = std::count_if(out.begin(), out.end(), [](const CB::Root& e) { return e.certified && std::abs(e.u - 0.3) < 2e-3; });
            CHECK(n == 2);
        } else {
            CHECK(std::any_of(out.begin(), out.end(), [](const CB::Root& e) { return !e.certified && e.ua <= 0.3 && 0.3 <= e.ub; }));
        }
    }
}

TEST_CASE("A double root is reported as an uncertified interval without a sign change", "[cheb_bernstein]") {
    std::vector<CB::Root> out;
    CB::real_roots<17>(from_roots<17>({0.2, 0.2, -0.5}), 0.0, out);
    CHECK(covered(out, 0.2));
    CHECK(covered(out, -0.5));
    bool saw_double = false;
    for (const auto& e : out) {
        if (e.ua <= 0.2 && 0.2 <= e.ub) {
            CHECK_FALSE(e.certified);
            saw_double = true;
        }
        if (e.ua <= -0.5 && -0.5 <= e.ub) CHECK(e.certified);
    }
    CHECK(saw_double);
}

TEST_CASE("Roots at the interval ends are reported", "[cheb_bernstein]") {
    std::vector<CB::Root> out;
    CB::real_roots<17>(from_roots<17>({1.0, -0.3}), 0.0, out);
    CHECK(covered(out, 1.0));
    CHECK(covered(out, -0.3));
    CB::real_roots<17>(from_roots<17>({-1.0, 0.6}), 0.0, out);
    CHECK(covered(out, -1.0));
    CHECK(covered(out, 0.6));
}

TEST_CASE("A root exactly at a subdivision point is reported once", "[cheb_bernstein]") {
    // u = 0 is the first split point; the three roots force subdivision
    std::vector<CB::Root> out;
    CB::real_roots<17>(from_roots<17>({-0.5, 0.0, 0.5}), 0.0, out);
    const auto n0 = std::count_if(out.begin(), out.end(), [](const CB::Root& e) { return e.ua <= 0.0 && 0.0 <= e.ub; });
    CHECK(n0 == 1);
    CHECK(covered(out, -0.5));
    CHECK(covered(out, 0.5));
}

TEST_CASE("Tolerance semantics: coefficients within tol have unknown sign", "[cheb_bernstein]") {
    // f = 1e-9 (u - 0.3): certified with tol = 0; with tol = 1e-6 every coefficient is ambiguous, so the
    // root can only be reported as unresolved -- but it must still be reported
    const auto c = from_roots<17>({0.3}, 1e-9);
    std::vector<CB::Root> out;
    CB::real_roots<17>(c, 0.0, out);
    REQUIRE(out.size() == 1);
    CHECK(out[0].certified);
    CHECK(out[0].u == Catch::Approx(0.3).margin(1e-12));

    const auto st = CB::real_roots<17>(c, 1e-6, out);
    CHECK(covered(out, 0.3));
    CHECK(std::none_of(out.begin(), out.end(), [](const CB::Root& e) { return e.certified; }));
    CHECK(st.unresolved >= 1);
}

TEST_CASE("Running out of node budget reports intervals instead of dropping roots", "[cheb_bernstein]") {
    const std::vector<double> roots = {-0.9, -0.62, -0.31, -0.05, 0.2, 0.41, 0.66, 0.88};
    CB::Options opt;
    opt.node_budget = 3;
    std::vector<CB::Root> out;
    const auto st = CB::real_roots<17>(from_roots<17>(roots), 0.0, out, opt);
    CHECK(st.budget_exhausted);
    for (double r : roots)
        CHECK(covered(out, r));
}

TEST_CASE("Randomized: every root is covered; certified intervals hold exactly one root", "[cheb_bernstein]") {
    std::mt19937_64 g(12345);
    std::uniform_real_distribution<double> U(-1, 1);
    std::uniform_int_distribution<int> K(1, 17);
    long certified = 0, total = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        const int k = K(g);
        std::vector<double> roots(k);
        for (auto& r : roots)
            r = U(g);
        std::sort(roots.begin(), roots.end());
        std::vector<CB::Root> out;
        CB::real_roots<17>(from_roots<17>(roots), 0.0, out);
        for (double r : roots) {
            CAPTURE(trial, r);
            CHECK(covered(out, r));
        }
        for (const auto& e : out) {
            ++total;
            if (!e.certified) continue;
            ++certified;
            const auto n = std::count_if(roots.begin(), roots.end(), [&](double r) { return e.ua <= r && r <= e.ub; });
            CAPTURE(trial, e.ua, e.ub);
            CHECK(n == 1);
        }
    }
    CHECK(certified > 0.99 * total);  // random roots are almost always separable
}

TEST_CASE("Invalid tolerances and non-finite input fail closed", "[cheb_bernstein]") {
    const auto c = from_roots<17>({0.3});
    std::vector<CB::Root> out;
    // a negative tolerance is treated as 0: the root is still found (it used to be excluded silently)
    CB::real_roots<17>(c, -1.5, out);
    CHECK(covered(out, 0.3));
    CB::real_roots<17>(c, -1e-15, out);
    CHECK(covered(out, 0.3));
    CB::real_roots<17>(c, -std::numeric_limits<double>::infinity(), out);
    CHECK(covered(out, 0.3));
    // NaN / infinite tolerance or a non-finite coefficient: nothing can be decided -> [-1, 1], uncertified
    for (double t : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        const auto st = CB::real_roots<17>(c, t, out);
        REQUIRE(out.size() == 1);
        CHECK_FALSE(out[0].certified);
        CHECK((out[0].ua == -1.0 && out[0].ub == 1.0));
        CHECK(st.nodes == 0);
    }
    auto cn = c;
    cn[5] = std::numeric_limits<double>::quiet_NaN();
    CB::real_roots<17>(cn, 0.0, out);
    REQUIRE(out.size() == 1);
    CHECK((out[0].ua == -1.0 && out[0].ub == 1.0 && !out[0].certified));
    // a root just inside the end, found with tol = 0 but lost with a tiny negative tol before the fix
    CB::Coeffs<1> c1{0.5 - std::ldexp(1.0, -53), 0.5};
    CB::real_roots<1>(c1, -1e-15, out);
    CHECK(covered(out, -(0.5 - std::ldexp(1.0, -53)) / 0.5));
}

TEST_CASE("The zero series is reported as unresolved and the statistics match the output", "[cheb_bernstein]") {
    std::vector<CB::Root> out;
    const auto st = CB::real_roots<17>(CB::Coeffs<17>{}, 0.0, out);
    CHECK(covered(out, -1.0));
    CHECK(covered(out, 0.0));
    CHECK(covered(out, 1.0));
    CHECK(std::none_of(out.begin(), out.end(), [](const CB::Root& e) { return e.certified; }));
    CHECK(st.unresolved == static_cast<long>(out.size()));
}

TEST_CASE("sign_change is only claimed where the end values are clearly of opposite sign", "[cheb_bernstein]") {
    // a simple root exactly at u = -1: the series is zero there, so no sign change can be certified -- but the
    // root must still be covered
    std::vector<CB::Root> out;
    CB::real_roots<3>(CB::Coeffs<3>{0.5, 0.5, 0, 0}, 0.0, out);  // 0.5 (u + 1)
    CHECK(covered(out, -1.0));
    for (const auto& e : out) {
        if (!e.sign_change) continue;
        const double fa = CB::clenshaw<3>(CB::Coeffs<3>{0.5, 0.5, 0, 0}, e.ua), fb = CB::clenshaw<3>(CB::Coeffs<3>{0.5, 0.5, 0, 0}, e.ub);
        CHECK((fa < 0) != (fb < 0));
    }
    // every entry claiming a sign change, here and for random series, really has opposite end signs
    std::mt19937_64 g(99);
    std::uniform_real_distribution<double> U(-1, 1);
    for (int trial = 0; trial < 300; ++trial) {
        std::vector<double> roots(6);
        for (auto& r : roots)
            r = U(g);
        roots[1] = roots[0] + 1e-12;  // force some unresolved intervals
        const auto c = from_roots<17>(roots);
        CB::real_roots<17>(c, 0.0, out);
        for (const auto& e : out)
            if (e.sign_change) CHECK((CB::clenshaw<17>(c, e.ua) < 0) != (CB::clenshaw<17>(c, e.ub) < 0));
    }
}

namespace {
template <int N>
void randomized_coverage(std::uint64_t seed, int trials) {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> U(-1, 1);
    std::uniform_int_distribution<int> K(1, N);
    for (int trial = 0; trial < trials; ++trial) {
        std::vector<double> roots(K(g));
        for (auto& r : roots)
            r = U(g);
        std::vector<CB::Root> out;
        CB::real_roots<N>(from_roots<N>(roots), 0.0, out);
        for (double r : roots) {
            CAPTURE(N, trial, r);
            CHECK(covered(out, r));
        }
        for (const auto& e : out)
            if (e.certified) CHECK(std::count_if(roots.begin(), roots.end(), [&](double r) { return e.ua <= r && r <= e.ub; }) == 1);
    }
}
}  // namespace

TEST_CASE("Randomized coverage at several degrees", "[cheb_bernstein]") {
    randomized_coverage<1>(1, 2000);
    randomized_coverage<2>(2, 2000);
    randomized_coverage<3>(3, 2000);
    randomized_coverage<5>(5, 2000);
    randomized_coverage<9>(9, 1000);
    randomized_coverage<17>(17, 500);
}

TEST_CASE("With tol > 0, every root of a function within tol is covered", "[cheb_bernstein]") {
    // g = f + tol * sin(K u) stays within tol of the series f; find g's roots by dense sampling and check that
    // each lies in a reported interval, and that g changes sign across every certified interval
    std::mt19937_64 gen(4242);
    std::uniform_real_distribution<double> U(-1, 1);
    const double tol = 1e-4;
    long n_certified = 0;
    for (int trial = 0; trial < 200; ++trial) {
        std::vector<double> roots(5);
        for (auto& r : roots)
            r = U(gen);
        const auto c = from_roots<17>(roots, 0.05);
        const double K = 50 + 500 * (U(gen) + 1);
        auto gfun = [&](double u) { return CB::clenshaw<17>(c, u) + tol * std::sin(K * u); };
        std::vector<CB::Root> out;
        CB::real_roots<17>(c, tol, out);
        const int M = 200000;
        double prev = gfun(-1.0);
        for (int i = 1; i <= M; ++i) {
            const double u = -1 + 2.0 * i / M, v = gfun(u);
            if ((v < 0) != (prev < 0)) {
                CAPTURE(trial, u);  // a root of g lies in [u - 2/M, u]
                CHECK(std::any_of(out.begin(), out.end(),
                                  [&](const CB::Root& e) { return e.ua <= u && u - 2.0 / M <= e.ub; }));  // bracket [u - 2/M, u] meets the interval
            }
            prev = v;
        }
        for (const auto& e : out)
            if (e.certified) {
                ++n_certified;
                CHECK((gfun(e.ua) < 0) != (gfun(e.ub) < 0));
            }
    }
    CHECK(n_certified > 50);  // the test must exercise certified intervals, not only coverage by wide unresolved ones
}

#endif  // ENABLE_CATCH
