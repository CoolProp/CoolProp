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
            CHECK(err < 1e-10 * std::max(1.0, sa));  // the bound itself is small
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

#endif  // ENABLE_CATCH
