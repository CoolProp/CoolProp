#pragma once
/**
 * Certified real-root isolation for Chebyshev series on [-1, 1].
 *
 * A degree-N Chebyshev series c_0 T_0(u) + ... + c_N T_N(u) is converted to the Bernstein basis on [-1, 1]
 * with an exact (rational, rounded once) conversion matrix.  Each interval is then tested on its Bernstein
 * coefficients b_j, whose signs are only trusted where |b_j| exceeds a tolerance (the caller's tolerance plus
 * a bound on the conversion and subdivision roundoff):
 *
 *  - convex hull: all b_j clearly of one sign -> the interval has no root;
 *  - Descartes' rule in the Bernstein basis: the number of sign changes V of (b_0, ..., b_N) bounds the number
 *    of roots and has the same parity; V = 0 -> no root;
 *  - V = 1 with clear signs throughout -> exactly one root, certified, refined by safeguarded Newton;
 *  - otherwise the interval is split at its midpoint by de Casteljau, which yields the Bernstein coefficients
 *    of both halves directly, and each half is tested again.
 *
 * Recursion stops at a depth cap (lower when a coefficient's sign is unknown, since subdividing cannot resolve
 * that) or when a node budget runs out.  Such an interval is never dropped: it is reported as uncertified,
 * with a flag telling whether the series changes sign across it.  So every root of the series in [-1, 1] lies
 * in a reported interval -- certified or not -- and discarding uncertified intervals is the caller's
 * decision, made explicitly.
 *
 * Costs: O(N^2) for the conversion and for each split.
 */

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <vector>

#include "CoolProp/numerics/cheb2bern_tables.h"

namespace CoolProp {
namespace ChebyshevBernstein {

/// Highest degree for which a conversion matrix is tabulated
constexpr int MAX_DEGREE = 17;

/// Chebyshev (or Bernstein) coefficients of a degree-N series, lowest order first
template <int N>
using Coeffs = std::array<double, N + 1>;

/// Value of the Chebyshev series at u in [-1, 1] (Clenshaw)
template <int N>
double clenshaw(const Coeffs<N>& c, double u) {
    double b1 = 0, b2 = 0;
    for (int k = N; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    return c[0] + u * b1 - b2;
}

/// Value f and derivative df/du of the Chebyshev series at u (Clenshaw with the differentiated recurrence)
template <int N>
void clenshaw_fd(const Coeffs<N>& c, double u, double& f, double& df) {
    double b1 = 0, b2 = 0, d1 = 0, d2 = 0;
    for (int k = N; k >= 1; --k) {
        const double b0 = c[k] + 2 * u * b1 - b2, d0 = 2 * b1 + 2 * u * d1 - d2;
        b2 = b1;
        b1 = b0;
        d2 = d1;
        d1 = d0;
    }
    f = c[0] + u * b1 - b2;
    df = b1 + u * d1 - d2;
}

/// sum_{k >= 1} |c_k|: |f(u) - c_0| <= l1_tail(c) on [-1, 1], so |c_0| > l1_tail + tol excludes a root cheaply
template <int N>
double l1_tail(const Coeffs<N>& c) {
    double s = 0;
    for (int k = 1; k <= N; ++k)
        s += std::abs(c[k]);
    return s;
}

/// Chebyshev coefficients of df/du (degree N - 1; the last entry is zero)
template <int N>
Coeffs<N> derivative(const Coeffs<N>& c) {
    Coeffs<N> d{};
    if (N == 0) return d;
    d[N - 1] = 2 * N * c[N];
    if (N >= 2) d[N - 2] = 2 * (N - 1) * c[N - 1];
    for (int k = N - 3; k >= 0; --k)
        d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
    d[0] *= 0.5;
    return d;
}

namespace detail {
template <int N>
struct Table;
#define COOLPROP_CHEB2BERN_TABLE(n)                     \
    template <>                                         \
    struct Table<n>                                     \
    {                                                   \
        static constexpr const auto& M = CHEB2BERN_##n; \
    };
COOLPROP_CHEB2BERN_TABLE(1)
COOLPROP_CHEB2BERN_TABLE(2)
COOLPROP_CHEB2BERN_TABLE(3)
COOLPROP_CHEB2BERN_TABLE(4)
COOLPROP_CHEB2BERN_TABLE(5)
COOLPROP_CHEB2BERN_TABLE(6)
COOLPROP_CHEB2BERN_TABLE(7)
COOLPROP_CHEB2BERN_TABLE(8)
COOLPROP_CHEB2BERN_TABLE(9)
COOLPROP_CHEB2BERN_TABLE(10)
COOLPROP_CHEB2BERN_TABLE(11)
COOLPROP_CHEB2BERN_TABLE(12)
COOLPROP_CHEB2BERN_TABLE(13)
COOLPROP_CHEB2BERN_TABLE(14)
COOLPROP_CHEB2BERN_TABLE(15)
COOLPROP_CHEB2BERN_TABLE(16)
COOLPROP_CHEB2BERN_TABLE(17)
#undef COOLPROP_CHEB2BERN_TABLE
}  // namespace detail

/// The exact (rounded once) Chebyshev-to-Bernstein conversion matrix of degree N: b = M c
template <int N>
constexpr const double (&cheb2bern_matrix())[N + 1][N + 1] {
    static_assert(N >= 1 && N <= MAX_DEGREE, "Chebyshev-to-Bernstein matrices are tabulated for degrees 1..17");
    return detail::Table<N>::M;
}

/// Bernstein coefficients (degree N, on [-1, 1]) of a Chebyshev series, and a bound on the rounding error of
/// each one: |b_j - exact_j| <= err for every j.
template <int N>
Coeffs<N> to_bernstein(const Coeffs<N>& c, double& err) {
    const auto& M = cheb2bern_matrix<N>();
    Coeffs<N> b{};
    err = 0;
    for (int j = 0; j <= N; ++j) {
        double s = 0, sa = 0;
        for (int k = 0; k <= N; ++k) {
            s += M[j][k] * c[k];
            sa += std::abs(M[j][k] * c[k]);
        }
        b[j] = s;
        // rounding of the N + 1 products and sums, plus the once-rounded matrix entries
        err = std::max(err, (N + 2) * DBL_EPSILON * sa);
    }
    return b;
}

/// Safeguarded Newton on the Chebyshev series inside a sign-change bracket [a, b]; fa = f(a).  Every step
/// that leaves the (shrinking) bracket is replaced by bisection.
template <int N>
double refine(const Coeffs<N>& c, double a, double b, double fa, double xtol) {
    double x = 0.5 * (a + b);
    for (int it = 0; it < 100; ++it) {
        double f, df;
        clenshaw_fd<N>(c, x, f, df);
        if (f == 0) return x;
        if ((f < 0) == (fa < 0))
            a = x;
        else
            b = x;
        double xn = x - f / df;
        if (!(xn > a && xn < b)) xn = 0.5 * (a + b);
        const double tolx = xtol * std::max(1.0, std::abs(x));
        if (std::abs(xn - x) <= tolx || b - a <= tolx) return xn;
        x = xn;
    }
    return x;
}

/// One root or root-containing interval in u
struct Root
{
    double u;          ///< refined root (certified, or uncertified with a sign change); interval midpoint otherwise
    double ua, ub;     ///< the interval known to contain it
    bool certified;    ///< exactly one simple root in [ua, ub], given the tolerance
    bool sign_change;  ///< the series changes sign across [ua, ub] (always true when certified)
};

struct Options
{
    int max_depth_ambiguous = 16;  ///< depth cap once a coefficient's sign is unknown (subdividing cannot resolve it)
    int max_depth = 48;            ///< depth cap for separating close roots (2^-48 of the interval is ~ulp)
    long node_budget = 2000;       ///< total subdivision nodes per call
    double xtol = 1e-12;           ///< relative tolerance of the Newton refinement in u
};

struct Stats
{
    long nodes = 0;                 ///< subdivision nodes visited
    long unresolved = 0;            ///< intervals reported uncertified
    bool budget_exhausted = false;  ///< the node budget ran out (some intervals were reported unresolved for that reason)
};

namespace detail {
template <int N>
struct Isolator
{
    const Coeffs<N>& c;
    double tol0;  // caller tolerance + conversion roundoff
    const Options& opt;
    std::vector<Root>& out;
    Stats st;
    long budget;

    void report_unresolved(double ua, double ub) {
        ++st.unresolved;
        const double fa = clenshaw<N>(c, ua), fb = clenshaw<N>(c, ub);
        const bool sc = (fa < 0) != (fb < 0) && fa != 0 && fb != 0;
        const double u = sc ? refine<N>(c, ua, ub, fa, opt.xtol) : 0.5 * (ua + ub);
        // merge with a touching unresolved interval (e.g. a root at a split point seen from both sides)
        if (!out.empty() && !out.back().certified && out.back().ub == ua) {
            Root& r = out.back();
            r.ub = ub;
            const double fl = clenshaw<N>(c, r.ua);
            r.sign_change = (fl < 0) != (fb < 0) && fl != 0 && fb != 0;
            r.u = r.sign_change ? refine<N>(c, r.ua, r.ub, fl, opt.xtol) : 0.5 * (r.ua + r.ub);
            return;
        }
        out.push_back({u, ua, ub, false, sc});
    }

    void rec(const Coeffs<N>& b, double ua, double ub, int depth) {
        ++st.nodes;
        // De Casteljau averages are convex combinations: each level adds at most ~eps * max|b| of rounding
        double bmax = 0;
        for (double v : b)
            bmax = std::max(bmax, std::abs(v));
        const double tol = tol0 + depth * DBL_EPSILON * bmax;
        bool amb = false, anypos = false, anyneg = false;
        int V = 0, last = 0;
        for (int j = 0; j <= N; ++j) {
            const int sgn = b[j] > tol ? 1 : b[j] < -tol ? -1 : 0;
            if (sgn == 0) {
                amb = true;
                continue;
            }
            (sgn > 0 ? anypos : anyneg) = true;
            if (last != 0 && sgn != last) ++V;
            last = sgn;
        }
        if (!amb && !(anypos && anyneg)) return;  // convex hull excludes zero
        if (!amb && V == 0) return;               // Descartes: no root
        if (!amb && V == 1) {                     // Descartes: exactly one root (ends then have opposite signs)
            out.push_back({refine<N>(c, ua, ub, clenshaw<N>(c, ua), opt.xtol), ua, ub, true, true});
            return;
        }
        if ((amb && depth >= opt.max_depth_ambiguous) || depth >= opt.max_depth) {
            report_unresolved(ua, ub);
            return;
        }
        if (--budget < 0) {
            st.budget_exhausted = true;
            report_unresolved(ua, ub);
            return;
        }
        Coeffs<N> L{}, R{}, t = b;
        L[0] = b[0];
        R[N] = b[N];
        for (int r = 1; r <= N; ++r) {
            for (int i = 0; i <= N - r; ++i)
                t[i] = 0.5 * (t[i] + t[i + 1]);
            L[r] = t[0];
            R[N - r] = t[N - r];
        }
        const double um = 0.5 * (ua + ub);
        rec(L, ua, um, depth + 1);
        rec(R, um, ub, depth + 1);
    }
};
}  // namespace detail

/**
 * All real roots of the degree-N Chebyshev series c on [-1, 1], in increasing order.
 *
 * @param c    Chebyshev coefficients
 * @param tol  absolute tolerance on the series values: Bernstein coefficients within tol (plus the roundoff
 *             bound) of zero have unknown sign.  Use the bound on |series - function| when the series
 *             approximates another function and the roots must be certified for that function; 0 otherwise.
 * @param out  receives the roots / root intervals (cleared first)
 * @return     the statistics of the call
 *
 * Guarantee: every root of the series in [-1, 1] lies in [ua, ub] of some reported entry.  Certified entries
 * contain exactly one simple root.
 */
template <int N>
Stats real_roots(const Coeffs<N>& c, double tol, std::vector<Root>& out, const Options& opt = Options()) {
    static_assert(N >= 1 && N <= MAX_DEGREE, "Chebyshev-to-Bernstein matrices are tabulated for degrees 1..17");
    out.clear();
    double err = 0;
    const Coeffs<N> b = to_bernstein<N>(c, err);
    detail::Isolator<N> iso{c, tol + err, opt, out, Stats{}, opt.node_budget};
    iso.rec(b, -1.0, 1.0, 0);
    return iso.st;
}

}  // namespace ChebyshevBernstein
}  // namespace CoolProp
