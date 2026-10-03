#pragma once
/**
 * Certified real-root isolation for Chebyshev series on [-1, 1].
 *
 * A degree-N Chebyshev series f(u) = c_0 T_0(u) + ... + c_N T_N(u) is converted to the Bernstein basis on
 * [-1, 1] with an exact (rational, rounded once) conversion matrix.  Each interval is then tested on its
 * Bernstein coefficients b_j, whose signs are only trusted where |b_j| exceeds a tolerance: the caller's
 * tolerance plus a rigorous bound on the rounding error accumulated in the conversion and the subdivisions.
 *
 *  - convex hull: all b_j clearly of one sign -> f has no root on the interval;
 *  - Descartes' rule in the Bernstein basis: the number of sign changes V of (b_0, ..., b_N) bounds the number
 *    of roots and has the same parity; V = 0 -> no root;
 *  - V = 1 with every sign clear -> exactly one root, certified, refined by safeguarded Newton;
 *  - otherwise the interval is split at its midpoint by de Casteljau, which yields the Bernstein coefficients
 *    of both halves directly, and each half is tested again.
 *
 * Recursion stops at a depth cap (lower when a coefficient's sign is unknown, since subdividing cannot resolve
 * that) or when a node budget runs out.  Such an interval is never dropped: it is reported as uncertified.
 * So every root of f in [-1, 1] lies in a reported interval, certified or not, and discarding uncertified
 * intervals is the caller's decision, made explicitly.
 *
 * With a tolerance tol > 0, the same statements hold for any function g with |g - f| <= tol on [-1, 1]: every
 * root of g lies in a reported interval, and across a certified interval g changes sign (at least one root).
 * Exactly one root is certified for f itself only -- g can cross several times where |f| < tol.
 *
 * Costs: O(N^2) for the conversion and for each split.
 */

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>
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
    if constexpr (N >= 1) {
        d[N - 1] = 2 * N * c[N];
    }
    if constexpr (N >= 2) {
        d[N - 2] = 2 * (N - 1) * c[N - 1];
        for (int k = N - 3; k >= 0; --k)
            d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1];
    }
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
/// each one: |b_j - exact_j| <= err for every j.  The bound covers the once-rounded matrix entries, the
/// products and the sums (each at most a unit roundoff, eps/2, relative), with a factor 2 to spare, plus an
/// absolute floor for coefficients in the subnormal range, where rounding is absolute rather than relative.
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
        err = std::max(err, (N + 2) * DBL_EPSILON * sa);
    }
    err += (N + 2) * std::numeric_limits<double>::denorm_min();
    return b;
}

/// Safeguarded Newton on the Chebyshev series inside a bracket [a, b] across which it changes sign;
/// negative_at_a says which side is which.  Every step that leaves the (shrinking) bracket is replaced by
/// bisection, so the result is always in [a, b].
template <int N>
double refine(const Coeffs<N>& c, double a, double b, bool negative_at_a, double xtol) {
    double x = 0.5 * (a + b);
    for (int it = 0; it < 100; ++it) {
        double f, df;
        clenshaw_fd<N>(c, x, f, df);
        if (f == 0) return x;
        if ((f < 0) == negative_at_a)
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
    double u;          ///< refined root when sign_change is true; the interval midpoint otherwise
    double ua, ub;     ///< the interval known to contain it
    bool certified;    ///< exactly one simple root of the series in [ua, ub] (given the tolerance, see real_roots)
    bool sign_change;  ///< the series has clearly opposite signs at ua and ub (so at least one root); always true
                       ///< when certified.  False means no such guarantee: an even number of roots (e.g. a double
                       ///< root), a root at an end of the interval (also at u = -1 or 1), or signs within tol
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
    long unresolved = 0;            ///< entries reported uncertified (after merging touching intervals)
    bool budget_exhausted = false;  ///< the node budget ran out (some intervals were reported unresolved for that reason)
};

namespace detail {
/// Per-call state of the isolation (values only; the series and the output are passed alongside)
struct IsolationState
{
    double tol = 0;  // caller tolerance (validated)
    Options opt;
    Stats st;
    long budget = 0;
    // the left-end Bernstein coefficient and tolerance of the last reported unresolved interval, to recompute
    // the sign-change flag when a touching interval is merged into it
    double last_b0 = 0, last_tol0 = 0;
};

/// Report [ua, ub] as unresolved, given its Bernstein coefficients b and their current tolerance
template <int N>
void report_unresolved(const Coeffs<N>& c, const Coeffs<N>& b, double tol, IsolationState& s, std::vector<Root>& out, double ua, double ub) {
    const bool merge = !out.empty() && !out.back().certified && out.back().ub == ua;  // touching (split point seen from both sides)
    const double b0 = merge ? s.last_b0 : b[0], tol0 = merge ? s.last_tol0 : tol, bN = b[N];
    const double a = merge ? out.back().ua : ua;
    // sign change only where both end values are clearly signed: the end Bernstein coefficients ARE the series
    // values there (up to the rounding the tolerance covers)
    const bool sc = ((b0 > tol0 && bN < -tol) || (b0 < -tol0 && bN > tol));
    const double u = sc ? refine<N>(c, a, ub, b0 < 0, s.opt.xtol) : 0.5 * (a + ub);
    if (merge) {
        Root& r = out.back();
        r.ub = ub;
        r.sign_change = sc;
        r.u = u;
        return;
    }
    ++s.st.unresolved;
    s.last_b0 = b0;
    s.last_tol0 = tol0;
    out.push_back({u, ua, ub, false, sc});
}

/// Isolate the roots on [ua, ub], whose Bernstein coefficients b carry an absolute rounding error <= e
template <int N>
void isolate(const Coeffs<N>& c, const Coeffs<N>& b, double ua, double ub, int depth, double e, IsolationState& s, std::vector<Root>& out) {
    ++s.st.nodes;
    // rounded up, so that e cannot be absorbed when it is below half an ulp of s.tol
    const double tol = (s.tol + e) * (1 + DBL_EPSILON);
    bool amb = false, anypos = false, anyneg = false;
    int V = 0, last = 0;
    double bmax = 0;
    for (int j = 0; j <= N; ++j) {
        bmax = std::max(bmax, std::abs(b[j]));
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
    if (!amb && V == 1) {                     // Descartes: exactly one root (the ends then have opposite signs)
        out.push_back({refine<N>(c, ua, ub, b[0] < 0, s.opt.xtol), ua, ub, true, true});
        return;
    }
    if ((amb && depth >= s.opt.max_depth_ambiguous) || depth >= s.opt.max_depth) {
        report_unresolved<N>(c, b, tol, s, out, ua, ub);
        return;
    }
    if (--s.budget < 0) {
        s.st.budget_exhausted = true;
        report_unresolved<N>(c, b, tol, s, out, ua, ub);
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
    // Each of the N rounds forms 0.5*(x + y) of values bounded by bmax + e: the sum rounds by at most
    // eps/2 * 2 (bmax + e), halved exactly -- or, in the subnormal range where rounding is absolute, by at most
    // denorm_min/2.  Averaging does not amplify the errors carried in, so the children carry
    // e + N * (eps/2 * (bmax + e) + denorm_min/2).
    // (N * denorm_min rather than N/2: an exact multiple, so it cannot round down when the eps term underflows)
    const double e_child = e + N * 0.5 * DBL_EPSILON * (bmax + e) + N * std::numeric_limits<double>::denorm_min();
    const double um = 0.5 * (ua + ub);
    isolate<N>(c, L, ua, um, depth + 1, e_child, s, out);
    isolate<N>(c, R, um, ub, depth + 1, e_child, s, out);
}
}  // namespace detail

/**
 * All real roots of the degree-N Chebyshev series c on [-1, 1], in increasing order.
 *
 * @param c    Chebyshev coefficients
 * @param tol  absolute tolerance on the series values: Bernstein coefficients within tol (plus the rounding
 *             bound) of zero have unknown sign.  Pass a bound on |g - series| to make the coverage guarantee
 *             hold for a function g that the series approximates; 0 otherwise.  A negative tol is treated as 0;
 *             a NaN or infinite tol, like a non-finite coefficient, reports [-1, 1] as one unresolved entry.
 * @param out  receives the roots / root intervals (cleared first)
 * @return     the statistics of the call
 *
 * Guarantee: every root of the series in [-1, 1] -- and, for tol > 0, of any g with |g - series| <= tol --
 * lies in [ua, ub] of some reported entry.  A certified entry contains exactly one simple root of the series,
 * and g changes sign across it.
 */
template <int N>
Stats real_roots(const Coeffs<N>& c, double tol, std::vector<Root>& out, const Options& opt = Options()) {
    static_assert(N >= 1 && N <= MAX_DEGREE, "Chebyshev-to-Bernstein matrices are tabulated for degrees 1..17");
    out.clear();
    detail::IsolationState s;
    s.opt = opt;
    s.budget = opt.node_budget;
    bool finite = std::isfinite(tol);
    for (double v : c)
        finite = finite && std::isfinite(v);
    if (!finite) {  // nothing can be decided: fail closed, without spending the budget on NaN comparisons
        s.st.unresolved = 1;
        out.push_back({0.0, -1.0, 1.0, false, false});
        return s.st;
    }
    s.tol = std::max(tol, 0.0);
    double err = 0;
    const Coeffs<N> b = to_bernstein<N>(c, err);
    detail::isolate<N>(c, b, -1.0, 1.0, 0, err, s, out);
    return s.st;
}

}  // namespace ChebyshevBernstein
}  // namespace CoolProp
