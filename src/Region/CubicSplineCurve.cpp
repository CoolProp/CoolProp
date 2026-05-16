#include "CoolProp/region/CubicSplineCurve.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>  // std::to_string + std::invalid_argument(string) — MSVC needs explicit include

namespace CoolProp {
namespace region {

namespace {

// Solve a tridiagonal system A x = d where A has subdiagonal a, diagonal
// b, superdiagonal c.  Standard Thomas algorithm — O(n) with one pass
// forward, one pass back.  Writes the solution into d in-place.
void thomas(std::vector<double>& sub, std::vector<double>& diag, std::vector<double>& sup, std::vector<double>& d) {
    const std::size_t n = diag.size();
    // Defense in depth: every CubicSplineCurve build() call feeds the
    // natural-cubic system whose diagonal is (h_lo + h_hi)/3 with both
    // h_i strictly positive (we validate strictly-increasing a above),
    // so the diagonal is always positive.  But if someone routes a
    // future call into this helper with a pathological matrix, the
    // forward sweep silently divides by zero -- check explicitly so we
    // surface a clear error instead of producing NaN coefficients.
    constexpr double kTinyDiag = 1e-300;
    for (std::size_t i = 0; i < n; ++i) {
        if (std::abs(diag[i]) < kTinyDiag) {
            throw std::invalid_argument("thomas: zero (or sub-normal) diagonal entry; tridiagonal system is singular");
        }
    }
    // Forward sweep.
    for (std::size_t i = 1; i < n; ++i) {
        const double m = sub[i] / diag[i - 1];
        diag[i] -= m * sup[i - 1];
        if (std::abs(diag[i]) < kTinyDiag) {
            throw std::invalid_argument("thomas: pivot vanished during forward sweep (system is singular at row " + std::to_string(i) + ")");
        }
        d[i] -= m * d[i - 1];
    }
    // Back substitution.
    d[n - 1] /= diag[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) {
        d[i] = (d[i] - sup[i] * d[i + 1]) / diag[i];
    }
}

// Evaluate the cubic on segment [a_i, a_{i+1}] at parameter t in [0,1]
// given knot values b_i, b_{i+1} and second derivatives M_i, M_{i+1}
// and segment width h = a_{i+1} - a_i.
//
//   p(t) = (1 - t) b_i + t b_{i+1}
//        + (h^2 / 6) [ ((1-t)^3 - (1-t)) M_i + (t^3 - t) M_{i+1} ]
double segment_eval(double t, double h, double b_i, double b_ip1, double M_i, double M_ip1) noexcept {
    const double s = 1.0 - t;
    const double h2_6 = h * h / 6.0;
    return s * b_i + t * b_ip1 + h2_6 * ((s * s * s - s) * M_i + (t * t * t - t) * M_ip1);
}

// dp/da on segment [a_i, a_{i+1}] at parameter t in [0,1]:
//
//   dp/dt = (b_{i+1} - b_i) + (h^2 / 6) [ (-3(1-t)^2 + 1) M_i
//                                          + ( 3 t^2  - 1) M_{i+1} ]
//   dp/da = dp/dt / h
//         = (b_{i+1} - b_i)/h + (h/6) [ (1 - 3(1-t)^2) M_i
//                                        + (3 t^2 - 1) M_{i+1} ]
double segment_eval_da(double t, double h, double b_i, double b_ip1, double M_i, double M_ip1) noexcept {
    const double s = 1.0 - t;
    const double h_6 = h / 6.0;
    return (b_ip1 - b_i) / h + h_6 * ((1.0 - 3.0 * s * s) * M_i + (3.0 * t * t - 1.0) * M_ip1);
}

// Tight min/max of the spline on its full build interval.
//
// Sweep every segment, locate the (up to two) critical points of the
// cubic from dp/dt = 0, and track the extremes against the knot values.
// dp/dt is a quadratic in t — root-find with the standard formula,
// clamped to [0, 1].
std::pair<double, double> compute_bounds(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& M) {
    double lo = b.front();
    double hi = b.front();
    for (const double v : b) {
        if (v < lo) {
            lo = v;
        }
        if (v > hi) {
            hi = v;
        }
    }
    const std::size_t n_seg = a.size() - 1;
    for (std::size_t i = 0; i < n_seg; ++i) {
        const double h = a[i + 1] - a[i];
        // dp/dt = (b_{i+1}-b_i) + (h^2/6)[(-3(1-t)^2+1) M_i + (3 t^2 - 1) M_{i+1}]
        // Expand to A t^2 + B t + C = 0:
        //   coefficient of t^2 : (h^2/6) * 3 * (M_{i+1} - M_i)
        //   coefficient of t   : (h^2/6) * 6 * M_i
        //   constant           : (b_{i+1}-b_i) + (h^2/6)*(-3 M_i*1) + (h^2/6)*(-M_{i+1})
        //                      = (b_{i+1}-b_i) - (h^2/6)*(3 M_i + M_{i+1})
        // (Algebra: -3(1-t)^2 + 1 = -3 + 6 t - 3 t^2 + 1 = -2 + 6 t - 3 t^2.)
        // Re-derive cleanly:
        //   -3(1-t)^2 + 1 = -3 + 6 t - 3 t^2 + 1 = -3 t^2 + 6 t - 2
        //   3 t^2 - 1
        // So dp/dt = D + (h^2/6) [(-3 t^2 + 6 t - 2) M_i + (3 t^2 - 1) M_{i+1}]
        //         = D + (h^2/6) [ 3 t^2 (M_{i+1} - M_i)
        //                        + 6 t  M_i
        //                        + (-2 M_i - M_{i+1}) ]
        // where D = b_{i+1} - b_i.
        const double D = b[i + 1] - b[i];
        const double k = h * h / 6.0;
        const double A = k * 3.0 * (M[i + 1] - M[i]);
        const double B = k * 6.0 * M[i];
        const double C = D + k * (-2.0 * M[i] - M[i + 1]);
        const auto try_t = [&](double t) {
            if (t > 0.0 && t < 1.0) {
                const double v = segment_eval(t, h, b[i], b[i + 1], M[i], M[i + 1]);
                if (v < lo) {
                    lo = v;
                }
                if (v > hi) {
                    hi = v;
                }
            }
        };
        if (std::abs(A) < 1e-300) {
            // Linear (degenerate quadratic): B t + C = 0.
            if (std::abs(B) > 1e-300) {
                try_t(-C / B);
            }
        } else {
            const double disc = B * B - 4.0 * A * C;
            if (disc >= 0.0) {
                const double sq = std::sqrt(disc);
                try_t((-B + sq) / (2.0 * A));
                try_t((-B - sq) / (2.0 * A));
            }
        }
    }
    return {lo, hi};
}

}  // namespace

std::unique_ptr<CubicSplineCurve> CubicSplineCurve::build(std::vector<double> a, std::vector<double> b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("CubicSplineCurve::build: a and b must have the same size");
    }
    const std::size_t n = a.size();
    if (n < 2) {
        throw std::invalid_argument("CubicSplineCurve::build: need at least 2 knots");
    }
    for (std::size_t i = 1; i < n; ++i) {
        if (!(a[i] > a[i - 1])) {
            throw std::invalid_argument("CubicSplineCurve::build: a must be strictly increasing");
        }
    }

    std::vector<double> M(n, 0.0);
    if (n >= 3) {
        // Build tridiagonal system for interior unknowns M[1..n-2] with
        // natural BCs M[0] = M[n-1] = 0.
        const std::size_t m = n - 2;
        std::vector<double> sub(m, 0.0);
        std::vector<double> diag(m, 0.0);
        std::vector<double> sup(m, 0.0);
        std::vector<double> rhs(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            const std::size_t k = i + 1;  // global index
            const double h_lo = a[k] - a[k - 1];
            const double h_hi = a[k + 1] - a[k];
            sub[i] = h_lo / 6.0;
            diag[i] = (h_lo + h_hi) / 3.0;
            sup[i] = h_hi / 6.0;
            rhs[i] = (b[k + 1] - b[k]) / h_hi - (b[k] - b[k - 1]) / h_lo;
        }
        // Natural BC eliminates the off-diagonal terms touching M[0],
        // M[n-1]: their coefficients in the rhs are zero, so sub[0] and
        // sup[m-1] feed into nothing (they multiply unknowns at i=-1 and
        // i=m, which don't exist in the reduced system).
        sub[0] = 0.0;
        sup[m - 1] = 0.0;
        thomas(sub, diag, sup, rhs);
        for (std::size_t i = 0; i < m; ++i) {
            M[i + 1] = rhs[i];
        }
    }
    // n == 2: M stays zero, the spline reduces to the straight line.

    auto bnds = compute_bounds(a, b, M);
    return std::unique_ptr<CubicSplineCurve>(new CubicSplineCurve(std::move(a), std::move(b), std::move(M), bnds.first, bnds.second));
}

CubicSplineCurve::CubicSplineCurve(std::vector<double> a, std::vector<double> b, std::vector<double> M, double b_min, double b_max)
  : a_(std::move(a)), b_(std::move(b)), M_(std::move(M)), b_min_(b_min), b_max_(b_max) {
    build_bucket_table_();
}

void CubicSplineCurve::build_bucket_table_() noexcept {
    const std::size_t n = a_.size();
    const double a_lo = a_.front();
    const double a_hi = a_.back();
    const double span = a_hi - a_lo;
    // Guard against degenerate single-segment curves where span could
    // be zero (build() already rejects this, but be defensive).
    a_lo_inv_step_ = (span > 0.0) ? (static_cast<double>(kBuckets) / span) : 0.0;
    // Sweep buckets left-to-right; maintain a single advancing knot
    // index `i` so the build is O(kBuckets + n) total.
    std::size_t i = 0;
    const auto max_seg_idx = static_cast<std::uint16_t>(n - 2);
    for (std::size_t k = 0; k < kBuckets; ++k) {
        const double a_bucket = a_lo + static_cast<double>(k) * span / static_cast<double>(kBuckets);
        // Advance i while a_[i+1] is still <= a_bucket.  Loop invariant
        // after this step: a_[i] <= a_bucket < a_[i+1] (modulo clamping).
        while (i + 1 < n - 1 && a_[i + 1] <= a_bucket) {
            ++i;
        }
        bucket_to_knot_[k] = (i > max_seg_idx) ? max_seg_idx : static_cast<std::uint16_t>(i);
    }
}

CubicSplineCurve::State CubicSplineCurve::state() const {
    return State{a_, b_, M_, b_min_, b_max_};
}

std::unique_ptr<CubicSplineCurve> CubicSplineCurve::from_state(State s) {
    const std::size_t n = s.a.size();
    if (s.b.size() != n || s.M.size() != n) {
        throw std::invalid_argument("CubicSplineCurve::from_state: a/b/M size mismatch");
    }
    if (n < 2) {
        throw std::invalid_argument("CubicSplineCurve::from_state: need at least 2 knots");
    }
    for (std::size_t i = 1; i < n; ++i) {
        if (!(s.a[i] > s.a[i - 1])) {
            throw std::invalid_argument("CubicSplineCurve::from_state: a must be strictly increasing");
        }
    }
    // Re-derive the tight (b_min, b_max) from the deserialized (a, b, M)
    // rather than trusting the values shipped in `s`.  A tampered or
    // mis-versioned blob could otherwise drop a bogus AABB into the
    // RegionAtlas and silently corrupt region dispatch.  The recompute
    // is O(n) — negligible vs the rest of deserialization.  The private
    // ctor rebuilds the bucket table from the knot array as part of
    // construction, so no extra work needed here.
    const auto bnds = compute_bounds(s.a, s.b, s.M);
    return std::unique_ptr<CubicSplineCurve>(new CubicSplineCurve(std::move(s.a), std::move(s.b), std::move(s.M), bnds.first, bnds.second));
}

std::size_t CubicSplineCurve::locate(double a) const noexcept {
    // O(1) amortized: hash into a precomputed bucket, then walk
    // forward at most a few knots.  Non-uniform knot spacing is fine
    // — the build_bucket_table_ pass already absorbed it.
    const std::size_t n = a_.size();
    const double a_lo = a_.front();
    // NaN/inf guard: the comparisons below are both false on NaN,
    // and the cast to ptrdiff_t a few lines down is UB on non-finite
    // values.  Return the leftmost segment for any non-finite input
    // — the eventual eval() math will propagate NaN to the caller
    // through the standard return path.
    if (!std::isfinite(a)) {
        return 0;
    }
    if (a <= a_lo) {
        return 0;
    }
    if (a >= a_.back()) {
        return n - 2;
    }
    // Bucket index in [0, kBuckets).  Clamp defensively in case
    // a_lo_inv_step_ is 0 or rounding pushes us out.
    auto k = static_cast<std::ptrdiff_t>((a - a_lo) * a_lo_inv_step_);
    if (k < 0) {
        k = 0;
    } else if (k >= static_cast<std::ptrdiff_t>(kBuckets)) {
        k = static_cast<std::ptrdiff_t>(kBuckets) - 1;
    }
    std::size_t i = bucket_to_knot_[static_cast<std::size_t>(k)];
    // Refine forward.  Loop bound matches the worst-case "knots per
    // bucket" which is ~ceil(n / kBuckets) for any reasonable knot
    // distribution (1 for n <= kBuckets, larger for dense knot tails
    // on log-uniform sampling).  Practical iterations: 0–1.
    while (i + 1 < n - 1 && a_[i + 1] <= a) {
        ++i;
    }
    return i;
}

double CubicSplineCurve::eval(double a) const noexcept {
    const std::size_t i = locate(a);
    const double h = a_[i + 1] - a_[i];
    const double t = (a - a_[i]) / h;
    return segment_eval(t, h, b_[i], b_[i + 1], M_[i], M_[i + 1]);
}

double CubicSplineCurve::eval_da(double a) const noexcept {
    const std::size_t i = locate(a);
    const double h = a_[i + 1] - a_[i];
    const double t = (a - a_[i]) / h;
    return segment_eval_da(t, h, b_[i], b_[i + 1], M_[i], M_[i + 1]);
}

std::pair<double, double> CubicSplineCurve::bounds() const noexcept {
    return {b_min_, b_max_};
}

std::pair<double, double> CubicSplineCurve::a_range() const noexcept {
    return {a_.front(), a_.back()};
}

}  // namespace region
}  // namespace CoolProp
