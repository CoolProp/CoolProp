#include "CoolProp/region/CubicSplineCurve.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace CoolProp {
namespace region {

namespace {

// Solve a tridiagonal system A x = d where A has subdiagonal a, diagonal
// b, superdiagonal c.  Standard Thomas algorithm — O(n) with one pass
// forward, one pass back.  Writes the solution into d in-place.
void thomas(std::vector<double>& sub, std::vector<double>& diag, std::vector<double>& sup, std::vector<double>& d) {
    const std::size_t n = diag.size();
    // Forward sweep.
    for (std::size_t i = 1; i < n; ++i) {
        const double m = sub[i] / diag[i - 1];
        diag[i] -= m * sup[i - 1];
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
  : a_(std::move(a)), b_(std::move(b)), M_(std::move(M)), b_min_(b_min), b_max_(b_max) {}

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
    return std::unique_ptr<CubicSplineCurve>(new CubicSplineCurve(std::move(s.a), std::move(s.b), std::move(s.M), s.b_min, s.b_max));
}

std::size_t CubicSplineCurve::locate(double a) const noexcept {
    // Binary search; clamp to a valid segment index.
    const std::size_t n = a_.size();
    if (a <= a_.front()) {
        return 0;
    }
    if (a >= a_.back()) {
        return n - 2;
    }
    std::size_t lo = 0;
    std::size_t hi = n - 1;
    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;
        if (a_[mid] <= a) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
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
