#include "CoolProp/region/PiecewiseChebyshevCurve.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace CoolProp {
namespace region {

namespace {

// Local portability constant.  M_PI is non-standard and MSVC does not
// define it without _USE_MATH_DEFINES; rather than reach for the
// project-wide shim in CPnumerics.h (which pulls in unrelated math
// utilities), keep this file self-contained.
constexpr double kPi = 3.141592653589793238462643383279502884;

// Compute Chebyshev coefficients c_0..c_N of f(s), s in [-1, 1], from
// samples taken at the Chebyshev-Lobatto nodes s_j = cos(j*pi/N),
// j = 0..N.  Returns N+1 coefficients such that
//
//   f(s) ~= sum_{k=0}^{N} c_k T_k(s)
//
// using the DCT-I formula
//
//   c_k = (2/N) * w_k * sum_{j=0..N} w_j * f(s_j) * cos(k j pi / N)
//
// where w_j = 1/2 for j in {0, N}, else 1.
std::vector<double> chebyshev_coeffs_from_samples(const std::vector<double>& fj) {
    const std::size_t Np1 = fj.size();
    if (Np1 < 2) {
        throw std::invalid_argument("chebyshev_coeffs_from_samples: need at least 2 samples");
    }
    const std::size_t N = Np1 - 1;
    std::vector<double> c(Np1, 0.0);
    const double pi_N = kPi / static_cast<double>(N);
    for (std::size_t k = 0; k <= N; ++k) {
        double s = 0.0;
        for (std::size_t j = 0; j <= N; ++j) {
            const double wj = (j == 0 || j == N) ? 0.5 : 1.0;
            s += wj * fj[j] * std::cos(static_cast<double>(k * j) * pi_N);
        }
        const double wk = (k == 0 || k == N) ? 0.5 : 1.0;
        c[k] = (2.0 / static_cast<double>(N)) * wk * s;
    }
    return c;
}

// Clenshaw recurrence: evaluate sum_{k=0}^{N} c_k T_k(s) at s.
double clenshaw_eval(const std::vector<double>& c, double s) noexcept {
    if (c.empty()) {
        return 0.0;
    }
    const std::size_t N = c.size() - 1;
    double b_kp1 = 0.0;
    double b_kp2 = 0.0;
    // Iterate k = N, N-1, ..., 1 (inclusive).  Using `k > 0` instead of
    // `k >= 1` for the condition is equivalent for unsigned k and lets
    // CodeQL see that the loop terminates -- the previous form (which
    // wrote `k >= 1` with an inner `break` for k == 1) tripped a "comparison
    // always true" alert because std::size_t is non-negative by
    // construction.
    for (std::size_t k = N; k > 0; --k) {
        const double b_k = 2.0 * s * b_kp1 - b_kp2 + c[k];
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    return s * b_kp1 - b_kp2 + c[0];
}

// Compute the coefficient vector of the *derivative* expansion.
//
// If f(s) = sum_{k=0}^{N} c_k T_k(s), then df/ds = sum_{k=0}^{N-1} c'_k T_k(s)
// with c' defined by the standard backward recurrence:
//
//   c'_N = 0
//   c'_{k-1} = c'_{k+1} + 2 k c_k     for k = N, N-1, ..., 1
//   c'_0 /= 2
//
// The returned vector has size N (NOT N+1) — the trivially-zero c'_N
// coefficient is dropped so Clenshaw doesn't iterate over it.
//
// Note this gives df/ds, not df/da; the chain rule via dt/da and ds/dt
// is handled at the call site.
std::vector<double> chebyshev_derivative_coeffs(const std::vector<double>& c) {
    const std::size_t Np1 = c.size();
    if (Np1 < 2) {
        return {};
    }
    const std::size_t N = Np1 - 1;
    // d has size N (degree-1 expansion); indices 0..N-1.
    std::vector<double> d(N, 0.0);
    // Backward recurrence: c'_{k-1} = c'_{k+1} + 2k c_k, k = N..1, with
    // c'_N implicitly zero (we never write to index N).
    for (std::size_t kp = N; kp > 0; --kp) {
        const double d_km1 = (kp + 1 < N ? d[kp + 1] : 0.0) + 2.0 * static_cast<double>(kp) * c[kp];
        d[kp - 1] = d_km1;
    }
    d[0] *= 0.5;
    return d;
}

// Tight extremum of the piece f(s) = sum c_k T_k(s) on s in [-1, 1].
// Drives the AABB so the cheap RegionAtlas filter rejects probes
// accurately — the conservative `c_0 ± sum |c_k|` bound can inflate
// the AABB 30–50% beyond the true range and defeat the filter.
//
// Strategy: f is smooth, so its extrema on [-1,1] are at s = ±1 or at
// the roots of f'(s) = 0 in (-1, 1).  Evaluate f at both endpoints,
// then root-find f' in (-1, 1) and evaluate f at each root.  Track
// running min/max.
//
// f' has degree N-1, with N+1 = c.size().  We sample f' on a fine
// grid (8 * N points), look for sign changes, and bisect each
// bracket to ~1e-10.  Cheap and robust for the modest degrees we use
// (build-time cost only).
std::pair<double, double> tight_piece_bounds(const std::vector<double>& c, const std::vector<double>& d) noexcept {
    auto eval_s = [&](double s) { return clenshaw_eval(c, s); };
    auto eval_deriv = [&](double s) { return clenshaw_eval(d, s); };

    double lo = eval_s(-1.0);
    double hi = lo;
    const double at_p1 = eval_s(1.0);
    if (at_p1 < lo) {
        lo = at_p1;
    }
    if (at_p1 > hi) {
        hi = at_p1;
    }
    // 8x oversample on derivative — sufficient to bracket all roots
    // of a polynomial of degree (N-1) where N is the curve's degree.
    const std::size_t N = c.size() - 1;
    const std::size_t n_samples = 8 * std::max<std::size_t>(N, 4);
    double prev_s = -1.0;
    double prev_d = eval_deriv(prev_s);
    for (std::size_t k = 1; k <= n_samples; ++k) {
        const double cur_s = -1.0 + 2.0 * static_cast<double>(k) / static_cast<double>(n_samples);
        const double cur_d = eval_deriv(cur_s);
        if (prev_d == 0.0 || (prev_d * cur_d < 0.0)) {
            // Bracket [prev_s, cur_s]; bisect.
            double a_lo = prev_s, a_hi = cur_s;
            double f_lo = prev_d, f_hi = cur_d;
            for (int it = 0; it < 60; ++it) {
                const double mid = 0.5 * (a_lo + a_hi);
                const double f_mid = eval_deriv(mid);
                if (f_lo * f_mid <= 0.0) {
                    a_hi = mid;
                    f_hi = f_mid;
                } else {
                    a_lo = mid;
                    f_lo = f_mid;
                }
                if (a_hi - a_lo < 1e-12) {
                    break;
                }
            }
            const double root_val = eval_s(0.5 * (a_lo + a_hi));
            if (root_val < lo) {
                lo = root_val;
            }
            if (root_val > hi) {
                hi = root_val;
            }
        }
        prev_s = cur_s;
        prev_d = cur_d;
    }
    return {lo, hi};
}

}  // namespace

std::unique_ptr<PiecewiseChebyshevCurve> PiecewiseChebyshevCurve::build(double a_lo, double a_hi, std::size_t n_pieces, std::size_t degree,
                                                                        ParamScale scale, const std::function<double(double)>& f) {
    if (!(a_hi > a_lo)) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::build: a_hi must exceed a_lo");
    }
    if (scale == ParamScale::LOG && !(a_lo > 0.0)) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::build: LOG requires a_lo > 0");
    }
    if (n_pieces < 1) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::build: need at least 1 piece");
    }
    if (degree < 2) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::build: degree must be >= 2");
    }

    const double t_lo = (scale == ParamScale::LOG) ? std::log(a_lo) : a_lo;
    const double t_hi = (scale == ParamScale::LOG) ? std::log(a_hi) : a_hi;
    const double dt = (t_hi - t_lo) / static_cast<double>(n_pieces);

    std::vector<Piece> pieces;
    pieces.reserve(n_pieces);
    double b_min = std::numeric_limits<double>::infinity();
    double b_max = -std::numeric_limits<double>::infinity();
    const std::size_t Np1 = degree + 1;
    const double pi_over_deg = kPi / static_cast<double>(degree);

    for (std::size_t p = 0; p < n_pieces; ++p) {
        Piece piece;
        piece.t_lo = t_lo + static_cast<double>(p) * dt;
        piece.t_hi = piece.t_lo + dt;
        piece.t_mid = 0.5 * (piece.t_lo + piece.t_hi);
        piece.inv_half_span = 2.0 / (piece.t_hi - piece.t_lo);
        // Sample f at Lobatto nodes s_j = cos(j*pi/N), j = 0..N.
        // Reject non-finite samples up front — without this check the
        // build silently produces a curve with NaN coefficients, and
        // every downstream eval / aabb_contains then returns NaN/false
        // and the atlas silently mis-dispatches.
        std::vector<double> fj(Np1);
        for (std::size_t j = 0; j <= degree; ++j) {
            const double s_j = std::cos(static_cast<double>(j) * pi_over_deg);
            const double t_j = piece.t_mid + 0.5 * dt * s_j;
            const double a_j = (scale == ParamScale::LOG) ? std::exp(t_j) : t_j;
            const double v = f(a_j);
            if (!std::isfinite(v)) {
                throw std::invalid_argument("PiecewiseChebyshevCurve::build: sample function returned non-finite value");
            }
            fj[j] = v;
        }
        piece.coeffs = chebyshev_coeffs_from_samples(fj);
        piece.deriv_coeffs = chebyshev_derivative_coeffs(piece.coeffs);
        // Tight bound on the piece via cubic-derivative root finding
        // (see tight_piece_bounds): a conservative |c_0| ± Σ|c_k| bound
        // is correct but pessimistic enough to inflate the AABB 30–50%
        // beyond the true range, which defeats the cheap RegionAtlas
        // filter.  Tight bounds keep curve_contains rare.
        const auto [piece_lo, piece_hi] = tight_piece_bounds(piece.coeffs, piece.deriv_coeffs);
        if (piece_lo < b_min) {
            b_min = piece_lo;
        }
        if (piece_hi > b_max) {
            b_max = piece_hi;
        }
        pieces.push_back(std::move(piece));
    }

    return std::unique_ptr<PiecewiseChebyshevCurve>(new PiecewiseChebyshevCurve(a_lo, a_hi, scale, std::move(pieces), b_min, b_max));
}

PiecewiseChebyshevCurve::PiecewiseChebyshevCurve(double a_lo, double a_hi, ParamScale scale, std::vector<Piece> pieces, double b_min, double b_max)
  : a_lo_(a_lo), a_hi_(a_hi), scale_(scale), pieces_(std::move(pieces)), b_min_(b_min), b_max_(b_max) {}

PiecewiseChebyshevCurve::State PiecewiseChebyshevCurve::state() const {
    State s;
    s.a_lo = a_lo_;
    s.a_hi = a_hi_;
    s.scale = scale_;
    s.b_min = b_min_;
    s.b_max = b_max_;
    s.pieces.reserve(pieces_.size());
    for (const auto& p : pieces_) {
        s.pieces.push_back(PieceState{p.t_lo, p.t_hi, p.inv_half_span, p.t_mid, p.coeffs, p.deriv_coeffs});
    }
    return s;
}

std::unique_ptr<PiecewiseChebyshevCurve> PiecewiseChebyshevCurve::from_state(State s) {
    if (s.pieces.empty()) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: pieces is empty");
    }
    if (!(s.a_hi > s.a_lo)) {
        throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: a_hi must exceed a_lo");
    }
    std::vector<Piece> pieces;
    pieces.reserve(s.pieces.size());
    for (auto& ps : s.pieces) {
        if (ps.coeffs.size() < 2) {
            throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: each piece needs >= 2 coefficients");
        }
        if (ps.deriv_coeffs.size() + 1 != ps.coeffs.size()) {
            throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: deriv_coeffs size must be coeffs.size() - 1");
        }
        if (!(ps.t_hi > ps.t_lo)) {
            throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: each piece must have t_hi > t_lo");
        }
        Piece p;
        p.t_lo = ps.t_lo;
        p.t_hi = ps.t_hi;
        p.inv_half_span = ps.inv_half_span;
        p.t_mid = ps.t_mid;
        p.coeffs = std::move(ps.coeffs);
        p.deriv_coeffs = std::move(ps.deriv_coeffs);
        pieces.push_back(std::move(p));
    }
    // `locate_piece` does an O(1) hit assuming pieces are sorted,
    // contiguous, and uniformly sized in t -- invariants build()
    // guarantees, but a malformed blob could break any of them and
    // silently route queries to the wrong piece.  Validate up-front so
    // the O(1) lookup math is always safe.
    const double dt0 = pieces.front().t_hi - pieces.front().t_lo;
    const double eps = 1e-9 * std::max(1.0, std::abs(dt0));
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (i > 0 && std::abs(pieces[i].t_lo - pieces[i - 1].t_hi) > eps) {
            throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: pieces must be contiguous in t (gap or overlap detected)");
        }
        const double dt = pieces[i].t_hi - pieces[i].t_lo;
        if (std::abs(dt - dt0) > eps) {
            throw std::invalid_argument("PiecewiseChebyshevCurve::from_state: pieces must be uniformly spaced in t (locate_piece assumes this)");
        }
    }
    return std::unique_ptr<PiecewiseChebyshevCurve>(new PiecewiseChebyshevCurve(s.a_lo, s.a_hi, s.scale, std::move(pieces), s.b_min, s.b_max));
}

double PiecewiseChebyshevCurve::to_t(double a) const noexcept {
    return (scale_ == ParamScale::LOG) ? std::log(a) : a;
}

double PiecewiseChebyshevCurve::dt_da(double a) const noexcept {
    return (scale_ == ParamScale::LOG) ? 1.0 / a : 1.0;
}

std::size_t PiecewiseChebyshevCurve::locate_piece(double t) const noexcept {
    const std::size_t n = pieces_.size();
    if (n == 1) {
        return 0;
    }
    // Pieces are uniformly spaced in t (the build guarantees this), so
    // an O(1) hit is available.
    const double t0 = pieces_.front().t_lo;
    const double tn = pieces_.back().t_hi;
    if (t <= t0) {
        return 0;
    }
    if (t >= tn) {
        return n - 1;
    }
    const double dt = pieces_.front().t_hi - pieces_.front().t_lo;
    auto idx = static_cast<std::ptrdiff_t>((t - t0) / dt);
    if (idx < 0) {
        return 0;
    }
    if (static_cast<std::size_t>(idx) >= n) {
        return n - 1;
    }
    return static_cast<std::size_t>(idx);
}

double PiecewiseChebyshevCurve::eval(double a) const noexcept {
    const double t = to_t(a);
    const std::size_t p = locate_piece(t);
    const Piece& pc = pieces_[p];
    const double s = (t - pc.t_mid) * pc.inv_half_span;
    return clenshaw_eval(pc.coeffs, s);
}

double PiecewiseChebyshevCurve::eval_da(double a) const noexcept {
    const double t = to_t(a);
    const std::size_t p = locate_piece(t);
    const Piece& pc = pieces_[p];
    const double s = (t - pc.t_mid) * pc.inv_half_span;
    // df/ds via the precomputed derivative coefficients (built once at
    // construction; see Piece::deriv_coeffs).  Chain rule via ds/dt
    // and dt/da at the call site.
    const double df_ds = clenshaw_eval(pc.deriv_coeffs, s);
    return df_ds * pc.inv_half_span * dt_da(a);
}

std::pair<double, double> PiecewiseChebyshevCurve::bounds() const noexcept {
    return {b_min_, b_max_};
}

std::pair<double, double> PiecewiseChebyshevCurve::a_range() const noexcept {
    return {a_lo_, a_hi_};
}

}  // namespace region
}  // namespace CoolProp
