#include "CoolProp/region/PiecewiseChebyshevCurve.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace CoolProp {
namespace region {

namespace {

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
    const double pi_N = M_PI / static_cast<double>(N);
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
    for (std::size_t k = N; k >= 1; --k) {
        const double b_k = 2.0 * s * b_kp1 - b_kp2 + c[k];
        b_kp2 = b_kp1;
        b_kp1 = b_k;
        if (k == 1) {
            break;
        }  // unsigned guard against k-- underflow
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
// Note this gives df/ds, not df/da; the chain rule via dt/da and ds/dt
// is handled at the call site.
std::vector<double> chebyshev_derivative_coeffs(const std::vector<double>& c) {
    const std::size_t Np1 = c.size();
    if (Np1 < 2) {
        return std::vector<double>(Np1, 0.0);
    }
    const std::size_t N = Np1 - 1;
    std::vector<double> d(Np1, 0.0);  // d[N] = 0
    for (std::size_t k = N; k >= 1; --k) {
        const double d_km1 = (k + 1 <= N ? d[k + 1] : 0.0) + 2.0 * static_cast<double>(k) * c[k];
        d[k - 1] = d_km1;
        if (k == 1) {
            break;
        }
    }
    d[0] *= 0.5;
    // The N-th coefficient of the derivative is zero (degree drops by 1).
    d[N] = 0.0;
    return d;
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
    const double pi_over_deg = M_PI / static_cast<double>(degree);

    for (std::size_t p = 0; p < n_pieces; ++p) {
        Piece piece;
        piece.t_lo = t_lo + static_cast<double>(p) * dt;
        piece.t_hi = piece.t_lo + dt;
        piece.t_mid = 0.5 * (piece.t_lo + piece.t_hi);
        piece.inv_half_span = 2.0 / (piece.t_hi - piece.t_lo);
        // Sample f at Lobatto nodes s_j = cos(j*pi/N), j = 0..N.
        std::vector<double> fj(Np1);
        for (std::size_t j = 0; j <= degree; ++j) {
            const double s_j = std::cos(static_cast<double>(j) * pi_over_deg);
            const double t_j = piece.t_mid + 0.5 * dt * s_j;
            const double a_j = (scale == ParamScale::LOG) ? std::exp(t_j) : t_j;
            const double v = f(a_j);
            fj[j] = v;
            if (v < b_min) {
                b_min = v;
            }
            if (v > b_max) {
                b_max = v;
            }
        }
        piece.coeffs = chebyshev_coeffs_from_samples(fj);
        pieces.push_back(std::move(piece));
    }

    return std::unique_ptr<PiecewiseChebyshevCurve>(new PiecewiseChebyshevCurve(a_lo, a_hi, scale, std::move(pieces), b_min, b_max));
}

PiecewiseChebyshevCurve::PiecewiseChebyshevCurve(double a_lo, double a_hi, ParamScale scale, std::vector<Piece> pieces, double b_min, double b_max)
  : a_lo_(a_lo), a_hi_(a_hi), scale_(scale), pieces_(std::move(pieces)), b_min_(b_min), b_max_(b_max) {}

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
    // df/ds via the derivative coefficient recurrence, then chain rule.
    const std::vector<double> d_coeffs = chebyshev_derivative_coeffs(pc.coeffs);
    const double df_ds = clenshaw_eval(d_coeffs, s);
    // ds/dt = inv_half_span; dt/da from to_t.
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
