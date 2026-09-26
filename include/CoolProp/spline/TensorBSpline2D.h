#ifndef COOLPROP_SPLINE_TENSOR_BSPLINE_2D_H
#define COOLPROP_SPLINE_TENSOR_BSPLINE_2D_H

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "CoolProp/Exceptions.h"

namespace CoolProp {
namespace spline {

// Tensor-product B-spline surface on a knot grid.
//
// Clamped and unclamped knot vectors are both supported.  The valid
// domain is the true support interval [knots[order-1], knots[n]], which
// coincides with [knots.front(), knots.back()] only when clamped.
//
//   f(x, y) = sum_i sum_j  c_ij  B_i,kx(x)  B_j,ky(y)
//
// ORDER CONVENTION.  `order` here is the B-spline ORDER, i.e. one more
// than the polynomial degree, matching the MATLAB spline-toolbox B-form
// and the Bollengier et al. (2019) water coefficient file (order 6 =
// quintic).  Note scipy.interpolate.BSpline takes the DEGREE instead;
// the two differ by one.  Getting this wrong is the likeliest source of
// a silently wrong answer, so it is stated in both directions here and
// in dev/scripts/gen_tensor_bspline_ref.py.
//
// The number of coefficients along an axis is fixed by the knot vector:
//   n = knots.size() - order
// and `coefs` is row-major with shape (nx, ny).
class TensorBSpline2D
{
   public:
    // Largest supported order per axis.  Evaluation uses fixed-size
    // stack scratch, so this is a hard bound, enforced by the
    // constructor rather than trusted from the caller.
    static constexpr std::size_t kMaxOrder = 16;

    TensorBSpline2D(std::vector<double> knots_x, std::vector<double> knots_y, std::size_t order_x, std::size_t order_y, std::vector<double> coefs);

    // Mixed partial d^(dx+dy) f / dx^dx dy^dy at (x, y); (0, 0) is the
    // value itself.  Any order is accepted: orders above the polynomial
    // degree of an axis are identically zero, as they should be.
    [[nodiscard]] double eval(double x, double y, unsigned dx = 0, unsigned dy = 0) const;

   private:
    // Throws ValueError unless (order, knots) describe a usable axis:
    // order in [1, kMaxOrder], enough knots that n = size - order is at
    // least 1, and a finite non-decreasing knot sequence.
    static void validate_axis(const std::vector<double>& knots, std::size_t order, const char* axis);

    // Throws ValueError unless `v` is finite and inside [knots[p],
    // knots[n]] -- the actual support of a clamped spline, which is
    // narrower than the knot vector's own extent.
    static void check_in_domain(const std::vector<double>& knots, std::size_t order, std::size_t n, double v, const char* axis);

    // Index of the knot span containing `v`, clamped to [p, n-1] so that
    // the right-hand endpoint belongs to the last span rather than
    // falling off the end.
    //
    // PRECONDITION: check_in_domain() has passed for `v`.  The bisection
    // below does not terminate for out-of-domain input -- it does not
    // merely return a wrong span, it hangs -- so this must never be
    // reached without the guard in front of it.
    [[nodiscard]] static std::size_t find_span(const std::vector<double>& knots, std::size_t order, std::size_t n, double v);

    // The `order` non-zero basis functions on span `span`, written to
    // out[0 .. order-1] and corresponding to coefficients
    // span-p .. span  (p = order - 1).  Cox-de Boor recurrence in the
    // triangular form of Piegl & Tiller, The NURBS Book, Alg. A2.2.
    static void basis_funs(const std::vector<double>& knots, std::size_t order, std::size_t span, double v, double* out);

    // The `nd`-th derivatives of those same `order` basis functions,
    // written to out[0 .. order-1].  Piegl & Tiller, Alg. A2.3
    // (DersBasisFuns), keeping only the row we asked for.
    //
    // Signed arithmetic throughout: the book's recurrence indexes
    // r - k, which goes negative, and doing that in std::size_t would
    // wrap to a huge value and read out of bounds.
    static void basis_ders(const std::vector<double>& knots, std::size_t order, std::size_t span, double v, unsigned nd, double* out);

    std::vector<double> kx_;
    std::vector<double> ky_;
    std::size_t ox_;
    std::size_t oy_;
    std::size_t nx_;
    std::size_t ny_;
    std::vector<double> coefs_;  // row-major (nx_, ny_)
};

inline void TensorBSpline2D::validate_axis(const std::vector<double>& knots, std::size_t order, const char* axis) {
    if (order == 0) {
        throw ValueError(std::string("TensorBSpline2D: ") + axis + " order must be at least 1 (order = degree + 1)");
    }
    if (order > kMaxOrder) {
        // Evaluation uses fixed-size stack scratch dimensioned by
        // kMaxOrder; accepting a larger order would overrun it.  Checked
        // before 2 * order below, so that product cannot overflow.
        throw ValueError(std::string("TensorBSpline2D: ") + axis + " order " + std::to_string(order) + " exceeds kMaxOrder "
                         + std::to_string(kMaxOrder));
    }
    if (knots.size() < 2 * order) {
        // n = knots.size() - order must be at least `order`, not merely at
        // least 1.  When n < order, find_span's right-endpoint early return
        // yields a span below p and eval's `sx - px + a` wraps in
        // std::size_t -- an out-of-bounds coefficient read, not merely a
        // wrong number.
        throw ValueError(std::string("TensorBSpline2D: ") + axis + " needs at least " + std::to_string(2 * order) + " knots for order "
                         + std::to_string(order) + " (so that n >= order), got " + std::to_string(knots.size()));
    }
    for (std::size_t i = 0; i < knots.size(); ++i) {
        if (!std::isfinite(knots[i])) {
            throw ValueError(std::string("TensorBSpline2D: non-finite ") + axis + " knot at index " + std::to_string(i));
        }
        if (i > 0 && knots[i] < knots[i - 1]) {
            throw ValueError(std::string("TensorBSpline2D: ") + axis + " knots must be non-decreasing, but knot " + std::to_string(i) + " decreases");
        }
    }
    // The last span must be non-degenerate.
    //
    // This is a VALUE condition, and it is the one that matters: for any
    // span s reached by the bisection the Cox-de Boor denominator is
    // knots[s+r+1] - knots[s+r+1-j], whose index pair always straddles
    // s / s+1, so knots[s] < knots[s+1] makes a zero denominator
    // impossible at ANY multiplicity.  The only way to reach a zero is
    // find_span's `v >= knots[n]` early return handing back span n-1 when
    // knots[n-1] == knots[n].  Two sub-cases, both admitted by every size
    // and multiplicity check: the domain collapsing to a point
    // (knots[order-1] == knots[n]), and a proper domain whose declared
    // right endpoint is the degenerate span.  Both previously returned a
    // silent NaN.
    //
    // Implied by this: knots[order-1] <= knots[n-1] < knots[n], so the
    // domain is always a non-empty interval.
    //
    // An earlier version of this guard rejected knot runs longer than
    // `order` instead.  That was wrong in both directions -- it admitted
    // both sub-cases above (a multiplicity of 2 at index n is enough) and
    // it refused legitimate vectors such as order 3 {0,0,0,0,1,2,3,4},
    // which scipy accepts and evaluates.
    const std::size_t n = knots.size() - order;
    if (!(knots[n - 1] < knots[n])) {
        throw ValueError(std::string("TensorBSpline2D: ") + axis + " has a degenerate last span: knots[" + std::to_string(n - 1) + "] == knots["
                         + std::to_string(n) + "] == " + std::to_string(knots[n])
                         + ", which makes the basis undefined at the right end of the domain");
    }
}

inline TensorBSpline2D::TensorBSpline2D(std::vector<double> knots_x, std::vector<double> knots_y, std::size_t order_x, std::size_t order_y,
                                        std::vector<double> coefs)
  : kx_(std::move(knots_x)), ky_(std::move(knots_y)), ox_(order_x), oy_(order_y), nx_(0), ny_(0), coefs_(std::move(coefs)) {
    // Both axes are validated before nx_/ny_ are derived, because the
    // subtraction that derives them is unsigned and would wrap on a
    // short knot vector instead of reporting the real problem.
    validate_axis(kx_, ox_, "x");
    validate_axis(ky_, oy_, "y");
    nx_ = kx_.size() - ox_;
    ny_ = ky_.size() - oy_;
    if (coefs_.size() != nx_ * ny_) {
        throw ValueError("TensorBSpline2D: expected " + std::to_string(nx_) + " x " + std::to_string(ny_) + " = " + std::to_string(nx_ * ny_)
                         + " coefficients, got " + std::to_string(coefs_.size()));
    }
}

inline void TensorBSpline2D::check_in_domain(const std::vector<double>& knots, std::size_t order, std::size_t n, double v, const char* axis) {
    // Tested first on purpose: every range comparison below is false for
    // a NaN, so a NaN would otherwise fall straight through the guard.
    if (!std::isfinite(v)) {
        throw ValueError(std::string("TensorBSpline2D: non-finite ") + axis + " input");
    }
    const double lo = knots[order - 1];
    const double hi = knots[n];
    if (v < lo || v > hi) {
        throw ValueError(std::string("TensorBSpline2D: ") + axis + " = " + std::to_string(v) + " is outside the spline domain [" + std::to_string(lo)
                         + ", " + std::to_string(hi) + "]");
    }
}

inline std::size_t TensorBSpline2D::find_span(const std::vector<double>& knots, std::size_t order, std::size_t n, double v) {
    const std::size_t p = order - 1;
    if (v >= knots[n]) {
        return n - 1;
    }
    std::size_t lo = p;
    std::size_t hi = n;
    std::size_t mid = (lo + hi) / 2;
    while (v < knots[mid] || v >= knots[mid + 1]) {
        if (v < knots[mid]) {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    return mid;
}

inline void TensorBSpline2D::basis_funs(const std::vector<double>& knots, std::size_t order, std::size_t span, double v, double* out) {
    const std::size_t p = order - 1;
    double left[kMaxOrder];
    double right[kMaxOrder];
    out[0] = 1.0;
    for (std::size_t j = 1; j <= p; ++j) {
        left[j] = v - knots[span + 1 - j];
        right[j] = knots[span + j] - v;
        double saved = 0.0;
        for (std::size_t r = 0; r < j; ++r) {
            const double temp = out[r] / (right[r + 1] + left[j - r]);
            out[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        out[j] = saved;
    }
}

inline void TensorBSpline2D::basis_ders(const std::vector<double>& knots, std::size_t order, std::size_t span, double v, unsigned nd, double* out) {
    const int p = static_cast<int>(order) - 1;
    const int s = static_cast<int>(span);

    // Above the polynomial degree every derivative vanishes identically.
    //
    // Compared BEFORE narrowing to int: static_cast<int>(nd) is negative
    // for nd >= 2^31, which slipped past this guard and then indexed
    // ders[-1].  p >= 0 always, the constructor having rejected order 0.
    if (nd > static_cast<unsigned>(p)) {
        for (int j = 0; j <= p; ++j) {
            out[j] = 0.0;
        }
        return;
    }
    const int k_want = static_cast<int>(nd);

    double ndu[kMaxOrder][kMaxOrder];
    double a[2][kMaxOrder];
    double ders[kMaxOrder][kMaxOrder];
    double left[kMaxOrder];
    double right[kMaxOrder];

    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = v - knots[s + 1 - j];
        right[j] = knots[s + j] - v;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }
    for (int j = 0; j <= p; ++j) {
        ders[0][j] = ndu[j][p];
    }

    for (int r = 0; r <= p; ++r) {
        int s1 = 0;
        int s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= k_want; ++k) {
            double d = 0.0;
            const int rk = r - k;
            const int pk = p - k;
            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = (rk >= -1) ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            const int tmp = s1;
            s1 = s2;
            s2 = tmp;
        }
    }

    // The recurrence above omits the falling-factorial factor
    // p (p-1) ... (p-k+1) on the k-th derivative row; apply it now.
    //
    // Accumulated in double, not int.  k_want can reach p, so the largest
    // value used is p! -- 13! = 6.2e9 at order 14, already past INT_MAX,
    // rising to 15! = 1.3e12 at kMaxOrder.  In int this silently returned
    // derivatives wrong by O(1) RELATIVE error across part of the range
    // kMaxOrder advertises.  Every partial product p!/(p-k)! is an exact
    // integer below 2^53, so double represents all of them exactly.
    double factor = p;
    for (int k = 1; k <= k_want; ++k) {
        for (int j = 0; j <= p; ++j) {
            ders[k][j] *= factor;
        }
        factor *= static_cast<double>(p - k);
    }

    for (int j = 0; j <= p; ++j) {
        out[j] = ders[k_want][j];
    }
}

inline double TensorBSpline2D::eval(double x, double y, unsigned dx, unsigned dy) const {
    check_in_domain(kx_, ox_, nx_, x, "x");
    check_in_domain(ky_, oy_, ny_, y, "y");
    const std::size_t sx = find_span(kx_, ox_, nx_, x);
    const std::size_t sy = find_span(ky_, oy_, ny_, y);
    double bx[kMaxOrder];
    double by[kMaxOrder];
    if (dx == 0) {
        basis_funs(kx_, ox_, sx, x, bx);
    } else {
        basis_ders(kx_, ox_, sx, x, dx, bx);
    }
    if (dy == 0) {
        basis_funs(ky_, oy_, sy, y, by);
    } else {
        basis_ders(ky_, oy_, sy, y, dy, by);
    }

    const std::size_t px = ox_ - 1;
    const std::size_t py = oy_ - 1;
    double acc = 0.0;
    for (std::size_t a = 0; a < ox_; ++a) {
        const std::size_t i = sx - px + a;
        double inner = 0.0;
        for (std::size_t b = 0; b < oy_; ++b) {
            const std::size_t j = sy - py + b;
            inner += by[b] * coefs_[i * ny_ + j];
        }
        acc += bx[a] * inner;
    }
    // Defence in depth.  The degenerate-last-span guard in validate_axis
    // is believed to make a non-finite result unreachable from a
    // constructible surface via a zero denominator, so this is not the
    // primary protection -- it is the backstop for anything that analysis
    // missed, and for genuine floating-point overflow, which IS reachable:
    // finite but enormous coefficients (order 15, |c| ~ 1e300) overflow to
    // inf in a high derivative and throw here rather than returning inf.
    // That is a deliberate behaviour choice: this class's contract is that
    // a caller never receives a silent non-finite value, and a
    // thermodynamic property built on one would propagate it invisibly.
    // One predictable branch against ~100 flops of spline evaluation.
    if (!std::isfinite(acc)) {
        throw ValueError("TensorBSpline2D: evaluation produced a non-finite result");
    }
    return acc;
}

}  // namespace spline
}  // namespace CoolProp

#endif  // COOLPROP_SPLINE_TENSOR_BSPLINE_2D_H
