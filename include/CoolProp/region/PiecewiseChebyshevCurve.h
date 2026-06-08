#ifndef COOLPROP_REGION_PIECEWISE_CHEBYSHEV_CURVE_H
#define COOLPROP_REGION_PIECEWISE_CHEBYSHEV_CURVE_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "CoolProp/region/BoundaryCurve.h"

namespace CoolProp {
namespace region {

// 1D curve b = f(a) represented as a piecewise Chebyshev expansion.
//
// The primary axis can be parametrised in two ways:
//   - LINEAR : pieces are uniform in a; Chebyshev variable s in [-1,1]
//              maps to a via a = a_mid + (a_hi - a_lo)/2 * s.
//   - LOG    : pieces are uniform in log(a); s maps to log(a).
//
// Each piece carries its own Chebyshev coefficient vector c_0..c_N for
// the expansion f(a) = sum_{k=0}^{N} c_k T_k(s) on its sub-interval.
//
// Analytic derivative df/da is computed via the standard Chebyshev
// recurrence on the coefficients (no finite differences) — see
// eval_da().
//
// CONTINUITY CAVEAT.  The curve is value-continuous across piece
// boundaries (adjacent Chebyshev fits share their endpoint samples)
// but the *derivative* is generally discontinuous — each piece's
// expansion is independent.  When this curve drives the secondary
// axis of a Region, the kink at each piece boundary propagates
// through η normalization into the SVD reconstruction and appears as
// a thin band of elevated error at the corresponding axis value.
// Phase 2a's e2e validation surfaced this empirically for Water:
// SAT_PIECES=6 placed a visible "green band" at 4/6 of the log-p
// range (~770 kPa).  When a smoother boundary is required, use
// CubicSplineCurve (C² continuous everywhere).
//
// Extracted and generalised from the ihb/SBTL Cheb1D / Cheb1DPiece
// classes (originally parametrised on log(p) only).
class PiecewiseChebyshevCurve final : public BoundaryCurve
{
   public:
    enum class ParamScale : std::uint8_t
    {
        LINEAR,
        LOG
    };

    // Build by sampling f at Chebyshev-Lobatto nodes inside each of
    // n_pieces sub-intervals of [a_lo, a_hi], with each piece carrying a
    // degree-N expansion (so N+1 nodes per piece including the endpoints).
    //
    // f is called O(n_pieces * (N+1)) times.  Adjacent pieces share an
    // endpoint sample; the build computes them twice for simplicity (the
    // extra cost is in the wash and keeps the code clean).
    static std::unique_ptr<PiecewiseChebyshevCurve> build(double a_lo, double a_hi, std::size_t n_pieces, std::size_t degree, ParamScale scale,
                                                          const std::function<double(double)>& f);

    // Plain-data snapshot of a single piece's coefficient state.
    struct PieceState
    {
        double t_lo;
        double t_hi;
        double inv_half_span;
        double t_mid;
        std::vector<double> coeffs;        // size degree+1
        std::vector<double> deriv_coeffs;  // size degree (drops trivially-zero last coeff)
    };
    // Snapshot of the whole curve.
    struct State
    {
        double a_lo;
        double a_hi;
        ParamScale scale;
        std::vector<PieceState> pieces;
        double b_min;
        double b_max;
    };
    [[nodiscard]] State state() const;

    // Trusted factory: reconstructs from a previously-captured State
    // without re-sampling f.  Used by the SBTL deserializer.
    static std::unique_ptr<PiecewiseChebyshevCurve> from_state(State s);

    [[nodiscard]] double eval(double a) const noexcept override;
    [[nodiscard]] double eval_da(double a) const noexcept override;
    [[nodiscard]] std::pair<double, double> bounds() const noexcept override;
    [[nodiscard]] std::pair<double, double> a_range() const noexcept override;

   private:
    struct Piece
    {
        double t_lo;                 // transformed left bound (a or log a)
        double t_hi;                 // transformed right bound
        double inv_half_span;        // 2 / (t_hi - t_lo); for s = inv_half_span*(t - t_mid)
        double t_mid;                // (t_lo + t_hi) / 2
        std::vector<double> coeffs;  // size degree+1; Chebyshev coefficients
        // Pre-computed derivative coefficients on this piece — size
        // degree (the trivially-zero c_N coefficient of the derivative
        // expansion is dropped).  eval_da used to recompute these on
        // every call (heap allocation + recurrence); precomputing here
        // saves the malloc and one wasted Clenshaw step per call.
        std::vector<double> deriv_coeffs;
    };

    PiecewiseChebyshevCurve(double a_lo, double a_hi, ParamScale scale, std::vector<Piece> pieces, double b_min, double b_max);

    [[nodiscard]] std::size_t locate_piece(double t) const noexcept;

    // Transform a in [a_lo, a_hi] into the Chebyshev primary coordinate
    // (a itself for LINEAR, log(a) for LOG).
    [[nodiscard]] double to_t(double a) const noexcept;

    // d t / d a — the Jacobian of the primary-axis transform.
    //   LINEAR : 1
    //   LOG    : 1 / a
    [[nodiscard]] double dt_da(double a) const noexcept;

    double a_lo_;
    double a_hi_;
    ParamScale scale_;
    std::vector<Piece> pieces_;
    double b_min_;
    double b_max_;
};

}  // namespace region
}  // namespace CoolProp

#endif  // COOLPROP_REGION_PIECEWISE_CHEBYSHEV_CURVE_H
