#ifndef COOLPROP_SVD_SVD_EVALUATOR_H
#define COOLPROP_SVD_SVD_EVALUATOR_H

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>

#include "CoolProp/svd/Hermite1D.h"
#include "CoolProp/svd/SVDDecomposition.h"

namespace CoolProp {
namespace svd {

// Shared per-point setup that all SVDEvaluators over the same (x_grid,
// y_grid) reuse.  In an SVDSurface every per-property evaluator in a
// region shares the same grids — so the locate() and basis-weight
// computation can be done ONCE and amortized across N property evals
// at the same (x, y).  See SVDSurface::eval_with_region_multi for the
// batched-output entry point that uses this.
struct SVDEvalContext
{
    std::size_t i = 0;  // x_grid cell index, [0, NX-2]
    std::size_t j = 0;  // y_grid cell index, [0, NY-2]
    double hx = 0.0;
    double hy = 0.0;
    double tx = 0.0;
    double ty = 0.0;
    HermiteBasis bx{};
    HermiteBasis by{};
};

// Hot-path 2D rank-r SVD evaluator.
//
// Holds a non-owning pointer to an SVDDecomposition (built offline by
// SVDBuilder) and evaluates the truncated reconstruction at arbitrary
// (x, y).  The cubic-Hermite kernel on each axis is C^1 at minimum and
// C^2 when the decomposition's slopes come from a natural cubic spline
// (the default).
//
// Inner loop, per call:
//   - locate i  in x_grid (binary search, ~log NX comparisons)
//   - locate j  in y_grid (binary search, ~log NY comparisons)
//   - compute Hermite basis at t_x, t_y (8 mults + 6 adds total)
//   - dot product:  acc = sum_{k=0..r-1} u_k(x) * v_s_k(y)
//   - (optional) exp(acc)
//
// At rank 20 this is ~40 FMAs in the hot loop; the binary searches and
// Hermite basis computation dominate the per-call cost.
class SVDEvaluator
{
   public:
    // Construct from a borrowed decomposition.  The caller must keep the
    // SVDDecomposition alive *at a stable address* for at least as long
    // as this evaluator.  Rvalue-ref construction is deleted to prevent
    // the obvious temporary foot-gun:
    //
    //   SVDEvaluator e(build_svd(...));  // <-- temporary, dangling
    //
    // The less-obvious foot-gun involves moving the decomposition out
    // from under the evaluator via a container.  This pattern looks
    // safe but ISN'T, because vector::push_back moves the inline
    // `decomp` and the evaluator's borrowed pointer is left dangling:
    //
    //   struct RegionData {
    //       SVDDecomposition decomp;                  // INLINE — gets moved!
    //       std::unique_ptr<SVDEvaluator> eval;
    //   };
    //   std::vector<RegionData> regions;
    //   RegionData rd;
    //   rd.decomp = build_svd(...);
    //   rd.eval = std::make_unique<SVDEvaluator>(rd.decomp);
    //   regions.push_back(std::move(rd));   // ⚠  eval->d_ now dangles
    //
    // The right pattern is to give the decomposition a stable heap
    // address — either by storing the decomposition itself behind
    // std::unique_ptr / std::shared_ptr, or by constructing the
    // evaluator AFTER push_back from the in-vector decomp:
    //
    //   struct RegionData {
    //       std::unique_ptr<SVDDecomposition> decomp;  // heap = stable
    //       std::unique_ptr<SVDEvaluator> eval;
    //   };
    //
    // See dev/svd_sbtl_e2e.cpp for the canonical multi-region setup.
    //
    // Invariants on the decomposition (NX, NY, rank, and the sizes of
    // U/V_S/dU_dx/dV_S_dy/x_grid/y_grid) are checked once at
    // construction so the eval hot path can assume them.
    explicit SVDEvaluator(const SVDDecomposition& decomp) : d_(&decomp) {
        const auto nx = static_cast<std::size_t>(decomp.NX);
        const auto ny = static_cast<std::size_t>(decomp.NY);
        const auto r = static_cast<std::size_t>(decomp.rank);
        if (decomp.NX < 2 || decomp.NY < 2 || decomp.rank <= 0 || decomp.x_grid.size() != nx || decomp.y_grid.size() != ny
            || decomp.U.size() != nx * r || decomp.dU_dx.size() != nx * r || decomp.V_S.size() != ny * r || decomp.dV_S_dy.size() != ny * r) {
            throw std::invalid_argument("SVDEvaluator: SVDDecomposition has inconsistent dimensions");
        }
        // Grids must be strictly increasing: locate() does a binary
        // search assuming a sorted axis, and the hx / hy divisions in
        // eval() would blow up on duplicate or descending knots.
        // SVDBuilder enforces this at build time, but an
        // SVDDecomposition can be constructed by hand or loaded from a
        // file, so re-check on the evaluator side.
        for (std::size_t i = 1; i < nx; ++i) {
            if (!(decomp.x_grid[i - 1] < decomp.x_grid[i])) {
                throw std::invalid_argument("SVDEvaluator: x_grid must be strictly increasing");
            }
        }
        for (std::size_t j = 1; j < ny; ++j) {
            if (!(decomp.y_grid[j - 1] < decomp.y_grid[j])) {
                throw std::invalid_argument("SVDEvaluator: y_grid must be strictly increasing");
            }
        }
    }
    // Forbid construction from a temporary — eval() reads from the
    // borrowed decomposition after the constructor returns.
    SVDEvaluator(SVDDecomposition&&) = delete;
    explicit SVDEvaluator(const SVDDecomposition&&) = delete;

    // Owning constructor: take a shared_ptr to the decomposition and
    // hold it for the evaluator's lifetime.  This is the safer pattern
    // for callers that can't easily guarantee a stable address (the
    // dangling-pointer foot-gun documented above).  The performance
    // cost relative to the borrowed-ref form is a single pointer hop
    // at construction; the eval hot path is identical.
    explicit SVDEvaluator(std::shared_ptr<const SVDDecomposition> decomp)
      : SVDEvaluator(decomp ? *decomp
                            : throw std::invalid_argument(
                                "SVDEvaluator: shared_ptr<SVDDecomposition> must not be null"))  // delegate to the ref ctor for the invariant check
    {
        owned_ = std::move(decomp);
    }

    // Compute the per-point context (cell indices + Hermite basis
    // weights) once for an (x, y) input.  Any number of SVDEvaluators
    // sharing the same grids can apply eval_with_context() against the
    // returned struct — the SVD-multi-eval batched path in SVDSurface
    // is built on this pair.
    [[nodiscard]] inline SVDEvalContext make_context(double x, double y) const noexcept {
        const SVDDecomposition& d = *d_;
        SVDEvalContext c;
        c.i = locate(d.x_grid, x);
        c.j = locate(d.y_grid, y);
        c.hx = d.x_grid[c.i + 1] - d.x_grid[c.i];
        c.hy = d.y_grid[c.j + 1] - d.y_grid[c.j];
        c.tx = (x - d.x_grid[c.i]) / c.hx;
        c.ty = (y - d.y_grid[c.j]) / c.hy;
        c.bx = hermite_basis(c.tx);
        c.by = hermite_basis(c.ty);
        return c;
    }

    // Evaluate the rank-r reconstruction given a pre-computed context.
    // Callers that need K outputs at the same (x, y) should reuse one
    // context across the K eval_with_context calls — the locate +
    // Hermite-basis setup is then paid once instead of K times.
    [[nodiscard]] inline double eval_with_context(const SVDEvalContext& c) const noexcept {
        const SVDDecomposition& d = *d_;
        const auto r = static_cast<std::size_t>(d.rank);
        // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        const double* U0 = d.U.data() + c.i * r;
        const double* U1 = d.U.data() + (c.i + 1) * r;
        const double* dU0 = d.dU_dx.data() + c.i * r;
        const double* dU1 = d.dU_dx.data() + (c.i + 1) * r;
        const double* V0 = d.V_S.data() + c.j * r;
        const double* V1 = d.V_S.data() + (c.j + 1) * r;
        const double* dV0 = d.dV_S_dy.data() + c.j * r;
        const double* dV1 = d.dV_S_dy.data() + (c.j + 1) * r;
        // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        const double bx10_hx = c.bx.h10 * c.hx;
        const double bx11_hx = c.bx.h11 * c.hx;
        const double by10_hy = c.by.h10 * c.hy;
        const double by11_hy = c.by.h11 * c.hy;
        double acc = 0.0;
        // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        for (std::size_t k = 0; k < r; ++k) {
            const double u_k = c.bx.h00 * U0[k] + bx10_hx * dU0[k] + c.bx.h01 * U1[k] + bx11_hx * dU1[k];
            const double v_k = c.by.h00 * V0[k] + by10_hy * dV0[k] + c.by.h01 * V1[k] + by11_hy * dV1[k];
            acc += u_k * v_k;
        }
        // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
        return (d.out_transform == OutputTransform::EXP) ? std::exp(acc) : acc;
    }

    // Evaluate the rank-r reconstruction at (x, y).  No clamping; if
    // (x, y) lies outside the grid the Hermite kernel extrapolates from
    // the boundary cell, which can lose accuracy quickly — callers
    // should normally gate on the Region they came from.  Equivalent
    // to make_context(x, y) followed by eval_with_context().
    [[nodiscard]] inline double eval(double x, double y) const noexcept {
        return eval_with_context(make_context(x, y));
    }

   private:
    // Binary search: returns i such that grid[i] <= x <= grid[i+1],
    // clamped to [0, n-2] so that the i+1 access is always in range.
    // Closed on both ends at the boundaries — x == grid.back() returns
    // n-2, so the Hermite kernel sees t = 1 on the last cell rather
    // than t = 0 on a non-existent next cell.
    [[nodiscard]] static inline std::size_t locate(const std::vector<double>& g, double x) noexcept {
        const std::size_t n = g.size();
        if (x <= g.front()) {
            return 0;
        }
        if (x >= g.back()) {
            return n - 2;
        }
        std::size_t lo = 0;
        std::size_t hi = n - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (g[mid] <= x) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

    const SVDDecomposition* d_;
    // Set only when constructed via the shared_ptr overload; the
    // borrowed-ref constructor leaves this null and d_ aliases the
    // caller-owned object.
    std::shared_ptr<const SVDDecomposition> owned_;
};

}  // namespace svd
}  // namespace CoolProp

#endif  // COOLPROP_SVD_SVD_EVALUATOR_H
