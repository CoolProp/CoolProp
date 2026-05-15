#ifndef COOLPROP_SVD_SVD_EVALUATOR_H
#define COOLPROP_SVD_SVD_EVALUATOR_H

#include <cmath>
#include <cstddef>
#include <utility>

#include "CoolProp/svd/Hermite1D.h"
#include "CoolProp/svd/SVDDecomposition.h"

namespace CoolProp {
namespace svd {

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
    explicit SVDEvaluator(const SVDDecomposition& decomp) noexcept : d_(&decomp) {}

    // Evaluate the rank-r reconstruction at (x, y).  No clamping; if
    // (x, y) lies outside the grid the Hermite kernel extrapolates from
    // the boundary cell, which can lose accuracy quickly — callers
    // should normally gate on the Region they came from.
    [[nodiscard]] inline double eval(double x, double y) const noexcept {
        const SVDDecomposition& d = *d_;
        const std::size_t i = locate(d.x_grid, x);
        const std::size_t j = locate(d.y_grid, y);
        const double hx = d.x_grid[i + 1] - d.x_grid[i];
        const double hy = d.y_grid[j + 1] - d.y_grid[j];
        const double tx = (x - d.x_grid[i]) / hx;
        const double ty = (y - d.y_grid[j]) / hy;
        const HermiteBasis bx = hermite_basis(tx);
        const HermiteBasis by = hermite_basis(ty);
        const std::size_t r = static_cast<std::size_t>(d.rank);
        const double* U0 = d.U.data() + i * r;
        const double* U1 = d.U.data() + (i + 1) * r;
        const double* dU0 = d.dU_dx.data() + i * r;
        const double* dU1 = d.dU_dx.data() + (i + 1) * r;
        const double* V0 = d.V_S.data() + j * r;
        const double* V1 = d.V_S.data() + (j + 1) * r;
        const double* dV0 = d.dV_S_dy.data() + j * r;
        const double* dV1 = d.dV_S_dy.data() + (j + 1) * r;
        const double bx10_hx = bx.h10 * hx;
        const double bx11_hx = bx.h11 * hx;
        const double by10_hy = by.h10 * hy;
        const double by11_hy = by.h11 * hy;
        double acc = 0.0;
        for (std::size_t k = 0; k < r; ++k) {
            const double u_k = bx.h00 * U0[k] + bx10_hx * dU0[k] + bx.h01 * U1[k] + bx11_hx * dU1[k];
            const double v_k = by.h00 * V0[k] + by10_hy * dV0[k] + by.h01 * V1[k] + by11_hy * dV1[k];
            acc += u_k * v_k;
        }
        return (d.out_transform == OutputTransform::EXP) ? std::exp(acc) : acc;
    }

   private:
    // Binary search: returns i such that grid[i] <= x < grid[i+1],
    // clamped to [0, n-2] so that the i+1 access is always in range.
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
};

}  // namespace svd
}  // namespace CoolProp

#endif  // COOLPROP_SVD_SVD_EVALUATOR_H
