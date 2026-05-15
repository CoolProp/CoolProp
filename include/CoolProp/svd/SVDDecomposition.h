#ifndef COOLPROP_SVD_SVD_DECOMPOSITION_H
#define COOLPROP_SVD_SVD_DECOMPOSITION_H

#include <cstdint>
#include <vector>

namespace CoolProp {
namespace svd {

// Output transform applied to the SVD result.  We store the SVD of a
// matrix in some transformed space (typically log of a strictly positive
// property), and undo the transform at evaluation time.
enum class OutputTransform : std::uint8_t
{
    IDENTITY,  // value = sum_k S_k u_k v_k   (no wrap)
    EXP        // value = exp(sum_k S_k u_k v_k)  (matrix was log of property)
};

// Slope source used to build the per-mode 1D slopes (dU_dx, dV_dy).
// Affects build only — the eval kernel reads slopes from arrays and is
// agnostic to how they were computed.
enum class SlopeSource : std::uint8_t
{
    NATURAL_CUBIC_SPLINE,  // default; C^2 continuous interpolant of each mode
    HERMITE_FD,            // central finite difference (svd_bench.cpp parity)
    PCHIP                  // monotone cubic; useful when monotonicity matters
};

// Plain-old-data container ferrying a 2D rank-r SVD from the offline
// builder to the runtime evaluator.
//
// Layout choices deliberately diverge from `dev/svd_bench.cpp` for
// performance:
//   - V is stored transposed relative to svd_bench.cpp's Vt: here
//     `V_S[j*rank + k]` is the k-th mode at y-index j, so the per-mode
//     scan in eval is unit-stride.
//   - The singular values S are *folded into* V_S — V_S[j,k] = V[j,k] * S[k].
//     The eval kernel becomes a plain dot product.  Storage saving:
//     rank doubles (negligible).  CPU saving: one multiply per mode per
//     call (small but free).
//   - dU_dx and dV_dy carry the precomputed slopes used by the cubic
//     Hermite kernel; one slope per (axis grid point, mode).
//
// The whole struct is trivially copyable except for the std::vectors;
// callers wanting GPU portability can `cudaMemcpy` the raw buffers
// pointed at by .data() and read them with an Eigen::Map.
struct SVDDecomposition
{
    std::int32_t NX = 0;    // number of grid points on x axis
    std::int32_t NY = 0;    // number of grid points on y axis
    std::int32_t rank = 0;  // truncation rank r (rank <= min(NX, NY))

    OutputTransform out_transform = OutputTransform::IDENTITY;
    SlopeSource slope_source = SlopeSource::NATURAL_CUBIC_SPLINE;

    std::vector<double> x_grid;  // size NX
    std::vector<double> y_grid;  // size NY

    // U: NX rows, rank columns, row-major.  U[i*rank + k] is mode k at x_i.
    std::vector<double> U;
    // dU/dx slopes co-located with U.
    std::vector<double> dU_dx;

    // V with S folded in: NY rows, rank columns, row-major.
    // V_S[j*rank + k] is sigma_k * V[j, k].
    std::vector<double> V_S;
    // dV_S/dy slopes co-located with V_S.
    std::vector<double> dV_S_dy;

    // Raw singular values kept for diagnostics / debug output.  Not used
    // by the hot path; the values are already folded into V_S.
    std::vector<double> S;
};

}  // namespace svd
}  // namespace CoolProp

#endif  // COOLPROP_SVD_SVD_DECOMPOSITION_H
