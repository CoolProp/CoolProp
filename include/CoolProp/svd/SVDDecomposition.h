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

// Slope source used to build the per-mode 1D slopes (dU_dx, dV_S_dy).
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
//   - dU_dx and dV_S_dy carry the precomputed slopes used by the cubic
//     Hermite kernel; one slope per (axis grid point, mode).
//
// This struct is NOT trivially copyable — the std::vector members own
// heap storage, so a memcpy of the SVDDecomposition object itself
// produces garbage.  Callers wanting GPU portability should `cudaMemcpy`
// the underlying contiguous buffers (U.data(), V_S.data(), dU_dx.data(),
// dV_S_dy.data(), x_grid.data(), y_grid.data()) into device buffers
// individually and read them with an Eigen::Map on the device side.
// Device-side evaluation will agree with host evaluation within a few
// ulp — not byte-exact — because nvcc fuses multiply-and-accumulate by
// default (`--fmad=true`) and the slopes themselves were computed on
// the host with whatever rounding the spline solver chose.
struct SVDDecomposition
{
    std::int32_t NX = 0;    // number of grid points on x axis
    std::int32_t NY = 0;    // number of grid points on y axis
    std::int32_t rank = 0;  // truncation rank r (rank <= min(NX, NY))

    OutputTransform out_transform = OutputTransform::IDENTITY;
    // Provenance only — the eval kernel reads slopes from dU_dx /
    // dV_S_dy and is agnostic to which strategy filled them.  Kept on
    // the struct so a deserialised decomposition is self-describing.
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
