#ifndef COOLPROP_SVD_SVD_BUILDER_H
#define COOLPROP_SVD_SVD_BUILDER_H

#include <cstddef>
#include <vector>

#include "CoolProp/svd/SVDDecomposition.h"

namespace CoolProp {
namespace svd {

// Offline construction of a rank-r SVD decomposition from a 2D grid of
// values.
//
// The input matrix M is assumed to be *already transformed* — i.e. if
// the caller wants exp() applied at evaluation time, they pass log(M)
// here and set out_transform = OutputTransform::EXP.  The builder
// itself never logs or exps; it just SVD-truncates and fits 1D slopes
// per mode.
struct SVDBuildOptions
{
    std::int32_t rank = 0;  // r; required, > 0, <= min(NX, NY)
    OutputTransform out_transform = OutputTransform::IDENTITY;
    SlopeSource slope_source = SlopeSource::NATURAL_CUBIC_SPLINE;
};

// Build a rank-r decomposition.
//
// M is a flattened NX-by-NY matrix in row-major order: M[i*NY + j] is
// the value at (x_grid[i], y_grid[j]).
//
// x_grid and y_grid must be strictly increasing.  No NaN entries are
// tolerated; callers responsible for filling holes (e.g. dome cells)
// before calling.
SVDDecomposition build_svd(const std::vector<double>& x_grid, const std::vector<double>& y_grid, const std::vector<double>& M,
                           const SVDBuildOptions& opts);

}  // namespace svd
}  // namespace CoolProp

#endif  // COOLPROP_SVD_SVD_BUILDER_H
