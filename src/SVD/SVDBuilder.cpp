#include "CoolProp/svd/SVDBuilder.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// Use Eigen's BDC SVD: divide-and-conquer is more efficient than Jacobi
// for the modest-sized (NX*NY ~ a few * 10^5) matrices we deal with.
#include <Eigen/Core>
#include <Eigen/SVD>

#include "CoolProp/region/CubicSplineCurve.h"

namespace CoolProp {
namespace svd {

namespace {

// Central-difference slopes on a (possibly non-uniform) grid.
//
//   m_0       = (y_1 - y_0) / (x_1 - x_0)
//   m_{N-1}   = (y_{N-1} - y_{N-2}) / (x_{N-1} - x_{N-2})
//   m_i       = (y_{i+1} - y_{i-1}) / (x_{i+1} - x_{i-1})  for 1 <= i <= N-2
std::vector<double> fd_slopes(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    std::vector<double> m(n, 0.0);
    if (n < 2) {
        return m;
    }
    m[0] = (y[1] - y[0]) / (x[1] - x[0]);
    m[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
    for (std::size_t i = 1; i + 1 < n; ++i) {
        m[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1]);
    }
    return m;
}

// Slopes from a natural cubic spline through (x, y).  Reuses the
// general-purpose CubicSplineCurve so we don't keep two copies of the
// tridiagonal solver around.
std::vector<double> spline_slopes(const std::vector<double>& x, const std::vector<double>& y) {
    auto spline = CoolProp::region::CubicSplineCurve::build(x, y);
    const std::size_t n = x.size();
    std::vector<double> m(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        m[i] = spline->eval_da(x[i]);
    }
    return m;
}

// Monotone cubic (PCHIP / Fritsch-Carlson) slopes.  Each interior slope
// is constructed so the resulting Hermite cubic is monotone on every
// segment when the input data is monotone.  At smooth maxima/minima the
// slope is forced to zero to prevent overshoot.
//
// References:
//   Fritsch & Carlson (1980), SIAM J. Numer. Anal. 17(2): 238-246.
std::vector<double> pchip_slopes(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    std::vector<double> m(n, 0.0);
    if (n < 2) {
        return m;
    }
    if (n == 2) {
        const double slope = (y[1] - y[0]) / (x[1] - x[0]);
        m[0] = slope;
        m[1] = slope;
        return m;
    }
    std::vector<double> h(n - 1, 0.0);
    std::vector<double> delta(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
        delta[i] = (y[i + 1] - y[i]) / h[i];
    }
    // Interior slopes via Fritsch-Carlson weighted harmonic mean.
    for (std::size_t i = 1; i + 1 < n; ++i) {
        if (delta[i - 1] * delta[i] <= 0.0) {
            m[i] = 0.0;
        } else {
            const double w1 = 2.0 * h[i] + h[i - 1];
            const double w2 = h[i] + 2.0 * h[i - 1];
            m[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
        }
    }
    // Endpoint slopes via three-point Steffen-style formulas, clamped to
    // preserve monotonicity at the boundary.
    m[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
    m[n - 1] = ((2.0 * h[n - 2] + h[n - 3]) * delta[n - 2] - h[n - 2] * delta[n - 3]) / (h[n - 2] + h[n - 3]);
    // Clamp endpoint slopes so they don't introduce non-monotonicity.
    if (m[0] * delta[0] <= 0.0) {
        m[0] = 0.0;
    } else if (std::abs(m[0]) > 3.0 * std::abs(delta[0])) {
        m[0] = 3.0 * delta[0];
    }
    if (m[n - 1] * delta[n - 2] <= 0.0) {
        m[n - 1] = 0.0;
    } else if (std::abs(m[n - 1]) > 3.0 * std::abs(delta[n - 2])) {
        m[n - 1] = 3.0 * delta[n - 2];
    }
    return m;
}

std::vector<double> compute_slopes(const std::vector<double>& x, const std::vector<double>& y, SlopeSource src) {
    switch (src) {
        case SlopeSource::HERMITE_FD:
            return fd_slopes(x, y);
        case SlopeSource::NATURAL_CUBIC_SPLINE:
            return spline_slopes(x, y);
        case SlopeSource::PCHIP:
            return pchip_slopes(x, y);
    }
    return fd_slopes(x, y);  // unreachable
}

}  // namespace

SVDDecomposition build_svd(const std::vector<double>& x_grid, const std::vector<double>& y_grid, const std::vector<double>& M,
                           const SVDBuildOptions& opts) {
    const std::size_t NX = x_grid.size();
    const std::size_t NY = y_grid.size();
    if (NX < 2 || NY < 2) {
        throw std::invalid_argument("build_svd: x_grid and y_grid must each have >= 2 entries");
    }
    if (M.size() != NX * NY) {
        throw std::invalid_argument("build_svd: M must be NX * NY entries in row-major order");
    }
    for (std::size_t i = 1; i < NX; ++i) {
        if (!(x_grid[i] > x_grid[i - 1])) {
            throw std::invalid_argument("build_svd: x_grid must be strictly increasing");
        }
    }
    for (std::size_t j = 1; j < NY; ++j) {
        if (!(y_grid[j] > y_grid[j - 1])) {
            throw std::invalid_argument("build_svd: y_grid must be strictly increasing");
        }
    }
    const std::size_t full_rank = std::min(NX, NY);
    if (opts.rank <= 0 || static_cast<std::size_t>(opts.rank) > full_rank) {
        throw std::invalid_argument("build_svd: rank must be in (0, min(NX, NY)]");
    }
    for (const double v : M) {
        if (!std::isfinite(v)) {
            throw std::invalid_argument("build_svd: M contains a non-finite entry");
        }
    }

    const auto r = static_cast<std::size_t>(opts.rank);

    // Eigen matrix view over M (row-major).
    using RowMajorMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::Map<const RowMajorMat> M_map(M.data(), static_cast<Eigen::Index>(NX), static_cast<Eigen::Index>(NY));
    // Eigen 3.5 / 5.0+ promoted ComputeThinU/V from a runtime constructor
    // argument to a template parameter on BDCSVD; the older form is
    // EIGEN_DEPRECATED in 3.5+ and slated for removal in the next major.
    // Use the template form when the installed Eigen supports it and
    // fall back to the runtime form on older installs.
#if defined(EIGEN_WORLD_VERSION) && (EIGEN_WORLD_VERSION > 3 || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION >= 5))
    Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M_map);
#else
    Eigen::BDCSVD<Eigen::MatrixXd> svd(M_map, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
    if (svd.info() != Eigen::Success) {
        throw std::runtime_error("build_svd: BDCSVD failed to converge");
    }

    const Eigen::MatrixXd& U_full = svd.matrixU();         // NX x min(NX,NY)
    const Eigen::MatrixXd& V_full = svd.matrixV();         // NY x min(NX,NY)
    const Eigen::VectorXd& s_full = svd.singularValues();  // min(NX,NY)

    SVDDecomposition d;
    d.NX = static_cast<std::int32_t>(NX);
    d.NY = static_cast<std::int32_t>(NY);
    d.rank = opts.rank;
    d.out_transform = opts.out_transform;
    d.slope_source = opts.slope_source;
    d.x_grid = x_grid;
    d.y_grid = y_grid;
    d.S.resize(r);
    for (std::size_t k = 0; k < r; ++k) {
        d.S[k] = s_full(static_cast<Eigen::Index>(k));
    }
    // U: row-major (NX, r).
    d.U.assign(NX * r, 0.0);
    for (std::size_t i = 0; i < NX; ++i) {
        for (std::size_t k = 0; k < r; ++k) {
            d.U[i * r + k] = U_full(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k));
        }
    }
    // V_S: row-major (NY, r) with sigma_k folded in.
    d.V_S.assign(NY * r, 0.0);
    for (std::size_t j = 0; j < NY; ++j) {
        for (std::size_t k = 0; k < r; ++k) {
            d.V_S[j * r + k] = V_full(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) * d.S[k];
        }
    }

    // Per-mode 1D slopes for the Hermite kernel.
    d.dU_dx.assign(NX * r, 0.0);
    d.dV_S_dy.assign(NY * r, 0.0);
    std::vector<double> u_col(NX);
    std::vector<double> v_col(NY);
    for (std::size_t k = 0; k < r; ++k) {
        for (std::size_t i = 0; i < NX; ++i) {
            u_col[i] = d.U[i * r + k];
        }
        const std::vector<double> u_slopes = compute_slopes(x_grid, u_col, opts.slope_source);
        for (std::size_t i = 0; i < NX; ++i) {
            d.dU_dx[i * r + k] = u_slopes[i];
        }
        for (std::size_t j = 0; j < NY; ++j) {
            v_col[j] = d.V_S[j * r + k];
        }
        const std::vector<double> v_slopes = compute_slopes(y_grid, v_col, opts.slope_source);
        for (std::size_t j = 0; j < NY; ++j) {
            d.dV_S_dy[j * r + k] = v_slopes[j];
        }
    }

    return d;
}

}  // namespace svd
}  // namespace CoolProp
