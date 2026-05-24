#if !defined(NO_TABULAR_BACKENDS)

#    include "BicubicBackend.h"

#    include <cmath>
#    include <limits>
#    include <optional>
#    include "MatrixMath.h"
#    include "DataStructures.h"
#    include "Backends/Helmholtz/PhaseEnvelopeRoutines.h"

namespace {
// Return the cubic root that lies in [0, 1] (the cell's normalized
// coordinate), or std::nullopt if no root is in-cell.
//
// solve_cubic returns up to three real roots in arbitrary order. The
// historical code picked the smallest-absolute-value root, which can be
// a far-negative root with no physical meaning and produces wildly wrong
// values once unscaled to the cell width (#1301: P = -760 MPa for CO2 at
// T=315K, rho=1 kg/m^3).
//
// Strict cell bound here is intentional: the bicubic interpolant has
// physical meaning only inside the cell. If the bisection landed in the
// wrong cell, the caller walks to a neighbour and retries; if no
// neighbouring cell contains the root either, the requested state is
// outside the table's domain and the caller throws.
//
// When multiple roots are in [0, 1] (can happen for cells near a
// saturation curve where the interpolant is non-monotonic), prefer the
// one closest to 0.5 — that's furthest from the cell edges and least
// sensitive to coefficient noise.
constexpr double kCellRootTol = 1e-6;

inline std::optional<double> find_in_cell_root(int N, double r0, double r1, double r2) {
    const double roots[3] = {r0, r1, r2};
    int best = -1;
    double best_score = std::numeric_limits<double>::infinity();
    for (int k = 0; k < N; ++k) {
        const double r = roots[k];
        if (r >= -kCellRootTol && r <= 1.0 + kCellRootTol) {
            const double score = std::abs(r - 0.5);
            if (score < best_score) {
                best_score = score;
                best = k;
            }
        }
    }
    if (best < 0) {
        return std::nullopt;
    }
    return roots[best];
}

// Maximum number of cells to walk away from the bisection's initial
// guess.  find_nearest_neighbor (TabularBackends.h) now bilinearly
// interpolates the otherkey column at the actual query position before
// bisecting (xkey branch), so the chosen cell is usually correct.  The
// remaining off-by-one cases come from the bilinear estimator being a
// linear approximation of the bicubic interpolant — the cubic in the
// neighbouring cell can hold the true root.  The immediate neighbour
// (one cell either side) is the only physically meaningful source:
// that's a real cell with its own valid cubic, not an extrapolation of
// the original cell's interpolant.  Walking further would just be
// extrapolation in disguise and would mask genuinely out-of-table
// queries — those are caught here by throwing.  #1301's catastrophic
// case lands ~50+ cells away, well outside this bound.
constexpr int kMaxCellWalk = 1;
}  // namespace

void CoolProp::BicubicBackend::find_native_nearest_good_indices(SinglePhaseGriddedTableData& table,
                                                                const std::vector<std::vector<CellCoeffs>>& coeffs, double x, double y,
                                                                std::size_t& i, std::size_t& j) {
    table.find_native_nearest_good_cell(x, y, i, j);
    const CellCoeffs& cell = coeffs[i][j];
    if (!cell.valid()) {
        if (auto alt = cell.get_alternate()) {
            // Get new good neighbor
            auto [ai, aj] = *alt;
            i = ai;
            j = aj;
        } else {
            if (!cell.valid()) {
                throw ValueError(format("Cell is invalid and has no good neighbors for x = %g, y= %g", x, y));
            }
        }
    }
}

/// Ask the derived class to find the nearest neighbor (pure virtual)
void CoolProp::BicubicBackend::find_nearest_neighbor(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                     const parameters variable1, const double value1, const parameters otherkey,
                                                     const double otherval, std::size_t& i, std::size_t& j) {
    table.find_nearest_neighbor(variable1, value1, otherkey, otherval, i, j);
    const CellCoeffs& cell = coeffs[i][j];
    if (!cell.valid()) {
        if (auto alt = cell.get_alternate()) {
            // Get new good neighbor
            auto [ai, aj] = *alt;
            i = ai;
            j = aj;
        } else {
            if (!cell.valid()) {
                throw ValueError(format("Cell is invalid and has no good neighbors for x = %g, y = %g", value1, otherval));
            }
        }
    }
}

/** Use the single_phase table to evaluate an output for a transport property
 *
 * Here we use linear interpolation because we don't have any information about the derivatives with respect to the
 * independent variables and it is too computationally expensive to build the derivatives numerically
 *
 * See also http://en.wikipedia.org/wiki/Bilinear_interpolation#Nonlinear
 */
double CoolProp::BicubicBackend::evaluate_single_phase_transport(SinglePhaseGriddedTableData& table, parameters output, double x, double y,
                                                                 std::size_t i, std::size_t j) {
    // By definition i,i+1,j,j+1 are all in range and valid
    std::vector<std::vector<double>>* f = nullptr;
    switch (output) {
        case iconductivity:
            f = &table.cond;
            break;
        case iviscosity:
            f = &table.visc;
            break;
        default:
            throw ValueError(format("invalid output variable to BicubicBackend::evaluate_single_phase_transport"));
    }
    double x1 = table.xvec[i], x2 = table.xvec[i + 1], y1 = table.yvec[j], y2 = table.yvec[j + 1];
    double f11 = (*f)[i][j], f12 = (*f)[i][j + 1], f21 = (*f)[i + 1][j], f22 = (*f)[i + 1][j + 1];
    double val =
      1 / ((x2 - x1) * (y2 - y1)) * (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y) + f12 * (x2 - x) * (y - y1) + f22 * (x - x1) * (y - y1));

    // Cache the output value calculated
    switch (output) {
        case iconductivity:
            _conductivity = val;
            break;
        case iviscosity:
            _viscosity = val;
            break;
        default:
            throw ValueError("Invalid output variable in evaluate_single_phase_transport");
    }
    return val;
}
// Use the single_phase table to evaluate an output
double CoolProp::BicubicBackend::evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                       const parameters output, const double x, const double y, const std::size_t i,
                                                       const std::size_t j) {
    // Get the cell
    const CellCoeffs& cell = coeffs[i][j];

    // Defense-in-depth: cells in the table's two-phase notch carry no bicubic
    // coefficients. Caller should have routed through find_native_nearest_good_indices
    // (or the saturation-curve cell-bump logic in TabularBackend::update) which
    // already resolve to a good neighbour. This guard catches the residual case
    // (e.g. #1950: PT_INPUTS at sub-saturation) so an invalid cell raises a
    // ValueError instead of dereferencing an empty alpha[] vector and segfaulting.
    if (!cell.valid()) {
        throw ValueError(format("evaluate_single_phase called on cell (%zu, %zu) with no bicubic coefficients (x=%g, y=%g)", i, j, x, y));
    }

    // Get the alpha coefficients
    const std::vector<double>& alpha = cell.get(output);

    // Normalized value in the range (0, 1)
    double xhat = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    double yhat = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);

    // Calculate the output value desired
    // Term multiplying x^0 using Horner's method
    double B0 = ((((0) + alpha[3 * 4 + 0]) * yhat + alpha[2 * 4 + 0]) * yhat + alpha[1 * 4 + 0]) * yhat + alpha[0 * 4 + 0];
    // Term multiplying x^1 using Horner's method
    double B1 = ((((0) + alpha[3 * 4 + 1]) * yhat + alpha[2 * 4 + 1]) * yhat + alpha[1 * 4 + 1]) * yhat + alpha[0 * 4 + 1];
    // Term multiplying x^2 using Horner's method
    double B2 = ((((0) + alpha[3 * 4 + 2]) * yhat + alpha[2 * 4 + 2]) * yhat + alpha[1 * 4 + 2]) * yhat + alpha[0 * 4 + 2];
    // Term multiplying x^3 using Horner's method
    double B3 = ((((0) + alpha[3 * 4 + 3]) * yhat + alpha[2 * 4 + 3]) * yhat + alpha[1 * 4 + 3]) * yhat + alpha[0 * 4 + 3];

    double val = ((((0) + B3) * xhat + B2) * xhat + B1) * xhat + B0;

    // Cache the output value calculated
    switch (output) {
        case iT:
            _T = val;
            break;
        case iDmolar:
            _rhomolar = val;
            break;
        case iSmolar:
            _smolar = val;
            break;
        case iHmolar:
            _hmolar = val;
            break;
        case iUmolar:
            _umolar = val;
            break;
        default:
            throw ValueError("Invalid output variable in evaluate_single_phase");
    }
    return val;
}
/// Use the single_phase table to evaluate an output
double CoolProp::BicubicBackend::evaluate_single_phase_derivative(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs,
                                                                  parameters output, double x, double y, std::size_t i, std::size_t j, std::size_t Nx,
                                                                  std::size_t Ny) {

    // Get the cell
    CellCoeffs& cell = coeffs[i][j];

    // Get the alpha coefficients
    const std::vector<double>& alpha = cell.get(output);

    // Normalized value in the range (0, 1)
    double xhat = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    double yhat = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);
    double dxhatdx = 1 / (table.xvec[i + 1] - table.xvec[i]);
    double dyhatdy = 1 / (table.yvec[j + 1] - table.yvec[j]);

    // Calculate the output value desired
    double val = 0;
    if (Nx == 1 && Ny == 0) {
        if (output == table.xkey) {
            return 1.0;
        }
        if (output == table.ykey) {
            return 0.0;
        }
        for (std::size_t l = 1; l < 4; ++l) {
            for (std::size_t m = 0; m < 4; ++m) {
                val += alpha[m * 4 + l] * l * pow(xhat, static_cast<int>(l - 1)) * pow(yhat, static_cast<int>(m));
            }
        }
        // val is now dz/dxhat|yhat
        return val * dxhatdx;
    } else if (Ny == 1 && Nx == 0) {
        if (output == table.ykey) {
            return 1.0;
        }
        if (output == table.xkey) {
            return 0.0;
        }
        for (std::size_t l = 0; l < 4; ++l) {
            for (std::size_t m = 1; m < 4; ++m) {
                val += alpha[m * 4 + l] * pow(xhat, static_cast<int>(l)) * m * pow(yhat, static_cast<int>(m - 1));
            }
        }
        // val is now dz/dyhat|xhat
        return val * dyhatdy;
    } else {
        throw ValueError("Invalid input");
    }
}

/// Use the single_phase table to invert for x given a y
///
/// Iterate cells in the i (xkey-axis) direction looking for one whose cubic
/// has a root inside the unit interval. find_nearest_neighbor's i bisection
/// (via bisect_segmented_vector_slice on the otherkey matrix at the chosen
/// j) can be off-by-one when the y-dependence is significant; walking ±1..N
/// cells recovers the correct cell so the final root lies in [0, 1] rather
/// than relying on out-of-cell cubic extrapolation. If no neighbour holds
/// the root either, the state is outside the table's domain — throw
/// instead of returning a wildly wrong (often negative) value (#1301).
void CoolProp::BicubicBackend::invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                     parameters other_key, double other, double y, std::size_t i, std::size_t j) {
    // yhat is fixed (same for every cell we try in the i direction)
    const double yhat = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);
    const double y_0 = 1, y_1 = yhat, y_2 = yhat * yhat, y_3 = yhat * yhat * yhat;

    auto try_cell = [&](std::size_t i_try) -> std::optional<double> {
        const CellCoeffs& cell = coeffs[i_try][j];
        if (!cell.valid()) {
            return std::nullopt;
        }
        const std::vector<double>& alpha = cell.get(other_key);
        const double a = alpha[3 + 0 * 4] * y_0 + alpha[3 + 1 * 4] * y_1 + alpha[3 + 2 * 4] * y_2 + alpha[3 + 3 * 4] * y_3;
        const double b = alpha[2 + 0 * 4] * y_0 + alpha[2 + 1 * 4] * y_1 + alpha[2 + 2 * 4] * y_2 + alpha[2 + 3 * 4] * y_3;
        const double c = alpha[1 + 0 * 4] * y_0 + alpha[1 + 1 * 4] * y_1 + alpha[1 + 2 * 4] * y_2 + alpha[1 + 3 * 4] * y_3;
        const double d = alpha[0 + 0 * 4] * y_0 + alpha[0 + 1 * 4] * y_1 + alpha[0 + 2 * 4] * y_2 + alpha[0 + 3 * 4] * y_3 - other;
        int N = 0;
        double r0 = NAN, r1 = NAN, r2 = NAN;
        solve_cubic(a, b, c, d, N, r0, r1, r2);
        if (N == 0) {
            return std::nullopt;
        }
        return find_in_cell_root(N, r0, r1, r2);
    };

    std::optional<double> xhat_opt = try_cell(i);
    std::size_t i_found = i;
    if (!xhat_opt) {
        for (int delta = 1; delta <= kMaxCellWalk; ++delta) {
            // Bound checks: cell coefficients are indexed by [i][j], cells go up to
            // table.xvec.size() - 2 (inclusive); the last index has no [i+1] neighbour.
            if (i + delta + 1 < table.xvec.size()) {
                xhat_opt = try_cell(i + delta);
                if (xhat_opt) {
                    i_found = i + delta;
                    break;
                }
            }
            if (i >= static_cast<std::size_t>(delta)) {
                xhat_opt = try_cell(i - delta);
                if (xhat_opt) {
                    i_found = i - delta;
                    break;
                }
            }
        }
    }

    if (!xhat_opt) {
        throw ValueError(format("BICUBIC: invert_single_phase_x could not find an in-cell cubic root for cell "
                                "(i=%zu, j=%zu) spanning xvec=[%g,%g], yvec=[%g,%g] or any neighbour within +/-%d "
                                "cells in the i direction (other_key=%d, other=%g, y=%g); requested state likely "
                                "lies outside the BICUBIC table's domain",
                                i, j, table.xvec[i], table.xvec[i + 1], table.yvec[j], table.yvec[j + 1], kMaxCellWalk, static_cast<int>(other_key),
                                other, y));
    }

    // Cache the cell we actually used so downstream evaluations see the
    // corrected position rather than the original bisection's mistake.
    i = i_found;
    cached_single_phase_i = i;
    const double xhat = *xhat_opt;
    const double val = xhat * (table.xvec[i + 1] - table.xvec[i]) + table.xvec[i];

    switch (table.xkey) {
        case iHmolar:
            _hmolar = val;
            break;
        case iT:
            _T = val;
            break;
        default:
            throw ValueError("Invalid output variable in invert_single_phase_x");
    }
}

/// Use the single_phase table to solve for y given an x
///
/// Iterate cells in the j (ykey-axis) direction looking for one whose cubic
/// has a root inside the unit interval. find_nearest_neighbor's j bisection
/// (via bisect_vector on the otherkey column at the lower-xkey edge of the
/// cell) can be off-by-one when the xkey-dependence is significant — e.g.
/// for the CO2 PT-table at T=320K, density along the supercritical
/// isotherm differs meaningfully from the lower edge of the i-cell, so the
/// bisection systematically misses by one j. Walking ±1..N cells recovers
/// the correct cell. If no neighbour holds the root either, the state is
/// outside the table's domain — throw instead of returning a wildly wrong
/// (often negative) value (#1301).
void CoolProp::BicubicBackend::invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                     parameters other_key, double other, double x, std::size_t i, std::size_t j) {
    // xhat is fixed (same for every cell we try in the j direction)
    const double xhat = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    const double x_0 = 1, x_1 = xhat, x_2 = xhat * xhat, x_3 = xhat * xhat * xhat;

    auto try_cell = [&](std::size_t j_try) -> std::optional<double> {
        const CellCoeffs& cell = coeffs[i][j_try];
        if (!cell.valid()) {
            return std::nullopt;
        }
        const std::vector<double>& alpha = cell.get(other_key);
        const double a = alpha[0 + 3 * 4] * x_0 + alpha[1 + 3 * 4] * x_1 + alpha[2 + 3 * 4] * x_2 + alpha[3 + 3 * 4] * x_3;
        const double b = alpha[0 + 2 * 4] * x_0 + alpha[1 + 2 * 4] * x_1 + alpha[2 + 2 * 4] * x_2 + alpha[3 + 2 * 4] * x_3;
        const double c = alpha[0 + 1 * 4] * x_0 + alpha[1 + 1 * 4] * x_1 + alpha[2 + 1 * 4] * x_2 + alpha[3 + 1 * 4] * x_3;
        const double d = alpha[0 + 0 * 4] * x_0 + alpha[1 + 0 * 4] * x_1 + alpha[2 + 0 * 4] * x_2 + alpha[3 + 0 * 4] * x_3 - other;
        int N = 0;
        double r0 = NAN, r1 = NAN, r2 = NAN;
        solve_cubic(a, b, c, d, N, r0, r1, r2);
        if (N == 0) {
            return std::nullopt;
        }
        return find_in_cell_root(N, r0, r1, r2);
    };

    std::optional<double> yhat_opt = try_cell(j);
    std::size_t j_found = j;
    if (!yhat_opt) {
        for (int delta = 1; delta <= kMaxCellWalk; ++delta) {
            if (j + delta + 1 < table.yvec.size()) {
                yhat_opt = try_cell(j + delta);
                if (yhat_opt) {
                    j_found = j + delta;
                    break;
                }
            }
            if (j >= static_cast<std::size_t>(delta)) {
                yhat_opt = try_cell(j - delta);
                if (yhat_opt) {
                    j_found = j - delta;
                    break;
                }
            }
        }
    }

    if (!yhat_opt) {
        throw ValueError(format("BICUBIC: invert_single_phase_y could not find an in-cell cubic root for cell "
                                "(i=%zu, j=%zu) spanning xvec=[%g,%g], yvec=[%g,%g] or any neighbour within +/-%d "
                                "cells in the j direction (other_key=%d, other=%g, x=%g); requested state likely "
                                "lies outside the BICUBIC table's domain",
                                i, j, table.xvec[i], table.xvec[i + 1], table.yvec[j], table.yvec[j + 1], kMaxCellWalk, static_cast<int>(other_key),
                                other, x));
    }

    j = j_found;
    cached_single_phase_j = j;
    const double yhat = *yhat_opt;
    const double val = yhat * (table.yvec[j + 1] - table.yvec[j]) + table.yvec[j];

    // Cache the output value calculated
    switch (table.ykey) {
        case iP:
            _p = val;
            break;
        default:
            throw ValueError("Invalid output variable in invert_single_phase_x");
    }
}

#endif  // !defined(NO_TABULAR_BACKENDS)
