#if !defined(NO_TABULAR_BACKENDS)

#    include "BicubicBackend.h"

#    include <cmath>
#    include <limits>
#    include "MatrixMath.h"
#    include "DataStructures.h"
#    include "Backends/Helmholtz/PhaseEnvelopeRoutines.h"

namespace {
// Pick the cubic root that lies in the unit interval [0, 1] (the cell's
// normalized coordinate). solve_cubic returns up to three real roots in
// arbitrary order; the historical code took the smallest-absolute-value
// root, which can be a large-negative root that has nothing to do with the
// physical solution and produces wildly wrong values once unscaled to the
// cell width (e.g. #1301: P = -760 MPa for CO2 at T=315K, rho=1 kg/m^3).
//
// Strategy:
//   1. Among the N roots, prefer those with r in [-tol, 1+tol] (in-cell
//      or within a small slop band). Pick the one closest to the cell
//      centre 0.5 (most numerically stable).
//   2. If no root is in-range, fall back to the closest-to-[0,1] root —
//      i.e. minimize max(0, -r) + max(0, r-1) — so we still emit something
//      bounded, and the caller can detect "outside cell" by the unscaled
//      value falling outside [vec[k], vec[k+1]].
inline double pick_unit_interval_root(int N, double r0, double r1, double r2) {
    constexpr double tol = 1e-6;
    const double roots[3] = {r0, r1, r2};
    int best_in = -1;
    double best_in_score = std::numeric_limits<double>::infinity();
    int best_any = 0;
    double best_any_dist = std::numeric_limits<double>::infinity();
    for (int k = 0; k < N; ++k) {
        const double r = roots[k];
        if (r >= -tol && r <= 1.0 + tol) {
            const double score = std::abs(r - 0.5);
            if (score < best_in_score) {
                best_in_score = score;
                best_in = k;
            }
        }
        const double dist = std::max(0.0, -r) + std::max(0.0, r - 1.0);
        if (dist < best_any_dist) {
            best_any_dist = dist;
            best_any = k;
        }
    }
    return (best_in >= 0) ? roots[best_in] : roots[best_any];
}
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
void CoolProp::BicubicBackend::invert_single_phase_x(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                     parameters other_key, double other, double y, std::size_t i, std::size_t j) {
    // Get the cell
    const CellCoeffs& cell = coeffs[i][j];

    // Get the alpha coefficients
    const std::vector<double>& alpha = cell.get(other_key);

    // Normalized value in the range (0, 1)
    double yhat = (y - table.yvec[j]) / (table.yvec[j + 1] - table.yvec[j]);

    double y_0 = 1, y_1 = yhat, y_2 = yhat * yhat, y_3 = yhat * yhat * yhat;

    double a = alpha[3 + 0 * 4] * y_0 + alpha[3 + 1 * 4] * y_1 + alpha[3 + 2 * 4] * y_2 + alpha[3 + 3 * 4] * y_3;          // factors of xhat^3
    double b = alpha[2 + 0 * 4] * y_0 + alpha[2 + 1 * 4] * y_1 + alpha[2 + 2 * 4] * y_2 + alpha[2 + 3 * 4] * y_3;          // factors of xhar^2
    double c = alpha[1 + 0 * 4] * y_0 + alpha[1 + 1 * 4] * y_1 + alpha[1 + 2 * 4] * y_2 + alpha[1 + 3 * 4] * y_3;          // factors of xhat
    double d = alpha[0 + 0 * 4] * y_0 + alpha[0 + 1 * 4] * y_1 + alpha[0 + 2 * 4] * y_2 + alpha[0 + 3 * 4] * y_3 - other;  // constant factors
    int N = 0;
    double xhat0 = NAN, xhat1 = NAN, xhat2 = NAN, val = NAN, xhat = _HUGE;
    solve_cubic(a, b, c, d, N, xhat0, xhat1, xhat2);
    if (N == 0) {
        throw ValueError("Could not find a solution in invert_single_phase_x");
    }
    // Pick the root in [0, 1] (the cell coordinate), not the abs-min root —
    // see pick_unit_interval_root() above and #1301 for context.
    xhat = pick_unit_interval_root(N, xhat0, xhat1, xhat2);

    // Unpack xhat into actual value
    // xhat = (x-x_{i})/(x_{i+1}-x_{i})
    val = xhat * (table.xvec[i + 1] - table.xvec[i]) + table.xvec[i];

    // Cache the output value calculated
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
void CoolProp::BicubicBackend::invert_single_phase_y(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                                     parameters other_key, double other, double x, std::size_t i, std::size_t j) {
    // Get the cell
    const CellCoeffs& cell = coeffs[i][j];

    // Get the alpha coefficients
    const std::vector<double>& alpha = cell.get(other_key);

    // Normalized value in the range (0, 1)
    double xhat = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);

    double x_0 = 1, x_1 = xhat, x_2 = xhat * xhat, x_3 = xhat * xhat * xhat;

    double a = alpha[0 + 3 * 4] * x_0 + alpha[1 + 3 * 4] * x_1 + alpha[2 + 3 * 4] * x_2 + alpha[3 + 3 * 4] * x_3;          // factors of yhat^3 (m= 3)
    double b = alpha[0 + 2 * 4] * x_0 + alpha[1 + 2 * 4] * x_1 + alpha[2 + 2 * 4] * x_2 + alpha[3 + 2 * 4] * x_3;          // factors of yhat^2
    double c = alpha[0 + 1 * 4] * x_0 + alpha[1 + 1 * 4] * x_1 + alpha[2 + 1 * 4] * x_2 + alpha[3 + 1 * 4] * x_3;          // factors of yhat
    double d = alpha[0 + 0 * 4] * x_0 + alpha[1 + 0 * 4] * x_1 + alpha[2 + 0 * 4] * x_2 + alpha[3 + 0 * 4] * x_3 - other;  // constant factors
    int N = 0;
    double yhat0 = NAN, yhat1 = NAN, yhat2 = NAN, val = NAN, yhat = _HUGE;
    solve_cubic(a, b, c, d, N, yhat0, yhat1, yhat2);
    if (N == 0) {
        throw ValueError("Could not find a solution in invert_single_phase_y");
    }
    // Pick the root in [0, 1] (the cell coordinate), not the abs-min root —
    // see pick_unit_interval_root() above and #1301 for context.
    yhat = pick_unit_interval_root(N, yhat0, yhat1, yhat2);

    // Unpack xhat into actual value
    // yhat = (y-y_{j})/(y_{j+1}-y_{j})
    val = yhat * (table.yvec[j + 1] - table.yvec[j]) + table.yvec[j];

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
