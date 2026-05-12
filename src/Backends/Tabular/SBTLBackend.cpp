#if !defined(NO_TABULAR_BACKENDS)

#    include "SBTLBackend.h"
#    include "DataStructures.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include <iostream>
#    include <map>
#    include <memory>
#    include <Eigen/Sparse>
#    include <Eigen/SparseLU>

namespace CoolProp {

namespace {  // Cubic B-spline math (Kunick-style C^2 splines).

// 1D cubic B-spline coefficients with NOT-A-KNOT boundary condition.
//
// Standard interpolation conditions: (C[i-1] + 4 C[i] + C[i+1]) / 6 = f[i].
// Not-a-knot BC at the second-from-edge knot: f'''(x_1^-) = f'''(x_1^+).
// In B-spline-coefficient form this gives:
//   C[-1] = 4 C[0] - 6 C[1] + 4 C[2] - C[3]    (derived from f''' continuity)
// Substituted into the i=0 interpolation row yields:
//   8 C[0] - 5 C[1] + 4 C[2] - C[3] = 6 f[0]
// Symmetric at the other end.  Result: a banded N×N system that's tridiag
// in the interior and has 4-band rows at row 0 and row N-1.  Solved here
// via Eigen::SparseLU; for uniform N=200 grid this is ~ms one-time per row
// per property.
//
// Not-a-knot is the SciPy default for cubic splines and produces much less
// boundary error than the f''=0 natural BC for data with non-zero curvature
// at the table edges (e.g. cold-isobar enthalpy on the LIQUID region's
// xnorm=0 boundary, which has finite ∂²h/∂P²).
std::vector<double> bspline_coeffs_1d(const std::vector<double>& f) {
    const std::size_t N = f.size();
    std::vector<double> C(N, 0.0);
    if (N == 0) return C;
    if (N <= 4) {
        // Tiny grids: fall back to natural BC for stability (not-a-knot
        // requires N ≥ 5 since the 4-band boundary rows touch C[0..3]).
        for (std::size_t i = 0; i < N; ++i)
            C[i] = f[i];
        if (N == 3) C[1] = (6.0 * f[1] - f[0] - f[2]) / 4.0;
        return C;
    }

    Eigen::SparseMatrix<double> A(static_cast<int>(N), static_cast<int>(N));
    Eigen::VectorXd b(N);
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(N * 4);

    // Row 0: not-a-knot at left edge.
    trip.emplace_back(0, 0, 8.0);
    trip.emplace_back(0, 1, -5.0);
    trip.emplace_back(0, 2, 4.0);
    trip.emplace_back(0, 3, -1.0);
    b[0] = 6.0 * f[0];

    // Interior rows: standard tridiag (1, 4, 1).
    for (std::size_t i = 1; i + 1 < N; ++i) {
        const int ii = static_cast<int>(i);
        trip.emplace_back(ii, ii - 1, 1.0);
        trip.emplace_back(ii, ii, 4.0);
        trip.emplace_back(ii, ii + 1, 1.0);
        b[ii] = 6.0 * f[i];
    }

    // Row N-1: not-a-knot at right edge.
    const int n = static_cast<int>(N) - 1;
    trip.emplace_back(n, n - 3, -1.0);
    trip.emplace_back(n, n - 2, 4.0);
    trip.emplace_back(n, n - 1, -5.0);
    trip.emplace_back(n, n, 8.0);
    b[n] = 6.0 * f[N - 1];

    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    if (solver.info() != Eigen::Success) {
        // Numerical failure: fall back to a passable default.
        for (std::size_t i = 0; i < N; ++i)
            C[i] = f[i];
        return C;
    }
    Eigen::VectorXd x = solver.solve(b);
    for (std::size_t i = 0; i < N; ++i)
        C[i] = x[static_cast<int>(i)];
    return C;
}

// Fill non-finite entries in a 1D vector by nearest-valid extrapolation.
// Returns false if there's no valid entry at all (caller writes the row
// off as NaN).
bool fill_holes_1d(std::vector<double>& v) {
    const std::size_t N = v.size();
    bool any = false;
    for (std::size_t i = 0; i < N; ++i)
        if (ValidNumber(v[i])) {
            any = true;
            break;
        }
    if (!any) return false;
    // Forward pass: copy nearest left-side valid value into NaN slots.
    double last = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < N; ++i) {
        if (ValidNumber(v[i]))
            last = v[i];
        else if (ValidNumber(last))
            v[i] = last;
    }
    // Reverse pass for any slots that came before the first valid entry.
    last = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < N; ++k) {
        std::size_t i = N - 1 - k;
        if (ValidNumber(v[i]))
            last = v[i];
        else if (ValidNumber(last))
            v[i] = last;
    }
    return true;
}

// Build the not-a-knot cubic B-spline matrix for a given N and cache the
// LU factorization. The matrix only depends on N, so we factor it once per
// N and reuse for every row/column (Nx + Ny solves per property).
class BsplineMatrixCache
{
   public:
    using Solver = Eigen::SparseLU<Eigen::SparseMatrix<double>>;
    Solver& solver_for(std::size_t N) {
        auto it = solvers_.find(N);
        if (it != solvers_.end()) return *it->second;
        Eigen::SparseMatrix<double> A(static_cast<int>(N), static_cast<int>(N));
        std::vector<Eigen::Triplet<double>> trip;
        trip.reserve(N * 4);
        // Not-a-knot row 0
        trip.emplace_back(0, 0, 8.0);
        trip.emplace_back(0, 1, -5.0);
        trip.emplace_back(0, 2, 4.0);
        trip.emplace_back(0, 3, -1.0);
        for (std::size_t i = 1; i + 1 < N; ++i) {
            const int ii = static_cast<int>(i);
            trip.emplace_back(ii, ii - 1, 1.0);
            trip.emplace_back(ii, ii, 4.0);
            trip.emplace_back(ii, ii + 1, 1.0);
        }
        const int n = static_cast<int>(N) - 1;
        trip.emplace_back(n, n - 3, -1.0);
        trip.emplace_back(n, n - 2, 4.0);
        trip.emplace_back(n, n - 1, -5.0);
        trip.emplace_back(n, n, 8.0);
        A.setFromTriplets(trip.begin(), trip.end());
        A.makeCompressed();
        auto solver = std::make_unique<Solver>();
        solver->analyzePattern(A);
        solver->factorize(A);
        Solver& ref = *solver;
        solvers_[N] = std::move(solver);
        return ref;
    }

   private:
    std::map<std::size_t, std::unique_ptr<Solver>> solvers_;
};

// File-local matrix cache — initialized on first call, freed at program exit.
BsplineMatrixCache& bspline_cache() {
    static BsplineMatrixCache cache;
    return cache;
}

// Fast 1D B-spline solve using the cached LU factorization.  Equivalent to
// bspline_coeffs_1d but ~100x faster when called many times with the same N.
std::vector<double> bspline_coeffs_1d_fast(const std::vector<double>& f) {
    const std::size_t N = f.size();
    if (N <= 4) return bspline_coeffs_1d(f);
    Eigen::VectorXd b(N);
    for (std::size_t i = 0; i < N; ++i)
        b[static_cast<int>(i)] = 6.0 * f[i];
    auto& solver = bspline_cache().solver_for(N);
    if (solver.info() != Eigen::Success) return bspline_coeffs_1d(f);
    Eigen::VectorXd x = solver.solve(b);
    std::vector<double> C(N);
    for (std::size_t i = 0; i < N; ++i)
        C[i] = x[static_cast<int>(i)];
    return C;
}

// 2D tensor-product natural cubic B-spline coefficient grid.  Two passes:
// first along x for each j, then along y for each i.  Non-finite entries in
// f are filled by nearest-valid extrapolation in each direction so a
// scattered hole doesn't contaminate the entire C grid through the two
// successive 1D solves.  Uses the cached LU factorization for ~100x build
// speedup over the per-row factor-and-solve.
std::vector<std::vector<double>> bspline_coeffs_2d(const std::vector<std::vector<double>>& f) {
    const std::size_t Nx = f.size();
    const std::size_t Ny = (Nx > 0) ? f[0].size() : 0;
    std::vector<std::vector<double>> g(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> C(Nx, std::vector<double>(Ny, 0.0));
    if (Nx == 0 || Ny == 0) return C;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t j = 0; j < Ny; ++j) {
        std::vector<double> col(Nx);
        for (std::size_t i = 0; i < Nx; ++i)
            col[i] = f[i][j];
        if (!fill_holes_1d(col)) {
            for (std::size_t i = 0; i < Nx; ++i)
                g[i][j] = nan;
            continue;
        }
        auto coeffs = bspline_coeffs_1d_fast(col);
        for (std::size_t i = 0; i < Nx; ++i)
            g[i][j] = coeffs[i];
    }
    for (std::size_t i = 0; i < Nx; ++i) {
        if (!fill_holes_1d(g[i])) {
            for (std::size_t j = 0; j < Ny; ++j)
                C[i][j] = nan;
            continue;
        }
        auto coeffs = bspline_coeffs_1d_fast(g[i]);
        C[i] = std::move(coeffs);
    }
    return C;
}

// B-spline basis polynomials in a unit cell.  Index a maps to physical
// stencil offset {-1, 0, 1, 2}.  [a][m] = coeff of t^m in B_a(t).
const double Bspline_poly[4][4] = {
  {1.0 / 6.0, -1.0 / 2.0, 1.0 / 2.0, -1.0 / 6.0},
  {2.0 / 3.0, 0.0, -1.0, 1.0 / 2.0},
  {1.0 / 6.0, 1.0 / 2.0, 1.0 / 2.0, -1.0 / 2.0},
  {0.0, 0.0, 0.0, 1.0 / 6.0},
};

// Not-a-knot ghost-cell extension.  At the left edge:
//   C[-1] = 4 C[0] - 6 C[1] + 4 C[2] - C[3]
// (derived from f'''(x_1^-) = f'''(x_1^+) in B-spline coefficient form.)
// Symmetric at the right edge.  This MUST match the BC used in
// bspline_coeffs_1d for consistency.
double C_ghost(const std::vector<std::vector<double>>& C, int ic, int jc) {
    const int Nx = static_cast<int>(C.size());
    const int Ny = (Nx > 0) ? static_cast<int>(C[0].size()) : 0;
    if (Nx == 0 || Ny == 0) return 0.0;
    auto inx = [&](int k) { return k >= 0 && k < Nx; };
    auto iny = [&](int k) { return k >= 0 && k < Ny; };
    auto get_x = [&](int ix, int iy) -> double {
        if (inx(ix) && iny(iy)) return C[ix][iy];
        if (iny(iy)) {
            if (ix == -1) {
                if (Nx >= 4) return 4.0 * C[0][iy] - 6.0 * C[1][iy] + 4.0 * C[2][iy] - C[3][iy];
                return 2.0 * C[0][iy] - C[1][iy];  // small-N fallback
            }
            if (ix == Nx) {
                if (Nx >= 4) return 4.0 * C[Nx - 1][iy] - 6.0 * C[Nx - 2][iy] + 4.0 * C[Nx - 3][iy] - C[Nx - 4][iy];
                return 2.0 * C[Nx - 1][iy] - C[Nx - 2][iy];
            }
        }
        return 0.0;
    };
    if (iny(jc)) return get_x(ic, jc);
    if (jc == -1) {
        if (Ny >= 4) return 4.0 * get_x(ic, 0) - 6.0 * get_x(ic, 1) + 4.0 * get_x(ic, 2) - get_x(ic, 3);
        return 2.0 * get_x(ic, 0) - get_x(ic, 1);
    }
    if (jc == Ny) {
        if (Ny >= 4) return 4.0 * get_x(ic, Ny - 1) - 6.0 * get_x(ic, Ny - 2) + 4.0 * get_x(ic, Ny - 3) - get_x(ic, Ny - 4);
        return 2.0 * get_x(ic, Ny - 1) - get_x(ic, Ny - 2);
    }
    return 0.0;
}

std::vector<double> bspline_polynomial_coeffs(const std::vector<std::vector<double>>& C, std::size_t i, std::size_t j) {
    std::vector<double> alpha(16, 0.0);
    for (int a = 0; a < 4; ++a) {
        const int ic = static_cast<int>(i) + a - 1;
        for (int b = 0; b < 4; ++b) {
            const int jc = static_cast<int>(j) + b - 1;
            const double Cab = C_ghost(C, ic, jc);
            if (!ValidNumber(Cab)) continue;
            for (int m = 0; m < 4; ++m) {
                if (Bspline_poly[a][m] == 0.0) continue;
                for (int n = 0; n < 4; ++n) {
                    alpha[m * 4 + n] += Cab * Bspline_poly[a][m] * Bspline_poly[b][n];
                }
            }
        }
    }
    return alpha;
}

std::vector<double> hermite_bicubic_polynomial_coeffs_impl(double f00, double f10, double f01, double f11, double fx00, double fx10, double fx01,
                                                           double fx11, double fy00, double fy10, double fy01, double fy11, double fxy00,
                                                           double fxy10, double fxy01, double fxy11) {
    std::vector<double> alpha(16, 0.0);
    // c_{0,0..3}: value, fy, value+fy terms at xi=0
    alpha[0 * 4 + 0] = f00;
    alpha[0 * 4 + 1] = fy00;
    alpha[0 * 4 + 2] = -3.0 * f00 + 3.0 * f01 - 2.0 * fy00 - fy01;
    alpha[0 * 4 + 3] = 2.0 * f00 - 2.0 * f01 + fy00 + fy01;
    // c_{1,0..3}: fx, fxy, mixed at xi=0
    alpha[1 * 4 + 0] = fx00;
    alpha[1 * 4 + 1] = fxy00;
    alpha[1 * 4 + 2] = -3.0 * fx00 + 3.0 * fx01 - 2.0 * fxy00 - fxy01;
    alpha[1 * 4 + 3] = 2.0 * fx00 - 2.0 * fx01 + fxy00 + fxy01;
    // c_{2,0..3}: order-2 in xi
    alpha[2 * 4 + 0] = -3.0 * f00 + 3.0 * f10 - 2.0 * fx00 - fx10;
    alpha[2 * 4 + 1] = -3.0 * fy00 + 3.0 * fy10 - 2.0 * fxy00 - fxy10;
    alpha[2 * 4 + 2] = 9.0 * f00 - 9.0 * f10 - 9.0 * f01 + 9.0 * f11 + 6.0 * fx00 + 3.0 * fx10 - 6.0 * fx01 - 3.0 * fx11 + 6.0 * fy00 - 6.0 * fy10
                       + 3.0 * fy01 - 3.0 * fy11 + 4.0 * fxy00 + 2.0 * fxy10 + 2.0 * fxy01 + fxy11;
    alpha[2 * 4 + 3] = -6.0 * f00 + 6.0 * f10 + 6.0 * f01 - 6.0 * f11 - 4.0 * fx00 - 2.0 * fx10 + 4.0 * fx01 + 2.0 * fx11 - 3.0 * fy00 + 3.0 * fy10
                       - 3.0 * fy01 + 3.0 * fy11 - 2.0 * fxy00 - fxy10 - 2.0 * fxy01 - fxy11;
    // c_{3,0..3}: order-3 in xi
    alpha[3 * 4 + 0] = 2.0 * f00 - 2.0 * f10 + fx00 + fx10;
    alpha[3 * 4 + 1] = 2.0 * fy00 - 2.0 * fy10 + fxy00 + fxy10;
    alpha[3 * 4 + 2] = -6.0 * f00 + 6.0 * f10 + 6.0 * f01 - 6.0 * f11 - 3.0 * fx00 - 3.0 * fx10 + 3.0 * fx01 + 3.0 * fx11 - 4.0 * fy00 + 4.0 * fy10
                       - 2.0 * fy01 + 2.0 * fy11 - 2.0 * fxy00 - 2.0 * fxy10 - fxy01 - fxy11;
    alpha[3 * 4 + 3] = 4.0 * f00 - 4.0 * f10 - 4.0 * f01 + 4.0 * f11 + 2.0 * fx00 + 2.0 * fx10 - 2.0 * fx01 - 2.0 * fx11 + 2.0 * fy00 - 2.0 * fy10
                       + 2.0 * fy01 - 2.0 * fy11 + fxy00 + fxy10 + fxy01 + fxy11;
    return alpha;
}

double eval_alpha(const std::vector<double>& alpha, double tx, double ty) {
    const double tx2 = tx * tx, tx3 = tx2 * tx;
    const double ty2 = ty * ty, ty3 = ty2 * ty;
    const double B0 = alpha[0] + alpha[1] * ty + alpha[2] * ty2 + alpha[3] * ty3;
    const double B1 = alpha[4] + alpha[5] * ty + alpha[6] * ty2 + alpha[7] * ty3;
    const double B2 = alpha[8] + alpha[9] * ty + alpha[10] * ty2 + alpha[11] * ty3;
    const double B3 = alpha[12] + alpha[13] * ty + alpha[14] * ty2 + alpha[15] * ty3;
    return B0 + B1 * tx + B2 * tx2 + B3 * tx3;
}

}  // namespace

// Hermite bicubic on a unit cell.  Takes 16 corner data (4 values, 4
// ∂f/∂xi, 4 ∂f/∂eta, 4 ∂²f/∂xi∂eta — derivatives must already be
// scaled to the unit-cell coordinates xi, eta ∈ [0,1]) and returns the
// 16-element alpha vector in SBTL convention: alpha[4m + n] = c_{m,n}
// such that f(xi, eta) = Σ_{m,n} c_{m,n} xi^m eta^n.
//
// Corner indexing: 00 = (xi=0,eta=0), 10 = (1,0), 01 = (0,1), 11 = (1,1).
// This is the standard bicubic interpolation closed form (see
// en.wikipedia.org/wiki/Bicubic_interpolation; coefficients match the
// Ainv matrix used by BicubicBackend in TabularBackends.cpp, transposed
// to SBTL convention so the same evaluate_single_phase polynomial can
// consume the alpha vector unchanged).
//
// First building block toward CoolProp-foi.1 (Hermite bicubic on the
// normph backbone).  Not wired into the build path yet; callable from
// tests and the upcoming iapws_conformance_mode path.
std::vector<double> hermite_bicubic_polynomial_coeffs(double f00, double f10, double f01, double f11, double fx00, double fx10, double fx01,
                                                      double fx11, double fy00, double fy10, double fy01, double fy11, double fxy00, double fxy10,
                                                      double fxy01, double fxy11) {
    return hermite_bicubic_polynomial_coeffs_impl(f00, f10, f01, f11, fx00, fx10, fx01, fx11, fy00, fy10, fy01, fy11, fxy00, fxy10, fxy01, fxy11);
}

// ---------------------------------------------------------------------------
// Helper: get pointer to a table field by parameter key
// ---------------------------------------------------------------------------
static std::vector<std::vector<double>>* sbtl_get_field(SinglePhaseGriddedTableData& table, parameters param) {
    switch (param) {
        case iT:
            return &table.T;
        case iP:
            return &table.p;
        case iDmolar:
            return &table.rhomolar;
        case iSmolar:
            return &table.smolar;
        case iHmolar:
            return &table.hmolar;
        case iUmolar:
            return &table.umolar;
        case ispeed_sound:
            return &table.speed_sound;
        case iviscosity:
            return &table.visc;
        case iconductivity:
            return &table.cond;
        default:
            throw ValueError("Invalid param in sbtl_get_field");
    }
}

// ---------------------------------------------------------------------------
// NormalizedPHTable: coordinate-aligned PH table.
// ---------------------------------------------------------------------------

void NormalizedPHTable::set_limits() {
    if (!AS) throw ValueError("NormalizedPHTable: AS is not set");

    // x-axis is xnorm ∈ [0, 1] regardless of region.
    xmin = 0.0;
    xmax = 1.0;

    // y-axis: pressure range depends on region.
    const double p_crit = static_cast<double>(AS->p_critical());
    const double p_min = static_cast<double>(AS->p_triple());
    const double p_max = static_cast<double>(AS->pmax());
    if (region_ == SUPER) {
        ymin = p_crit;
        ymax = p_max;
    } else {
        // LIQUID and VAPOR: from p_triple up to (just below) p_crit.
        ymin = p_min;
        // Stay 0.1% off p_crit so saturation_hmolar_*(p) stays in the well-
        // conditioned region of the cached sat curve / superancillary.
        ymax = p_crit * 0.999;
    }
}

double NormalizedPHTable::xnorm_from_h(double h, double P, double h_sat) const {
    // Source of truth per region:
    //   LIQUID: h_lo from Cheb (h(T_min, p)),     h_hi from H-superanc (h_sat,L)
    //   VAPOR : h_lo from H-superanc (h_sat,V),    h_hi from Cheb (h(T_max, p))
    //   SUPER : both from Cheb
    // Caller passes h_sat (machine-precision via H-superancillary) for the
    // sat-boundary; if NaN, fall back to the Cheb (less accurate near the
    // critical cusp, but still works for smooth-monotonic stretches).
    auto cheb_or_isobar = [&](const Cheb1D& cb, const std::vector<double>& fallback) -> double {
        if (cb.valid()) return cb.eval(P);
        if (fallback.empty()) return std::numeric_limits<double>::quiet_NaN();
        const std::size_t Nrows = yvec.size();
        std::size_t j = 0;
        bisect_vector(yvec, P, j);
        if (j >= Nrows - 1) j = Nrows - 2;
        const double frac = (std::log(P) - std::log(yvec[j])) / (std::log(yvec[j + 1]) - std::log(yvec[j]));
        return fallback[j] + frac * (fallback[j + 1] - fallback[j]);
    };
    double h_lo, h_hi;
    switch (region_) {
        case LIQUID:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = std::isfinite(h_sat) ? h_sat : cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case VAPOR:
            h_lo = std::isfinite(h_sat) ? h_sat : cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case SUPER:
        default:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
    }
    return (h - h_lo) / (h_hi - h_lo);
}

void SBTLBackend::populate_normph_bounds(NormalizedPHTable& table) {
    // Sample yvec must already be laid down (set_limits has run).  We compute
    // h_lo(P) and h_hi(P) at each row and store; subsequent xnorm_from_h
    // queries linearly interpolate in log-P.
    if (table.yvec.empty()) {
        throw ValueError("populate_normph_bounds: yvec not yet populated; call set_limits / resize first");
    }
    const std::size_t Ny = table.yvec.size();
    table.h_lo_isobar.assign(Ny, _HUGE);
    table.h_hi_isobar.assign(Ny, _HUGE);

    const double T_min = std::max(static_cast<double>(this->AS->Ttriple()), static_cast<double>(this->AS->Tmin()));
    const double T_max_ext = static_cast<double>(this->AS->Tmax()) * 1.499;

    auto safe_h_at_PT = [this](double P, double T) -> double {
        try {
            this->AS->update(PT_INPUTS, P, T);
            return static_cast<double>(this->AS->hmolar());
        } catch (...) {
            return _HUGE;
        }
    };
    // Melting-line lower bound on T per isobar.  Returns T_melt(P) directly
    // (no safety margin) when the fluid has a melting ancillary; -inf
    // otherwise.  Using exactly T_melt keeps the maximum useful liquid
    // range — important for cryogens like helium where Tc - Ttriple ~ 3K
    // and even a 1K margin would lose 30% of the range.  If the EOS
    // happens to throw exactly at T = T_melt, fill_holes_1d will substitute
    // the nearest interior valid value during spline build.
    const bool fluid_has_melt = this->AS && this->AS->has_melting_line();
    auto T_melt_floor = [this, fluid_has_melt](double P) -> double {
        if (!fluid_has_melt) return -std::numeric_limits<double>::infinity();
        try {
            return static_cast<double>(this->AS->melting_line(iT, iP, P));
        } catch (...) {
            return -std::numeric_limits<double>::infinity();
        }
    };
    auto safe_h_sat_L = [this](double P) -> double {
        try {
            return saturation_hmolar_liquid(P);
        } catch (...) {
            return _HUGE;
        }
    };
    auto safe_h_sat_V = [this](double P) -> double {
        try {
            return saturation_hmolar_vapor(P);
        } catch (...) {
            return _HUGE;
        }
    };

    for (std::size_t j = 0; j < Ny; ++j) {
        const double P = table.yvec[j];
        switch (table.region()) {
            case NormalizedPHTable::LIQUID: {
                const double h_sat_L = safe_h_sat_L(P);
                table.h_hi_isobar[j] = h_sat_L;
                // T_lo = T_melt(P) when the fluid has a melting ancillary,
                // else T_min.  The Helmholtz EOS is valid down to the
                // melting curve — including water's negative-slope ice I
                // band where T_melt < T_triple between ~5 and ~200 MPa.
                const double T_melt = T_melt_floor(P);
                const double T_lo = std::isfinite(T_melt) ? T_melt : T_min;
                table.h_lo_isobar[j] = safe_h_at_PT(P, T_lo);
                break;
            }
            case NormalizedPHTable::VAPOR: {
                table.h_lo_isobar[j] = safe_h_sat_V(P);
                table.h_hi_isobar[j] = safe_h_at_PT(P, T_max_ext);
                break;
            }
            case NormalizedPHTable::SUPER: {
                // SUPER: P >= p_crit, T from Tcrit to T_max_ext.  Using
                // T_lo = T_crit makes h_lo = h(Tcrit, P) — the critical
                // isotherm — which is a smooth function of P for P >
                // pcrit and avoids both (a) the freezing curve (T_min
                // T_lo = T_melt(P) when the fluid has a melting ancillary,
                // else T_min.  Uses the melting curve as the authoritative
                // lower-T bound — including when it dips BELOW T_triple
                // (water's negative-slope ice I region, ~5-200 MPa, where
                // T_melt < T_triple).  The earlier (1.05*Tcrit) excluded
                // the entire supercritical-liquid wing (T<T_crit, P>p_crit)
                // — for CO2 ~100K of useful range below T_crit at high P.
                const double T_melt = T_melt_floor(P);
                const double T_lo = std::isfinite(T_melt) ? T_melt : T_min;
                table.h_lo_isobar[j] = safe_h_at_PT(P, T_lo);
                table.h_hi_isobar[j] = safe_h_at_PT(P, T_max_ext);
                break;
            }
        }
    }

    // Build Chebyshev expansions of h_lo(p) and h_hi(p) over the table's
    // pressure range, sampled at Chebyshev–Lobatto nodes via direct EOS
    // queries.  Spectrally-accurate for the smooth boundary functions
    // and eliminates the linear-in-log-p xnorm_from_h interpolation
    // error that previously bottlenecked rho accuracy near steep features.
    //
    // The Cheb1D becomes the source of truth: after building it, we
    // OVERWRITE h_lo_isobar/h_hi_isobar with Cheb1D-evaluated values at
    // the yvec grid points.  This keeps cell-build (which samples h at
    // h_lo_isobar + xnorm * (h_hi_isobar - h_lo_isobar)) consistent with
    // cell-lookup (which uses Cheb1D to derive xnorm from h).  Without
    // this consistency, even 1e-12 discrepancy in xnorm at yvec nodes
    // makes cell polynomials overshoot.
    // Piecewise Chebyshev breakpoints.  Subdivision concentrates
    // resolution where the function changes shape — for LIQUID/VAPOR
    // h_sat boundaries this means a fine piece near p_crit where the
    // sqrt-cusp lives.  For SUPER, one piece is sufficient (smooth
    // throughout, h_lo/h_hi just functions of T_min/T_max isotherms).
    // Nudge p_lo upward by a small relative amount so the lowest
    // Lobatto node — exp(log(p_lo)) plus floating-point roundoff —
    // doesn't slip below the EOS / melting-line ancillary lower bound
    // (water's melting curve has p_triple as its strict lower bound;
    // queries at p < p_triple throw).  1e-4 relative is well below
    // tabular discretization error and keeps the Cheb fit clean.
    const double p_lo = table.yvec.front() * 1.0001;
    const double p_hi = table.yvec.back() * 0.9999;
    // Build Chebyshev expansions only for the NON-saturation boundaries:
    //   LIQUID: h_lo = h(T_min, p)        — Cheb
    //   LIQUID: h_hi = h_sat,L(p)         — H-superancillary at lookup
    //   VAPOR : h_lo = h_sat,V(p)         — H-superancillary at lookup
    //   VAPOR : h_hi = h(T_max_ext, p)    — Cheb
    //   SUPER : both                      — Cheb
    // Saturation enthalpies are already machine-precision via the
    // H-superancillaries (Bell et al.); re-fitting with Chebyshev would
    // just inherit any error there with no upside.
    std::vector<double> breakpoints;
    if (table.region() == NormalizedPHTable::SUPER) {
        breakpoints = {p_lo, p_hi};
    } else if (this->AS) {
        const double p_crit = static_cast<double>(this->AS->p_critical());
        breakpoints = {p_lo, 0.5 * p_crit, 0.85 * p_crit, 0.97 * p_crit, p_hi};
    } else {
        breakpoints = {p_lo, p_hi};
    }
    const std::size_t N_cheb = 32;
    auto h_lo_fn = [&](double P) -> double {
        // Only LIQUID and SUPER need a Cheb for h_lo (h(T_min, p)).
        const double T_melt_p = T_melt_floor(P);
        const double T_lo_p = std::isfinite(T_melt_p) ? T_melt_p : T_min;
        return safe_h_at_PT(P, T_lo_p);
    };
    auto h_hi_fn = [&](double P) -> double {
        // Only VAPOR and SUPER need a Cheb for h_hi (h(T_max_ext, p)).
        return safe_h_at_PT(P, T_max_ext);
    };
    try {
        if (table.region() != NormalizedPHTable::VAPOR) {
            table.h_lo_cheb = Cheb1D::build(breakpoints, N_cheb, h_lo_fn);
        }
        if (table.region() != NormalizedPHTable::LIQUID) {
            table.h_hi_cheb = Cheb1D::build(breakpoints, N_cheb, h_hi_fn);
        }
        // Reconcile the per-row arrays with the per-region source-of-truth
        // so cell-build at yvec[j] and cell-lookup agree on h_lo/h_hi:
        // sat boundaries via H-superancillary, others via Cheb1D.
        for (std::size_t j = 0; j < Ny; ++j) {
            const double Pj = table.yvec[j];
            switch (table.region()) {
                case NormalizedPHTable::LIQUID:
                    table.h_lo_isobar[j] = table.h_lo_cheb.eval(Pj);
                    table.h_hi_isobar[j] = safe_h_sat_L(Pj);
                    break;
                case NormalizedPHTable::VAPOR:
                    table.h_lo_isobar[j] = safe_h_sat_V(Pj);
                    table.h_hi_isobar[j] = table.h_hi_cheb.eval(Pj);
                    break;
                case NormalizedPHTable::SUPER:
                    table.h_lo_isobar[j] = table.h_lo_cheb.eval(Pj);
                    table.h_hi_isobar[j] = table.h_hi_cheb.eval(Pj);
                    break;
            }
        }
    } catch (...) {
        table.h_lo_cheb = Cheb1D{};
        table.h_hi_cheb = Cheb1D{};
    }
}

void SBTLBackend::build_normph_table(NormalizedPHTable& table) {
    if (!this->AS) throw ValueError("build_normph_table: AS is not set");
    table.AS = this->AS;
    table.set_limits();
    table.resize(table.Nx, table.Ny);  // allocates LIST_OF_MATRICES storage; fills xvec, yvec.
    populate_normph_bounds(table);

    // For LIQUID/VAPOR regions some isobars may sit at pressures where the
    // H-superancillary or sat-cache isn't well-defined (e.g. very close to
    // p_crit on the dome).  Skip rows where bounds are non-finite or
    // collapsed.
    for (std::size_t j = 0; j < table.Ny; ++j) {
        const double h_lo = table.h_lo_isobar[j];
        const double h_hi = table.h_hi_isobar[j];
        if (!ValidNumber(h_lo) || !ValidNumber(h_hi) || !(h_hi > h_lo)) {
            // Leave whole isobar as holes.
            continue;
        }
        const double P = table.yvec[j];
        for (std::size_t i = 0; i < table.Nx; ++i) {
            const double xnorm = table.xvec[i];  // ∈ [0, 1] linearly
            const double h = h_lo + xnorm * (h_hi - h_lo);
            // Saturation boundaries: xnorm=1 of LIQUID (h=h_sat,L) and
            // xnorm=0 of VAPOR (h=h_sat,V).  HmolarP_INPUTS at exactly the
            // sat enthalpy is numerically ambiguous (Q can be 0.999.. or
            // 1.0001 depending on flash residual), and the resulting state
            // may flip between sat-side and slightly-2phase.  Force the
            // saturated-side properties via PQ_INPUTS at these boundaries
            // — that gives us rho_sat, T_sat etc. exactly, which is what
            // the routing-time xnorm = 0/1 query SHOULD see.
            const bool sat_boundary =
              (table.region() == NormalizedPHTable::LIQUID && i == table.Nx - 1) || (table.region() == NormalizedPHTable::VAPOR && i == 0);
            try {
                if (sat_boundary) {
                    const double Q_target = (table.region() == NormalizedPHTable::VAPOR) ? 1.0 : 0.0;
                    this->AS->update(PQ_INPUTS, P, Q_target);
                } else {
                    this->AS->update(HmolarP_INPUTS, h, P);
                }
                if (!ValidNumber(this->AS->rhomolar())) {
                    continue;
                }
            } catch (...) {
                continue;
            }
            // Skip two-phase interior states.  Allow Q=0/1 (saturated side
            // accessed via PQ at the boundaries); skip strictly-interior 2-phase.
            const double q = static_cast<double>(this->AS->Q());
            if (!sat_boundary && q > 0.0 && q < 1.0) {
                continue;
            }
            table.T[i][j] = static_cast<double>(this->AS->T());
            table.p[i][j] = static_cast<double>(this->AS->p());
            table.rhomolar[i][j] = static_cast<double>(this->AS->rhomolar());
            table.hmolar[i][j] = static_cast<double>(this->AS->hmolar());
            table.smolar[i][j] = static_cast<double>(this->AS->smolar());
            table.umolar[i][j] = static_cast<double>(this->AS->umolar());
            try {
                table.speed_sound[i][j] = static_cast<double>(this->AS->speed_sound());
            } catch (...) {
                // speed_sound may not be defined in some regions (e.g. inside
                // the metastable spinodal); leave the cell as a _HUGE hole
                // and build_bspline_coeffs will mark cells with non-finite
                // corners as invalid.
            }
            try {
                table.visc[i][j] = static_cast<double>(this->AS->viscosity());
                table.cond[i][j] = static_cast<double>(this->AS->conductivity());
            } catch (...) {
                // Transport may fail in some regions; leave as holes.
            }
        }
    }
}

void SBTLBackend::build_normph_tables() {
    if (normph_tables_built) return;
    build_normph_table(_normph_liquid);
    build_normph_table(_normph_vapor);
    build_normph_table(_normph_super);
    build_bspline_coeffs(_normph_liquid, _coeffs_normph_liquid);
    build_bspline_coeffs(_normph_vapor, _coeffs_normph_vapor);
    build_bspline_coeffs(_normph_super, _coeffs_normph_super);
    normph_tables_built = true;
}

// ---------------------------------------------------------------------------
// NormalizedPTTable: coordinate-aligned PT table.  Mirror of NormalizedPHTable
// with iT in place of iHmolar.  Cell boundary at xnorm=1 (LIQUID) is exactly
// the saturation curve, so cells never straddle the dome — this is the fix
// for the 5 % compressed-liquid PT residue documented in
// dev/sbtl_pt_outstanding_work.md.
// ---------------------------------------------------------------------------

void NormalizedPTTable::set_limits() {
    if (!AS) throw ValueError("NormalizedPTTable: AS is not set");
    xmin = 0.0;
    xmax = 1.0;
    const double p_crit = static_cast<double>(AS->p_critical());
    const double p_min = static_cast<double>(AS->p_triple());
    const double p_max = static_cast<double>(AS->pmax());
    if (region_ == SUPER) {
        ymin = p_crit;
        ymax = p_max;
    } else {
        ymin = p_min;
        ymax = p_crit * 0.999;
    }
}

double NormalizedPTTable::xnorm_from_T(double T, double P, double T_sat) const {
    auto cheb_or_isobar = [&](const Cheb1D& cb, const std::vector<double>& fallback) -> double {
        if (cb.valid()) return cb.eval(P);
        if (fallback.empty()) return std::numeric_limits<double>::quiet_NaN();
        const std::size_t Nrows = yvec.size();
        std::size_t j = 0;
        bisect_vector(yvec, P, j);
        if (j >= Nrows - 1) j = Nrows - 2;
        const double frac = (std::log(P) - std::log(yvec[j])) / (std::log(yvec[j + 1]) - std::log(yvec[j]));
        return fallback[j] + frac * (fallback[j + 1] - fallback[j]);
    };
    double T_lo = NAN, T_hi = NAN;
    switch (region_) {
        case LIQUID:
            T_lo = cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = std::isfinite(T_sat) ? T_sat : cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case VAPOR:
            T_lo = std::isfinite(T_sat) ? T_sat : cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case SUPER:
        default:
            T_lo = cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
    }
    return (T - T_lo) / (T_hi - T_lo);
}

double NormalizedPTTable::T_from_xnorm(double xnorm, double P, double T_sat) const {
    auto cheb_or_isobar = [&](const Cheb1D& cb, const std::vector<double>& fallback) -> double {
        if (cb.valid()) return cb.eval(P);
        if (fallback.empty()) return std::numeric_limits<double>::quiet_NaN();
        const std::size_t Nrows = yvec.size();
        std::size_t j = 0;
        bisect_vector(yvec, P, j);
        if (j >= Nrows - 1) j = Nrows - 2;
        const double frac = (std::log(P) - std::log(yvec[j])) / (std::log(yvec[j + 1]) - std::log(yvec[j]));
        return fallback[j] + frac * (fallback[j + 1] - fallback[j]);
    };
    double T_lo = NAN, T_hi = NAN;
    switch (region_) {
        case LIQUID:
            T_lo = cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = std::isfinite(T_sat) ? T_sat : cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case VAPOR:
            T_lo = std::isfinite(T_sat) ? T_sat : cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case SUPER:
        default:
            T_lo = cheb_or_isobar(T_lo_cheb, T_lo_isobar);
            T_hi = cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
    }
    return T_lo + xnorm * (T_hi - T_lo);
}

void SBTLBackend::saturation_T_LV(double p, double& T_L, double& T_V) {
    if (_Tsat_LV_cache_p == p) {
        T_L = _Tsat_LV_cache_TL;
        T_V = _Tsat_LV_cache_TV;
        return;
    }
    // Bound check BEFORE calling AS->update: PQ_INPUTS extrapolates the sat
    // curve along the metastable branch for p < p_triple instead of throwing.
    // That leaves AS in a non-physical state and pollutes the cache with a
    // bogus T_sat.  Throw cleanly here so the caller (and any sibling backend
    // sharing AS) sees no side effect.
    const double p_min_pure = static_cast<double>(this->AS->p_triple());
    const double p_max_pure = static_cast<double>(this->AS->p_critical());
    if (!(p >= p_min_pure && p <= p_max_pure)) {
        throw ValueError(format("SBTL saturation_T_LV: P=%g Pa outside the sat range [%g, %g] Pa.", p, p_min_pure, p_max_pure));
    }
    // Pure fluids: T_sat,L = T_sat,V = T_sat(P).  A single PQ_INPUTS at Q=0
    // suffices.  Pseudo-pure / true mixtures are rejected at construction so
    // this path always runs against a single-component fluid.
    this->AS->update(PQ_INPUTS, p, 0.0);
    const double T_sat = static_cast<double>(this->AS->T());
    T_L = T_sat;
    T_V = T_sat;
    _Tsat_LV_cache_p = p;
    _Tsat_LV_cache_TL = T_L;
    _Tsat_LV_cache_TV = T_V;
}

void SBTLBackend::populate_normpt_bounds(NormalizedPTTable& table) {
    if (table.yvec.empty()) {
        throw ValueError("populate_normpt_bounds: yvec not yet populated; call set_limits / resize first");
    }
    const std::size_t Ny = table.yvec.size();
    table.T_lo_isobar.assign(Ny, _HUGE);
    table.T_hi_isobar.assign(Ny, _HUGE);

    const double T_min = std::max(static_cast<double>(this->AS->Ttriple()), static_cast<double>(this->AS->Tmin()));
    // T_max_ext: cap at the published EOS upper bound (no extrapolation).
    // The Helmholtz EOS is only guaranteed valid up to Tmax; extrapolating
    // beyond gives ~5–50 % silent rho errors at xnorm near 1 in cells whose
    // sampled corners fall outside validity.  Stay inside the EOS envelope.
    const double T_max_ext = static_cast<double>(this->AS->Tmax());
    const bool fluid_has_melt = this->AS && this->AS->has_melting_line();

    auto T_melt_floor = [this, fluid_has_melt](double P) -> double {
        if (!fluid_has_melt) return -std::numeric_limits<double>::infinity();
        try {
            return static_cast<double>(this->AS->melting_line(iT, iP, P));
        } catch (...) {
            return -std::numeric_limits<double>::infinity();
        }
    };
    auto safe_T_sat_L = [this](double P) -> double {
        try {
            this->AS->update(PQ_INPUTS, P, 0.0);
            return static_cast<double>(this->AS->T());
        } catch (...) {
            return _HUGE;
        }
    };
    auto safe_T_sat_V = [this](double P) -> double {
        try {
            this->AS->update(PQ_INPUTS, P, 1.0);
            return static_cast<double>(this->AS->T());
        } catch (...) {
            return _HUGE;
        }
    };

    for (std::size_t j = 0; j < Ny; ++j) {
        const double P = table.yvec[j];
        switch (table.region()) {
            case NormalizedPTTable::LIQUID: {
                const double T_melt = T_melt_floor(P);
                const double T_lo = std::isfinite(T_melt) ? T_melt : T_min;
                table.T_lo_isobar[j] = T_lo;
                table.T_hi_isobar[j] = safe_T_sat_L(P);
                break;
            }
            case NormalizedPTTable::VAPOR: {
                table.T_lo_isobar[j] = safe_T_sat_V(P);
                table.T_hi_isobar[j] = T_max_ext;
                break;
            }
            case NormalizedPTTable::SUPER: {
                const double T_melt = T_melt_floor(P);
                const double T_lo = std::isfinite(T_melt) ? T_melt : T_min;
                table.T_lo_isobar[j] = T_lo;
                table.T_hi_isobar[j] = T_max_ext;
                break;
            }
        }
    }

    // Cheb1D fits for the non-saturation boundaries (LIQUID T_lo, VAPOR T_hi,
    // SUPER both).  Saturation boundary T_sat,L/V comes straight from
    // PQ_INPUTS at lookup time — already machine-precision.
    const double p_lo = table.yvec.front() * 1.0001;
    const double p_hi = table.yvec.back() * 0.9999;
    std::vector<double> breakpoints;
    if (table.region() == NormalizedPTTable::SUPER) {
        breakpoints = {p_lo, p_hi};
    } else if (this->AS) {
        const double p_crit = static_cast<double>(this->AS->p_critical());
        breakpoints = {p_lo, 0.5 * p_crit, 0.85 * p_crit, 0.97 * p_crit, p_hi};
    } else {
        breakpoints = {p_lo, p_hi};
    }
    const std::size_t N_cheb = 32;
    auto T_lo_fn = [&](double P) -> double {
        const double T_melt = T_melt_floor(P);
        return std::isfinite(T_melt) ? T_melt : T_min;
    };
    auto T_hi_fn = [&](double /*P*/) -> double { return T_max_ext; };
    try {
        if (table.region() != NormalizedPTTable::VAPOR) {
            table.T_lo_cheb = Cheb1D::build(breakpoints, N_cheb, T_lo_fn);
        }
        if (table.region() != NormalizedPTTable::LIQUID) {
            table.T_hi_cheb = Cheb1D::build(breakpoints, N_cheb, T_hi_fn);
        }
        // Reconcile: cell-build samples T at T_lo_isobar + xnorm * (T_hi - T_lo);
        // cell-lookup uses Cheb1D for the same.  Source-of-truth at row j must
        // be the Cheb evaluation for the non-sat side, sat HEOS for the sat side.
        for (std::size_t j = 0; j < Ny; ++j) {
            const double Pj = table.yvec[j];
            switch (table.region()) {
                case NormalizedPTTable::LIQUID:
                    table.T_lo_isobar[j] = table.T_lo_cheb.eval(Pj);
                    // T_hi_isobar = T_sat,L stays as-is from PQ_INPUTS.
                    break;
                case NormalizedPTTable::VAPOR:
                    // T_lo_isobar = T_sat,V stays as-is.
                    table.T_hi_isobar[j] = table.T_hi_cheb.eval(Pj);
                    break;
                case NormalizedPTTable::SUPER:
                    table.T_lo_isobar[j] = table.T_lo_cheb.eval(Pj);
                    table.T_hi_isobar[j] = table.T_hi_cheb.eval(Pj);
                    break;
            }
        }
    } catch (...) {
        table.T_lo_cheb = Cheb1D{};
        table.T_hi_cheb = Cheb1D{};
    }
}

void SBTLBackend::build_normpt_table(NormalizedPTTable& table) {
    if (!this->AS) throw ValueError("build_normpt_table: AS is not set");
    table.AS = this->AS;
    table.set_limits();
    table.resize(table.Nx, table.Ny);
    populate_normpt_bounds(table);

    for (std::size_t j = 0; j < table.Ny; ++j) {
        const double T_lo = table.T_lo_isobar[j];
        const double T_hi = table.T_hi_isobar[j];
        if (!ValidNumber(T_lo) || !ValidNumber(T_hi) || !(T_hi > T_lo)) {
            continue;
        }
        const double P = table.yvec[j];
        for (std::size_t i = 0; i < table.Nx; ++i) {
            const double xnorm = table.xvec[i];
            const double T = T_lo + xnorm * (T_hi - T_lo);
            const bool sat_boundary =
              (table.region() == NormalizedPTTable::LIQUID && i == table.Nx - 1) || (table.region() == NormalizedPTTable::VAPOR && i == 0);
            try {
                if (sat_boundary) {
                    const double Q_target = (table.region() == NormalizedPTTable::VAPOR) ? 1.0 : 0.0;
                    this->AS->update(PQ_INPUTS, P, Q_target);
                } else {
                    this->AS->update(PT_INPUTS, P, T);
                }
                if (!ValidNumber(this->AS->rhomolar())) continue;
            } catch (...) {
                continue;
            }
            const double q = static_cast<double>(this->AS->Q());
            if (!sat_boundary && q > 0.0 && q < 1.0) continue;
            table.T[i][j] = static_cast<double>(this->AS->T());
            table.p[i][j] = static_cast<double>(this->AS->p());
            table.rhomolar[i][j] = static_cast<double>(this->AS->rhomolar());
            table.hmolar[i][j] = static_cast<double>(this->AS->hmolar());
            table.smolar[i][j] = static_cast<double>(this->AS->smolar());
            table.umolar[i][j] = static_cast<double>(this->AS->umolar());
            try {
                table.visc[i][j] = static_cast<double>(this->AS->viscosity());
                table.cond[i][j] = static_cast<double>(this->AS->conductivity());
            } catch (...) {
                // Transport optional.
            }
        }
    }
}

void SBTLBackend::build_normpt_tables() {
    if (normpt_tables_built) return;
    build_normpt_table(_normpt_liquid);
    build_normpt_table(_normpt_vapor);
    build_normpt_table(_normpt_super);
    build_bspline_coeffs(_normpt_liquid, _coeffs_normpt_liquid);
    build_bspline_coeffs(_normpt_vapor, _coeffs_normpt_vapor);
    build_bspline_coeffs(_normpt_super, _coeffs_normpt_super);
    normpt_tables_built = true;
}

void SBTLBackend::build_bspline_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs) {
    if (!coeffs.empty()) return;

    // Core thermodynamic params drive cell validity.  Auxiliary params
    // (speed_sound / viscosity / conductivity) build their own polynomial
    // alongside but DON'T affect cell validity — at sat-curve boundaries
    // and metastable spinodals the EOS may throw for speed_sound or
    // transport, and we don't want one missing corner to disqualify the
    // entire cell from the thermodynamic-property lookups.
    const int core_count = 6;
    const parameters core_params[core_count] = {iDmolar, iT, iSmolar, iHmolar, iP, iUmolar};
    const int aux_count = 3;
    const parameters aux_params[aux_count] = {ispeed_sound, iviscosity, iconductivity};

    coeffs.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));

    // First pass: cell validity by 4 finite corners across core params only.
    for (std::size_t i = 0; i < table.Nx - 1; ++i) {
        for (std::size_t j = 0; j < table.Ny - 1; ++j) {
            bool ok = true;
            for (int k = 0; k < core_count && ok; ++k) {
                parameters param = core_params[k];
                if (param == table.xkey || param == table.ykey) continue;
                std::vector<std::vector<double>>* fp = sbtl_get_field(table, param);
                if (!ValidNumber((*fp)[i][j]) || !ValidNumber((*fp)[i + 1][j]) || !ValidNumber((*fp)[i][j + 1])
                    || !ValidNumber((*fp)[i + 1][j + 1])) {
                    ok = false;
                }
            }
            if (ok) {
                coeffs[i][j].set_valid();
                coeffs[i][j].dx_dxhat = table.xvec[i + 1] - table.xvec[i];
                coeffs[i][j].dy_dyhat = std::log(table.yvec[j + 1]) - std::log(table.yvec[j]);
            } else {
                coeffs[i][j].set_invalid();
            }
        }
    }

    // Second pass: per-property B-spline + per-cell polynomial expansion.
    // Self-check: after computing alpha, evaluate the polynomial at each
    // cell corner (t_x ∈ {0, 1}, t_y ∈ {0, 1}) and verify the value matches
    // the corresponding stored data value.  Mismatches indicate a bug in the
    // B-spline construction or the polynomial expansion — print and abort
    // the property's per-cell fill so the test loop can find them.
    int self_check_failures = 0;
    const int total_count = core_count + aux_count;
    for (int k = 0; k < total_count; ++k) {
        const bool is_aux = (k >= core_count);
        parameters param = is_aux ? aux_params[k - core_count] : core_params[k];
        if (param == table.xkey || param == table.ykey) continue;
        std::vector<std::vector<double>>* fp = sbtl_get_field(table, param);

        // For aux params we replace any non-finite corner with 0.0 before
        // building the global B-spline so the spline operator stays well
        // defined; per-cell we then skip writing the polynomial into cells
        // that touch a non-finite corner.  Core params are guaranteed
        // finite for cells flagged valid by the first pass.
        std::vector<std::vector<double>> aux_scratch;
        if (is_aux) {
            aux_scratch = *fp;  // copy
            for (auto& row : aux_scratch) {
                for (double& v : row) {
                    if (!ValidNumber(v)) v = 0.0;
                }
            }
        }
        const auto& source = is_aux ? aux_scratch : *fp;

        const auto C = bspline_coeffs_2d(source);

        for (std::size_t i = 0; i < table.Nx - 1; ++i) {
            for (std::size_t j = 0; j < table.Ny - 1; ++j) {
                if (!coeffs[i][j].valid()) continue;

                if (is_aux) {
                    // Skip cells touching a non-finite corner — they get no
                    // aux polynomial and calc_speed_sound etc. will fall back
                    // to HEOS (or throw) at lookup time.
                    if (!ValidNumber((*fp)[i][j]) || !ValidNumber((*fp)[i + 1][j]) || !ValidNumber((*fp)[i][j + 1])
                        || !ValidNumber((*fp)[i + 1][j + 1])) {
                        continue;
                    }
                }

                std::vector<double> alpha = bspline_polynomial_coeffs(C, i, j);

                // Self-check: polynomial value at (0,0) corner must equal f[i][j].
                const double v00 = eval_alpha(alpha, 0.0, 0.0);
                const double v10 = eval_alpha(alpha, 1.0, 0.0);
                const double v01 = eval_alpha(alpha, 0.0, 1.0);
                const double v11 = eval_alpha(alpha, 1.0, 1.0);
                const double f00 = source[i][j];
                const double f10 = source[i + 1][j];
                const double f01 = source[i][j + 1];
                const double f11 = source[i + 1][j + 1];
                const double scale = std::max({std::abs(f00), std::abs(f10), std::abs(f01), std::abs(f11), 1e-30});
                const double tol = 1e-8 * scale;
                if (self_check_failures < 5
                    && (std::abs(v00 - f00) > tol || std::abs(v10 - f10) > tol || std::abs(v01 - f01) > tol || std::abs(v11 - f11) > tol)) {
                    std::cout << "SBTL bspline self-check fail param=" << static_cast<int>(param) << " cell=(" << i << "," << j << ")  f=(" << f00
                              << "," << f10 << "," << f01 << "," << f11 << ")  v=(" << v00 << "," << v10 << "," << v01 << "," << v11 << ")\n";
                    ++self_check_failures;
                }
                coeffs[i][j].set(param, alpha);
            }
        }
    }

    // Third pass: alternates for invalid cells.
    for (std::size_t i = 0; i < table.Nx - 1; ++i) {
        for (std::size_t j = 0; j < table.Ny - 1; ++j) {
            if (coeffs[i][j].valid()) continue;
            for (int r = 1; r < static_cast<int>(std::max(table.Nx, table.Ny)); ++r) {
                bool found = false;
                for (int di = -r; di <= r && !found; ++di) {
                    for (int dj = -r; dj <= r && !found; ++dj) {
                        const int ii = static_cast<int>(i) + di, jj = static_cast<int>(j) + dj;
                        if (ii < 0 || jj < 0 || ii >= static_cast<int>(table.Nx - 1) || jj >= static_cast<int>(table.Ny - 1)) continue;
                        if (!coeffs[ii][jj].valid()) continue;
                        coeffs[i][j].set_alternate(static_cast<std::size_t>(ii), static_cast<std::size_t>(jj));
                        found = true;
                    }
                }
                if (found) break;
            }
        }
    }
}

double NormalizedPHTable::h_from_xnorm(double xnorm, double P, double h_sat) const {
    auto cheb_or_isobar = [&](const Cheb1D& cb, const std::vector<double>& fallback) -> double {
        if (cb.valid()) return cb.eval(P);
        if (fallback.empty()) return std::numeric_limits<double>::quiet_NaN();
        const std::size_t Nrows = yvec.size();
        std::size_t j = 0;
        bisect_vector(yvec, P, j);
        if (j >= Nrows - 1) j = Nrows - 2;
        const double frac = (std::log(P) - std::log(yvec[j])) / (std::log(yvec[j + 1]) - std::log(yvec[j]));
        return fallback[j] + frac * (fallback[j + 1] - fallback[j]);
    };
    double h_lo, h_hi;
    switch (region_) {
        case LIQUID:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = std::isfinite(h_sat) ? h_sat : cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case VAPOR:
            h_lo = std::isfinite(h_sat) ? h_sat : cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case SUPER:
        default:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
    }
    return h_lo + xnorm * (h_hi - h_lo);
}

// ---------------------------------------------------------------------------
// H-superancillary fast path: returns true and writes h_sat,Q(p) to `out`.
// Triggers HelmholtzEOSMixtureBackend::ensure_caloric_superancillaries() on
// first call so the H/S expansions are built lazily via add_variable() — this
// requires the SuperAncillary to be instantiated with Eigen::ArrayXd (see
// CoolPropFluid.h's SuperAncillary_t alias).  Returns false for mixtures or
// fluids without an underlying SUPERANCILLARY block.
// ---------------------------------------------------------------------------

bool SBTLBackend::try_h_superanc(double p, short Q, double& out) const {
    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (!helm) return false;
    std::shared_ptr<EquationOfState::SuperAncillary_t> super_anc;
    try {
        // ensure_caloric_superancillaries() is a no-op if H/S are already
        // built; on first call it pays the ~10–50 ms one-time build cost.
        helm->ensure_caloric_superancillaries();
        super_anc = helm->get_components()[0].EOS().get_superanc();
    } catch (...) {
        return false;
    }
    if (!super_anc) return false;
    try {
        const double T_sat = super_anc->get_invlnp().value().eval(std::log(p));
        out = super_anc->get_approx1d('H', Q).eval(T_sat);
        return true;
    } catch (...) {
        return false;
    }
}

void SBTLBackend::saturation_hmolar_LV(double p, double& h_L, double& h_V) {
    // 1-deep cache: subsequent queries at the same p reuse the result.
    if (p == _hsat_LV_cache_p) {
        h_L = _hsat_LV_cache_hL;
        h_V = _hsat_LV_cache_hV;
        return;
    }
    // Try the H-superancillary fast path with both Q=0 and Q=1, sharing
    // the T_sat(p) computation between the two phase queries.
    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (helm) {
        try {
            helm->ensure_caloric_superancillaries();
            auto super_anc = helm->get_components()[0].EOS().get_superanc();
            if (super_anc) {
                const double T_sat = super_anc->get_invlnp().value().eval(std::log(p));
                h_L = super_anc->get_approx1d('H', 0).eval(T_sat);
                h_V = super_anc->get_approx1d('H', 1).eval(T_sat);
                _hsat_LV_cache_p = p;
                _hsat_LV_cache_hL = h_L;
                _hsat_LV_cache_hV = h_V;
                return;
            }
        } catch (...) {
            // Fall through to per-Q HEOS update calls.
        }
    }
    // Fallback: separate PQ_INPUTS calls for fluids without an
    // H-superancillary block.
    h_L = saturation_hmolar_liquid(p);
    h_V = saturation_hmolar_vapor(p);
    _hsat_LV_cache_p = p;
    _hsat_LV_cache_hL = h_L;
    _hsat_LV_cache_hV = h_V;
}

// ---------------------------------------------------------------------------
// Cell finding — trivial bisect on the normalized table's xvec/yvec.
//
// All cells on a NormalizedPHTable / NormalizedPTTable are valid by
// construction (build_bspline_coeffs marks a cell valid iff its four
// data corners are finite, and the table-fill skips entire isobars
// where the bounds aren't well-defined).  The only failure mode is OOB
// (caller passed an xnorm outside [0, 1] or a P outside [yvec.front,
// yvec.back]) — guarded upstream in update().
// ---------------------------------------------------------------------------
void SBTLBackend::find_native_nearest_good_indices(SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs, double x,
                                                   double y, std::size_t& i, std::size_t& j) {
    bisect_vector(table.xvec, x, i);
    bisect_vector(table.yvec, y, j);
    if (!coeffs[i][j].valid()) {
        throw ValueError(format("SBTLBackend: cell (%zu, %zu) is invalid at xnorm=%g, P=%g", i, j, x, y));
    }
}

void SBTLBackend::find_nearest_neighbor(SinglePhaseGriddedTableData& /*table*/, const std::vector<std::vector<CellCoeffs>>& /*coeffs*/,
                                        const parameters /*variable1*/, const double /*value1*/, const parameters /*otherkey*/,
                                        const double /*otherval*/, std::size_t& /*i*/, std::size_t& /*j*/) {
    throw ValueError("SBTL backend (pure-fluid PT/PH only) does not support find_nearest_neighbor "
                     "(used for non-PT/PH input pairs).  Use PT_INPUTS or HmassP_INPUTS, or use HEOS for other inputs.");
}

// ---------------------------------------------------------------------------
// Transport properties: bilinear interpolation
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase_transport(SinglePhaseGriddedTableData& /*table*/, parameters /*output*/, double /*x*/, double /*y*/,
                                                    std::size_t /*i*/, std::size_t /*j*/) {
    throw ValueError("SBTL: tabular transport not supported; use HEOS direct via viscosity()/conductivity() accessors");
}

// ---------------------------------------------------------------------------
// Polynomial evaluation: bi-cubic (degree-3, 16 coeffs) per-cell polynomial
// at (xnorm, log P) in cell (i, j) of an active NormalizedPHTable /
// NormalizedPTTable.  Built against xnorm ∈ [0, 1] × uniform-log-P, so
// xi is linear in xnorm and eta is linear in log P.
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase(const SinglePhaseGriddedTableData& table, const std::vector<std::vector<CellCoeffs>>& coeffs,
                                          const parameters output, const double x, const double y, const std::size_t i, const std::size_t j) {
    const CellCoeffs& cell = coeffs[i][j];
    const std::vector<double>& alpha = cell.get(output);

    const double xi = (x - table.xvec[i]) / (table.xvec[i + 1] - table.xvec[i]);
    const double eta = (std::log(y) - std::log(table.yvec[j])) / (std::log(table.yvec[j + 1]) - std::log(table.yvec[j]));

    const double xi2 = xi * xi, xi3 = xi2 * xi;
    const double eta2 = eta * eta, eta3 = eta2 * eta;
    const double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2 + alpha[3] * eta3;
    const double B1 = alpha[4] + alpha[5] * eta + alpha[6] * eta2 + alpha[7] * eta3;
    const double B2 = alpha[8] + alpha[9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
    const double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
    const double val = B0 + B1 * xi + B2 * xi2 + B3 * xi3;

    switch (output) {
        case iT:
            _T = val;
            break;
        case iP:
            _p = val;
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
        case ispeed_sound:
            _speed_sound = val;
            break;
        case iviscosity:
            _viscosity = val;
            break;
        case iconductivity:
            _conductivity = val;
            break;
        default:
            throw ValueError("Invalid output variable in SBTLBackend::evaluate_single_phase");
    }
    return val;
}

// ---------------------------------------------------------------------------
// Derivative and inversion paths — not supported.  Derivative-driven
// properties (cp, cv, speed of sound, viscosity, conductivity, first_partial_deriv)
// route through HEOS direct; inversion is only used by the base-class machinery
// for non-PT/PH input pairs, which we reject upstream in update().
// ---------------------------------------------------------------------------

double SBTLBackend::evaluate_single_phase_derivative(SinglePhaseGriddedTableData& /*table*/, std::vector<std::vector<CellCoeffs>>& /*coeffs*/,
                                                     parameters /*output*/, double /*x*/, double /*y*/, std::size_t /*i*/, std::size_t /*j*/,
                                                     std::size_t /*Nx*/, std::size_t /*Ny*/) {
    throw ValueError("SBTL: polynomial derivatives not supported; use HEOS direct via calc_first_partial_deriv or the cp/cv/c_sound accessors");
}

void SBTLBackend::invert_single_phase_x(const SinglePhaseGriddedTableData& /*table*/, const std::vector<std::vector<CellCoeffs>>& /*coeffs*/,
                                        parameters /*other_key*/, double /*other*/, double /*y*/, std::size_t /*i*/, std::size_t /*j*/) {
    throw ValueError("SBTL backend (pure-fluid PT/PH only) does not support input-pair inversion — "
                     "use PT_INPUTS or HmassP_INPUTS directly, or HEOS for other inputs");
}

void SBTLBackend::invert_single_phase_y(const SinglePhaseGriddedTableData& /*table*/, const std::vector<std::vector<CellCoeffs>>& /*coeffs*/,
                                        parameters /*other_key*/, double /*other*/, double /*x*/, std::size_t /*i*/, std::size_t /*j*/) {
    throw ValueError("SBTL backend (pure-fluid PT/PH only) does not support input-pair inversion — "
                     "use PT_INPUTS or HmassP_INPUTS directly, or HEOS for other inputs");
}

// ---------------------------------------------------------------------------
// update() override: routes HmolarP/HmassP and PT through the
// coordinate-aligned normalized tables.  Anything else (other input
// pairs, imposed two-phase) throws cleanly — SBTL is pure-fluid PT/PH only.
//
// State-machine invariant for the property accessors below: at the top
// of every update() call we (a) invalidate the AbstractState cached
// scalars via clear(), (b) zero ALL active-table pointers and the
// critbox flag.  Each fast-path branch that succeeds then sets EXACTLY
// ONE of {_critbox_active, _normph_active_table, _normpt_active_table}.
// Without the unconditional reset, a PH update followed by a PT update
// would leave _normph_active_table set and the property accessors would
// silently route through the wrong (stale-h, stale-xnorm) table.
// ---------------------------------------------------------------------------
void SBTLBackend::update(CoolProp::input_pairs input_pair, double val1, double val2) {
    // Invariant for downstream property accessors (calc_hmolar, calc_smolar,
    // calc_cpmolar, calc_viscosity, ...): they read from the active
    // normalized table (via _normph_active_table / _normpt_active_table) or
    // the critbox-cached values.  All caches must be invalidated at the
    // entry to update() so a successful fast path can leave its own state
    // behind without inheriting stale state from a prior call.  The
    // AS->clear() invalidates the AbstractState's memoized scalar cache
    // (hmolar/smolar/umolar/cpmolar/viscosity/...) — without it, the
    // second update() in
    //   AS->update(PH, h1, p1); _ = AS->smolar();
    //   AS->update(PT, p2, T2); s = AS->smolar();     // returns cached s1
    // would silently return the previous value.
    this->clear();
    _critbox_active = false;
    _normph_active_table = nullptr;
    _normph_active_coeffs = nullptr;
    _normpt_active_table = nullptr;
    _normpt_active_coeffs = nullptr;

    // imposed_phase_index handling: the normalized tables come in L/V/SUPER
    // flavours that map directly onto the imposed-phase enum, so a user's
    // specify_phase() can be honoured WITHOUT falling through to the legacy
    // tables.  This matters in particular on the saturation curve: a
    // PT_INPUTS query at exactly T = T_sat(P) is formally ambiguous between
    // saturated-liquid (rho_L) and saturated-vapor (rho_V) — the user
    // disambiguates via specify_phase(iphase_liquid|iphase_gas) before the
    // update, and we route to the corresponding table.  The mapping:
    //   iphase_not_imposed                  -> auto: route by T vs T_sat (or h vs h_sat)
    //   iphase_liquid                       -> LIQUID table (xnorm=1 boundary is sat-L)
    //   iphase_gas                          -> VAPOR  table (xnorm=0 boundary is sat-V)
    //   iphase_supercritical                -> SUPER  table
    //   iphase_supercritical_liquid|gas     -> SUPER  table
    //   iphase_twophase                     -> rejected (PT in dome is ambiguous)
    //   iphase_critical_point               -> SUPER table (degenerate; routes through critbox)
    // Pure-fluid only: mixtures (real or pseudo-pure) gate the fast path off
    // entirely, so they still fall through to the base-class imposed-phase
    // logic against the legacy tables.

    // Helper: derive a region tag from imposed_phase_index for use inside the
    // PH and PT fast paths.  Returns -1 (auto), 0 (LIQUID), 1 (VAPOR), 2 (SUPER),
    // or throws for iphase_twophase.  iphase_critical_point routes through the
    // SUPER table; queries near (T_c, p_c) are caught upstream by the critbox.
    enum class ForcedRegion
    {
        Auto = -1,
        Liquid = 0,
        Vapor = 1,
        Super = 2,
        TwoPhase = 3
    };
    const ForcedRegion forced = [&]() -> ForcedRegion {
        switch (imposed_phase_index) {
            case iphase_not_imposed:
                return ForcedRegion::Auto;
            case iphase_liquid:
                return ForcedRegion::Liquid;
            case iphase_gas:
                return ForcedRegion::Vapor;
            case iphase_supercritical:
            case iphase_supercritical_liquid:
            case iphase_supercritical_gas:
            case iphase_critical_point:
                return ForcedRegion::Super;
            case iphase_twophase:
                return ForcedRegion::TwoPhase;
            default:
                return ForcedRegion::Auto;
        }
    }();

    // Two-phase imposed: PT_INPUTS in the dome is ambiguous (rho_L vs rho_V),
    // and HmassP_INPUTS at exactly h_sat is also ambiguous.  Reject cleanly
    // so the caller switches to PQ_INPUTS or specify_phase(iphase_liquid|gas).
    if (forced == ForcedRegion::TwoPhase) {
        throw ValueError("SBTL: imposed two-phase is not supported (PT_INPUTS / HmassP_INPUTS are ambiguous in the dome); "
                         "use PQ_INPUTS or specify_phase(iphase_liquid|iphase_gas) to disambiguate");
    }

    // Coordinate-aligned PH fast path.
    if (normph_tables_built && (input_pair == HmolarP_INPUTS || input_pair == HmassP_INPUTS)) {
        double h_molar = NAN, p = NAN;
        if (input_pair == HmolarP_INPUTS) {
            h_molar = val1;
            p = val2;
        } else {
            // HmassP -> molar: multiply mass-h by molar mass (h_molar = h_mass * M)
            h_molar = val1 * static_cast<double>(this->AS->molar_mass());
            p = val2;
        }

        const double p_crit = static_cast<double>(this->AS->p_critical());

        // Critical-region HEOS-fallback box.  Within a small (P, h)
        // neighbourhood of the critical point, no smooth interpolation
        // (B-spline, Lagrange, or any polynomial) can faithfully represent
        // the cusp in ρ(h, P) that arises because cp diverges at (Tc, pc)
        // and the saturation curve terminates with a vertical tangent there.
        // Mirroring Kunick's IAPWS-IF97 design, the iterative region is a
        // SMALL 2D box defined on both axes — much tighter than a 1D band
        // on P alone, since queries far from h_crit but at p ≈ pc (e.g.
        // hot supercritical gas at p = pc + 5 bar) are well-behaved and
        // shouldn't pay the ~100 us EOS-direct cost.  Box dimensions
        // chosen as ~1.5 % of (h_crit, p_crit) — comfortably wider than
        // the polynomial-interpolation breakdown zone but tight enough
        // to leave most queries on the fast SBTL path.
        // 2D HEOS-fallback box around the critical point.  Now that the
        // SUPER normph table covers down to T_melt(P) at all P > p_crit,
        // both sides of the critical pressure use the SAME tight box —
        // only queries with h close to h_crit (where the polynomial
        // backbones can't reproduce the cusp) actually need EOS-direct.
        const double dp_below = 0.10e6;    // ±1 bar below pcrit
        const double dp_above = 0.10e6;    // ±1 bar above pcrit
        const double dh_crit_frac = 0.30;  // ±30 % of |h_crit|
        // Lazy-compute h_crit at the critical point on the first query
        // that lands near pcrit.  std::optional<double> keeps the unset
        // state explicit instead of relying on a NaN/Inf sentinel.
        if (!_hmolar_crit_cached.has_value() && std::abs(p - p_crit) < 2.0 * dp_above) {
            try {
                this->AS->update(DmolarT_INPUTS, this->AS->rhomolar_critical(), this->AS->T_critical());
                _hmolar_crit_cached = static_cast<double>(this->AS->hmolar());
            } catch (...) {
                _hmolar_crit_cached = std::numeric_limits<double>::infinity();
            }
        }
        const double h_crit = _hmolar_crit_cached.value_or(std::numeric_limits<double>::infinity());
        const bool in_box = std::isfinite(h_crit) && std::abs(p - p_crit) < ((p < p_crit) ? dp_below : dp_above)
                            && std::abs(h_molar - h_crit) < dh_crit_frac * std::abs(h_crit);
        if (in_box) {
            try {
                this->AS->update(HmolarP_INPUTS, h_molar, p);
                _hmolar = h_molar;
                _p = p;
                _T = static_cast<double>(this->AS->T());
                _rhomolar = static_cast<double>(this->AS->rhomolar());
                _Q = -1000;
                // Honor imposed_phase if the user specified one — same
                // contract as the normph fast path below.  iphase_critical_point
                // is excluded so a user with imposed_critical_point at a state
                // far from the critical point doesn't get mislabeled.
                if (imposed_phase_index != iphase_not_imposed && imposed_phase_index != iphase_critical_point) {
                    _phase = imposed_phase_index;
                } else {
                    _phase = this->AS->phase();
                }
                _critbox_active = true;
                _critbox_smolar = static_cast<double>(this->AS->smolar());
                _critbox_umolar = static_cast<double>(this->AS->umolar());
                _critbox_cpmolar = static_cast<double>(this->AS->cpmolar());
                _critbox_cvmolar = static_cast<double>(this->AS->cvmolar());
                _critbox_speed_sound = static_cast<double>(this->AS->speed_sound());
                try {
                    _critbox_viscosity = static_cast<double>(this->AS->viscosity());
                } catch (...) {
                    _critbox_viscosity = _HUGE;
                }
                try {
                    _critbox_conductivity = static_cast<double>(this->AS->conductivity());
                } catch (...) {
                    _critbox_conductivity = _HUGE;
                }
                using_single_phase_table = true;
                selected_table = SELECTED_PH_TABLE;
                return;
            } catch (...) {
                // EOS rejected the state (out of range, etc.) — fall through.
                _critbox_active = false;
            }
        }

        NormalizedPHTable* tbl = nullptr;
        const std::vector<std::vector<CellCoeffs>>* coeffs = nullptr;
        // h_sat_active is the sat-side boundary value to pass down into
        // xnorm_from_h: routing + xnorm use the SAME h_sat (machine-
        // precision via H-superancillary), no boundary disagreement.
        // NaN means no sat boundary applies (SUPER region).
        double h_sat_active = std::numeric_limits<double>::quiet_NaN();
        // Region selection: imposed phase (when set) wins over the automatic
        // sat-curve comparison.  This is how a user disambiguates a query
        // at exactly h = h_sat,L (or h_sat,V) — specify_phase(iphase_liquid)
        // routes to LIQUID, specify_phase(iphase_gas) to VAPOR.
        if (forced == ForcedRegion::Super || (forced == ForcedRegion::Auto && p >= p_crit)) {
            tbl = &_normph_super;
            coeffs = &_coeffs_normph_super;
        } else if (forced == ForcedRegion::Liquid) {
            // Sub- or supercritical: user wants liquid.  Below p_crit use the
            // LIQUID table and the H-superancillary sat-L boundary; above
            // p_crit there's no sat curve to anchor to, so fall back to SUPER.
            if (p < p_crit) {
                tbl = &_normph_liquid;
                coeffs = &_coeffs_normph_liquid;
                try {
                    double h_sat_L = NAN, h_sat_V = NAN;
                    saturation_hmolar_LV(p, h_sat_L, h_sat_V);
                    h_sat_active = h_sat_L;
                } catch (...) {
                }
            } else {
                tbl = &_normph_super;
                coeffs = &_coeffs_normph_super;
            }
        } else if (forced == ForcedRegion::Vapor) {
            if (p < p_crit) {
                tbl = &_normph_vapor;
                coeffs = &_coeffs_normph_vapor;
                try {
                    double h_sat_L = NAN, h_sat_V = NAN;
                    saturation_hmolar_LV(p, h_sat_L, h_sat_V);
                    h_sat_active = h_sat_V;
                } catch (...) {
                }
            } else {
                tbl = &_normph_super;
                coeffs = &_coeffs_normph_super;
            }
        } else {
            // forced == ForcedRegion::Auto and p < p_crit: use h vs sat to pick.
            double h_sat_L = std::numeric_limits<double>::quiet_NaN();
            double h_sat_V = std::numeric_limits<double>::quiet_NaN();
            try {
                saturation_hmolar_LV(p, h_sat_L, h_sat_V);
            } catch (...) {
            }
            if (std::isfinite(h_sat_L) && std::isfinite(h_sat_V)) {
                if (h_molar <= h_sat_L) {
                    tbl = &_normph_liquid;
                    coeffs = &_coeffs_normph_liquid;
                    h_sat_active = h_sat_L;  // LIQUID's h_hi
                } else if (h_molar >= h_sat_V) {
                    tbl = &_normph_vapor;
                    coeffs = &_coeffs_normph_vapor;
                    h_sat_active = h_sat_V;  // VAPOR's h_lo
                }
                // else: two-phase; fall through to base class which has
                // saturation-blend logic.
            }
        }

        if (tbl) {
            try {
                const double xnorm = tbl->xnorm_from_h(h_molar, p, h_sat_active);
                // build_bspline_coeffs marks boundary cells valid (i=0 and
                // i=Nx-2) as long as their 4 data corners are finite, and the
                // resulting cubic B-spline polynomial interpolates exactly
                // at those corners.  Don't reject boundary queries — the
                // boundary cell handles its own range cleanly with proper
                // not-a-knot ghost extension.
                if (std::isfinite(xnorm) && xnorm >= 0.0 && xnorm <= 1.0) {
                    std::size_t i = 0, j = 0;
                    find_native_nearest_good_indices(*tbl, *coeffs, xnorm, p, i, j);
                    const CellCoeffs& cell = (*coeffs)[i][j];
                    if (cell.valid()) {
                        // Set state: evaluate the table directly (xnorm × P) to
                        // get T, rho, etc.  These hide the base class's cache
                        // because we're not going through it for this path.
                        _hmolar = h_molar;
                        _p = p;
                        const double T_val = evaluate_single_phase(*tbl, *coeffs, iT, xnorm, p, i, j);
                        const double rho_val = evaluate_single_phase(*tbl, *coeffs, iDmolar, xnorm, p, i, j);
                        // Sanity check: cubic B-spline through hole-filled
                        // data near table edges can extrapolate to
                        // non-physical values (e.g. negative density at
                        // cold compressed liquid where the y-row was
                        // skipped).  Reject and fall through to base class.
                        if (!std::isfinite(rho_val) || rho_val <= 0.0 || !std::isfinite(T_val) || T_val <= 0.0) {
                            throw ValueError("SBTL normph polynomial gave non-physical T or rho");
                        }
                        _T = T_val;
                        _rhomolar = rho_val;
                        // Q-tag at saturation boundaries (LIQUID xnorm=1 ↔ sat-L,
                        // VAPOR xnorm=0 ↔ sat-V).  Matches the PT path so both
                        // input pairs return consistent Q at the dome boundary.
                        // Tolerance is loose enough to absorb 1e-10 round-off
                        // in xnorm from the (h - h_sat)/(h_hi - h_sat) algebra.
                        if (tbl == &_normph_liquid && std::abs(xnorm - 1.0) < 1e-10) {
                            _Q = 0;
                        } else if (tbl == &_normph_vapor && std::abs(xnorm) < 1e-10) {
                            _Q = 1;
                        } else {
                            _Q = -1000;
                        }
                        // Phase tag: imposed phase wins if set (caller intent
                        // is documented); otherwise derive from active table.
                        if (imposed_phase_index != iphase_not_imposed && imposed_phase_index != iphase_critical_point) {
                            _phase = imposed_phase_index;
                        } else if (tbl == &_normph_super) {
                            _phase = (T_val > static_cast<double>(this->AS->T_critical())) ? iphase_supercritical_gas : iphase_supercritical_liquid;
                        } else if (tbl == &_normph_liquid) {
                            _phase = iphase_liquid;
                        } else {
                            _phase = iphase_gas;
                        }
                        using_single_phase_table = true;
                        selected_table = SELECTED_PH_TABLE;
                        // Stash the indices in the existing cache slots so
                        // SBTL property accessors (which evaluate via
                        // evaluate_single_phase_phmolar) get routed through the
                        // normph path via the active-table override below.
                        cached_single_phase_i = i;
                        cached_single_phase_j = j;
                        _normph_active_table = tbl;
                        _normph_active_coeffs = coeffs;
                        _normph_xnorm = xnorm;
                        return;
                    }
                }
            } catch (...) {
                // Fall through to base class on any failure.
            }
        }
    }

    // Coordinate-aligned PT fast path (pure fluids only; built in ctor).
    // Mirror of the PH normph block above — selects LIQUID / VAPOR / SUPER
    // based on (T, P) vs sat curve, then evaluates via xnorm = (T - T_lo) / (T_hi - T_lo)
    // on the appropriate normpt table.  The cell boundary at xnorm=1 (LIQUID)
    // is exactly the saturation curve, so cells never straddle the dome.
    if (normpt_tables_built && input_pair == PT_INPUTS) {
        const double p = val1;
        const double T_query = val2;
        const double p_crit = static_cast<double>(this->AS->p_critical());
        const double T_crit = static_cast<double>(this->AS->T_critical());

        // Critical-region HEOS-fallback box (mirror of the PH path's box,
        // widened for PT).  Two effects to capture:
        //   (a) Cusp at the critical point itself — narrow, ±0.1 MPa,
        //       ±5% T_crit suffices on the subcritical side.
        //   (b) Supercritical shoulder where ρ(T, P) drops steeply above
        //       p_crit and just above T_crit — extends out to ~p_crit*1.3
        //       and ~T_crit*1.2.  Polynomial backbones can't reproduce
        //       this shoulder at the resolution of a 200-cell SUPER grid;
        //       routing through HEOS direct here recovers BICUBIC-grade
        //       accuracy at ~10 µs per query.
        // Critical-region HEOS-fallback box dimensions.  The cusp at
        // (T_c, p_c) plus the supercritical-shoulder where ρ(T, p) drops
        // steeply between T_c and T_c*1.25 at p just above p_c — and the
        // mirror of that shoulder on the subcritical side just below the
        // critical point — are all routed to HEOS direct (~10 µs / call)
        // instead of through the polynomial backbone.  Numbers chosen
        // empirically to bring the random-PT max relative error from
        // ~3 % down to ~1e-3 across a wide-window CO2 / Water sweep.
        const double dp_below_frac = 0.25;  // -25 % p_crit below
        const double dp_above_frac = 0.75;  // +75 % p_crit above
        const double dT_below_frac = 0.10;  // -10 % T_crit below
        const double dT_above_frac = 0.30;  // +30 % T_crit above
        const bool in_crit_box = (p < p_crit ? (p_crit - p) < dp_below_frac * p_crit : (p - p_crit) < dp_above_frac * p_crit)
                                 && (T_query < T_crit ? (T_crit - T_query) < dT_below_frac * T_crit : (T_query - T_crit) < dT_above_frac * T_crit);
        if (in_crit_box) {
            try {
                this->AS->update(PT_INPUTS, p, T_query);
                _T = T_query;
                _p = p;
                _rhomolar = static_cast<double>(this->AS->rhomolar());
                _hmolar = static_cast<double>(this->AS->hmolar());
                _Q = -1000;
                if (imposed_phase_index != iphase_not_imposed && imposed_phase_index != iphase_critical_point) {
                    _phase = imposed_phase_index;
                } else {
                    _phase = this->AS->phase();
                }
                _critbox_active = true;
                _critbox_smolar = static_cast<double>(this->AS->smolar());
                _critbox_umolar = static_cast<double>(this->AS->umolar());
                _critbox_cpmolar = static_cast<double>(this->AS->cpmolar());
                _critbox_cvmolar = static_cast<double>(this->AS->cvmolar());
                _critbox_speed_sound = static_cast<double>(this->AS->speed_sound());
                try {
                    _critbox_viscosity = static_cast<double>(this->AS->viscosity());
                } catch (...) {
                    _critbox_viscosity = _HUGE;
                }
                try {
                    _critbox_conductivity = static_cast<double>(this->AS->conductivity());
                } catch (...) {
                    _critbox_conductivity = _HUGE;
                }
                using_single_phase_table = true;
                selected_table = SELECTED_PT_TABLE;
                return;
            } catch (...) {
                _critbox_active = false;
                // Fall through to normpt path.
            }
        }

        NormalizedPTTable* tbl = nullptr;
        const std::vector<std::vector<CellCoeffs>>* coeffs = nullptr;
        double T_sat_active = std::numeric_limits<double>::quiet_NaN();

        // OOB pressure guard runs BEFORE saturation_T_LV: the PQ_INPUTS
        // inside saturation_T_LV silently extrapolates the sat curve
        // along the metastable branch when p < p_triple, mutating AS into
        // a non-physical state.  Detect the OOB range upfront against the
        // appropriate table's yvec.front/back so we throw cleanly without
        // touching the AS state.
        const auto pt_range_for = [&](NormalizedPTTable* t) -> std::pair<double, double> { return {t->yvec.front(), t->yvec.back()}; };
        // Region selection: imposed phase wins over the automatic
        // sat-curve comparison.  This is essential at PT queries on the
        // saturation curve: T = T_sat,L(P) is ambiguous between sat-L and
        // sat-V, and the user disambiguates via
        // specify_phase(iphase_liquid|iphase_gas) before the update.
        //   - forced=Liquid / Vapor : route to the corresponding subcritical
        //     table at p<p_crit; fall back to SUPER at p≥p_crit (where no
        //     sat curve exists to anchor to).
        //   - forced=Super          : route to SUPER table regardless of p
        //     (caller asked for the supercritical phase, e.g.
        //     iphase_supercritical_liquid for p>p_c, T<T_c).
        //   - forced=Auto           : pick LIQUID/VAPOR/SUPER from p vs p_c
        //     and T vs T_sat.
        if (forced == ForcedRegion::Super) {
            // User-imposed supercritical, regardless of p vs p_crit.
            auto [pmin, pmax] = pt_range_for(&_normpt_super);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the SUPER normpt "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            tbl = &_normpt_super;
            coeffs = &_coeffs_normpt_super;
        } else if (forced == ForcedRegion::Liquid && p < p_crit) {
            // User-imposed liquid (sub-critical p).
            auto [pmin, pmax] = pt_range_for(&_normpt_liquid);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the LIQUID normpt "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            // Read T_sat,L for the xnorm=1 boundary; xnorm_from_T uses it.
            // Sat lookup may throw if p is outside the sat range; let
            // xnorm_from_T fall back to the Cheb-fitted boundary.
            try {
                double T_sat_L = NAN, T_sat_V = NAN;
                saturation_T_LV(p, T_sat_L, T_sat_V);
                T_sat_active = T_sat_L;
            } catch (...) {
            }
            tbl = &_normpt_liquid;
            coeffs = &_coeffs_normpt_liquid;
        } else if (forced == ForcedRegion::Vapor && p < p_crit) {
            // User-imposed vapor (sub-critical p).
            auto [pmin, pmax] = pt_range_for(&_normpt_vapor);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the VAPOR normpt "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            try {
                double T_sat_L = NAN, T_sat_V = NAN;
                saturation_T_LV(p, T_sat_L, T_sat_V);
                T_sat_active = T_sat_V;
            } catch (...) {
            }
            tbl = &_normpt_vapor;
            coeffs = &_coeffs_normpt_vapor;
        } else if ((forced == ForcedRegion::Liquid || forced == ForcedRegion::Vapor) && p >= p_crit) {
            // User imposed liquid/gas but p is supercritical.  Route through
            // SUPER (which covers both sides of the T_c isotherm at p≥p_c).
            auto [pmin, pmax] = pt_range_for(&_normpt_super);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the SUPER normpt "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            tbl = &_normpt_super;
            coeffs = &_coeffs_normpt_super;
        } else if (p >= p_crit) {
            auto [pmin, pmax] = pt_range_for(&_normpt_super);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the SUPER normpt "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            tbl = &_normpt_super;
            coeffs = &_coeffs_normpt_super;
        } else if (p > 0.0) {
            // forced=Auto, subcritical: pick LIQUID/VAPOR from T vs T_sat.
            auto [pmin, pmax] = pt_range_for(&_normpt_liquid);
            if (p < pmin || p > pmax) {
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa is outside the subcritical "
                                        "normpt table range [%g, %g] Pa — use HEOS directly.",
                                        p, pmin, pmax));
            }
            double T_sat_L = std::numeric_limits<double>::quiet_NaN();
            double T_sat_V = std::numeric_limits<double>::quiet_NaN();
            try {
                saturation_T_LV(p, T_sat_L, T_sat_V);
            } catch (...) {
            }
            if (std::isfinite(T_sat_L) && std::isfinite(T_sat_V)) {
                // T = T_sat (exact saturation): convention is liquid side
                // (sat-L properties at the LIQUID xnorm=1 boundary).
                // Caller wanting vapor-side must specify_phase(iphase_gas).
                if (T_query <= T_sat_L) {
                    tbl = &_normpt_liquid;
                    coeffs = &_coeffs_normpt_liquid;
                    T_sat_active = T_sat_L;
                } else if (T_query >= T_sat_V) {
                    tbl = &_normpt_vapor;
                    coeffs = &_coeffs_normpt_vapor;
                    T_sat_active = T_sat_V;
                } else {
                    // Strictly inside the dome — invalid PT_INPUTS.
                    throw NotImplementedError(format("SBTL PT_INPUTS query (T=%g K, P=%g Pa) is inside the "
                                                     "two-phase dome (T_sat=%g K).  PT_INPUTS in the dome "
                                                     "is ambiguous; use PQ_INPUTS or HmassP_INPUTS instead, "
                                                     "or specify_phase(iphase_liquid|iphase_gas) before "
                                                     "the update to disambiguate.",
                                                     T_query, p, T_sat_L));
                }
            }
        }

        if (tbl) {
            try {
                const double xnorm = tbl->xnorm_from_T(T_query, p, T_sat_active);
                if (std::isfinite(xnorm) && xnorm >= 0.0 && xnorm <= 1.0) {
                    std::size_t i = 0, j = 0;
                    find_native_nearest_good_indices(*tbl, *coeffs, xnorm, p, i, j);
                    const CellCoeffs& cell = (*coeffs)[i][j];
                    if (cell.valid()) {
                        _T = T_query;
                        _p = p;
                        const double rho_val = evaluate_single_phase(*tbl, *coeffs, iDmolar, xnorm, p, i, j);
                        if (!ValidNumber(rho_val) || rho_val <= 0.0) {
                            throw ValueError("SBTL normpt polynomial gave non-physical rho");
                        }
                        _rhomolar = rho_val;
                        _hmolar = evaluate_single_phase(*tbl, *coeffs, iHmolar, xnorm, p, i, j);
                        _smolar = evaluate_single_phase(*tbl, *coeffs, iSmolar, xnorm, p, i, j);
                        _umolar = evaluate_single_phase(*tbl, *coeffs, iUmolar, xnorm, p, i, j);
                        // _Q tag: 0 at LIQUID xnorm=1 boundary (saturated liquid),
                        // 1 at VAPOR xnorm=0 boundary (saturated vapor),
                        // -1000 elsewhere (single-phase interior or supercritical).
                        // Tolerance is 1e-10 to absorb routine f.p. noise from
                        // the (T - T_sat) subtraction in xnorm_from_T; tighter
                        // tolerances miss the boundary case where the user
                        // computes T_sat slightly differently than we did.
                        if (tbl == &_normpt_liquid && std::abs(xnorm - 1.0) < 1e-10) {
                            _Q = 0;
                        } else if (tbl == &_normpt_vapor && std::abs(xnorm) < 1e-10) {
                            _Q = 1;
                        } else {
                            _Q = -1000;
                        }
                        // Phase tag: prefer the user's imposed phase when they
                        // imposed one (their intent is documented).  Otherwise
                        // derive from the active table.
                        if (imposed_phase_index != iphase_not_imposed && imposed_phase_index != iphase_critical_point) {
                            _phase = imposed_phase_index;
                        } else if (tbl == &_normpt_super) {
                            _phase = (T_query > static_cast<double>(this->AS->T_critical())) ? iphase_supercritical_gas : iphase_supercritical_liquid;
                        } else if (tbl == &_normpt_liquid) {
                            _phase = iphase_liquid;
                        } else {
                            _phase = iphase_gas;
                        }
                        using_single_phase_table = true;
                        selected_table = SELECTED_PT_TABLE;
                        cached_single_phase_i = i;
                        cached_single_phase_j = j;
                        _normpt_active_table = tbl;
                        _normpt_active_coeffs = coeffs;
                        _normpt_xnorm = xnorm;
                        return;
                    }
                }
            } catch (const NotImplementedError&) {
                throw;  // Re-throw user-facing errors (in-dome PT_INPUTS).
            } catch (...) {
                // Internal failure: re-throw as a clean error so the caller
                // sees a failure rather than silent fall-through (SBTL is
                // pure-fluid PT/PH only — no fall-back path).
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa, T=%g K could not be evaluated through "
                                        "the normalized PT path; the table may not cover this state.  "
                                        "Use HEOS directly.",
                                        val1, val2));
            }
        }
    }

    // No table-driven fast path matched.  SBTL only supports PT_INPUTS,
    // HmolarP_INPUTS, and HmassP_INPUTS for pure fluids.  Anything else is
    // a user error in this configuration — surface a clean message.
    throw ValueError(format("SBTL: input pair %s not supported (only PT_INPUTS, HmassP_INPUTS, HmolarP_INPUTS).  "
                            "Use HEOS for other input pairs.",
                            get_input_pair_short_desc(input_pair).c_str()));
}

// ---------------------------------------------------------------------------
// Property accessors: return critbox-cached or normph-routed values
// ---------------------------------------------------------------------------
CoolPropDbl SBTLBackend::calc_hmolar(void) {
    if (_critbox_active) return static_cast<CoolPropDbl>(_hmolar);
    if (_normph_active_table) return static_cast<CoolPropDbl>(_hmolar);
    if (_normpt_active_table) return static_cast<CoolPropDbl>(_hmolar);
    return TabularBackend::calc_hmolar();
}
CoolPropDbl SBTLBackend::calc_smolar(void) {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_smolar);
    if (_normph_active_table) {
        return static_cast<CoolPropDbl>(evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, iSmolar, _normph_xnorm,
                                                              static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
    }
    if (_normpt_active_table) return static_cast<CoolPropDbl>(_smolar);
    return TabularBackend::calc_smolar();
}
CoolPropDbl SBTLBackend::calc_umolar(void) {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_umolar);
    if (_normph_active_table) {
        return static_cast<CoolPropDbl>(evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, iUmolar, _normph_xnorm,
                                                              static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
    }
    if (_normpt_active_table) return static_cast<CoolPropDbl>(_umolar);
    return TabularBackend::calc_umolar();
}
CoolPropDbl SBTLBackend::calc_rhomolar(void) {
    if (_critbox_active) return static_cast<CoolPropDbl>(_rhomolar);
    if (_normph_active_table) return static_cast<CoolPropDbl>(_rhomolar);
    if (_normpt_active_table) return static_cast<CoolPropDbl>(_rhomolar);
    return TabularBackend::calc_rhomolar();
}

// cpmolar / cvmolar / speed_sound: SBTL does not build polynomial
// derivatives over its normalized tables (the build path skips the
// xkey/ykey polynomials, and the xnorm chain-rule scaling
// dxi/dT = 1/(T_hi-T_lo) varies per isobar in ways the base derivative
// evaluator wouldn't capture).  Route everything through HEOS direct
// after priming AS at the requested state — ~10 µs/call, same envelope
// as the critbox path.
// ensure_heos_at_active_state: prime the underlying HEOS at the state the
// user actually requested.  For an active normpt table the user supplied
// (T, p) directly and we use PT_INPUTS.  For an active normph table the
// user supplied (h, p) and we use HmolarP_INPUTS — calling PT_INPUTS with
// the interpolated _T would re-evaluate HEOS at a slightly-off T (table
// interpolation error ~1e-5 relative), giving cp/cv/c_sound values that
// disagree with the EOS-direct answer at the *requested* (h, p) by up to
// several percent in regions where cp(T) is steep (e.g. supercritical
// shoulder).  Using the right input pair keeps the state consistent
// between the SBTL-evaluated h/s/u and the HEOS-evaluated cp/cv/c_sound.
void SBTLBackend::ensure_heos_at_active_state() {
    if (_normph_active_table != nullptr) {
        this->AS->update(HmolarP_INPUTS, static_cast<double>(_hmolar), static_cast<double>(_p));
    } else if (_normpt_active_table != nullptr) {
        this->AS->update(PT_INPUTS, static_cast<double>(_p), static_cast<double>(_T));
    } else {
        // Fallback: most accurate input pair available is PT (we have it).
        this->AS->update(PT_INPUTS, static_cast<double>(_p), static_cast<double>(_T));
    }
}
// First-partial-derivative routing: the base-class calc_first_partial_deriv
// fills the Jacobian by calling evaluate_single_phase_*_derivative(output, ...)
// for output ∈ {Of, Wrt, Constant}.  For an active normalized table this
// breaks when output == table.xkey or table.ykey because the build path
// SKIPS the xkey/ykey polynomials (they're recoverable from the boundary
// functions, not stored as cell coefficients), so the polynomial dispatch
// would dereference an empty alpha vector.  Routing the whole derivative
// through HEOS direct here uses ~10 µs total — same envelope as the
// individual cp/cv/c_sound overrides — and gives the exact answer at the
// active state.
CoolPropDbl SBTLBackend::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
    if (_critbox_active || _normph_active_table || _normpt_active_table) {
        ensure_heos_at_active_state();
        return this->AS->first_partial_deriv(Of, Wrt, Constant);
    }
    return TabularBackend::calc_first_partial_deriv(Of, Wrt, Constant);
}
CoolPropDbl SBTLBackend::calc_cpmolar() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_cpmolar);
    if (_normph_active_table || _normpt_active_table) {
        ensure_heos_at_active_state();
        return static_cast<CoolPropDbl>(this->AS->cpmolar());
    }
    return TabularBackend::calc_cpmolar();
}
CoolPropDbl SBTLBackend::calc_cvmolar() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_cvmolar);
    if (_normph_active_table || _normpt_active_table) {
        ensure_heos_at_active_state();
        return static_cast<CoolPropDbl>(this->AS->cvmolar());
    }
    return TabularBackend::calc_cvmolar();
}
CoolPropDbl SBTLBackend::calc_speed_sound() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_speed_sound);
    // SBTL native: evaluate the speed-of-sound polynomial directly from
    // the active normalized table (matches Kunick's IAPWS-IF97 design —
    // each spline function has its own coefficient block).  ~1 µs vs
    // ~10 µs for HEOS direct, and consistent with the rest of the SBTL
    // output (no formulation-comparison contamination when SBTL was
    // built from IF97 and the caller queries against IF97).  Cells whose
    // EOS corners returned NaN for w get no aux polynomial built — fall
    // back to the base class (which routes via HEOS direct) in that case.
    if (_normph_active_table != nullptr && _normph_active_coeffs != nullptr) {
        const auto& cell = (*_normph_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(ispeed_sound).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, ispeed_sound, _normph_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    if (_normpt_active_table != nullptr && _normpt_active_coeffs != nullptr) {
        const auto& cell = (*_normpt_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(ispeed_sound).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normpt_active_table, *_normpt_active_coeffs, ispeed_sound, _normpt_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    return TabularBackend::calc_speed_sound();
}
CoolPropDbl SBTLBackend::calc_viscosity() {
    if (_critbox_active) {
        if (!ValidNumber(_critbox_viscosity)) {
            throw ValueError("SBTL critbox viscosity not available for this state");
        }
        return static_cast<CoolPropDbl>(_critbox_viscosity);
    }
    if (_normph_active_table != nullptr && _normph_active_coeffs != nullptr) {
        const auto& cell = (*_normph_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(iviscosity).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, iviscosity, _normph_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    if (_normpt_active_table != nullptr && _normpt_active_coeffs != nullptr) {
        const auto& cell = (*_normpt_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(iviscosity).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normpt_active_table, *_normpt_active_coeffs, iviscosity, _normpt_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    return TabularBackend::calc_viscosity();
}
CoolPropDbl SBTLBackend::calc_conductivity() {
    if (_critbox_active) {
        if (!ValidNumber(_critbox_conductivity)) {
            throw ValueError("SBTL critbox conductivity not available for this state");
        }
        return static_cast<CoolPropDbl>(_critbox_conductivity);
    }
    if (_normph_active_table != nullptr && _normph_active_coeffs != nullptr) {
        const auto& cell = (*_normph_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(iconductivity).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normph_active_table, *_normph_active_coeffs, iconductivity, _normph_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    if (_normpt_active_table != nullptr && _normpt_active_coeffs != nullptr) {
        const auto& cell = (*_normpt_active_coeffs)[cached_single_phase_i][cached_single_phase_j];
        if (!cell.get(iconductivity).empty()) {
            return static_cast<CoolPropDbl>(evaluate_single_phase(*_normpt_active_table, *_normpt_active_coeffs, iconductivity, _normpt_xnorm,
                                                                  static_cast<double>(_p), cached_single_phase_i, cached_single_phase_j));
        }
    }
    return TabularBackend::calc_conductivity();
}

}  // namespace CoolProp

#endif  // !defined(NO_TABULAR_BACKENDS)
