#if !defined(NO_TABULAR_BACKENDS)

#    include "SBTLBackend.h"
#    include "DataStructures.h"
#    include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#    include <algorithm>
#    include <iostream>
#    include <fstream>
#    include <map>
#    include <memory>
#    include <mutex>
#    include <utility>
#    include <Eigen/Sparse>
#    include <msgpack.hpp>
#    include "miniz.h"
#    include "CPfilepaths.h"
#    include "Configuration.h"
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
    // Thread-safety: this is a process-wide cache (file-static via
    // bspline_cache()).  If two threads concurrently construct SBTL
    // backends with different grid sizes, both could race on solvers_.
    // Magic-statics gives us thread-safe initialization of the cache
    // object itself, but the map mutation here needs its own lock.
    // Contention is negligible — the lookup is hit once per (Nx, Ny)
    // pair across the program's lifetime, then becomes a const map
    // read.
    Solver& solver_for(std::size_t N) {
        std::lock_guard<std::mutex> lock(solvers_mutex_);
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
    std::mutex solvers_mutex_;
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
const std::array<std::array<double, 4>, 4> Bspline_poly = {{
  {1.0 / 6.0, -1.0 / 2.0, 1.0 / 2.0, -1.0 / 6.0},
  {2.0 / 3.0, 0.0, -1.0, 1.0 / 2.0},
  {1.0 / 6.0, 1.0 / 2.0, 1.0 / 2.0, -1.0 / 2.0},
  {0.0, 0.0, 0.0, 1.0 / 6.0},
}};

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
// normph backbone).  Wired into build_bspline_coeffs as the primary
// alpha source for the 4 core derived params on normph / normpt
// tables; cells where corner partials throw fall back to cubic alpha.
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
// Corner-state derivative collection for Hermite bicubic build path
// (CoolProp-foi.1).  Reads value + 4 partials per core property (rho, T,
// s, u) from a HEOS state primed by the caller.  Returns false if any
// EOS call throws or yields non-finite, so the build loop can skip
// unusable corners (sat boundary, metastable spinodal).
// ---------------------------------------------------------------------------
bool populate_corner_derivatives(AbstractState& AS, SBTLCornerDerivs& out) {
    auto safe = [&](double v) { return ValidNumber(v) ? v : std::numeric_limits<double>::quiet_NaN(); };
    try {
        out.rhomolar = safe(static_cast<double>(AS.rhomolar()));
        out.T = safe(static_cast<double>(AS.T()));
        out.smolar = safe(static_cast<double>(AS.smolar()));
        out.umolar = safe(static_cast<double>(AS.umolar()));
        if (!ValidNumber(out.rhomolar) || !ValidNumber(out.T) || !ValidNumber(out.smolar) || !ValidNumber(out.umolar)) {
            return false;
        }

        // First partials at constant p (∂f/∂h).
        out.drho_dh = safe(static_cast<double>(AS.first_partial_deriv(iDmolar, iHmolar, iP)));
        out.dT_dh = safe(static_cast<double>(AS.first_partial_deriv(iT, iHmolar, iP)));
        out.ds_dh = safe(static_cast<double>(AS.first_partial_deriv(iSmolar, iHmolar, iP)));
        out.du_dh = safe(static_cast<double>(AS.first_partial_deriv(iUmolar, iHmolar, iP)));

        // First partials at constant h (∂f/∂p).
        out.drho_dp = safe(static_cast<double>(AS.first_partial_deriv(iDmolar, iP, iHmolar)));
        out.dT_dp = safe(static_cast<double>(AS.first_partial_deriv(iT, iP, iHmolar)));
        out.ds_dp = safe(static_cast<double>(AS.first_partial_deriv(iSmolar, iP, iHmolar)));
        out.du_dp = safe(static_cast<double>(AS.first_partial_deriv(iUmolar, iP, iHmolar)));

        // Second partials at constant p (∂²f/∂h²).
        out.d2rho_dh2 = safe(static_cast<double>(AS.second_partial_deriv(iDmolar, iHmolar, iP, iHmolar, iP)));
        out.d2T_dh2 = safe(static_cast<double>(AS.second_partial_deriv(iT, iHmolar, iP, iHmolar, iP)));
        out.d2s_dh2 = safe(static_cast<double>(AS.second_partial_deriv(iSmolar, iHmolar, iP, iHmolar, iP)));
        out.d2u_dh2 = safe(static_cast<double>(AS.second_partial_deriv(iUmolar, iHmolar, iP, iHmolar, iP)));

        // Mixed second partials ∂²f/∂h∂p (d/dp of ∂f/∂h|_p, at constant h).
        out.d2rho_dhp = safe(static_cast<double>(AS.second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar)));
        out.d2T_dhp = safe(static_cast<double>(AS.second_partial_deriv(iT, iHmolar, iP, iP, iHmolar)));
        out.d2s_dhp = safe(static_cast<double>(AS.second_partial_deriv(iSmolar, iHmolar, iP, iP, iHmolar)));
        out.d2u_dhp = safe(static_cast<double>(AS.second_partial_deriv(iUmolar, iHmolar, iP, iP, iHmolar)));
    } catch (...) {
        return false;
    }
    out.valid = ValidNumber(out.drho_dh) && ValidNumber(out.drho_dp) && ValidNumber(out.d2rho_dh2) && ValidNumber(out.d2rho_dhp);
    return out.valid;
}

// PT-table analogue.  Reads ∂f/∂T, ∂f/∂p, ∂²f/∂T², ∂²f/∂T∂p for
// f in {rho, h, s, u} from a primed HEOS state.  Same usage contract
// as populate_corner_derivatives — returns false on any throw or NaN.
bool populate_corner_derivatives_pt(AbstractState& AS, SBTLCornerDerivsPT& out) {
    auto safe = [&](double v) { return ValidNumber(v) ? v : std::numeric_limits<double>::quiet_NaN(); };
    try {
        out.rhomolar = safe(static_cast<double>(AS.rhomolar()));
        out.hmolar = safe(static_cast<double>(AS.hmolar()));
        out.smolar = safe(static_cast<double>(AS.smolar()));
        out.umolar = safe(static_cast<double>(AS.umolar()));
        if (!ValidNumber(out.rhomolar) || !ValidNumber(out.hmolar) || !ValidNumber(out.smolar) || !ValidNumber(out.umolar)) {
            return false;
        }
        out.drho_dT = safe(static_cast<double>(AS.first_partial_deriv(iDmolar, iT, iP)));
        out.dh_dT = safe(static_cast<double>(AS.first_partial_deriv(iHmolar, iT, iP)));
        out.ds_dT = safe(static_cast<double>(AS.first_partial_deriv(iSmolar, iT, iP)));
        out.du_dT = safe(static_cast<double>(AS.first_partial_deriv(iUmolar, iT, iP)));

        out.drho_dp = safe(static_cast<double>(AS.first_partial_deriv(iDmolar, iP, iT)));
        out.dh_dp = safe(static_cast<double>(AS.first_partial_deriv(iHmolar, iP, iT)));
        out.ds_dp = safe(static_cast<double>(AS.first_partial_deriv(iSmolar, iP, iT)));
        out.du_dp = safe(static_cast<double>(AS.first_partial_deriv(iUmolar, iP, iT)));

        out.d2rho_dT2 = safe(static_cast<double>(AS.second_partial_deriv(iDmolar, iT, iP, iT, iP)));
        out.d2h_dT2 = safe(static_cast<double>(AS.second_partial_deriv(iHmolar, iT, iP, iT, iP)));
        out.d2s_dT2 = safe(static_cast<double>(AS.second_partial_deriv(iSmolar, iT, iP, iT, iP)));
        out.d2u_dT2 = safe(static_cast<double>(AS.second_partial_deriv(iUmolar, iT, iP, iT, iP)));

        out.d2rho_dTp = safe(static_cast<double>(AS.second_partial_deriv(iDmolar, iT, iP, iP, iT)));
        out.d2h_dTp = safe(static_cast<double>(AS.second_partial_deriv(iHmolar, iT, iP, iP, iT)));
        out.d2s_dTp = safe(static_cast<double>(AS.second_partial_deriv(iSmolar, iT, iP, iP, iT)));
        out.d2u_dTp = safe(static_cast<double>(AS.second_partial_deriv(iUmolar, iT, iP, iP, iT)));
    } catch (...) {
        return false;
    }
    out.valid = ValidNumber(out.drho_dT) && ValidNumber(out.drho_dp) && ValidNumber(out.d2rho_dT2) && ValidNumber(out.d2rho_dTp);
    return out.valid;
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

    // For LIQUID/VAPOR only: probe whether HEOS PT_INPUTS evaluates cleanly
    // at the corners of the table's pressure envelope.  Some cryogens
    // (Argon, Helium) have an HEOS melting-line ancillary whose lower-p
    // bound sits slightly above p_triple — AS->update at p = p_triple then
    // throws "unable to calculate melting line T(p)".  Left in place the
    // failing row contaminates Cheb1D sample values (via safe_h_at_PT
    // returning _HUGE) and all subsequent table cells become invalid.
    // SUPER region is intentionally exempted: at p = p_crit and T = T_min
    // the fluid is below the melting line, so the probe would always fail
    // for any cryogen and we'd walk ymin far above p_crit, opening a
    // [p_crit, ymin_walked] hole in supercritical coverage.
    if (region_ != SUPER) {
        const double T_min_HEOS = std::max(static_cast<double>(AS->Ttriple()), static_cast<double>(AS->Tmin()));
        const double T_max_HEOS = static_cast<double>(AS->Tmax());
        auto probe_ok = [&](double P) -> bool {
            for (double T : {T_min_HEOS, T_max_HEOS}) {
                try {
                    AS->update(PT_INPUTS, P, T);
                } catch (...) {
                    return false;
                }
            }
            return true;
        };
        for (int trial = 0; trial < 30 && ymin < ymax; ++trial) {
            if (probe_ok(ymin)) break;
            ymin *= 1.05;
        }
    }
}

double NormalizedPHTable::xnorm_from_h(double h, double P, std::optional<double> h_sat) const {
    // Source of truth per region:
    //   LIQUID: h_lo from Cheb (h(T_min, p)),     h_hi from H-superanc (h_sat,L)
    //   VAPOR : h_lo from H-superanc (h_sat,V),    h_hi from Cheb (h(T_max, p))
    //   SUPER : both from Cheb
    // Caller passes h_sat (machine-precision via H-superancillary) for the
    // sat-boundary; if nullopt, fall back to the Cheb (less accurate near
    // the critical cusp, but still works for smooth-monotonic stretches).
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
    double h_lo = 0.0, h_hi = 0.0;
    switch (region_) {
        case LIQUID:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = h_sat ? *h_sat : cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case VAPOR:
            h_lo = h_sat ? *h_sat : cheb_or_isobar(h_lo_cheb, h_lo_isobar);
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
    //
    // For cryogens (Argon, Helium, Neon) the HEOS melting-line ancillary's
    // lower-p bound sits just above the triple-point pressure — sampling
    // Cheb-Lobatto nodes there returns _HUGE through safe_h_at_PT and the
    // resulting polynomial evaluates to NaN at intermediate points
    // (huge - huge cancellation).  Walk p_lo upward until both h_lo_fn
    // (h at T_lo) and h_hi_fn (h at T_max_ext) return finite values, so
    // every Lobatto node feeds a clean number into the build.
    double p_lo = table.yvec.front() * 1.0001;
    const double p_hi = table.yvec.back() * 0.9999;
    {
        const bool need_h_lo = (table.region() != NormalizedPHTable::VAPOR);
        const bool need_h_hi = (table.region() != NormalizedPHTable::LIQUID);
        for (int trial = 0; trial < 40; ++trial) {
            const double T_melt_p = T_melt_floor(p_lo);
            const double T_lo_p = std::isfinite(T_melt_p) ? T_melt_p : T_min;
            const double v_lo = need_h_lo ? safe_h_at_PT(p_lo, T_lo_p) : 0.0;
            const double v_hi = need_h_hi ? safe_h_at_PT(p_lo, T_max_ext) : 0.0;
            const bool ok_lo = std::isfinite(v_lo) && std::abs(v_lo) < 0.5 * _HUGE;
            const bool ok_hi = std::isfinite(v_hi) && std::abs(v_hi) < 0.5 * _HUGE;
            if (ok_lo && ok_hi) break;
            p_lo *= 1.05;
            if (p_lo >= p_hi) {
                p_lo = table.yvec.front() * 1.0001;
                break;
            }
        }
    }
    // Build Chebyshev expansions only for the NON-saturation boundaries:
    //   LIQUID: h_lo = h(T_min, p)        — Cheb
    //   LIQUID: h_hi = h_sat,L(p)         — H-superancillary at lookup
    //   VAPOR : h_lo = h_sat,V(p)         — H-superancillary at lookup
    //   VAPOR : h_hi = h(T_max_ext, p)    — Cheb
    //   SUPER : both                      — Cheb
    // Saturation enthalpies are already machine-precision via the
    // H-superancillaries (Bell et al.); re-fitting with Chebyshev would
    // just inherit any error there with no upside.
    // Piecewise Chebyshev with low degree per piece (8) and more pieces near
    // p_crit.  Replaces an earlier 4-piece × degree-32 layout that was both
    // slow to evaluate (33-coeff Clenshaw on every property access) and
    // poorly conditioned for the sqrt-cusp in h_sat as p -> p_crit.
    // Fewer coeffs per piece + cusp-concentrated breakpoints gives:
    //   - faster Cheb1D.eval (9-coeff Clenshaw, ~3x less arithmetic),
    //   - smaller stored polynomials (10*9=90 vs 4*33=132 doubles),
    //   - better near-critical resolution where the cusp dominates.
    std::vector<double> breakpoints;
    if (table.region() == NormalizedPHTable::SUPER) {
        breakpoints = {p_lo, p_hi};  // smooth functions, single piece suffices
    } else if (this->AS) {
        const double p_crit = static_cast<double>(this->AS->p_critical());
        breakpoints = {p_lo,          0.10 * p_crit, 0.30 * p_crit,  0.50 * p_crit,  0.70 * p_crit, 0.85 * p_crit,
                       0.92 * p_crit, 0.96 * p_crit, 0.985 * p_crit, 0.995 * p_crit, p_hi};
        // Drop any breakpoint that doesn't fit between p_lo and p_hi
        // (can happen for narrow tables or if p_lo got walked up to be
        // above one of the fixed fractions).
        breakpoints.erase(std::remove_if(breakpoints.begin() + 1, breakpoints.end() - 1, [&](double v) { return v <= p_lo || v >= p_hi; }),
                          breakpoints.end() - 1);
    } else {
        breakpoints = {p_lo, p_hi};
    }
    const std::size_t N_cheb = 8;
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

std::pair<std::vector<double>, std::vector<double>>
  SBTLBackend::build_adaptive_normph_grid(int region_int, double ymin, double ymax, std::size_t Nx_max, std::size_t Ny_max, double tol_relerr) const {
    // Joint xvec+yvec adaptive grid by alternating row/column bisection.
    //
    // Algorithm:
    //   1. Seed with coarse uniform xvec ([0.01, 0.99]) and log-uniform
    //      yvec ([ymin, ymax]).
    //   2. Loop, alternating two passes:
    //      a) Row pass: for each yvec[j..j+1] pair, evaluate the actual
    //         Hermite bicubic at every cell midpoint along that row
    //         (using current xvec).  If ANY cell exceeds tol_relerr,
    //         insert the geometric midpoint of yvec[j..j+1] into yvec.
    //      b) Column pass: same, transposed — for each xvec[i..i+1] pair,
    //         scan all cells in the column against current yvec and
    //         bisect the column if any fail.
    //   3. Stop when a full row+column pass produces no new splits, or
    //      when budgets Nx_max / Ny_max are hit.
    //
    // The acceptance criterion is the production cell residual at the
    // cell midpoint: build the same Hermite bicubic the cell will store
    // (corner ρ + HEOS first partials + cross deriv, chain-ruled into
    // (xnorm, log p) cell coords and through the active output transform),
    // evaluate at (ξ=0.5, η=0.5), invert the transform, compare to HEOS
    // truth at the physical midpoint.  This is the residual the user
    // sees at lookup — no proxy.
    //
    // Falls back to log-uniform on any setup failure.
    auto log_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        const double la = std::log(a), lb = std::log(b);
        for (std::size_t k = 0; k < n; ++k) {
            const double frac = static_cast<double>(k) / static_cast<double>(n - 1);
            v.push_back(std::exp(la + frac * (lb - la)));
        }
        return v;
    };
    auto linear_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            const double frac = static_cast<double>(k) / static_cast<double>(n - 1);
            v.push_back(a + frac * (b - a));
        }
        return v;
    };
    constexpr std::size_t Nx_seed = 16;
    constexpr std::size_t Ny_seed = 16;
    // Inset slightly inside [0, 1] so HmolarP_INPUTS at the seed xnorm=0/1
    // doesn't sit exactly on the sat boundary (which the dome flash can
    // wobble around).  Final xvec stretches to [0, 1] via the inset
    // anchors getting bisected outward in the early refinement passes
    // (not strictly necessary — we anchor inside).
    std::vector<double> xvec = linear_uniform(0.0, 1.0, Nx_seed);
    std::vector<double> yvec = log_uniform(ymin, ymax, Ny_seed);

    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (!helm) return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};
    double rho_crit_local = 0.0;
    try {
        helm->ensure_caloric_superancillaries();
        rho_crit_local = static_cast<double>(helm->rhomolar_critical());
    } catch (...) {
        return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};
    }
    if (rho_crit_local <= 0.0) return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};

    std::shared_ptr<CoolProp::AbstractState> probe_AS;
    try {
        const auto names = helm->fluid_names();
        if (names.empty()) return {std::move(xvec), std::move(yvec)};
        probe_AS.reset(CoolProp::AbstractState::factory("HEOS", names[0]));
    } catch (...) {
        return {std::move(xvec), std::move(yvec)};
    }

    const double T_min_eos = std::max(static_cast<double>(probe_AS->Ttriple()), static_cast<double>(probe_AS->Tmin()));
    const double T_max_ext = static_cast<double>(probe_AS->Tmax()) * 1.499;
    const bool is_vapor = (region_int == 1);
    const bool is_liquid = (region_int == 0);
    const bool is_super = (region_int == 2);
    if (!is_vapor && !is_liquid && !is_super) return {std::move(xvec), std::move(yvec)};

    // Region-specific h bounds & their dp derivatives at a given p.
    // SUPER uses T_lo = T_min_eos and T_hi = T_max_ext (both constants in T)
    // as an approximation for the cold/hot isotherm boundaries.  This
    // ignores the melting-curve floor T_melt(P) that production uses for
    // LIQUID-side T_lo in SUPER — acceptable for shallow-melting fluids
    // (Argon, CO2, R245fa, D6).  For aggressive-melting fluids (Water at
    // very high P), the cold corner sits below the melting curve in HEOS
    // extrapolation territory; the resulting adaptive grid may have
    // slightly different cell layout than what production stores there,
    // but the worst-case states near the critical pillar (the actual
    // accuracy bottleneck) are unaffected.
    auto bounds_at_p = [&probe_AS, T_min_eos, T_max_ext, is_vapor, is_liquid](double p, double& h_lo, double& dh_lo_dp, double& h_hi,
                                                                              double& dh_hi_dp) -> bool {
        try {
            if (is_vapor) {
                probe_AS->update(PQ_INPUTS, p, 1.0);
                h_lo = static_cast<double>(probe_AS->hmolar());
                dh_lo_dp = probe_AS->first_saturation_deriv(iHmolar, iP);
                probe_AS->update(PT_INPUTS, p, T_max_ext);
                h_hi = static_cast<double>(probe_AS->hmolar());
                dh_hi_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
            } else if (is_liquid) {
                probe_AS->update(PT_INPUTS, p, T_min_eos);
                h_lo = static_cast<double>(probe_AS->hmolar());
                dh_lo_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
                probe_AS->update(PQ_INPUTS, p, 0.0);
                h_hi = static_cast<double>(probe_AS->hmolar());
                dh_hi_dp = probe_AS->first_saturation_deriv(iHmolar, iP);
            } else {  // SUPER
                probe_AS->update(PT_INPUTS, p, T_min_eos);
                h_lo = static_cast<double>(probe_AS->hmolar());
                dh_lo_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
                probe_AS->update(PT_INPUTS, p, T_max_ext);
                h_hi = static_cast<double>(probe_AS->hmolar());
                dh_hi_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
            }
            return std::isfinite(h_lo) && std::isfinite(h_hi) && std::isfinite(dh_lo_dp) && std::isfinite(dh_hi_dp);
        } catch (...) {
            return false;
        }
    };

    // Full HEOS corner data at (xnorm, p): ρ + h-partials + bounds + their dp derivatives.
    struct Corner
    {
        double f;     // ρ
        double f_h;   // ∂ρ/∂h|p
        double f_p;   // ∂ρ/∂p|h
        double f_hh;  // ∂²ρ/∂h²
        double f_hp;  // ∂²ρ/∂h∂p
        double h_lo;  // table low-h boundary at this p
        double h_hi;  // table high-h boundary at this p
        double dh_lo_dp;
        double dh_hi_dp;
        double p;
        bool ok;
    };
    auto compute_corner = [&probe_AS, &bounds_at_p](double xnorm, double p) -> Corner {
        Corner c{};
        c.p = p;
        if (!bounds_at_p(p, c.h_lo, c.dh_lo_dp, c.h_hi, c.dh_hi_dp)) {
            c.ok = false;
            return c;
        }
        const double h = c.h_lo + xnorm * (c.h_hi - c.h_lo);
        try {
            probe_AS->update(HmolarP_INPUTS, h, p);
            c.f = static_cast<double>(probe_AS->rhomolar());
            c.f_h = probe_AS->first_partial_deriv(iDmolar, iHmolar, iP);
            c.f_p = probe_AS->first_partial_deriv(iDmolar, iP, iHmolar);
            c.f_hh = probe_AS->second_partial_deriv(iDmolar, iHmolar, iP, iHmolar, iP);
            c.f_hp = probe_AS->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);
            c.ok = std::isfinite(c.f) && std::isfinite(c.f_h) && std::isfinite(c.f_p) && std::isfinite(c.f_hh) && std::isfinite(c.f_hp);
        } catch (...) {
            c.ok = false;
        }
        return c;
    };

    // HEOS truth ρ at a query (xnorm, p) for the midpoint comparison.
    auto rho_truth = [&probe_AS, &bounds_at_p](double xnorm, double p) -> double {
        double h_lo = 0.0, dh_lo_dp = 0.0, h_hi = 0.0, dh_hi_dp = 0.0;
        if (!bounds_at_p(p, h_lo, dh_lo_dp, h_hi, dh_hi_dp)) return std::numeric_limits<double>::quiet_NaN();
        const double h = h_lo + xnorm * (h_hi - h_lo);
        try {
            probe_AS->update(HmolarP_INPUTS, h, p);
            return static_cast<double>(probe_AS->rhomolar());
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    };

    // Cell midpoint residual: build the production-style Hermite bicubic,
    // evaluate at (ξ=0.5, η=0.5), invert the transform, compare to HEOS.
    const SBTLTransform tform = sbtl_transform_for_property(region_int, iDmolar);
    const double y_ref = sbtl_transform_y_ref(region_int, iDmolar, rho_crit_local);
    // Same cusp-distance guard production uses: if any corner sits within
    // 1e-6·y_ref of the reference (ρ_crit for CUSP_DIST), the chain-rule
    // derivative g'(f) = 0.5/√d blows up; production falls back to B-spline
    // for that cell, so the adaptive metric should not try to drive
    // bisection there (return 0, treating the cell as "not Hermite's
    // responsibility").
    const double cusp_min_d =
      (tform == SBTLTransform::CUSP_DIST_ABOVE_REF || tform == SBTLTransform::CUSP_DIST_BELOW_REF) ? 1e-6 * std::abs(y_ref) : 0.0;
    auto cell_relerr = [&rho_truth, rho_crit_local, tform, y_ref, cusp_min_d](const Corner& c00, const Corner& c10, const Corner& c01,
                                                                              const Corner& c11, double xn_a, double xn_b, double p_a,
                                                                              double p_b) -> double {
        const double Delta_xi = xn_b - xn_a;
        const double dlog_p = std::log(p_b) - std::log(p_a);
        const std::array<const Corner*, 4> cs = {&c00, &c10, &c01, &c11};
        const std::array<double, 4> xn = {xn_a, xn_b, xn_a, xn_b};
        std::array<double, 4> f{}, fxi{}, feta{}, fxieta{};
        for (int k = 0; k < 4; ++k) {
            const Corner& c = *cs[k];
            if (!c.ok) return std::numeric_limits<double>::infinity();
            if (tform == SBTLTransform::CUSP_DIST_ABOVE_REF && (c.f - y_ref) <= cusp_min_d) return 0.0;
            if (tform == SBTLTransform::CUSP_DIST_BELOW_REF && (y_ref - c.f) <= cusp_min_d) return 0.0;
            const double Delta_h_p = c.h_hi - c.h_lo;
            const double dDelta_h_dlogp = c.p * (c.dh_hi_dp - c.dh_lo_dp);
            const double dh_dlogp_at_xn = c.p * (c.dh_lo_dp + xn[k] * (c.dh_hi_dp - c.dh_lo_dp));
            double df_dxnorm = c.f_h * Delta_h_p;
            double df_dlogp = c.p * c.f_p + c.f_h * dh_dlogp_at_xn;
            double d2f = Delta_h_p * (c.f_hh * dh_dlogp_at_xn + c.p * c.f_hp) + c.f_h * dDelta_h_dlogp;
            double fval = c.f;
            if (tform != SBTLTransform::IDENTITY) {
                if (!(fval > 0.0)) return std::numeric_limits<double>::infinity();
                const double gp = sbtl_transform_deriv(tform, fval, y_ref);
                const double gpp = sbtl_transform_second_deriv(tform, fval, y_ref);
                const double new_d2 = gpp * df_dxnorm * df_dlogp + gp * d2f;
                df_dxnorm *= gp;
                df_dlogp *= gp;
                d2f = new_d2;
                fval = sbtl_transform_apply(tform, fval, y_ref);
            }
            f[k] = fval;
            fxi[k] = df_dxnorm * Delta_xi;
            feta[k] = df_dlogp * dlog_p;
            fxieta[k] = d2f * Delta_xi * dlog_p;
        }
        const auto alpha = hermite_bicubic_polynomial_coeffs(f[0], f[1], f[2], f[3], fxi[0], fxi[1], fxi[2], fxi[3], feta[0], feta[1], feta[2],
                                                             feta[3], fxieta[0], fxieta[1], fxieta[2], fxieta[3]);
        constexpr int N_test = 3;
        constexpr std::array<double, N_test> ts = {0.25, 0.5, 0.75};
        double worst_relerr = 0.0;
        for (int ti = 0; ti < N_test; ++ti) {
            for (int tj = 0; tj < N_test; ++tj) {
                const double xi = ts[ti], eta = ts[tj];
                const double xi2 = xi * xi, xi3 = xi2 * xi;
                const double eta2 = eta * eta, eta3 = eta2 * eta;
                const double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2 + alpha[3] * eta3;
                const double B1 = alpha[4] + alpha[5] * eta + alpha[6] * eta2 + alpha[7] * eta3;
                const double B2 = alpha[8] + alpha[9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
                const double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
                const double F_test = B0 + B1 * xi + B2 * xi2 + B3 * xi3;
                const double r_pred = sbtl_transform_inverse(tform, F_test, y_ref);
                const double xn_test = xn_a + xi * (xn_b - xn_a);
                const double p_test = std::exp(std::log(p_a) + eta * (std::log(p_b) - std::log(p_a)));
                const double r_heos = rho_truth(xn_test, p_test);
                if (!std::isfinite(r_pred) || !std::isfinite(r_heos)) return std::numeric_limits<double>::infinity();
                const double denom = std::max(std::abs(r_heos), rho_crit_local * 0.01);
                const double relerr = std::abs(r_pred - r_heos) / denom;
                if (relerr > worst_relerr) worst_relerr = relerr;
            }
        }
        return worst_relerr;
    };

    auto recompute_corners = [&compute_corner](const std::vector<double>& xv, const std::vector<double>& yv) -> std::vector<std::vector<Corner>> {
        std::vector<std::vector<Corner>> m(xv.size(), std::vector<Corner>(yv.size()));
        for (std::size_t i = 0; i < xv.size(); ++i) {
            for (std::size_t j = 0; j < yv.size(); ++j) {
                m[i][j] = compute_corner(xv[i], yv[j]);
            }
        }
        return m;
    };

    constexpr int max_iter = 14;
    for (int iter = 0; iter < max_iter; ++iter) {
        auto corners = recompute_corners(xvec, yvec);
        bool any_split = false;
        if (yvec.size() < Ny_max) {
            std::vector<double> new_yvec;
            new_yvec.reserve(yvec.size() * 2);
            for (std::size_t j = 0; j + 1 < yvec.size(); ++j) {
                new_yvec.push_back(yvec[j]);
                bool row_bad = false;
                for (std::size_t i = 0; i + 1 < xvec.size(); ++i) {
                    if (cell_relerr(corners[i][j], corners[i + 1][j], corners[i][j + 1], corners[i + 1][j + 1], xvec[i], xvec[i + 1], yvec[j],
                                    yvec[j + 1])
                        > tol_relerr) {
                        row_bad = true;
                        break;
                    }
                }
                if (row_bad && new_yvec.size() + (yvec.size() - j - 1) < Ny_max) {
                    new_yvec.push_back(std::exp(0.5 * (std::log(yvec[j]) + std::log(yvec[j + 1]))));
                    any_split = true;
                }
            }
            new_yvec.push_back(yvec.back());
            yvec = std::move(new_yvec);
        }
        if (xvec.size() < Nx_max) {
            corners = recompute_corners(xvec, yvec);
            std::vector<double> new_xvec;
            new_xvec.reserve(xvec.size() * 2);
            for (std::size_t i = 0; i + 1 < xvec.size(); ++i) {
                new_xvec.push_back(xvec[i]);
                bool col_bad = false;
                for (std::size_t j = 0; j + 1 < yvec.size(); ++j) {
                    if (cell_relerr(corners[i][j], corners[i + 1][j], corners[i][j + 1], corners[i + 1][j + 1], xvec[i], xvec[i + 1], yvec[j],
                                    yvec[j + 1])
                        > tol_relerr) {
                        col_bad = true;
                        break;
                    }
                }
                if (col_bad && new_xvec.size() + (xvec.size() - i - 1) < Nx_max) {
                    new_xvec.push_back(0.5 * (xvec[i] + xvec[i + 1]));
                    any_split = true;
                }
            }
            new_xvec.push_back(xvec.back());
            xvec = std::move(new_xvec);
        }
        if (!any_split) break;
    }

    return {std::move(xvec), std::move(yvec)};
}

std::pair<std::vector<double>, std::vector<double>>
  SBTLBackend::build_adaptive_normpt_grid(int region_int, double ymin, double ymax, std::size_t Nx_max, std::size_t Ny_max, double tol_relerr) const {
    // PT analogue of build_adaptive_normph_grid.  Same algorithm; T
    // replaces h as the independent variable on the xnorm axis.
    auto log_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        const double la = std::log(a), lb = std::log(b);
        for (std::size_t k = 0; k < n; ++k) {
            v.push_back(std::exp(la + (lb - la) * static_cast<double>(k) / static_cast<double>(n - 1)));
        }
        return v;
    };
    auto linear_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        for (std::size_t k = 0; k < n; ++k)
            v.push_back(a + (b - a) * static_cast<double>(k) / static_cast<double>(n - 1));
        return v;
    };
    constexpr std::size_t Nx_seed = 16, Ny_seed = 16;
    std::vector<double> xvec = linear_uniform(0.0, 1.0, Nx_seed);
    std::vector<double> yvec = log_uniform(ymin, ymax, Ny_seed);

    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (!helm) return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};
    double rho_crit_local = 0.0;
    try {
        helm->ensure_caloric_superancillaries();
        rho_crit_local = static_cast<double>(helm->rhomolar_critical());
    } catch (...) {
        return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};
    }
    if (rho_crit_local <= 0.0) return {linear_uniform(0.0, 1.0, 200), log_uniform(ymin, ymax, 200)};

    std::shared_ptr<CoolProp::AbstractState> probe_AS;
    try {
        const auto names = helm->fluid_names();
        if (names.empty()) return {std::move(xvec), std::move(yvec)};
        probe_AS.reset(CoolProp::AbstractState::factory("HEOS", names[0]));
    } catch (...) {
        return {std::move(xvec), std::move(yvec)};
    }

    const double T_min_eos = std::max(static_cast<double>(probe_AS->Ttriple()), static_cast<double>(probe_AS->Tmin()));
    const double T_max_ext = static_cast<double>(probe_AS->Tmax()) * 1.499;
    const bool is_vapor = (region_int == 1);
    const bool is_liquid = (region_int == 0);
    const bool is_super = (region_int == 2);
    if (!is_vapor && !is_liquid && !is_super) return {std::move(xvec), std::move(yvec)};

    // Region-specific T bounds & their dp derivatives at a given p.
    // SUPER uses T_lo = T_min_eos and T_hi = T_max_ext (both constants).
    // Production's SUPER T_lo follows the melting curve T_melt(p); using
    // a constant approximation here ignores that curvature but is exact
    // for Argon/CO2 (and any shallow-melting fluid) and accepts a small
    // build-vs-production cell-layout mismatch near the cold isotherm
    // for steep-melting fluids (water at very high p).  The near-critical
    // pillar — the actual accuracy bottleneck — is far from the cold
    // boundary so the approximation doesn't affect worst-case states.
    auto bounds_at_p = [&probe_AS, T_min_eos, T_max_ext, is_vapor, is_liquid](double p, double& T_lo, double& dT_lo_dp, double& T_hi,
                                                                              double& dT_hi_dp) -> bool {
        try {
            if (is_vapor) {
                probe_AS->update(PQ_INPUTS, p, 1.0);
                T_lo = static_cast<double>(probe_AS->T());
                dT_lo_dp = probe_AS->first_saturation_deriv(iT, iP);
                T_hi = T_max_ext;
                dT_hi_dp = 0.0;
            } else if (is_liquid) {
                T_lo = T_min_eos;
                dT_lo_dp = 0.0;
                probe_AS->update(PQ_INPUTS, p, 0.0);
                T_hi = static_cast<double>(probe_AS->T());
                dT_hi_dp = probe_AS->first_saturation_deriv(iT, iP);
            } else {  // SUPER
                T_lo = T_min_eos;
                dT_lo_dp = 0.0;
                T_hi = T_max_ext;
                dT_hi_dp = 0.0;
            }
            return std::isfinite(T_lo) && std::isfinite(T_hi) && std::isfinite(dT_lo_dp) && std::isfinite(dT_hi_dp);
        } catch (...) {
            return false;
        }
    };

    struct Corner
    {
        double f;     // ρ
        double f_T;   // ∂ρ/∂T|p
        double f_p;   // ∂ρ/∂p|T
        double f_TT;  // ∂²ρ/∂T²
        double f_Tp;  // ∂²ρ/∂T∂p
        double T_lo;
        double T_hi;
        double dT_lo_dp;
        double dT_hi_dp;
        double p;
        bool ok;
    };
    auto compute_corner = [&probe_AS, &bounds_at_p](double xnorm, double p) -> Corner {
        Corner c{};
        c.p = p;
        if (!bounds_at_p(p, c.T_lo, c.dT_lo_dp, c.T_hi, c.dT_hi_dp)) {
            c.ok = false;
            return c;
        }
        const double T = c.T_lo + xnorm * (c.T_hi - c.T_lo);
        try {
            probe_AS->update(PT_INPUTS, p, T);
            c.f = static_cast<double>(probe_AS->rhomolar());
            c.f_T = probe_AS->first_partial_deriv(iDmolar, iT, iP);
            c.f_p = probe_AS->first_partial_deriv(iDmolar, iP, iT);
            c.f_TT = probe_AS->second_partial_deriv(iDmolar, iT, iP, iT, iP);
            c.f_Tp = probe_AS->second_partial_deriv(iDmolar, iT, iP, iP, iT);
            c.ok = std::isfinite(c.f) && std::isfinite(c.f_T) && std::isfinite(c.f_p) && std::isfinite(c.f_TT) && std::isfinite(c.f_Tp);
        } catch (...) {
            c.ok = false;
        }
        return c;
    };
    auto rho_truth = [&probe_AS, &bounds_at_p](double xnorm, double p) -> double {
        double T_lo = 0.0, dT_lo_dp = 0.0, T_hi = 0.0, dT_hi_dp = 0.0;
        if (!bounds_at_p(p, T_lo, dT_lo_dp, T_hi, dT_hi_dp)) return std::numeric_limits<double>::quiet_NaN();
        const double T = T_lo + xnorm * (T_hi - T_lo);
        try {
            probe_AS->update(PT_INPUTS, p, T);
            return static_cast<double>(probe_AS->rhomolar());
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    };

    const SBTLTransform tform = sbtl_transform_for_property(region_int, iDmolar);
    const double y_ref = sbtl_transform_y_ref(region_int, iDmolar, rho_crit_local);
    const double cusp_min_d =
      (tform == SBTLTransform::CUSP_DIST_ABOVE_REF || tform == SBTLTransform::CUSP_DIST_BELOW_REF) ? 1e-6 * std::abs(y_ref) : 0.0;

    auto cell_relerr = [&rho_truth, rho_crit_local, tform, y_ref, cusp_min_d](const Corner& c00, const Corner& c10, const Corner& c01,
                                                                              const Corner& c11, double xn_a, double xn_b, double p_a,
                                                                              double p_b) -> double {
        const double Delta_xi = xn_b - xn_a;
        const double dlog_p = std::log(p_b) - std::log(p_a);
        const std::array<const Corner*, 4> cs = {&c00, &c10, &c01, &c11};
        const std::array<double, 4> xn = {xn_a, xn_b, xn_a, xn_b};
        std::array<double, 4> f{}, fxi{}, feta{}, fxieta{};
        for (int k = 0; k < 4; ++k) {
            const Corner& c = *cs[k];
            if (!c.ok) return std::numeric_limits<double>::infinity();
            if (tform == SBTLTransform::CUSP_DIST_ABOVE_REF && (c.f - y_ref) <= cusp_min_d) return 0.0;
            if (tform == SBTLTransform::CUSP_DIST_BELOW_REF && (y_ref - c.f) <= cusp_min_d) return 0.0;
            const double Delta_T_p = c.T_hi - c.T_lo;
            const double dDelta_T_dlogp = c.p * (c.dT_hi_dp - c.dT_lo_dp);
            const double dT_dlogp_at_xn = c.p * (c.dT_lo_dp + xn[k] * (c.dT_hi_dp - c.dT_lo_dp));
            double df_dxnorm = c.f_T * Delta_T_p;
            double df_dlogp = c.p * c.f_p + c.f_T * dT_dlogp_at_xn;
            double d2f = Delta_T_p * (c.f_TT * dT_dlogp_at_xn + c.p * c.f_Tp) + c.f_T * dDelta_T_dlogp;
            double fval = c.f;
            if (tform != SBTLTransform::IDENTITY) {
                if (!(fval > 0.0)) return std::numeric_limits<double>::infinity();
                const double gp = sbtl_transform_deriv(tform, fval, y_ref);
                const double gpp = sbtl_transform_second_deriv(tform, fval, y_ref);
                const double new_d2 = gpp * df_dxnorm * df_dlogp + gp * d2f;
                df_dxnorm *= gp;
                df_dlogp *= gp;
                d2f = new_d2;
                fval = sbtl_transform_apply(tform, fval, y_ref);
            }
            f[k] = fval;
            fxi[k] = df_dxnorm * Delta_xi;
            feta[k] = df_dlogp * dlog_p;
            fxieta[k] = d2f * Delta_xi * dlog_p;
        }
        const auto alpha = hermite_bicubic_polynomial_coeffs(f[0], f[1], f[2], f[3], fxi[0], fxi[1], fxi[2], fxi[3], feta[0], feta[1], feta[2],
                                                             feta[3], fxieta[0], fxieta[1], fxieta[2], fxieta[3]);
        constexpr int N_test = 3;
        constexpr std::array<double, N_test> ts = {0.25, 0.5, 0.75};
        double worst_relerr = 0.0;
        for (int ti = 0; ti < N_test; ++ti) {
            for (int tj = 0; tj < N_test; ++tj) {
                const double xi = ts[ti], eta = ts[tj];
                const double xi2 = xi * xi, xi3 = xi2 * xi;
                const double eta2 = eta * eta, eta3 = eta2 * eta;
                const double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2 + alpha[3] * eta3;
                const double B1 = alpha[4] + alpha[5] * eta + alpha[6] * eta2 + alpha[7] * eta3;
                const double B2 = alpha[8] + alpha[9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
                const double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
                const double F_test = B0 + B1 * xi + B2 * xi2 + B3 * xi3;
                const double r_pred = sbtl_transform_inverse(tform, F_test, y_ref);
                const double xn_test = xn_a + xi * (xn_b - xn_a);
                const double p_test = std::exp(std::log(p_a) + eta * (std::log(p_b) - std::log(p_a)));
                const double r_heos = rho_truth(xn_test, p_test);
                if (!std::isfinite(r_pred) || !std::isfinite(r_heos)) return std::numeric_limits<double>::infinity();
                const double denom = std::max(std::abs(r_heos), rho_crit_local * 0.01);
                const double relerr = std::abs(r_pred - r_heos) / denom;
                if (relerr > worst_relerr) worst_relerr = relerr;
            }
        }
        return worst_relerr;
    };

    auto recompute_corners = [&compute_corner](const std::vector<double>& xv, const std::vector<double>& yv) -> std::vector<std::vector<Corner>> {
        std::vector<std::vector<Corner>> m(xv.size(), std::vector<Corner>(yv.size()));
        for (std::size_t i = 0; i < xv.size(); ++i) {
            for (std::size_t j = 0; j < yv.size(); ++j) {
                m[i][j] = compute_corner(xv[i], yv[j]);
            }
        }
        return m;
    };

    constexpr int max_iter = 14;
    for (int iter = 0; iter < max_iter; ++iter) {
        auto corners = recompute_corners(xvec, yvec);
        bool any_split = false;
        if (yvec.size() < Ny_max) {
            std::vector<double> new_yvec;
            new_yvec.reserve(yvec.size() * 2);
            for (std::size_t j = 0; j + 1 < yvec.size(); ++j) {
                new_yvec.push_back(yvec[j]);
                bool row_bad = false;
                for (std::size_t i = 0; i + 1 < xvec.size(); ++i) {
                    if (cell_relerr(corners[i][j], corners[i + 1][j], corners[i][j + 1], corners[i + 1][j + 1], xvec[i], xvec[i + 1], yvec[j],
                                    yvec[j + 1])
                        > tol_relerr) {
                        row_bad = true;
                        break;
                    }
                }
                if (row_bad && new_yvec.size() + (yvec.size() - j - 1) < Ny_max) {
                    new_yvec.push_back(std::exp(0.5 * (std::log(yvec[j]) + std::log(yvec[j + 1]))));
                    any_split = true;
                }
            }
            new_yvec.push_back(yvec.back());
            yvec = std::move(new_yvec);
        }
        if (xvec.size() < Nx_max) {
            corners = recompute_corners(xvec, yvec);
            std::vector<double> new_xvec;
            new_xvec.reserve(xvec.size() * 2);
            for (std::size_t i = 0; i + 1 < xvec.size(); ++i) {
                new_xvec.push_back(xvec[i]);
                bool col_bad = false;
                for (std::size_t j = 0; j + 1 < yvec.size(); ++j) {
                    if (cell_relerr(corners[i][j], corners[i + 1][j], corners[i][j + 1], corners[i + 1][j + 1], xvec[i], xvec[i + 1], yvec[j],
                                    yvec[j + 1])
                        > tol_relerr) {
                        col_bad = true;
                        break;
                    }
                }
                if (col_bad && new_xvec.size() + (xvec.size() - i - 1) < Nx_max) {
                    new_xvec.push_back(0.5 * (xvec[i] + xvec[i + 1]));
                    any_split = true;
                }
            }
            new_xvec.push_back(xvec.back());
            xvec = std::move(new_xvec);
        }
        if (!any_split) break;
    }

    return {std::move(xvec), std::move(yvec)};
}

std::vector<double> SBTLBackend::build_adaptive_yvec(double ymin, double ymax, std::size_t Ny_target, std::size_t Ny_max,
                                                     const std::function<double(double)>& /*prop_at_p*/, double /*tol_rel*/) const {
    // Adaptive subdivision driven by the H-superancillary's piece boundaries.
    //
    // Seed the grid with the superancillary's interval edges that fall
    // inside (ymin, ymax) — these are placed by the superancillary's own
    // dyadic-splitting build right where T_sat (and therefore h_sat /
    // ρ_sat) changes shape, so they're already-fluid-tuned anchors.
    //
    // Refinement criterion: ρ_sat (the property whose sqrt-cusp drives
    // SBTL's cell-level error near p_crit) linear-interpolation error
    // at the log-p midpoint of each interval, normalised by ρ_crit.
    // Bisect intervals where this exceeds tol_rho.  ρ_sat itself is what
    // SBTL has to represent — using it directly as the refinement
    // sensitivity captures the cusp without proxy artifacts.
    //
    // For subcritical tables ymax = 0.999·p_crit, so ρ_sat,V is well-
    // defined throughout (ymin, ymax) via the H-superancillary.  For
    // SUPER (above p_crit) ρ_sat doesn't apply — caller skips the
    // adaptive path for that region.
    //
    // Falls back to log-uniform on any access failure (mixtures, fluids
    // lacking the superancillary, etc.).
    const std::size_t Ny_floor = std::max(Ny_target, static_cast<std::size_t>(40));
    auto log_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        const double la = std::log(a), lb = std::log(b);
        for (std::size_t k = 0; k < n; ++k) {
            const double frac = static_cast<double>(k) / static_cast<double>(n - 1);
            v.push_back(std::exp(la + frac * (lb - la)));
        }
        return v;
    };
    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (!helm) return log_uniform(ymin, ymax, Ny_floor);

    // Step 1: seed anchors as a SMALL log-uniform grid (no fluid-specific
    // breakpoints).  Using the H-superancillary piece edges as hard anchors
    // was tried — it produced discrete bands of high-error cells at each
    // anchor location because the bisection-density-across-anchor changed
    // abruptly, creating the same kind of cell-spacing discontinuity that
    // the earlier two-zone scheme had at 0.5·p_crit (just at multiple
    // locations instead of one).  Starting from a uniform log grid and
    // letting the bisection drive density wherever ρ_sat curvature
    // demands it produces a smooth cell-spacing profile.
    constexpr std::size_t Ny_seed = 16;
    std::vector<double> anchors = log_uniform(ymin, ymax, Ny_seed);
    double rho_crit = 0.0;
    try {
        helm->ensure_caloric_superancillaries();  // probe AS path below needs it
        rho_crit = static_cast<double>(helm->rhomolar_critical());
    } catch (...) {
        return log_uniform(ymin, ymax, Ny_floor);
    }
    if (rho_crit <= 0.0) return log_uniform(ymin, ymax, Ny_floor);

    // Step 2: bisection-to-tolerance on ρ_sat,V linear-interp error at log midpoint.
    // ρ_sat,V at the supercritical-shoulder edge (p just below p_crit) varies as
    //   ρ_sat,V ~ ρ_crit − C (p_crit − p)^(1/2)
    // so an interval whose midpoint shows large lerp-error in ρ_sat is one
    // the bicubic can't bridge with a single cell — bisect it.  Tolerance
    // tol_rho is in units of ρ_crit; 1 % of ρ_crit at each refinement step
    // converges in ~6 levels of bisection near p_crit, ending with cells
    // whose log-p width ~ 0.0005·p_crit there.
    // Use a SEPARATE AS instance for the ρ_sat probes so the SBTL build's
    // AS->update sequence isn't disrupted (the cell-fill self-checks at
    // 1e-12 relative tolerance are sensitive to HEOS internal-cache state
    // carrying over from PQ_INPUTS calls into subsequent HmolarP_INPUTS
    // calls).  The shared AS would otherwise alternate between PQ + PH
    // state and pick up tiny initial-guess differences.
    std::shared_ptr<CoolProp::AbstractState> probe_AS;
    try {
        const std::vector<std::string> names = helm->fluid_names();
        if (names.empty()) return log_uniform(ymin, ymax, Ny_floor);
        probe_AS.reset(CoolProp::AbstractState::factory("HEOS", names[0]));
    } catch (...) {
        return log_uniform(ymin, ymax, Ny_floor);
    }
    // Cell-midpoint ρ-error metric.  For each candidate interval (p_a, p_b),
    // probe ρ at the same η coordinate at three pressures (p_a, p_mid, p_b)
    // across η ∈ {0.02, 0.25, 0.5, 0.75, 0.98} (avoiding the sat boundaries
    // exactly so the HmolarP_INPUTS probe doesn't land in the dome's
    // ambiguity).  Compare HEOS(η, p_mid) to the linear interpolation of
    // ρ between (η, p_a) and (η, p_b) — this is the lerp-residual the
    // cell would carry, which upper-bounds the actual Hermite-bicubic
    // residual (cubic is more accurate than linear, given the same
    // corners).  Take the max over η probes as the interval's "needs
    // refinement?" score.
    //
    // Compared to the earlier ρ_sat,V-at-η=0 lerp criterion, this measures
    // bulk-cell error too, not just the dome boundary.  Bands of bad
    // states away from the dome (~Argon at 0.5·p_crit, Water at multiple
    // bisection-level boundaries) should drop because bisection now sees
    // those bulk cells when they're poorly approximated.
    //
    // η-coordinate is region-dependent.  For LIQUID: η spans [0, 1] from
    // h_min(p)=h(T_min,p) to h_sat,L(p).  For VAPOR: η spans from
    // h_sat,V(p) to h_hi(p)=h(T_max_ext,p).
    const double T_min_eos = std::max(static_cast<double>(probe_AS->Ttriple()), static_cast<double>(probe_AS->Tmin()));
    const double T_max_ext = static_cast<double>(probe_AS->Tmax()) * 1.499;
    auto h_bounds_at_p_LIQ = [&probe_AS, T_min_eos](double p, double& h_lo, double& h_hi) -> bool {
        try {
            probe_AS->update(PT_INPUTS, p, T_min_eos);
            h_lo = static_cast<double>(probe_AS->hmolar());
            probe_AS->update(PQ_INPUTS, p, 0.0);
            h_hi = static_cast<double>(probe_AS->hmolar());
            return std::isfinite(h_lo) && std::isfinite(h_hi);
        } catch (...) {
            return false;
        }
    };
    auto h_bounds_at_p_VAP = [&probe_AS, T_max_ext](double p, double& h_lo, double& h_hi) -> bool {
        try {
            probe_AS->update(PQ_INPUTS, p, 1.0);
            h_lo = static_cast<double>(probe_AS->hmolar());
            probe_AS->update(PT_INPUTS, p, T_max_ext);
            h_hi = static_cast<double>(probe_AS->hmolar());
            return std::isfinite(h_lo) && std::isfinite(h_hi);
        } catch (...) {
            return false;
        }
    };
    // Determine region by probing: if VAPOR-style h bounds work and produce
    // a positive span, treat as VAPOR; same for LIQUID; if both work pick
    // by checking that h_hi - h_lo is wider in one (the more-relevant
    // region for our table).  Heuristic: probe at p_mid in the middle of
    // the yvec.
    bool is_vapor = false, is_liquid = false;
    {
        const double p_probe = std::exp(0.5 * (std::log(ymin) + std::log(ymax)));
        double h_lo, h_hi;
        is_vapor = h_bounds_at_p_VAP(p_probe, h_lo, h_hi) && (h_hi > h_lo);
        is_liquid = h_bounds_at_p_LIQ(p_probe, h_lo, h_hi) && (h_hi > h_lo);
    }
    auto h_bounds_at_p = [is_vapor, is_liquid, &h_bounds_at_p_LIQ, &h_bounds_at_p_VAP](double p, double& h_lo, double& h_hi) -> bool {
        // Heuristic resolved at p_mid: prefer VAPOR if it works (most
        // problematic for the sqrt-cusp), fall back to LIQUID.
        if (is_vapor) return h_bounds_at_p_VAP(p, h_lo, h_hi);
        if (is_liquid) return h_bounds_at_p_LIQ(p, h_lo, h_hi);
        return false;
    };

    auto rho_at_eta_p = [&probe_AS, &h_bounds_at_p](double eta, double p) -> double {
        double h_lo = 0.0, h_hi = 0.0;
        if (!h_bounds_at_p(p, h_lo, h_hi)) return std::numeric_limits<double>::quiet_NaN();
        const double h = h_lo + eta * (h_hi - h_lo);
        try {
            probe_AS->update(HmolarP_INPUTS, h, p);
            return static_cast<double>(probe_AS->rhomolar());
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    };

    // True cell-vs-HEOS metric.  For each candidate interval (p_a, p_b)
    // and each η probe, build the cubic-Hermite-in-log-p that the actual
    // SBTL cell would use along that constant-η slice — using HEOS-supplied
    // values AND first derivatives at the corners — and compare its value
    // at the log-p midpoint to HEOS truth.  This is the actual cell-fit
    // residual the user will see, not the lerp upper bound.
    //
    // The Hermite cubic basis at t=0.5 evaluates to:
    //   h00(0.5) = h01(0.5) = 0.5
    //   h10(0.5) = +0.125,  h11(0.5) = −0.125
    // so the midpoint value is
    //   ρ_mid_hermite = 0.5·(ρ_a + ρ_b) + 0.125·Δlogp·(ρ'_a − ρ'_b)
    // where ρ'_x = dρ/d(log p) at p_x, computed by central finite difference.
    //
    // ANALYTIC total derivative dρ/d(log p) at constant η — what the cell
    // actually stores along its constant-η log-p slice.  Two sources of
    // p-dependence at fixed η:
    //   (1) ρ varies with p at fixed h:         ∂ρ/∂p|_h
    //   (2) h itself shifts with p at fixed η:  dh/dp|_η = dh_lo/dp + η·(dh_hi/dp − dh_lo/dp)
    // Chain rule:
    //   dρ/dp|_η = ∂ρ/∂p|_h + ∂ρ/∂h|_p · dh/dp|_η
    //   dρ/d(log p)|_η = p · dρ/dp|_η
    //
    // For each region the h_lo and h_hi derivatives come from different
    // HEOS paths:
    //   VAPOR : h_lo = h_sat,V(p) via PQ + first_saturation_deriv(H, P);
    //           h_hi = h(T_max_ext, p) at fixed T via first_partial_deriv(H, P, T)
    //   LIQUID: h_lo = h(T_min, p);  h_hi = h_sat,L(p) via PQ + first_saturation_deriv
    auto bounds_and_derivs_at_p = [&probe_AS, T_min_eos, T_max_ext, is_vapor, is_liquid](double p, double& h_lo, double& dh_lo_dp, double& h_hi,
                                                                                         double& dh_hi_dp) -> bool {
        try {
            if (is_vapor) {
                probe_AS->update(PQ_INPUTS, p, 1.0);
                h_lo = static_cast<double>(probe_AS->hmolar());
                dh_lo_dp = probe_AS->first_saturation_deriv(iHmolar, iP);
                probe_AS->update(PT_INPUTS, p, T_max_ext);
                h_hi = static_cast<double>(probe_AS->hmolar());
                dh_hi_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
            } else if (is_liquid) {
                probe_AS->update(PT_INPUTS, p, T_min_eos);
                h_lo = static_cast<double>(probe_AS->hmolar());
                dh_lo_dp = probe_AS->first_partial_deriv(iHmolar, iP, iT);
                probe_AS->update(PQ_INPUTS, p, 0.0);
                h_hi = static_cast<double>(probe_AS->hmolar());
                dh_hi_dp = probe_AS->first_saturation_deriv(iHmolar, iP);
            } else {
                return false;
            }
            return std::isfinite(h_lo) && std::isfinite(h_hi) && std::isfinite(dh_lo_dp) && std::isfinite(dh_hi_dp);
        } catch (...) {
            return false;
        }
    };

    // Full corner data at (xnorm, p): ρ and the four HEOS partials needed
    // to chain-rule into (xnorm, log p) cell-coords and apply the same
    // output transform that production uses.  Matches build_normph_hermite_alphas.
    auto corner_full_data = [&probe_AS, &bounds_and_derivs_at_p](double xnorm, double p, double& f, double& f_h, double& f_p, double& f_hh,
                                                                 double& f_hp, double& dh_lo_dp, double& dh_hi_dp) -> bool {
        double h_lo = 0.0, h_hi = 0.0;
        if (!bounds_and_derivs_at_p(p, h_lo, dh_lo_dp, h_hi, dh_hi_dp)) return false;
        const double h = h_lo + xnorm * (h_hi - h_lo);
        try {
            probe_AS->update(HmolarP_INPUTS, h, p);
            f = static_cast<double>(probe_AS->rhomolar());
            f_h = probe_AS->first_partial_deriv(iDmolar, iHmolar, iP);
            f_p = probe_AS->first_partial_deriv(iDmolar, iP, iHmolar);
            f_hh = probe_AS->second_partial_deriv(iDmolar, iHmolar, iP, iHmolar, iP);
            f_hp = probe_AS->second_partial_deriv(iDmolar, iHmolar, iP, iP, iHmolar);
        } catch (...) {
            return false;
        }
        return std::isfinite(f) && std::isfinite(f_h) && std::isfinite(f_p) && std::isfinite(f_hh) && std::isfinite(f_hp);
    };

    // Which output transform will the production cell apply?  Must match
    // build_normph_hermite_alphas, otherwise the metric and the stored
    // polynomial disagree.
    const int region_int = is_vapor ? 1 : (is_liquid ? 0 : -1);
    const SBTLTransform tform = (region_int >= 0) ? sbtl_transform_for_property(region_int, iDmolar) : SBTLTransform::IDENTITY;
    const double y_ref = (region_int >= 0) ? sbtl_transform_y_ref(region_int, iDmolar, rho_crit) : 0.0;

    // True 2D-bicubic-vs-HEOS metric.  For a candidate yvec interval
    // (p_a, p_b) and each xnorm test cell, build the actual production-
    // style Hermite bicubic (corner ρ + chain-ruled first derivs +
    // cross deriv, all carried through the active output transform),
    // evaluate at (ξ=0.5, η=0.5), invert the transform, and compare to
    // HEOS truth at the same physical (xnorm_mid, p_mid).
    //
    // Why: the prior 1D-Hermite-in-ρ proxy underestimates the residual
    // the user sees, because it (a) ignores the LOG/CUSP transform the
    // bicubic actually stores in, and (b) misses the xi-direction
    // interpolation component of the 2D fit.  Empirically: production
    // cells passed the proxy at <1e-7 residual but exhibited >1e-3 user-
    // visible error.  Replacing the proxy with the actual bicubic
    // evaluation removes that blind spot.
    // Tile the xnorm direction with TEST cells at production xnorm spacing
    // (Δxi ≈ 1/Nx_test).  Pre-fetch corner data at Nx_test+1 xnorm anchors
    // × {p_a, p_b}, then for each adjacent pair build the bicubic and
    // check its midpoint.  Tiling (no gaps between test cells) ensures
    // every production cell is represented; the previous strided-probe
    // version missed cells falling between probes.
    constexpr std::size_t Nx_test = 200;  // matches production Nx baseline
    constexpr double Delta_xi = 1.0 / static_cast<double>(Nx_test);
    auto cell_max_relerr = [&corner_full_data, &rho_at_eta_p, &bounds_and_derivs_at_p, rho_crit, tform, y_ref, Delta_xi](double p_a,
                                                                                                                         double p_b) -> double {
        const double dlog_p = std::log(p_b) - std::log(p_a);
        const double p_mid = std::exp(0.5 * (std::log(p_a) + std::log(p_b)));
        struct Corner
        {
            double f, f_h, f_p, f_hh, f_hp, dh_lo_dp, dh_hi_dp;
            bool ok;
        };
        // Stay in [0.01, 0.99] to avoid HmolarP_INPUTS landing exactly on
        // the dome at xnorm=0 (sat) or hitting the T_max_ext extrapolation
        // cliff at xnorm=1.
        const double xn_lo_global = 0.01, xn_hi_global = 0.99;
        const std::size_t n_anchors = Nx_test + 1;
        std::vector<Corner> col_a(n_anchors), col_b(n_anchors);
        std::vector<double> xn_grid(n_anchors);
        for (std::size_t k = 0; k < n_anchors; ++k) {
            xn_grid[k] = xn_lo_global + (xn_hi_global - xn_lo_global) * static_cast<double>(k) / static_cast<double>(Nx_test);
            col_a[k].ok = corner_full_data(xn_grid[k], p_a, col_a[k].f, col_a[k].f_h, col_a[k].f_p, col_a[k].f_hh, col_a[k].f_hp, col_a[k].dh_lo_dp,
                                           col_a[k].dh_hi_dp);
            col_b[k].ok = corner_full_data(xn_grid[k], p_b, col_b[k].f, col_b[k].f_h, col_b[k].f_p, col_b[k].f_hh, col_b[k].f_hp, col_b[k].dh_lo_dp,
                                           col_b[k].dh_hi_dp);
        }
        double max_relerr = 0.0;
        for (std::size_t k = 0; k + 1 < n_anchors; ++k) {
            const double xn_a = xn_grid[k];
            const double xn_b = xn_grid[k + 1];
            const double xn_mid = 0.5 * (xn_a + xn_b);
            if (!col_a[k].ok || !col_a[k + 1].ok || !col_b[k].ok || !col_b[k + 1].ok) continue;
            const std::array<const Corner*, 4> cds = {&col_a[k], &col_a[k + 1], &col_b[k], &col_b[k + 1]};
            const std::array<double, 4> xn = {xn_a, xn_b, xn_a, xn_b};
            const std::array<double, 4> p_at = {p_a, p_a, p_b, p_b};
            std::array<double, 4> f{}, fxi{}, feta{}, fxieta{};
            bool ok = true;
            for (int c = 0; c < 4; ++c) {
                if (!cds[c]->ok) {
                    ok = false;
                    break;
                }
                const double Delta_h = 1.0;  // we scale by Δh below via f_h·(h_hi-h_lo); but we need it for chain rule
                // Recompute Δh and dh-derivatives at this corner's p (cheaper to refetch from corner):
                const double p_k = p_at[c];
                const double h_lo_p_unused = 0.0;
                (void)h_lo_p_unused;
                // Reconstruct Δh and dDelta_h_dp from the dh_lo_dp / dh_hi_dp we cached:
                // Need h_lo, h_hi for this p — easiest: get from bounds_and_derivs_at_p again,
                // but that's a duplicate call.  Instead, since rho_at_eta_p was already called
                // with the same h_lo/h_hi internally, we don't have direct access here.  Re-query.
                double h_lo_loc = 0.0, h_hi_loc = 0.0, dh_lo_dp_unused = 0.0, dh_hi_dp_unused = 0.0;
                if (!bounds_and_derivs_at_p(p_k, h_lo_loc, dh_lo_dp_unused, h_hi_loc, dh_hi_dp_unused)) {
                    ok = false;
                    break;
                }
                const double Delta_h_p = h_hi_loc - h_lo_loc;
                const double dDelta_h_dlogp = p_k * (cds[c]->dh_hi_dp - cds[c]->dh_lo_dp);
                const double dh_dlogp_at_xn = p_k * (cds[c]->dh_lo_dp + xn[c] * (cds[c]->dh_hi_dp - cds[c]->dh_lo_dp));
                (void)Delta_h;
                // Partials in (xnorm-global, log p) coords:
                double df_dxnorm = cds[c]->f_h * Delta_h_p;
                double df_dlogp = p_k * cds[c]->f_p + cds[c]->f_h * dh_dlogp_at_xn;
                double d2f = Delta_h_p * (cds[c]->f_hh * dh_dlogp_at_xn + p_k * cds[c]->f_hp) + cds[c]->f_h * dDelta_h_dlogp;
                double fval = cds[c]->f;
                // Apply the active output transform — same chain rule as build_normph_hermite_alphas.
                if (tform != SBTLTransform::IDENTITY) {
                    if (!(fval > 0.0)) {
                        ok = false;
                        break;
                    }
                    const double gp = sbtl_transform_deriv(tform, fval, y_ref);
                    const double gpp = sbtl_transform_second_deriv(tform, fval, y_ref);
                    const double new_d2 = gpp * df_dxnorm * df_dlogp + gp * d2f;
                    df_dxnorm *= gp;
                    df_dlogp *= gp;
                    d2f = new_d2;
                    fval = sbtl_transform_apply(tform, fval, y_ref);
                }
                f[c] = fval;
                fxi[c] = df_dxnorm * Delta_xi;
                feta[c] = df_dlogp * dlog_p;
                fxieta[c] = d2f * Delta_xi * dlog_p;
            }
            if (!ok) continue;
            const auto alpha = hermite_bicubic_polynomial_coeffs(f[0], f[1], f[2], f[3], fxi[0], fxi[1], fxi[2], fxi[3], feta[0], feta[1], feta[2],
                                                                 feta[3], fxieta[0], fxieta[1], fxieta[2], fxieta[3]);
            // Bicubic eval at (ξ=0.5, η=0.5):
            const double xi = 0.5, eta = 0.5;
            const double xi2 = xi * xi, xi3 = xi2 * xi;
            const double eta2 = eta * eta, eta3 = eta2 * eta;
            const double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2 + alpha[3] * eta3;
            const double B1 = alpha[4] + alpha[5] * eta + alpha[6] * eta2 + alpha[7] * eta3;
            const double B2 = alpha[8] + alpha[9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
            const double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
            const double F_mid = B0 + B1 * xi + B2 * xi2 + B3 * xi3;
            const double r_pred = sbtl_transform_inverse(tform, F_mid, y_ref);
            const double r_heos_mid = rho_at_eta_p(xn_mid, p_mid);
            if (!std::isfinite(r_pred) || !std::isfinite(r_heos_mid)) continue;
            const double denom = std::max(std::abs(r_heos_mid), rho_crit * 0.01);
            const double relerr = std::abs(r_pred - r_heos_mid) / denom;
            if (relerr > max_relerr) max_relerr = relerr;
        }
        return max_relerr;
    };

    // Tolerance: 1e-4 relative — what G13-15 demands for ρ.  Bisection
    // refines until the Hermite cubic in log p matches HEOS at the cell
    // midpoint to this tolerance, for the worst η probe in the cell.
    const double tol_relerr = 1e-4;
    int iter = 0;
    constexpr int max_iter = 12;
    while (iter < max_iter && anchors.size() < Ny_max) {
        std::vector<double> next;
        next.reserve(anchors.size() * 2);
        bool progressing = false;
        for (std::size_t j = 0; j + 1 < anchors.size(); ++j) {
            const double a = anchors[j], b = anchors[j + 1];
            next.push_back(a);
            if (anchors.size() + (next.size() - j - 1) < Ny_max) {
                const double mid = std::exp(0.5 * (std::log(a) + std::log(b)));
                if (cell_max_relerr(a, b) > tol_relerr) {
                    next.push_back(mid);
                    progressing = true;
                }
            }
        }
        next.push_back(anchors.back());
        anchors = std::move(next);
        if (!progressing) break;
        ++iter;
    }

    // Step 3: pad with log-uniform fillers between sparse anchor pairs to bring
    // total count up to Ny_target.  The bisection above places anchors based
    // on local sensitivity; fillers distribute remaining rows proportional to
    // each interval's log-p span so smooth bulk regions still get adequate
    // resolution (without those fillers a fluid with few cusp features
    // could end up with only ~30 rows total).
    if (anchors.size() >= Ny_target) return anchors;
    const std::size_t deficit = Ny_target - anchors.size();
    const double total_log_span = std::log(ymax) - std::log(ymin);
    std::vector<double> yvec;
    yvec.reserve(Ny_target);
    for (std::size_t j = 0; j + 1 < anchors.size(); ++j) {
        const double a = anchors[j], b = anchors[j + 1];
        const double log_span = std::log(b) - std::log(a);
        const std::size_t n_fill = static_cast<std::size_t>(std::round(deficit * log_span / total_log_span));
        yvec.push_back(a);
        for (std::size_t i = 1; i <= n_fill; ++i) {
            const double frac = static_cast<double>(i) / static_cast<double>(n_fill + 1);
            yvec.push_back(std::exp(std::log(a) + frac * log_span));
        }
    }
    yvec.push_back(anchors.back());
    return yvec;
}

std::vector<double> SBTLBackend::build_adaptive_xvec(int region_int, std::size_t Nx_target, std::size_t Nx_max, double ymin, double ymax) const {
    // Adaptive bisection on η ∈ [0, 1] driven by Hermite-cubic-vs-HEOS
    // error at η midpoints.  Probes at the log-midpoint pressure
    // p_anchor as a representative state (single-p probe trades some
    // optimality for build cost).  Same mechanics as build_adaptive_yvec
    // along the orthogonal axis.
    //
    // Hermite cubic at η midpoint with corner derivatives dρ/dη|_p:
    //   dρ/dη|_p = (∂ρ/∂h)|_p · (h_hi - h_lo)
    //
    // The (h_hi - h_lo) factor is the cell's η-stretch at p_anchor — it
    // scales how rapidly ρ changes per unit η.  All-analytic from HEOS
    // partials; no FD.
    const std::size_t Nx_floor = std::max(Nx_target, static_cast<std::size_t>(40));
    auto linear_uniform = [](double a, double b, std::size_t n) {
        std::vector<double> v;
        v.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            const double frac = static_cast<double>(k) / static_cast<double>(n - 1);
            v.push_back(a + frac * (b - a));
        }
        return v;
    };

    auto helm = std::dynamic_pointer_cast<HelmholtzEOSMixtureBackend>(this->AS);
    if (!helm) return linear_uniform(0.0, 1.0, Nx_floor);

    double rho_crit = 0.0;
    std::shared_ptr<CoolProp::AbstractState> probe_AS;
    try {
        helm->ensure_caloric_superancillaries();
        rho_crit = static_cast<double>(helm->rhomolar_critical());
        const std::vector<std::string> names = helm->fluid_names();
        if (names.empty()) return linear_uniform(0.0, 1.0, Nx_floor);
        probe_AS.reset(CoolProp::AbstractState::factory("HEOS", names[0]));
    } catch (...) {
        return linear_uniform(0.0, 1.0, Nx_floor);
    }
    if (rho_crit <= 0.0) return linear_uniform(0.0, 1.0, Nx_floor);

    // For an η interval, probe at ALL p values across the yvec seed
    // (not a single p_anchor) — covers the full row of cells stacked
    // above each η interval, not a slice.  Per Ian: "the check should
    // be for all 16×16 cells, not one per row".  Pressure probes are
    // 16 log-uniform values in (ymin, ymax) — matches the yvec seed.
    constexpr std::size_t Np_probes = 16;
    std::vector<double> p_probes(Np_probes);
    for (std::size_t k = 0; k < Np_probes; ++k) {
        const double frac = static_cast<double>(k) / static_cast<double>(Np_probes - 1);
        p_probes[k] = std::exp(std::log(ymin) + frac * (std::log(ymax) - std::log(ymin)));
    }
    const double T_min_eos = std::max(static_cast<double>(probe_AS->Ttriple()), static_cast<double>(probe_AS->Tmin()));
    const double T_max_ext = static_cast<double>(probe_AS->Tmax()) * 1.499;
    const bool is_liquid = (region_int == 0);
    const bool is_vapor = (region_int == 1);

    // For each probe pressure, get (h_lo, h_hi) so we can map η→h.
    std::vector<double> h_lo_at_p(Np_probes), h_hi_at_p(Np_probes);
    std::vector<bool> p_probe_valid(Np_probes, false);
    for (std::size_t k = 0; k < Np_probes; ++k) {
        try {
            if (is_vapor) {
                probe_AS->update(PQ_INPUTS, p_probes[k], 1.0);
                h_lo_at_p[k] = static_cast<double>(probe_AS->hmolar());
                probe_AS->update(PT_INPUTS, p_probes[k], T_max_ext);
                h_hi_at_p[k] = static_cast<double>(probe_AS->hmolar());
            } else if (is_liquid) {
                probe_AS->update(PT_INPUTS, p_probes[k], T_min_eos);
                h_lo_at_p[k] = static_cast<double>(probe_AS->hmolar());
                probe_AS->update(PQ_INPUTS, p_probes[k], 0.0);
                h_hi_at_p[k] = static_cast<double>(probe_AS->hmolar());
            }
            p_probe_valid[k] = std::isfinite(h_lo_at_p[k]) && std::isfinite(h_hi_at_p[k]) && (h_hi_at_p[k] > h_lo_at_p[k]);
        } catch (...) {
            p_probe_valid[k] = false;
        }
    }

    auto rho_and_drho_deta = [&probe_AS](double eta, double p, double h_lo_p, double h_hi_p, double& rho, double& drho_deta) -> bool {
        const double delta_h_p = h_hi_p - h_lo_p;
        const double h = h_lo_p + eta * delta_h_p;
        try {
            probe_AS->update(HmolarP_INPUTS, h, p);
            rho = static_cast<double>(probe_AS->rhomolar());
            const double drho_dh = probe_AS->first_partial_deriv(iDmolar, iHmolar, iP);
            drho_deta = drho_dh * delta_h_p;
        } catch (...) {
            return false;
        }
        return std::isfinite(rho) && std::isfinite(drho_deta);
    };

    auto cell_eta_relerr = [&rho_and_drho_deta, &p_probes, &h_lo_at_p, &h_hi_at_p, &p_probe_valid, rho_crit, Np_probes](double eta_a,
                                                                                                                        double eta_b) -> double {
        const double eta_mid = 0.5 * (eta_a + eta_b);
        const double d_eta = eta_b - eta_a;
        double max_relerr = 0.0;
        // Take the worst-case relerr across all pressure probes — covers
        // every cell in the column above this η interval.
        for (std::size_t k = 0; k < Np_probes; ++k) {
            if (!p_probe_valid[k]) continue;
            const double p_k = p_probes[k], h_lo_k = h_lo_at_p[k], h_hi_k = h_hi_at_p[k];
            double r_a = 0.0, dr_a = 0.0, r_b = 0.0, dr_b = 0.0, r_mid_heos = 0.0, dr_mid_unused = 0.0;
            if (!rho_and_drho_deta(eta_a, p_k, h_lo_k, h_hi_k, r_a, dr_a)) continue;
            if (!rho_and_drho_deta(eta_b, p_k, h_lo_k, h_hi_k, r_b, dr_b)) continue;
            if (!rho_and_drho_deta(eta_mid, p_k, h_lo_k, h_hi_k, r_mid_heos, dr_mid_unused)) continue;
            const double r_mid_hermite = 0.5 * (r_a + r_b) + 0.125 * d_eta * (dr_a - dr_b);
            const double denom = std::max(std::abs(r_mid_heos), rho_crit * 0.01);
            const double relerr = std::abs(r_mid_hermite - r_mid_heos) / denom;
            if (relerr > max_relerr) max_relerr = relerr;
        }
        return max_relerr;
    };

    // Seed: uniform 16-point η grid in [0, 1].
    constexpr std::size_t Nx_seed = 16;
    std::vector<double> anchors = linear_uniform(0.0, 1.0, Nx_seed);

    const double tol_relerr = 1e-4;
    int iter = 0;
    constexpr int max_iter = 12;
    while (iter < max_iter && anchors.size() < Nx_max) {
        std::vector<double> next;
        next.reserve(anchors.size() * 2);
        bool progressing = false;
        for (std::size_t j = 0; j + 1 < anchors.size(); ++j) {
            const double a = anchors[j], b = anchors[j + 1];
            next.push_back(a);
            if (anchors.size() + (next.size() - j - 1) < Nx_max) {
                if (cell_eta_relerr(a, b) > tol_relerr) {
                    next.push_back(0.5 * (a + b));
                    progressing = true;
                }
            }
        }
        next.push_back(anchors.back());
        anchors = std::move(next);
        if (!progressing) break;
        ++iter;
    }

    // Pad with linear-uniform fillers up to Nx_target so smooth-η bulk
    // intervals don't end up with only the seed-grid density.
    if (anchors.size() >= Nx_target) return anchors;
    const std::size_t deficit = Nx_target - anchors.size();
    std::vector<double> xvec;
    xvec.reserve(Nx_target);
    for (std::size_t j = 0; j + 1 < anchors.size(); ++j) {
        const double a = anchors[j], b = anchors[j + 1];
        const double span = b - a;
        const std::size_t n_fill = static_cast<std::size_t>(std::round(deficit * span));
        xvec.push_back(a);
        for (std::size_t i = 1; i <= n_fill; ++i) {
            const double frac = static_cast<double>(i) / static_cast<double>(n_fill + 1);
            xvec.push_back(a + frac * span);
        }
    }
    xvec.push_back(anchors.back());
    return xvec;
}

void SBTLBackend::build_normph_table(NormalizedPHTable& table) {
    if (!this->AS) throw ValueError("build_normph_table: AS is not set");
    table.AS = this->AS;
    table.set_limits();

    // Determine final Ny BEFORE the single resize() call.  The base-class
    // resize(Nx, Ny) only grows the outer dimension; inner vectors keep
    // their old Ny size when Nx is unchanged.  Calling resize twice with
    // different Ny therefore corrupts the LIST_OF_MATRICES storage.
    // Compute the adaptive yvec first, set table.Ny to its size, then
    // resize exactly once.
    std::vector<double> adaptive_yvec;
    if (table.region() != NormalizedPHTable::SUPER) {
        const short Q = (table.region() == NormalizedPHTable::LIQUID) ? 0 : 1;
        auto h_sat_at_p = [this, Q](double p) -> double {
            double out = std::numeric_limits<double>::quiet_NaN();
            if (this->try_h_superanc(p, Q, out)) return out;
            return (Q == 0) ? this->saturation_hmolar_liquid(p) : this->saturation_hmolar_vapor(p);
        };
        // Subcritical regions get an adaptive yvec seeded from the
        // H-superancillary's piece breakpoints (already concentrated near
        // p_crit by the superancillary's own dyadic-splitting build),
        // then bisected where the seed leaves rows that miss accuracy
        // spec at the cell midpoint.  SUPER stays log-uniform — its
        // property surface is smooth.  Cap Ny growth at 2 × baseline so
        // storage stays bounded.
        adaptive_yvec = build_adaptive_yvec(table.ymin, table.ymax, table.Ny, /*Ny_max=*/10 * table.Ny, h_sat_at_p,
                                            /*tol_rel=*/0.005);
        if (adaptive_yvec.size() >= 4) {
            table.Ny = adaptive_yvec.size();
        }
    }

    // Joint xvec+yvec adaptive grid for LIQUID/VAPOR — replaces the
    // separate single-axis yvec / xvec bisection above and below.  The
    // joint refinement uses the actual production cell residual (Hermite
    // bicubic at midpoint vs HEOS) so the cell layout matches what the
    // user actually sees, not a 1D proxy.
    std::vector<double> adaptive_xvec;
    // SUPER region needs the melting curve T_melt(P) as the cold-T
    // boundary; using T_min_eos as a flat approximation in the adaptive
    // grid samples HEOS in the solid extrapolation zone for high-p
    // queries on steep-melting fluids (water).  Until the melting-curve
    // path is wired through bounds_at_p / corner data, restrict adaptive
    // refinement to LIQUID and VAPOR.
    if (table.region() == NormalizedPHTable::LIQUID || table.region() == NormalizedPHTable::VAPOR) {
        auto grid = build_adaptive_normph_grid(static_cast<int>(table.region()), table.ymin, table.ymax, /*Nx_max=*/2 * table.Nx,
                                               /*Ny_max=*/4 * table.Ny, /*tol_relerr=*/1e-4);
        if (grid.first.size() >= 4 && grid.second.size() >= 4) {
            adaptive_xvec = std::move(grid.first);
            adaptive_yvec = std::move(grid.second);
            table.Nx = adaptive_xvec.size();
            table.Ny = adaptive_yvec.size();
        }
    }

    table.resize(table.Nx, table.Ny);  // allocates LIST_OF_MATRICES storage; fills xvec, yvec.
    if (!adaptive_yvec.empty() && adaptive_yvec.size() == table.yvec.size()) {
        table.yvec = std::move(adaptive_yvec);  // override log-uniform with adaptive layout
    }
    if (!adaptive_xvec.empty() && adaptive_xvec.size() == table.xvec.size()) {
        table.xvec = std::move(adaptive_xvec);
    }

    // Same cusp-concentration trick on the xnorm (η) axis for subcritical
    // regions: ρ_sat,L(p) has its sqrt-cusp at η=1 in LIQUID; ρ_sat,V(p)
    // has it at η=0 in VAPOR.  Concentrating ~60 % of cells in the half
    // of [0, 1] adjacent to the sat boundary shrinks the cell η-span
    // there and brings the cell-level Hermite bicubic closer to the cusp
    // without overshooting elsewhere.
    //
    // SUPER xvec stays uniform.  Compressed-liquid (low η) and
    // supercritical-gas extension (high η) span fundamentally different
    // output-variable ranges: ρ varies by ~3× across the L half but
    // ~1000× across the G half.  A single Hermite bicubic on a uniform
    // xvec cannot represent both halves to G13-15 precision without per-
    // half output scaling (linear-ρ vs log-ρ).  Closing the residual
    // SUPER error needs Kunick's L/G split + per-property transforms —
    // out of scope here, queued in dev/sbtl_pt_outstanding_work.md.
    // xvec is now set by the joint build_adaptive_normph_grid above; the
    // legacy single-axis build_adaptive_xvec call is no longer needed
    // for LIQUID/VAPOR regions.
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

// ---------------------------------------------------------------------------
// Disk persistence: read/write one (table, coeffs) pair per file.
// Mirrors BICUBIC/TTSE's pattern (msgpack + miniz zlib) but bundled in a
// single file per region so the corner data and per-cell alpha vectors
// can be round-tripped together — SBTL's coeffs aren't reconstructible
// cheaply from corner data alone (the Hermite chain rule calls HEOS
// second_partial_deriv at every grid node), so they're worth persisting.
// ---------------------------------------------------------------------------
namespace {
template <typename T, typename Coeffs>
void sbtl_write_region(const std::string& path_to_tables, const std::string& filename, const T& table, const Coeffs& coeffs) {
    // Pack table and coeffs back-to-back in the same sbuffer; read side
    // unpacks them in the same order.  Avoids needing std::pair<T, Coeffs>
    // to be default-constructible / copy-constructible.
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, table);
    msgpack::pack(sbuf, coeffs);
    const std::string tab_path = path_to_tables + "/" + filename + ".bin";
    const std::string z_path = tab_path + ".z";
    std::vector<char> buffer(sbuf.size() + 32);
    uLong out_size = static_cast<uLong>(buffer.size());
    compress(reinterpret_cast<unsigned char*>(buffer.data()), &out_size, reinterpret_cast<const unsigned char*>(sbuf.data()),
             static_cast<mz_ulong>(sbuf.size()));
    std::ofstream ofs(z_path.c_str(), std::ofstream::binary);
    ofs.write(buffer.data(), out_size);
    if (CoolProp::get_config_bool(SAVE_RAW_TABLES)) {
        std::ofstream raw(tab_path.c_str(), std::ofstream::binary);
        raw.write(sbuf.data(), sbuf.size());
    }
}

template <typename T, typename Coeffs>
bool sbtl_try_read_region(const std::string& path_to_tables, const std::string& filename, T& table, Coeffs& coeffs) {
    const std::string z_path = path_to_tables + "/" + filename + ".bin.z";
    std::vector<char> raw;
    try {
        raw = get_binary_file_contents(z_path.c_str());
    } catch (...) {
        return false;
    }
    std::vector<unsigned char> buf(raw.size() * 5);
    uLong buf_size = static_cast<uLong>(buf.size());
    int code = Z_BUF_ERROR;
    while (code == Z_BUF_ERROR) {
        code = uncompress(buf.data(), &buf_size, reinterpret_cast<unsigned char*>(raw.data()), static_cast<mz_ulong>(raw.size()));
        if (code == Z_BUF_ERROR) {
            buf.resize(buf.size() * 5);
            buf_size = static_cast<uLong>(buf.size());
        }
    }
    if (code != 0) return false;
    try {
        // Unpack the two msgpack objects back-to-back: first the table,
        // then the coeffs grid.  msgpack::unpacker handles the streaming.
        msgpack::unpacker pac;
        pac.reserve_buffer(buf_size);
        std::memcpy(pac.buffer(), buf.data(), buf_size);
        pac.buffer_consumed(buf_size);
        msgpack::object_handle oh;
        if (!pac.next(oh)) return false;
        T loaded_table;
        oh.get().convert(loaded_table);
        if (!pac.next(oh)) return false;
        Coeffs loaded_coeffs;
        oh.get().convert(loaded_coeffs);
        // Validate grid dimensions match the in-memory tables (could
        // change if TABULAR_NX/NY config differs from cached build).
        if (loaded_table.Nx != table.Nx || loaded_table.Ny != table.Ny) {
            return false;
        }
        table = std::move(loaded_table);
        coeffs = std::move(loaded_coeffs);
        // Reinstate axis vectors and X-macro matrix members from the loaded
        // `matrices` map.
        table.unpack();
        return true;
    } catch (...) {
        return false;
    }
}
}  // namespace

bool SBTLBackend::try_load_normph_tables() {
    if (normph_tables_built) return true;
    const std::string path = path_to_tables();
    if (!sbtl_try_read_region(path, "sbtl_normph_liquid", _normph_liquid, _coeffs_normph_liquid)) return false;
    if (!sbtl_try_read_region(path, "sbtl_normph_vapor", _normph_vapor, _coeffs_normph_vapor)) return false;
    if (!sbtl_try_read_region(path, "sbtl_normph_super", _normph_super, _coeffs_normph_super)) return false;
    _normph_liquid.AS = this->AS;
    _normph_vapor.AS = this->AS;
    _normph_super.AS = this->AS;
    normph_tables_built = true;
    return true;
}

bool SBTLBackend::try_load_normpt_tables() {
    if (normpt_tables_built) return true;
    const std::string path = path_to_tables();
    if (!sbtl_try_read_region(path, "sbtl_normpt_liquid", _normpt_liquid, _coeffs_normpt_liquid)) return false;
    if (!sbtl_try_read_region(path, "sbtl_normpt_vapor", _normpt_vapor, _coeffs_normpt_vapor)) return false;
    if (!sbtl_try_read_region(path, "sbtl_normpt_super", _normpt_super, _coeffs_normpt_super)) return false;
    _normpt_liquid.AS = this->AS;
    _normpt_vapor.AS = this->AS;
    _normpt_super.AS = this->AS;
    normpt_tables_built = true;
    return true;
}

void SBTLBackend::write_normph_tables() {
    const std::string path = path_to_tables();
    make_dirs(path);
    // pack() copies the X-macro matrix members into the inherited `matrices`
    // map so they get serialised; mirror what BICUBIC/TTSE do.
    _normph_liquid.pack();
    _normph_vapor.pack();
    _normph_super.pack();
    sbtl_write_region(path, "sbtl_normph_liquid", _normph_liquid, _coeffs_normph_liquid);
    sbtl_write_region(path, "sbtl_normph_vapor", _normph_vapor, _coeffs_normph_vapor);
    sbtl_write_region(path, "sbtl_normph_super", _normph_super, _coeffs_normph_super);
}

void SBTLBackend::write_normpt_tables() {
    const std::string path = path_to_tables();
    make_dirs(path);
    _normpt_liquid.pack();
    _normpt_vapor.pack();
    _normpt_super.pack();
    sbtl_write_region(path, "sbtl_normpt_liquid", _normpt_liquid, _coeffs_normpt_liquid);
    sbtl_write_region(path, "sbtl_normpt_vapor", _normpt_vapor, _coeffs_normpt_vapor);
    sbtl_write_region(path, "sbtl_normpt_super", _normpt_super, _coeffs_normpt_super);
}

void SBTLBackend::build_normph_tables() {
    if (normph_tables_built) return;
    if (try_load_normph_tables()) return;
    build_normph_table(_normph_liquid);
    build_normph_table(_normph_vapor);
    build_normph_table(_normph_super);
    build_bspline_coeffs(_normph_liquid, _coeffs_normph_liquid);
    build_bspline_coeffs(_normph_vapor, _coeffs_normph_vapor);
    build_bspline_coeffs(_normph_super, _coeffs_normph_super);
    normph_tables_built = true;
    try {
        write_normph_tables();
    } catch (...) {
        // Persistence failures aren't fatal — the in-memory tables work.
    }
}

void SBTLBackend::build_normph_hermite_alphas(NormalizedPHTable& table, std::vector<std::vector<CellCoeffs>>& coeffs) {
    if (coeffs.empty()) return;
    if (!this->AS) return;

    // Per-row saturation envelope data: h_lo(p_j), h_hi(p_j), and
    // their d/d(log p) values used to chain-rule (h, p) HEOS partials
    // into the (xnorm, log p) coordinate.
    struct RowSat
    {
        double h_lo{0.0};
        double h_hi{0.0};
        double dh_lo_dlogp{0.0};
        double dh_hi_dlogp{0.0};
        bool valid{false};
    };
    std::vector<RowSat> row_sat(table.Ny);
    for (std::size_t j = 0; j < table.Ny; ++j) {
        const double p = table.yvec[j];
        const double h_lo = table.h_lo_isobar[j];
        const double h_hi = table.h_hi_isobar[j];
        double dh_lo_dlogp = std::numeric_limits<double>::quiet_NaN();
        double dh_hi_dlogp = std::numeric_limits<double>::quiet_NaN();
        // For saturation boundaries (LIQUID h_hi = h_sat,L; VAPOR h_lo =
        // h_sat,V) use the HEOS analytic d/dp via first_saturation_deriv,
        // converted to d/dlogp by multiplication by p.  For non-sat
        // boundaries (LIQUID h_lo = h(T_min, p), VAPOR h_hi = h(T_max_ext, p),
        // SUPER both edges) use the Cheb1D fit's analytic derivative.
        //
        // Pre-fix: only the Cheb path was wired in, leaving sat-side
        // boundaries with NaN derivatives.  Every VAPOR row then failed
        // the Hermite-upgrade validity check, so VAPOR cells kept their
        // (lower-accuracy) plain B-spline polynomial.  Same on the
        // LIQUID sat-L side.  Filling those derivatives unlocks the
        // Hermite upgrade for the entire subcritical table and aligns
        // production with what the adaptive metric tests.
        try {
            if (table.region() == NormalizedPHTable::LIQUID) {
                if (table.h_lo_cheb.valid()) dh_lo_dlogp = table.h_lo_cheb.eval_dlogp(p);
                this->AS->update(PQ_INPUTS, p, 0.0);
                dh_hi_dlogp = p * static_cast<double>(this->AS->first_saturation_deriv(iHmolar, iP));
            } else if (table.region() == NormalizedPHTable::VAPOR) {
                this->AS->update(PQ_INPUTS, p, 1.0);
                dh_lo_dlogp = p * static_cast<double>(this->AS->first_saturation_deriv(iHmolar, iP));
                if (table.h_hi_cheb.valid()) dh_hi_dlogp = table.h_hi_cheb.eval_dlogp(p);
            } else {  // SUPER — both edges are Cheb-fit (no sat curve)
                if (table.h_lo_cheb.valid()) dh_lo_dlogp = table.h_lo_cheb.eval_dlogp(p);
                if (table.h_hi_cheb.valid()) dh_hi_dlogp = table.h_hi_cheb.eval_dlogp(p);
            }
        } catch (...) {
            // Leave NaN; row will be marked invalid below.
        }
        row_sat[j].h_lo = h_lo;
        row_sat[j].h_hi = h_hi;
        row_sat[j].dh_lo_dlogp = dh_lo_dlogp;
        row_sat[j].dh_hi_dlogp = dh_hi_dlogp;
        row_sat[j].valid = ValidNumber(h_lo) && ValidNumber(h_hi) && ValidNumber(dh_lo_dlogp) && ValidNumber(dh_hi_dlogp) && h_hi > h_lo;
    }

    // Per-node corner-derivative cache (lazy).  Each node is shared by up
    // to 4 cells; populate once and reuse to avoid 4× redundant HEOS updates.
    std::vector<std::vector<SBTLCornerDerivs>> corner(table.Nx, std::vector<SBTLCornerDerivs>(table.Ny));
    std::vector<std::vector<bool>> visited(table.Nx, std::vector<bool>(table.Ny, false));

    auto get_corner = [&](std::size_t i, std::size_t j) -> SBTLCornerDerivs& {
        if (visited[i][j]) return corner[i][j];
        visited[i][j] = true;
        if (!row_sat[j].valid) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        const double xnorm = table.xvec[i];
        const double p = table.yvec[j];
        const double h_lo = row_sat[j].h_lo;
        const double h_hi = row_sat[j].h_hi;
        const double h = h_lo + xnorm * (h_hi - h_lo);
        const bool sat_boundary =
          (table.region() == NormalizedPHTable::LIQUID && i == table.Nx - 1) || (table.region() == NormalizedPHTable::VAPOR && i == 0);
        try {
            if (sat_boundary) {
                const double Q = (table.region() == NormalizedPHTable::VAPOR) ? 1.0 : 0.0;
                this->AS->update(PQ_INPUTS, p, Q);
            } else {
                this->AS->update(HmolarP_INPUTS, h, p);
            }
        } catch (...) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        const double q = static_cast<double>(this->AS->Q());
        if (!sat_boundary && q > 0.0 && q < 1.0) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        populate_corner_derivatives(*(this->AS), corner[i][j]);
        return corner[i][j];
    };

    // Pull (f, f_h, f_p, f_hh, f_hp) for a given property from corner derivs.
    auto extract = [](const SBTLCornerDerivs& cd, parameters prop, double& f, double& f_h, double& f_p, double& f_hh, double& f_hp) -> bool {
        switch (prop) {
            case iDmolar:
                f = cd.rhomolar;
                f_h = cd.drho_dh;
                f_p = cd.drho_dp;
                f_hh = cd.d2rho_dh2;
                f_hp = cd.d2rho_dhp;
                break;
            case iT:
                f = cd.T;
                f_h = cd.dT_dh;
                f_p = cd.dT_dp;
                f_hh = cd.d2T_dh2;
                f_hp = cd.d2T_dhp;
                break;
            case iSmolar:
                f = cd.smolar;
                f_h = cd.ds_dh;
                f_p = cd.ds_dp;
                f_hh = cd.d2s_dh2;
                f_hp = cd.d2s_dhp;
                break;
            case iUmolar:
                f = cd.umolar;
                f_h = cd.du_dh;
                f_p = cd.du_dp;
                f_hh = cd.d2u_dh2;
                f_hp = cd.d2u_dhp;
                break;
            default:
                return false;
        }
        return ValidNumber(f) && ValidNumber(f_h) && ValidNumber(f_p) && ValidNumber(f_hh) && ValidNumber(f_hp);
    };

    const std::array<parameters, 4> core_props = {iDmolar, iT, iSmolar, iUmolar};
    const double rho_crit_molar = static_cast<double>(this->AS->rhomolar_critical());
    int hermite_cell_count = 0;
    for (std::size_t i = 0; i + 1 < table.Nx; ++i) {
        for (std::size_t j = 0; j + 1 < table.Ny; ++j) {
            if (!coeffs[i][j].valid()) continue;
            const auto& c00 = get_corner(i, j);
            const auto& c10 = get_corner(i + 1, j);
            const auto& c01 = get_corner(i, j + 1);
            const auto& c11 = get_corner(i + 1, j + 1);
            if (!c00.valid || !c10.valid || !c01.valid || !c11.valid) continue;

            const double Delta_xi = table.xvec[i + 1] - table.xvec[i];
            const double Delta_eta = std::log(table.yvec[j + 1]) - std::log(table.yvec[j]);

            const std::array<double, 4> xn = {table.xvec[i], table.xvec[i + 1], table.xvec[i], table.xvec[i + 1]};
            const std::array<std::size_t, 4> jr = {j, j, j + 1, j + 1};
            const std::array<const SBTLCornerDerivs*, 4> cds = {&c00, &c10, &c01, &c11};

            bool any_property_upgraded = false;
            for (parameters prop : core_props) {
                const SBTLTransform tform = sbtl_transform_for_property(static_cast<int>(table.region()), prop);
                const double y_ref = sbtl_transform_y_ref(static_cast<int>(table.region()), prop, rho_crit_molar);
                // CUSP_DIST derivatives blow up as f → y_ref.  Require a minimum
                // distance of 1e-6·y_ref at every corner so the chain-rule terms
                // stay finite; cells that brush ρ_crit too closely fall back to
                // the existing bilinear coefficients for this property.
                const double cusp_min_d =
                  (tform == SBTLTransform::CUSP_DIST_ABOVE_REF || tform == SBTLTransform::CUSP_DIST_BELOW_REF) ? 1e-6 * std::abs(y_ref) : 0.0;
                std::array<double, 4> f{};
                std::array<double, 4> fxi{};
                std::array<double, 4> feta{};
                std::array<double, 4> fxieta{};
                bool ok = true;
                for (int k = 0; k < 4; ++k) {
                    double f_h = 0.0, f_p = 0.0, f_hh = 0.0, f_hp = 0.0;
                    if (!extract(*cds[k], prop, f[k], f_h, f_p, f_hh, f_hp)) {
                        ok = false;
                        break;
                    }
                    const std::size_t row = jr[k];
                    const double xnorm_k = xn[k];
                    const double p_k = table.yvec[row];
                    const double Delta_h = row_sat[row].h_hi - row_sat[row].h_lo;
                    const double dDelta_h_dlogp = row_sat[row].dh_hi_dlogp - row_sat[row].dh_lo_dlogp;
                    const double dh_dlogp = row_sat[row].dh_lo_dlogp + xnorm_k * dDelta_h_dlogp;

                    // ∂f̃/∂xnorm = f_h · Δh
                    double df_dxnorm = f_h * Delta_h;
                    // ∂f̃/∂(log p) at fixed xnorm = p · f_p + f_h · (dh_lo' + xnorm · dΔh')
                    double df_dlogp = p_k * f_p + f_h * dh_dlogp;
                    // ∂²f̃/∂xnorm∂(log p) = Δh · [f_hh · dh_dlogp + p · f_hp] + f_h · dΔh'
                    double d2f = Delta_h * (f_hh * dh_dlogp + p_k * f_hp) + f_h * dDelta_h_dlogp;

                    // Chain-rule the corner derivatives through the output
                    // transform g.  For F = g(f):
                    //   F_x   = g'(f) · f_x
                    //   F_y   = g'(f) · f_y
                    //   F_xy  = g''(f) · f_x · f_y + g'(f) · f_xy
                    if (tform != SBTLTransform::IDENTITY) {
                        if (!(f[k] > 0.0)) {
                            ok = false;
                            break;
                        }  // log/cusp require y>0
                        if (tform == SBTLTransform::CUSP_DIST_ABOVE_REF && (f[k] - y_ref) <= cusp_min_d) {
                            ok = false;
                            break;
                        }
                        if (tform == SBTLTransform::CUSP_DIST_BELOW_REF && (y_ref - f[k]) <= cusp_min_d) {
                            ok = false;
                            break;
                        }
                        const double gp = sbtl_transform_deriv(tform, f[k], y_ref);
                        const double gpp = sbtl_transform_second_deriv(tform, f[k], y_ref);
                        const double new_d2 = gpp * df_dxnorm * df_dlogp + gp * d2f;
                        df_dxnorm *= gp;
                        df_dlogp *= gp;
                        d2f = new_d2;
                        f[k] = sbtl_transform_apply(tform, f[k], y_ref);
                    }

                    fxi[k] = df_dxnorm * Delta_xi;
                    feta[k] = df_dlogp * Delta_eta;
                    fxieta[k] = d2f * Delta_xi * Delta_eta;
                }
                if (!ok) continue;
                const auto alpha = hermite_bicubic_polynomial_coeffs(f[0], f[1], f[2], f[3], fxi[0], fxi[1], fxi[2], fxi[3], feta[0], feta[1],
                                                                     feta[2], feta[3], fxieta[0], fxieta[1], fxieta[2], fxieta[3]);
                coeffs[i][j].set(prop, alpha);
                any_property_upgraded = true;
            }
            if (any_property_upgraded) ++hermite_cell_count;
        }
    }
    if (get_debug_level() > 0) {
        std::cout << "SBTL: upgraded " << hermite_cell_count << " cells to Hermite bicubic for normph region " << static_cast<int>(table.region())
                  << "\n";
    }
}

void SBTLBackend::build_normpt_hermite_alphas(NormalizedPTTable& table, std::vector<std::vector<CellCoeffs>>& coeffs) {
    if (coeffs.empty()) return;
    if (!this->AS) return;

    // Per-row sat envelope data on the T-axis.  All derivatives are analytic:
    //   * sat boundaries  (LIQUID T_hi = T_sat,L, VAPOR T_lo = T_sat,V)
    //     → HEOS-analytic via first_saturation_deriv(iT, iP), scaled by p.
    //   * non-sat boundaries (LIQUID T_lo, VAPOR T_hi, SUPER both)
    //     → analytic derivative of the Cheb1D fit through HEOS-sampled
    //       values (eval_dlogp).  Pairs with h_*_cheb.eval at lookup time
    //       so the build and lookup paths see the same boundary function.
    //
    // No FD fallback.  If neither source supplies a finite derivative the
    // row is marked invalid and falls back to the cubic-B-spline backbone
    // (which is the same behaviour as any other row with non-finite corner
    // data).
    struct RowSat
    {
        double T_lo{0.0};
        double T_hi{0.0};
        double dT_lo_dlogp{0.0};
        double dT_hi_dlogp{0.0};
        bool valid{false};
    };
    std::vector<RowSat> row_sat(table.Ny);
    for (std::size_t j = 0; j < table.Ny; ++j) {
        const double p = table.yvec[j];
        const double T_lo = table.T_lo_isobar[j];
        const double T_hi = table.T_hi_isobar[j];
        const bool lo_is_sat = (table.region() == NormalizedPTTable::VAPOR);
        const bool hi_is_sat = (table.region() == NormalizedPTTable::LIQUID);
        double dT_lo_dlogp = std::numeric_limits<double>::quiet_NaN();
        double dT_hi_dlogp = std::numeric_limits<double>::quiet_NaN();
        try {
            if (lo_is_sat) {
                this->AS->update(PQ_INPUTS, p, 1.0);
                dT_lo_dlogp = p * static_cast<double>(this->AS->first_saturation_deriv(iT, iP));
            } else if (table.T_lo_cheb.valid()) {
                dT_lo_dlogp = table.T_lo_cheb.eval_dlogp(p);
            }
            if (hi_is_sat) {
                this->AS->update(PQ_INPUTS, p, 0.0);
                dT_hi_dlogp = p * static_cast<double>(this->AS->first_saturation_deriv(iT, iP));
            } else if (table.T_hi_cheb.valid()) {
                dT_hi_dlogp = table.T_hi_cheb.eval_dlogp(p);
            }
        } catch (...) {
            // Leave NaN; row will be marked invalid below.
        }
        row_sat[j].T_lo = T_lo;
        row_sat[j].T_hi = T_hi;
        row_sat[j].dT_lo_dlogp = dT_lo_dlogp;
        row_sat[j].dT_hi_dlogp = dT_hi_dlogp;
        row_sat[j].valid = ValidNumber(T_lo) && ValidNumber(T_hi) && ValidNumber(dT_lo_dlogp) && ValidNumber(dT_hi_dlogp) && T_hi > T_lo;
    }

    std::vector<std::vector<SBTLCornerDerivsPT>> corner(table.Nx, std::vector<SBTLCornerDerivsPT>(table.Ny));
    std::vector<std::vector<bool>> visited(table.Nx, std::vector<bool>(table.Ny, false));

    auto get_corner = [&](std::size_t i, std::size_t j) -> SBTLCornerDerivsPT& {
        if (visited[i][j]) return corner[i][j];
        visited[i][j] = true;
        if (!row_sat[j].valid) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        const double xnorm = table.xvec[i];
        const double p = table.yvec[j];
        const double T_lo = row_sat[j].T_lo;
        const double T_hi = row_sat[j].T_hi;
        const double T = T_lo + xnorm * (T_hi - T_lo);
        const bool sat_boundary =
          (table.region() == NormalizedPTTable::LIQUID && i == table.Nx - 1) || (table.region() == NormalizedPTTable::VAPOR && i == 0);
        try {
            if (sat_boundary) {
                const double Q = (table.region() == NormalizedPTTable::VAPOR) ? 1.0 : 0.0;
                this->AS->update(PQ_INPUTS, p, Q);
            } else {
                this->AS->update(PT_INPUTS, p, T);
            }
        } catch (...) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        const double q = static_cast<double>(this->AS->Q());
        if (!sat_boundary && q > 0.0 && q < 1.0) {
            corner[i][j].valid = false;
            return corner[i][j];
        }
        populate_corner_derivatives_pt(*(this->AS), corner[i][j]);
        return corner[i][j];
    };

    auto extract = [](const SBTLCornerDerivsPT& cd, parameters prop, double& f, double& f_T, double& f_p, double& f_TT, double& f_Tp) -> bool {
        switch (prop) {
            case iDmolar:
                f = cd.rhomolar;
                f_T = cd.drho_dT;
                f_p = cd.drho_dp;
                f_TT = cd.d2rho_dT2;
                f_Tp = cd.d2rho_dTp;
                break;
            case iHmolar:
                f = cd.hmolar;
                f_T = cd.dh_dT;
                f_p = cd.dh_dp;
                f_TT = cd.d2h_dT2;
                f_Tp = cd.d2h_dTp;
                break;
            case iSmolar:
                f = cd.smolar;
                f_T = cd.ds_dT;
                f_p = cd.ds_dp;
                f_TT = cd.d2s_dT2;
                f_Tp = cd.d2s_dTp;
                break;
            case iUmolar:
                f = cd.umolar;
                f_T = cd.du_dT;
                f_p = cd.du_dp;
                f_TT = cd.d2u_dT2;
                f_Tp = cd.d2u_dTp;
                break;
            default:
                return false;
        }
        return ValidNumber(f) && ValidNumber(f_T) && ValidNumber(f_p) && ValidNumber(f_TT) && ValidNumber(f_Tp);
    };

    const std::array<parameters, 4> core_props = {iDmolar, iHmolar, iSmolar, iUmolar};
    const double rho_crit_molar = static_cast<double>(this->AS->rhomolar_critical());
    int hermite_cell_count = 0;
    for (std::size_t i = 0; i + 1 < table.Nx; ++i) {
        for (std::size_t j = 0; j + 1 < table.Ny; ++j) {
            if (!coeffs[i][j].valid()) continue;
            const auto& c00 = get_corner(i, j);
            const auto& c10 = get_corner(i + 1, j);
            const auto& c01 = get_corner(i, j + 1);
            const auto& c11 = get_corner(i + 1, j + 1);
            if (!c00.valid || !c10.valid || !c01.valid || !c11.valid) continue;

            const double Delta_xi = table.xvec[i + 1] - table.xvec[i];
            const double Delta_eta = std::log(table.yvec[j + 1]) - std::log(table.yvec[j]);

            const std::array<double, 4> xn = {table.xvec[i], table.xvec[i + 1], table.xvec[i], table.xvec[i + 1]};
            const std::array<std::size_t, 4> jr = {j, j, j + 1, j + 1};
            const std::array<const SBTLCornerDerivsPT*, 4> cds = {&c00, &c10, &c01, &c11};

            bool any_property_upgraded = false;
            for (parameters prop : core_props) {
                const SBTLTransform tform = sbtl_transform_for_property(static_cast<int>(table.region()), prop);
                const double y_ref = sbtl_transform_y_ref(static_cast<int>(table.region()), prop, rho_crit_molar);
                const double cusp_min_d =
                  (tform == SBTLTransform::CUSP_DIST_ABOVE_REF || tform == SBTLTransform::CUSP_DIST_BELOW_REF) ? 1e-6 * std::abs(y_ref) : 0.0;
                std::array<double, 4> f{};
                std::array<double, 4> fxi{};
                std::array<double, 4> feta{};
                std::array<double, 4> fxieta{};
                bool ok = true;
                for (int k = 0; k < 4; ++k) {
                    double f_T = 0.0, f_p = 0.0, f_TT = 0.0, f_Tp = 0.0;
                    if (!extract(*cds[k], prop, f[k], f_T, f_p, f_TT, f_Tp)) {
                        ok = false;
                        break;
                    }
                    const std::size_t row = jr[k];
                    const double xnorm_k = xn[k];
                    const double p_k = table.yvec[row];
                    const double Delta_T = row_sat[row].T_hi - row_sat[row].T_lo;
                    const double dDelta_T_dlogp = row_sat[row].dT_hi_dlogp - row_sat[row].dT_lo_dlogp;
                    const double dT_dlogp = row_sat[row].dT_lo_dlogp + xnorm_k * dDelta_T_dlogp;

                    double df_dxnorm = f_T * Delta_T;
                    double df_dlogp = p_k * f_p + f_T * dT_dlogp;
                    double d2f = Delta_T * (f_TT * dT_dlogp + p_k * f_Tp) + f_T * dDelta_T_dlogp;

                    // Chain-rule corner derivatives through the output
                    // transform g (same form as PH path; see comment there).
                    if (tform != SBTLTransform::IDENTITY) {
                        if (!(f[k] > 0.0)) {
                            ok = false;
                            break;
                        }
                        if (tform == SBTLTransform::CUSP_DIST_ABOVE_REF && (f[k] - y_ref) <= cusp_min_d) {
                            ok = false;
                            break;
                        }
                        if (tform == SBTLTransform::CUSP_DIST_BELOW_REF && (y_ref - f[k]) <= cusp_min_d) {
                            ok = false;
                            break;
                        }
                        const double gp = sbtl_transform_deriv(tform, f[k], y_ref);
                        const double gpp = sbtl_transform_second_deriv(tform, f[k], y_ref);
                        const double new_d2 = gpp * df_dxnorm * df_dlogp + gp * d2f;
                        df_dxnorm *= gp;
                        df_dlogp *= gp;
                        d2f = new_d2;
                        f[k] = sbtl_transform_apply(tform, f[k], y_ref);
                    }

                    fxi[k] = df_dxnorm * Delta_xi;
                    feta[k] = df_dlogp * Delta_eta;
                    fxieta[k] = d2f * Delta_xi * Delta_eta;
                }
                if (!ok) continue;
                const auto alpha = hermite_bicubic_polynomial_coeffs(f[0], f[1], f[2], f[3], fxi[0], fxi[1], fxi[2], fxi[3], feta[0], feta[1],
                                                                     feta[2], feta[3], fxieta[0], fxieta[1], fxieta[2], fxieta[3]);
                coeffs[i][j].set(prop, alpha);
                any_property_upgraded = true;
            }
            if (any_property_upgraded) ++hermite_cell_count;
        }
    }
    if (get_debug_level() > 0) {
        std::cout << "SBTL: upgraded " << hermite_cell_count << " cells to Hermite bicubic for normpt region " << static_cast<int>(table.region())
                  << "\n";
    }
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

double NormalizedPTTable::xnorm_from_T(double T, double P, std::optional<double> T_sat) const {
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
            T_hi = T_sat ? *T_sat : cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case VAPOR:
            T_lo = T_sat ? *T_sat : cheb_or_isobar(T_lo_cheb, T_lo_isobar);
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

double NormalizedPTTable::T_from_xnorm(double xnorm, double P, std::optional<double> T_sat) const {
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
            T_hi = T_sat ? *T_sat : cheb_or_isobar(T_hi_cheb, T_hi_isobar);
            break;
        case VAPOR:
            T_lo = T_sat ? *T_sat : cheb_or_isobar(T_lo_cheb, T_lo_isobar);
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
    if (tsat_LV_cache_p == p) {
        T_L = tsat_LV_cache_TL;
        T_V = tsat_LV_cache_TV;
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
    tsat_LV_cache_p = p;
    tsat_LV_cache_TL = T_L;
    tsat_LV_cache_TV = T_V;
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
    // Piecewise Chebyshev with low degree per piece (8) and more pieces
    // near p_crit.  Mirrors the PH refactor; same rationale (cheaper eval,
    // better cusp resolution).  T_hi is constant (T_max_ext), so for SUPER
    // a single piece is sufficient; T_melt(p) is smooth so few pieces also
    // suffice there.  LIQUID/VAPOR's T_sat boundary uses PQ_INPUTS at
    // lookup, not this Cheb, so the cusp resolution here matters less than
    // in the PH path — but uniform structure is easier to maintain.
    std::vector<double> breakpoints;
    if (table.region() == NormalizedPTTable::SUPER) {
        breakpoints = {p_lo, p_hi};
    } else if (this->AS) {
        const double p_crit = static_cast<double>(this->AS->p_critical());
        breakpoints = {p_lo,          0.10 * p_crit, 0.30 * p_crit,  0.50 * p_crit,  0.70 * p_crit, 0.85 * p_crit,
                       0.92 * p_crit, 0.96 * p_crit, 0.985 * p_crit, 0.995 * p_crit, p_hi};
        breakpoints.erase(std::remove_if(breakpoints.begin() + 1, breakpoints.end() - 1, [&](double v) { return v <= p_lo || v >= p_hi; }),
                          breakpoints.end() - 1);
    } else {
        breakpoints = {p_lo, p_hi};
    }
    const std::size_t N_cheb = 8;
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

    // Mirror of the PH adaptive-yvec block above.  PT's saturation
    // boundary is T_sat(p) (smooth where ρ_sat has its cusp), so refinement
    // here matters less than in PH — but for symmetry and so the persisted
    // cache layout matches across coord systems, do the same dance.
    std::vector<double> adaptive_yvec;
    if (table.region() != NormalizedPTTable::SUPER) {
        auto T_sat_at_p = [this](double p) -> double {
            double T_L = NAN, T_V = NAN;
            this->saturation_T_LV(p, T_L, T_V);
            return T_L;  // L == V for pure fluids
        };
        adaptive_yvec = build_adaptive_yvec(table.ymin, table.ymax, table.Ny, /*Ny_max=*/10 * table.Ny, T_sat_at_p,
                                            /*tol_rel=*/0.005);
        if (adaptive_yvec.size() >= 4) {
            table.Ny = adaptive_yvec.size();
        }
    }

    // SUPER region excluded for the same reason as PH path — the constant
    // T_lo approximation breaks for high-p queries on water (samples in
    // the solid extrapolation zone), making the adaptive grid worse than
    // the default log-uniform.  Restrict to LIQUID and VAPOR.
    std::vector<double> adaptive_xvec;
    if (table.region() == NormalizedPTTable::LIQUID || table.region() == NormalizedPTTable::VAPOR) {
        auto grid = build_adaptive_normpt_grid(static_cast<int>(table.region()), table.ymin, table.ymax, /*Nx_max=*/2 * table.Nx,
                                               /*Ny_max=*/4 * table.Ny, /*tol_relerr=*/1e-4);
        if (grid.first.size() >= 4 && grid.second.size() >= 4) {
            adaptive_xvec = std::move(grid.first);
            adaptive_yvec = std::move(grid.second);
            table.Nx = adaptive_xvec.size();
            table.Ny = adaptive_yvec.size();
        }
    }

    table.resize(table.Nx, table.Ny);
    if (!adaptive_yvec.empty() && adaptive_yvec.size() == table.yvec.size()) {
        table.yvec = std::move(adaptive_yvec);
    }
    if (!adaptive_xvec.empty() && adaptive_xvec.size() == table.xvec.size()) {
        table.xvec = std::move(adaptive_xvec);
    }
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
    if (try_load_normpt_tables()) return;
    build_normpt_table(_normpt_liquid);
    build_normpt_table(_normpt_vapor);
    build_normpt_table(_normpt_super);
    build_bspline_coeffs(_normpt_liquid, _coeffs_normpt_liquid);
    build_bspline_coeffs(_normpt_vapor, _coeffs_normpt_vapor);
    build_bspline_coeffs(_normpt_super, _coeffs_normpt_super);
    normpt_tables_built = true;
    try {
        write_normpt_tables();
    } catch (...) {
        // Persistence failures aren't fatal — the in-memory tables work.
    }
}

void SBTLBackend::build_bspline_coeffs(SinglePhaseGriddedTableData& table, std::vector<std::vector<CellCoeffs>>& coeffs) {
    if (!coeffs.empty()) return;

    // Core thermodynamic params drive cell validity.  Auxiliary params
    // (speed_sound / viscosity / conductivity) build their own polynomial
    // alongside but DON'T affect cell validity — at sat-curve boundaries
    // and metastable spinodals the EOS may throw for speed_sound or
    // transport, and we don't want one missing corner to disqualify the
    // entire cell from the thermodynamic-property lookups.
    const std::array<parameters, 6> core_params = {iDmolar, iT, iSmolar, iHmolar, iP, iUmolar};
    const std::array<parameters, 3> aux_params = {ispeed_sound, iviscosity, iconductivity};
    const int core_count = static_cast<int>(core_params.size());
    const int aux_count = static_cast<int>(aux_params.size());

    coeffs.resize(table.Nx - 1, std::vector<CellCoeffs>(table.Ny - 1));

    // Hermite bicubic overlay: for a NormalizedPHTable we route the 4 core
    // derived params (rho, T, s, u) through a Hermite bicubic path that
    // uses EOS-supplied first and cross derivatives chain-ruled into
    // (xnorm, log p); for a NormalizedPTTable, the analogous 4 props
    // (rho, h, s, u) route through (xnorm_T, log p).  Aux + coordinate-
    // aligned params (h/p or T/p, w, η, λ) keep their cubic-B-spline alpha.
    // Cells where Hermite chain-rule throws fall back to cubic alpha
    // computed in the loop below — single-cell fallback so cell.valid()
    // always implies a populated polynomial.
    auto* normph_for_hermite = dynamic_cast<NormalizedPHTable*>(&table);
    auto* normpt_for_hermite = dynamic_cast<NormalizedPTTable*>(&table);

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
        // Cubic B-spline runs for every cell × every param.  The Hermite
        // pass below overlays its alpha for cells where corner-derivative
        // data is available; cells where Hermite fails (e.g. EOS partials
        // throw) keep their cubic alpha so cell `valid()` still implies a
        // populated polynomial — no eval-time empty-vector crashes.
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

        // Output-variable scaling: for properties whose region demands a
        // non-identity transform (ρ in VAPOR / SUPER => log), build the
        // spline on the TRANSFORMED matrix.  The stored alpha coefficients
        // then encode the polynomial fit to g(value) — log(ρ) in this
        // phase — and the lookup path applies g⁻¹ to recover the physical
        // value.  Skips properties whose transform is IDENTITY (no-op).
        const int region_int_for_transform = [&]() -> int {
            if (auto* ph = dynamic_cast<NormalizedPHTable*>(&table)) return static_cast<int>(ph->region());
            if (auto* pt = dynamic_cast<NormalizedPTTable*>(&table)) return static_cast<int>(pt->region());
            return -1;  // unknown — caller is some other table type; identity
        }();
        const SBTLTransform tform =
          (region_int_for_transform >= 0) ? sbtl_transform_for_property(region_int_for_transform, param) : SBTLTransform::IDENTITY;
        const double y_ref_src = (region_int_for_transform >= 0)
                                   ? sbtl_transform_y_ref(region_int_for_transform, param, static_cast<double>(this->AS->rhomolar_critical()))
                                   : 0.0;
        std::vector<std::vector<double>> transform_scratch;
        if (tform != SBTLTransform::IDENTITY) {
            transform_scratch = is_aux ? aux_scratch : *fp;
            for (auto& row : transform_scratch) {
                for (double& v : row) {
                    bool on_domain = ValidNumber(v) && v > 0.0;
                    if (on_domain && tform == SBTLTransform::CUSP_DIST_ABOVE_REF && !(v > y_ref_src)) on_domain = false;
                    if (on_domain && tform == SBTLTransform::CUSP_DIST_BELOW_REF && !(v < y_ref_src)) on_domain = false;
                    if (on_domain) {
                        v = sbtl_transform_apply(tform, v, y_ref_src);
                    } else {
                        v = _HUGE;  // mark so per-cell B-spline skips
                    }
                }
            }
        }
        const auto& source = (tform != SBTLTransform::IDENTITY) ? transform_scratch : (is_aux ? aux_scratch : *fp);

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

    // Hermite bicubic pass.  Cells with valid corner derivatives get
    // Hermite alpha overlaid on top of the cubic alpha already populated
    // by the loop above; cells where EOS partials fail (sat boundary,
    // metastable spinodal, near-critical cusp) keep cubic so cell.valid()
    // always implies a populated polynomial.
    if (normph_for_hermite != nullptr) {
        build_normph_hermite_alphas(*normph_for_hermite, coeffs);
    } else if (normpt_for_hermite != nullptr) {
        build_normpt_hermite_alphas(*normpt_for_hermite, coeffs);
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

double NormalizedPHTable::h_from_xnorm(double xnorm, double P, std::optional<double> h_sat) const {
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
    double h_lo = 0.0, h_hi = 0.0;
    switch (region_) {
        case LIQUID:
            h_lo = cheb_or_isobar(h_lo_cheb, h_lo_isobar);
            h_hi = h_sat ? *h_sat : cheb_or_isobar(h_hi_cheb, h_hi_isobar);
            break;
        case VAPOR:
            h_lo = h_sat ? *h_sat : cheb_or_isobar(h_lo_cheb, h_lo_isobar);
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

double SBTLBackend::evaluate_single_phase_pre(const std::vector<std::vector<CellCoeffs>>& coeffs, const parameters output, const double xi,
                                              const double eta, const std::size_t i, const std::size_t j) {
    const std::vector<double>& alpha = coeffs[i][j].get(output);
    const double xi2 = xi * xi, xi3 = xi2 * xi;
    const double eta2 = eta * eta, eta3 = eta2 * eta;
    const double B0 = alpha[0] + alpha[1] * eta + alpha[2] * eta2 + alpha[3] * eta3;
    const double B1 = alpha[4] + alpha[5] * eta + alpha[6] * eta2 + alpha[7] * eta3;
    const double B2 = alpha[8] + alpha[9] * eta + alpha[10] * eta2 + alpha[11] * eta3;
    const double B3 = alpha[12] + alpha[13] * eta + alpha[14] * eta2 + alpha[15] * eta3;
    return B0 + B1 * xi + B2 * xi2 + B3 * xi3;
}

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
                } else {
                    // h_sat_L < h_molar < h_sat_V → strictly inside the dome.
                    // PT_INPUTS at the analogous state already throws a clean
                    // NotImplementedError (see normpt branch); mirror that
                    // here so callers don't get the misleading "input pair
                    // HmolarP_INPUTS not supported" terminal throw at the
                    // bottom of update().
                    throw NotImplementedError(
                      format("SBTL HmolarP/HmassP_INPUTS at (h=%g J/mol, P=%g Pa) is strictly inside the two-phase dome "
                             "(h_sat,L=%g, h_sat,V=%g).  Use PQ_INPUTS or specify_phase(iphase_liquid|iphase_gas) to disambiguate.",
                             h_molar, p, h_sat_L, h_sat_V));
                }
            }
        }

        if (tbl) {
            // Per-table OOB pressure guard — mirror of the PT path's check.
            // Without this, an OOB p with a coincidentally in-range xnorm
            // would silently extrapolate the bicubic, and an OOB p whose
            // sat-lookup throws would reach the unconditional final-throw
            // below with the misleading "input pair not supported" message.
            if (p < tbl->yvec.front() || p > tbl->yvec.back()) {
                throw ValueError(format("SBTL HmolarP/HmassP_INPUTS at P=%g Pa is outside the normph "
                                        "table range [%g, %g] Pa — use HEOS directly.",
                                        p, tbl->yvec.front(), tbl->yvec.back()));
            }
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
                        // Compute in-cell (xi, eta) once and reuse across
                        // every property evaluation.  Caches feed both the
                        // eager T/rho here and the lazy s/u/μ/k accessors.
                        const double xi = (xnorm - tbl->xvec[i]) / (tbl->xvec[i + 1] - tbl->xvec[i]);
                        const double eta = (std::log(p) - std::log(tbl->yvec[j])) / (std::log(tbl->yvec[j + 1]) - std::log(tbl->yvec[j]));
                        _normph_xi = xi;
                        _normph_eta = eta;
                        const int region_int = static_cast<int>(tbl->region());
                        const double rho_crit_inv = static_cast<double>(this->AS->rhomolar_critical());
                        const SBTLTransform tform_T = sbtl_transform_for_property(region_int, iT);
                        const SBTLTransform tform_rho = sbtl_transform_for_property(region_int, iDmolar);
                        const double y_ref_T = sbtl_transform_y_ref(region_int, iT, rho_crit_inv);
                        const double y_ref_rho = sbtl_transform_y_ref(region_int, iDmolar, rho_crit_inv);
                        const double T_val = sbtl_transform_inverse(tform_T, evaluate_single_phase_pre(*coeffs, iT, xi, eta, i, j), y_ref_T);
                        const double rho_val =
                          sbtl_transform_inverse(tform_rho, evaluate_single_phase_pre(*coeffs, iDmolar, xi, eta, i, j), y_ref_rho);
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
            } catch (const NotImplementedError&) {
                throw;  // Re-throw user-facing errors (in-dome HmassP).
            } catch (const std::exception& e) {
                // Internal failure: re-throw with context so the caller sees
                // the real problem rather than the misleading "input pair
                // not supported" terminal throw at the end of update().
                throw ValueError(format("SBTL HmolarP/HmassP_INPUTS at (h=%g J/mol, P=%g Pa) could not be "
                                        "evaluated through the normalized PH path: %s.  Use HEOS directly.",
                                        h_molar, p, e.what()));
            }
            // Reached only when xnorm fell outside [0,1] or the cell at the
            // resolved (i,j) was marked invalid (e.g. table edge with no
            // neighbour-recovery).  Surface a clean error rather than
            // dropping to the misleading "input pair not supported" throw.
            // Re-eval xnorm in a try/catch so a sat-lookup race here doesn't
            // replace the diagnostic with whatever xnorm_from_h raises.
            {
                double xnorm_d = std::numeric_limits<double>::quiet_NaN();
                try {
                    xnorm_d = tbl->xnorm_from_h(h_molar, p, h_sat_active);
                } catch (...) {
                    // leave xnorm_d as NaN; message still useful
                }
                throw ValueError(format("SBTL HmolarP/HmassP_INPUTS at (h=%g J/mol, P=%g Pa) is outside the "
                                        "normph table coverage (xnorm=%g, p_range=[%g,%g]).  "
                                        "Use HEOS directly.",
                                        h_molar, p, xnorm_d, tbl->yvec.front(), tbl->yvec.back()));
            }
        } else {
            // tbl == nullptr: routing fell through without picking a table.
            // Most likely cause: saturation_hmolar_LV(p) threw (sat lookup
            // race), and forced==Auto with p<p_crit so all three regions
            // were skipped.  Surface a clean error instead of falling
            // through to the terminal "input pair not supported" message.
            throw ValueError(format("SBTL HmolarP/HmassP_INPUTS at (h=%g J/mol, P=%g Pa) could not be routed "
                                    "to a normalized table — saturation lookup likely failed at this pressure. "
                                    "Use HEOS directly.",
                                    h_molar, p));
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
            } catch (const std::exception& e) {
                // PQ_INPUTS solver failed inside saturation_T_LV (rare:
                // numerical convergence wobble near p_crit, or stored
                // ancillary issue).  Don't fall through to the misleading
                // "input pair not supported" terminal throw; surface what
                // actually failed.
                throw ValueError(format("SBTL PT_INPUTS at P=%g Pa: could not resolve T_sat via PQ_INPUTS — %s.  "
                                        "Use HEOS directly for this state.",
                                        p, e.what()));
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
                        // Precompute (xi, eta) once; cache for lazy accessors.
                        const double xi = (xnorm - tbl->xvec[i]) / (tbl->xvec[i + 1] - tbl->xvec[i]);
                        const double eta = (std::log(p) - std::log(tbl->yvec[j])) / (std::log(tbl->yvec[j + 1]) - std::log(tbl->yvec[j]));
                        _normpt_xi = xi;
                        _normpt_eta = eta;
                        const int region_int_pt = static_cast<int>(tbl->region());
                        const SBTLTransform tform_rho_pt = sbtl_transform_for_property(region_int_pt, iDmolar);
                        const double y_ref_rho_pt = sbtl_transform_y_ref(region_int_pt, iDmolar, static_cast<double>(this->AS->rhomolar_critical()));
                        const double rho_val =
                          sbtl_transform_inverse(tform_rho_pt, evaluate_single_phase_pre(*coeffs, iDmolar, xi, eta, i, j), y_ref_rho_pt);
                        if (!ValidNumber(rho_val) || rho_val <= 0.0) {
                            throw ValueError("SBTL normpt polynomial gave non-physical rho");
                        }
                        _rhomolar = rho_val;
                        // h, s, u deferred to lazy accessors (calc_hmolar /
                        // calc_smolar / calc_umolar via the _normpt_active_table
                        // override).  Eager evaluation here was burning ~450 ns
                        // per update on three Hermite-bicubic polynomial evals
                        // that are usually unused — most PT consumers only ask
                        // for rho.  Matches the PH path's lazy contract.
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
CoolPropDbl SBTLBackend::calc_hmolar() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_hmolar);
    if (_normph_active_table) return static_cast<CoolPropDbl>(_hmolar);
    if (_normpt_active_table) {
        return static_cast<CoolPropDbl>(
          evaluate_single_phase_pre(*_normpt_active_coeffs, iHmolar, _normpt_xi, _normpt_eta, cached_single_phase_i, cached_single_phase_j));
    }
    return TabularBackend::calc_hmolar();
}
CoolPropDbl SBTLBackend::calc_smolar() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_smolar);
    if (_normph_active_table) {
        return static_cast<CoolPropDbl>(
          evaluate_single_phase_pre(*_normph_active_coeffs, iSmolar, _normph_xi, _normph_eta, cached_single_phase_i, cached_single_phase_j));
    }
    if (_normpt_active_table) {
        return static_cast<CoolPropDbl>(
          evaluate_single_phase_pre(*_normpt_active_coeffs, iSmolar, _normpt_xi, _normpt_eta, cached_single_phase_i, cached_single_phase_j));
    }
    return TabularBackend::calc_smolar();
}
CoolPropDbl SBTLBackend::calc_umolar() {
    if (_critbox_active) return static_cast<CoolPropDbl>(_critbox_umolar);
    if (_normph_active_table) {
        return static_cast<CoolPropDbl>(
          evaluate_single_phase_pre(*_normph_active_coeffs, iUmolar, _normph_xi, _normph_eta, cached_single_phase_i, cached_single_phase_j));
    }
    if (_normpt_active_table) {
        return static_cast<CoolPropDbl>(
          evaluate_single_phase_pre(*_normpt_active_coeffs, iUmolar, _normpt_xi, _normpt_eta, cached_single_phase_i, cached_single_phase_j));
    }
    return TabularBackend::calc_umolar();
}
CoolPropDbl SBTLBackend::calc_rhomolar() {
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
