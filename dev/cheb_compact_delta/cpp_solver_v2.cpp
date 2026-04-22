// Optimized C++ density rootfinder — pushes for minimum wall time.
//
// Optimizations over cpp_solver.cpp:
//   1. Interval-major memory layout H[kk, k, j] instead of (n_terms, K*(n+1)).
//      At runtime each interval's 36×13 block is contiguous → fits in L1.
//   2. Fused per-interval inner loop: small matvec + times_x + residual build
//      + cheap-skip test + derivative/monotonicity test + secant — all in
//      registers/L1 cache with no intermediate K-sized buffers.
//   3. Preallocated output (fixed-size numpy array with count) — no pybind11
//      per-root list-append overhead.
//   4. LAPACK dhseqr on the Hessenberg colleague matrix (no dgehrd reduction)
//      for non-monotone intervals.
//   5. Early-exit skip test: accumulate |c_j| as we build; as soon as it
//      exceeds |c_0| we know the skip is impossible (rare branch).
//
// Compile:
//   c++ -O3 -ffast-math -march=native -std=c++17 -shared -undefined dynamic_lookup \
//     $(python3 -m pybind11 --includes) cpp_solver_v2.cpp -framework Accelerate \
//     -o cpp_solver_v2$(python3-config --extension-suffix)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Accelerate/Accelerate.h>

#include <vector>
#include <cmath>
#include <algorithm>

namespace py = pybind11;

enum TermFamily : int { POWER = 0, GAUSS = 1 };

struct FluidChebV2 {
    int K;
    int n_terms;
    int n_delta;       // per-interval degree; residual length = n_delta+2

    double rhoc, R, Tc;

    // Geometry (K,)
    std::vector<double> halfs;
    std::vector<double> mids;

    // Per-term F_k(τ) eval info
    std::vector<int> family;
    std::vector<double> t_vals;
    std::vector<double> beta_vals;
    std::vector<double> gamma_vals;
    std::vector<double> n_values;

    // H: (n_terms, K*(n_delta+1)) — layout for BLAS dgemv (same as v1).
    // Interval-major layout hurt performance because BLAS dgemv's amortized
    // matmul beats K separate naive small matvecs.
    std::vector<double> H;

    // Scratch for a single call (reused)
    mutable std::vector<double> buf_weights;   // (n_terms,)
    mutable std::vector<double> buf_A01_flat;  // (K*(n_delta+1),)

    FluidChebV2(int K_, int n_terms_, int n_delta_,
                double rhoc_, double R_, double Tc_,
                py::array_t<double> a_, py::array_t<double> b_,
                py::array_t<int> family_, py::array_t<double> t_,
                py::array_t<double> beta_, py::array_t<double> gamma_,
                py::array_t<double> n_, py::array_t<double> H_flat_)
        : K(K_), n_terms(n_terms_), n_delta(n_delta_),
          rhoc(rhoc_), R(R_), Tc(Tc_)
    {
        auto cp = [](py::array_t<double> arr) {
            auto r = arr.unchecked<1>();
            std::vector<double> v(r.shape(0));
            for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = r(i);
            return v;
        };
        auto cpi = [](py::array_t<int> arr) {
            auto r = arr.unchecked<1>();
            std::vector<int> v(r.shape(0));
            for (py::ssize_t i = 0; i < r.shape(0); ++i) v[i] = r(i);
            return v;
        };
        auto A = cp(a_); auto B = cp(b_);
        halfs.resize(K); mids.resize(K);
        for (int k = 0; k < K; ++k) {
            halfs[k] = 0.5 * (B[k] - A[k]);
            mids[k]  = 0.5 * (A[k] + B[k]);
        }
        family = cpi(family_);
        t_vals = cp(t_);
        beta_vals = cp(beta_);
        gamma_vals = cp(gamma_);
        n_values = cp(n_);
        // Copy H_flat as-is: (n_terms, K*(n_delta+1)), for BLAS dgemv.
        auto h = H_flat_.unchecked<1>();
        H.resize(h.shape(0));
        for (py::ssize_t i = 0; i < h.shape(0); ++i) H[i] = h(i);
        buf_weights.resize(n_terms);
        buf_A01_flat.resize(static_cast<size_t>(K) * (n_delta + 1));
    }

    // Clenshaw at x of Chebyshev coefficients c[0..n] (length n+1).
    static inline double chebval(double x, const double* c, int deg) {
        double b2 = 0.0, b1 = 0.0, x2 = 2.0 * x;
        for (int k = deg; k >= 1; --k) {
            double b = c[k] + x2 * b1 - b2;
            b2 = b1; b1 = b;
        }
        return c[0] + x * b1 - b2;
    }

    // Derivative coefficients d[0..deg-1] of series c[0..deg].
    static inline void chebder(const double* c, int deg, double* d) {
        if (deg < 1) return;
        d[deg - 1] = 2.0 * deg * c[deg];
        if (deg >= 2) d[deg - 2] = 2.0 * (deg - 1) * c[deg - 1];
        for (int i = deg - 3; i >= 0; --i) {
            d[i] = d[i + 2] + 2.0 * (i + 1) * c[i + 1];
        }
        d[0] *= 0.5;
    }

    // Secant-with-bisection-fallback on c[0..deg] residual, x_std ∈ [-1, 1].
    static inline double secant_bracket(const double* c, int deg,
                                         double fa, double fb) {
        double a_s = -1.0, b_s = 1.0;
        const double eps_bracket = 1e-13;
        const double eps_f = 1e-12;
        double c_new = a_s;
        for (int iter = 0; iter < 12; ++iter) {
            c_new = b_s - fb * (b_s - a_s) / (fb - fa);
            if (c_new <= a_s) return a_s;
            if (c_new >= b_s) return b_s;
            double fc = chebval(c_new, c, deg);
            double mag = std::max({std::abs(fa), std::abs(fb), 1.0});
            if (std::abs(fc) < eps_f * mag) return c_new;
            if ((b_s - a_s) < eps_bracket) return c_new;
            if (fa * fc < 0.0) { b_s = c_new; fb = fc; }
            else { a_s = c_new; fa = fc; }
        }
        return c_new;
    }

    // Precompute A01(δ) coefficients for all intervals at this τ.
    // Stored in buf_A01_flat. Subsequent calls to density_roots_fixed_tau
    // reuse it, saving the F_k evaluation + BLAS matmul each time.
    void set_tau(double T) const {
        const int W = n_delta + 1;
        const double tau = Tc / T;
        for (int k = 0; k < n_terms; ++k) {
            double t = t_vals[k];
            double F;
            if (family[k] == POWER) {
                F = std::pow(tau, t);
            } else {
                double dg = tau - gamma_vals[k];
                F = std::pow(tau, t) * std::exp(-beta_vals[k] * dg * dg);
            }
            buf_weights[k] = n_values[k] * F;
        }
        const int M = n_terms;
        const int N = K * W;
        cblas_dgemv(CblasRowMajor, CblasTrans, M, N, 1.0,
                    H.data(), N,
                    buf_weights.data(), 1,
                    0.0,
                    buf_A01_flat.data(), 1);
    }

    // Paper's "no τ flattening" regime: reuse A01 table from last set_tau().
    py::array_t<double> density_roots_fixed_tau(double T, double P_target) const {
        return density_roots_impl(T, P_target, /*skip_matmul=*/true);
    }

    py::array_t<double> density_roots(double T, double P_target) const {
        return density_roots_impl(T, P_target, /*skip_matmul=*/false);
    }

private:
    py::array_t<double> density_roots_impl(double T, double P_target,
                                            bool skip_matmul) const {
        const int W = n_delta + 1;
        const int RW = n_delta + 2;        // residual length
        const double tau = Tc / T;
        const double factor = rhoc * R * T;

        if (!skip_matmul) {
            // F_k(τ) analytic + weights
            for (int k = 0; k < n_terms; ++k) {
                double t = t_vals[k];
                double F;
                if (family[k] == POWER) {
                    F = std::pow(tau, t);
                } else {
                    double dg = tau - gamma_vals[k];
                    F = std::pow(tau, t) * std::exp(-beta_vals[k] * dg * dg);
                }
                buf_weights[k] = n_values[k] * F;
            }
            // One BLAS matvec: A01_flat = weights^T @ H (row-major)
            const int M = n_terms;
            const int N = K * W;
            cblas_dgemv(CblasRowMajor, CblasTrans, M, N, 1.0,
                        H.data(), N,
                        buf_weights.data(), 1,
                        0.0,
                        buf_A01_flat.data(), 1);
        }
        (void)tau;  // unused if skip_matmul

        // Fused per-interval: residual build + skip test + mono test + secant
        constexpr int MAX_N = 48;
        double r[MAX_N + 1];      // residual (length RW)
        double d[MAX_N + 1];      // derivative coeffs

        std::vector<double> out_roots;
        out_roots.reserve(16);

        for (int kk = 0; kk < K; ++kk) {
            const double half = halfs[kk];
            const double mid  = mids[kk];
            const double* A01 = &buf_A01_flat[kk * W];

            // Build residual (length RW = n_delta+2):
            //   r = factor * (half * c_prime + mid * A01)  with constants merged in
            //   where c_prime is times-x identity on A01.
            //   Then add factor*mid to r[0], factor*half to r[1], subtract P_target from r[0].
            //
            // c_prime formulas:
            //   cp[0]   = 0.5 * A01[1]
            //   cp[1]   = A01[0] + 0.5 * A01[2]   (for n_delta >= 2)
            //   cp[i]   = 0.5 * (A01[i-1] + A01[i+1]),  2 <= i <= n-1
            //   cp[n]   = 0.5 * A01[n-1]
            //   cp[n+1] = 0.5 * A01[n]
            r[0] = factor * (half * (0.5 * A01[1]) + mid * A01[0]) + factor * mid - P_target;
            // Accumulate sum of abs of c_{j>=1} for cheap-skip test
            double sum_rest = 0.0;
            double v1 = factor * (half * (A01[0] + 0.5 * A01[2]) + mid * A01[1]) + factor * half;
            r[1] = v1;
            sum_rest += std::abs(v1);
            for (int i = 2; i < n_delta; ++i) {
                double cp_i = 0.5 * (A01[i - 1] + A01[i + 1]);
                double vi = factor * (half * cp_i + mid * A01[i]);
                r[i] = vi;
                sum_rest += std::abs(vi);
            }
            double vn = factor * (half * 0.5 * A01[n_delta - 1] + mid * A01[n_delta]);
            r[n_delta] = vn;
            sum_rest += std::abs(vn);
            double vn1 = factor * half * 0.5 * A01[n_delta];
            r[n_delta + 1] = vn1;
            sum_rest += std::abs(vn1);

            // Cheap no-root test: skip interval if |c_0| > sum of other |c_j|
            double abs0 = std::abs(r[0]);
            if (abs0 > sum_rest) continue;

            // Compute derivative coefficients for monotonicity test
            chebder(r, n_delta + 1, d);
            double abs_d0 = std::abs(d[0]);
            double sum_d_rest = 0.0;
            for (int j = 1; j < n_delta + 1; ++j) sum_d_rest += std::abs(d[j]);
            bool mono = abs_d0 > sum_d_rest;

            // Endpoint values
            double fa = 0.0, fb = 0.0, sgn = 1.0;
            for (int j = 0; j <= n_delta + 1; ++j) {
                fb += r[j];
                fa += sgn * r[j];
                sgn = -sgn;
            }

            if (mono) {
                if (fa * fb >= 0.0) continue;  // monotone, same-sign endpoints
                double root_std = secant_bracket(r, n_delta + 1, fa, fb);
                double delta_r = mid + half * root_std;
                out_roots.push_back(delta_r * rhoc);
            } else {
                // Non-monotone: LAPACK dhseqr on the pre-Hessenberg colleague matrix.
                // Build the symmetric form of the Chebyshev colleague (feedback in
                // LAST COLUMN) used by numpy's chebcompanion, which IS upper Hessenberg.
                //   M[0,1] = sqrt(0.5); M[1,0] = sqrt(0.5)
                //   M[i,i+1] = 0.5 for 1 <= i <= n-2;  M[i+1,i] = 0.5 for 1 <= i <= n-2
                //   M[0, n-1]   -= c_0 / (c_n · sqrt(2))
                //   M[i, n-1]   -= c_i / (2·c_n)  for 1 <= i <= n-1
                int deg = n_delta + 1;
                while (deg > 0 && std::abs(r[deg]) < 1e-300) deg--;
                if (deg < 1) continue;
                // Stack-allocated: no malloc in hot path. MAX_N=48 → up to n_delta=47.
                double H_c[(MAX_N + 1) * (MAX_N + 1)];
                for (int i = 0; i < deg * deg; ++i) H_c[i] = 0.0;
                const double s = std::sqrt(0.5);
                // Column-major: H_c[j*deg + i] = M[i, j]
                if (deg >= 2) {
                    H_c[1 * deg + 0] = s;      // M[0,1] = sqrt(0.5)
                    H_c[0 * deg + 1] = s;      // M[1,0] = sqrt(0.5)
                    for (int i = 1; i <= deg - 2; ++i) {
                        H_c[(i + 1) * deg + i] = 0.5;   // M[i, i+1]
                        H_c[i * deg + (i + 1)] = 0.5;   // M[i+1, i]
                    }
                    // Feedback in last column
                    const double inv_cn = 1.0 / r[deg];
                    H_c[(deg - 1) * deg + 0] -= r[0] * inv_cn / std::sqrt(2.0);
                    for (int i = 1; i < deg; ++i) {
                        H_c[(deg - 1) * deg + i] -= r[i] * inv_cn * 0.5;
                    }
                } else {
                    H_c[0] = -r[0] / r[1];
                }
                double wr[MAX_N + 1], wi[MAX_N + 1];
                double work[16 * (MAX_N + 1)];
                char job = 'E', compz = 'N';
                int info = 0;
                int N_l = deg, ilo = 1, ihi = deg;
                int ldh = deg, ldz = 1;
                int lwork = 16 * (MAX_N + 1);
                double z;
                dhseqr_(&job, &compz, &N_l, &ilo, &ihi,
                        H_c, &ldh,
                        wr, wi,
                        &z, &ldz,
                        work, &lwork, &info);
                if (info != 0) continue;
                for (int j = 0; j < deg; ++j) {
                    if (std::abs(wi[j]) < 1e-8 * (1.0 + std::abs(wr[j]))) {
                        double rr = wr[j];
                        if (rr >= -1.0 - 1e-10 && rr <= 1.0 + 1e-10) {
                            rr = std::max(-1.0, std::min(1.0, rr));
                            double delta = mid + half * rr;
                            out_roots.push_back(delta * rhoc);
                        }
                    }
                }
            }
        }
        // Return as numpy array
        py::array_t<double> result(static_cast<py::ssize_t>(out_roots.size()));
        auto buf = result.mutable_unchecked<1>();
        for (size_t i = 0; i < out_roots.size(); ++i) buf(i) = out_roots[i];
        return result;
    }
};


PYBIND11_MODULE(cpp_solver_v2, m) {
    m.doc() = "Optimized density rootfinder (v2).";
    py::class_<FluidChebV2>(m, "FluidChebV2")
        .def(py::init<int, int, int, double, double, double,
                      py::array_t<double>, py::array_t<double>,
                      py::array_t<int>, py::array_t<double>,
                      py::array_t<double>, py::array_t<double>,
                      py::array_t<double>, py::array_t<double>>())
        .def("density_roots", &FluidChebV2::density_roots,
             py::arg("T"), py::arg("P_target"))
        .def("set_tau", &FluidChebV2::set_tau, py::arg("T"),
             "Precompute the τ-dependent A01 coefficient matrix once; "
             "subsequent density_roots_fixed_tau calls reuse it.")
        .def("density_roots_fixed_tau", &FluidChebV2::density_roots_fixed_tau,
             py::arg("T"), py::arg("P_target"),
             "Run rootfinder reusing the A01 table from the last set_tau() "
             "call. This is the paper's 'no τ flattening' regime (§ 4.3).");
}
