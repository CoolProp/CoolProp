// C++ port of the hot path from fast_solver.py. Benchmarks the single-digit-μs
// target claimed by the 2018 Bell & Alpert paper for small-component mixtures.
//
// Architecture mirrors the Python code:
//   - H_flat:     (n_terms, K*(n_delta+1)) flattened coefficient matrix
//   - intervals:  per-piece [a, b] in δ
//   - At runtime: compute F_k(τ), matmul with H_flat via BLAS, add analytic
//     times-x, cheap-skip + monotonicity test, inlined secant on mono+bracket
//     intervals. Non-monotone intervals call LAPACK dgeev on the colleague
//     matrix (fallback).
//
// Compile (macOS, Accelerate BLAS/LAPACK):
//   c++ -O3 -ffast-math -march=native -std=c++17 -shared -undefined dynamic_lookup \
//     $(python3 -m pybind11 --includes) cpp_solver.cpp -framework Accelerate \
//     -o cpp_solver$(python3-config --extension-suffix)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Accelerate/Accelerate.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>

namespace py = pybind11;

enum TermFamily : int { POWER = 0, GAUSS = 1 };

struct FluidChebCpp {
    // Dimensions
    int K;
    int n_terms;
    int n_delta;  // per-interval degree; coeff-vector length = n_delta+1

    // Constants
    double rhoc, R, Tc;

    // Geometry
    std::vector<double> intervals_a;   // (K,)
    std::vector<double> intervals_b;   // (K,)

    // Per-term info for F_k(τ) closed-form evaluation
    std::vector<int> family;           // (n_terms,)  POWER or GAUSS
    std::vector<double> t_vals;        // (n_terms,)  τ exponent
    std::vector<double> beta_vals;     // (n_terms,)  Gaussian τ damping (0 for POWER)
    std::vector<double> gamma_vals;    // (n_terms,)  Gaussian τ center (ignored for POWER)
    std::vector<double> n_values;      // (n_terms,)

    // H_flat: layout (n_terms, K*(n_delta+1)), row-major
    std::vector<double> H_flat;

    // Scratch buffers reused per call to avoid per-call allocations
    mutable std::vector<double> buf_weights;   // (n_terms,)
    mutable std::vector<double> buf_A01_flat;  // (K*(n_delta+1),)
    mutable std::vector<double> buf_resid;     // (K*(n_delta+2),)  residual per interval (longer by 1)

    FluidChebCpp(int K_, int n_terms_, int n_delta_,
                 double rhoc_, double R_, double Tc_,
                 py::array_t<double> a_, py::array_t<double> b_,
                 py::array_t<int> family_, py::array_t<double> t_,
                 py::array_t<double> beta_, py::array_t<double> gamma_,
                 py::array_t<double> n_, py::array_t<double> H_flat_)
        : K(K_), n_terms(n_terms_), n_delta(n_delta_),
          rhoc(rhoc_), R(R_), Tc(Tc_) {
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
        intervals_a = cp(a_);
        intervals_b = cp(b_);
        family = cpi(family_);
        t_vals = cp(t_);
        beta_vals = cp(beta_);
        gamma_vals = cp(gamma_);
        n_values = cp(n_);
        auto h = H_flat_.unchecked<1>();
        H_flat.resize(static_cast<size_t>(n_terms) *
                      static_cast<size_t>(K) *
                      static_cast<size_t>(n_delta + 1));
        for (size_t i = 0; i < H_flat.size(); ++i) H_flat[i] = h(i);
        buf_weights.resize(n_terms);
        buf_A01_flat.resize(static_cast<size_t>(K) * (n_delta + 1));
        buf_resid.resize(static_cast<size_t>(K) * (n_delta + 2));
    }

    // Evaluate F_k(τ) for all terms (analytic, scalar τ).
    void eval_F_into_weights(double tau) const {
        for (int k = 0; k < n_terms; ++k) {
            double F;
            double t = t_vals[k];
            if (family[k] == POWER) {
                F = std::pow(tau, t);
            } else {  // GAUSS
                double dg = tau - gamma_vals[k];
                F = std::pow(tau, t) * std::exp(-beta_vals[k] * dg * dg);
            }
            buf_weights[k] = n_values[k] * F;
        }
    }

    // Chebyshev value at x for coefficients c[0..n] (Clenshaw).
    static inline double chebval(double x, const double* c, int n) {
        if (n == 0) return c[0];
        double b2 = 0.0, b1 = 0.0, x2 = 2.0 * x;
        for (int k = n; k >= 1; --k) {
            double b = c[k] + x2 * b1 - b2;
            b2 = b1; b1 = b;
        }
        return c[0] + x * b1 - b2;
    }

    // Derivative coefficients d[0..n-1] of Chebyshev series c[0..n].
    static void chebder(const double* c, int n, double* d) {
        if (n < 1) return;
        d[n - 1] = 2.0 * n * c[n];
        if (n >= 2) d[n - 2] = 2.0 * (n - 1) * c[n - 1];
        for (int i = n - 3; i >= 0; --i) {
            d[i] = d[i + 2] + 2.0 * (i + 1) * c[i + 1];
        }
        d[0] *= 0.5;
    }

    // Inline secant-with-bisection-fallback for bracketed monotone case in x_std ∈ [-1, 1].
    // c is length n+1 (residual coefficients on the interval).
    static double secant_bracket(const double* c, int n, double fa, double fb) {
        double a_s = -1.0, b_s = 1.0;
        const double eps_bracket = 1e-13;
        const double eps_f = 1e-12;
        double c_new = a_s;
        for (int iter = 0; iter < 12; ++iter) {
            c_new = b_s - fb * (b_s - a_s) / (fb - fa);
            if (c_new <= a_s) return a_s;
            if (c_new >= b_s) return b_s;
            double fc = chebval(c_new, c, n);
            double mag = std::max({std::abs(fa), std::abs(fb), 1.0});
            if (std::abs(fc) < eps_f * mag) return c_new;
            if ((b_s - a_s) < eps_bracket) return c_new;
            if (fa * fc < 0.0) { b_s = c_new; fb = fc; }
            else { a_s = c_new; fa = fc; }
        }
        return c_new;
    }

    // Main entry: solve P(T, ρ) = P_target; return all density roots in mol/m^3.
    py::list density_roots(double T, double P_target) const {
        double tau = Tc / T;
        eval_F_into_weights(tau);

        // BLAS matvec: A01_flat[j] = Σ_k weights[k] * H_flat[k, j]
        // H_flat is (n_terms, K*(n_delta+1)) row-major; with CBlasRowMajor that's
        // H viewed as (rows=n_terms, cols=K*(n_delta+1)); we want H^T @ w = res
        // which cblas_dgemv handles as CblasTrans if we pass the matrix as-is.
        const int M = n_terms;
        const int N = K * (n_delta + 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, M, N, 1.0,
                    H_flat.data(), N,
                    buf_weights.data(), 1,
                    0.0,
                    buf_A01_flat.data(), 1);

        // Build the pressure residual per interval into buf_resid[k, 0..n_delta+1]
        // (length n_delta+2). Residual = ρ_r R T [δ + δ·A01] − P_target as a single
        // Chebyshev series in x_std on [a_k, b_k]. The "δ" part contributes
        //   ρ_r R T (mid + half·x_std) which adds mid·(ρ_r R T) to c_0 and
        //   half·(ρ_r R T) to c_1.
        // The "δ·A01" part is computed via the multiply-by-x identity on the
        // interval: c_new[0] = 0.5 a_1; c_new[1] = a_0 + 0.5 a_2; c_new[i] =
        // 0.5 (a_{i-1}+a_{i+1}); c_new[n] = 0.5 a_{n-1}; c_new[n+1] = 0.5 a_n.
        // Then multiply by half and add mid * padded A01.
        const int n = n_delta;
        const double factor = rhoc * R * T;
        for (int kk = 0; kk < K; ++kk) {
            const double half = 0.5 * (intervals_b[kk] - intervals_a[kk]);
            const double mid  = 0.5 * (intervals_b[kk] + intervals_a[kk]);
            const double* A01 = &buf_A01_flat[kk * (n + 1)];
            double* r = &buf_resid[kk * (n + 2)];
            // times_x result, then scale by half, add mid*A01, then by factor.
            // c_prime lives implicitly in r[0..n+1].
            double c0 = 0.5 * A01[1];
            double c1 = (n >= 2) ? (A01[0] + 0.5 * A01[2]) : A01[0];
            r[0] = factor * (half * c0 + mid * A01[0]);
            r[1] = factor * (half * c1 + ((n + 1 > 1) ? mid * A01[1] : 0.0));
            for (int i = 2; i < n; ++i) {
                double cp_i = 0.5 * (A01[i - 1] + A01[i + 1]);
                r[i] = factor * (half * cp_i + mid * A01[i]);
            }
            if (n >= 1) {
                double cp_n = 0.5 * A01[n - 1];
                r[n] = factor * (half * cp_n + mid * A01[n]);
                double cp_np1 = 0.5 * A01[n];
                r[n + 1] = factor * half * cp_np1;
            }
            // Add ρ_r R T · δ contribution: factor*mid to c0 and factor*half to c1
            r[0] += factor * mid;
            r[1] += factor * half;
            // Subtract P_target (constant goes into c0)
            r[0] -= P_target;
        }

        py::list roots;
        std::vector<double> d_tmp(n + 1);
        for (int kk = 0; kk < K; ++kk) {
            const double* r = &buf_resid[kk * (n + 2)];
            // Cheap no-root test: |c_0| > Σ |c_{j>=1}|
            double abs0 = std::abs(r[0]);
            double sum_rest = 0.0;
            for (int j = 1; j <= n + 1; ++j) sum_rest += std::abs(r[j]);
            if (abs0 > sum_rest) continue;

            // Monotonicity test on the derivative coefficients of r (length n+1)
            chebder(r, n + 1, d_tmp.data());
            double abs_d0 = std::abs(d_tmp[0]);
            double sum_d_rest = 0.0;
            for (int j = 1; j < n + 1; ++j) sum_d_rest += std::abs(d_tmp[j]);
            bool mono = abs_d0 > sum_d_rest;

            // Endpoint values
            double fa = 0.0, fb = 0.0;
            double sgn = 1.0;
            for (int j = 0; j <= n + 1; ++j) {
                fb += r[j];
                fa += sgn * r[j];
                sgn = -sgn;
            }
            if (mono) {
                if (fa * fb >= 0.0) continue;  // monotone, no sign change → no root
                double root_std = secant_bracket(r, n + 1, fa, fb);
                double half = 0.5 * (intervals_b[kk] - intervals_a[kk]);
                double mid  = 0.5 * (intervals_b[kk] + intervals_a[kk]);
                double delta = mid + half * root_std;
                roots.append(delta * rhoc);
            } else {
                // Non-monotone: fall back to LAPACK eig on colleague matrix.
                // Build a colleague matrix of size (n+1)×(n+1) per numpy convention.
                // For simplicity use dgeev (LAPACK), finding all eigenvalues.
                // (n+1 is degree, so matrix is (n+1) x (n+1))
                const int deg = n + 1;
                std::vector<double> coeffs(deg + 1);
                for (int j = 0; j <= deg; ++j) coeffs[j] = r[j];
                // Drop trailing near-zeros to get polynomial of effective degree
                int eff_deg = deg;
                while (eff_deg > 0 && std::abs(coeffs[eff_deg]) < 1e-300) eff_deg--;
                if (eff_deg < 1) continue;
                // Colleague matrix for Chebyshev-1st-kind: eff_deg × eff_deg, upper Hessenberg
                std::vector<double> M(eff_deg * eff_deg, 0.0);
                if (eff_deg >= 2) {
                    M[0 * eff_deg + 1] = 1.0;
                    for (int i = 1; i < eff_deg - 1; ++i) {
                        M[i * eff_deg + (i - 1)] = 0.5;
                        M[i * eff_deg + (i + 1)] = 0.5;
                    }
                    M[(eff_deg - 1) * eff_deg + (eff_deg - 2)] = 0.5;
                }
                // Last row: coefficient feedback
                // a_j = -c_j / (2 c_n);  correction on last col
                if (eff_deg >= 2) {
                    for (int j = 0; j < eff_deg; ++j)
                        M[(eff_deg - 1) * eff_deg + j] = -coeffs[j] / (2.0 * coeffs[eff_deg]);
                    M[(eff_deg - 1) * eff_deg + (eff_deg - 1)] += 0.5;
                } else {
                    M[0] = -coeffs[0] / coeffs[1];
                }
                // Call dgeev (column-major). Our M is row-major; transpose for column-major.
                std::vector<double> Mt(eff_deg * eff_deg);
                for (int i = 0; i < eff_deg; ++i)
                    for (int j = 0; j < eff_deg; ++j)
                        Mt[j * eff_deg + i] = M[i * eff_deg + j];
                std::vector<double> wr(eff_deg), wi(eff_deg);
                char jobvl = 'N', jobvr = 'N';
                int info = 0;
                int lwork = 4 * eff_deg;
                std::vector<double> work(lwork);
                int lda = eff_deg;
                int ldvl = 1, ldvr = 1;
                double vl, vr;
                int N_l = eff_deg;
                dgeev_(&jobvl, &jobvr, &N_l, Mt.data(), &lda,
                       wr.data(), wi.data(), &vl, &ldvl, &vr, &ldvr,
                       work.data(), &lwork, &info);
                if (info != 0) continue;
                for (int j = 0; j < eff_deg; ++j) {
                    if (std::abs(wi[j]) < 1e-8 * (1.0 + std::abs(wr[j]))) {
                        double rr = wr[j];
                        if (rr >= -1.0 - 1e-10 && rr <= 1.0 + 1e-10) {
                            rr = std::max(-1.0, std::min(1.0, rr));
                            double half = 0.5 * (intervals_b[kk] - intervals_a[kk]);
                            double mid  = 0.5 * (intervals_b[kk] + intervals_a[kk]);
                            double delta = mid + half * rr;
                            roots.append(delta * rhoc);
                        }
                    }
                }
            }
        }
        return roots;
    }
};


PYBIND11_MODULE(cpp_solver, m) {
    m.doc() = "C++ density rootfinder, mirrors fast_solver.density_roots_vectorized";
    py::class_<FluidChebCpp>(m, "FluidChebCpp")
        .def(py::init<int, int, int, double, double, double,
                      py::array_t<double>, py::array_t<double>,
                      py::array_t<int>, py::array_t<double>,
                      py::array_t<double>, py::array_t<double>,
                      py::array_t<double>, py::array_t<double>>())
        .def("density_roots", &FluidChebCpp::density_roots,
             py::arg("T"), py::arg("P_target"));
}
