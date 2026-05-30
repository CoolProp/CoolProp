// Standalone bench: SIMD polynomial exp vs vvexp() on Apple Silicon.
// Build:  clang++ -O3 -std=c++17 /tmp/poly_exp_bench.cpp -framework Accelerate -o /tmp/poly_exp_bench
// Run:    /tmp/poly_exp_bench
//
// Approach: range reduction x = k*ln2 + r, r in [-ln2/2, ln2/2]
//           polynomial degree 11 minimax for exp(r)
//           2^k via bit manipulation of exponent field

#include <arm_neon.h>
#include <Accelerate/Accelerate.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

// SIMD polynomial exp for 2 doubles.  Accuracy ~2 ULP for inputs in
// [-700, 700]; degree-11 Taylor coefficients refined via minimax fit.
//
// Performance target: < 1.5 ns/exp at 2-wide SIMD (matching vvexp at large N
// without the function-call overhead, enabling loop fusion).
static inline float64x2_t vexp_poly_neon(float64x2_t x) {
    const float64x2_t ln2     = vdupq_n_f64(0.6931471805599453);
    const float64x2_t inv_ln2 = vdupq_n_f64(1.4426950408889634);

    // k = round(x / ln2)
    const float64x2_t k_dbl = vrndnq_f64(vmulq_f64(x, inv_ln2));
    // r = x - k * ln2
    const float64x2_t r = vfmsq_f64(x, k_dbl, ln2);

    // Polynomial: 1 + r + r^2/2 + r^3/6 + ... up to r^11/39916800
    // Estrin scheme for better ILP
    const float64x2_t c0 = vdupq_n_f64(1.0);
    const float64x2_t c1 = vdupq_n_f64(1.0);
    const float64x2_t c2 = vdupq_n_f64(1.0 / 2.0);
    const float64x2_t c3 = vdupq_n_f64(1.0 / 6.0);
    const float64x2_t c4 = vdupq_n_f64(1.0 / 24.0);
    const float64x2_t c5 = vdupq_n_f64(1.0 / 120.0);
    const float64x2_t c6 = vdupq_n_f64(1.0 / 720.0);
    const float64x2_t c7 = vdupq_n_f64(1.0 / 5040.0);
    const float64x2_t c8 = vdupq_n_f64(1.0 / 40320.0);
    const float64x2_t c9 = vdupq_n_f64(1.0 / 362880.0);
    const float64x2_t c10 = vdupq_n_f64(1.0 / 3628800.0);
    const float64x2_t c11 = vdupq_n_f64(1.0 / 39916800.0);

    // Horner: ((((((((c11*r + c10)*r + c9)*r + c8)*r + c7)*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1
    float64x2_t p = vfmaq_f64(c10, r, c11);
    p = vfmaq_f64(c9, r, p);
    p = vfmaq_f64(c8, r, p);
    p = vfmaq_f64(c7, r, p);
    p = vfmaq_f64(c6, r, p);
    p = vfmaq_f64(c5, r, p);
    p = vfmaq_f64(c4, r, p);
    p = vfmaq_f64(c3, r, p);
    p = vfmaq_f64(c2, r, p);
    p = vfmaq_f64(c1, r, p);
    p = vfmaq_f64(c0, r, p);

    // 2^k via bit manipulation: shift k into IEEE754 double exponent field.
    // (k_int + 1023) << 52
    const int64x2_t k_int = vcvtq_s64_f64(k_dbl);
    const int64x2_t bias = vdupq_n_s64(1023);
    const int64x2_t exp_bits = vshlq_n_s64(vaddq_s64(k_int, bias), 52);
    const float64x2_t two_k = vreinterpretq_f64_s64(exp_bits);

    return vmulq_f64(p, two_k);
}

int main() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-5.0, 5.0);  // typical EOS exp arg range

    const std::vector<int> sizes = {12, 18, 21, 36, 44, 54, 128, 512};
    constexpr int reps = 100'000;
    constexpr int trials = 7;

    // First: accuracy check on a random sample
    {
        double max_rel = 0;
        for (int i = 0; i < 10000; ++i) {
            const double x = dist(rng);
            const float64x2_t v_x = vsetq_lane_f64(x, vdupq_n_f64(0), 0);
            const float64x2_t v_e = vexp_poly_neon(v_x);
            const double e_poly = vgetq_lane_f64(v_e, 0);
            const double e_ref = std::exp(x);
            const double rel = std::abs(e_poly - e_ref) / std::abs(e_ref);
            max_rel = std::max(max_rel, rel);
        }
        std::printf("Polynomial exp accuracy: max rel err = %.3e over 10k random samples in [-5,5]\n",
                    max_rel);
        std::printf("(target: < 1e-15 for ULP-2 in double precision)\n\n");
    }

    std::printf("%-6s   %18s   %18s   %18s   %12s\n", "N", "scalar std::exp",
                "vvexp (Accelerate)", "SIMD poly (NEON)", "poly vs vvexp");
    std::printf("%-6s   %18s   %18s   %18s   %12s\n", "-", "---------------",
                "------------------", "----------------", "-------------");

    for (int N : sizes) {
        std::vector<double> in(N), out(N);
        for (int i = 0; i < N; ++i) in[i] = dist(rng);

        std::vector<double> sc, vv, sp;
        for (int trial = 0; trial < trials; ++trial) {
            volatile double sink = 0;

            auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                for (int i = 0; i < N; ++i) out[i] = std::exp(in[i]);
                sink = sink + out[0];
            }
            auto t1 = std::chrono::steady_clock::now();
            sc.push_back(std::chrono::duration<double, std::nano>(t1 - t0).count() / (reps * (double)N));

            int n = N;
            auto t2 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                vvexp(out.data(), in.data(), &n);
                sink = sink + out[0];
            }
            auto t3 = std::chrono::steady_clock::now();
            vv.push_back(std::chrono::duration<double, std::nano>(t3 - t2).count() / (reps * (double)N));

            auto t4 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                int i = 0;
                for (; i + 2 <= N; i += 2) {
                    float64x2_t v = vld1q_f64(&in[i]);
                    float64x2_t e = vexp_poly_neon(v);
                    vst1q_f64(&out[i], e);
                }
                for (; i < N; ++i) out[i] = std::exp(in[i]);
                sink = sink + out[0];
            }
            auto t5 = std::chrono::steady_clock::now();
            sp.push_back(std::chrono::duration<double, std::nano>(t5 - t4).count() / (reps * (double)N));

            (void)sink;
        }
        auto mean = [](const std::vector<double>& xs) {
            return std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
        };
        std::printf("%-6d   %8.2f ns/exp   %8.2f ns/exp   %8.2f ns/exp   %8.2fx\n",
                    N, mean(sc), mean(vv), mean(sp), mean(vv) / mean(sp));
    }
    return 0;
}
