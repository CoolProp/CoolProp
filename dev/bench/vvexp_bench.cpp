// Standalone: Apple Accelerate vvexp() vs scalar std::exp() loop.
// Build:  clang++ -O3 -std=c++17 /tmp/vvexp_bench.cpp -framework Accelerate -o /tmp/vvexp_bench
// Run:    /tmp/vvexp_bench
//
// Measures per-exp wall time for arrays of varying size, comparable to the
// CoolProp residual-Helmholtz term count (12-54 typical).

#include <Accelerate/Accelerate.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <vector>

int main() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-2.0, 2.0);

    const std::vector<int> sizes = {12, 18, 21, 36, 44, 54, 128, 512, 4096};
    constexpr int reps = 200'000;
    constexpr int trials = 5;

    std::printf("%-6s   scalar (ns/exp)        vvexp (ns/exp)     scalar/vvexp\n", "N");
    std::printf("%-6s   ---------------        --------------     ------------\n", "-");

    for (int N : sizes) {
        std::vector<double> in(N), out(N);
        for (int i = 0; i < N; ++i) in[i] = dist(rng);

        std::vector<double> scalar_ns_per_exp, vvexp_ns_per_exp;

        for (int trial = 0; trial < trials; ++trial) {
            // Scalar: looped std::exp
            volatile double scalar_sink = 0;
            auto t0 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                for (int i = 0; i < N; ++i) out[i] = std::exp(in[i]);
                scalar_sink = scalar_sink + out[0];
            }
            auto t1 = std::chrono::steady_clock::now();
            scalar_ns_per_exp.push_back(
                std::chrono::duration<double, std::nano>(t1 - t0).count() / (reps * (double)N));
            (void)scalar_sink;

            // vvexp from Accelerate
            volatile double vvexp_sink = 0;
            int n = N;
            auto t2 = std::chrono::steady_clock::now();
            for (int r = 0; r < reps; ++r) {
                vvexp(out.data(), in.data(), &n);
                vvexp_sink = vvexp_sink + out[0];
            }
            auto t3 = std::chrono::steady_clock::now();
            vvexp_ns_per_exp.push_back(
                std::chrono::duration<double, std::nano>(t3 - t2).count() / (reps * (double)N));
            (void)vvexp_sink;
        }

        auto mean = [](const std::vector<double>& xs) {
            return std::accumulate(xs.begin(), xs.end(), 0.0) / xs.size();
        };
        auto std_ = [&mean](const std::vector<double>& xs) {
            double m = mean(xs), v = 0;
            for (double x : xs) v += (x - m) * (x - m);
            return std::sqrt(v / xs.size());
        };
        double sm = mean(scalar_ns_per_exp), ss = std_(scalar_ns_per_exp);
        double vm = mean(vvexp_ns_per_exp), vs = std_(vvexp_ns_per_exp);
        std::printf("%-6d   %6.2f ± %5.2f         %6.2f ± %5.2f      %5.2f×\n",
                    N, sm, ss, vm, vs, sm / vm);
    }
    return 0;
}
