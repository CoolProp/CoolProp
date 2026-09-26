// Whole-triangle derivatives for N = 1..4, per method: values + min-of-repeats timings.
// Also times two "only what you need" baselines from the first experiment.
#include <chrono>
#include <cstdio>
#include <vector>
#include "api_all.hpp"

using namespace spike;
MethodAll methodall_taylor2(int), methodall_polar(int), methodall_teqp(int);
MethodTable method_double(int), method_numdual(int);

template <class F>
double time_ns(F&& f) {
    const int reps = 7, iters = 5000;
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        for (int it = 0; it < iters; ++it)
            f(it);
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / iters);
    }
    return best;
}

int main() {
    PCSAFTParams propane{1, {2.0020}, {3.6184}, {208.11}};
    PCSAFTParams c123{3, {1.0000, 1.6069, 2.0020}, {3.7039, 3.5206, 3.6184}, {150.03, 191.42, 208.11}};
    struct S
    {
        const char* name;
        PCSAFTParams p;
        double T, rho;
        XArr x;
    };
    const std::vector<S> states = {
      {"propane_liq", propane, 300.0, 11000.0, {1.0}},
      {"propane_gas", propane, 300.0, 100.0, {1.0}},
      {"C1C2C3_dense", c123, 250.0, 8000.0, {0.5, 0.3, 0.2}},
    };
    std::printf("kind,method,N,state,i,j,value\n");
    volatile double sink = 0;
    for (auto mk : {methodall_taylor2, methodall_polar, methodall_teqp}) {
        const MethodAll m = mk(0);
        for (int N = 1; N <= NMAXORDER; ++N) {
            if (!m.byN[N]) continue;
            for (const auto& s : states) {
                double o[NSLOTS];
                m.byN[N](s.p, s.T, s.rho, s.x, o);
                for (int k = 0; k <= N; ++k)
                    for (int j = 0; j <= k; ++j)
                        std::printf("val,%s,%d,%s,%d,%d,%.17g\n", m.name, N, s.name, k - j, j, o[idx2(k - j, j)]);
                const double ns = time_ns([&](int it) {
                    m.byN[N](s.p, s.T, s.rho * (1.0 + 1e-12 * (it & 7)), s.x, o);
                    sink = sink + o[0];
                });
                std::printf("ns,%s,%d,%s,0,0,%.6g\n", m.name, N, s.name, ns);
            }
        }
    }
    // baselines: one double evaluation; numdual Dual3 in rho (A00..A03); numdual HyperDual A11
    for (const auto& s : states) {
        double o[4];
        const MethodTable d = method_double(0), n = method_numdual(0);
        std::printf("ns,value_only,0,%s,0,0,%.6g\n", s.name, time_ns([&](int it) {
                        d.Ar0n(s.p, s.T, s.rho * (1.0 + 1e-12 * (it & 7)), s.x, o);
                        sink = sink + o[0];
                    }));
        std::printf("ns,numdual_Ar0n,3,%s,0,0,%.6g\n", s.name, time_ns([&](int it) {
                        n.Ar0n(s.p, s.T, s.rho * (1.0 + 1e-12 * (it & 7)), s.x, o);
                        sink = sink + o[0];
                    }));
        std::printf("ns,numdual_Ar11,2,%s,0,0,%.6g\n", s.name, time_ns([&](int it) {
                        n.Ar11(s.p, s.T, s.rho * (1.0 + 1e-12 * (it & 7)), s.x, o);
                        sink = sink + o[0];
                    }));
    }
}
