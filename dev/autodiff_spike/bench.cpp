// Runs every method on a few states, prints values (for the accuracy check against the
// mpmath reference) and min-of-repeats timings.  Methods live in separate TUs and are
// reached through function pointers, so nothing is hoisted across calls.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>
#include "api.hpp"

using namespace spike;
MethodTable method_double(int), method_ad_real(int), method_ad_dual(int), method_numdual(int), method_mcx(int), method_cstep(int), method_fd(int);

struct State
{
    const char* name;
    PCSAFTParams p;
    double T, rho;
    XArr x;
};

int main() {
    // Gross & Sadowski 2001 parameters
    PCSAFTParams propane{1, {2.0020}, {3.6184}, {208.11}};
    PCSAFTParams c123{3, {1.0000, 1.6069, 2.0020}, {3.7039, 3.5206, 3.6184}, {150.03, 191.42, 208.11}};
    const std::vector<State> states = {
      {"propane_liq", propane, 300.0, 11000.0, {1.0}},
      {"propane_gas", propane, 300.0, 100.0, {1.0}},
      {"C1C2C3_dense", c123, 250.0, 8000.0, {0.5, 0.3, 0.2}},
    };
    const std::vector<MethodTable> methods = {method_double(0), method_ad_real(0), method_ad_dual(0), method_numdual(0),
                                              method_mcx(0),    method_cstep(0),   method_fd(0)};
    const char* qnames[4] = {"Ar0n", "Arn0", "Ar11", "gradx"};

    std::printf("kind,method,state,quantity,index,value\n");
    for (const auto& m : methods) {
        const Fn fns[4] = {m.Ar0n, m.Arn0, m.Ar11, m.gradx};
        for (const auto& s : states) {
            for (int q = 0; q < 4; ++q) {
                double out[4];
                fns[q](s.p, s.T, s.rho, s.x, out);
                const int n = (q == 0) ? 4 : (q == 1) ? 3 : (q == 2) ? 1 : int(s.p.N);
                for (int i = 0; i < n; ++i)
                    std::printf("val,%s,%s,%s,%d,%.17g\n", m.name, s.name, qnames[q], i, out[i]);

                // timing: min over repeats of the mean per call
                const int reps = 7, iters = 20000;
                double best = 1e300, sink = 0;
                for (int r = 0; r < reps; ++r) {
                    const auto t0 = std::chrono::steady_clock::now();
                    for (int it = 0; it < iters; ++it) {
                        fns[q](s.p, s.T, s.rho * (1.0 + 1e-12 * (it & 7)), s.x, out);
                        sink += out[0];
                    }
                    const auto t1 = std::chrono::steady_clock::now();
                    best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / iters);
                }
                std::printf("ns,%s,%s,%s,%d,%.6g\n", m.name, s.name, qnames[q], (sink == 12345.0), best);
            }
        }
    }
}
