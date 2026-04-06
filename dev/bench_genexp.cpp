// Benchmark ResidualHelmholtzGeneralizedExponential::all() directly.
//
// Build (from repo root):
//   cmake build -DCOOLPROP_MY_MAIN=../dev/bench_genexp.cpp && cmake --build build --target Main -j
// Run:
//   ./build/Main
//
// Tests the GenExp.all() inner loop for a variety of fluids covering
// the main term combinations:
//   Propane : power + Gaussian (many terms, no non-analytic)
//   Water   : power + Gaussian + non-analytic (tests around GenExp)
//   CO2     : power + Gaussian (fewer terms than propane)
//   Methane : pure power terms (l_double path only)
//
// Reports ns/call for GenExp.all() and alphar.all() (full container).
// Use this to measure before/after optimizations.

#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolPropFluid.h"
#include "Helmholtz.h"

#include <chrono>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace CoolProp;
using Clock = std::chrono::high_resolution_clock;

static volatile double sink = 0.0;

static double ns_per_call(Clock::time_point t0, Clock::time_point t1, int N) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / (double)N;
}

static void bench_fluid(const char* name, int N) {
    auto heos = std::make_unique<HelmholtzEOSMixtureBackend>(std::vector<std::string>{name});

    ResidualHelmholtzContainer& alphar = heos->get_components()[0].EOS().alphar;
    ResidualHelmholtzGeneralizedExponential& GenExp = alphar.GenExp;

    const double Tr   = heos->get_components()[0].EOS().reduce.T;
    const double rhor = heos->get_components()[0].EOS().reduce.rhomolar;
    const double tau  = Tr / 300.0;
    const double rho  = rhor * 0.5;

    // Print term counts
    printf("\n── %s ──\n", name);
    printf("  GenExp terms : %zu\n", GenExp.elements.size());
    printf("  flags        : delta_li=%d tau_mi=%d eta1=%d eta2=%d beta1=%d beta2=%d\n",
           (int)GenExp.delta_li_in_u, (int)GenExp.tau_mi_in_u,
           (int)GenExp.eta1_in_u,    (int)GenExp.eta2_in_u,
           (int)GenExp.beta1_in_u,   (int)GenExp.beta2_in_u);

    // ── Warm up ───────────────────────────────────────────────────────────────
    for (int i = 0; i < 500; ++i) {
        double delta = rho / rhor + i * 1e-9;
        HelmholtzDerivatives d;
        GenExp.all(tau, delta, d);
        sink += d.alphar;
    }

    // ── Benchmark GenExp.all() ────────────────────────────────────────────────
    auto t0 = Clock::now();
    for (int i = 0; i < N; ++i) {
        double delta = rho / rhor + i * 1e-9;
        HelmholtzDerivatives d;
        GenExp.all(tau, delta, d);
        sink += d.alphar + d.dalphar_ddelta + d.d2alphar_ddelta2;
    }
    double ns_genexp = ns_per_call(t0, Clock::now(), N);

    // ── Benchmark alphar.all() (full container incl. NonAnalytic etc.) ────────
    // Warm up
    for (int i = 0; i < 500; ++i) {
        double delta = rho / rhor + i * 1e-9;
        HelmholtzDerivatives d = alphar.all(tau, delta, false);
        sink += d.alphar;
    }
    t0 = Clock::now();
    for (int i = 0; i < N; ++i) {
        double delta = rho / rhor + i * 1e-9;
        HelmholtzDerivatives d = alphar.all(tau, delta, false);
        sink += d.alphar + d.dalphar_ddelta + d.d2alphar_ddelta2;
    }
    double ns_alphar = ns_per_call(t0, Clock::now(), N);

    printf("  GenExp.all()  : %7.1f ns/call\n", ns_genexp);
    printf("  alphar.all()  : %7.1f ns/call\n", ns_alphar);
    printf("  NonAnalytic+other overhead: %.1f ns\n", ns_alphar - ns_genexp);
}

int main() {
    const int N = 500'000;
    printf("N = %d per measurement\n", N);

    bench_fluid("Propane",  N);
    bench_fluid("Methane",  N);
    bench_fluid("CO2",      N);
    bench_fluid("Water",    N);

    (void)sink;
    return 0;
}
