// Microbenchmark: cost of evaluating internal energy u(T, rho) for a
// cubic EOS (SRK / PR), compared against the full Helmholtz reference
// (HEOS).  Motivation: when designing UV tabular inputs for SVDSBTL we
// want to know whether a tabular forward eval is actually faster than
// just calling a cubic EOS directly.
//
// Each timed loop runs `update(DmolarT_INPUTS, rho, T)` followed by
// `umolar()`.  update() clears the property cache every call, so
// umolar() genuinely recomputes each iteration -- this is the true
// forward-eval cost, not a cache hit.
//
// Build:
//   cmake -B build -DCOOLPROP_BUILD_CUBIC_U_BENCH=ON -DCMAKE_BUILD_TYPE=Release
//   cmake --build build --target bench_cubic_u -j
//
// Run:
//   ./build/bench_cubic_u                 # Methane at (300 K, 5000 mol/m^3)
//   ./build/bench_cubic_u Water 350 20000 # custom fluid, T (K), rho (mol/m^3)
//
// Env:
//   BENCH_REPEATS (default 2,000,000)
//
// Dev-only microbenchmark: the lints below (varargs printf, argv/atof handling, the
// deliberate bit-exact memcmp on doubles, the forwarding-ref timing helper, main-may-throw)
// are all appropriate for throwaway benchmark tooling and match the sibling dev benchmarks.
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,cppcoreguidelines-missing-std-forward,cert-exp42-c,bugprone-suspicious-memory-comparison,cert-flp37-c,bugprone-exception-escape,cert-err34-c,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays)

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "CoolProp/AbstractState.h"
#include "CoolProp/DataStructures.h"
#include "Backends/Cubics/CubicBackend.h"
#include "Backends/Cubics/GeneralizedCubic.h"

namespace {

std::size_t env_size_t(const char* name, std::size_t fallback) {
    const char* v = std::getenv(name);
    if (v == nullptr || v[0] == 0) {
        return fallback;
    }
    try {
        return static_cast<std::size_t>(std::stoul(v));
    } catch (...) {
        return fallback;
    }
}

// Time a callable that returns a double; returns mean ns/call.
template <typename F>
double timeit_double(F&& f, std::size_t repeats) {
    volatile double sink = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < repeats; ++i) {
        sink += f();
    }
    const auto t1 = std::chrono::steady_clock::now();
    (void)sink;
    return std::chrono::duration<double, std::nano>(t1 - t0).count() / static_cast<double>(repeats);
}

// Time one backend's u(T, rho) forward eval; prints ns/call and the
// value it computed (so we can confirm the backends agree on the state).
// When impose_phase is true the single-phase index is imposed before the
// timed loop, which lets the cubic backend skip the superancillary
// phase-determination work in update_DmolarT().
void bench_backend(const std::string& backend, const std::string& fluid, double T, double rho_molar, std::size_t repeats, bool impose_phase) {
    const std::string label = backend + (impose_phase ? " (phase)" : "");
    std::shared_ptr<CoolProp::AbstractState> AS;
    try {
        AS = std::shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory(backend, fluid));
    } catch (const std::exception& e) {
        std::printf("  %-18s  FAILED to construct: %s\n", label.c_str(), e.what());
        return;
    }

    if (impose_phase) {
        // Any imposed index sends update_DmolarT() down the direct-pressure
        // branch; the value only sets the _phase label, not the u result.
        AS->specify_phase(CoolProp::iphase_gas);
    }

    // Warm-up + correctness probe (also pays one-time lazy init costs).
    double u_probe = 0.0;
    try {
        AS->update(CoolProp::DmolarT_INPUTS, rho_molar, T);
        u_probe = AS->umolar();
    } catch (const std::exception& e) {
        std::printf("  %-18s  FAILED to evaluate: %s\n", label.c_str(), e.what());
        return;
    }

    const double ns = timeit_double(
      [&]() {
          AS->update(CoolProp::DmolarT_INPUTS, rho_molar, T);
          return AS->umolar();
      },
      repeats);

    std::printf("  %-18s  %8.1f ns/call   u = %.6g J/mol\n", label.c_str(), ns, u_probe);
}

// Time just the single residual derivative umolar actually needs:
// alphar(tau, delta, z, itau=1, idelta=0), called directly on the cubic
// object.  This is the floor cost if the backend computed only the one
// derivative instead of all 15 via CubicResidualHelmholtz::all().
void bench_raw_dalphar_dtau(const std::string& backend, const std::string& fluid, double T, double rho_molar, std::size_t repeats) {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, fluid));
    auto* cb = dynamic_cast<CoolProp::AbstractCubicBackend*>(AS.get());
    if (cb == nullptr) {
        std::printf("  %-18s  not a cubic backend\n", backend.c_str());
        return;
    }
    std::shared_ptr<AbstractCubic>& cubic = cb->get_cubic();
    const std::vector<double> z = {1.0};
    const double tau = cubic->get_Tr() / T;
    const double delta = rho_molar / cubic->get_rhor();

    const double ns = timeit_double([&]() { return cubic->alphar(tau, delta, z, 1, 0); }, repeats);
    std::printf("  %-18s  %8.1f ns/call   (raw alphar itau=1 only)\n", (backend + " d/dtau").c_str(), ns);
}

// Verify the shared-intermediate assembly in CubicResidualHelmholtz::all() is bit-for-bit
// identical to 15 independent alphar() calls (alphar() itself is unchanged).  Returns the
// number of bitwise mismatches found across a small (tau, delta) grid.
int verify_bitexact(const std::string& backend, const std::string& fluid) {
    std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory(backend, fluid));
    auto* cb = dynamic_cast<CoolProp::AbstractCubicBackend*>(AS.get());
    if (cb == nullptr) return 0;
    std::shared_ptr<AbstractCubic>& cubic = cb->get_cubic();
    const std::vector<double> z = {1.0};
    const std::size_t combos[15][2] = {{0, 0}, {1, 0}, {0, 1}, {2, 0}, {1, 1}, {0, 2}, {3, 0}, {2, 1},
                                       {1, 2}, {0, 3}, {4, 0}, {3, 1}, {2, 2}, {1, 3}, {0, 4}};
    int mismatches = 0;
    for (int it_tau = 0; it_tau <= 10; ++it_tau) {
        const double tau = 0.7 + 0.13 * it_tau;
        for (int it_del = 0; it_del <= 10; ++it_del) {
            const double delta = 0.2 + 0.27 * it_del;
            // shared-intermediate assembly, mirroring CubicResidualHelmholtz::all()
            std::array<double, 5> pm{}, pp{}, ta{};
            for (std::size_t n = 0; n <= 4; ++n) {
                pm[n] = cubic->psi_minus(delta, z, 0, n);
                pp[n] = cubic->psi_plus(delta, z, n);
                ta[n] = cubic->tau_times_a(tau, z, n);
            }
            const double den = cubic->get_R_u() * cubic->get_Tr();
            for (auto& c : combos) {
                const std::size_t it = c[0], id = c[1];
                const double assembled = (it == 0 ? pm[id] : 0.0) - ta[it] / den * pp[id];
                const double direct = cubic->alphar(tau, delta, z, it, id);
                if (std::memcmp(&assembled, &direct, sizeof(double)) != 0) ++mismatches;
            }
        }
    }
    std::printf("  %-10s  bit-exact check: %d mismatch(es)\n", backend.c_str(), mismatches);
    return mismatches;
}

}  // namespace

int main(int argc, char** argv) {
    std::string fluid = (argc > 1) ? argv[1] : "Methane";
    const double T = (argc > 2) ? std::atof(argv[2]) : 300.0;           // K
    const double rho_molar = (argc > 3) ? std::atof(argv[3]) : 5000.0;  // mol/m^3
    const std::size_t repeats = env_size_t("BENCH_REPEATS", 2000000);

    std::printf("u(T, rho) forward-eval benchmark\n");
    std::printf("  fluid   = %s\n", fluid.c_str());
    std::printf("  T       = %g K\n", T);
    std::printf("  rho     = %g mol/m^3\n", rho_molar);
    std::printf("  repeats = %zu\n\n", repeats);

    // Cubic EOS first (the thing of interest), HEOS as a reference.
    // Each cubic is timed twice: default (lets update_DmolarT do
    // superancillary phase determination) and with the phase imposed.
    bench_backend("SRK", fluid, T, rho_molar, repeats, false);
    bench_backend("SRK", fluid, T, rho_molar, repeats, true);
    bench_backend("PR", fluid, T, rho_molar, repeats, false);
    bench_backend("PR", fluid, T, rho_molar, repeats, true);
    bench_backend("HEOS", fluid, T, rho_molar, repeats, false);

    std::printf("\n  -- floor: single needed residual derivative --\n");
    bench_raw_dalphar_dtau("SRK", fluid, T, rho_molar, repeats);
    bench_raw_dalphar_dtau("PR", fluid, T, rho_molar, repeats);

    std::printf("\n  -- bit-exactness of shared-intermediate all() vs 15 alphar() calls --\n");
    verify_bitexact("SRK", fluid);
    verify_bitexact("PR", fluid);

    return 0;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg,cppcoreguidelines-missing-std-forward,cert-exp42-c,bugprone-suspicious-memory-comparison,cert-flp37-c,bugprone-exception-escape,cert-err34-c,cppcoreguidelines-pro-bounds-pointer-arithmetic,modernize-avoid-c-arrays)
