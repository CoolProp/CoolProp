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
// NOLINTBEGIN(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg)

#include <chrono>
#include <cstdio>
#include <cstdlib>
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
void bench_backend(const std::string& backend, const std::string& fluid, double T, double rho_molar, std::size_t repeats,
                   bool impose_phase) {
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
void bench_raw_dalphar_dtau(const std::string& backend, const std::string& fluid, double T, double rho_molar,
                            std::size_t repeats) {
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

}  // namespace

int main(int argc, char** argv) {
    std::string fluid = (argc > 1) ? argv[1] : "Methane";
    const double T = (argc > 2) ? std::atof(argv[2]) : 300.0;          // K
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

    return 0;
}

// NOLINTEND(cppcoreguidelines-pro-type-vararg,cert-err33-c,hicpp-vararg)
