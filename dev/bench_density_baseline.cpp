// Baseline characterization of HelmholtzEOSMixtureBackend::solver_rho_Tp
// across four mixture systems. For each (T,p) cell we record:
//   - full PT_flash time (ms, best-of-N)        : user-visible cost
//   - phase classification (CoolProp phase index)
//   - post-flash density (mol/m^3)
//   - isolated solver_rho_Tp time (μs, best-of-M): marginal density-solve cost
//   - solver_rho_Tp success flag
//
// NOTE: tier capture (which fallback resolved the flash) is not available — the
// `get_last_pt_flash_tier` helper referenced in dev/bench_ch4_h2s_grid.cpp was
// never committed. Adding it is a separate follow-up issue.
//
// Issue: CoolProp-aor (baseline characterization)
//
// Build (from build/):
//   cmake .. -DCOOLPROP_MY_MAIN=../dev/bench_density_baseline.cpp \
//            -DCOOLPROP_STATIC_LIBRARY=ON
//   ninja Main
//
// Run:
//   ./Main amarillo > dev/bench_density_amarillo.csv
//   ./Main ch4h2s   > dev/bench_density_ch4h2s.csv
//   ./Main co2n2    > dev/bench_density_co2n2.csv
//   ./Main c10c1    > dev/bench_density_c10c1.csv

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "AbstractState.h"
#include "Backends/Helmholtz/FlashRoutines.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "CoolProp.h"

using namespace CoolProp;
using clk = std::chrono::steady_clock;

struct System
{
    const char* key;
    const char* fluids;
    std::vector<double> z;
    int NT, NP;
    double T_min, T_max;  // K
    double p_min, p_max;  // Pa
    int flash_repeat;     // outer best-of for full flash
    int solver_inner;     // inner loop count for isolated solver_rho_Tp
};

static const std::vector<System> SYSTEMS = {
  // Amarillo AGA-8 reference natural gas
  {"amarillo",
   "Methane&Nitrogen&CarbonDioxide&Ethane&Propane&IsoButane&n-Butane&Isopentane&n-Pentane&n-Hexane",
   {0.906724, 0.031284, 0.004676, 0.045279, 0.008280, 0.001037, 0.001563, 0.000321, 0.000443, 0.000393},
   60,
   60,
   140.0,
   260.0,
   0.5e5,
   80e5,
   5,
   200},

  // CH4/H2S — sour gas, H2S critical ~373 K
  {"ch4h2s", "Methane&HydrogenSulfide", {0.70, 0.30}, 60, 60, 200.0, 420.0, 1e5, 150e5, 5, 200},

  // CO2/N2 — cuts CO2 critical (304 K, 73.8 bar)
  {"co2n2", "CarbonDioxide&Nitrogen", {0.80, 0.20}, 60, 60, 250.0, 400.0, 10e5, 250e5, 5, 200},

  // n-Decane/Methane — heavy + light, retrograde region
  {"c10c1", "n-Decane&Methane", {0.30, 0.70}, 60, 60, 250.0, 500.0, 1e5, 500e5, 5, 200},
};

static const System* find_system(const std::string& key) {
    for (const auto& s : SYSTEMS) {
        if (key == s.key) return &s;
    }
    return nullptr;
}

int main(int argc, char** argv) {
    std::string key = (argc > 1) ? argv[1] : "amarillo";
    const System* sys = find_system(key);
    if (!sys) {
        std::fprintf(stderr, "Unknown system: %s\nAvailable: ", key.c_str());
        for (const auto& s : SYSTEMS)
            std::fprintf(stderr, "%s ", s.key);
        std::fprintf(stderr, "\n");
        return 1;
    }
    std::fprintf(stderr, "System: %s  fluids=%s  grid=%dx%d  T=[%.0f,%.0f]K  p=[%.2f,%.1f]bar\n", sys->key, sys->fluids, sys->NT, sys->NP, sys->T_min,
                 sys->T_max, sys->p_min / 1e5, sys->p_max / 1e5);

    auto AS = shared_ptr<AbstractState>(AbstractState::factory("HEOS", sys->fluids));
    AS->set_mole_fractions(sys->z);
    auto* HEOS = static_cast<HelmholtzEOSMixtureBackend*>(AS.get());

    // Warm-up
    try {
        AS->update(PT_INPUTS, 1e5, 0.5 * (sys->T_min + sys->T_max));
    } catch (...) {
    }

    std::printf("T_K,p_Pa,phase,flash_ms,rho_molm3,solver_us,solver_ok\n");

    for (int j = 0; j < sys->NP; ++j) {
        double log_p = std::log(sys->p_min) + (std::log(sys->p_max) - std::log(sys->p_min)) * j / (sys->NP - 1);
        double p = std::exp(log_p);

        for (int i = 0; i < sys->NT; ++i) {
            double T = sys->T_min + (sys->T_max - sys->T_min) * i / (sys->NT - 1);

            // ---- Full PT flash: best-of-flash_repeat ----
            double best_ms = 1e9;
            int phase_out = -1;
            double rho_out = std::nan("");
            bool flash_ok = false;

            for (int r = 0; r < sys->flash_repeat; ++r) {
                auto t0 = clk::now();
                int phase = -1;
                double rho = std::nan("");
                bool ok = false;
                try {
                    AS->update(PT_INPUTS, p, T);
                    phase = static_cast<int>(AS->phase());
                    rho = AS->rhomolar();
                    ok = true;
                } catch (...) {
                    phase = -1;
                }
                double ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
                if (ms < best_ms) {
                    best_ms = ms;
                    phase_out = phase;
                    rho_out = rho;
                    flash_ok = ok;
                }
            }

            // ---- Isolated solver_rho_Tp: inner-loop best-of ----
            // Only meaningful if the flash succeeded and gave us a known density.
            // We pass rhomolar_guess=-1 to exercise the SRK-guess path (cold start
            // for the rootfinder, warm for the backend).
            double best_us = std::nan("");
            int solver_ok = 0;
            if (flash_ok && std::isfinite(rho_out)) {
                // First, verify the solver converges at all on this point.
                try {
                    double rho_check = HEOS->solver_rho_Tp(T, p, -1.0);
                    if (std::isfinite(rho_check)) solver_ok = 1;
                } catch (...) {
                    solver_ok = 0;
                }

                if (solver_ok) {
                    // Repeat M times, total time / M = per-call.
                    auto t0 = clk::now();
                    for (int k = 0; k < sys->solver_inner; ++k) {
                        try {
                            (void)HEOS->solver_rho_Tp(T, p, -1.0);
                        } catch (...) {
                        }
                    }
                    double total_us = std::chrono::duration<double, std::micro>(clk::now() - t0).count();
                    best_us = total_us / sys->solver_inner;
                }
            }

            std::printf("%.4f,%.2f,%d,%.4f,%.6g,%.4f,%d\n", T, p, phase_out, best_ms, rho_out, best_us, solver_ok);
        }
    }

    std::fprintf(stderr, "Done: %s\n", sys->key);
    return 0;
}
