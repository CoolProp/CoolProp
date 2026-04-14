// Compare build() vs build_isochoric() phase envelope tracers.
// Prints: method, mixture, N_pts, time_ms, built(0/1), T_max_K, p_max_bar
//
// Build:
//   cmake .. -DCOOLPROP_MY_MAIN=../dev/bench_phase_envelope.cpp -DCOOLPROP_STATIC_LIBRARY=ON
//   ninja Main
// Run:
//   ./Main

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>
#include "AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/Helmholtz/PhaseEnvelopeRoutines.h"
#include "CoolProp.h"

using namespace CoolProp;
using clk = std::chrono::steady_clock;

struct MixtureSpec {
    const char* label;
    const char* fluids;
    std::vector<double> z;
};

// Type-I and Type-II binary and ternary mixtures used in natural-gas applications.
// build() is known to succeed for all of these.
static const MixtureSpec MIXTURES[] = {
    {"CH4/C2H6 50/50",   "Methane&Ethane",             {0.5, 0.5}},
    {"CH4/nC4 80/20",    "Methane&n-Butane",            {0.8, 0.2}},
    {"CH4/C3H8/N2",      "Methane&Propane&Nitrogen",    {0.70, 0.20, 0.10}},
    {"CH4/C2H6/C3H8",    "Methane&Ethane&Propane",      {0.70, 0.20, 0.10}},
    {"CH4/CO2/N2",       "Methane&CarbonDioxide&Nitrogen", {0.80, 0.10, 0.10}},
};

struct EnvResult {
    const char* method;
    const char* label;
    std::size_t n_pts;
    double time_ms;
    bool built;
    double T_max;
    double p_max;
};

static EnvResult run(const char* method, const MixtureSpec& mix) {
    auto AS = std::shared_ptr<AbstractState>(AbstractState::factory("HEOS", mix.fluids));
    auto* heos = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    heos->set_mole_fractions(mix.z);

    auto t0 = clk::now();
    if (std::string(method) == "build") {
        PhaseEnvelopeRoutines::build(*heos);
    } else {
        PhaseEnvelopeRoutines::build_isochoric(*heos);
    }
    double ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();

    const auto& env = heos->PhaseEnvelope;
    double T_max = env.T.empty() ? 0.0 : *std::max_element(env.T.begin(), env.T.end());
    double p_max = env.p.empty() ? 0.0 : *std::max_element(env.p.begin(), env.p.end());

    return {method, mix.label, env.T.size(), ms, env.built, T_max, p_max / 1e5};
}

int main() {
    static const char* methods[] = {"build", "build_isochoric"};

    printf("%-22s  %-16s  %6s  %10s  %5s  %8s  %9s\n",
           "mixture", "method", "N_pts", "time_ms", "built", "T_max_K", "p_max_bar");
    printf("%-22s  %-16s  %6s  %10s  %5s  %8s  %9s\n",
           "----------------------", "----------------", "------", "----------",
           "-----", "--------", "---------");
    fflush(stdout);

    for (const auto& mix : MIXTURES) {
        for (int rep = 0; rep < 3; ++rep) {
            for (const char* method : methods) {
                try {
                    auto r = run(method, mix);
                    printf("%-22s  %-16s  %6zu  %10.1f  %5d  %8.2f  %9.2f\n",
                           r.label, r.method, r.n_pts, r.time_ms,
                           (int)r.built, r.T_max, r.p_max);
                } catch (const std::exception& e) {
                    printf("%-22s  %-16s  ERROR: %s\n", mix.label, method, e.what());
                }
                fflush(stdout);
            }
        }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
