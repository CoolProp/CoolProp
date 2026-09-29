// Verdict stability under ~1e-10 relative input perturbations.  For each state in the input list (exact doubles regenerated from
// bench_ratio's RNG for the N2/C1/C2/nC4/nC5 run), flash NPERT perturbed copies (fresh object each) and count how many are two-phase.
#include "CoolProp/AbstractState.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
using namespace CoolProp;
static bool twophase(double T, double p) {
    std::unique_ptr<AbstractState> AS(AbstractState::factory("GERG2008", "Nitrogen&Methane&Ethane&n-Butane&n-Pentane"));
    AS->set_mole_fractions({0.3797, 0.3225, 0.278, 0.0014, 0.0184});
    try {
        AS->update(PT_INPUTS, p, T);
        return AS->Q() > 0 && AS->Q() < 1;
    } catch (...) {
        return false;
    }
}
int main(int argc, char** argv) {
    const int NPERT = argc > 1 ? std::atoi(argv[1]) : 40;
    const int NS = 10000;
    std::mt19937_64 g(42);
    std::vector<double> Ts(NS), ps(NS);
    for (int i = 0; i < NS; ++i) {
        Ts[i] = 100 + 250 * std::uniform_real_distribution<double>(0, 1)(g);
        ps[i] = std::exp(std::log(1e4) + (std::log(3e7) - std::log(1e4)) * std::uniform_real_distribution<double>(0, 1)(g));
    }
    std::mt19937_64 h(1);
    std::uniform_real_distribution<double> u(-1, 1);
    int idx, nstates = 0, nflaky = 0;
    long ntot = 0, n2 = 0;
    while (std::scanf("%d", &idx) == 1) {
        int k2 = 0;
        for (int k = 0; k < NPERT; ++k)
            k2 += twophase(Ts[idx] * (1 + 3e-10 * u(h)), ps[idx] * (1 + 3e-10 * u(h)));
        ++nstates;
        ntot += NPERT;
        n2 += k2;
        const bool flaky = k2 > 0 && k2 < NPERT;
        nflaky += flaky;
        if (flaky) std::printf("  flaky i=%d T=%.3f p=%.4g: two-phase in %d/%d perturbations\n", idx, Ts[idx], ps[idx], k2, NPERT);
    }
    std::printf("states: %d, flaky (verdict flips under 3e-10 perturbation): %d; overall two-phase fraction %.3f\n", nstates, nflaky,
                double(n2) / ntot);
}
