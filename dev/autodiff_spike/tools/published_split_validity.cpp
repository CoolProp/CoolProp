// Sequence of PT flashes on ONE object (public API only); flag published splits whose phases are not in equilibrium.
#include "CoolProp/AbstractState.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <memory>
using namespace CoolProp;
int main() {
    const std::string f = "Nitrogen&Methane&Ethane&n-Butane&n-Pentane";
    const std::vector<double> z = {0.3797, 0.3225, 0.278, 0.0014, 0.0184};
    std::unique_ptr<AbstractState> AS(AbstractState::factory("GERG2008", f)), L(AbstractState::factory("GERG2008", f)),
      V(AbstractState::factory("GERG2008", f)), Fr(AbstractState::factory("GERG2008", f));
    AS->set_mole_fractions(z);
    double T, p, a, b, c;
    int n = 0, bad = 0, idx = 0;
    double prevT = NAN, prevp = NAN;
    while (std::scanf("%lf %lf %lf %lf %lf", &T, &p, &a, &b, &c) == 5) {
        ++idx;
        if (std::getenv("SET_Z_EACH")) AS->set_mole_fractions(z);
        try {
            AS->update(PT_INPUTS, p, T);
        } catch (...) {
            prevT = T;
            prevp = p;
            continue;
        }
        if (AS->phase() == iphase_twophase) {
            ++n;
            const auto x = AS->mole_fractions_liquid_double(), y = AS->mole_fractions_vapor_double();
            L->set_mole_fractions(x);
            V->set_mole_fractions(y);
            L->specify_phase(iphase_liquid);
            V->specify_phase(iphase_gas);  // imposed phase: direct (rho, T) evaluation, no flash
            try {
                L->update(DmolarT_INPUTS, AS->saturated_liquid_keyed_output(iDmolar), T);
                V->update(DmolarT_INPUTS, AS->saturated_vapor_keyed_output(iDmolar), T);
            } catch (std::exception& e) {
                std::printf("#%d check failed: %s\n", idx, e.what());
                prevT = T;
                prevp = p;
                continue;
            }
            double dl = 0;
            for (std::size_t i = 0; i < z.size(); ++i)
                dl = std::max(dl, std::abs(std::log(V->fugacity(i) / L->fugacity(i))));
            double mb = 0;
            const double Q = AS->Q();
            for (std::size_t i = 0; i < z.size(); ++i)
                mb = std::max(mb, std::abs(z[i] - ((1 - Q) * x[i] + Q * y[i])));
            if (dl > 1e-6 || mb > 1e-8) {
                ++bad;
                Fr->set_mole_fractions(z);
                Fr->update(PT_INPUTS, p, T);
                std::printf(
                  "#%d T=%.3f p=%.6g Q=%.4f max|dlnf|=%.2e mass-bal=%.1e xN2 L/V=%.4f/%.4f | fresh object: Q=%.4f | previous state T=%.3f p=%.6g\n",
                  idx, T, p, Q, dl, mb, x[0], y[0], Fr->Q(), prevT, prevp);
            }
        }
        prevT = T;
        prevp = p;
    }
    std::printf("published two-phase states: %d, non-equilibrium or mass-balance-violating: %d\n", n, bad);
}
