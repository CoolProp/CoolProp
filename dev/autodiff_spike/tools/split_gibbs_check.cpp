// In-model check: is the published split's Gibbs energy below the best single-phase root at the feed?
#include "CoolProp/AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/Helmholtz/MixtureDerivatives.h"
#include "gerg_cheb_solver.hpp"
#include <cmath>
#include <cstdio>
#include <memory>
using namespace CoolProp;
// molar Gibbs / RT with the mixture's ideal-gas part up to terms common to all phases at the same T,p:
// g/RT = sum_i x_i (ln x_i + ln phi_i)  (+ ln p + mu0_i(T)/RT terms that cancel between split and single phase)
static double gRT(HelmholtzEOSMixtureBackend& H) {
    double g = 0; const auto x = H.get_mole_fractions();
    for (std::size_t i = 0; i < x.size(); ++i) if (x[i] > 0) g += x[i] * (std::log(x[i]) + MixtureDerivatives::ln_fugacity_coefficient(H, i, XN_INDEPENDENT));
    return g;
}
int main(int argc, char** argv) {
    const bool air = argc > 1;
    const std::string FL = air ? "Nitrogen&Oxygen&Argon&CarbonDioxide&Water" : "Methane&HydrogenSulfide";
    const std::vector<double> Z = air ? std::vector<double>{0.7654, 0.2053, 0.0090, 0.0003, 0.02} : std::vector<double>{0.5, 0.5};
    double T, p; int idx;
    std::unique_ptr<AbstractState> AS(AbstractState::factory("GERG2008", FL));
    std::unique_ptr<AbstractState> F(AbstractState::factory("GERG2008", FL));
    auto* H = dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    auto* HF = dynamic_cast<HelmholtzEOSMixtureBackend*>(F.get());
    std::vector<double> z = Z;
    gergcheb::Solver sv; sv.build(HF, 4.0, 1.3 * HF->Reducing->Tr(z) / (air ? 250.0 : 180.0), 1e-6);
    while (std::scanf("%d %lf %lf", &idx, &T, &p) == 3) {
        AS->set_mole_fractions(z);
        AS->update(PT_INPUTS, p, T);
        const double Q = AS->Q();
        // best single phase: min g over mechanically stable roots at the feed
        F->set_mole_fractions(z);
        gergcheb::Solver::State S; sv.assemble(T, z, S); gergcheb::Solver::Root r[gergcheb::MAXROOTS];
        const int n = sv.roots(S, p, r, gergcheb::Solver::Polish::All);
        // reference single phase: the kernel's spinodal-branch selection (NOT min-g over all roots: the
        // alpha^r-well roots near delta ~ 1 have spuriously low g)
        const int ks = sv.select(S, p, r, n);
        double rs = ks >= 0 ? r[ks].rho : NAN, gs = 1e300;
        if (ks >= 0) { HF->update_DmolarT_direct(rs, T); gs = gRT(*HF); }
        if (Q > 0 && Q < 1) {
            auto& L = *H->SatL; auto& V = *H->SatV;
            const double g2 = (1 - Q) * gRT(L) + Q * gRT(V);
            double spread = 0; for (std::size_t i = 0; i < z.size(); ++i) spread = std::max(spread, std::abs(L.get_mole_fractions()[i] - V.get_mole_fractions()[i]));
            std::printf("%5d T=%7.2f p=%6.2f MPa  split Q=%.3f rhoL=%8.1f rhoV=%8.1f x0 L/V=%.4f/%.4f xLast L/V=%.4g/%.4g | g_split-g_single = %+.3e  (single rho=%.1f)\n", idx, T, p / 1e6, Q,
                        L.rhomolar(), V.rhomolar(), L.get_mole_fractions()[0], V.get_mole_fractions()[0], L.get_mole_fractions().back(), V.get_mole_fractions().back(), g2 - gs, rs);
        } else {
            std::printf("%5d T=%7.2f p=%6.2f MPa  single rho=%.1f (best single-phase root %.1f)\n", idx, T, p / 1e6, AS->rhomolar(), rs);
        }
    }
}
