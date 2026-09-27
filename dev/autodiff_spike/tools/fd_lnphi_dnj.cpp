// Validate n*dln(phi_i)/dn_j at const T,p against central finite differences in mole numbers.
#include "CoolProp/AbstractState.h"
#include "Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"
#include "Backends/Helmholtz/MixtureDerivatives.h"
#include <cmath>
#include <cstdio>
#include <memory>
using namespace CoolProp;
using MD = MixtureDerivatives;
static std::vector<double> lnphi_at(HelmholtzEOSMixtureBackend& H, const std::vector<double>& n, double T, double p, double rho0) {
    double s = 0; for (double v : n) s += v;
    std::vector<double> x(n.size()); for (std::size_t i = 0; i < n.size(); ++i) x[i] = n[i] / s;
    H.set_mole_fractions(x);
    H.update_TP_guessrho(T, p, rho0);
    std::vector<double> out(n.size());
    for (std::size_t i = 0; i < n.size(); ++i) out[i] = MD::ln_fugacity_coefficient(H, i, XN_INDEPENDENT);
    return out;
}
static void check(const char* fl, std::vector<double> x, double T, double p) {
    std::unique_ptr<AbstractState> AS(AbstractState::factory("HEOS", fl));
    auto& H = *dynamic_cast<HelmholtzEOSMixtureBackend*>(AS.get());
    AS->set_mole_fractions(x);
    AS->update(PT_INPUTS, p, T);
    const double rho0 = AS->rhomolar();
    const std::size_t N = x.size();
    std::printf("%s T=%g p=%g rho=%g phase=%d\n", fl, T, p, rho0, (int)AS->phase());
    std::vector<std::vector<double>> A(N, std::vector<double>(N)), Bd(N, std::vector<double>(N)), C(N, std::vector<double>(N)), F(N, std::vector<double>(N));
    H.update_DmolarT_direct(rho0, T);
    for (std::size_t i = 0; i < N; ++i) {
        double sx = 0;
        for (std::size_t k = 0; k < N; ++k) sx += x[k] * MD::dln_fugacity_coefficient_dxj__constT_p_xi(H, i, k, XN_INDEPENDENT);
        for (std::size_t j = 0; j < N; ++j) {
            A[i][j] = MD::ndln_fugacity_coefficient_dnj__constT_p(H, i, j, XN_INDEPENDENT);
            Bd[i][j] = MD::ndln_fugacity_coefficient_dnj__constT_p(H, i, j, XN_DEPENDENT);
            C[i][j] = MD::dln_fugacity_coefficient_dxj__constT_p_xi(H, i, j, XN_INDEPENDENT) - sx;
        }
    }
    for (std::size_t j = 0; j < N; ++j) {
        const double h = 1e-5;
        std::vector<double> np = x, nm = x; np[j] += h; nm[j] -= h;
        auto lp = lnphi_at(H, np, T, p, rho0), lm = lnphi_at(H, nm, T, p, rho0);
        for (std::size_t i = 0; i < N; ++i) F[i][j] = (lp[i] - lm[i]) / (2 * h);  // n_total = 1 at the base point
    }
    double eA = 0, eB = 0, eC = 0, sc = 0;
    for (std::size_t i = 0; i < N; ++i) for (std::size_t j = 0; j < N; ++j) {
        sc = std::max(sc, std::abs(F[i][j]));
        eA = std::max(eA, std::abs(A[i][j] - F[i][j])); eB = std::max(eB, std::abs(Bd[i][j] - F[i][j])); eC = std::max(eC, std::abs(C[i][j] - F[i][j]));
    }
    std::printf("  max|FD|=%.3g  err: ndlnphi_dnj(XN_INDEP)=%.3g  ndlnphi_dnj(XN_DEP)=%.3g  dlnphi/dx-sum(XN_INDEP)=%.3g\n", sc, eA, eB, eC);
    H.set_mole_fractions(x);
}
int main() {
    check("Methane&Ethane", {0.5, 0.5}, 293.643, 197959);
    check("Methane&Ethane", {0.3, 0.7}, 200, 5e6);
    check("Methane&Ethane&Propane", {0.5, 0.3, 0.2}, 250, 3e6);
    check("Methane&Nitrogen&CarbonDioxide&Ethane&Propane", {0.8, 0.05, 0.05, 0.07, 0.03}, 180, 4e6);
}
