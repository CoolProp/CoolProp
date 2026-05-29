#ifndef MELTING_CALORIC_H
#define MELTING_CALORIC_H

#include <vector>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include "superancillary/superancillary.h"
#include <Eigen/Dense>

namespace CoolProp {

class HelmholtzEOSMixtureBackend;  // forward decl

/// Caloric properties along the solid-liquid melting curve, parametrized by
/// ln(p) (the monotone curve variable; T is double-valued on water's curve).
class MeltingCaloric
{
   public:
    MeltingCaloric() = default;

    /// Walk the melting curve over [pmin, pmax] into raw (lnp, T, rho, h, s)
    /// samples. n>=8. EOS density/caloric come from update(PT_INPUTS,...).
    void sample(HelmholtzEOSMixtureBackend& H, std::size_t n);

    std::size_t n_samples() const { return m_lnp.size(); }
    double sample_lnp(std::size_t i) const { return m_lnp[i]; }
    double sample_T(std::size_t i) const { return m_T[i]; }
    double sample_rho(std::size_t i) const { return m_rho[i]; }
    double sample_h(std::size_t i) const { return m_h[i]; }
    double sample_s(std::size_t i) const { return m_s[i]; }

    /// sample() probe + per-segment Chebyshev fits vs ln(p). Idempotent; sets built()==true on success.
    void build(HelmholtzEOSMixtureBackend& H);
    bool built() const { return m_built; }

    double lnp_min() const { return m_T_approx ? m_T_approx->xmin() : 0.0; }
    double lnp_max() const { return m_T_approx ? m_T_approx->xmax() : 0.0; }

    double eval_T(double lnp) const { return m_T_approx->eval(lnp); }
    double eval_rho(double lnp) const { return m_rho_approx->eval(lnp); }
    double eval_h(double lnp) const { return m_h_approx->eval(lnp); }
    double eval_s(double lnp) const { return m_s_approx->eval(lnp); }

   protected:
    std::vector<double> m_lnp, m_T, m_rho, m_h, m_s;

    using Approx = CoolProp::superancillary::ChebyshevApproximation1D<Eigen::ArrayXd>;
    bool m_built = false;
    std::optional<Approx> m_T_approx, m_rho_approx, m_h_approx, m_s_approx;
};

}  // namespace CoolProp
#endif
