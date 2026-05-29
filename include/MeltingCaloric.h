#ifndef MELTING_CALORIC_H
#define MELTING_CALORIC_H

#include <vector>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>

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

   protected:
    std::vector<double> m_lnp, m_T, m_rho, m_h, m_s;
};

}  // namespace CoolProp
#endif
