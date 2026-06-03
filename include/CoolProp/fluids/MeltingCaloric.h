#ifndef MELTING_CALORIC_H
#define MELTING_CALORIC_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>
#include "CoolProp/superancillary/superancillary.h"
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

    [[nodiscard]] std::size_t n_samples() const {
        return m_lnp.size();
    }
    [[nodiscard]] double sample_lnp(std::size_t i) const {
        return m_lnp[i];
    }
    [[nodiscard]] double sample_T(std::size_t i) const {
        return m_T[i];
    }
    [[nodiscard]] double sample_rho(std::size_t i) const {
        return m_rho[i];
    }
    [[nodiscard]] double sample_h(std::size_t i) const {
        return m_h[i];
    }
    [[nodiscard]] double sample_s(std::size_t i) const {
        return m_s[i];
    }

    /// sample() probe + per-segment Chebyshev fits vs ln(p). Idempotent; sets built()==true on success.
    void build(HelmholtzEOSMixtureBackend& H);
    [[nodiscard]] bool built() const {
        return m_built;
    }

    [[nodiscard]] double lnp_min() const {
        return m_T_approx ? m_T_approx->xmin() : 0.0;
    }
    [[nodiscard]] double lnp_max() const {
        return m_T_approx ? m_T_approx->xmax() : 0.0;
    }

    [[nodiscard]] double eval_T(double lnp) const {
        return m_T_approx ? m_T_approx->eval(lnp) : std::numeric_limits<double>::quiet_NaN();
    }
    [[nodiscard]] double eval_rho(double lnp) const {
        return m_rho_approx ? m_rho_approx->eval(lnp) : std::numeric_limits<double>::quiet_NaN();
    }
    [[nodiscard]] double eval_h(double lnp) const {
        return m_h_approx ? m_h_approx->eval(lnp) : std::numeric_limits<double>::quiet_NaN();
    }
    [[nodiscard]] double eval_s(double lnp) const {
        return m_s_approx ? m_s_approx->eval(lnp) : std::numeric_limits<double>::quiet_NaN();
    }

    /// Find a (T0, rho0) seed for a target whose caloric values are expressed in
    /// THIS object's build frame (s_cache, h_cache). Returns false if no melting-
    /// line entropy intersection exists. Disambiguates multiple intersections by
    /// closeness in enthalpy.
    bool seed_for_hs(double s_cache, double h_cache, double& T0, double& rho0) const;

    /// The (a1, a2) alpha0 offset pair that was active when build() was called,
    /// or nullopt if build() has not completed successfully.
    [[nodiscard]] std::optional<std::pair<double, double>> stamp() const {
        return m_stamp;
    }

    /// Minimum temperature on the melting curve (may be below the triple-point T
    /// for fluids like water that fold back, e.g. ~251 K for water).
    /// Returns 0.0 if build() has not completed successfully.
    [[nodiscard]] double curve_Tmin() const {
        return m_curve_Tmin;
    }

   protected:
    std::vector<double> m_lnp, m_T, m_rho, m_h, m_s;

    using Approx = CoolProp::superancillary::ChebyshevApproximation1D<Eigen::ArrayXd>;
    bool m_built = false;
    std::optional<Approx> m_T_approx, m_rho_approx, m_h_approx, m_s_approx;
    std::optional<std::pair<double, double>> m_stamp;
    double m_curve_Tmin = 0.0;
};

/// Return a process-global, lazily-built MeltingCaloric for H's (pure) fluid.
/// Returns nullptr if the fluid is not pure or has no melting line. Built once
/// per fluid name; subsequent calls return the cached instance.
std::shared_ptr<MeltingCaloric> get_melting_caloric_cached(HelmholtzEOSMixtureBackend& H);

}  // namespace CoolProp
#endif
