#ifndef TRANSPORTROUTINES_H
#define TRANSPORTROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

class TransportRoutines
{
public:
    /**
    \brief The general dilute gas viscosity from used for ECS

    \f[
    \eta^0 = \displaystyle\frac{26.692\times 10^{-9}\sqrt{MT}}{\sigma^2\Omega^{(2,2)}(T^*)}
    \f]
    \f[
    \Omega^{(2,2)}(T^*)=1.16145(T^*)^{-0.14874}+0.52487\exp(-0.77320T^*)+2.16178\exp(-2.43787T^*)
    \f]
    with \f$T^* = \frac{T}{\varepsilon/k}\f$ and \f$\sigma\f$ in nm, M is in kg/kmol. Yields viscosity in Pa-s.
    */
    static long double general_dilute_gas_viscosity(HelmholtzEOSMixtureBackend &HEOS);
    
    /**
    \brief A dilute gas term that has a form like

    \f[
    \eta^0 = \displaystyle\frac{A\sqrt{MT}}{\sigma^2\mathfrak{S}(T^*)}
    \f]
    \f[
    \mathfrak{S}(T^*)=\exp\left(\sum_ia_i[\ln T^*]^i\right)
    \f]
    with \f$T^* = \frac{T}{\varepsilon/k}\f$ and \f$\sigma\f$ in nm, M is in kg/kmol. Yields viscosity in Pa-s.
    */
    static long double dilute_gas_viscosity_collision_integral(HelmholtzEOSMixtureBackend &HEOS);
};

}; /* namespace CoolProp */
#endif