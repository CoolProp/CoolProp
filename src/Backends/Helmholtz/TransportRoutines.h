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
    \brief The dilute gas viscosity term that is based on collision integral or effective cross section

    \f[
    \eta^0 = \displaystyle\frac{A\sqrt{MT}}{\sigma^2\mathfrak{S}(T^*)}
    \f]
    \f[
    \mathfrak{S}(T^*)=\exp\left(\sum_ia_i[\ln T^*]^{t_i}\right)
    \f]
    with \f$T^* = \frac{T}{\varepsilon/k}\f$ and \f$\sigma\f$ in nm, M is in kg/kmol. Yields viscosity in Pa-s.
    */
    static long double dilute_gas_viscosity(HelmholtzEOSMixtureBackend &HEOS);

    /** 
    \brief The modified Batschinski-Hildebrand contribution to the viscosity

    \f[
    \Delta\eta = \displaystyle\sum_{i}a_{i}\frac{\delta^{d_i}}{\tau_{t_j}}+\left(\displaystyle\sum_{k}\displaystyle\frac{f_k}{\tau^{k-1}}\right)\delta^c\left(\frac{1}{\delta_0(\tau)-\delta}-\frac{1}{\delta_0(\tau)}\right)
    \f]
    where \f$\tau = T_c/T\f$ and \f$\delta = \rho/\rho_c\f$ 
    \f[
    \delta_0(\tau) = \displaystyle\frac{\displaystyle\sum_{i}g_i\tau^{h_i}}{\displaystyle\sum_{i}p_i\tau^{q_i}}
    \f]
    The more general form of \f$\delta_0(\tau)\f$ is selected in order to be able to handle all the forms in the literature 
    */
    static long double modified_Batschinski_Hildebrand_viscosity_term(HelmholtzEOSMixtureBackend &HEOS);

    /** 
    \brief The initial density dependence term \f$B_{\eta}\f$ from Rainwater-Friend theory

    The total contribution from this term is given by
    \f[
    \eta_{RF} = \eta_0B_{\eta}\rho
    \f]
    where \f$\eta_0\f$ is the dilute gas viscosity in Pa-s and \f$\rho\f$ is the molar density in mol/m\f$^3\f$ and \f$B_{\eta}\f$ is in m^3/mol.

    \f[
    B_{\eta}(T) = B_{\eta}^*(T^*)N_A\sigma_{\eta}^3
    \f]
    where \f$N_A\f$ is Avogadros number \f$6.022\times 10^{23}\f$ mol\f$^{-1}\f$ and \f$\sigma_{\eta}\f$ is in m.

    \f[
    B_{\eta}^*(T^*) = \displaystyle\sum_ib_i(T^*)^{c_i}
    \f]

    IMPORTANT: This function returns \f$B_{\eta}\f$, not \f$\eta_{RF}\f$
    */
    static long double initial_density_dependence_viscosity_term(HelmholtzEOSMixtureBackend &HEOS);
};

}; /* namespace CoolProp */
#endif