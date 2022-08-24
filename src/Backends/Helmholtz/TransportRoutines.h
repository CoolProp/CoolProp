#ifndef TRANSPORTROUTINES_H
#define TRANSPORTROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

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
    static CoolPropDbl viscosity_dilute_kinetic_theory(HelmholtzEOSMixtureBackend& HEOS);

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
    static CoolPropDbl viscosity_dilute_collision_integral(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief A dilute gas viscosity term formed of summation of power terms

    \f[
    \eta^0 = \displaystyle\sum_ia_iT^{t_i}
    \f]
    with T in K, \f$\eta^0\f$ in Pa-s
    */
    static CoolPropDbl viscosity_dilute_powers_of_T(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief A dilute gas viscosity term formed of summation of power terms of the reduced temperature

    \f[
    \eta^0 = \displaystyle\sum_ia_i(T/T_c)^{t_i}
    \f]
    with T in K, \f$\eta^0\f$ in Pa-s
    */
    static CoolPropDbl viscosity_dilute_powers_of_Tr(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl viscosity_dilute_collision_integral_powers_of_T(HelmholtzEOSMixtureBackend& HEOS);

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
    B_{\eta}^*(T^*) = \displaystyle\sum_ib_i(T^*)^{t_i}
    \f]

    IMPORTANT: This function returns \f$B_{\eta}\f$, not \f$\eta_{RF}\f$
    */
    static CoolPropDbl viscosity_initial_density_dependence_Rainwater_Friend(HelmholtzEOSMixtureBackend& HEOS);

    /**
     * \brief An empirical form for the initial density dependence
     *
     * Given by the polynomial-like form
     * \f[
     *  \eta^1 = \sum_i n_i\delta^{d_i}\tau^{t_i}
     * \f]
     * where the output is in Pa-s
     */
    static CoolPropDbl viscosity_initial_density_dependence_empirical(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief The modified Batschinski-Hildebrand contribution to the viscosity

    \f[
    \Delta\eta = \displaystyle\sum_{i}a_{i}\delta^{d1_i}\tau^{t1_j}\exp(\gamma_i\delta^{l_i})+\left(\displaystyle\sum_{i}f_i\delta^{d2_i}\tau^{t2_i}\right)\left(\frac{1}{\delta_0(\tau)-\delta}-\frac{1}{\delta_0(\tau)}\right)
    \f]
    where \f$\tau = T_c/T\f$ and \f$\delta = \rho/\rho_c\f$
    \f[
    \delta_0(\tau) = \displaystyle\frac{\displaystyle\sum_{i}g_i\tau^{h_i}}{\displaystyle\sum_{i}p_i\tau^{q_i}}
    \f]
    The more general form of \f$\delta_0(\tau)\f$ is selected in order to be able to handle all the forms in the literature
    */
    static CoolPropDbl viscosity_higher_order_modified_Batschinski_Hildebrand(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl viscosity_dilute_ethane(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_dilute_cyclohexane(HelmholtzEOSMixtureBackend& HEOS);

    /** \brief Viscosity hardcoded for Methanol
     *
     * From Xiang et al., A New Reference Correlation for the Viscosity of Methanol, J. Phys. Chem. Ref. Data, Vol. 35, No. 4, 2006
     */
    static CoolPropDbl viscosity_methanol_hardcoded(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl viscosity_heavywater_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_water_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_helium_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_R23_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_m_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_o_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_p_xylene_hardcoded(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl viscosity_toluene_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_ethane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_hydrogen_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_benzene_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_hexane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl viscosity_heptane_higher_order_hardcoded(HelmholtzEOSMixtureBackend& HEOS);

    /**
     * @brief Higher-order viscosity term from friction theory of Sergio Quinones-Cisneros
     *
     * Several functional forms have been proposed and this function attempts to handle all of them
     * \f$ \eta_{HO} = \kappa_ap_a + \kappa_r\Delta p_r + \kappa_i p_{id}+\kappa_{aa}p_a^2 + \kappa_{drdr}\Delta p_r^2 + \kappa_{rr}p_{r}^2 + \kappa_{ii}p_{id}^2 +\kappa_{rrr}p_r^3 + \kappa_{aaa}p_a^3
     *
     * Watch out that sometimes it is \f$\Delta p_r\f$ and other times it is \f$p_r\f$!
     *
	 * 1e5 for conversion from Pa -> bar
	 *
     * \f[ p_r = T \frac{\partial p}{\partial T}\right|_{\rho}/1e5 \f]
     * \f[ p_a = p - p_r \f]
     * \f[ p_{id} = \rho R T \f] / 1e5 \f]
     * \f[ \Delta p_r = p_r - p_{id} \f]
     * \f[ \psi_1 = \exp(\tau)-c_1 \f]
     * \f[ \psi_2 = \exp(\tau^2)-c_2 \f]
     * \f[ \kappa_i = (A_{i,0} + A_{i,1}\psi_1 + A_{i,2}\psi_2)\tau \f]
     * \f[ \kappa_a = (A_{a,0} + A_{a,1}\psi_1 + A_{a,2}\psi_2)\tau^{N_a} \f]
     * \f[ \kappa_{aa} = (A_{aa,0} + A_{aa,1}\psi_1 + A_{aa,2}\psi_2)\tau^{N_{aa}} \f]
     * \f[ \kappa_r = (A_{r,0} + A_{r,1}\psi_1 + A_{r,2}\psi_2)\tau^{N_r} \f]
     * \f[ \kappa_{rr} = (A_{rr,0} + A_{rr,1}\psi_1 + A_{rr,2}\psi_2)\tau^{N_{rr}} \f]
     * \f[ \kappa_{drdr} = (A_{drdr,0} + A_{drdr,1}\psi_1 + A_{drdr,2}\psi_2)\tau^{N_{drdr}} \f]
     * \f[ \kappa_{aa} = (A_{aa,0} + F_{Aaa,1}\psi_1 + F.Aaa[2]\psi_2)\tau^{N_{aa}} \f]
     * \f[ \kappa_{rrr} = (A_{rrr,0} + A_{rrr,1}\psi_1 + A_{rrr,2}\psi_2)\tau^{N_{rrr}} \f]
     * \f[ \kappa_{aaa} = (A_{aaa,0} + A_{aaa,1}\psi_1 + A_{aaa,2}\psi_2)\tau^{N_{aaa}} \f]
     *
     * @param HEOS The instance to use
     * @return
     */
    static CoolPropDbl viscosity_higher_order_friction_theory(HelmholtzEOSMixtureBackend& HEOS);

    /**
     * Implement the method of:
     *
     * Chung, Ting Horng, et al. "Generalized multiparameter correlation for nonpolar and polar fluid transport properties."
     * Industrial & engineering chemistry research 27(4) (1988): 671-679.
     */
    static CoolPropDbl viscosity_Chung(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief The general dilute gas conductivity term formed of a ratio of polynomial like terms

    \f[
    \lambda^0 = \frac{A_i\displaystyle\sum_iT_r^{n_i}}{B_i\displaystyle\sum_iT_r^{m_i}}
    \f]
    with \f$\lambda^0\f$ in W/m/K, T_r is the reduced temperature \f$T_{r} = T/T_{red}\f$
    */
    static CoolPropDbl conductivity_dilute_ratio_polynomials(HelmholtzEOSMixtureBackend& HEOS);

    /**

    This term is given by
    \f[
    \Delta\lambda(\rho,T) = \displaystyle\sum_iA_i\tau^{t,i}\delta^{d_i}
    \f]

    As used by Assael, Perkins, Huber, etc., the residual term is given by
    \f[
    \Delta\lambda(\rho,T) = \displaystyle\sum_i(B_{1,i}+B_{2,i}(T/T_c))(\rho/\rho_c)^i
    \f]
    which can be easily converted by noting that \f$\tau=Tc/T\f$ and \f$\delta=\rho/\rho_c\f$
    */
    static CoolPropDbl conductivity_residual_polynomial(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief The simplified critical conductivity term of Olchowy and Sengers

    Olchowy, G. A. & Sengers, J. V. (1989), "A Simplified Representation for the Thermal Conductivity of Fluids in the Critical Region", International Journal of Thermophysics, 10, (2), 417-426

    \f[
        \lambda^{(c)} = \frac{\rho c_p R_DkT}{6\pi\eta\zeta}(\Omega-\Omega_0)
    \f]
    \f[
        \Omega = \frac{2}{\pi}\left[ \left( \frac{c_p-c_v}{c_p}\right)\arctan(q_d\zeta)+\frac{c_v}{c_p}q_d\zeta \right]
    \f]
    \f[
        \Omega_0 = \frac{2}{\pi}\left[1-\exp\left(-\frac{1}{(q_d\zeta)^{-1}+(q_d\zeta\rho_c/\rho)^2/3} \right) \right]
    \f]
    \f[
        \zeta = \zeta_0\left(\frac{p_c\rho}{\Gamma\rho_c^2}\right)^{\nu/\gamma}\left[\left.\frac{\partial \rho(T,\rho)}{\partial p} \right|_{T}- \frac{T_R}{T}\left.\frac{\partial \rho(T_R,\rho)}{\partial p} \right|_{T}  \right]^{\nu/\gamma},
    \f]
    where \f$\lambda^{(c)}\f$ is in W\f$\cdot\f$m\f$^{-1}\f$\f$\cdot\f$K\f$^{-1}\f$, \f$\zeta\f$ is in m,
    \f$c_p\f$ and \f$c_v\f$ are in J\f$\cdot\f$kg\f$^{-1}\cdot\f$K\f$^{-1}\f$, \f$p\f$ and \f$p_c\f$ are in Pa,
    \f$\rho\f$ and \f$\rho_c\f$ are in mol\f$\cdot\f$m\f$^{-3}\f$, \f$\eta\f$ is the viscosity in Pa\f$\cdot\f$s,
    and the remaining parameters are defined in the following tables.

    It should be noted that some authors use slightly different values for the "universal" constants

    Coefficients for use in the simplified Olchowy-Sengers critical term
    Parameter             | Variable     | Value
    ---------             | --------     | ------
    Boltzmann constant    | \f$k\f$      | \f$1.3806488\times 10^{-23}\f$ J\f$\cdot\f$K\f$^{-1}\f$
    Universal amplitude   | \f$R_D\f$    | 1.03
    Critical exponent     | \f$\nu\f$    | 0.63
    Critical exponent     | \f$\gamma\f$ | 1.239
    Reference temperature | \f$T_R\f$    | 1.5\f$T_c\f$

    Recommended default constants (see Huber (I&ECR, 2003))
    Parameter        | Variable     | Value
    ---------        | --------     | ------
    Amplitude        | \f$\Gamma\f$ | 0.0496
    Amplitude        |\f$\zeta_0\f$ | 1.94 \f$\times\f$ 10\f$^{-10}\f$ m
    Effective cutoff | \f$q_d\f$    | 2 \f$\times\f$ 10\f$^{9}\f$ m

    */
    static CoolPropDbl conductivity_critical_simplified_Olchowy_Sengers(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl conductivity_critical_hardcoded_CO2_ScalabrinJPCRD2006(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_critical_hardcoded_R123(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_dilute_hardcoded_CO2(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_dilute_hardcoded_ethane(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl conductivity_dilute_eta0_and_poly(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_residual_polynomial_and_exponential(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl conductivity_hardcoded_heavywater(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_hardcoded_water(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_hardcoded_R23(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_hardcoded_helium(HelmholtzEOSMixtureBackend& HEOS);
    static CoolPropDbl conductivity_hardcoded_methane(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl conductivity_critical_hardcoded_ammonia(HelmholtzEOSMixtureBackend& HEOS);

    /**
    \brief Calculate the viscosity using the extended corresponding states method

    This method is covered in depth in

    Bell, I. H.; Wronski, J.; Quoilin, S. & Lemort, V. (2014), Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp, Industrial & Engineering Chemistry Research, 53, (6), 2498-2508

    which is originally based on the methods presented in

    Huber, M. L., Laesecke, A. and Perkins, R. A., (2003), Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a, Industrial & Engineering Chemistry Research, v. 42, pp. 3163-3178

    and

    McLinden, M. O.; Klein, S. A. & Perkins, R. A. (2000), An extended corresponding states model for the thermal conductivity of refrigerants and refrigerant mixtures, Int. J. Refrig., 23, 43-63


    */
    static CoolPropDbl viscosity_ECS(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference);

    static CoolPropDbl viscosity_rhosr(HelmholtzEOSMixtureBackend& HEOS);

    static CoolPropDbl conductivity_ECS(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference);

    /* \brief Solver for the conformal state for ECS model
     *
     */
    static void conformal_state_solver(HelmholtzEOSMixtureBackend& HEOS, HelmholtzEOSMixtureBackend& HEOS_Reference, CoolPropDbl& T0,
                                       CoolPropDbl& rhomolar0);

}; /* class TransportRoutines */

}; /* namespace CoolProp */
#endif
