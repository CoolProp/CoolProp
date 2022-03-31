#include "MixtureDerivatives.h"
#include "Backends/Cubics/CubicBackend.h"

namespace CoolProp {

CoolPropDbl MixtureDerivatives::dln_fugacity_i_dT__constp_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return dln_fugacity_coefficient_dT__constp_n(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dp__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return dln_fugacity_coefficient_dp__constT_n(HEOS, i, xN_flag) + 1 / HEOS.p();
}
CoolPropDbl MixtureDerivatives::fugacity_i(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.mole_fractions[i] * HEOS.rhomolar() * HEOS.gas_constant() * HEOS.T() * exp(dnalphar_dni__constT_V_nj(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::ln_fugacity_coefficient(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag) - log(1 + HEOS._delta.pt() * HEOS.dalphar_dDelta());
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dT__constrho_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return 1 / HEOS.T() * (1 - HEOS.tau() * HEOS.dalphar_dTau() - HEOS.tau() * d_ndalphardni_dTau(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_drho__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return 1 / HEOS.rhomolar() * (1 + HEOS.delta() * HEOS.dalphar_dDelta() + HEOS.delta() * d_ndalphardni_dDelta(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dnalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // GERG Equation 7.42
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2nalphar_dni_dT(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return -HEOS._tau.pt() / HEOS._T * (HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double T = HEOS._reducing.T / HEOS._tau.pt();
    CoolPropDbl R_u = HEOS.gas_constant();
    return d2nalphar_dni_dT(HEOS, i, xN_flag) + 1 / T - partial_molar_volume(HEOS, i, xN_flag) / (R_u * T) * dpdT__constV_n(HEOS);
}
CoolPropDbl MixtureDerivatives::partial_molar_volume(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return -ndpdni__constT_V_nj(HEOS, i, xN_flag) / ndpdV__constT_n(HEOS);
}

CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // GERG equation 7.30
    CoolPropDbl R_u = HEOS.gas_constant();
    double partial_molar_volumeval = partial_molar_volume(HEOS, i, xN_flag);  // [m^3/mol]
    double term1 = partial_molar_volumeval / (R_u * HEOS._T);                 // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
    double term2 = 1.0 / HEOS.p();
    return term1 - term2;
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dtau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return -1 / HEOS.tau() + HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_ddelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return 1 + HEOS.delta() * HEOS.dalphar_dDelta() + HEOS.delta() * d_ndalphardni_dDelta(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                              x_N_dependency_flag xN_flag) {
    // This is a term to which some more might be added depending on i and j
    CoolPropDbl val = dln_fugacity_coefficient_dxj__constT_p_xi(HEOS, i, j, xN_flag);
    const std::vector<CoolPropDbl>& x = HEOS.get_mole_fractions();
    std::size_t N = x.size();
    if (i == N - 1) {
        val += -1 / x[N - 1];
    } else if (i == j) {
        val += 1 / x[j];
    }
    return val;
}
CoolPropDbl MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                x_N_dependency_flag xN_flag) {
    if (xN_flag == XN_INDEPENDENT) {
        throw ValueError("dln_fugacity_dxj__constT_rho_xi only valid for xN_DEPENDENT for now");
    }
    CoolPropDbl rhor = HEOS.Reducing->rhormolar(HEOS.get_mole_fractions());
    CoolPropDbl Tr = HEOS.Reducing->Tr(HEOS.get_mole_fractions());
    CoolPropDbl dTrdxj = HEOS.Reducing->dTrdxi__constxj(HEOS.get_mole_fractions(), j, xN_flag);
    CoolPropDbl drhordxj = HEOS.Reducing->drhormolardxi__constxj(HEOS.get_mole_fractions(), j, xN_flag);

    // These lines are all the same
    CoolPropDbl line1 = dln_fugacity_i_dtau__constdelta_x(HEOS, i, xN_flag) * 1 / HEOS.T() * dTrdxj;
    CoolPropDbl line2 = -dln_fugacity_i_ddelta__consttau_x(HEOS, i, xN_flag) * 1 / rhor * drhordxj;
    CoolPropDbl line4 = HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag) + d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag);

    const std::vector<CoolPropDbl>& x = HEOS.get_mole_fractions();
    std::size_t N = x.size();

    CoolPropDbl line3 = 1 / rhor * HEOS.Reducing->drhormolardxi__constxj(x, j, xN_flag) + 1 / Tr * HEOS.Reducing->dTrdxi__constxj(x, j, xN_flag);
    ;

    // This is a term to which some more might be added depending on i and j
    if (i == N - 1) {
        line3 += -1 / x[N - 1];
    } else if (i == j) {
        line3 += 1 / x[j];
    } else {
    }

    return line1 + line2 + line3 + line4;
}

CoolPropDbl MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag) {
    double s = (HEOS.mole_fractions[i] > DBL_EPSILON) ? Kronecker_delta(i, j) / HEOS.mole_fractions[i] : 0;
    return s + nd2nalphardnidnj__constT_V(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                         x_N_dependency_flag xN_flag) {
    return d_ndalphardni_dTau(HEOS, j, xN_flag) + d_nd_ndalphardni_dnj_dTau__constdelta_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           x_N_dependency_flag xN_flag) {
    return d2_ndalphardni_dTau2(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dTau2__constdelta_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_ddelta2__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           x_N_dependency_flag xN_flag) {
    return d2_ndalphardni_dDelta2(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dDelta2__consttau_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           x_N_dependency_flag xN_flag) {
    return d2_ndalphardni_dDelta_dTau(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dDelta_dTau__constx(HEOS, i, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                         x_N_dependency_flag xN_flag) {
    return d_ndalphardni_dDelta(HEOS, j, xN_flag) + d_nd_ndalphardni_dnj_dDelta__consttau_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           std::size_t k, x_N_dependency_flag xN_flag) {
    CoolPropDbl s = (HEOS.mole_fractions[i] > DBL_EPSILON) ? -Kronecker_delta(i, j) * Kronecker_delta(i, k) / pow(HEOS.mole_fractions[i], 2) : 0;
    return s + d_ndalphardni_dxj__constdelta_tau_xi(HEOS, j, k, xN_flag) + d_nd_ndalphardni_dnj_dxk__consttau_delta(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                            std::size_t k, x_N_dependency_flag xN_flag) {
    return d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, j, k, xN_flag) + d2_nd_ndalphardni_dnj_dxk_dTau__constdelta(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                            std::size_t k, x_N_dependency_flag xN_flag) {
    return d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, j, k, xN_flag) + d2_nd_ndalphardni_dnj_dxk_dDelta__consttau(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::nd_ndln_fugacity_i_dnj_dnk__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                                        x_N_dependency_flag xN_flag) {
    double sum = d_ndln_fugacity_i_dnj_dtau__constdelta_x(HEOS, i, j, xN_flag) * ndtaudni__constT_V_nj(HEOS, k, xN_flag)
                 + d_ndln_fugacity_i_dnj_ddelta__consttau_x(HEOS, i, j, xN_flag) * nddeltadni__constT_V_nj(HEOS, k, xN_flag)
                 + d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HEOS, i, j, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; ++m) {
        sum -= HEOS.mole_fractions[m] * d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HEOS, i, j, m, xN_flag);
    }
    return sum;
}

CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                          x_N_dependency_flag xN_flag) {
    // Gernert 3.115
    CoolPropDbl R_u = HEOS.gas_constant();
    // partial molar volume is -dpdn/dpdV, so need to flip the sign here
    return d2nalphar_dxj_dni__constT_V(HEOS, j, i, xN_flag)
           - partial_molar_volume(HEOS, i, xN_flag) / (R_u * HEOS._T) * dpdxj__constT_V_xi(HEOS, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::dpdxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag) {
    // Gernert 3.130
    CoolPropDbl R_u = HEOS.gas_constant();
    return HEOS._rhomolar * R_u * HEOS._T
           * (ddelta_dxj__constT_V_xi(HEOS, j, xN_flag) * HEOS.dalphar_dDelta()
              + HEOS._delta.pt() * d_dalpharddelta_dxj__constT_V_xi(HEOS, j, xN_flag));
}

CoolPropDbl MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag) {
    // Gernert Equation 3.134 (Catch test provided)
    return HEOS.d2alphar_dDelta2() * ddelta_dxj__constT_V_xi(HEOS, j, xN_flag) + HEOS.d2alphar_dDelta_dTau() * dtau_dxj__constT_V_xi(HEOS, j, xN_flag)
           + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag) {
    //Gernert 3.119 (Catch test provided)
    return HEOS.dalphar_dDelta() * ddelta_dxj__constT_V_xi(HEOS, j, xN_flag) + HEOS.dalphar_dTau() * dtau_dxj__constT_V_xi(HEOS, j, xN_flag)
           + HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag) {
    // Gernert 3.118
    return d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag)
           + ddelta_dxj__constT_V_xi(HEOS, j, xN_flag) * d_ndalphardni_dDelta(HEOS, i, xN_flag)
           + dtau_dxj__constT_V_xi(HEOS, j, xN_flag) * d_ndalphardni_dTau(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::ddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag) {
    // Gernert 3.121 (Catch test provided)
    return -HEOS._delta.pt() / HEOS._reducing.rhomolar * HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::dtau_dxj__constT_V_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t j, x_N_dependency_flag xN_flag) {
    // Gernert 3.122 (Catch test provided)
    return 1 / HEOS._T * HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::dpdT__constV_n(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl R_u = HEOS.gas_constant();
    return HEOS._rhomolar * R_u * (1 + HEOS._delta.pt() * HEOS.dalphar_dDelta() - HEOS._delta.pt() * HEOS._tau.pt() * HEOS.d2alphar_dDelta_dTau());
}
CoolPropDbl MixtureDerivatives::dpdrho__constT_n(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl R_u = HEOS.gas_constant();
    return R_u * HEOS._T * (1 + 2 * HEOS._delta.pt() * HEOS.dalphar_dDelta() + pow(HEOS._delta.pt(), 2) * HEOS.d2alphar_dDelta2());
}
CoolPropDbl MixtureDerivatives::ndpdV__constT_n(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl R_u = HEOS.gas_constant();
    return -pow(HEOS._rhomolar, 2) * R_u * HEOS._T
           * (1 + 2 * HEOS._delta.pt() * HEOS.dalphar_dDelta() + pow(HEOS._delta.pt(), 2) * HEOS.d2alphar_dDelta2());
}
CoolPropDbl MixtureDerivatives::ndpdni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // Eqn 7.64 and 7.63
    CoolPropDbl R_u = HEOS.gas_constant();
    double ndrhorbar_dni__constnj = HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
    double ndTr_dni__constnj = HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; ++k) {
        summer += HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    double nd2alphar_dni_dDelta = HEOS._delta.pt() * HEOS.d2alphar_dDelta2() * (1 - 1 / HEOS._reducing.rhomolar * ndrhorbar_dni__constnj)
                                  + HEOS._tau.pt() * HEOS.d2alphar_dDelta_dTau() / HEOS._reducing.T * ndTr_dni__constnj
                                  + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, i, xN_flag) - summer;
    return HEOS._rhomolar * R_u * HEOS._T
           * (1 + HEOS._delta.pt() * HEOS.dalphar_dDelta() * (2 - 1 / HEOS._reducing.rhomolar * ndrhorbar_dni__constnj)
              + HEOS._delta.pt() * nd2alphar_dni_dDelta);
}

CoolPropDbl MixtureDerivatives::ndalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 = HEOS._delta.pt() * HEOS.dalphar_dDelta()
                   * (1 - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double term2 = HEOS._tau.pt() * HEOS.dalphar_dTau() * (1 / HEOS._reducing.T) * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        s += HEOS.mole_fractions[k] * HEOS.residual_helmholtz->dalphar_dxi(HEOS, k, xN_flag);
    }
    double term3 = HEOS.residual_helmholtz->dalphar_dxi(HEOS, i, xN_flag);
    return term1 + term2 + term3 - s;
}
CoolPropDbl MixtureDerivatives::ndln_fugacity_coefficient_dnj__constT_p(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    CoolPropDbl R_u = HEOS.gas_constant();
    return nd2nalphardnidnj__constT_V(HEOS, j, i, xN_flag) + 1
           - partial_molar_volume(HEOS, j, xN_flag) / (R_u * HEOS._T) * ndpdni__constT_V_nj(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::nddeltadni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS._delta.pt() - HEOS._delta.pt() / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_nddeltadni_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return 1 - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_nddeltadni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag) {
    double rhor = HEOS._reducing.rhomolar;
    return -HEOS.delta() / rhor
           * (HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)
              - 1 / rhor * HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag)
                  * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::d2_nddeltadni_dxj_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                   x_N_dependency_flag xN_flag) {
    return d_nddeltadni_dxj__constdelta_tau(HEOS, i, j, xN_flag) / HEOS.delta();
}

CoolPropDbl MixtureDerivatives::ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS._tau.pt() / HEOS._reducing.T * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}

CoolPropDbl MixtureDerivatives::d_ndtaudni_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return 1 / HEOS._reducing.T * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndtaudni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                               x_N_dependency_flag xN_flag) {
    double Tr = HEOS._reducing.T;
    return HEOS.tau() / Tr
           * (HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)
              - 1 / Tr * HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)
                  * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::d2_ndtaudni_dxj_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                 x_N_dependency_flag xN_flag) {
    return d_ndtaudni_dxj__constdelta_tau(HEOS, i, j, xN_flag) / HEOS.tau();
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                     x_N_dependency_flag xN_flag) {
    double line1 = HEOS._delta.pt() * HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag)
                   * (1 - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double line3 = HEOS._tau.pt() * HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag) * (1 / HEOS._reducing.T)
                   * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
    double line2 = -HEOS._delta.pt() * HEOS.dalphar_dDelta() * (1 / HEOS._reducing.rhomolar)
                   * (HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)
                      - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag)
                          * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double line4 = HEOS._tau.pt() * HEOS.dalphar_dTau() * (1 / HEOS._reducing.T)
                   * (HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)
                      - 1 / HEOS._reducing.T * HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)
                          * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));

    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        s += HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d2alphardxidxj(HEOS, j, k, xN_flag);
    }
    double line5 = HEOS.residual_helmholtz->d2alphardxidxj(HEOS, i, j, xN_flag) - HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag) - s;
    return line1 + line2 + line3 + line4 + line5;
}

CoolPropDbl MixtureDerivatives::nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                           x_N_dependency_flag xN_flag) {
    return ndalphar_dni__constT_V_nj(HEOS, j, xN_flag)           // First term from 7.46
           + nd_ndalphardni_dnj__constT_V(HEOS, i, j, xN_flag);  // 7.47
}

CoolPropDbl MixtureDerivatives::nd_ndalphardni_dnj__constT_V(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                             x_N_dependency_flag xN_flag) {
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double line3 = d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        line3 -= HEOS.mole_fractions[k] * d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dTau__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    double line1 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dTau2(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag)
                   + d_ndalphardni_dTau(HEOS, i, xN_flag) * d_ndtaudni_dTau(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        summer += HEOS.mole_fractions[k] * d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag);
    }
    double line3 = d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, j, xN_flag) - summer;
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dTau2__constdelta_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                          x_N_dependency_flag xN_flag) {
    double line1 = d3_ndalphardni_dDelta_dTau2(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = 2 * d2_ndalphardni_dTau2(HEOS, i, xN_flag) * d_ndtaudni_dTau(HEOS, j, xN_flag);
    double line3 = d3_ndalphardni_dTau3(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        summer += HEOS.mole_fractions[k] * d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, k, xN_flag);
    }
    double line4 = d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, j, xN_flag) - summer;
    return line1 + line2 + line3 + line4;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dDelta2__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                          x_N_dependency_flag xN_flag) {
    double line1 = d3_ndalphardni_dDelta3(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = 2 * d2_ndalphardni_dDelta2(HEOS, i, xN_flag) * d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line3 = d3_ndalphardni_dDelta2_dTau(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        summer += HEOS.mole_fractions[k] * d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, k, xN_flag);
    }
    double line4 = d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, j, xN_flag) - summer;
    return line1 + line2 + line3 + line4;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dDelta_dTau__constx(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                          x_N_dependency_flag xN_flag) {
    double line1 = d3_ndalphardni_dDelta2_dTau(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line3 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * d_ndtaudni_dTau(HEOS, j, xN_flag);
    double line4 = d3_ndalphardni_dDelta_dTau2(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        summer += HEOS.mole_fractions[k] * d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag);
    }
    double line5 = d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, j, xN_flag) - summer;
    return line1 + line2 + line3 + line4 + line5;
}
CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dDelta__consttau_x(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    double line1 = d2_ndalphardni_dDelta2(HEOS, i, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag)
                   + d_ndalphardni_dDelta(HEOS, i, xN_flag) * d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        summer += HEOS.mole_fractions[k] * d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag);
    }
    double line3 = d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, j, xN_flag) - summer;
    return line1 + line2 + line3;
}

CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dxk__consttau_delta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                         std::size_t k, x_N_dependency_flag xN_flag) {
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag) * d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);

    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag) * d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HEOS, i, j, k, xN_flag) - d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k, xN_flag);

    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        line3 -= HEOS.mole_fractions[m] * d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           std::size_t k, x_N_dependency_flag xN_flag) {
    double line1 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);

    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag) * d2_ndtaudni_dxj_dTau__constdelta(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dTau2(HEOS, i, xN_flag) * d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag) * d_ndtaudni_dTau(HEOS, j, xN_flag)
                   + d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, k, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HEOS, i, j, k, xN_flag) - d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag);

    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        line3 -= HEOS.mole_fractions[m] * d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           std::size_t k, x_N_dependency_flag xN_flag) {
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag) * d2_nddeltadni_dxj_dDelta__consttau(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dDelta2(HEOS, i, xN_flag) * d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag) * d_nddeltadni_dDelta(HEOS, j, xN_flag)
                   + d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, k, xN_flag) * nddeltadni__constT_V_nj(HEOS, j, xN_flag);

    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag) * d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                   + d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag) * ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HEOS, i, j, k, xN_flag) - d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        line3 -= HEOS.mole_fractions[m] * d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // The first line
    double term1 = (HEOS._delta.pt() * HEOS.d2alphar_dDelta2() + HEOS.dalphar_dDelta())
                   * (1 - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 =
      HEOS._tau.pt() * HEOS.d2alphar_dDelta_dTau() * (1 / HEOS._reducing.T) * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 = (2 * HEOS.d2alphar_dDelta2() + HEOS.delta() * HEOS.d3alphar_dDelta3()) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = HEOS.tau() * HEOS.d3alphar_dDelta2_dTau() * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 = (3 * HEOS.d3alphar_dDelta3() + HEOS.delta() * HEOS.d4alphar_dDelta4()) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = HEOS.tau() * HEOS.d4alphar_dDelta3_dTau() * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 =
      (HEOS.d2alphar_dDelta_dTau() + HEOS.delta() * HEOS.d3alphar_dDelta2_dTau()) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (HEOS.tau() * HEOS.d3alphar_dDelta_dTau2() + HEOS.d2alphar_dDelta_dTau()) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta2_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 =
      (2 * HEOS.d3alphar_dDelta2_dTau() + HEOS.delta() * HEOS.d4alphar_dDelta3_dTau()) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 =
      (HEOS.tau() * HEOS.d4alphar_dDelta2_dTau2() + HEOS.d3alphar_dDelta2_dTau()) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 =
      (HEOS.d3alphar_dDelta_dTau2() + HEOS.delta() * HEOS.d4alphar_dDelta2_dTau2()) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 =
      (HEOS.tau() * HEOS.d4alphar_dDelta_dTau3() + 2 * HEOS.d3alphar_dDelta_dTau2()) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 = HEOS.delta() * HEOS.d3alphar_dDelta_dTau2() * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (2 * HEOS.d2alphar_dTau2() + HEOS.tau() * HEOS.d3alphar_dTau3()) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dTau3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    double term1 = HEOS.delta() * HEOS.d4alphar_dDelta_dTau3() * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (3 * HEOS.d3alphar_dTau3() + HEOS.tau() * HEOS.d4alphar_dTau4()) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dDelta__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                       x_N_dependency_flag xN_flag) {
    double term1 =
      (HEOS.dalphar_dDelta() + HEOS.delta() * HEOS.d2alphar_dDelta2()) * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag)
                    + HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag))
                   * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.tau() * HEOS.d2alphar_dDelta_dTau() * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 =
      HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term5 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, k, j, xN_flag)
                 + Kronecker_delta(k, j) * HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dTau__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                       x_N_dependency_flag xN_flag) {
    double term1 = HEOS.delta() * HEOS.d2alphar_dDelta_dTau() * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 =
      HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = (HEOS.tau() * HEOS.d2alphar_dTau2() + HEOS.dalphar_dTau()) * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 =
      (HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag))
      * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term5 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, k, j, xN_flag)
                 + Kronecker_delta(k, j) * HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dTau2__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    double term1 = HEOS.delta() * HEOS.d3alphar_dDelta_dTau2() * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 =
      HEOS.delta() * HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, j, xN_flag) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = (HEOS.tau() * HEOS.d3alphar_dTau3() + 2 * HEOS.d2alphar_dTau2()) * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 =
      (HEOS.tau() * HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, j, xN_flag) + 2 * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag))
      * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term5 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, k, j, xN_flag)
                 + Kronecker_delta(k, j) * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dDelta2__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    double term1 =
      (2 * HEOS.d2alphar_dDelta2() + HEOS.delta() * HEOS.d3alphar_dDelta3()) * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (2 * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag)
                    + HEOS.delta() * HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, j, xN_flag))
                   * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.tau() * HEOS.d3alphar_dDelta2_dTau() * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 =
      HEOS.tau() * HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, j, xN_flag) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term5 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, k, j, xN_flag)
                 + Kronecker_delta(k, j) * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dDelta_dTau__constxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                        x_N_dependency_flag xN_flag) {
    double term1 =
      (HEOS.d2alphar_dDelta_dTau() + HEOS.delta() * HEOS.d3alphar_dDelta2_dTau()) * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)
                    + HEOS.delta() * HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, j, xN_flag))
                   * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 =
      (HEOS.tau() * HEOS.d3alphar_dDelta_dTau2() + HEOS.d2alphar_dDelta_dTau()) * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = (HEOS.tau() * HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, j, xN_flag)
                    + HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag))
                   * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term5 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, k, j, xN_flag)
                 + Kronecker_delta(k, j) * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}

CoolPropDbl MixtureDerivatives::d_ndalphardni_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // The first line
    double term1 = HEOS._delta.pt() * HEOS.d2alphar_dDelta_dTau()
                   * (1 - 1 / HEOS._reducing.rhomolar * HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = (HEOS._tau.pt() * HEOS.d2alphar_dTau2() + HEOS.dalphar_dTau()) * (1 / HEOS._reducing.T)
                   * HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (unsigned int k = 0; k < kmax; k++) {
        term3 -= HEOS.mole_fractions[k] * HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                          std::size_t k, x_N_dependency_flag xN_flag) {
    double term1 =
      HEOS.delta()
      * (HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag) * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag)
         + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag) * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag));
    double term2 =
      HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag) * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.delta() * HEOS.dalphar_dDelta() * HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);

    double term4 =
      HEOS.tau()
      * (HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag) * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag)
         + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag) * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag));
    double term5 =
      HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term6 = HEOS.tau() * HEOS.dalphar_dTau() * HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);

    double term7 =
      HEOS.residual_helmholtz->d3alphardxidxjdxk(HEOS, i, j, k, xN_flag) - 2 * HEOS.residual_helmholtz->d2alphardxidxj(HEOS, j, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        term7 -= HEOS.mole_fractions[m] * HEOS.residual_helmholtz->d3alphardxidxjdxk(HEOS, j, k, m, xN_flag);
    }

    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           std::size_t k, x_N_dependency_flag xN_flag) {
    double term1a = HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)
                    * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term1b = HEOS.delta() * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, j, k, xN_flag)
                    * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term1c = HEOS.delta() * HEOS.d2alphar_dDelta_dTau() * HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term1d = HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag)
                    * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term1 = term1a + term1b + term1c + term1d;

    double term2a =
      (HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag))
      * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term2b = (HEOS.tau() * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, j, k, xN_flag)
                     + HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag))
                    * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term2c =
      (HEOS.tau() * HEOS.d2alphar_dTau2() + HEOS.dalphar_dTau()) * HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term2d =
      (HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag))
      * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = term2a + term2b + term2c + term2d;

    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dxk_dTau(HEOS, i, j, k, xN_flag)
                   - 2 * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        term3 -= HEOS.mole_fractions[m] * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dxk_dTau(HEOS, j, k, m, xN_flag);
    }

    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                                                           std::size_t k, x_N_dependency_flag xN_flag) {
    double term1a = (HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag)
                     + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag))
                    * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term1b = (HEOS.delta() * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, j, k, xN_flag)
                     + HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag))
                    * HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term1c =
      (HEOS.delta() * HEOS.d2alphar_dDelta2() + HEOS.dalphar_dDelta()) * HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term1d = (HEOS.delta() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag)
                     + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag))
                    * HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term1 = term1a + term1b + term1c + term1d;

    double term2a = HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)
                    * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term2b =
      HEOS.tau() * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, j, k, xN_flag) * HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term2c = HEOS.tau() * HEOS.d2alphar_dDelta_dTau() * HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term2d = HEOS.tau() * HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag)
                    * HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = term2a + term2b + term2c + term2d;

    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dxk_dDelta(HEOS, i, j, k, xN_flag)
                   - 2 * HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        mmax--;
    }
    for (unsigned int m = 0; m < mmax; m++) {
        term3 -= HEOS.mole_fractions[m] * HEOS.residual_helmholtz->d4alphar_dxi_dxj_dxk_dDelta(HEOS, j, k, m, xN_flag);
    }
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::dalpha0_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // Reducing values are constant for all components under consideration
    double Tr = HEOS.T_reducing();
    double rhor = HEOS.rhomolar_reducing();

    // Values for the i-th component
    double Tci = HEOS.get_fluid_constant(i, iT_critical);
    double rhoci = HEOS.get_fluid_constant(i, irhomolar_critical);
    double tau_oi = HEOS.tau() * Tci / Tr;
    double delta_oi = HEOS.delta() * rhor / rhoci;
    double Rratioi = 1;  //HEOS.gas_constant()/HEOS.components[i].EOS().R_u;

    double logxi = (std::abs(HEOS.mole_fractions[i]) > DBL_EPSILON) ? (log(HEOS.mole_fractions[i])) : 0;
    double term = Rratioi * HEOS.components[i].EOS().alpha0.base(tau_oi, delta_oi) + logxi + 1;

    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (std::size_t k = 0; k < kmax; ++k) {
        double xk = HEOS.mole_fractions[k];
        double Tck = HEOS.get_fluid_constant(k, iT_critical);
        double rhock = HEOS.get_fluid_constant(k, irhomolar_critical);
        double tau_ok = HEOS.tau() * Tck / Tr;
        double delta_ok = HEOS.delta() * rhor / rhock;
        double dtauok_dxi = -tau_ok / Tr * HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag);              // (Gernert, supp, B.19)
        double ddeltaok_dxi = delta_ok / rhor * HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag);  // (Gernert, supp. B.20)

        double Rratiok = 1;  //HEOS.gas_constant()/HEOS.components[k].EOS().R_u;
        HelmholtzDerivatives alpha0kterms = HEOS.components[k].EOS().alpha0.all(tau_ok, delta_ok);
        double dalpha0_ok_dxi = alpha0kterms.dalphar_dtau * dtauok_dxi + alpha0kterms.dalphar_ddelta * ddeltaok_dxi;
        term += xk * (Rratiok * dalpha0_ok_dxi);
    }
    return term;
}

CoolPropDbl MixtureDerivatives::d2alpha0_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // Reducing values are constant for all components under consideration
    double Tr = HEOS.T_reducing();
    double rhor = HEOS.rhomolar_reducing();

    // Values for the i-th component
    double Tci = HEOS.get_fluid_constant(i, iT_critical);
    double rhoci = HEOS.get_fluid_constant(i, irhomolar_critical);
    double tau_oi = HEOS.tau() * Tci / Tr;
    double delta_oi = HEOS.delta() * rhor / rhoci;
    double Rratioi = 1;  //HEOS.gas_constant()/HEOS.components[i].EOS().R_u;

    double term = rhor / rhoci * Rratioi * HEOS.components[i].EOS().alpha0.dDelta(tau_oi, delta_oi);

    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (std::size_t k = 0; k < kmax; ++k) {
        double xk = HEOS.mole_fractions[k];
        double Tck = HEOS.get_fluid_constant(k, iT_critical);
        double rhock = HEOS.get_fluid_constant(k, irhomolar_critical);
        double tau_ok = HEOS.tau() * Tck / Tr;
        double delta_ok = HEOS.delta() * rhor / rhock;
        double dtauok_dxi = -tau_ok / Tr * HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag);  // (Gernert, supp, B.19)
        double drhor_dxi = HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag);
        double ddeltaok_dxi = delta_ok / rhor * drhor_dxi;  // (Gernert, supp. B.20)

        //double Rratiok = 1;//HEOS.gas_constant()/HEOS.components[k].EOS().R_u;
        HelmholtzDerivatives alpha0kterms = HEOS.components[k].EOS().alpha0.all(tau_ok, delta_ok);
        double dalpha0ok_ddeltaok = alpha0kterms.dalphar_ddelta;

        double d_dalpha0ok_ddeltaok_dxi = alpha0kterms.d2alphar_ddelta_dtau * dtauok_dxi + alpha0kterms.d2alphar_ddelta2 * ddeltaok_dxi;
        term += xk / rhock * (rhor * d_dalpha0ok_ddeltaok_dxi + drhor_dxi * dalpha0ok_ddeltaok);
    }
    return term;
}

CoolPropDbl MixtureDerivatives::d2alpha0_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    // Reducing values are constant for all components under consideration
    double Tr = HEOS.T_reducing();
    double rhor = HEOS.rhomolar_reducing();

    // Values for the i-th component
    double Tci = HEOS.get_fluid_constant(i, iT_critical);
    double rhoci = HEOS.get_fluid_constant(i, irhomolar_critical);
    double tau_oi = HEOS.tau() * Tci / Tr;
    double delta_oi = HEOS.delta() * rhor / rhoci;
    double Rratioi = 1;  //HEOS.gas_constant()/HEOS.components[i].EOS().R_u;

    double term = Tci / Tr * Rratioi * HEOS.components[i].EOS().alpha0.dTau(tau_oi, delta_oi);

    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (std::size_t k = 0; k < kmax; ++k) {
        double xk = HEOS.mole_fractions[k];
        double Tck = HEOS.get_fluid_constant(k, iT_critical);
        double rhock = HEOS.get_fluid_constant(k, irhomolar_critical);
        double tau_ok = HEOS.tau() * Tck / Tr;
        double delta_ok = HEOS.delta() * rhor / rhock;
        double dTr_dxi = HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag);
        double dtauok_dxi = -tau_ok / Tr * dTr_dxi;  // (Gernert, supp, B.19)
        double drhor_dxi = HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag);
        double ddeltaok_dxi = delta_ok / rhor * drhor_dxi;  // (Gernert, supp. B.20)

        //double Rratiok = 1;//HEOS.gas_constant()/HEOS.components[k].EOS().R_u;
        HelmholtzDerivatives alpha0kterms = HEOS.components[k].EOS().alpha0.all(tau_ok, delta_ok);
        double dalpha0ok_dtauok = alpha0kterms.dalphar_dtau;
        double d_dalpha0ok_dTauok_dxi = alpha0kterms.d2alphar_dtau2 * dtauok_dxi + alpha0kterms.d2alphar_ddelta_dtau * ddeltaok_dxi;
        term += xk * Tck * (1 / Tr * d_dalpha0ok_dTauok_dxi + -1 / POW2(Tr) * dTr_dxi * dalpha0ok_dtauok);
    }
    return term;
}

CoolPropDbl MixtureDerivatives::d2alpha0dxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
    // Reducing values are constant for all components under consideration
    double Tr = HEOS.T_reducing();
    double rhor = HEOS.rhomolar_reducing();

    // Values for the i-th component
    double Tci = HEOS.get_fluid_constant(i, iT_critical);
    double rhoci = HEOS.get_fluid_constant(i, irhomolar_critical);
    double tau_oi = HEOS.tau() * Tci / Tr;
    double delta_oi = HEOS.delta() * rhor / rhoci;
    double dTr_dxi = HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag);
    double drhor_dxi = HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag);

    // Values for the j-th component
    double Tcj = HEOS.get_fluid_constant(j, iT_critical);
    double rhocj = HEOS.get_fluid_constant(j, irhomolar_critical);
    double tau_oj = HEOS.tau() * Tcj / Tr;
    double delta_oj = HEOS.delta() * rhor / rhocj;
    double dTr_dxj = HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag);
    double drhor_dxj = HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag);

    // Cross-terms with both i & j
    double dtauoi_dxj = -tau_oi / Tr * dTr_dxj;         // (Gernert, supp, B.19)
    double ddeltaoi_dxj = delta_oi / rhor * drhor_dxj;  // (Gernert, supp. B.20)
    double dtauoj_dxi = -tau_oj / Tr * dTr_dxi;         // (Gernert, supp, B.19)
    double ddeltaoj_dxi = delta_oj / rhor * drhor_dxi;  // (Gernert, supp. B.20)
    double d2Tr_dxidxj = HEOS.Reducing->d2Trdxidxj(HEOS.mole_fractions, i, j, xN_flag);
    double d2rhor_dxidxj = HEOS.Reducing->d2rhormolardxidxj(HEOS.mole_fractions, i, j, xN_flag);

    //double Rratioi = 1;//HEOS.gas_constant()/HEOS.components[i].EOS().R_u;
    HelmholtzDerivatives alpha0iterms = HEOS.components[i].EOS().alpha0.all(tau_oi, delta_oi),
                         alpha0jterms = HEOS.components[j].EOS().alpha0.all(tau_oj, delta_oj);

    double d_dalpha0oi_dxj = alpha0iterms.dalphar_dtau * dtauoi_dxj + alpha0iterms.dalphar_ddelta * ddeltaoi_dxj;
    double d_dalpha0oj_dxi = alpha0jterms.dalphar_dtau * dtauoj_dxi + alpha0jterms.dalphar_ddelta * ddeltaoj_dxi;

    double xi = HEOS.mole_fractions[i], xj = HEOS.mole_fractions[j];
    double Kronecker_delta_over_xi = (xj > DBL_EPSILON && xi > DBL_EPSILON) ? Kronecker_delta(i, j) / xi : 0;
    double term = d_dalpha0oi_dxj + d_dalpha0oj_dxi + Kronecker_delta_over_xi;

    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT) {
        kmax--;
    }
    for (std::size_t k = 0; k < kmax; ++k) {
        // Values for the k-th component
        double xk = HEOS.mole_fractions[k];
        double Tck = HEOS.get_fluid_constant(k, iT_critical);
        double rhock = HEOS.get_fluid_constant(k, irhomolar_critical);
        double tau_ok = HEOS.tau() * Tck / Tr;
        double delta_ok = HEOS.delta() * rhor / rhock;

        double dtauok_dxj = -tau_ok / Tr * dTr_dxj;         // (Gernert, supp, B.19)
        double ddeltaok_dxj = delta_ok / rhor * drhor_dxj;  // (Gernert, supp. B.20)
        double dtauok_dxi = -tau_ok / Tr * dTr_dxi;         // (Gernert, supp, B.19)
        double ddeltaok_dxi = delta_ok / rhor * drhor_dxi;  // (Gernert, supp. B.20)

        HelmholtzDerivatives alpha0kterms = HEOS.components[k].EOS().alpha0.all(tau_ok, delta_ok);
        double dalpha0ok_dtauok = alpha0kterms.dalphar_dtau;
        double d2tauok_dxidxj = -Tck * HEOS.tau() * (POW2(Tr) * d2Tr_dxidxj - dTr_dxi * (2 * Tr * dTr_dxj)) / POW4(Tr);
        double d_dalpha0ok_dtauok_dxj = alpha0kterms.d2alphar_dtau2 * dtauok_dxj + alpha0kterms.d2alphar_ddelta_dtau * ddeltaok_dxj;

        double dalpha0ok_ddeltaok = alpha0kterms.dalphar_ddelta;
        double d2deltaok_dxidxj = HEOS.delta() / rhock * d2rhor_dxidxj;
        double d_dalpha0ok_ddeltaok_dxj = alpha0kterms.d2alphar_ddelta_dtau * dtauok_dxj + alpha0kterms.d2alphar_ddelta2 * ddeltaok_dxj;

        term += xk
                * (dalpha0ok_dtauok * d2tauok_dxidxj + d_dalpha0ok_dtauok_dxj * dtauok_dxi + dalpha0ok_ddeltaok * d2deltaok_dxidxj
                   + d_dalpha0ok_ddeltaok_dxj * ddeltaok_dxi);
    }
    return term;
}

/// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\delta\f$
CoolPropDbl MixtureDerivatives::dpsi_dDelta(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.gas_constant() * HEOS.T() * (HEOS.delta() * dalpha_dDelta(HEOS) + alpha(HEOS));
}

/// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\tau\f$
CoolPropDbl MixtureDerivatives::dpsi_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.delta() * HEOS.gas_constant() * HEOS.T() / HEOS.tau() * (HEOS.tau() * dalpha_dTau(HEOS) - alpha(HEOS));
}

/// Return the first partial of Helmholtz energy density \f$\psi = \rho a\f$ with respect to \f$\tau\f$
CoolPropDbl MixtureDerivatives::dpsir_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.delta() * HEOS.gas_constant() * HEOS.T() / HEOS.tau() * (HEOS.tau() * dalphar_dTau(HEOS) - alphar(HEOS));
}

CoolPropDbl MixtureDerivatives::dpsi_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / HEOS.tau()
           * (alpha(HEOS, xN_flag) * d_rhorTr_dxi(HEOS, i, xN_flag) + HEOS.rhomolar_reducing() * HEOS.T_reducing() * dalpha_dxi(HEOS, i, xN_flag));
}

CoolPropDbl MixtureDerivatives::dpsir_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / HEOS.tau()
           * (alphar(HEOS, xN_flag) * d_rhorTr_dxi(HEOS, i, xN_flag) + HEOS.rhomolar_reducing() * HEOS.T_reducing() * dalphar_dxi(HEOS, i, xN_flag));
}

CoolPropDbl MixtureDerivatives::d_rhorTr_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    GERG2008ReducingFunction* GERG = static_cast<GERG2008ReducingFunction*>(HEOS.Reducing.get());
    return HEOS.rhomolar_reducing() * GERG->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag)
           + HEOS.T_reducing() * GERG->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag);
}

CoolPropDbl MixtureDerivatives::d2_rhorTr_dxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
    GERG2008ReducingFunction* GERG = static_cast<GERG2008ReducingFunction*>(HEOS.Reducing.get());
    return (HEOS.rhomolar_reducing() * GERG->d2Trdxidxj(HEOS.mole_fractions, i, j, xN_flag)
            + GERG->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag) * GERG->dTrdxi__constxj(HEOS.mole_fractions, i, xN_flag)
            + HEOS.T_reducing() * GERG->d2rhormolardxidxj(HEOS.mole_fractions, i, j, xN_flag)
            + GERG->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag) * GERG->drhormolardxi__constxj(HEOS.mole_fractions, i, xN_flag));
}

CoolPropDbl MixtureDerivatives::d2psi_dDelta2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.gas_constant() * HEOS.T() * (HEOS.delta() * d2alpha_dDelta2(HEOS) + 2 * dalpha_dDelta(HEOS));
}

CoolPropDbl MixtureDerivatives::d2psi_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.gas_constant() * HEOS.T() / HEOS.tau()
           * (HEOS.tau() * dalpha_dTau(HEOS) - alpha(HEOS) - HEOS.delta() * dalpha_dDelta(HEOS)
              + HEOS.tau() * HEOS.delta() * d2alpha_dDelta_dTau(HEOS));
}
CoolPropDbl MixtureDerivatives::d2psir_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    return HEOS.rhomolar_reducing() * HEOS.gas_constant() * HEOS.T() / HEOS.tau()
           * (HEOS.tau() * dalphar_dTau(HEOS) - alphar(HEOS) - HEOS.delta() * dalphar_dDelta(HEOS)
              + HEOS.tau() * HEOS.delta() * d2alphar_dDelta_dTau(HEOS));
}
CoolPropDbl MixtureDerivatives::d2psi_dTau2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    double tau = HEOS.tau();
    return HEOS.rhomolar_reducing() * HEOS.delta() * HEOS.gas_constant() * HEOS.T() / POW2(tau)
           * (POW2(tau) * d2alpha_dTau2(HEOS) - 2 * tau * dalpha_dTau(HEOS) + 2 * alpha(HEOS));
}
CoolPropDbl MixtureDerivatives::d2psir_dTau2(HelmholtzEOSMixtureBackend& HEOS, x_N_dependency_flag xN_flag) {
    double tau = HEOS.tau();
    return HEOS.rhomolar_reducing() * HEOS.delta() * HEOS.gas_constant() * HEOS.T() / POW2(tau)
           * (POW2(tau) * d2alphar_dTau2(HEOS) - 2 * tau * dalphar_dTau(HEOS) + 2 * alphar(HEOS));
}
CoolPropDbl MixtureDerivatives::d2psi_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.gas_constant() / HEOS.tau()
           * (d_rhorTr_dxi(HEOS, i, xN_flag) * (HEOS.delta() * dalpha_dDelta(HEOS) + alpha(HEOS))
              + HEOS.rhomolar_reducing() * HEOS.T_reducing() * (HEOS.delta() * d2alpha_dxi_dDelta(HEOS, i, xN_flag) + dalpha_dxi(HEOS, i, xN_flag)));
}
CoolPropDbl MixtureDerivatives::d2psi_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / POW2(HEOS.tau())
           * (d_rhorTr_dxi(HEOS, i, xN_flag) * (HEOS.tau() * dalpha_dTau(HEOS) - alpha(HEOS))
              + HEOS.rhomolar_reducing() * HEOS.T_reducing() * (HEOS.tau() * d2alpha_dxi_dTau(HEOS, i, xN_flag) - dalpha_dxi(HEOS, i, xN_flag)));
}
CoolPropDbl MixtureDerivatives::d2psir_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / POW2(HEOS.tau())
           * (d_rhorTr_dxi(HEOS, i, xN_flag) * (HEOS.tau() * dalphar_dTau(HEOS) - alphar(HEOS))
              + HEOS.rhomolar_reducing() * HEOS.T_reducing() * (HEOS.tau() * d2alphar_dxi_dTau(HEOS, i, xN_flag) - dalphar_dxi(HEOS, i, xN_flag)));
}
CoolPropDbl MixtureDerivatives::d2psi_dxi_dxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / HEOS.tau()
           * (alpha(HEOS) * d2_rhorTr_dxidxj(HEOS, i, j, xN_flag) + dalpha_dxi(HEOS, i, xN_flag) * d_rhorTr_dxi(HEOS, j, xN_flag)
              + dalpha_dxi(HEOS, j, xN_flag) * d_rhorTr_dxi(HEOS, i, xN_flag)
              + HEOS.rhomolar_reducing() * HEOS.T_reducing() * d2alphadxidxj(HEOS, i, j, xN_flag));
}
CoolPropDbl MixtureDerivatives::d2psir_dxi_dxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
    return HEOS.delta() * HEOS.gas_constant() / HEOS.tau()
           * (alphar(HEOS) * d2_rhorTr_dxidxj(HEOS, i, j, xN_flag) + dalphar_dxi(HEOS, i, xN_flag) * d_rhorTr_dxi(HEOS, j, xN_flag)
              + dalphar_dxi(HEOS, j, xN_flag) * d_rhorTr_dxi(HEOS, i, xN_flag)
              + HEOS.rhomolar_reducing() * HEOS.T_reducing() * d2alphardxidxj(HEOS, i, j, xN_flag));
}

} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>

#    include "Backends/Cubics/CubicBackend.h"
#    include "Backends/Cubics/VTPRBackend.h"

using namespace CoolProp;

double mix_deriv_err_func(double numeric, double analytic) {
    if (std::abs(analytic) < DBL_EPSILON) {
        return std::abs(numeric - analytic);
    } else {
        return std::abs(numeric / analytic - 1);
    }
}

typedef CoolProp::MixtureDerivatives MD;

enum derivative
{
    NO_DERIV = 0,
    TAU,
    DELTA,
    XI,
    XJ,
    XK,
    T_CONSTP,
    P_CONSTT,
    T_CONSTRHO,
    RHO_CONSTT,
    CONST_TAU_DELTA,
    CONST_TRHO
};

typedef double (*zero_mole_fraction_pointer)(CoolProp::HelmholtzEOSMixtureBackend& HEOS, CoolProp::x_N_dependency_flag xN_flag);
typedef double (*one_mole_fraction_pointer)(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, CoolProp::x_N_dependency_flag xN_flag);
typedef double (*two_mole_fraction_pointer)(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j,
                                            CoolProp::x_N_dependency_flag xN_flag);
typedef double (*three_mole_fraction_pointer)(CoolProp::HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                              CoolProp::x_N_dependency_flag xN_flag);

template <class backend>
class DerivativeFixture
{
   public:
    shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS, HEOS_plus_tau, HEOS_minus_tau, HEOS_plus_delta, HEOS_minus_delta, HEOS_plus_T__constp,
      HEOS_minus_T__constp, HEOS_plus_p__constT, HEOS_minus_p__constT, HEOS_plus_T__constrho, HEOS_minus_T__constrho, HEOS_plus_rho__constT,
      HEOS_minus_rho__constT;
    std::vector<shared_ptr<CoolProp::HelmholtzEOSMixtureBackend>> HEOS_plus_z, HEOS_minus_z, HEOS_plus_z__constTrho, HEOS_minus_z__constTrho,
      HEOS_plus_n, HEOS_minus_n, HEOS_plus_2z, HEOS_minus_2z, HEOS_plus_2z__constTrho, HEOS_minus_2z__constTrho;
    CoolProp::x_N_dependency_flag xN;
    double dtau, ddelta, dz, dn, tol, dT, drho, dp;
    DerivativeFixture() : xN(XN_INDEPENDENT) {
        dtau = 1e-6;
        ddelta = 1e-6;
        dz = 1e-6;
        dn = 1e-6;
        dT = 1e-3;
        drho = 1e-3;
        dp = 1;
        tol = 5e-6;
        std::vector<std::string> names;
        names.push_back("n-Pentane");
        names.push_back("Ethane");
        names.push_back("n-Propane");
        names.push_back("n-Butane");
        std::vector<CoolPropDbl> mole_fractions;
        mole_fractions.push_back(0.1);
        mole_fractions.push_back(0.0);
        mole_fractions.push_back(0.3);
        mole_fractions.push_back(0.6);
        HEOS.reset(new backend(names));
        HEOS->set_mole_fractions(mole_fractions);
        HEOS->specify_phase(CoolProp::iphase_gas);
        HEOS->update_DmolarT_direct(300, 300);

        init_state(HEOS_plus_tau);
        init_state(HEOS_minus_tau);
        init_state(HEOS_plus_delta);
        init_state(HEOS_minus_delta);
        init_state(HEOS_plus_T__constp);
        init_state(HEOS_minus_T__constp);
        init_state(HEOS_plus_p__constT);
        init_state(HEOS_minus_p__constT);
        init_state(HEOS_plus_T__constrho);
        init_state(HEOS_minus_T__constrho);
        init_state(HEOS_plus_rho__constT);
        init_state(HEOS_minus_rho__constT);
        // Constant composition, varying state
        HEOS_plus_tau->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS->rhomolar_reducing(), HEOS->T_reducing() / (HEOS->tau() + dtau));
        HEOS_minus_tau->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS->rhomolar_reducing(), HEOS->T_reducing() / (HEOS->tau() - dtau));
        HEOS_plus_delta->update(CoolProp::DmolarT_INPUTS, (HEOS->delta() + ddelta) * HEOS->rhomolar_reducing(), HEOS->T_reducing() / HEOS->tau());
        HEOS_minus_delta->update(CoolProp::DmolarT_INPUTS, (HEOS->delta() - ddelta) * HEOS->rhomolar_reducing(), HEOS->T_reducing() / HEOS->tau());
        HEOS_plus_T__constp->update(CoolProp::PT_INPUTS, HEOS->p(), HEOS->T() + dT);
        HEOS_minus_T__constp->update(CoolProp::PT_INPUTS, HEOS->p(), HEOS->T() - dT);
        HEOS_plus_p__constT->update(CoolProp::PT_INPUTS, HEOS->p() + dp, HEOS->T());
        HEOS_minus_p__constT->update(CoolProp::PT_INPUTS, HEOS->p() - dp, HEOS->T());
        HEOS_plus_T__constrho->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar(), HEOS->T() + dT);
        HEOS_minus_T__constrho->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar(), HEOS->T() - dT);
        HEOS_plus_rho__constT->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar() + drho, HEOS->T());
        HEOS_minus_rho__constT->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar() - drho, HEOS->T());

        // Varying mole fractions
        HEOS_plus_z.resize(4);
        HEOS_minus_z.resize(4);
        HEOS_plus_z__constTrho.resize(4);
        HEOS_minus_z__constTrho.resize(4);
        HEOS_plus_2z.resize(4);
        HEOS_minus_2z.resize(4);
        HEOS_plus_2z__constTrho.resize(4);
        HEOS_minus_2z__constTrho.resize(4);
        for (int i = 0; i < HEOS_plus_z.size(); ++i) {
            init_state(HEOS_plus_z[i]);
            init_state(HEOS_plus_2z[i]);
            init_state(HEOS_plus_z__constTrho[i]);
            std::vector<double> zz = HEOS->get_mole_fractions(), zz2 = HEOS->get_mole_fractions();
            zz[i] += dz;
            zz2[i] += 2 * dz;
            if (xN == CoolProp::XN_DEPENDENT) {
                zz[zz.size() - 1] -= dz;
                zz2[zz2.size() - 1] -= 2 * dz;
            }
            HEOS_plus_z[i]->set_mole_fractions(zz);
            HEOS_plus_z[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_plus_z[i]->rhomolar_reducing(),
                                   HEOS_plus_z[i]->T_reducing() / HEOS->tau());
            HEOS_plus_2z[i]->set_mole_fractions(zz2);
            HEOS_plus_2z[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_plus_2z[i]->rhomolar_reducing(),
                                    HEOS_plus_2z[i]->T_reducing() / HEOS->tau());
            HEOS_plus_z__constTrho[i]->set_mole_fractions(zz);
            HEOS_plus_z__constTrho[i]->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar(), HEOS->T());
        }
        for (int i = 0; i < HEOS_minus_z.size(); ++i) {
            init_state(HEOS_minus_z[i]);
            init_state(HEOS_minus_2z[i]);
            init_state(HEOS_minus_z__constTrho[i]);
            std::vector<double> zz = HEOS->get_mole_fractions(), zz2 = HEOS->get_mole_fractions();
            zz[i] -= dz;
            zz2[i] -= 2 * dz;
            if (xN == CoolProp::XN_DEPENDENT) {
                zz[zz.size() - 1] += dz;
                zz2[zz2.size() - 1] += 2 * dz;
            }
            HEOS_minus_z[i]->set_mole_fractions(zz);
            HEOS_minus_z[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_minus_z[i]->rhomolar_reducing(),
                                    HEOS_minus_z[i]->T_reducing() / HEOS->tau());
            HEOS_minus_2z[i]->set_mole_fractions(zz2);
            HEOS_minus_2z[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_minus_2z[i]->rhomolar_reducing(),
                                     HEOS_minus_2z[i]->T_reducing() / HEOS->tau());
            HEOS_minus_z__constTrho[i]->set_mole_fractions(zz);
            HEOS_minus_z__constTrho[i]->update(CoolProp::DmolarT_INPUTS, HEOS->rhomolar(), HEOS->T());
        }

        // Varying mole numbers
        HEOS_plus_n.resize(4);
        HEOS_minus_n.resize(4);
        for (int i = 0; i < HEOS_plus_n.size(); ++i) {
            init_state(HEOS_plus_n[i]);
            std::vector<double> zz = HEOS->get_mole_fractions();
            zz[i] += dn;
            for (int j = 0; j < HEOS_minus_n.size(); ++j) {
                zz[i] /= 1 + dn;
            }
            HEOS_plus_n[i]->set_mole_fractions(zz);
            HEOS_plus_n[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_plus_n[i]->rhomolar_reducing(),
                                   HEOS_plus_n[i]->T_reducing() / HEOS->tau());
        }
        for (int i = 0; i < HEOS_plus_z.size(); ++i) {
            init_state(HEOS_minus_n[i]);
            std::vector<double> zz = HEOS->get_mole_fractions();
            zz[i] -= dn;
            for (int j = 0; j < HEOS_minus_n.size(); ++j) {
                zz[i] /= 1 - dn;
            }
            HEOS_minus_n[i]->set_mole_fractions(zz);
            HEOS_minus_n[i]->update(CoolProp::DmolarT_INPUTS, HEOS->delta() * HEOS_minus_n[i]->rhomolar_reducing(),
                                    HEOS_minus_n[i]->T_reducing() / HEOS->tau());
        }
    };
    void init_state(shared_ptr<CoolProp::HelmholtzEOSMixtureBackend>& other) {
        other.reset(HEOS->get_copy());
        other->set_mole_fractions(HEOS->get_mole_fractions());
        other->specify_phase(CoolProp::iphase_gas);  // Something homogeneous
    }
    void zero(const std::string& name, zero_mole_fraction_pointer f, zero_mole_fraction_pointer g = NULL, derivative wrt = NO_DERIV) {
        double analytic = f(*HEOS, xN);
        double numeric = 0;
        if (wrt == TAU) {
            numeric = (g(*HEOS_plus_tau, xN) - g(*HEOS_minus_tau, xN)) / (2 * dtau);
        } else if (wrt == DELTA) {
            numeric = (g(*HEOS_plus_delta, xN) - g(*HEOS_minus_delta, xN)) / (2 * ddelta);
        }
        CAPTURE(name);
        CAPTURE(analytic);
        CAPTURE(numeric);
        CAPTURE(xN);
        double error = mix_deriv_err_func(numeric, analytic);
        CAPTURE(error);
        CHECK(error < tol);
    }
    void one(const std::string& name, one_mole_fraction_pointer f, one_mole_fraction_pointer g = NULL, derivative wrt = NO_DERIV) {
        for (int i = 0; i < 4; ++i) {
            double analytic = f(*HEOS, i, xN);
            double numeric = 0;
            if (wrt == TAU) {
                numeric = (g(*HEOS_plus_tau, i, xN) - g(*HEOS_minus_tau, i, xN)) / (2 * dtau);
            } else if (wrt == DELTA) {
                numeric = (g(*HEOS_plus_delta, i, xN) - g(*HEOS_minus_delta, i, xN)) / (2 * ddelta);
            } else if (wrt == T_CONSTP) {
                numeric = (g(*HEOS_plus_T__constp, i, xN) - g(*HEOS_minus_T__constp, i, xN)) / (2 * dT);
            } else if (wrt == P_CONSTT) {
                numeric = (g(*HEOS_plus_p__constT, i, xN) - g(*HEOS_minus_p__constT, i, xN)) / (2 * dp);
            } else if (wrt == T_CONSTRHO) {
                double g1 = g(*HEOS_plus_T__constrho, i, xN), g2 = g(*HEOS_minus_T__constrho, i, xN);
                numeric = (g1 - g2) / (2 * dT);
            } else if (wrt == RHO_CONSTT) {
                numeric = (g(*HEOS_plus_rho__constT, i, xN) - g(*HEOS_minus_rho__constT, i, xN)) / (2 * drho);
            }
            CAPTURE(name);
            CAPTURE(i);
            CAPTURE(analytic);
            CAPTURE(numeric);
            CAPTURE(xN);
            double error = mix_deriv_err_func(numeric, analytic);
            CAPTURE(error);
            CHECK(error < tol);
        }
    };
    void one_comp(const std::string& name, one_mole_fraction_pointer f, zero_mole_fraction_pointer g, derivative wrt = CONST_TAU_DELTA) {
        for (int i = 0; i < 4; ++i) {
            double analytic;
            CHECK_NOTHROW(analytic = f(*HEOS, i, xN));
            double numeric = -10000;
            if (wrt == CONST_TAU_DELTA) {
                if (HEOS->get_mole_fractions()[i] > dz) {
                    CHECK_NOTHROW(numeric = (g(*HEOS_plus_z[i], xN) - g(*HEOS_minus_z[i], xN)) / (2 * dz));
                } else {
                    CHECK_NOTHROW(numeric = (-3 * g(*HEOS, xN) + 4 * g(*HEOS_plus_z[i], xN) - g(*HEOS_plus_2z[i], xN)) / (2 * dz));
                }
            } else if (wrt == CONST_TRHO) {
                CHECK_NOTHROW(numeric = (g(*HEOS_plus_z__constTrho[i], xN) - g(*HEOS_minus_z__constTrho[i], xN)) / (2 * dz));
            }

            CAPTURE(name);
            CAPTURE(i);
            CAPTURE(analytic);
            CAPTURE(numeric);
            CAPTURE(xN);
            double error = mix_deriv_err_func(numeric, analytic);
            CAPTURE(error);
            CHECK(error < tol);
        }
    }
    void two(const std::string& name, two_mole_fraction_pointer f, two_mole_fraction_pointer g = NULL, derivative wrt = NO_DERIV) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double analytic = f(*HEOS, i, j, xN);
                bool is_other = false;
                double numeric = 0;
                if (wrt == TAU) {
                    numeric = (g(*HEOS_plus_tau, i, j, xN) - g(*HEOS_minus_tau, i, j, xN)) / (2 * dtau);
                    is_other = true;
                } else if (wrt == DELTA) {
                    numeric = (g(*HEOS_plus_delta, i, j, xN) - g(*HEOS_minus_delta, i, j, xN)) / (2 * ddelta);
                    is_other = true;
                }
                CAPTURE(name);
                CAPTURE(i);
                CAPTURE(j);
                CAPTURE(analytic);
                CAPTURE(numeric);
                CAPTURE(xN);
                double error = mix_deriv_err_func(numeric, analytic);
                CAPTURE(error);
                CHECK(error < tol);
            }
        }
    }
    void two_comp(const std::string& name, two_mole_fraction_pointer f, one_mole_fraction_pointer g, derivative wrt = NO_DERIV) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double analytic = f(*HEOS, i, j, xN);
                double numeric = 500;
                if (HEOS->get_mole_fractions()[j] > 2 * dz) {
                    // Second order centered difference in composition
                    CHECK_NOTHROW(numeric = (g(*HEOS_plus_z[j], i, xN) - g(*HEOS_minus_z[j], i, xN)) / (2 * dz));
                } else {
                    // Forward difference in composition
                    CHECK_NOTHROW(numeric = (-3 * g(*HEOS, i, xN) + 4 * g(*HEOS_plus_z[j], i, xN) - g(*HEOS_plus_2z[j], i, xN)) / (2 * dz));
                }
                CAPTURE(name);
                CAPTURE(i);
                CAPTURE(j);
                CAPTURE(analytic);
                CAPTURE(numeric);
                CAPTURE(xN);
                double error = mix_deriv_err_func(numeric, analytic);
                CAPTURE(error);
                CHECK(error < tol);
            }
        }
    }
    void three(const std::string& name, three_mole_fraction_pointer f, three_mole_fraction_pointer g = NULL, derivative wrt = NO_DERIV) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    double analytic = f(*HEOS, i, j, k, xN);
                    double numeric = 0;
                    if (wrt == TAU) {
                        numeric = (g(*HEOS_plus_tau, i, j, k, xN) - g(*HEOS_minus_tau, i, j, k, xN)) / (2 * dtau);
                    } else if (wrt == DELTA) {
                        numeric = (g(*HEOS_plus_delta, i, j, k, xN) - g(*HEOS_minus_delta, i, j, k, xN)) / (2 * ddelta);
                    }
                    CAPTURE(name);
                    CAPTURE(i);
                    CAPTURE(j);
                    CAPTURE(analytic);
                    CAPTURE(numeric);
                    CAPTURE(xN);
                    double error = mix_deriv_err_func(numeric, analytic);
                    CAPTURE(error);
                    CHECK(error < tol);
                }
            }
        }
    }
    void three_comp(const std::string& name, three_mole_fraction_pointer f, two_mole_fraction_pointer g, derivative wrt = NO_DERIV) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    double analytic = f(*HEOS, i, j, k, xN);
                    double numeric;
                    if (HEOS->get_mole_fractions()[i] > 2 * dz) {
                        CHECK_NOTHROW(numeric = (g(*HEOS_plus_z[k], i, j, xN) - g(*HEOS_minus_z[k], i, j, xN)) / (2 * dz));
                    } else {
                        CHECK_NOTHROW(numeric =
                                        (-3 * g(*HEOS, i, j, xN) + 4 * g(*HEOS_plus_z[k], i, j, xN) - g(*HEOS_plus_2z[k], i, j, xN)) / (2 * dz));
                    }
                    CAPTURE(name);
                    CAPTURE(i);
                    CAPTURE(j);
                    CAPTURE(k);
                    CAPTURE(analytic);
                    CAPTURE(numeric);
                    CAPTURE(xN);
                    double error = mix_deriv_err_func(numeric, analytic);
                    CAPTURE(error);
                    CHECK(error < tol);
                }
            }
        }
    }
    Eigen::MatrixXd get_matrix(CoolProp::HelmholtzEOSMixtureBackend& HEOS, const std::string name) {
        Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(HEOS, xN);
        if (name == "Lstar") {
            return Lstar;
        } else if (name == "Mstar") {
            return MixtureDerivatives::Mstar(HEOS, xN, Lstar);
        } else {
            throw ValueError("Must be one of Lstar or Mstar");
        }
    }
    void matrix_derivative(const std::string name, const std::string& wrt) {
        CAPTURE(name);
        CAPTURE(wrt);
        double err = 10000;
        Eigen::MatrixXd analytic, numeric, Lstar, dLstar_dTau, dLstar_dDelta;
        if (name == "Mstar") {
            Lstar = MixtureDerivatives::Lstar(*HEOS, xN);
            dLstar_dDelta = MixtureDerivatives::dLstar_dX(*HEOS, xN, CoolProp::iDelta);
            dLstar_dTau = MixtureDerivatives::dLstar_dX(*HEOS, xN, CoolProp::iTau);
        }
        if (wrt == "tau") {
            if (name == "Lstar") {
                analytic = MixtureDerivatives::dLstar_dX(*HEOS, xN, CoolProp::iTau);
            } else if (name == "Mstar") {
                analytic = MixtureDerivatives::dMstar_dX(*HEOS, xN, CoolProp::iTau, Lstar, dLstar_dTau);
            } else {
                throw ValueError("argument not understood");
            }
            numeric = (get_matrix(*HEOS_plus_tau, name) - get_matrix(*HEOS_minus_tau, name)) / (2 * dtau);
        } else if (wrt == "delta") {
            if (name == "Lstar") {
                analytic = MixtureDerivatives::dLstar_dX(*HEOS, xN, CoolProp::iDelta);
            } else if (name == "Mstar") {
                analytic = MixtureDerivatives::dMstar_dX(*HEOS, xN, CoolProp::iDelta, Lstar, dLstar_dDelta);
            } else {
                throw ValueError("argument not understood");
            }
            numeric = (get_matrix(*HEOS_plus_delta, name) - get_matrix(*HEOS_minus_delta, name)) / (2 * ddelta);
        } else {
            throw ValueError("argument not understood");
        }
        CAPTURE(analytic);
        CAPTURE(numeric);
        Eigen::MatrixXd rel_error = ((analytic - numeric).cwiseQuotient(analytic));
        CAPTURE(rel_error);
        err = rel_error.squaredNorm();
        CAPTURE(err);
        CHECK(err < 1e-8);
    }

    void run_checks() {

        two_comp("d_PSI_rho_dxj", MD::d_PSI_rho_dxj, MD::PSI_rho);
        two_comp("d_PSI_T_dxj", MD::d_PSI_T_dxj, MD::PSI_T);

        one("d_ndalphardni_dDelta", MD::d_ndalphardni_dDelta, MD::ndalphar_dni__constT_V_nj, DELTA);
        one("d2_ndalphardni_dDelta2", MD::d2_ndalphardni_dDelta2, MD::d_ndalphardni_dDelta, DELTA);
        one("d3_ndalphardni_dDelta3", MD::d3_ndalphardni_dDelta3, MD::d2_ndalphardni_dDelta2, DELTA);
        one("d_ndalphardni_dTau", MD::d_ndalphardni_dTau, MD::ndalphar_dni__constT_V_nj, TAU);
        one("d2_ndalphardni_dTau2", MD::d2_ndalphardni_dTau2, MD::d_ndalphardni_dTau, TAU);
        one("d3_ndalphardni_dTau3", MD::d3_ndalphardni_dTau3, MD::d2_ndalphardni_dTau2, TAU);
        one("d2_ndalphardni_dDelta_dTau", MD::d2_ndalphardni_dDelta_dTau, MD::d_ndalphardni_dDelta, TAU);
        one("d3_ndalphardni_dDelta2_dTau", MD::d3_ndalphardni_dDelta2_dTau, MD::d2_ndalphardni_dDelta2, TAU);
        one("d3_ndalphardni_dDelta_dTau2", MD::d3_ndalphardni_dDelta_dTau2, MD::d2_ndalphardni_dDelta_dTau, TAU);

        zero("dalphar_dDelta", MD::dalphar_dDelta, MD::alphar, DELTA);
        zero("d2alphar_dDelta2", MD::d2alphar_dDelta2, MD::dalphar_dDelta, DELTA);
        zero("dalphar_dTau", MD::dalphar_dTau, MD::alphar, TAU);
        zero("d2alphar_dTau2", MD::d2alphar_dTau2, MD::dalphar_dTau, TAU);

        zero("dalpha0_dDelta", MD::dalpha0_dDelta, MD::alpha0, DELTA);
        zero("d2alpha0_dDelta2", MD::d2alpha0_dDelta2, MD::dalpha0_dDelta, DELTA);
        zero("dalpha0_dTau", MD::dalpha0_dTau, MD::alpha0, TAU);
        zero("d2alpha0_dTau2", MD::d2alpha0_dTau2, MD::dalpha0_dTau, TAU);

        one_comp("dalpha0_dxi", MD::dalpha0_dxi, MD::alpha0);
        one("d2alpha0_dxi_dDelta", MD::d2alpha0_dxi_dDelta, MD::dalpha0_dxi, DELTA);
        one("d2alpha0_dxi_dTau", MD::d2alpha0_dxi_dTau, MD::dalpha0_dxi, TAU);
        two_comp("d2alpha0dxidxj", MD::d2alpha0dxidxj, MD::dalpha0_dxi);

        one_comp("dalpha_dxi", MD::dalpha_dxi, MD::alpha);
        one("d2alpha_dxi_dDelta", MD::d2alpha_dxi_dDelta, MD::dalpha_dxi, DELTA);
        one("d2alpha_dxi_dTau", MD::d2alpha_dxi_dTau, MD::dalpha_dxi, TAU);
        two_comp("d2alphadxidxj", MD::d2alphadxidxj, MD::dalpha_dxi);

        zero("dpsi_dDelta", MD::dpsi_dDelta, MD::psi, DELTA);
        zero("dpsi_dTau", MD::dpsi_dTau, MD::psi, TAU);
        zero("d2psi_dDelta2", MD::d2psi_dDelta2, MD::dpsi_dDelta, DELTA);
        zero("d2psi_dDelta_dTau", MD::d2psi_dDelta_dTau, MD::dpsi_dDelta, TAU);
        zero("d2psi_dTau2", MD::d2psi_dTau2, MD::dpsi_dTau, TAU);
        one_comp("dpsi_dxi", MD::dpsi_dxi, MD::psi);
        one_comp("d2psi_dxi_dDelta", MD::d2psi_dxi_dDelta, MD::dpsi_dDelta);
        one_comp("d2psi_dxi_dTau", MD::d2psi_dxi_dTau, MD::dpsi_dTau);
        two_comp("d2psi_dxi_dxj", MD::d2psi_dxi_dxj, MD::dpsi_dxi);

        //two_comp("d_ndalphardni_dxj__constT_V_xi", MD::d_ndalphardni_dxj__constT_V_xi, MD::ndalphar_dni__constT_V_nj);

        one_comp("dalphar_dxi", MD::dalphar_dxi, MD::alphar);
        two_comp("d2alphardxidxj", MD::d2alphardxidxj, MD::dalphar_dxi);
        three_comp("d3alphardxidxjdxk", MD::d3alphardxidxjdxk, MD::d2alphardxidxj);
        one("d2alphar_dxi_dTau", MD::d2alphar_dxi_dTau, MD::dalphar_dxi, TAU);
        one("d2alphar_dxi_dDelta", MD::d2alphar_dxi_dDelta, MD::dalphar_dxi, DELTA);
        one("d3alphar_dxi_dDelta2", MD::d3alphar_dxi_dDelta2, MD::d2alphar_dxi_dDelta, DELTA);
        one("d3alphar_dxi_dTau2", MD::d3alphar_dxi_dTau2, MD::d2alphar_dxi_dTau, TAU);
        one("d4alphar_dxi_dTau3", MD::d4alphar_dxi_dTau3, MD::d3alphar_dxi_dTau2, TAU);
        one("d3alphar_dxi_dDelta_dTau", MD::d3alphar_dxi_dDelta_dTau, MD::d2alphar_dxi_dDelta, TAU);
        one("d4alphar_dxi_dDelta_dTau2", MD::d4alphar_dxi_dDelta_dTau2, MD::d3alphar_dxi_dDelta_dTau, TAU);
        one("d4alphar_dxi_dDelta2_dTau", MD::d4alphar_dxi_dDelta2_dTau, MD::d3alphar_dxi_dDelta2, TAU);
        two("d3alphar_dxi_dxj_dDelta", MD::d3alphar_dxi_dxj_dDelta, MD::d2alphardxidxj, DELTA);
        two("d4alphar_dxi_dxj_dDelta2", MD::d4alphar_dxi_dxj_dDelta2, MD::d3alphar_dxi_dxj_dDelta, DELTA);
        two("d4alphar_dxi_dxj_dDelta_dTau", MD::d4alphar_dxi_dxj_dDelta_dTau, MD::d3alphar_dxi_dxj_dDelta, TAU);
        two("d3alphar_dxi_dxj_dTau", MD::d3alphar_dxi_dxj_dTau, MD::d2alphardxidxj, TAU);
        two("d4alphar_dxi_dxj_dTau2", MD::d4alphar_dxi_dxj_dTau2, MD::d3alphar_dxi_dxj_dTau, TAU);
        one_comp("d_dalpharddelta_dxj__constT_V_xi", MD::d_dalpharddelta_dxj__constT_V_xi, MD::dalphar_dDelta, CONST_TRHO);

        two_comp("d_ndalphardni_dxj__constdelta_tau_xi", MD::d_ndalphardni_dxj__constdelta_tau_xi, MD::ndalphar_dni__constT_V_nj);
        two("d2_ndalphardni_dxj_dDelta__consttau_xi", MD::d2_ndalphardni_dxj_dDelta__consttau_xi, MD::d_ndalphardni_dxj__constdelta_tau_xi, DELTA);
        two("d3_ndalphardni_dxj_dDelta2__consttau_xi", MD::d3_ndalphardni_dxj_dDelta2__consttau_xi, MD::d2_ndalphardni_dxj_dDelta__consttau_xi,
            DELTA);
        two("d2_ndalphardni_dxj_dTau__constdelta_xi", MD::d2_ndalphardni_dxj_dTau__constdelta_xi, MD::d_ndalphardni_dxj__constdelta_tau_xi, TAU);
        two("d3_ndalphardni_dxj_dTau2__constdelta_xi", MD::d3_ndalphardni_dxj_dTau2__constdelta_xi, MD::d2_ndalphardni_dxj_dTau__constdelta_xi, TAU);
        two("d3_ndalphardni_dxj_dDelta_dTau__constxi", MD::d3_ndalphardni_dxj_dDelta_dTau__constxi, MD::d2_ndalphardni_dxj_dDelta__consttau_xi, TAU);
        three_comp("d2_ndalphardni_dxj_dxk__constdelta_tau_xi", MD::d2_ndalphardni_dxj_dxk__constdelta_tau_xi,
                   MD::d_ndalphardni_dxj__constdelta_tau_xi);
        three("d3_ndalphardni_dxj_dxk_dTau__constdelta_xi", MD::d3_ndalphardni_dxj_dxk_dTau__constdelta_xi,
              MD::d2_ndalphardni_dxj_dxk__constdelta_tau_xi, TAU);
        three("d3_ndalphardni_dxj_dxk_dDelta__consttau_xi", MD::d3_ndalphardni_dxj_dxk_dDelta__consttau_xi,
              MD::d2_ndalphardni_dxj_dxk__constdelta_tau_xi, DELTA);

        //        two("nd_ndalphardni_dnj__constT_V", MD::nd_ndalphardni_dnj__constT_V);
        two("d_nd_ndalphardni_dnj_dTau__constdelta_x", MD::d_nd_ndalphardni_dnj_dTau__constdelta_x, MD::nd_ndalphardni_dnj__constT_V, TAU);
        two("d2_nd_ndalphardni_dnj_dTau2__constdelta_x", MD::d2_nd_ndalphardni_dnj_dTau2__constdelta_x, MD::d_nd_ndalphardni_dnj_dTau__constdelta_x,
            TAU);
        two("d_nd_ndalphardni_dnj_dDelta__consttau_x", MD::d_nd_ndalphardni_dnj_dDelta__consttau_x, MD::nd_ndalphardni_dnj__constT_V, DELTA);
        two("d2_nd_ndalphardni_dnj_dDelta_dTau__constx", MD::d2_nd_ndalphardni_dnj_dDelta_dTau__constx, MD::d_nd_ndalphardni_dnj_dDelta__consttau_x,
            TAU);
        two("d2_nd_ndalphardni_dnj_dDelta2__consttau_x", MD::d2_nd_ndalphardni_dnj_dDelta2__consttau_x, MD::d_nd_ndalphardni_dnj_dDelta__consttau_x,
            DELTA);
        three_comp("d_nd_ndalphardni_dnj_dxk__consttau_delta", MD::d_nd_ndalphardni_dnj_dxk__consttau_delta, MD::nd_ndalphardni_dnj__constT_V);
        three("d2_nd_ndalphardni_dnj_dxk_dTau__constdelta", MD::d2_nd_ndalphardni_dnj_dxk_dTau__constdelta,
              MD::d_nd_ndalphardni_dnj_dxk__consttau_delta, TAU);
        three("d2_nd_ndalphardni_dnj_dxk_dDelta__consttau", MD::d2_nd_ndalphardni_dnj_dxk_dDelta__consttau,
              MD::d_nd_ndalphardni_dnj_dxk__consttau_delta, DELTA);

        // Xn-dep only        two_comp("dln_fugacity_dxj__constT_rho_xi", MD::dln_fugacity_dxj__constT_rho_xi, MD::ln_fugacity);
        three("d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau", MD::d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau,
              MD::d_ndln_fugacity_i_dnj_ddxk__consttau_delta, DELTA);
        three("d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta", MD::d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta,
              MD::d_ndln_fugacity_i_dnj_ddxk__consttau_delta, TAU);
        two("d2_ndln_fugacity_i_dnj_ddelta_dtau__constx", MD::d2_ndln_fugacity_i_dnj_ddelta_dtau__constx,
            MD::d_ndln_fugacity_i_dnj_ddelta__consttau_x, TAU);
        two("d_ndln_fugacity_i_dnj_ddelta__consttau_x", MD::d_ndln_fugacity_i_dnj_ddelta__consttau_x, MD::ndln_fugacity_i_dnj__constT_V_xi, DELTA);
        two("d_ndln_fugacity_i_dnj_dtau__constdelta_x", MD::d_ndln_fugacity_i_dnj_dtau__constdelta_x, MD::ndln_fugacity_i_dnj__constT_V_xi, TAU);
        three_comp("d_ndln_fugacity_i_dnj_ddxk__consttau_delta", MD::d_ndln_fugacity_i_dnj_ddxk__consttau_delta, MD::ndln_fugacity_i_dnj__constT_V_xi,
                   TAU);
        one("dln_fugacity_i_dT__constrho_n", MD::dln_fugacity_i_dT__constrho_n, MD::ln_fugacity, T_CONSTRHO);
        one("dln_fugacity_i_drho__constT_n", MD::dln_fugacity_i_drho__constT_n, MD::ln_fugacity, RHO_CONSTT);
        one("dln_fugacity_i_dT__constp_n", MD::dln_fugacity_i_dT__constp_n, MD::ln_fugacity, T_CONSTP);
        one("dln_fugacity_i_dp__constT_n", MD::dln_fugacity_i_dp__constT_n, MD::ln_fugacity, P_CONSTT);
        one("dln_fugacity_coefficient_dT__constp_n", MD::dln_fugacity_coefficient_dT__constp_n, MD::ln_fugacity_coefficient, T_CONSTP);
        one("dln_fugacity_coefficient_dp__constT_n", MD::dln_fugacity_coefficient_dp__constT_n, MD::ln_fugacity_coefficient, P_CONSTT);

        three_comp("d2_PSI_T_dxj_dxk", MD::d2_PSI_T_dxj_dxk, MD::d_PSI_T_dxj);
        three_comp("d2_PSI_rho_dxj_dxk", MD::d2_PSI_rho_dxj_dxk, MD::d_PSI_rho_dxj);

        three("d_n2Aijk_dTau", MD::d_n2Aijk_dTau, MD::n2Aijk, TAU);
        three("d_n2Aijk_dDelta", MD::d_n2Aijk_dDelta, MD::n2Aijk, DELTA);
        two("d_nAij_dTau", MD::d_nAij_dTau, MD::nAij, TAU);
        two("d_nAij_dDelta", MD::d_nAij_dDelta, MD::nAij, DELTA);

        two_comp("d_nddeltadni_dxj__constdelta_tau", MD::d_nddeltadni_dxj__constdelta_tau, MD::nddeltadni__constT_V_nj);
        two_comp("d_ndtaudni_dxj__constdelta_tau", MD::d_ndtaudni_dxj__constdelta_tau, MD::ndtaudni__constT_V_nj);
        two("d2_ndtaudni_dxj_dTau__constdelta", MD::d2_ndtaudni_dxj_dTau__constdelta, MD::d_ndtaudni_dxj__constdelta_tau, TAU);

        one_comp("dTrdxi__constxj", MD::dTrdxi__constxj, MD::Tr);
        one_comp("d2Tr_dxidbetaT", MD::d2Tr_dxidbetaT, MD::dTr_dbetaT);
        one_comp("d2Tr_dxidgammaT", MD::d2Tr_dxidgammaT, MD::dTr_dgammaT);
        // (??)two_comp("d2Trdxi2__constxj", MD::d2Trdxi2__constxj, MD::dTrdxi__constxj);
        two_comp("d2Trdxidxj", MD::d2Trdxidxj, MD::dTrdxi__constxj);
        three_comp("d3Trdxidxjdxk", MD::d3Trdxidxjdxk, MD::d2Trdxidxj);
        one_comp("dtaudxj__constT_V_xi", MD::dtau_dxj__constT_V_xi, MD::tau, CONST_TRHO);
        two_comp("d_ndTrdni_dxj__constxi", MD::d_ndTrdni_dxj__constxi, MD::ndTrdni__constnj);
        three_comp("d2_ndTrdni_dxj_dxk__constxi", MD::d2_ndTrdni_dxj_dxk__constxi, MD::d_ndTrdni_dxj__constxi);

        one_comp("drhormolardxi__constxj", MD::drhormolardxi__constxj, MD::rhormolar);
        one_comp("d2rhormolar_dxidbetaV", MD::d2rhormolar_dxidbetaV, MD::drhormolar_dbetaV);
        one_comp("d2rhormolar_dxidgammaV", MD::d2rhormolar_dxidgammaV, MD::drhormolar_dgammaV);
        // (??) two_comp("d2rhormolardxi2__constxj", MD::d2rhormolardxi2__constxj, MD::drhormolardxi__constxj);
        two_comp("d2rhormolardxidxj", MD::d2rhormolardxidxj, MD::drhormolardxi__constxj);
        three_comp("d3rhormolardxidxjdxk", MD::d3rhormolardxidxjdxk, MD::d2rhormolardxidxj);
        one_comp("ddelta_dxj__constT_V_xi", MD::ddelta_dxj__constT_V_xi, MD::delta, CONST_TRHO);
        two_comp("d_ndrhorbardni_dxj__constxi", MD::d_ndrhorbardni_dxj__constxi, MD::ndrhorbardni__constnj);
        three_comp("d2_ndrhorbardni_dxj_dxk__constxi", MD::d2_ndrhorbardni_dxj_dxk__constxi, MD::d_ndrhorbardni_dxj__constxi);

        one_comp("dpdxj__constT_V_xi", MD::dpdxj__constT_V_xi, MD::p, CONST_TRHO);

        matrix_derivative("Lstar", "tau");
        matrix_derivative("Lstar", "delta");
        matrix_derivative("Mstar", "tau");
        matrix_derivative("Mstar", "delta");
    }
};

TEST_CASE_METHOD(DerivativeFixture<HelmholtzEOSMixtureBackend>, "Check derivatives for HEOS", "[mixture_derivs2]") {
    run_checks();
};
TEST_CASE_METHOD(DerivativeFixture<PengRobinsonBackend>, "Check derivatives for Peng-Robinson", "[mixture_derivs2]") {
    tol = 1e-4;  // Relax the tolerance a bit
    run_checks();
};
TEST_CASE_METHOD(DerivativeFixture<SRKBackend>, "Check derivatives for SRK", "[mixture_derivs2]") {
    tol = 1e-4;  // Relax the tolerance a bit
    run_checks();
};
    // Make sure you set the VTPR UNIFAC path with something like set_config_string(VTPR_UNIFAC_PATH, "/Users/ian/Code/CUBAC/dev/unifaq/");
    //TEST_CASE_METHOD(DerivativeFixture<VTPRBackend>, "Check derivatives for VTPR", "[mixture_derivs2]")
    //{
    //    tol = 1e-4; // Relax the tolerance a bit
    //    run_checks();
    //};

#endif
