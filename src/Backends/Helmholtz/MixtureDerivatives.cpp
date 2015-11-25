#include "MixtureDerivatives.h"

namespace CoolProp{

CoolPropDbl MixtureDerivatives::dln_fugacity_i_dT__constp_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return dln_fugacity_coefficient_dT__constp_n(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dp__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return dln_fugacity_coefficient_dp__constT_n(HEOS, i, xN_flag) + 1/HEOS.p();
}
CoolPropDbl MixtureDerivatives::fugacity_i(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
    return HEOS.mole_fractions[i]*HEOS.rhomolar()*HEOS.gas_constant()*HEOS.T()*exp( dnalphar_dni__constT_V_nj(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::ln_fugacity_coefficient(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag)-log(1+HEOS._delta.pt()*HEOS.dalphar_dDelta());
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dT__constrho_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return 1/HEOS.T()*(1-HEOS.tau()*HEOS.dalphar_dTau()-HEOS.tau()*d_ndalphardni_dTau(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_drho__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return 1/HEOS.rhomolar()*(1+HEOS.delta()*HEOS.dalphar_dDelta() + HEOS.delta()*d_ndalphardni_dDelta(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dnalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // GERG Equation 7.42
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2nalphar_dni_dT(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -HEOS._tau.pt()/HEOS._T*(HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double T = HEOS._reducing.T/HEOS._tau.pt();
    CoolPropDbl R_u = HEOS.gas_constant();
    return d2nalphar_dni_dT(HEOS, i, xN_flag) + 1/T-partial_molar_volume(HEOS, i, xN_flag)/(R_u*T)*dpdT__constV_n(HEOS);
}
CoolPropDbl MixtureDerivatives::partial_molar_volume(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -ndpdni__constT_V_nj(HEOS, i, xN_flag)/ndpdV__constT_n(HEOS);
}

CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // GERG equation 7.30
    CoolPropDbl R_u = HEOS.gas_constant();
    double partial_molar_volumeval = partial_molar_volume(HEOS, i, xN_flag); // [m^3/mol]
    double term1 = partial_molar_volumeval/(R_u*HEOS._T); // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
    double term2 = 1.0/HEOS.p();
    return term1 - term2;
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_dtau__constdelta_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -1/HEOS.tau() + HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_i_ddelta__consttau_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return 1 + HEOS.delta()*HEOS.dalphar_dDelta() + HEOS.delta()*d_ndalphardni_dDelta(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    // This is a term to which some more might be added depending on i and j
    CoolPropDbl val = dln_fugacity_coefficient_dxj__constT_p_xi(HEOS, i, j, xN_flag);
    const std::vector<CoolPropDbl> &x = HEOS.get_mole_fractions();
    std::size_t N = x.size();
    if (i == N-1){
        val += -1/x[N-1];
    }
    else if (i == j){
        val += 1/x[j];
    }
    return val;
}
CoolPropDbl MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){throw ValueError("dln_fugacity_dxj__constT_rho_xi only valid for xN_DEPENDENT for now");}
    CoolPropDbl rhor = HEOS.Reducing->rhormolar(HEOS.get_mole_fractions());
    CoolPropDbl Tr = HEOS.Reducing->Tr(HEOS.get_mole_fractions());
    CoolPropDbl dTrdxj = HEOS.Reducing->dTrdxi__constxj(HEOS.get_mole_fractions(),j,xN_flag);
    CoolPropDbl drhordxj = HEOS.Reducing->drhormolardxi__constxj(HEOS.get_mole_fractions(),j,xN_flag);
    
    // These lines are all the same
    CoolPropDbl line1 = dln_fugacity_i_dtau__constdelta_x(HEOS, i, xN_flag)*1/HEOS.T()*dTrdxj;
    CoolPropDbl line2 = -dln_fugacity_i_ddelta__consttau_x(HEOS, i, xN_flag)*1/rhor*drhordxj;
    CoolPropDbl line4 = HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag) + d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag);
    
    const std::vector<CoolPropDbl> &x = HEOS.get_mole_fractions();
    std::size_t N = x.size();
    
    CoolPropDbl line3 = 1/rhor*HEOS.Reducing->drhormolardxi__constxj(x, j, xN_flag) + 1/Tr*HEOS.Reducing->dTrdxi__constxj(x, j, xN_flag);;
    
    // This is a term to which some more might be added depending on i and j
    if (i == N-1){
        line3 += -1/x[N-1];
    }
    else if (i == j){
        line3 += 1/x[j];
    }
    else{}

    return line1 + line2 + line3 + line4;
}

CoolPropDbl MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return Kronecker_delta(i, j)/HEOS.mole_fractions[i]+nd2nalphardnidnj__constT_V(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return d_ndalphardni_dTau(HEOS, j, xN_flag) + d_nd_ndalphardni_dnj_dTau__constdelta_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return d2_ndalphardni_dTau2(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dTau2__constdelta_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_ddelta2__consttau_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return d2_ndalphardni_dDelta2(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dDelta2__consttau_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return d2_ndalphardni_dDelta_dTau(HEOS, j, xN_flag) + d2_nd_ndalphardni_dnj_dDelta_dTau__constx(HEOS, i, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    return d_ndalphardni_dDelta(HEOS, j, xN_flag) + d_nd_ndalphardni_dnj_dDelta__consttau_x(HEOS, i, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){
    CoolPropDbl s = -Kronecker_delta(i, j)*Kronecker_delta(i, k)/pow(HEOS.mole_fractions[i], 2);
    return s + d_ndalphardni_dxj__constdelta_tau_xi(HEOS, j, k, xN_flag) + d_nd_ndalphardni_dnj_dxk__consttau_delta(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){
    return d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, j, k, xN_flag) + d2_nd_ndalphardni_dnj_dxk_dTau__constdelta(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){
    return d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, j, k, xN_flag) + d2_nd_ndalphardni_dnj_dxk_dDelta__consttau(HEOS, i, j, k, xN_flag);
}
CoolPropDbl MixtureDerivatives::nd_ndln_fugacity_i_dnj_dnk__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag){
    double sum = d_ndln_fugacity_i_dnj_dtau__constdelta_x(HEOS, i, j, xN_flag)*ndtaudni__constT_V_nj(HEOS, k, xN_flag)
        + d_ndln_fugacity_i_dnj_ddelta__consttau_x(HEOS, i, j, xN_flag)*nddeltadni__constT_V_nj(HEOS, k, xN_flag)
        + d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HEOS, i, j, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ mmax--; }
    for (unsigned int m = 0; m < mmax; ++m)
    {
        sum -= HEOS.mole_fractions[m]*d_ndln_fugacity_i_dnj_ddxk__consttau_delta(HEOS, i, j, m, xN_flag);
    }
    return sum;
}

CoolPropDbl MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.115
    CoolPropDbl R_u = HEOS.gas_constant();
    // partial molar volume is -dpdn/dpdV, so need to flip the sign here
    return d2nalphar_dxj_dni__constT_V(HEOS, j, i, xN_flag) - partial_molar_volume(HEOS, i, xN_flag)/(R_u*HEOS._T)*dpdxj__constT_V_xi(HEOS, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::dpdxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.130
    CoolPropDbl R_u = HEOS.gas_constant();
    return HEOS._rhomolar*R_u*HEOS._T*(ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)*HEOS.dalphar_dDelta()+HEOS._delta.pt()*d_dalpharddelta_dxj__constT_V_xi(HEOS, j, xN_flag));
}

CoolPropDbl MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert Equation 3.134 (Catch test provided)
    return HEOS.d2alphar_dDelta2()*ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)
         + HEOS.d2alphar_dDelta_dTau()*dtau_dxj__constT_V_xi(HEOS, j, xN_flag)
         + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    //Gernert 3.119 (Catch test provided)
    return HEOS.dalphar_dDelta()*ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)+HEOS.dalphar_dTau()*dtau_dxj__constT_V_xi(HEOS, j, xN_flag)+HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.118
    return d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i,j, xN_flag)
          + ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)*d_ndalphardni_dDelta(HEOS, i, xN_flag)
          + dtau_dxj__constT_V_xi(HEOS, j, xN_flag)*d_ndalphardni_dTau(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::ddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.121 (Catch test provided)
    return -HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions,j, xN_flag);
}
CoolPropDbl MixtureDerivatives::dtau_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.122 (Catch test provided)
    return 1/HEOS._T*HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions,j, xN_flag);
}

CoolPropDbl MixtureDerivatives::dpdT__constV_n(HelmholtzEOSMixtureBackend &HEOS)
{
    CoolPropDbl R_u = HEOS.gas_constant();
    return HEOS._rhomolar*R_u*(1+HEOS._delta.pt()*HEOS.dalphar_dDelta()-HEOS._delta.pt()*HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau());
}
CoolPropDbl MixtureDerivatives::dpdrho__constT_n(HelmholtzEOSMixtureBackend &HEOS)
{
    CoolPropDbl R_u = HEOS.gas_constant();
    return R_u*HEOS._T*(1+2*HEOS._delta.pt()*HEOS.dalphar_dDelta()+pow(HEOS._delta.pt(),2)*HEOS.d2alphar_dDelta2());
}
CoolPropDbl MixtureDerivatives::ndpdV__constT_n(HelmholtzEOSMixtureBackend &HEOS)
{
    CoolPropDbl R_u = HEOS.gas_constant();
    return -pow(HEOS._rhomolar,2)*R_u*HEOS._T*(1+2*HEOS._delta.pt()*HEOS.dalphar_dDelta()+pow(HEOS._delta.pt(),2)*HEOS.d2alphar_dDelta2());
}
CoolPropDbl MixtureDerivatives::ndpdni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // Eqn 7.64 and 7.63
    CoolPropDbl R_u = HEOS.gas_constant();
    double ndrhorbar_dni__constnj = HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag);
    double ndTr_dni__constnj = HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions,i, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; ++k)
    {
        summer += HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    double nd2alphar_dni_dDelta = HEOS._delta.pt()*HEOS.d2alphar_dDelta2()*(1-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()/HEOS._reducing.T*ndTr_dni__constnj+HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, i, xN_flag)-summer;
    return HEOS._rhomolar*R_u*HEOS._T*(1+HEOS._delta.pt()*HEOS.dalphar_dDelta()*(2-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._delta.pt()*nd2alphar_dni_dDelta);
}

CoolPropDbl MixtureDerivatives::ndalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag));
    double term2 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions,i, xN_flag);

    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        s += HEOS.mole_fractions[k]*HEOS.residual_helmholtz->dalphar_dxi(HEOS, k, xN_flag);
    }
    double term3 = HEOS.residual_helmholtz->dalphar_dxi(HEOS, i, xN_flag);
    return term1 + term2 + term3 - s;
}
CoolPropDbl MixtureDerivatives::ndln_fugacity_coefficient_dnj__constT_p(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    CoolPropDbl R_u = HEOS.gas_constant();
    return nd2nalphardnidnj__constT_V(HEOS, j, i, xN_flag) + 1 - partial_molar_volume(HEOS, j, xN_flag)/(R_u*HEOS._T)*ndpdni__constT_V_nj(HEOS, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::nddeltadni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS._delta.pt()-HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_nddeltadni_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return 1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_nddeltadni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double rhor = HEOS._reducing.rhomolar;
    return -HEOS.delta()/rhor*(HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)
        -1/rhor*HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions, j, xN_flag)*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::d2_nddeltadni_dxj_dDelta__consttau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    return d_nddeltadni_dxj__constdelta_tau(HEOS, i, j, xN_flag)/HEOS.delta();
}

CoolPropDbl MixtureDerivatives::ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS._tau.pt()/HEOS._reducing.T*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}

CoolPropDbl MixtureDerivatives::d_ndtaudni_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return 1/HEOS._reducing.T*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndtaudni_dxj__constdelta_tau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double Tr = HEOS._reducing.T;
    return HEOS.tau()/Tr*(HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag) 
                      -1/Tr*HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));
}
CoolPropDbl MixtureDerivatives::d2_ndtaudni_dxj_dTau__constdelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    return d_ndtaudni_dxj__constdelta_tau(HEOS, i, j, xN_flag)/HEOS.tau();
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = HEOS._delta.pt()*HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag)*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double line3 = HEOS._tau.pt()*HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag)*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
    double line2 = -HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1/HEOS._reducing.rhomolar)*(HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)-1/HEOS._reducing.rhomolar*HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions,j, xN_flag)*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag));
    double line4 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*(HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.mole_fractions,i,j, xN_flag)-1/HEOS._reducing.T*HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));
    
    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        s += HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d2alphardxidxj(HEOS, j, k, xN_flag);
    }
    double line5 = HEOS.residual_helmholtz->d2alphardxidxj(HEOS, i, j, xN_flag)-HEOS.residual_helmholtz->dalphar_dxi(HEOS, j, xN_flag)-s;
    return line1 + line2 + line3 + line4 + line5;
}

CoolPropDbl MixtureDerivatives::nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    return ndalphar_dni__constT_V_nj(HEOS, j, xN_flag) // First term from 7.46
    + nd_ndalphardni_dnj__constT_V(HEOS, i, j, xN_flag); // 7.47
}

CoolPropDbl MixtureDerivatives::nd_ndalphardni_dnj__constT_V(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double line3 = d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        line3 -= HEOS.mole_fractions[k]*d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dTau__constdelta_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dTau2(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag)
                   + d_ndalphardni_dTau(HEOS, i, xN_flag)*d_ndtaudni_dTau(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag);
    }
    double line3 = d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, j, xN_flag)-summer;
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dTau2__constdelta_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d3_ndalphardni_dDelta_dTau2(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = 2*d2_ndalphardni_dTau2(HEOS, i, xN_flag)*d_ndtaudni_dTau(HEOS, j, xN_flag);
    double line3 = d3_ndalphardni_dTau3(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, k, xN_flag);
    }
    double line4 = d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, j, xN_flag)-summer;
    return line1 + line2 + line3 + line4;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dDelta2__consttau_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d3_ndalphardni_dDelta3(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = 2*d2_ndalphardni_dDelta2(HEOS, i, xN_flag)*d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line3 = d3_ndalphardni_dDelta2_dTau(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, k, xN_flag);
    }
    double line4 = d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, j, xN_flag)-summer;
    return line1 + line2 + line3 + line4;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dDelta_dTau__constx(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d3_ndalphardni_dDelta2_dTau(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line3 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*d_ndtaudni_dTau(HEOS, j, xN_flag); 
    double line4 = d3_ndalphardni_dDelta_dTau2(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag);
    }
    double line5 = d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, j, xN_flag)-summer;
    return line1 + line2 + line3 + line4 + line5;
}
CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dDelta__consttau_x(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = d2_ndalphardni_dDelta2(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag)
                + d_ndalphardni_dDelta(HEOS, i, xN_flag)*d_nddeltadni_dDelta(HEOS, j, xN_flag);
    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag);
    }
    double line3 = d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, j, xN_flag)-summer;
    return line1 + line2 + line3;
}

CoolPropDbl MixtureDerivatives::d_nd_ndalphardni_dnj_dxk__consttau_delta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag)*d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                + d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    
    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag)*d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
                + d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HEOS, i, j, k, xN_flag) - d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k, xN_flag);

    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ mmax--; }
    for (unsigned int m = 0; m < mmax; m++)
    {
        line3 -= HEOS.mole_fractions[m]*d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dxk_dTau__constdelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double line1 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
        + d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);

    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag)*d2_ndtaudni_dxj_dTau__constdelta(HEOS, j, k, xN_flag)
        + d2_ndalphardni_dTau2(HEOS, i, xN_flag)*d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
        + d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag)*d_ndtaudni_dTau(HEOS, j, xN_flag)
        + d3_ndalphardni_dxj_dTau2__constdelta_xi(HEOS, i, k, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HEOS, i, j, k, xN_flag) - d2_ndalphardni_dxj_dTau__constdelta_xi(HEOS, i, k, xN_flag);

    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ mmax--; }
    for (unsigned int m = 0; m < mmax; m++)
    {
        line3 -= HEOS.mole_fractions[m]*d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d2_nd_ndalphardni_dnj_dxk_dDelta__consttau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag)*d2_nddeltadni_dxj_dDelta__consttau(HEOS, j, k, xN_flag)
        + d2_ndalphardni_dDelta2(HEOS, i, xN_flag)*d_nddeltadni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
        + d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag)*d_nddeltadni_dDelta(HEOS, j, xN_flag)
        + d3_ndalphardni_dxj_dDelta2__consttau_xi(HEOS, i, k, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag); 
    
    double line2 = d2_ndalphardni_dDelta_dTau(HEOS, i, xN_flag)*d_ndtaudni_dxj__constdelta_tau(HEOS, j, k, xN_flag)
        + d3_ndalphardni_dxj_dDelta_dTau__constxi(HEOS, i, k, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);

    double line3 = d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HEOS, i, j, k, xN_flag) - d2_ndalphardni_dxj_dDelta__consttau_xi(HEOS, i, k, xN_flag);
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ mmax--; }
    for (unsigned int m = 0; m < mmax; m++)
    {
        line3 -= HEOS.mole_fractions[m]*d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HEOS, i, m, k, xN_flag);
    }
    return line1 + line2 + line3;
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = (HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta())*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dDelta2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = (2*HEOS.d2alphar_dDelta2() + HEOS.delta()*HEOS.d3alphar_dDelta3())*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = HEOS.tau()*HEOS.d3alphar_dDelta2_dTau()*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta3(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = (3*HEOS.d3alphar_dDelta3() + HEOS.delta()*HEOS.d4alphar_dDelta4())*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = HEOS.tau()*HEOS.d4alphar_dDelta3_dTau()*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dDelta_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = (HEOS.d2alphar_dDelta_dTau() + HEOS.delta()*HEOS.d3alphar_dDelta2_dTau())*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (HEOS.tau()*HEOS.d3alphar_dDelta_dTau2() + HEOS.d2alphar_dDelta_dTau())*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta2_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = (2*HEOS.d3alphar_dDelta2_dTau() + HEOS.delta()*HEOS.d4alphar_dDelta3_dTau())*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (HEOS.tau()*HEOS.d4alphar_dDelta2_dTau2() + HEOS.d3alphar_dDelta2_dTau())*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dDelta_dTau2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = (HEOS.d3alphar_dDelta_dTau2() + HEOS.delta()*HEOS.d4alphar_dDelta2_dTau2())*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (HEOS.tau()*HEOS.d4alphar_dDelta_dTau3() + 2*HEOS.d3alphar_dDelta_dTau2())*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dTau2(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS.delta()*HEOS.d3alphar_dDelta_dTau2()*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (2*HEOS.d2alphar_dTau2() + HEOS.tau()*HEOS.d3alphar_dTau3())*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dTau3(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS.delta()*HEOS.d4alphar_dDelta_dTau3()*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term2 = (3*HEOS.d3alphar_dTau3() + HEOS.tau()*HEOS.d4alphar_dTau4())*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dDelta__consttau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double term1 = (HEOS.dalphar_dDelta() + HEOS.delta()*HEOS.d2alphar_dDelta2())*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag) + HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag))*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.tau()*HEOS.d2alphar_dDelta_dTau()*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term5 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, k, j, xN_flag) + Kronecker_delta(k, j)*HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dTau__constdelta_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS.delta()*HEOS.d2alphar_dDelta_dTau()*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = (HEOS.tau()*HEOS.d2alphar_dTau2()+HEOS.dalphar_dTau())*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = (HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag)+HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag))*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term5 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, k, j, xN_flag) + Kronecker_delta(k, j)*HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dTau2__constdelta_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS.delta()*HEOS.d3alphar_dDelta_dTau2()*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = HEOS.delta()*HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, j, xN_flag)*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = (HEOS.tau()*HEOS.d3alphar_dTau3()+2*HEOS.d2alphar_dTau2())*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = (HEOS.tau()*HEOS.residual_helmholtz->d4alphar_dxi_dTau3(HEOS, j, xN_flag)+2*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag))*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term5 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, k, j, xN_flag) + Kronecker_delta(k, j)*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dDelta2__consttau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    double term1 = (2*HEOS.d2alphar_dDelta2() + HEOS.delta()*HEOS.d3alphar_dDelta3())*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (2*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag)+HEOS.delta()*HEOS.residual_helmholtz->d4alphar_dxi_dDelta3(HEOS, j, xN_flag))*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.tau()*HEOS.d3alphar_dDelta2_dTau()*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = HEOS.tau()*HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, j, xN_flag)*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term5 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, k, j, xN_flag) + Kronecker_delta(k, j)*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}
CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dDelta_dTau__constxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag){
    double term1 = (HEOS.d2alphar_dDelta_dTau() + HEOS.delta()*HEOS.d3alphar_dDelta2_dTau())*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = (HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)+HEOS.delta()*HEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(HEOS, j, xN_flag))*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = (HEOS.tau()*HEOS.d3alphar_dDelta_dTau2()+HEOS.d2alphar_dDelta_dTau())*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term4 = (HEOS.tau()*HEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(HEOS, j, xN_flag)+HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag))*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term5 = HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, i, j, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term5 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, k, j, xN_flag) + Kronecker_delta(k, j)*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3 + term4 + term5;
}

CoolPropDbl MixtureDerivatives::d_ndalphardni_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = (HEOS._tau.pt()*HEOS.d2alphar_dTau2()+HEOS.dalphar_dTau())*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d2_ndalphardni_dxj_dxk__constdelta_tau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS.delta()*(HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag)*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag)+HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag)*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag));
    double term2 = HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag)*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term3 = HEOS.delta()*HEOS.dalphar_dDelta()*HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    
    double term4 = HEOS.tau()*(HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag)*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag)+HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag)*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag));
    double term5 = HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag)*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term6 = HEOS.tau()*HEOS.dalphar_dTau()*HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);

    /// All derivatives with dalphar/(dxi,dxj,dxk) are zero
    double term7 = -2*HEOS.residual_helmholtz->d2alphardxidxj(HEOS, j, k, xN_flag);

    return term1 + term2 + term3 + term4 + term5 + term6 + term7;
}

CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double term1a = HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term1b = HEOS.delta()*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, j, k, xN_flag)*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term1c = HEOS.delta()*HEOS.d2alphar_dDelta_dTau()*HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term1d = HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag)*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term1 = term1a + term1b + term1c + term1d;

    double term2a = (HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, j, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, j, xN_flag))*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term2b = (HEOS.tau()*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(HEOS, j, k, xN_flag) + HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag))*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term2c = (HEOS.tau()*HEOS.d2alphar_dTau2() + HEOS.dalphar_dTau())*HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term2d = (HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dTau2(HEOS, k, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dTau(HEOS, k, xN_flag))*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = term2a + term2b + term2c + term2d;

    /// All derivatives of dalphar/(dxi,dxj,dxk) are zero
    double term3 = -2*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(HEOS, j, k, xN_flag);
    return term1 + term2 + term3;
}

CoolPropDbl MixtureDerivatives::d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, std::size_t k, x_N_dependency_flag xN_flag)
{
    double term1a = (HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, j, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, j, xN_flag))*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term1b = (HEOS.delta()*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(HEOS, j, k, xN_flag) + HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag))*HEOS.Reducing->PSI_rho(HEOS.mole_fractions, i, xN_flag);
    double term1c = (HEOS.delta()*HEOS.d2alphar_dDelta2() + HEOS.dalphar_dDelta())*HEOS.Reducing->d2_PSI_rho_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term1d = (HEOS.delta()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta2(HEOS, k, xN_flag) + HEOS.residual_helmholtz->d2alphar_dxi_dDelta(HEOS, k, xN_flag))*HEOS.Reducing->d_PSI_rho_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term1 = term1a + term1b + term1c + term1d;

    double term2a = HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, j, xN_flag)*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, k, xN_flag);
    double term2b = HEOS.tau()*HEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(HEOS, j, k, xN_flag)*HEOS.Reducing->PSI_T(HEOS.mole_fractions, i, xN_flag);
    double term2c = HEOS.tau()*HEOS.d2alphar_dDelta_dTau()*HEOS.Reducing->d2_PSI_T_dxj_dxk(HEOS.mole_fractions, i, j, k, xN_flag);
    double term2d = HEOS.tau()*HEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(HEOS, k, xN_flag)*HEOS.Reducing->d_PSI_T_dxj(HEOS.mole_fractions, i, j, xN_flag);
    double term2 = term2a + term2b + term2c + term2d;

    /// All derivatives of dalphar)/(dxi,dxj,dxk) are zero
    double term3 = -2*HEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(HEOS, j, k, xN_flag);
    return term1 + term2 + term3;
}


} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#include "catch.hpp"

using namespace CoolProp;

double mix_deriv_err_func(double numeric, double analytic)
{
    if (std::abs(analytic) < DBL_EPSILON){
        return std::abs(numeric - analytic);
    }
    else{
        return std::abs(numeric/analytic-1);
    }
}

static bool fluids_set = false;
static const std::size_t Ncomp_max = 6;

// These are composition invariant
// ** Levels **
// 0: number of components in the mixture
// 1: component index
static std::vector<std::vector<shared_ptr<HelmholtzEOSMixtureBackend> > > HEOS, 
                                                                          HEOS_plusT_constrho, HEOS_minusT_constrho, 
                                                                          HEOS_plusrho_constT, HEOS_minusrho_constT,
                                                                          HEOS_plusz_xNindep, HEOS_minusz_xNindep,
                                                                          HEOS_plusz_xNdep, HEOS_minusz_xNdep,
                                                                          HEOS_plusz_consttaudelta_xNindep, HEOS_minusz_consttaudelta_xNindep,
                                                                          HEOS_plusz_consttaudelta_xNdep, HEOS_minusz_consttaudelta_xNdep;

static const double T1 = 319.325, rho1 = 13246.6, dT = 1e-3, drho = 1e-3, dz = 1e-6;

void setup_state(std::vector<shared_ptr<HelmholtzEOSMixtureBackend> > & HEOS, std::size_t Ncomp, double increment, x_N_dependency_flag xN_flag = XN_INDEPENDENT)
{
    std::vector<std::string> names(Ncomp);
    std::vector<CoolPropDbl> z(Ncomp);
    if (Ncomp == 2){
        names[0] = "Methane"; names[1] = "H2S";
        z[0] = 0.4; z[1] = 0.6;
    }
    else if (Ncomp == 3){
        names[0] = "Ethane"; names[1] = "Propane"; names[2] = "Methane";
        z[0] = 0.3; z[1] = 0.4; z[2] = 0.3;
    }
    else if (Ncomp == 4){
        names[0] = "Ethane"; names[1] = "Propane"; names[2] = "Methane"; names[3] = "n-Butane";
        z[0] = 0.3; z[1] = 0.4; z[2] = 0.2; z[3] = 0.1;
    }
    for (std::size_t i = 0; i < HEOS.size(); ++i){
        std::vector<CoolPropDbl> zn = z;
        zn[i] += increment;
        if (xN_flag == XN_DEPENDENT){ zn[zn.size()-1] -= increment; }
        HEOS[i].reset(new HelmholtzEOSMixtureBackend(names));
        HEOS[i]->specify_phase(iphase_gas);
        HEOS[i]->set_mole_fractions(zn);
    }
}

// Set up all the fluids
void connect_fluids(){
    if (!fluids_set){
        HEOS.resize(Ncomp_max);
        HEOS_plusT_constrho.resize(Ncomp_max);
        HEOS_minusT_constrho.resize(Ncomp_max);
        HEOS_plusrho_constT.resize(Ncomp_max);
        HEOS_minusrho_constT.resize(Ncomp_max);
        HEOS_plusz_xNindep.resize(Ncomp_max);
        HEOS_minusz_xNindep.resize(Ncomp_max);
        HEOS_plusz_consttaudelta_xNindep.resize(Ncomp_max);
        HEOS_minusz_consttaudelta_xNindep.resize(Ncomp_max);
        HEOS_plusz_xNdep.resize(Ncomp_max);
        HEOS_minusz_xNdep.resize(Ncomp_max);
        HEOS_plusz_consttaudelta_xNdep.resize(Ncomp_max);
        HEOS_minusz_consttaudelta_xNdep.resize(Ncomp_max);

        for (std::size_t Ncomp = 2; Ncomp <= 4; ++Ncomp){    
            HEOS[Ncomp].resize(1);
            HEOS_plusT_constrho[Ncomp].resize(1);
            HEOS_minusT_constrho[Ncomp].resize(1);
            HEOS_plusrho_constT[Ncomp].resize(1);
            HEOS_minusrho_constT[Ncomp].resize(1);
            HEOS_plusz_xNindep[Ncomp].resize(Ncomp);
            HEOS_minusz_xNindep[Ncomp].resize(Ncomp);
            HEOS_plusz_consttaudelta_xNindep[Ncomp].resize(Ncomp);
            HEOS_minusz_consttaudelta_xNindep[Ncomp].resize(Ncomp);
            HEOS_plusz_xNdep[Ncomp].resize(Ncomp);
            HEOS_minusz_xNdep[Ncomp].resize(Ncomp);
            HEOS_plusz_consttaudelta_xNdep[Ncomp].resize(Ncomp);
            HEOS_minusz_consttaudelta_xNdep[Ncomp].resize(Ncomp);

            setup_state(HEOS[Ncomp], Ncomp, 0);
            setup_state(HEOS_plusT_constrho[Ncomp], Ncomp, 0);
            setup_state(HEOS_minusT_constrho[Ncomp], Ncomp, 0);
            setup_state(HEOS_plusrho_constT[Ncomp], Ncomp, 0);
            setup_state(HEOS_minusrho_constT[Ncomp], Ncomp, 0);
            setup_state(HEOS_plusz_xNindep[Ncomp], Ncomp, dz);
            setup_state(HEOS_minusz_xNindep[Ncomp], Ncomp, -dz);
            setup_state(HEOS_plusz_consttaudelta_xNindep[Ncomp], Ncomp, dz);
            setup_state(HEOS_minusz_consttaudelta_xNindep[Ncomp], Ncomp, -dz);
            setup_state(HEOS_plusz_xNdep[Ncomp], Ncomp, dz, XN_DEPENDENT);
            setup_state(HEOS_minusz_xNdep[Ncomp], Ncomp, -dz, XN_DEPENDENT);
            setup_state(HEOS_plusz_consttaudelta_xNdep[Ncomp], Ncomp, dz, XN_DEPENDENT);
            setup_state(HEOS_minusz_consttaudelta_xNdep[Ncomp], Ncomp, -dz, XN_DEPENDENT);

            HEOS[Ncomp][0]->update(DmolarT_INPUTS, rho1, T1);
            HEOS_plusT_constrho[Ncomp][0]->update(DmolarT_INPUTS, rho1, T1 + dT);
            HEOS_minusT_constrho[Ncomp][0]->update(DmolarT_INPUTS, rho1, T1 - dT);
            HEOS_plusrho_constT[Ncomp][0]->update(DmolarT_INPUTS, rho1 + drho, T1);
            HEOS_minusrho_constT[Ncomp][0]->update(DmolarT_INPUTS, rho1 - drho, T1);

            for (std::size_t i = 0; i < Ncomp; ++i){
                
                HEOS_plusz_xNindep[Ncomp][i]->update(DmolarT_INPUTS, rho1, T1);
                HEOS_minusz_xNindep[Ncomp][i]->update(DmolarT_INPUTS, rho1, T1);
                HEOS_plusz_xNdep[Ncomp][i]->update(DmolarT_INPUTS, rho1, T1);
                HEOS_minusz_xNdep[Ncomp][i]->update(DmolarT_INPUTS, rho1, T1);

                HEOS_plusz_consttaudelta_xNindep[Ncomp][i]->calc_reducing_state();
                SimpleState red = HEOS_plusz_consttaudelta_xNindep[Ncomp][i]->get_reducing_state();
                HEOS_plusz_consttaudelta_xNindep[Ncomp][i]->update(DmolarT_INPUTS, red.rhomolar*HEOS[Ncomp][0]->delta(), red.T/HEOS[Ncomp][0]->tau());
                HEOS_plusz_consttaudelta_xNdep[Ncomp][i]->calc_reducing_state();
                red = HEOS_plusz_consttaudelta_xNdep[Ncomp][i]->get_reducing_state();
                HEOS_plusz_consttaudelta_xNdep[Ncomp][i]->update(DmolarT_INPUTS, red.rhomolar*HEOS[Ncomp][0]->delta(), red.T/HEOS[Ncomp][0]->tau());

                HEOS_minusz_consttaudelta_xNindep[Ncomp][i]->calc_reducing_state();
                red = HEOS_minusz_consttaudelta_xNindep[Ncomp][i]->get_reducing_state();
                HEOS_minusz_consttaudelta_xNindep[Ncomp][i]->update(DmolarT_INPUTS, red.rhomolar*HEOS[Ncomp][0]->delta(), red.T/HEOS[Ncomp][0]->tau());
                HEOS_minusz_consttaudelta_xNdep[Ncomp][i]->calc_reducing_state();
                red = HEOS_minusz_consttaudelta_xNdep[Ncomp][i]->get_reducing_state();
                HEOS_minusz_consttaudelta_xNdep[Ncomp][i]->update(DmolarT_INPUTS, red.rhomolar*HEOS[Ncomp][0]->delta(), red.T/HEOS[Ncomp][0]->tau());
            }
        }
        fluids_set = true;
    }
}

TEST_CASE("Mixture derivative checks", "[mixtures],[mixture_derivs]")
{
    connect_fluids();

    for (std::size_t Ncomp = 2; Ncomp <= 4; Ncomp++)
    {
        std::ostringstream ss00;
        ss00 << "Mixture with " << Ncomp << " components";
        SECTION(ss00.str(),"")
        {
            x_N_dependency_flag xN_flag;
            for (std::size_t xNxN = 1; xNxN > 0; xNxN--){
                std::ostringstream ss000;
                std::string xN_string;
                if (xNxN == 0){
                    xN_flag = XN_DEPENDENT;
                    xN_string = "Mole fractions dependent";
                }
                else{
                    xN_flag = XN_INDEPENDENT;
                    xN_string = "Mole fractions independent";
                }
            
                ss000 << xN_string;
                SECTION(ss000.str(),"")
                {                    
                    
                    HelmholtzEOSMixtureBackend &rHEOS = *(HEOS[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_plusT_constrho = *(HEOS_plusT_constrho[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_minusT_constrho = *(HEOS_minusT_constrho[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_plusrho_constT = *(HEOS_plusrho_constT[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_minusrho_constT = *(HEOS_minusrho_constT[Ncomp][0].get());

                    const std::vector<CoolPropDbl> &z = rHEOS.get_mole_fractions();
                    
                    SECTION("adj(Mstar)", "")
                    {
                        if (Ncomp == 2){
                        Eigen::MatrixXd Mstar = MixtureDerivatives::Mstar(rHEOS, xN_flag);
                        Eigen::MatrixXd analytic = adjugate(Mstar);
                        Eigen::MatrixXd numeric(2,2);
                        numeric << Mstar(1,1), -Mstar(0,1), -Mstar(1,0), Mstar(0,0);
                        double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CAPTURE(Mstar);
                        CHECK(err < 1e-8);
                        }
                    }
                    SECTION("d(adj(Lstar))/dDelta", "")
                    {
                        Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(rHEOS, xN_flag);
                        Eigen::MatrixXd dLstar_dDelta = MixtureDerivatives::dLstar_dX(rHEOS, xN_flag, CoolProp::iDelta);
                        Eigen::MatrixXd analytic = adjugate_derivative(Lstar, dLstar_dDelta);

                        Eigen::MatrixXd adj_Lstar_plus = adjugate(MixtureDerivatives::Lstar(rHEOS_plusrho_constT, xN_flag)); double delta1 = rHEOS_plusrho_constT.delta();
                        Eigen::MatrixXd adj_Lstar_minus = adjugate(MixtureDerivatives::Lstar(rHEOS_minusrho_constT, xN_flag)); double delta2 = rHEOS_minusrho_constT.delta();
                        
                        Eigen::MatrixXd numeric = (adj_Lstar_plus - adj_Lstar_minus)/(delta1-delta2);
                        double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CHECK(err < 1e-8);
                    }
                    SECTION("d(M1)/dTau", "")
                    {
                        Eigen::MatrixXd Mstar = MixtureDerivatives::Mstar(rHEOS, xN_flag);
                        Eigen::MatrixXd dMstar_dTau = MixtureDerivatives::dMstar_dX(rHEOS, xN_flag, CoolProp::iTau);
                        Eigen::MatrixXd adjM = adjugate(Mstar);
                        double analytic = (adjM*dMstar_dTau).trace();
                        
                        double detMstar_plus = MixtureDerivatives::Mstar(rHEOS_plusT_constrho, xN_flag).determinant(); double tau1 = rHEOS_plusT_constrho.tau();
                        double detMstar_minus = MixtureDerivatives::Mstar(rHEOS_minusT_constrho, xN_flag).determinant(); double tau2 = rHEOS_minusT_constrho.tau();

                        double numeric = (detMstar_plus - detMstar_minus)/(tau1-tau2);
                        double err = mix_deriv_err_func(numeric, analytic);
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CAPTURE(dMstar_dTau);
                        CAPTURE(adjM);
                        CHECK(err < 1e-8);
                    }
                    SECTION("d(M1)/dDelta", "")
                    {
                        Eigen::MatrixXd Mstar = MixtureDerivatives::Mstar(rHEOS, xN_flag);
                        Eigen::MatrixXd dMstar_dDelta = MixtureDerivatives::dMstar_dX(rHEOS, xN_flag, CoolProp::iDelta);
                        double analytic = (adjugate(Mstar)*dMstar_dDelta).trace();
                        
                        double detMstar_plus = MixtureDerivatives::Mstar(rHEOS_plusrho_constT, xN_flag).determinant(); double delta1 = rHEOS_plusrho_constT.delta();
                        double detMstar_minus = MixtureDerivatives::Mstar(rHEOS_minusrho_constT, xN_flag).determinant(); double delta2 = rHEOS_minusrho_constT.delta();
                        
                        double numeric = (detMstar_plus - detMstar_minus)/(delta1-delta2);
                        double err = mix_deriv_err_func(numeric, analytic);
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CHECK(err < 1e-8);
                    }
                    SECTION("d(L1)/dDelta", "")
                    {
                        Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(rHEOS, xN_flag);
                        Eigen::MatrixXd dLstar_dDelta = MixtureDerivatives::dLstar_dX(rHEOS, xN_flag, CoolProp::iDelta);
                        double analytic = (adjugate(Lstar)*dLstar_dDelta).trace();
                        
                        double detLstar_plus = MixtureDerivatives::Lstar(rHEOS_plusrho_constT, xN_flag).determinant(); double delta1 = rHEOS_plusrho_constT.delta();
                        double detLstar_minus = MixtureDerivatives::Lstar(rHEOS_minusrho_constT, xN_flag).determinant(); double delta2 = rHEOS_minusrho_constT.delta();
                        
                        double numeric = (detLstar_plus - detLstar_minus)/(delta1-delta2);
                        double err = mix_deriv_err_func(numeric, analytic);
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CHECK(err < 1e-8);
                    }

                    SECTION("adj(Lstar)", "")
                    {
                        if (Ncomp == 2){
                        Eigen::MatrixXd Lstar = MixtureDerivatives::Lstar(rHEOS, xN_flag);
                        Eigen::MatrixXd analytic = adjugate(Lstar);
                        Eigen::MatrixXd numeric(2,2);
                        numeric << Lstar(1,1), -Lstar(0,1), -Lstar(1,0), Lstar(0,0);
                        double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
                        CAPTURE(numeric);
                        CAPTURE(analytic);
                        CAPTURE(Lstar);
                        CHECK(err < 1e-8);
                        }
                    }
					SECTION("dLstar_Tau", "")
					{
						Eigen::MatrixXd analytic = MixtureDerivatives::dLstar_dX(rHEOS, xN_flag, CoolProp::iTau);
						Eigen::MatrixXd m1 = MixtureDerivatives::Lstar(rHEOS_plusT_constrho, xN_flag); double tau1 = rHEOS_plusT_constrho.tau();
						Eigen::MatrixXd m2 = MixtureDerivatives::Lstar(rHEOS_minusT_constrho, xN_flag); double tau2 = rHEOS_minusT_constrho.tau();
						Eigen::MatrixXd numeric = (m1 - m2)/(tau1 - tau2);
						double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
						CAPTURE(numeric);
                        CAPTURE(analytic);
						CHECK(err < 1e-8);
					}
					SECTION("dLstar_dDelta", "")
                    {
						Eigen::MatrixXd analytic = MixtureDerivatives::dLstar_dX(rHEOS, xN_flag, CoolProp::iDelta);
						Eigen::MatrixXd m1 = MixtureDerivatives::Lstar(rHEOS_plusrho_constT, xN_flag); double delta1 = rHEOS_plusrho_constT.delta();
						Eigen::MatrixXd m2 = MixtureDerivatives::Lstar(rHEOS_minusrho_constT, xN_flag); double delta2 = rHEOS_minusrho_constT.delta();
						Eigen::MatrixXd numeric = (m1 - m2)/(delta1 - delta2);
						double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
						CAPTURE(numeric);
						CAPTURE(analytic);
						CHECK(err < 1e-8);
					}
					SECTION("dMstar_dTau", "")
					{
						Eigen::MatrixXd analytic = MixtureDerivatives::dMstar_dX(rHEOS, xN_flag, CoolProp::iTau);
                        Eigen::MatrixXd m1 = MixtureDerivatives::Mstar(rHEOS_plusT_constrho, xN_flag); double tau1 = rHEOS_plusT_constrho.tau();
                        Eigen::MatrixXd m2 = MixtureDerivatives::Mstar(rHEOS_minusT_constrho, xN_flag); double tau2 = rHEOS_minusT_constrho.tau();
						Eigen::MatrixXd numeric = (m1 - m2)/(tau1 - tau2);
						double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
						CAPTURE(numeric);
						CAPTURE(analytic);
						CHECK(err < 1e-8);
					}
					std::ostringstream ss3o7; ss3o7 << "dMstar_dDelta";
					SECTION(ss3o7.str(), "")
					{
						Eigen::MatrixXd analytic = MixtureDerivatives::dMstar_dX(rHEOS, xN_flag, CoolProp::iDelta);
						Eigen::MatrixXd m1 = MixtureDerivatives::Mstar(rHEOS_plusrho_constT, xN_flag); double delta1 = rHEOS_plusrho_constT.delta();
						Eigen::MatrixXd m2 = MixtureDerivatives::Mstar(rHEOS_minusrho_constT, xN_flag); double delta2 = rHEOS_minusrho_constT.delta();
						Eigen::MatrixXd numeric = (m1 - m2)/(delta1 - delta2);
						double err = ((analytic-numeric).array()/analytic.array()).cwiseAbs().sum()/Ncomp/Ncomp;
						CAPTURE(numeric);
						CAPTURE(analytic);
						CHECK(err < 1e-8);
					}
                
                    // These ones only require the i index
                    for (std::size_t i = 0; i< Ncomp; ++i)
                    {
                        HelmholtzEOSMixtureBackend & rHEOS_pluszi = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_xNindep[Ncomp][i].get()) : *(HEOS_plusz_xNdep[Ncomp][i].get());
                        HelmholtzEOSMixtureBackend & rHEOS_minuszi = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_xNindep[Ncomp][i].get()) : *(HEOS_minusz_xNdep[Ncomp][i].get());

                        HelmholtzEOSMixtureBackend &rHEOS_pluszi_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_consttaudelta_xNindep[Ncomp][i].get()) : *(HEOS_plusz_consttaudelta_xNdep[Ncomp][i].get());
                        HelmholtzEOSMixtureBackend &rHEOS_minuszi_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_consttaudelta_xNindep[Ncomp][i].get()) : *(HEOS_minusz_consttaudelta_xNdep[Ncomp][i].get());
                        
                        std::ostringstream ss0;
                        ss0 << "dln_fugacity_i_dT__constrho_n, i=" << i;
                        SECTION(ss0.str(),"")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_i_dT__constrho_n(rHEOS, i, xN_flag);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_plusT_constrho, i, xN_flag));
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minusT_constrho, i, xN_flag));
                            double numeric = (v1 - v2)/(2*dT);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss0a;
                        ss0a << "dln_fugacity_i_dp__constT_n, i=" << i;
                        SECTION(ss0a.str(),"")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_i_dp__constT_n(rHEOS, i, xN_flag);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_plusrho_constT, i, xN_flag)), p1 = rHEOS_plusrho_constT.p();
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minusrho_constT, i, xN_flag)), p2 = rHEOS_minusrho_constT.p();
                            double numeric = (v1 - v2)/(p1 - p2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-6);
                        }
                        std::ostringstream ss0b;
                        ss0b << "dln_fugacity_i_drho__constT_n, i=" << i;
                        SECTION(ss0b.str(), "")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_i_drho__constT_n(rHEOS, i, xN_flag);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_plusrho_constT, i, xN_flag));
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minusrho_constT, i, xN_flag));
                            double numeric = (v1 - v2)/(2*dT);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-6);
                        }
                        std::ostringstream ss2;
                        ss2 << "dln_fugacity_coefficient_dp__constT_n, i=" << i;
                        SECTION(ss2.str(), "")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS_plusrho_constT, i, xN_flag), p1 = rHEOS_plusrho_constT.p();
                            double v2 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS_minusrho_constT, i, xN_flag), p2 = rHEOS_minusrho_constT.p();
                            double numeric = (v1 - v2)/(p1 - p2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CHECK(err < 1e-6);
                        }
                        /*std::ostringstream ss1;
                        ss1 << "dln_fugacity_coefficient_dT__constp_n, i=" << i;
                        SECTION(ss1.str(), "")
                        {
                            double T1 = 300, dT = 1e-3;
                            rHEOS.specify_phase(iphase_gas);
                            
                            rHEOS.update(PT_INPUTS, 101325, T1);
                            double analytic = MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(rHEOS, i, xN_flag);
                            
                            rHEOS.update(PT_INPUTS, 101325, T1 + dT);
                            double v1 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag);
                            rHEOS.update(PT_INPUTS, 101325, T1 - dT);
                            double v2 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag);
                            
                            double numeric = (v1 - v2)/(2*dT);
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-6);
                        }
                        std::ostringstream ss1a;
                        ss1a << "dln_fugacity_i_dT__constp_n, i=" << i;
                        SECTION(ss1a.str(), "")
                        {
                            double T1 = 300, dT = 1e-3;
                            rHEOS.specify_phase(iphase_gas);
                            
                            rHEOS.update(PT_INPUTS, 101325, T1);
                            double analytic = MixtureDerivatives::dln_fugacity_i_dT__constp_n(rHEOS, i, xN_flag);
                            
                            rHEOS.update(PT_INPUTS, 101325, T1 + dT);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS, i, xN_flag));
                            rHEOS.update(PT_INPUTS, 101325, T1 - dT);
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS, i, xN_flag));
                            
                            double numeric = (v1 - v2)/(2*dT);
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-6);
                        }
                        */
                        std::ostringstream ss3;
                        ss3 << "d_ndalphardni_dDelta, i=" << i;
                        SECTION(ss3.str(), "")
                        {
                            double analytic = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_plusrho_constT, i, xN_flag),  delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3a;
                        ss3a << "d2alphar_dxi_dDelta, i=" << i;
                        SECTION(ss3a.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusrho_constT.residual_helmholtz->dalphar_dxi(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = rHEOS_minusrho_constT.residual_helmholtz->dalphar_dxi(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CHECK(err < 1e-6);
                        }
                        std::ostringstream ss4a;
                        ss4a << "d2alphar_dxi_dTau, i=" << i;
                        SECTION(ss4a.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d2alphar_dxi_dTau(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->dalphar_dxi(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->dalphar_dxi(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss4;
                        ss4 << "d_ndalphardni_dTau, i=" << i;
                        SECTION(ss4.str(), "")
                        {
                            double analytic = MixtureDerivatives::d_ndalphardni_dTau(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss5;
                        ss5 << "dpdxj__constT_V_xi, i=" << i;
                        SECTION(ss5.str(), "")
                        {
                            double analytic = MixtureDerivatives::dpdxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.p();
                            double v2 = rHEOS_minuszi.p();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-6);
                        }       
                        std::ostringstream ss5a;
                        ss5a << "dtaudxj__constT_V_xi, i=" << i;
                        SECTION(ss5a.str(), "")
                        {
                            double analytic = MixtureDerivatives::dtau_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.tau();
                            double v2 = rHEOS_minuszi.tau();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }       
                        std::ostringstream ss5b;
                        ss5b << "ddeltadxj__constT_V_xi, i=" << i;
                        SECTION(ss5b.str(), "")
                        {
                            double analytic = MixtureDerivatives::ddelta_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.delta();
                            double v2 = rHEOS_minuszi.delta();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }       
                        std::ostringstream ss6;
                        ss6 << "d_dalpharddelta_dxj__constT_V_xi, i=" << i;
                        SECTION(ss6.str(), "")
                        {
                            double analytic = MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.dalphar_dDelta();
                            double v2 = rHEOS_minuszi.dalphar_dDelta();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }    
                        std::ostringstream ss7;
                        ss7 << "dTrdxi__constxj, i=" << i;
                        SECTION(ss7.str(), "")
                        {
                            double analytic = rHEOS.Reducing->dTrdxi__constxj(rHEOS.get_mole_fractions(), i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->Tr(rHEOS_pluszi.get_mole_fractions());
                            double v2 = rHEOS_minuszi.Reducing->Tr(rHEOS_minuszi.get_mole_fractions());
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss8;
                        ss8 << "drhormolardxi__constxj, i=" << i;
                        SECTION(ss8.str(), "")
                        {
                            double analytic = rHEOS.Reducing->drhormolardxi__constxj(rHEOS.get_mole_fractions(), i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->rhormolar(rHEOS_pluszi.get_mole_fractions());
                            double v2 = rHEOS_minuszi.Reducing->rhormolar(rHEOS_minuszi.get_mole_fractions());
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3c;
                        ss3c << "d2Trdxi2__constxj, i=" << i;
                        SECTION(ss3c.str(), "")
                        {
                            if (i == Ncomp-1){ break; }
                            double analytic = rHEOS.Reducing->d2Trdxi2__constxj(z, i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->dTrdxi__constxj(rHEOS_pluszi.get_mole_fractions(), i, xN_flag);
                            double v2 = rHEOS_minuszi.Reducing->dTrdxi__constxj(rHEOS_minuszi.get_mole_fractions(), i, xN_flag);
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3d;
                        ss3d << "dalphar_dxi, i=" << i;
                        SECTION(ss3d.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->dalphar_dxi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi_consttaudelta.alphar();
                            double v2 = rHEOS_minuszi_consttaudelta.alphar();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3e; ss3e << "d3alphar_dxi_dDelta2, i=" << i;
                        SECTION(ss3e.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d3alphar_dxi_dDelta2(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusrho_constT.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = rHEOS_minusrho_constT.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3f; ss3f << "d3alphar_dxi_dTau2, i=" << i;
                        SECTION(ss3f.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d3alphar_dxi_dTau2(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->d2alphar_dxi_dTau(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->d2alphar_dxi_dTau(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3ff; ss3ff << "d4alphar_dxi_dTau3, i=" << i;
                        SECTION(ss3ff.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dTau3(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->d3alphar_dxi_dTau2(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->d3alphar_dxi_dTau2(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }

                        std::ostringstream ss3g; ss3g << "d3alphar_dxi_dDelta_Tau, i=" << i;
                        SECTION(ss3g.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d3alphar_dxi_dDelta_dTau(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3gg; ss3gg << "d4alphar_dxi_dDelta_Tau2, i=" << i;
                        SECTION(ss3gg.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dDelta_dTau2(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->d3alphar_dxi_dDelta_dTau(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->d3alphar_dxi_dDelta_dTau(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3ggg; ss3ggg << "d4alphar_dxi_dDelta2_Tau, i=" << i;
                        SECTION(ss3ggg.str(), "")
                        {
                            double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dDelta2_dTau(rHEOS, i, xN_flag);
                            double v1 = rHEOS_plusT_constrho.residual_helmholtz->d3alphar_dxi_dDelta2(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = rHEOS_minusT_constrho.residual_helmholtz->d3alphar_dxi_dDelta2(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        
                        std::ostringstream ss3h; ss3h << "d2_ndalphardni_dTau2, i=" << i;
                        SECTION(ss3h.str(), "")
                        {
                            double analytic = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3hh; ss3hh << "d3_ndalphardni_dTau3, i=" << i;
                        SECTION(ss3hh.str(), "")
                        {
                            double analytic = MixtureDerivatives::d3_ndalphardni_dTau3(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3i; ss3i << "d2_ndalphardni_dDelta_dTau, i=" << i;
                        SECTION(ss3i.str(), "")
                        {
                            double analytic = MixtureDerivatives::d2_ndalphardni_dDelta_dTau(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3ii; ss3ii << "d3_ndalphardni_dDelta_dTau2, i=" << i;
                        SECTION(ss3ii.str(), "")
                        {
                            double analytic = MixtureDerivatives::d3_ndalphardni_dDelta_dTau2(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3iii; ss3iii << "d3_ndalphardni_dDelta2_dTau, i=" << i;
                        SECTION(ss3iii.str(), "")
                        {
                            double analytic = MixtureDerivatives::d3_ndalphardni_dDelta2_dTau(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d2_ndalphardni_dDelta_dTau(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::d2_ndalphardni_dDelta_dTau(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3j; ss3j << "d2_ndalphardni_dDelta2, i=" << i;
                        SECTION(ss3j.str(), "")
                        {
                            double analytic = MixtureDerivatives::d2_ndalphardni_dDelta2(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        std::ostringstream ss3k; ss3k << "d3_ndalphardni_dDelta3, i=" << i;
                        SECTION(ss3k.str(), "")
                        {
                            double analytic = MixtureDerivatives::d3_ndalphardni_dDelta3(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::d2_ndalphardni_dDelta2(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::d2_ndalphardni_dDelta2(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = mix_deriv_err_func(numeric, analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        
                        // These derivatives depend on both the i and j indices
                        for (std::size_t j = 0; j < Ncomp; ++j){
                            if (xN_flag == XN_DEPENDENT && j == Ncomp){ continue; }
                            HelmholtzEOSMixtureBackend & rHEOS_pluszj = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_xNindep[Ncomp][j].get()) : *(HEOS_plusz_xNdep[Ncomp][j].get());
                            HelmholtzEOSMixtureBackend & rHEOS_minuszj = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_xNindep[Ncomp][j].get()) : *(HEOS_minusz_xNdep[Ncomp][j].get());

                            HelmholtzEOSMixtureBackend &rHEOS_pluszj_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_consttaudelta_xNindep[Ncomp][j].get()) : *(HEOS_plusz_consttaudelta_xNdep[Ncomp][j].get());
                            HelmholtzEOSMixtureBackend &rHEOS_minuszj_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_consttaudelta_xNindep[Ncomp][j].get()) : *(HEOS_minusz_consttaudelta_xNdep[Ncomp][j].get());
                        
                            std::ostringstream ss1a;
                            ss1a << "dln_fugacity_dxj__constT_rho_xi, i=" << i << ", j=" << j;
                            SECTION(ss1a.str(), "")
                            {
                                if (xN_flag == XN_INDEPENDENT){continue;}
                                double analytic = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(rHEOS, i, j, xN_flag);
                                double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_pluszj, i, xN_flag));
                                double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minuszj, i, xN_flag));
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss2;
                            ss2 << "d_ndTrdni_dxj, i=" << i << ", j=" << j;
                            SECTION(ss2.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS.get_mole_fractions(), i, j, xN_flag);
                                double v1 = rHEOS_pluszj.Reducing->ndTrdni__constnj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS_minuszj.Reducing->ndTrdni__constnj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss4;
                            ss4 << "d_ndrhomolarrdni_dxj, i=" << i << ", j=" << j;
                            SECTION(ss4.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS.get_mole_fractions(), i, j, xN_flag);
                                double v1 = rHEOS_pluszj.Reducing->ndrhorbardni__constnj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS_minuszj.Reducing->ndrhorbardni__constnj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3;
                            ss3 << "d_ndalphardni_dxj__constT_V_xi, i=" << i << ", j=" << j;
                            SECTION(ss3.str(), "")
                            {
                                double analytic = MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_pluszj, i, xN_flag);
                                double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_minuszj, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3a;
                            ss3a << "d2alphardxidxj, i=" << i << ", j=" << j;
                            SECTION(ss3a.str(), "")
                            {
                                double analytic = rHEOS.residual_helmholtz->d2alphardxidxj(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->dalphar_dxi(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->dalphar_dxi(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){break;}
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-6);
                            }
                            std::ostringstream ss3b;
                            ss3b << "d2Trdxidxj, i=" << i << ", j=" << j;
                            SECTION(ss3b.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d2Trdxidxj(z, i, j, xN_flag);
                                double v1 = rHEOS.Reducing->dTrdxi__constxj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS.Reducing->dTrdxi__constxj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3c;
                            ss3c << "d_PSI_rho_dxj, i=" << i << ", j=" << j;
                            SECTION(ss3c.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d_PSI_rho_dxj(z, i, j, xN_flag);
                                double v1 = rHEOS.Reducing->PSI_rho(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS.Reducing->PSI_rho(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3d; ss3d << "d_PSI_T_dxj, i=" << i << ", j=" << j;
                            SECTION(ss3d.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d_PSI_T_dxj(z, i, j, xN_flag);
                                double v1 = rHEOS.Reducing->PSI_T(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS.Reducing->PSI_T(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3k; ss3k << "d2_ndalphardni_dxj_dDelta, i=" << i << ", j=" << j;
                            SECTION(ss3k.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d2_ndalphardni_dxj_dDelta__consttau_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3kk; ss3kk << "d3_ndalphardni_dxj_dDelta2__consttau_xi, i=" << i << ", j=" << j;
                            SECTION(ss3kk.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d3_ndalphardni_dxj_dDelta2__consttau_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d2_ndalphardni_dDelta2(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::d2_ndalphardni_dDelta2(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3l; ss3l << "d2_ndalphardni_dxj_dTau, i=" << i << ", j=" << j;
                            SECTION(ss3l.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d2_ndalphardni_dxj_dTau__constdelta_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::d_ndalphardni_dTau(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3ll; ss3ll << "d3_ndalphardni_dxj_dTau2__constdelta_xi, i=" << i << ", j=" << j;
                            SECTION(ss3ll.str(), "")
                            {
                                double analytic = MixtureDerivatives::d3_ndalphardni_dxj_dTau2__constdelta_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::d2_ndalphardni_dTau2(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3lll; ss3lll << "d3_ndalphardni_dxj_dDelta_dTau__constxi, i=" << i << ", j=" << j;
                            SECTION(ss3lll.str(), "")
                            {
                                double analytic = MixtureDerivatives::d3_ndalphardni_dxj_dDelta_dTau__constxi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d2_ndalphardni_dDelta_dTau(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::d2_ndalphardni_dDelta_dTau(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3m; ss3m << "d3alphar_dxi_dxj_dDelta, i=" << i << ", j=" << j;
                            SECTION(ss3m.str(), "")
                            {
                                double analytic = rHEOS.residual_helmholtz->d3alphar_dxi_dxj_dDelta(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->d2alphar_dxi_dDelta(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3mm; ss3mm << "d4alphar_dxi_dxj_dDelta2, i=" << i << ", j=" << j;
                            SECTION(ss3mm.str(), "")
                            {
                                double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta2(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dDelta2(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dDelta2(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3mmm; ss3mmm << "d4alphar_dxi_dxj_dDelta_dTau, i=" << i << ", j=" << j;
                            SECTION(ss3mmm.str(), "")
                            {
                                double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dxj_dDelta_dTau(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dDelta_dTau(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dDelta_dTau(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3n; ss3n << "d3alphar_dxi_dxj_dTau, i=" << i << ", j=" << j;
                            SECTION(ss3n.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = rHEOS.residual_helmholtz->d3alphar_dxi_dxj_dTau(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->d2alphar_dxi_dTau(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->d2alphar_dxi_dTau(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3nn; ss3nn << "d4alphar_dxi_dxj_dTau2, i=" << i << ", j=" << j;
                            SECTION(ss3nn.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = rHEOS.residual_helmholtz->d4alphar_dxi_dxj_dTau2(rHEOS, i, j, xN_flag);
                                double v1 = rHEOS_pluszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dTau2(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = rHEOS_minuszj_consttaudelta.residual_helmholtz->d3alphar_dxi_dTau2(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            
                            std::ostringstream ss3o; ss3o << "d_nd_ndalphardni_dnj_dDelta__consttau_x, i=" << i << ", j=" << j;
                            SECTION(ss3o.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_nd_ndalphardni_dnj_dDelta__consttau_x(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_plusrho_constT, i, j, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                                double v2 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_minusrho_constT, i, j, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                                double numeric = (v1 - v2)/(delta1 - delta2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3o1; ss3o1 << "d_ndln_fugacity_i_dnj_ddelta__consttau_x, i=" << i << ", j=" << j;
                            SECTION(ss3o1.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_plusrho_constT, i, j, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                                double v2 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_minusrho_constT, i, j, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                                double numeric = (v1 - v2)/(delta1 - delta2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3o2; ss3o2 << "d_ndln_fugacity_i_dnj_dtau__constdelta_x, i=" << i << ", j=" << j;
                            SECTION(ss3o2.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                                double v2 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                                double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }

                            std::ostringstream ss3p; ss3p << "d_nd_ndalphardni_dnj_dTau__constdelta_x, i=" << i << ", j=" << j;
                            SECTION(ss3p.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_nd_ndalphardni_dnj_dTau__constdelta_x(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                                double v2 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                                double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3q; ss3q << "d_nddeltadni_dxj__constdelta_tau, i=" << i << ", j=" << j;
                            SECTION(ss3q.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_nddeltadni_dxj__constdelta_tau(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::nddeltadni__constT_V_nj(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::nddeltadni__constT_V_nj(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3r; ss3r << "d_ndtaudni_dxj__constdelta_tau, i=" << i << ", j=" << j;
                            SECTION(ss3r.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d_ndtaudni_dxj__constdelta_tau(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::ndtaudni__constT_V_nj(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::ndtaudni__constT_V_nj(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-6);
                            }
							std::ostringstream ss3s; ss3s << "d2_ndtaudni_dxj_dTau__constdelta, i=" << i << ", j=" << j;
							SECTION(ss3s.str(), "")
							{
								double analytic = MixtureDerivatives::d2_ndtaudni_dxj_dTau__constdelta(rHEOS, i, j, xN_flag);
								double v1 = MixtureDerivatives::d_ndtaudni_dxj__constdelta_tau(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
								double v2 = MixtureDerivatives::d_ndtaudni_dxj__constdelta_tau(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
								double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
								CAPTURE(numeric);
								CAPTURE(analytic);
								CHECK(err < 1e-7);
							}
							std::ostringstream ss3t; ss3t << "d_nAij_dX, i=" << i << ", j=" << j;
							SECTION(ss3t.str(), "")
							{
								double analytic = MixtureDerivatives::d_nAij_dX(rHEOS, i, j, xN_flag, CoolProp::iTau);
								double v1 = MixtureDerivatives::nAij(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
								double v2 = MixtureDerivatives::nAij(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
								double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
								CAPTURE(numeric);
								CAPTURE(analytic);
								CHECK(err < 1e-7);
							}
							std::ostringstream ss3o3; ss3o3 << "d_ndln_fugacity_i_dnj_dtau__constdelta_x, i=" << i << ", j=" << j;
                            SECTION(ss3o3.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d2_ndln_fugacity_i_dnj_dtau2__constdelta_x(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                                double v2 = MixtureDerivatives::d_ndln_fugacity_i_dnj_dtau__constdelta_x(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                                double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
							std::ostringstream ss3o4; ss3o4 << "d2_ndln_fugacity_i_dnj_ddelta_dtau__constx, i=" << i << ", j=" << j;
                            SECTION(ss3o4.str(), "")
                            {
                                if (i == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d2_ndln_fugacity_i_dnj_ddelta_dtau__constx(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(rHEOS_plusT_constrho, i, j, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                                double v2 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddelta__consttau_x(rHEOS_minusT_constrho, i, j, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                                double numeric = (v1 - v2)/(tau1 - tau2);
                                double err = mix_deriv_err_func(numeric, analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            
                            // These derivatives depend on i,j, and k indices
                            for (std::size_t k = 0; k < Ncomp; ++k){
                                HelmholtzEOSMixtureBackend & rHEOS_pluszk = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_xNindep[Ncomp][k].get()) : *(HEOS_plusz_xNdep[Ncomp][k].get());
                                HelmholtzEOSMixtureBackend & rHEOS_minuszk = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_xNindep[Ncomp][k].get()) : *(HEOS_minusz_xNdep[Ncomp][k].get());

                                HelmholtzEOSMixtureBackend &rHEOS_pluszk_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_consttaudelta_xNindep[Ncomp][k].get()) : *(HEOS_plusz_consttaudelta_xNdep[Ncomp][k].get());
                                HelmholtzEOSMixtureBackend &rHEOS_minuszk_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_consttaudelta_xNindep[Ncomp][k].get()) : *(HEOS_minusz_consttaudelta_xNdep[Ncomp][k].get());

                                std::ostringstream ss1;    ss1 << "d3Trdxidxjdxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss1.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d3Trdxidxjdxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d2Trdxidxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d2Trdxidxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-6);
                                }
                                std::ostringstream ss2;    ss2 << "d3rhormolardxidxjdxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss2.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d3rhormolardxidxjdxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d2rhormolardxidxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d2rhormolardxidxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss3; ss3 << "d2_ndTrdni_dxj_dxk__constxi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss3.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d2_ndTrdni_dxj_dxk__constxi(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss4; ss4 << "d2_ndrhorbardni_dxj_dxk__constxi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss4.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d2_ndrhorbardni_dxj_dxk__constxi(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss5; ss5 << "d2_PSI_T_dxj_dxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss5.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d2_PSI_T_dxj_dxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_PSI_T_dxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_PSI_T_dxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss6; ss6 << "d2_PSI_rho_dxj_dxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss6.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.Reducing->d2_PSI_rho_dxj_dxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_PSI_rho_dxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_PSI_rho_dxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss7; ss7 << "d2_ndalphardni_dxj_dxk__constdelta_tau_xi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss7.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = MixtureDerivatives::d2_ndalphardni_dxj_dxk__constdelta_tau_xi(rHEOS, i, j, k, xN_flag);
                                    double v1 = MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-6);
                                }
                                std::ostringstream ss7a; ss7a << "d3_ndalphardni_dxj_dxk_dDelta__consttau_xi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss7a.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = MixtureDerivatives::d3_ndalphardni_dxj_dxk_dDelta__consttau_xi(rHEOS, i, j, k, xN_flag);
                                    double v1 = MixtureDerivatives::d2_ndalphardni_dxj_dDelta__consttau_xi(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = MixtureDerivatives::d2_ndalphardni_dxj_dDelta__consttau_xi(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-6);
                                }
                                std::ostringstream ss7b; ss7b << "d3_ndalphardni_dxj_dxk_dTau__constdelta_xi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss7b.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = MixtureDerivatives::d3_ndalphardni_dxj_dxk_dTau__constdelta_xi(rHEOS, i, j, k, xN_flag);
                                    double v1 = MixtureDerivatives::d2_ndalphardni_dxj_dTau__constdelta_xi(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = MixtureDerivatives::d2_ndalphardni_dxj_dTau__constdelta_xi(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss8; ss8 << "d3alphardxidxjdxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss8.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = rHEOS.residual_helmholtz->d3alphardxidxjdxk(rHEOS, i, j, k, xN_flag);
                                    double v1 = rHEOS_pluszk_consttaudelta.residual_helmholtz->d2alphardxidxj(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = rHEOS_minuszk_consttaudelta.residual_helmholtz->d2alphardxidxj(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){ break; }
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-8);
                                }
                                
                                std::ostringstream ss9; ss9 << "d_nd_ndalphardni_dnj_dxk__consttau_delta, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss9.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = MixtureDerivatives::d_nd_ndalphardni_dnj_dxk__consttau_delta(rHEOS, i, j, k, xN_flag);
                                    double v1 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = MixtureDerivatives::nd_ndalphardni_dnj__constT_V(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){ break; }
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-6);
                                }
                                std::ostringstream ss10; ss10 << "d_ndln_fugacity_i_dnj_ddxk__consttau_delta, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss10.str(), "")
                                {
                                    if ((xN_flag == XN_DEPENDENT) && (i == Ncomp-1 || j == Ncomp-1 || k == Ncomp-1)){ break; }
                                    double analytic = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(rHEOS, i, j, k, xN_flag);
                                    double v1 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_pluszk_consttaudelta, i, j, xN_flag);
                                    double v2 = MixtureDerivatives::ndln_fugacity_i_dnj__constT_V_xi(rHEOS_minuszk_consttaudelta, i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){ break; }
                                    double err = mix_deriv_err_func(numeric, analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-8);
                                }
								std::ostringstream ss3o4; ss3o4 << "d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta, i=" << i << ", j=" << j << ", k=" << k;
								SECTION(ss3o4.str(), "")
								{
									if (i == Ncomp-1){ break; }
									double analytic = MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dTau__constdelta(rHEOS, i, j, k, xN_flag);
									double v1 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(rHEOS_plusT_constrho, i, j, k, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
									double v2 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(rHEOS_minusT_constrho, i, j, k, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
									double numeric = (v1 - v2)/(tau1 - tau2);
                                    double err = mix_deriv_err_func(numeric, analytic);
									CAPTURE(numeric);
									CAPTURE(analytic);
									CHECK(err < 1e-7);
								}
								std::ostringstream ss3oo4; ss3oo4 << "d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau, i=" << i << ", j=" << j << ", k=" << k;
								SECTION(ss3oo4.str(), "")
								{
									if (i == Ncomp-1){ break; }
									double analytic = MixtureDerivatives::d2_ndln_fugacity_i_dnj_dxk_dDelta__consttau(rHEOS, i, j, k, xN_flag);
									double v1 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(rHEOS_plusrho_constT, i, j, k, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
									double v2 = MixtureDerivatives::d_ndln_fugacity_i_dnj_ddxk__consttau_delta(rHEOS_minusrho_constT, i, j, k, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
									double numeric = (v1 - v2)/(delta1 - delta2);
                                    double err = mix_deriv_err_func(numeric, analytic);
									CAPTURE(numeric);
									CAPTURE(analytic);
									CHECK(err < 1e-7);
								}
								std::ostringstream ss3o5; ss3o5 << "d_n2Aijk_dTau, i=" << i << ", j=" << j << ", k=" << k;
								SECTION(ss3o5.str(), "")
								{
									double analytic = MixtureDerivatives::d_n2Aijk_dX(rHEOS, i, j, k, xN_flag, CoolProp::iTau);
									double v1 = MixtureDerivatives::n2Aijk(rHEOS_plusT_constrho, i, j, k, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
									double v2 = MixtureDerivatives::n2Aijk(rHEOS_minusT_constrho, i, j, k, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
									double numeric = (v1 - v2)/(tau1 - tau2);
                                    double err = mix_deriv_err_func(numeric, analytic);
									CAPTURE(numeric);
									CAPTURE(analytic);
									CHECK(err < 1e-7);
								}
								std::ostringstream ss3o6; ss3o6 << "d_n2Aijk_dDelta, i=" << i << ", j=" << j << ", k=" << k;
								SECTION(ss3o6.str(), "")
								{
									double analytic = MixtureDerivatives::d_n2Aijk_dX(rHEOS, i, j, k, xN_flag, CoolProp::iDelta);
									double v1 = MixtureDerivatives::n2Aijk(rHEOS_plusrho_constT, i, j, k, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
									double v2 = MixtureDerivatives::n2Aijk(rHEOS_minusrho_constT, i, j, k, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
									double numeric = (v1 - v2)/(delta1 - delta2);
                                    double err = mix_deriv_err_func(numeric, analytic);
									CAPTURE(numeric);
									CAPTURE(analytic);
									CHECK(err < 1e-7);
								}
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif


