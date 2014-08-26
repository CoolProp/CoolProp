#include "MixtureDerivatives.h"

namespace CoolProp{
    
long double MixtureDerivatives::dalphar_dxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS.components[i]->pEOS->baser(HEOS._tau, HEOS._delta) + HEOS.Excess.dalphar_dxi(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
}
long double MixtureDerivatives::d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS.components[i]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dTau(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
}
long double MixtureDerivatives::d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS.components[i]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dDelta(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
}
long double MixtureDerivatives::d2alphardxidxj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    return 0                           + HEOS.Excess.d2alphardxidxj(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i, j);
}

long double MixtureDerivatives::fugacity_i(HelmholtzEOSMixtureBackend &HEOS, std::size_t i){
    return HEOS.mole_fractions[i]*HEOS.rhomolar()*HEOS.gas_constant()*HEOS.T()*exp( dnalphar_dni__constT_V_nj(HEOS, i));
}
long double MixtureDerivatives::ln_fugacity_coefficient(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i)-log(1+HEOS._delta.pt()*HEOS.dalphar_dDelta());
}
long double MixtureDerivatives::dln_fugacity_coefficient_dT__constrho_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    double dtau_dT = -HEOS._tau.pt()/HEOS._T; //[1/K]
    return (HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i)-1/(1+HEOS._delta.pt()*HEOS.dalphar_dDelta())*(HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()))*dtau_dT;
}
long double MixtureDerivatives::dln_fugacity_coefficient_drho__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    double ddelta_drho = 1/HEOS._reducing.rhomolar; //[m^3/mol]
    return (HEOS.dalphar_dDelta() + d_ndalphardni_dDelta(HEOS, i)-1/(1+HEOS._delta.pt()*HEOS.dalphar_dDelta())*(HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta()))*ddelta_drho;
}
long double MixtureDerivatives::dnalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    // GERG Equation 7.42
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i);
}
long double MixtureDerivatives::d2nalphar_dni_dT(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return -HEOS._tau.pt()/HEOS._T*(HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i));
}
long double MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    double T = HEOS._reducing.T/HEOS._tau.pt();
    long double R_u = HEOS.gas_constant();
    return d2nalphar_dni_dT(HEOS, i) + 1/T-partial_molar_volume(HEOS, i)/(R_u*T)*dpdT__constV_n(HEOS);
}
long double MixtureDerivatives::partial_molar_volume(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return -ndpdni__constT_V_nj(HEOS, i)/ndpdV__constT_n(HEOS);
}

long double MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    // GERG equation 7.30
    long double R_u = HEOS.gas_constant();
    double partial_molar_volumeval = partial_molar_volume(HEOS, i); // [m^3/mol]
    double term1 = partial_molar_volumeval/(R_u*HEOS._T); // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
    double term2 = 1.0/HEOS.p();
    return term1 - term2;
}

long double MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    // Gernert 3.115
    long double R_u = HEOS.gas_constant();
    // partial molar volume is -dpdn/dpdV, so need to flip the sign here
    return d2nalphar_dxj_dni__constT_V(HEOS, j, i) - partial_molar_volume(HEOS, i)/(R_u*HEOS._T)*dpdxj__constT_V_xi(HEOS, j);
}
long double MixtureDerivatives::dpdxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j)
{
    // Gernert 3.130
    long double R_u = HEOS.gas_constant();
    return HEOS._rhomolar*R_u*HEOS._T*(ddelta_dxj__constT_V_xi(HEOS, j)*HEOS.dalphar_dDelta()+HEOS._delta.pt()*d_dalpharddelta_dxj__constT_V_xi(HEOS, j));
}

long double MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j)
{
    // Gernert Equation 3.134 (Catch test provided)
    return HEOS.d2alphar_dDelta2()*ddelta_dxj__constT_V_xi(HEOS, j)
         + HEOS.d2alphar_dDelta_dTau()*dtau_dxj__constT_V_xi(HEOS, j)
         + d2alphar_dxi_dDelta(HEOS, j);
}

long double MixtureDerivatives::dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j)
{
    //Gernert 3.119 (Catch test provided)
    return HEOS.dalphar_dDelta()*ddelta_dxj__constT_V_xi(HEOS, j)+HEOS.dalphar_dTau()*dtau_dxj__constT_V_xi(HEOS, j)+dalphar_dxi(HEOS, j);
}
long double MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    // Gernert 3.118
    return d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i,j)
          + ddelta_dxj__constT_V_xi(HEOS, j)*d_ndalphardni_dDelta(HEOS, i)
          + dtau_dxj__constT_V_xi(HEOS, j)*d_ndalphardni_dTau(HEOS, i);
}
long double MixtureDerivatives::ddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j)
{
    // Gernert 3.121 (Catch test provided)
    return -HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing.p->drhormolardxi__constxj(HEOS.mole_fractions,j);
}
long double MixtureDerivatives::dtau_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j)
{
    // Gernert 3.122 (Catch test provided)
    return 1/HEOS._T*HEOS.Reducing.p->dTrdxi__constxj(HEOS.mole_fractions,j);
}

long double MixtureDerivatives::dpdT__constV_n(HelmholtzEOSMixtureBackend &HEOS)
{
    long double R_u = HEOS.gas_constant();
    return HEOS._rhomolar*R_u*(1+HEOS._delta.pt()*HEOS.dalphar_dDelta()-HEOS._delta.pt()*HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau());
}
long double MixtureDerivatives::dpdrho__constT_n(HelmholtzEOSMixtureBackend &HEOS)
{
    long double R_u = HEOS.gas_constant();
    return R_u*HEOS._T*(1+2*HEOS._delta.pt()*HEOS.dalphar_dDelta()+pow(HEOS._delta.pt(),2)*HEOS.d2alphar_dDelta2());
}
long double MixtureDerivatives::ndpdV__constT_n(HelmholtzEOSMixtureBackend &HEOS)
{
    long double R_u = HEOS.gas_constant();
    return -pow(HEOS._rhomolar,2)*R_u*HEOS._T*(1+2*HEOS._delta.pt()*HEOS.dalphar_dDelta()+pow(HEOS._delta.pt(),2)*HEOS.d2alphar_dDelta2());
}
long double MixtureDerivatives::ndpdni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    // Eqn 7.64 and 7.63
    long double R_u = HEOS.gas_constant();
    double ndrhorbar_dni__constnj = HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i);
    double ndTr_dni__constnj = HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions,i);
    double summer = 0;
    for (unsigned int k = 0; k < HEOS.mole_fractions.size(); ++k)
    {
        summer += HEOS.mole_fractions[k]*d2alphar_dxi_dDelta(HEOS, k);
    }
    double nd2alphar_dni_dDelta = HEOS._delta.pt()*HEOS.d2alphar_dDelta2()*(1-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()/HEOS._reducing.T*ndTr_dni__constnj+d2alphar_dxi_dDelta(HEOS, i)-summer;
    return HEOS._rhomolar*R_u*HEOS._T*(1+HEOS._delta.pt()*HEOS.dalphar_dDelta()*(2-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._delta.pt()*nd2alphar_dni_dDelta);
}

long double MixtureDerivatives::ndalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    double term1 = HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i));
    double term2 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions,i);

    double s = 0;
    for (unsigned int k = 0; k < HEOS.mole_fractions.size(); k++)
    {
        s += HEOS.mole_fractions[k]*dalphar_dxi(HEOS, k);
    }
    double term3 = dalphar_dxi(HEOS, i);
    return term1 + term2 + term3 - s;
}
long double MixtureDerivatives::ndln_fugacity_coefficient_dnj__constT_p(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    long double R_u = HEOS.gas_constant();
    return nd2nalphardnidnj__constT_V(HEOS, j, i) + 1 - partial_molar_volume(HEOS, j)/(R_u*HEOS._T)*ndpdni__constT_V_nj(HEOS, i);
}
long double MixtureDerivatives::nddeltadni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS._delta.pt()-HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i);
}
long double MixtureDerivatives::ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    return HEOS._tau.pt()/HEOS._reducing.T*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i);
}
long double MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    double line1 = HEOS._delta.pt()*d2alphar_dxi_dDelta(HEOS, j)*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i));
    double line2 = -HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1/HEOS._reducing.rhomolar)*(HEOS.Reducing.p->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j)-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->drhormolardxi__constxj(HEOS.mole_fractions,j)*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i));
    double line3 = HEOS._tau.pt()*d2alphar_dxi_dTau(HEOS, j)*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i);
    double line4 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*(HEOS.Reducing.p->d_ndTrdni_dxj__constxi(HEOS.mole_fractions,i,j)-1/HEOS._reducing.T*HEOS.Reducing.p->dTrdxi__constxj(HEOS.mole_fractions, j)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i));
    double s = 0;
    for (unsigned int m = 0; m < HEOS.mole_fractions.size(); m++)
    {
        s += HEOS.mole_fractions[m]*d2alphardxidxj(HEOS, j,m);
    }
    double line5 = d2alphardxidxj(HEOS, i,j)-dalphar_dxi(HEOS, j)-s;
    return line1+line2+line3+line4+line5;
}
long double MixtureDerivatives::nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j)
{
    double line0 = ndalphar_dni__constT_V_nj(HEOS, j); // First term from 7.46
    double line1 = d_ndalphardni_dDelta(HEOS, i)*nddeltadni__constT_V_nj(HEOS, j);
    double line2 = d_ndalphardni_dTau(HEOS, i)*ndtaudni__constT_V_nj(HEOS, j);
    double summer = 0;
    for (unsigned int k = 0; k < HEOS.mole_fractions.size(); k++)
    {
        summer += HEOS.mole_fractions[k]*d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k);
    }
    double line3 = d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j)-summer;
    return line0 + line1 + line2 + line3;
}
long double MixtureDerivatives::d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    // The first line
    double term1 = (HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta())*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i));

    // The second line
    double term2 = HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i);

    // The third line
    double term3 = d2alphar_dxi_dDelta(HEOS, i);
    for (unsigned int k = 0; k < HEOS.mole_fractions.size(); k++)
    {
        term3 -= HEOS.mole_fractions[k]*d2alphar_dxi_dDelta(HEOS, k);
    }
    return term1 + term2 + term3;
}

long double MixtureDerivatives::d_ndalphardni_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i)
{
    // The first line
    double term1 = HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i));

    // The second line
    double term2 = (HEOS._tau.pt()*HEOS.d2alphar_dTau2()+HEOS.dalphar_dTau())*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i);

    // The third line
    double term3 = d2alphar_dxi_dTau(HEOS, i);
    for (unsigned int k = 0; k < HEOS.mole_fractions.size(); k++)
    {
        term3 -= HEOS.mole_fractions[k]*d2alphar_dxi_dTau(HEOS, k);
    }
    return term1 + term2 + term3;
}


} /* namespace CoolProp */
