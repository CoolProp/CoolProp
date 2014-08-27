#include "MixtureDerivatives.h"

namespace CoolProp{
    
long double MixtureDerivatives::dalphar_dxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->baser(HEOS._tau, HEOS._delta) + HEOS.Excess.dalphar_dxi(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<long double> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return _HUGE;
        double dar_dxi = HEOS.components[i]->pEOS->baser(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->baser(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->alphar(HEOS._tau, HEOS._delta);
        dar_dxi += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i==k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->alphar(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->alphar(HEOS._tau, HEOS._delta);
            dar_dxi += x[i]*(Fikarik - FiNariN - FkNarkN);
        }
        return dar_dxi;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}
long double MixtureDerivatives::d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dTau(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<long double> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return _HUGE;
        double d2ar_dxi_dTau = HEOS.components[i]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->dalphar_dTau(HEOS._tau, HEOS._delta);
        d2ar_dxi_dTau += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i==k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->dalphar_dTau(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->dalphar_dTau(HEOS._tau, HEOS._delta);
            d2ar_dxi_dTau += x[i]*(Fikarik - FiNariN - FkNarkN);
        }
        return d2ar_dxi_dTau;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
        
}
long double MixtureDerivatives::d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dDelta(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
        }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<long double> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return _HUGE;
        double d2ar_dxi_dDelta = HEOS.components[i]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->dalphar_dDelta(HEOS._tau, HEOS._delta);
        d2ar_dxi_dDelta += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i==k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->dalphar_dDelta(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->dalphar_dDelta(HEOS._tau, HEOS._delta);
            d2ar_dxi_dDelta += x[i]*(Fikarik - FiNariN - FkNarkN);
        }
        return d2ar_dxi_dDelta;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}

long double MixtureDerivatives::d2alphardxidxj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
            return 0                           + HEOS.Excess.d2alphardxidxj(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i, j);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::size_t N = HEOS.mole_fractions.size();
        if (i == N-1){ return 0;}
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->alphar(HEOS._tau, HEOS._delta);
        if (i == j) { return FiNariN; }
        if (j == N-1){ return _HUGE; }
        
        double Fijarij = HEOS.Excess.F[i][j]*HEOS.Excess.DepartureFunctionMatrix[i][j]->alphar(HEOS._tau, HEOS._delta);
        double FjNarjN = HEOS.Excess.F[j][N-1]*HEOS.Excess.DepartureFunctionMatrix[j][N-1]->alphar(HEOS._tau, HEOS._delta);
        return Fijarij-FiNariN-FjNarjN;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}

long double MixtureDerivatives::fugacity_i(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag){
    return HEOS.mole_fractions[i]*HEOS.rhomolar()*HEOS.gas_constant()*HEOS.T()*exp( dnalphar_dni__constT_V_nj(HEOS, i, xN_flag));
}
long double MixtureDerivatives::ln_fugacity_coefficient(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag)-log(1+HEOS._delta.pt()*HEOS.dalphar_dDelta());
}
long double MixtureDerivatives::dln_fugacity_coefficient_dT__constrho_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double dtau_dT = -HEOS._tau.pt()/HEOS._T; //[1/K]
    return (HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag)-1/(1+HEOS._delta.pt()*HEOS.dalphar_dDelta())*(HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()))*dtau_dT;
}
long double MixtureDerivatives::dln_fugacity_coefficient_drho__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double ddelta_drho = 1/HEOS._reducing.rhomolar; //[m^3/mol]
    return (HEOS.dalphar_dDelta() + d_ndalphardni_dDelta(HEOS, i, xN_flag)-1/(1+HEOS._delta.pt()*HEOS.dalphar_dDelta())*(HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta()))*ddelta_drho;
}
long double MixtureDerivatives::dnalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // GERG Equation 7.42
    return HEOS.alphar() + ndalphar_dni__constT_V_nj(HEOS, i, xN_flag);
}
long double MixtureDerivatives::d2nalphar_dni_dT(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -HEOS._tau.pt()/HEOS._T*(HEOS.dalphar_dTau() + d_ndalphardni_dTau(HEOS, i, xN_flag));
}
long double MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double T = HEOS._reducing.T/HEOS._tau.pt();
    long double R_u = HEOS.gas_constant();
    return d2nalphar_dni_dT(HEOS, i, xN_flag) + 1/T-partial_molar_volume(HEOS, i, xN_flag)/(R_u*T)*dpdT__constV_n(HEOS);
}
long double MixtureDerivatives::partial_molar_volume(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return -ndpdni__constT_V_nj(HEOS, i, xN_flag)/ndpdV__constT_n(HEOS);
}

long double MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // GERG equation 7.30
    long double R_u = HEOS.gas_constant();
    double partial_molar_volumeval = partial_molar_volume(HEOS, i, xN_flag); // [m^3/mol]
    double term1 = partial_molar_volumeval/(R_u*HEOS._T); // m^3/mol/(N*m)*mol = m^2/N = 1/Pa
    double term2 = 1.0/HEOS.p();
    return term1 - term2;
}

long double MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.115
    long double R_u = HEOS.gas_constant();
    // partial molar volume is -dpdn/dpdV, so need to flip the sign here
    return d2nalphar_dxj_dni__constT_V(HEOS, j, i, xN_flag) - partial_molar_volume(HEOS, i, xN_flag)/(R_u*HEOS._T)*dpdxj__constT_V_xi(HEOS, j, xN_flag);
}
long double MixtureDerivatives::dpdxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.130
    long double R_u = HEOS.gas_constant();
    return HEOS._rhomolar*R_u*HEOS._T*(ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)*HEOS.dalphar_dDelta()+HEOS._delta.pt()*d_dalpharddelta_dxj__constT_V_xi(HEOS, j, xN_flag));
}

long double MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert Equation 3.134 (Catch test provided)
    return HEOS.d2alphar_dDelta2()*ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)
         + HEOS.d2alphar_dDelta_dTau()*dtau_dxj__constT_V_xi(HEOS, j, xN_flag)
         + d2alphar_dxi_dDelta(HEOS, j, xN_flag);
}

long double MixtureDerivatives::dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    //Gernert 3.119 (Catch test provided)
    return HEOS.dalphar_dDelta()*ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)+HEOS.dalphar_dTau()*dtau_dxj__constT_V_xi(HEOS, j, xN_flag)+dalphar_dxi(HEOS, j, xN_flag);
}
long double MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.118
    return d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i,j, xN_flag)
          + ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)*d_ndalphardni_dDelta(HEOS, i, xN_flag)
          + dtau_dxj__constT_V_xi(HEOS, j, xN_flag)*d_ndalphardni_dTau(HEOS, i, xN_flag);
}
long double MixtureDerivatives::ddelta_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.121 (Catch test provided)
    return -HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing.p->drhormolardxi__constxj(HEOS.mole_fractions,j, xN_flag);
}
long double MixtureDerivatives::dtau_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    // Gernert 3.122 (Catch test provided)
    return 1/HEOS._T*HEOS.Reducing.p->dTrdxi__constxj(HEOS.mole_fractions,j, xN_flag);
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
long double MixtureDerivatives::ndpdni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // Eqn 7.64 and 7.63
    long double R_u = HEOS.gas_constant();
    double ndrhorbar_dni__constnj = HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag);
    double ndTr_dni__constnj = HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions,i, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; ++k)
    {
        summer += HEOS.mole_fractions[k]*d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    double nd2alphar_dni_dDelta = HEOS._delta.pt()*HEOS.d2alphar_dDelta2()*(1-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()/HEOS._reducing.T*ndTr_dni__constnj+d2alphar_dxi_dDelta(HEOS, i, xN_flag)-summer;
    return HEOS._rhomolar*R_u*HEOS._T*(1+HEOS._delta.pt()*HEOS.dalphar_dDelta()*(2-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._delta.pt()*nd2alphar_dni_dDelta);
}

long double MixtureDerivatives::ndalphar_dni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    double term1 = HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag));
    double term2 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions,i, xN_flag);

    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        s += HEOS.mole_fractions[k]*dalphar_dxi(HEOS, k, xN_flag);
    }
    double term3 = dalphar_dxi(HEOS, i, xN_flag);
    return term1 + term2 + term3 - s;
}
long double MixtureDerivatives::ndln_fugacity_coefficient_dnj__constT_p(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    long double R_u = HEOS.gas_constant();
    return nd2nalphardnidnj__constT_V(HEOS, j, i, xN_flag) + 1 - partial_molar_volume(HEOS, j, xN_flag)/(R_u*HEOS._T)*ndpdni__constT_V_nj(HEOS, i, xN_flag);
}
long double MixtureDerivatives::nddeltadni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS._delta.pt()-HEOS._delta.pt()/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag);
}
long double MixtureDerivatives::ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS._tau.pt()/HEOS._reducing.T*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}
long double MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = HEOS._delta.pt()*d2alphar_dxi_dDelta(HEOS, j, xN_flag)*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double line2 = -HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1/HEOS._reducing.rhomolar)*(HEOS.Reducing.p->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->drhormolardxi__constxj(HEOS.mole_fractions,j, xN_flag)*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag));
    double line3 = HEOS._tau.pt()*d2alphar_dxi_dTau(HEOS, j, xN_flag)*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
    double line4 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*(HEOS.Reducing.p->d_ndTrdni_dxj__constxi(HEOS.mole_fractions,i,j, xN_flag)-1/HEOS._reducing.T*HEOS.Reducing.p->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));
    double s = 0;
    std::size_t mmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ mmax--; }
    for (unsigned int m = 0; m < mmax; m++)
    {
        s += HEOS.mole_fractions[m]*d2alphardxidxj(HEOS, j,m, xN_flag);
    }
    double line5 = d2alphardxidxj(HEOS, i,j, xN_flag)-dalphar_dxi(HEOS, j, xN_flag)-s;
    return line1+line2+line3+line4+line5;
}
long double MixtureDerivatives::nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line0 = ndalphar_dni__constT_V_nj(HEOS, j, xN_flag); // First term from 7.46
    double line1 = d_ndalphardni_dDelta(HEOS, i, xN_flag)*nddeltadni__constT_V_nj(HEOS, j, xN_flag);
    double line2 = d_ndalphardni_dTau(HEOS, i, xN_flag)*ndtaudni__constT_V_nj(HEOS, j, xN_flag);
    double summer = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        summer += HEOS.mole_fractions[k]*d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, k, xN_flag);
    }
    double line3 = d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag)-summer;
    return line0 + line1 + line2 + line3;
}
long double MixtureDerivatives::d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = (HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta())*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = d2alphar_dxi_dDelta(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}
long double MixtureDerivatives::d_ndalphardni_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing.p->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = (HEOS._tau.pt()*HEOS.d2alphar_dTau2()+HEOS.dalphar_dTau())*(1/HEOS._reducing.T)*HEOS.Reducing.p->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

    // The third line
    double term3 = d2alphar_dxi_dTau(HEOS, i, xN_flag);
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        term3 -= HEOS.mole_fractions[k]*d2alphar_dxi_dTau(HEOS, k, xN_flag);
    }
    return term1 + term2 + term3;
}

} /* namespace CoolProp */

#ifdef ENABLE_CATCH
#include "catch.hpp"

using namespace CoolProp;
TEST_CASE("Mixture derivative checks", "[mixtures],[mixture_derivs]")
{
    /* Set up a test class for the mixture tests */
    std::vector<std::string> names(2);
    names[0] = "Ethane"; names[1] = "Propane";
    std::vector<long double> z(2);
    z[0] = 0.25; z[1] = 1-z[0];
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS(new HelmholtzEOSMixtureBackend(names));
    HelmholtzEOSMixtureBackend &rHEOS = *(HEOS.get());
    rHEOS.set_mole_fractions(z);
    
    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    // These ones only require the i index
    for (std::size_t i = 0; i< z.size();++i)
    {
        std::ostringstream ss1;
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
            CHECK(err < 1e-8);
        }
        
        std::ostringstream ss2;
        ss2 << "dln_fugacity_coefficient_dp__constT_n, i=" << i;
        SECTION(ss2.str(), "")
        {
            double p0 = 101325, dP = 1e-4;
            rHEOS.specify_phase(iphase_gas);
            
            rHEOS.update(PT_INPUTS, p0, 300);
            double analytic = MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rHEOS, i, xN_flag);
            double rho1 = rHEOS.rhomolar();
            
            rHEOS.update(DmolarT_INPUTS, rho1 + dP, 300);
            double v1 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag), p1 = rHEOS.p();
            rHEOS.update(DmolarT_INPUTS, rho1 - dP, 300);
            double v2 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag), p2 = rHEOS.p();
            
            double numeric = (v1 - v2)/(p1 - p2);
            double err = std::abs((numeric-analytic)/analytic);
            CHECK(err < 1e-8);
        }
        std::ostringstream ss3;
        ss3 << "d_ndalphardni_dDelta, i=" << i;
        SECTION(ss3.str(), "")
        {
            double p1 = 101325, dP = 1e-1;
            rHEOS.specify_phase(iphase_gas);
            
            rHEOS.update(PT_INPUTS, p1, 300);
            double analytic = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS, i, xN_flag);
            
            rHEOS.update(PT_INPUTS, p1 + dP, 300);
            double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag), delta1 = rHEOS.delta();
            rHEOS.update(PT_INPUTS, p1 - dP, 300);
            double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag), delta2 = rHEOS.delta();
            
            double numeric = (v1 - v2)/(delta1 - delta2);
            double err = std::abs((numeric-analytic)/analytic);
            CHECK(err < 1e-8);
        }
        std::ostringstream ss4;
        ss4 << "d_ndalphardni_dTau, i=" << i;
        SECTION(ss4.str(), "")
        {
            double p1 = 101325, dT = 1e-2;
            rHEOS.specify_phase(iphase_gas);
            rHEOS.update(PT_INPUTS, 101325, 300);
            double rho1 = rHEOS.rhomolar();
            
            rHEOS.update(DmolarT_INPUTS, rho1, 300);
            double analytic = MixtureDerivatives::d_ndalphardni_dTau(rHEOS, i, xN_flag);
            
            rHEOS.update(DmolarT_INPUTS, rho1, 300 + dT);
            double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag), tau1 = rHEOS.tau();
            rHEOS.update(DmolarT_INPUTS, rho1, 300 - dT);
            double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag), tau2 = rHEOS.tau();
            
            double numeric = (v1 - v2)/(tau1 - tau2);
            double err = std::abs((numeric-analytic)/analytic);
            CHECK(err < 1e-8);
        }
        std::ostringstream ss5;
        ss5 << "dpdxj__constT_V_xi, i=" << i;
        SECTION(ss5.str(), "")
        {
            double dz = 1e-6;
            rHEOS.specify_phase(iphase_gas);
            rHEOS.update(DmolarT_INPUTS, 300, 300);
            
            double rho1 = rHEOS.rhomolar();
            double analytic = MixtureDerivatives::dpdxj__constT_V_xi(rHEOS, i, xN_flag);
            std::vector<long double> zp = z; /// Copy base composition
            zp[i] += dz;
            rHEOS.set_mole_fractions(zp);
            rHEOS.update(DmolarT_INPUTS, rho1, 300);
            double v1 = rHEOS.p();
            std::vector<long double> zm = z; /// Copy base composition
            zm[i] -= dz;
            rHEOS.set_mole_fractions(zm);
            rHEOS.update(DmolarT_INPUTS, rho1, 300);
            double v2 = rHEOS.p();
            
            double numeric = (v1 - v2)/(2*dz);
            double err = std::abs((numeric-analytic)/analytic);
            
            CAPTURE(numeric);
            CAPTURE(analytic);
            CHECK(err < 1e-80);
        }
        
        std::ostringstream ss6;
        ss6 << "d_dalpharddelta_dxj__constT_V_xi, i=" << i;
        SECTION(ss6.str(), "")
        {
            double dz = 1e-6;
            rHEOS.specify_phase(iphase_gas);
            rHEOS.update(DmolarT_INPUTS, 300, 300);
            
            double rho1 = rHEOS.rhomolar();
            double analytic = MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(rHEOS, i, xN_flag);
            
            // Increment mole fraction
            std::vector<long double> zp = z; /// Copy base composition
            zp[i] += dz;
            rHEOS.set_mole_fractions(zp);
            rHEOS.update(DmolarT_INPUTS, rho1, 300);
            double v1 = rHEOS.dalphar_dDelta();
            
            // Decrement mole fraction
            std::vector<long double> zm = z; /// Copy base composition
            zm[i] -= dz;
            rHEOS.set_mole_fractions(zm);
            rHEOS.update(DmolarT_INPUTS, rho1, 300);
            double v2 = rHEOS.dalphar_dDelta();
            
            double numeric = (v1 - v2)/(2*dz);
            double err = std::abs((numeric-analytic)/analytic);
            
            CAPTURE(numeric);
            CAPTURE(analytic);
            CHECK(err < 1e-60);
        }
            
        // These derivatives depend on both the i and j indices
        for (std::size_t j = 0; j < z.size(); ++j){
            std::ostringstream ss1;
            ss1 << "dln_fugacity_coefficient_dxj__constT_p_xi, i=" << i << ", j=" << j;
            SECTION(ss1.str(), "")
            {
                double dz = 1e-6;
                rHEOS.specify_phase(iphase_gas);
                rHEOS.update(DmolarT_INPUTS, 300, 300);
                double p = rHEOS.p();
                CAPTURE(p);
                
                double rho1 = rHEOS.rhomolar();
                double analytic = MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(rHEOS, i, j, xN_flag);
                std::vector<long double> zp = z; /// Copy base composition
                zp[j] += dz;
                rHEOS.set_mole_fractions(zp);
                rHEOS.update(PT_INPUTS, p, 300);
                double v1 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag);
                std::vector<long double> zm = z; /// Copy base composition
                zm[j] -= dz;
                rHEOS.set_mole_fractions(zm);
                rHEOS.update(PT_INPUTS, p, 300);
                double v2 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS, i, xN_flag);
                
                double numeric = (v1 - v2)/(2*dz);
                double err = std::abs((numeric-analytic)/analytic);
                
                CAPTURE(numeric);
                CAPTURE(analytic);
                CHECK(err < 1e-8);
            }
            std::ostringstream ss3;
            ss2 << "d_ndalphardni_dxj__constT_V_xi, i=" << i << ", j=" << j;
            SECTION(ss2.str(), "")
            {
                double dz = 1e-6;
                rHEOS.specify_phase(iphase_gas);
                rHEOS.update(DmolarT_INPUTS, 300, 300);
                
                double rho1 = rHEOS.rhomolar();
                double analytic = MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(rHEOS, i, j, xN_flag);
                
                // Increment mole fraction
                std::vector<long double> zp = z; /// Copy base composition
                zp[j] += dz;
                rHEOS.set_mole_fractions(zp);
                rHEOS.update(DmolarT_INPUTS, rho1, 300);
                double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag);
                
                // Decrement mole fraction
                std::vector<long double> zm = z; /// Copy base composition
                zm[j] -= dz;
                rHEOS.set_mole_fractions(zm);
                rHEOS.update(DmolarT_INPUTS, rho1, 300);
                double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS, i, xN_flag);
                
                double numeric = (v1 - v2)/(2*dz);
                double err = std::abs((numeric-analytic)/analytic);
                
                CAPTURE(numeric);
                CAPTURE(analytic);
                CHECK(err < 1e-60);
            }
            
        }
        
    }
}
#endif


