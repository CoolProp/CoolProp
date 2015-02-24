#include "MixtureDerivatives.h"

namespace CoolProp{
    
CoolPropDbl MixtureDerivatives::dalphar_dxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->baser(HEOS._tau, HEOS._delta) + HEOS.Excess.dalphar_dxi(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i == N-1) return 0;
        double dar_dxi = HEOS.components[i]->pEOS->baser(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->baser(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->alphar(HEOS._tau, HEOS._delta);
        dar_dxi += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i == k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->alphar(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->alphar(HEOS._tau, HEOS._delta);
            dar_dxi += x[k]*(Fikarik - FiNariN - FkNarkN);
        }
        return dar_dxi;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}
CoolPropDbl MixtureDerivatives::d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dTau(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return 0;
        double d2ar_dxi_dTau = HEOS.components[i]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->dalphar_dTau(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->dalphar_dTau(HEOS._tau, HEOS._delta);
        d2ar_dxi_dTau += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i==k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->dalphar_dTau(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->dalphar_dTau(HEOS._tau, HEOS._delta);
            d2ar_dxi_dTau += x[k]*(Fikarik - FiNariN - FkNarkN);
        }
        return d2ar_dxi_dTau;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
        
}
CoolPropDbl MixtureDerivatives::d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dDelta(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
        }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return 0;
        double d2ar_dxi_dDelta = HEOS.components[i]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta) - HEOS.components[N-1]->pEOS->dalphar_dDelta(HEOS._tau, HEOS._delta);
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->dalphar_dDelta(HEOS._tau, HEOS._delta);
        d2ar_dxi_dDelta += (1-2*x[i])*FiNariN;
        for (std::size_t k = 0; k < N-1; ++k){
            if (i==k) continue;
            double Fikarik = HEOS.Excess.F[i][k]*HEOS.Excess.DepartureFunctionMatrix[i][k]->dalphar_dDelta(HEOS._tau, HEOS._delta);
            double FkNarkN = HEOS.Excess.F[k][N-1]*HEOS.Excess.DepartureFunctionMatrix[k][N-1]->dalphar_dDelta(HEOS._tau, HEOS._delta);
            d2ar_dxi_dDelta += x[k]*(Fikarik - FiNariN - FkNarkN);
        }
        return d2ar_dxi_dDelta;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}

CoolPropDbl MixtureDerivatives::d2alphardxidxj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
            return 0                           + HEOS.Excess.d2alphardxidxj(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i, j);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::size_t N = HEOS.mole_fractions.size();
        if (i == N-1 || j == N-1){ return 0;}
        double FiNariN = HEOS.Excess.F[i][N-1]*HEOS.Excess.DepartureFunctionMatrix[i][N-1]->alphar(HEOS._tau, HEOS._delta);
        if (i == j) { return -2*FiNariN; }
        
        double Fijarij = HEOS.Excess.F[i][j]*HEOS.Excess.DepartureFunctionMatrix[i][j]->alphar(HEOS._tau, HEOS._delta);
        double FjNarjN = HEOS.Excess.F[j][N-1]*HEOS.Excess.DepartureFunctionMatrix[j][N-1]->alphar(HEOS._tau, HEOS._delta);
        return Fijarij - FiNariN - FjNarjN;
    }
    else{
        throw ValueError(format("xN_flag is invalid"));
    }
}

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
    CoolPropDbl line4 = dalphar_dxi(HEOS, j, xN_flag) + d_ndalphardni_dxj__constdelta_tau_xi(HEOS, i, j, xN_flag);
    
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
         + d2alphar_dxi_dDelta(HEOS, j, xN_flag);
}

CoolPropDbl MixtureDerivatives::dalphar_dxj__constT_V_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t j, x_N_dependency_flag xN_flag)
{
    //Gernert 3.119 (Catch test provided)
    return HEOS.dalphar_dDelta()*ddelta_dxj__constT_V_xi(HEOS, j, xN_flag)+HEOS.dalphar_dTau()*dtau_dxj__constT_V_xi(HEOS, j, xN_flag)+dalphar_dxi(HEOS, j, xN_flag);
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
        summer += HEOS.mole_fractions[k]*d2alphar_dxi_dDelta(HEOS, k, xN_flag);
    }
    double nd2alphar_dni_dDelta = HEOS._delta.pt()*HEOS.d2alphar_dDelta2()*(1-1/HEOS._reducing.rhomolar*ndrhorbar_dni__constnj)+HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()/HEOS._reducing.T*ndTr_dni__constnj+d2alphar_dxi_dDelta(HEOS, i, xN_flag)-summer;
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
        s += HEOS.mole_fractions[k]*dalphar_dxi(HEOS, k, xN_flag);
    }
    double term3 = dalphar_dxi(HEOS, i, xN_flag);
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
CoolPropDbl MixtureDerivatives::ndtaudni__constT_V_nj(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    return HEOS._tau.pt()/HEOS._reducing.T*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
}
CoolPropDbl MixtureDerivatives::d_ndalphardni_dxj__constdelta_tau_xi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
{
    double line1 = HEOS._delta.pt()*d2alphar_dxi_dDelta(HEOS, j, xN_flag)*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));
    double line3 = HEOS._tau.pt()*d2alphar_dxi_dTau(HEOS, j, xN_flag)*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);
    double line2 = -HEOS._delta.pt()*HEOS.dalphar_dDelta()*(1/HEOS._reducing.rhomolar)*(HEOS.Reducing->d_ndrhorbardni_dxj__constxi(HEOS.mole_fractions, i, j, xN_flag)-1/HEOS._reducing.rhomolar*HEOS.Reducing->drhormolardxi__constxj(HEOS.mole_fractions,j, xN_flag)*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions,i, xN_flag));
    double line4 = HEOS._tau.pt()*HEOS.dalphar_dTau()*(1/HEOS._reducing.T)*(HEOS.Reducing->d_ndTrdni_dxj__constxi(HEOS.mole_fractions,i,j, xN_flag)-1/HEOS._reducing.T*HEOS.Reducing->dTrdxi__constxj(HEOS.mole_fractions, j, xN_flag)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag));
    
    double s = 0;
    std::size_t kmax = HEOS.mole_fractions.size();
    if (xN_flag == XN_DEPENDENT){ kmax--; }
    for (unsigned int k = 0; k < kmax; k++)
    {
        s += HEOS.mole_fractions[k]*d2alphardxidxj(HEOS, j, k, xN_flag);
    }
    double line5 = d2alphardxidxj(HEOS, i, j, xN_flag)-dalphar_dxi(HEOS, j, xN_flag)-s;
    return line1 + line2 + line3 + line4 + line5;
}
CoolPropDbl MixtureDerivatives::nd2nalphardnidnj__constT_V(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag)
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
CoolPropDbl MixtureDerivatives::d_ndalphardni_dDelta(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = (HEOS._delta.pt()*HEOS.d2alphar_dDelta2()+HEOS.dalphar_dDelta())*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = HEOS._tau.pt()*HEOS.d2alphar_dDelta_dTau()*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

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
CoolPropDbl MixtureDerivatives::d_ndalphardni_dTau(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    // The first line
    double term1 = HEOS._delta.pt()*HEOS.d2alphar_dDelta_dTau()*(1-1/HEOS._reducing.rhomolar*HEOS.Reducing->ndrhorbardni__constnj(HEOS.mole_fractions, i, xN_flag));

    // The second line
    double term2 = (HEOS._tau.pt()*HEOS.d2alphar_dTau2()+HEOS.dalphar_dTau())*(1/HEOS._reducing.T)*HEOS.Reducing->ndTrdni__constnj(HEOS.mole_fractions, i, xN_flag);

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
    for (std::size_t Ncomp = 2; Ncomp <= 4; Ncomp++)
    {
        std::ostringstream ss00;
        ss00 << "Mixture with " << Ncomp << " components";
        SECTION(ss00.str(),"")
        {
            x_N_dependency_flag xN_flag;
            for (std::size_t xNxN = 0; xNxN <=1; xNxN++){
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
                    std::vector<std::string> names;
                    std::vector<CoolPropDbl> z;
                    shared_ptr<HelmholtzEOSMixtureBackend> HEOS, HEOS_plusT_constrho, HEOS_plusrho_constT, HEOS_minusT_constrho, HEOS_minusrho_constT;
                    names.resize(Ncomp);
                    z.resize(Ncomp);
                    
                    /* Set up a test class for the mixture tests */
                    if (Ncomp == 2)
                    {
                        names[0] = "Ethane"; names[1] = "Propane";
                        z[0] = 0.3; z[1] = 0.7;
                    }
                    else if (Ncomp == 3)
                    {
                        names[0] = "Ethane"; names[1] = "Propane"; names[2] = "Methane";
                        z[0] = 0.3; z[1] = 0.4; z[2] = 0.3;
                    }
                    else if (Ncomp == 4)
                    {
                        names[0] = "Ethane"; names[1] = "Propane"; names[2] = "Methane"; names[3] = "n-Butane";
                        z[0] = 0.3; z[1] = 0.4; z[2] = 0.2; z[3] = 0.1;
                    }
                    double T1 = 300, rho1 = 300, dT = 1e-3, drho = 1e-3, dz = 1e-6;
                    
                    HEOS.reset(new HelmholtzEOSMixtureBackend(names));
                    HEOS->specify_phase(iphase_gas);
                    HelmholtzEOSMixtureBackend &rHEOS = *(HEOS.get());
                    rHEOS.set_mole_fractions(z);
                    rHEOS.update(DmolarT_INPUTS, rho1, T1);
                    
                    HEOS_plusT_constrho.reset(new HelmholtzEOSMixtureBackend(names));
                    HEOS_plusT_constrho->specify_phase(iphase_gas); 
                    HelmholtzEOSMixtureBackend &rHEOS_plusT_constrho = *(HEOS_plusT_constrho.get());
                    rHEOS_plusT_constrho.set_mole_fractions(z);
                    rHEOS_plusT_constrho.update(DmolarT_INPUTS, rho1, T1 + dT);
                    
                    HEOS_minusT_constrho.reset(new HelmholtzEOSMixtureBackend(names));
                    HEOS_minusT_constrho->specify_phase(iphase_gas); 
                    HelmholtzEOSMixtureBackend &rHEOS_minusT_constrho = *(HEOS_minusT_constrho.get());
                    rHEOS_minusT_constrho.set_mole_fractions(z);
                    rHEOS_minusT_constrho.update(DmolarT_INPUTS, rho1, T1 - dT);
                    
                    HEOS_plusrho_constT.reset(new HelmholtzEOSMixtureBackend(names));
                    HEOS_plusrho_constT->specify_phase(iphase_gas);            
                    HelmholtzEOSMixtureBackend &rHEOS_plusrho_constT = *(HEOS_plusrho_constT.get());
                    rHEOS_plusrho_constT.set_mole_fractions(z);
                    rHEOS_plusrho_constT.update(DmolarT_INPUTS, rho1 + drho, T1);
                    
                    HEOS_minusrho_constT.reset(new HelmholtzEOSMixtureBackend(names));
                    HEOS_minusrho_constT->specify_phase(iphase_gas);
                    HelmholtzEOSMixtureBackend &rHEOS_minusrho_constT = *(HEOS_minusrho_constT.get());
                    rHEOS_minusrho_constT.set_mole_fractions(z);
                    rHEOS_minusrho_constT.update(DmolarT_INPUTS, rho1 - drho, T1);
                    
                    rHEOS.update(DmolarT_INPUTS, rho1, T1);
                
                    // These ones only require the i index
                    for (std::size_t i = 0; i< z.size();++i)
                    {
                        shared_ptr<HelmholtzEOSMixtureBackend> HEOS, HEOS_pluszi, HEOS_minuszi;
                        HEOS_pluszi.reset(new HelmholtzEOSMixtureBackend(names));
                        HEOS_pluszi->specify_phase(iphase_gas); 
                        HelmholtzEOSMixtureBackend &rHEOS_pluszi = *(HEOS_pluszi.get());
                        std::vector<CoolPropDbl> zp = z; /// Copy base composition
                        zp[i] += dz; 
                        if (xN_flag == XN_DEPENDENT)
                            zp[z.size()-1] -= dz;
                        rHEOS_pluszi.set_mole_fractions(zp);
                        rHEOS_pluszi.update(DmolarT_INPUTS, rho1, T1);
                    
                        HEOS_minuszi.reset(new HelmholtzEOSMixtureBackend(names));
                        HEOS_minuszi->specify_phase(iphase_gas); 
                        HelmholtzEOSMixtureBackend &rHEOS_minuszi = *(HEOS_minuszi.get());
                        std::vector<CoolPropDbl> zm = z; /// Copy base composition
                        zm[i] -= dz; 
                        if (xN_flag == XN_DEPENDENT)
                            zm[z.size()-1] += dz;
                        rHEOS_minuszi.set_mole_fractions(zm);
                        rHEOS_minuszi.update(DmolarT_INPUTS, rho1, T1);
                        
                        std::ostringstream ss0;
                        ss0 << "dln_fugacity_i_dT__constrho_n, i=" << i;
                        SECTION(ss0.str(),"")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_i_dT__constrho_n(rHEOS, i, xN_flag);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_plusT_constrho, i, xN_flag));
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minusT_constrho, i, xN_flag));
                            double numeric = (v1 - v2)/(2*dT);
                            double err = std::abs((numeric-analytic)/analytic);
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
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss0b;
                        ss0b << "dln_fugacity_i_drho__constT_n, i=" << i;
                        SECTION(ss0b.str(), "")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_i_drho__constT_n(rHEOS, i, xN_flag);
                            double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_plusrho_constT, i, xN_flag));
                            double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minusrho_constT, i, xN_flag));
                            double numeric = (v1 - v2)/(2*dT);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss2;
                        ss2 << "dln_fugacity_coefficient_dp__constT_n, i=" << i;
                        SECTION(ss2.str(), "")
                        {
                            double analytic = MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS_plusrho_constT, i, xN_flag), p1 = rHEOS_plusrho_constT.p();
                            double v2 = MixtureDerivatives::ln_fugacity_coefficient(rHEOS_minusrho_constT, i, xN_flag), p2 = rHEOS_minusrho_constT.p();
                            double numeric = (v1 - v2)/(p1 - p2);
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-8);
                        }
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
                        
                        std::ostringstream ss3;
                        ss3 << "d_ndalphardni_dDelta, i=" << i;
                        SECTION(ss3.str(), "")
                        {
                            double analytic = MixtureDerivatives::d_ndalphardni_dDelta(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_plusrho_constT, i, xN_flag),  delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3a;
                        ss3a << "d2alphar_dxi_dDelta, i=" << i;
                        SECTION(ss3a.str(), "")
                        {
                            if (i == z.size()-1){break;}
                            double analytic = MixtureDerivatives::d2alphar_dxi_dDelta(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::dalphar_dxi(rHEOS_plusrho_constT, i, xN_flag), delta1 = rHEOS_plusrho_constT.delta();
                            double v2 = MixtureDerivatives::dalphar_dxi(rHEOS_minusrho_constT, i, xN_flag), delta2 = rHEOS_minusrho_constT.delta();
                            double numeric = (v1 - v2)/(delta1 - delta2);
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss4a;
                        ss4a << "d2alphar_dxi_dTau, i=" << i;
                        SECTION(ss4a.str(), "")
                        {
                            if (i == z.size()-1){ break; }
                            double analytic = MixtureDerivatives::d2alphar_dxi_dTau(rHEOS, i, xN_flag);
                            double v1 = MixtureDerivatives::dalphar_dxi(rHEOS_plusT_constrho, i, xN_flag), tau1 = rHEOS_plusT_constrho.tau();
                            double v2 = MixtureDerivatives::dalphar_dxi(rHEOS_minusT_constrho, i, xN_flag), tau2 = rHEOS_minusT_constrho.tau();
                            double numeric = (v1 - v2)/(tau1 - tau2);
                            double err = std::abs((numeric-analytic)/analytic);
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
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss5;
                        ss5 << "dpdxj__constT_V_xi, i=" << i;
                        SECTION(ss5.str(), "")
                        {
                            if (i==z.size()-1){break;}
                            double analytic = MixtureDerivatives::dpdxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.p();
                            double v2 = rHEOS_minuszi.p();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }       
                        std::ostringstream ss5a;
                        ss5a << "dtaudxj__constT_V_xi, i=" << i;
                        SECTION(ss5a.str(), "")
                        {
                            if (i==z.size()-1){break;}
                            double analytic = MixtureDerivatives::dtau_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.tau();
                            double v2 = rHEOS_minuszi.tau();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }       
                        std::ostringstream ss5b;
                        ss5b << "ddeltadxj__constT_V_xi, i=" << i;
                        SECTION(ss5b.str(), "")
                        {
                            if (i==z.size()-1){break;}
                            double analytic = MixtureDerivatives::ddelta_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.delta();
                            double v2 = rHEOS_minuszi.delta();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }       
                        std::ostringstream ss6;
                        ss6 << "d_dalpharddelta_dxj__constT_V_xi, i=" << i;
                        SECTION(ss6.str(), "")
                        {
                            if (i==z.size()-1){break;}
                            double analytic = MixtureDerivatives::d_dalpharddelta_dxj__constT_V_xi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi.dalphar_dDelta();
                            double v2 = rHEOS_minuszi.dalphar_dDelta();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }    
                        std::ostringstream ss7;
                        ss7 << "dTrdxi__constxj, i=" << i;
                        SECTION(ss7.str(), "")
                        {
                            if (i == z.size()-1){break;}
                            double analytic = rHEOS.Reducing->dTrdxi__constxj(rHEOS.get_mole_fractions(), i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->Tr(rHEOS_pluszi.get_mole_fractions());
                            double v2 = rHEOS_minuszi.Reducing->Tr(rHEOS_minuszi.get_mole_fractions());
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss8;
                        ss8 << "drhormolardxi__constxj, i=" << i;
                        SECTION(ss8.str(), "")
                        {
                            if (i == z.size()-1){break;}
                            double analytic = rHEOS.Reducing->drhormolardxi__constxj(rHEOS.get_mole_fractions(), i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->rhormolar(rHEOS_pluszi.get_mole_fractions());
                            double v2 = rHEOS_minuszi.Reducing->rhormolar(rHEOS_minuszi.get_mole_fractions());
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3c;
                        ss3c << "d2Trdxi2__constxj, i=" << i;
                        SECTION(ss3c.str(), "")
                        {
                            if (i == z.size()-1){break;}
                            double analytic = rHEOS.Reducing->d2Trdxi2__constxj(z, i, xN_flag);
                            double v1 = rHEOS_pluszi.Reducing->dTrdxi__constxj(rHEOS_pluszi.get_mole_fractions(), i, xN_flag);
                            double v2 = rHEOS_minuszi.Reducing->dTrdxi__constxj(rHEOS_minuszi.get_mole_fractions(), i, xN_flag);
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3d;
                        ss3d << "dalphar_dxi, i=" << i;
                        SECTION(ss3d.str(), "")
                        {
                            if (i == z.size()-1){break;}
                            double analytic = MixtureDerivatives::dalphar_dxi(rHEOS, i, xN_flag);
                            
                            shared_ptr<HelmholtzEOSMixtureBackend> plus(new HelmholtzEOSMixtureBackend(names));
                            plus->specify_phase(iphase_gas);
                            plus->set_mole_fractions(zp);
                            plus->calc_reducing_state();
                            SimpleState red = plus->get_reducing_state();
                            plus->update(DmolarT_INPUTS, red.rhomolar*rHEOS.delta(), red.T/rHEOS.tau());
                            double v1 = plus->alphar();
                            shared_ptr<HelmholtzEOSMixtureBackend> minus(new HelmholtzEOSMixtureBackend(names));
                            minus->specify_phase(iphase_gas);
                            minus->set_mole_fractions(zm);
                            minus->calc_reducing_state();
                            red = minus->get_reducing_state();
                            minus->update(DmolarT_INPUTS, red.rhomolar*rHEOS.delta(), red.T/rHEOS.tau());
                            double v2 = minus->alphar();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
                            CAPTURE(numeric);
                            CAPTURE(analytic);
                            CHECK(err < 1e-7);
                        }
                        
                        // These derivatives depend on both the i and j indices
                        for (std::size_t j = 0; j < z.size(); ++j){
                            if (xN_flag == XN_DEPENDENT && j == z.size()){ continue; }
                            shared_ptr<HelmholtzEOSMixtureBackend> HEOS, HEOS_pluszj, HEOS_minuszj;
                            HEOS_pluszj.reset(new HelmholtzEOSMixtureBackend(names));
                            HEOS_pluszj->specify_phase(iphase_gas); 
                            HelmholtzEOSMixtureBackend &rHEOS_pluszj = *(HEOS_pluszj.get());
                            std::vector<CoolPropDbl> zp = z; /// Copy base composition
                            zp[j] += dz; 
                            if (xN_flag == XN_DEPENDENT)
                                zp[z.size()-1] -= dz;
                            rHEOS_pluszj.set_mole_fractions(zp);
                            rHEOS_pluszj.update(DmolarT_INPUTS, rho1, T1);
                        
                            HEOS_minuszj.reset(new HelmholtzEOSMixtureBackend(names));
                            HEOS_minuszj->specify_phase(iphase_gas); 
                            HelmholtzEOSMixtureBackend &rHEOS_minuszj = *(HEOS_minuszj.get());
                            std::vector<CoolPropDbl> zm = z; /// Copy base composition
                            zm[j] -= dz; 
                            if (xN_flag == XN_DEPENDENT)
                                zm[z.size()-1] += dz;
                            rHEOS_minuszj.set_mole_fractions(zm);
                            rHEOS_minuszj.update(DmolarT_INPUTS, rho1, T1);
                        
                            std::ostringstream ss1a;
                            ss1a << "dln_fugacity_dxj__constT_rho_xi, i=" << i << ", j=" << j;
                            SECTION(ss1a.str(), "")
                            {
                                if (xN_flag == XN_INDEPENDENT){continue;}
                                if (j == z.size()-1){break;}
                                double analytic = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(rHEOS, i, j, xN_flag);
                                double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_pluszj, i, xN_flag));
                                double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minuszj, i, xN_flag));
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss2;
                            ss2 << "d_ndTrdni_dxj, i=" << i << ", j=" << j;
                            SECTION(ss2.str(), "")
                            {
                                if (j == z.size()-1){break;}
                                double analytic = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS.get_mole_fractions(), i, j, xN_flag);
                                double v1 = rHEOS_pluszj.Reducing->ndTrdni__constnj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS_minuszj.Reducing->ndTrdni__constnj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss4;
                            ss4 << "d_ndrhomolarrdni_dxj, i=" << i << ", j=" << j;
                            SECTION(ss4.str(), "")
                            {
                                if (j == z.size()-1){break;}
                                double analytic = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS.get_mole_fractions(), i, j, xN_flag);
                                double v1 = rHEOS_pluszj.Reducing->ndrhorbardni__constnj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS_minuszj.Reducing->ndrhorbardni__constnj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3;
                            ss3 << "d_ndalphardni_dxj__constT_V_xi, i=" << i << ", j=" << j;
                            SECTION(ss3.str(), "")
                            {
                                if (j == z.size()-1){break;}
                                double analytic = MixtureDerivatives::d_ndalphardni_dxj__constT_V_xi(rHEOS, i, j, xN_flag);
                                double v1 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_pluszj, i, xN_flag);
                                double v2 = MixtureDerivatives::ndalphar_dni__constT_V_nj(rHEOS_minuszj, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss3a;
                            ss3a << "d2alphardxidxj, i=" << i << ", j=" << j;
                            SECTION(ss3a.str(), "")
                            {
                                if (j == z.size()-1){break;}
                                double analytic = MixtureDerivatives::d2alphardxidxj(rHEOS,i,j,xN_flag);
                                
                                shared_ptr<HelmholtzEOSMixtureBackend> plus(new HelmholtzEOSMixtureBackend(names));
                                plus->specify_phase(iphase_gas);
                                plus->set_mole_fractions(zp);
                                plus->calc_reducing_state();
                                SimpleState red = plus->get_reducing_state();
                                plus->update(DmolarT_INPUTS, red.rhomolar*rHEOS.delta(), red.T/rHEOS.tau());
                                double v1 = MixtureDerivatives::dalphar_dxi(*(plus.get()), i, xN_flag);
                                shared_ptr<HelmholtzEOSMixtureBackend> minus(new HelmholtzEOSMixtureBackend(names));
                                minus->specify_phase(iphase_gas);
                                minus->set_mole_fractions(zm);
                                minus->calc_reducing_state();
                                red = minus->get_reducing_state();
                                minus->update(DmolarT_INPUTS, red.rhomolar*rHEOS.delta(), red.T/rHEOS.tau());
                                double v2 = MixtureDerivatives::dalphar_dxi(*(minus.get()), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){break;}
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3b;
                            ss3b << "d2Trdxidxj, i=" << i << ", j=" << j;
                            SECTION(ss3b.str(), "")
                            {
                                if (j == z.size()-1 || i == j){break;}
                                double analytic = rHEOS.Reducing->d2Trdxidxj(z, i, j, xN_flag);
                                double v1 = rHEOS.Reducing->dTrdxi__constxj(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS.Reducing->dTrdxi__constxj(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){break;}
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                        }
                    }
                }
            }
        }
    }
}
#endif


