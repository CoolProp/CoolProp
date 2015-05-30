#include "MixtureDerivatives.h"

namespace CoolProp{
    
CoolPropDbl MixtureDerivatives::dalphar_dxi(HelmholtzEOSMixtureBackend &HEOS, std::size_t i, x_N_dependency_flag xN_flag)
{
    if (xN_flag == XN_INDEPENDENT){
        return HEOS.components[i].EOS().baser(HEOS._tau, HEOS._delta) + HEOS.Excess.dalphar_dxi(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i == N-1) return 0;
        double dar_dxi = HEOS.components[i].EOS().baser(HEOS._tau, HEOS._delta) - HEOS.components[N-1].EOS().baser(HEOS._tau, HEOS._delta);
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
        return HEOS.components[i].EOS().dalphar_dTau(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dTau(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
    }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return 0;
        double d2ar_dxi_dTau = HEOS.components[i].EOS().dalphar_dTau(HEOS._tau, HEOS._delta) - HEOS.components[N-1].EOS().dalphar_dTau(HEOS._tau, HEOS._delta);
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
        return HEOS.components[i].EOS().dalphar_dDelta(HEOS._tau, HEOS._delta) + HEOS.Excess.d2alphar_dxi_dDelta(HEOS._tau, HEOS._delta, HEOS.mole_fractions, i);
        }
    else if(xN_flag == XN_DEPENDENT){
        std::vector<CoolPropDbl> &x = HEOS.mole_fractions;
        std::size_t N = x.size();
        if (i==N-1) return 0;
        double d2ar_dxi_dDelta = HEOS.components[i].EOS().dalphar_dDelta(HEOS._tau, HEOS._delta) - HEOS.components[N-1].EOS().dalphar_dDelta(HEOS._tau, HEOS._delta);
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

static const double T1 = 300, rho1 = 300, dT = 1e-3, drho = 1e-3, dz = 1e-6;

void setup_state(std::vector<shared_ptr<HelmholtzEOSMixtureBackend> > & HEOS, std::size_t Ncomp, double increment, x_N_dependency_flag xN_flag = XN_INDEPENDENT)
{
    std::vector<std::string> names(Ncomp);
    std::vector<CoolPropDbl> z(Ncomp);
    if (Ncomp == 2){
        names[0] = "Ethane"; names[1] = "Propane";
        z[0] = 0.3; z[1] = 0.7;
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
                    
                    HelmholtzEOSMixtureBackend &rHEOS = *(HEOS[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_plusT_constrho = *(HEOS_plusT_constrho[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_minusT_constrho = *(HEOS_minusT_constrho[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_plusrho_constT = *(HEOS_plusrho_constT[Ncomp][0].get());
                    HelmholtzEOSMixtureBackend &rHEOS_minusrho_constT = *(HEOS_minusrho_constT[Ncomp][0].get());

                    const std::vector<CoolPropDbl> &z = rHEOS.get_mole_fractions();
                
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
                            double err = std::abs((numeric-analytic)/analytic);
                            CHECK(err < 1e-8);
                        }
                        std::ostringstream ss3a;
                        ss3a << "d2alphar_dxi_dDelta, i=" << i;
                        SECTION(ss3a.str(), "")
                        {
                            if (i == Ncomp-1){ break; }
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
                            if (i == Ncomp-1){ break; }
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
                            if (i==Ncomp-1){ break; }
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
                            if (i==Ncomp-1){ break; }
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
                            if (i==Ncomp-1){ break; }
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
                            if (i==Ncomp-1){ break; }
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
                            if (i == Ncomp-1){ break; }
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
                            if (i == Ncomp-1){ break; }
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
                            if (i == Ncomp-1){ break; }
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
                            if (i == Ncomp-1){ break; }
                            double analytic = MixtureDerivatives::dalphar_dxi(rHEOS, i, xN_flag);
                            double v1 = rHEOS_pluszi_consttaudelta.alphar();
                            double v2 = rHEOS_minuszi_consttaudelta.alphar();
                            double numeric = (v1 - v2)/(2*dz);
                            double err = std::abs((numeric-analytic)/analytic);
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
                                if (j == Ncomp-1){break;}
                                double analytic = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(rHEOS, i, j, xN_flag);
                                double v1 = log(MixtureDerivatives::fugacity_i(rHEOS_pluszj, i, xN_flag));
                                double v2 = log(MixtureDerivatives::fugacity_i(rHEOS_minuszj, i, xN_flag));
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-7);
                            }
                            std::ostringstream ss2;
                            ss2 << "d_ndTrdni_dxj, i=" << i << ", j=" << j;
                            SECTION(ss2.str(), "")
                            {
                                if (j == Ncomp-1){ break; }
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
                                if (j == Ncomp-1){ break; }
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
                                if (j == Ncomp-1){ break; }
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
                                if (j == Ncomp-1){ break; }
                                double analytic = MixtureDerivatives::d2alphardxidxj(rHEOS,i,j,xN_flag);
                                double v1 = MixtureDerivatives::dalphar_dxi(rHEOS_pluszj_consttaudelta, i, xN_flag);
                                double v2 = MixtureDerivatives::dalphar_dxi(rHEOS_minuszj_consttaudelta, i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){break;}
                                double err;
                                if (std::abs(analytic) > DBL_EPSILON){
                                    err = std::abs((numeric-analytic)/analytic);
                                }
                                else{
                                    err = numeric-analytic;
                                }
                                
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            std::ostringstream ss3b;
                            ss3b << "d2Trdxidxj, i=" << i << ", j=" << j;
                            SECTION(ss3b.str(), "")
                            {
                                if (j == Ncomp-1 || i == j){ break; }
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
                            std::ostringstream ss3c;
                            ss3c << "d_PSI_rho_dxj, i=" << i << ", j=" << j;
                            SECTION(ss3c.str(), "")
                            {
                                double analytic = rHEOS.Reducing->d_PSI_rho_dxj(z, i, j, xN_flag);
                                double v1 = rHEOS.Reducing->PSI_rho(rHEOS_pluszj.get_mole_fractions(), i, xN_flag);
                                double v2 = rHEOS.Reducing->PSI_rho(rHEOS_minuszj.get_mole_fractions(), i, xN_flag);
                                double numeric = (v1 - v2)/(2*dz);
                                double err = std::abs((numeric-analytic)/analytic);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){ break; }
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
                                double err = std::abs((numeric-analytic)/analytic);
                                if (std::abs(numeric) < DBL_EPSILON && std::abs(analytic) < DBL_EPSILON){ break; }
                                CAPTURE(numeric);
                                CAPTURE(analytic);
                                CHECK(err < 1e-8);
                            }
                            // These derivatives depend on i,j, and k indices
                            for (std::size_t k = 0; k < Ncomp; ++k){
                                if (xN_flag == XN_DEPENDENT && k == Ncomp-1){ continue; }
                                HelmholtzEOSMixtureBackend & rHEOS_pluszk = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_xNindep[Ncomp][k].get()) : *(HEOS_plusz_xNdep[Ncomp][k].get());
                                HelmholtzEOSMixtureBackend & rHEOS_minuszk = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_xNindep[Ncomp][k].get()) : *(HEOS_minusz_xNdep[Ncomp][k].get());

                                HelmholtzEOSMixtureBackend &rHEOS_pluszk_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_plusz_consttaudelta_xNindep[Ncomp][k].get()) : *(HEOS_plusz_consttaudelta_xNdep[Ncomp][k].get());
                                HelmholtzEOSMixtureBackend &rHEOS_minuszk_consttaudelta = (xN_flag == XN_INDEPENDENT) ? *(HEOS_minusz_consttaudelta_xNindep[Ncomp][k].get()) : *(HEOS_minusz_consttaudelta_xNdep[Ncomp][k].get());

                                std::ostringstream ss1;    ss1 << "d3Trdxidxjdxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss1.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d3Trdxidxjdxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d2Trdxidxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d2Trdxidxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss2;    ss2 << "d3rhormolardxidxjdxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss2.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d3rhormolardxidxjdxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d2rhormolardxidxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d2rhormolardxidxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss3; ss3 << "d2_ndTrdni_dxj_dxk__constxi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss3.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d2_ndTrdni_dxj_dxk__constxi(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_ndTrdni_dxj__constxi(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss4; ss4 << "d2_ndrhorbardni_dxj_dxk__constxi, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss4.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d2_ndrhorbardni_dxj_dxk__constxi(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_ndrhorbardni_dxj__constxi(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss5; ss5 << "d2_PSI_T_dxj_dxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss5.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d2_PSI_T_dxj_dxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_PSI_T_dxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_PSI_T_dxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
                                    CAPTURE(numeric);
                                    CAPTURE(analytic);
                                    CHECK(err < 1e-7);
                                }
                                std::ostringstream ss6; ss6 << "d2_PSI_rho_dxj_dxk, i=" << i << ", j=" << j << ", k=" << k;
                                SECTION(ss6.str(), "")
                                {
                                    if (xN_flag == XN_INDEPENDENT){ continue; }
                                    if (j == Ncomp-1){ break; }
                                    double analytic = rHEOS.Reducing->d2_PSI_rho_dxj_dxk(rHEOS.get_mole_fractions(), i, j, k, xN_flag);
                                    double v1 = rHEOS.Reducing->d_PSI_rho_dxj(rHEOS_pluszk.get_mole_fractions(), i, j, xN_flag);
                                    double v2 = rHEOS.Reducing->d_PSI_rho_dxj(rHEOS_minuszk.get_mole_fractions(), i, j, xN_flag);
                                    double numeric = (v1 - v2)/(2*dz);
                                    double err = std::abs((numeric-analytic)/analytic);
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


