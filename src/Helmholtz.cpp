#include <numeric>
#include "Helmholtz.h"

namespace CoolProp{

CoolPropDbl kahanSum(std::vector<CoolPropDbl> &x)
{
    CoolPropDbl sum = x[0], y, t;
    CoolPropDbl c = 0.0;          //A running compensation for lost low-order bits.
    for (std::size_t i = 1; i < x.size(); ++i)
    {
        y = x[i] - c;    //So far, so good: c is zero.
        t = sum + y;         //Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y;   //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
        sum = t;             //Algebraically, c should always be zero. Beware eagerly optimising compilers!
    }
    return sum;
}
bool wayToSort(CoolPropDbl i, CoolPropDbl j) { return std::abs(i) > std::abs(j); }

// define function to be applied coefficient-wise
double ramp(double x)
{
    if (x > 0)
        return x;
    else
        return 0;
}

/*
void ResidualHelmholtzGeneralizedExponential::allEigen(const CoolPropDbl &tau, const CoolPropDbl &delta, HelmholtzDerivatives &derivs) throw()
{
    double log_tau = log(tau), log_delta = log(delta), 
           one_over_delta = 1/delta, one_over_tau = 1/tau; // division is much slower than multiplication, so do one division here
    
    Eigen::Map<Eigen::ArrayXd> nE(&(n[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> dE(&(d[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> tE(&(t[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> cE(&(c[0]), elements.size());
    Eigen::Map<Eigen::ArrayXi> l_intE(&(l_int[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> l_doubleE(&(l_double[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> eta1E(&(eta1[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> eta2E(&(eta2[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> epsilon1E(&(epsilon1[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> epsilon2E(&(epsilon2[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> beta1E(&(beta1[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> beta2E(&(beta2[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> gamma1E(&(gamma1[0]), elements.size());
    Eigen::Map<Eigen::ArrayXd> gamma2E(&(gamma2[0]), elements.size());
    
    // ****************************************
    // The u part in exp(u) and its derivatives
    // ****************************************
    
    #if defined(EIGEN_VECTORIZE_SSE2)
        //std::cout << "EIGEN_VECTORIZE_SSE2" << std::endl;
    #endif

    // Set the u part of exp(u) to zero
    uE.fill(0);
    du_ddeltaE.fill(0);
    du_dtauE.fill(0);
    d2u_ddelta2E.fill(0);
    d2u_dtau2E.fill(0);
    d3u_ddelta3E.fill(0);
    d3u_dtau3E.fill(0);
    
    if (delta_li_in_u){
        Eigen::ArrayXd u_increment = -cE*(log_delta*l_doubleE).exp(); //pow(delta,L) -> exp(L*log(delta))
        uE += u_increment;
        du_ddeltaE += l_doubleE*u_increment*one_over_delta;
        d2u_ddelta2E += (l_doubleE-1)*l_doubleE*u_increment*one_over_delta*one_over_delta;
        d3u_ddelta3E += (l_doubleE-2)*(l_doubleE-1)*l_doubleE*u_increment*one_over_delta*one_over_delta*one_over_delta;
    }
    
//    if (tau_mi_in_u){
//        CoolPropDbl omegai = el.omega, m_double = el.m_double;
//        if (std::abs(m_double) > 0){
//            CoolPropDbl u_increment = -omegai*pow(tau, m_double);
//            CoolPropDbl du_dtau_increment = m_double*u_increment*one_over_tau;
//            CoolPropDbl d2u_dtau2_increment = (m_double-1)*du_dtau_increment*one_over_tau;
//            CoolPropDbl d3u_dtau3_increment = (m_double-2)*d2u_dtau2_increment*one_over_tau;
//            u += u_increment;
//            du_dtau += du_dtau_increment;
//            d2u_dtau2 += d2u_dtau2_increment;
//            d3u_dtau3 += d3u_dtau3_increment;
//        }
//    }
    if (eta1_in_u){
        uE += -eta1E*(delta-epsilon1E);
        du_ddeltaE += -eta1E;
    }
    if (eta2_in_u){
        uE += -eta2E*POW2(delta-epsilon2E);
        du_ddeltaE += -2*eta2E*(delta-epsilon2E);
        d2u_ddelta2E += -2*eta2E;
    }
    if (beta1_in_u){
        uE += -beta1E*(tau-gamma1E);
        du_dtauE += -beta1E;
    }
    if (beta2_in_u){
        uE += -beta2E*POW2(tau-gamma2E);
        du_dtauE += -2*beta2E*(tau-gamma2E);
        d2u_dtau2E += -2*beta2E;
    }
    
    Eigen::ArrayXd ndteuE = nE*exp(tE*log_tau + dE*log_delta + uE);
    Eigen::ArrayXd B_deltaE = delta*du_ddeltaE + dE;
    Eigen::ArrayXd B_tauE = tau*du_dtauE + tE;
    Eigen::ArrayXd B_delta2E = POW2(delta)*(d2u_ddelta2E + du_ddeltaE.square()) + 2*dE*delta*du_ddeltaE + dE*(dE-1);
    Eigen::ArrayXd B_tau2E = POW2(tau)*(d2u_dtau2E + du_dtauE.square()) + 2*tE*tau*du_dtauE + tE*(tE-1);
    Eigen::ArrayXd B_delta3E = POW3(delta)*d3u_ddelta3E + 3*dE*POW2(delta)*d2u_ddelta2E+3*POW3(delta)*d2u_ddelta2E*du_ddeltaE+3*dE*POW2(delta*du_ddeltaE)+3*dE*(dE-1)*delta*du_ddeltaE+dE*(dE-1)*(dE-2)+POW3(delta*du_ddeltaE);
    Eigen::ArrayXd B_tau3E = POW3(tau)*d3u_dtau3E + 3*tE*POW2(tau)*d2u_dtau2E+3*POW3(tau)*d2u_dtau2E*du_dtauE+3*tE*POW2(tau*du_dtauE)+3*tE*(tE-1)*tau*du_dtauE+tE*(tE-1)*(tE-2)+POW3(tau*du_dtauE);
    
    derivs.alphar                +=  ndteuE.sum();
    derivs.dalphar_ddelta        += (ndteuE*B_deltaE).sum()*one_over_delta;
    derivs.dalphar_dtau          += (ndteuE*B_tauE).sum()*one_over_tau;
    derivs.d2alphar_ddelta2      += (ndteuE*B_delta2E).sum()*POW2(one_over_delta);
    derivs.d2alphar_dtau2        += (ndteuE*B_tau2E).sum()*POW2(one_over_tau);
    derivs.d2alphar_ddelta_dtau  += (ndteuE*B_deltaE*B_tauE).sum()*one_over_delta*one_over_tau;
    
    derivs.d3alphar_ddelta3      += (ndteuE*B_delta3E).sum()*POW3(one_over_delta);
    derivs.d3alphar_dtau3        += (ndteuE*B_tau3E).sum()*POW3(one_over_tau);
    derivs.d3alphar_ddelta2_dtau += (ndteuE*B_delta2E*B_tauE).sum()*POW2(one_over_delta)*one_over_tau;
    derivs.d3alphar_ddelta_dtau2 += (ndteuE*B_deltaE*B_tau2E).sum()*one_over_delta*POW2(one_over_tau);
    
    return;
};
*/
void ResidualHelmholtzGeneralizedExponential::all(const CoolPropDbl &tau, const CoolPropDbl &delta, HelmholtzDerivatives &derivs) throw()
{
    CoolPropDbl log_tau = log(tau), log_delta = log(delta), ndteu, 
                one_over_delta = 1/delta, one_over_tau = 1/tau; // division is much slower than multiplication, so do one division here
    
    // Maybe split the construction of u and other parts into two separate loops?  
    // If both loops can get vectorized, could be worth it.
    const std::size_t N = elements.size();
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzGeneralizedExponentialElement &el = elements[i];
        CoolPropDbl ni = el.n, di = el.d, ti = el.t;
        
        // Set the u part of exp(u) to zero
        CoolPropDbl u = 0;
        CoolPropDbl du_ddelta = 0;
        CoolPropDbl du_dtau = 0;
        CoolPropDbl d2u_ddelta2 = 0;
        CoolPropDbl d2u_dtau2 = 0;
        CoolPropDbl d3u_ddelta3 = 0;
        CoolPropDbl d3u_dtau3 = 0;
        
        if (delta_li_in_u){
            CoolPropDbl  ci = el.c, l_double = el.l_double;
            int l_int = el.l_int;
            if (ValidNumber(l_double) && l_int > 0){
                CoolPropDbl u_increment = -ci*pow(delta, l_int);
                CoolPropDbl du_ddelta_increment = l_double*u_increment*one_over_delta;
                CoolPropDbl d2u_ddelta2_increment = (l_double-1)*du_ddelta_increment*one_over_delta;
                CoolPropDbl d3u_ddelta3_increment = (l_double-2)*d2u_ddelta2_increment*one_over_delta;
                u += u_increment;
                du_ddelta += du_ddelta_increment;
                d2u_ddelta2 += d2u_ddelta2_increment;
                d3u_ddelta3 += d3u_ddelta3_increment;
            }
        }
        if (tau_mi_in_u){
            CoolPropDbl omegai = el.omega, m_double = el.m_double;
            if (std::abs(m_double) > 0){
                CoolPropDbl u_increment = -omegai*pow(tau, m_double);
                CoolPropDbl du_dtau_increment = m_double*u_increment*one_over_tau;
                CoolPropDbl d2u_dtau2_increment = (m_double-1)*du_dtau_increment*one_over_tau;
                CoolPropDbl d3u_dtau3_increment = (m_double-2)*d2u_dtau2_increment*one_over_tau;
                u += u_increment;
                du_dtau += du_dtau_increment;
                d2u_dtau2 += d2u_dtau2_increment;
                d3u_dtau3 += d3u_dtau3_increment;
            }
        }
        if (eta1_in_u){
            CoolPropDbl eta1 = el.eta1, epsilon1 = el.epsilon1;
            if (ValidNumber(eta1)){
                u += -eta1*(delta-epsilon1);
                du_ddelta += -eta1;
            }
        }
        if (eta2_in_u){
            CoolPropDbl eta2 = el.eta2, epsilon2 = el.epsilon2;
            if (ValidNumber(eta2)){
                u += -eta2*POW2(delta-epsilon2);
                du_ddelta += -2*eta2*(delta-epsilon2);
                d2u_ddelta2 += -2*eta2;
            }
        }
        if (beta1_in_u){
            CoolPropDbl beta1 = el.beta1, gamma1 = el.gamma1;
            if (ValidNumber(beta1)){
                u += -beta1*(tau-gamma1);
                du_dtau += -beta1;
            }
        }
        if (beta2_in_u){
            CoolPropDbl beta2 = el.beta2, gamma2 = el.gamma2;
            if (ValidNumber(beta2)){
                u += -beta2*POW2(tau-gamma2);
                du_dtau += -2*beta2*(tau-gamma2);
                d2u_dtau2 += -2*beta2;
            }
        }
        
        ndteu = ni*exp(ti*log_tau+di*log_delta+u);
        
        CoolPropDbl B_delta = (delta*du_ddelta + di);
        CoolPropDbl B_tau = (tau*du_dtau + ti);
        CoolPropDbl B_delta2 = (POW2(delta)*(d2u_ddelta2 + POW2(du_ddelta)) + 2*di*delta*du_ddelta + di*(di-1));
        CoolPropDbl B_tau2 = (POW2(tau)*(d2u_dtau2 + POW2(du_dtau)) + 2*ti*tau*du_dtau + ti*(ti-1));
        CoolPropDbl B_delta3 = POW3(delta)*d3u_ddelta3 + 3*di*POW2(delta)*d2u_ddelta2+3*POW3(delta)*d2u_ddelta2*du_ddelta+3*di*POW2(delta*du_ddelta)+3*di*(di-1)*delta*du_ddelta+di*(di-1)*(di-2)+POW3(delta*du_ddelta);
        CoolPropDbl B_tau3 = POW3(tau)*d3u_dtau3 + 3*ti*POW2(tau)*d2u_dtau2+3*POW3(tau)*d2u_dtau2*du_dtau+3*ti*POW2(tau*du_dtau)+3*ti*(ti-1)*tau*du_dtau+ti*(ti-1)*(ti-2)+POW3(tau*du_dtau);
        
        derivs.alphar += ndteu;
        
        derivs.dalphar_ddelta += ndteu*B_delta;
        derivs.dalphar_dtau += ndteu*B_tau;
        
        derivs.d2alphar_ddelta2 += ndteu*B_delta2;
        derivs.d2alphar_ddelta_dtau += ndteu*B_delta*B_tau;
        derivs.d2alphar_dtau2 += ndteu*B_tau2;
        
        derivs.d3alphar_ddelta3 += ndteu*B_delta3;
        derivs.d3alphar_ddelta2_dtau += ndteu*B_delta2*B_tau;
        derivs.d3alphar_ddelta_dtau2 += ndteu*B_delta*B_tau2;
        derivs.d3alphar_dtau3 += ndteu*B_tau3;
    }
    derivs.dalphar_ddelta        *= one_over_delta;
    derivs.dalphar_dtau          *= one_over_tau;
    derivs.d2alphar_ddelta2      *= POW2(one_over_delta);
    derivs.d2alphar_dtau2        *= POW2(one_over_tau);
    derivs.d2alphar_ddelta_dtau  *= one_over_delta*one_over_tau;
    
    derivs.d3alphar_ddelta3      *= POW3(one_over_delta);
    derivs.d3alphar_dtau3        *= POW3(one_over_tau);
    derivs.d3alphar_ddelta2_dtau *= POW2(one_over_delta)*one_over_tau;
    derivs.d3alphar_ddelta_dtau2 *= one_over_delta*POW2(one_over_tau);
    
    
    
    return;
};

void ResidualHelmholtzNonAnalytic::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzNonAnalytic",doc.GetAllocator());

    rapidjson::Value _n(rapidjson::kArrayType), _a(rapidjson::kArrayType), _b(rapidjson::kArrayType),
                     _beta(rapidjson::kArrayType), _A(rapidjson::kArrayType), _B(rapidjson::kArrayType),
                     _C(rapidjson::kArrayType), _D(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &elem = elements[i];
        _n.PushBack((double)elem.n, doc.GetAllocator());
        _a.PushBack((double)elem.a, doc.GetAllocator());
        _b.PushBack((double)elem.b, doc.GetAllocator());
        _beta.PushBack((double)elem.beta, doc.GetAllocator());
        _A.PushBack((double)elem.A, doc.GetAllocator());
        _B.PushBack((double)elem.B, doc.GetAllocator());
        _C.PushBack((double)elem.C, doc.GetAllocator());
        _D.PushBack((double)elem.D, doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("a",_a,doc.GetAllocator());
    el.AddMember("b",_b,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("A",_A,doc.GetAllocator());
    el.AddMember("B",_B,doc.GetAllocator());
    el.AddMember("C",_C,doc.GetAllocator());
    el.AddMember("D",_D,doc.GetAllocator());
}

void ResidualHelmholtzNonAnalytic::all(const CoolPropDbl &tau, const CoolPropDbl &delta, HelmholtzDerivatives &derivs) throw()
{
    if (N==0){return;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl theta = (1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        
        CoolPropDbl PSI = exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta = -2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dPSI_dTau = -2.0*Di*(tau-1.0)*PSI;
        CoolPropDbl dPSI2_dDelta2 = (2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        CoolPropDbl dPSI2_dDelta_dTau = 4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        CoolPropDbl dPSI2_dTau2 = (2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        CoolPropDbl dPSI3_dDelta3 = 2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        CoolPropDbl dPSI3_dDelta2_dTau = (2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*dPSI_dTau;
        CoolPropDbl dPSI3_dDelta_dTau2 = 2*Di*(2*Di*pow(tau-1,2)-1)*dPSI_dDelta;
        CoolPropDbl dPSI3_dTau3 = 2.0*Di*PSI*(-4*Di*Di*pow(tau-1,3)+6*Di*(tau-1));
        
        CoolPropDbl DELTA = pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl dDELTA_dTau = -2*theta;
        CoolPropDbl dDELTA2_dDelta_dTau = 2.0*Ai/(betai)*pow(pow(delta-1,2),1.0/(2.0*betai)-0.5);
        CoolPropDbl dDELTA_dDelta = (delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        CoolPropDbl dDELTA3_dDelta2_dTau = 2.0*Ai*(betai-1)/(betai*betai)*pow(pow(delta-1,2),1/(2*betai)-1.0);
        
        CoolPropDbl dDELTAbi_dDelta, dDELTA2_dDelta2, dDELTAbi2_dDelta2, dDELTAbi3_dDelta3, dDELTA3_dDelta3;
        if (std::abs(delta-1) < 10*DBL_EPSILON){
            dDELTAbi_dDelta = 0;
            dDELTA2_dDelta2 = 0;
            dDELTA3_dDelta3 = 0;
            dDELTAbi2_dDelta2 = 0;
            dDELTAbi3_dDelta3 = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
            CoolPropDbl dDELTA_dDelta_over_delta_minus_1=(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
            dDELTA2_dDelta2 = dDELTA_dDelta_over_delta_minus_1+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
            CoolPropDbl PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
            CoolPropDbl dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
            dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
            dDELTAbi2_dDelta2 = bi*(pow(DELTA,bi-1)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
            dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));
        }
        
        CoolPropDbl dDELTAbi_dTau = -2.0*theta*bi*pow(DELTA,bi-1.0);
        
        CoolPropDbl dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;
        CoolPropDbl dDELTAbi2_dTau2 = 2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);
        CoolPropDbl dDELTAbi3_dTau3 = -12.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2)-8*pow(theta,3)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3);
        CoolPropDbl dDELTAbi3_dDelta_dTau2 = 2*bi*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+4*pow(theta,2)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta+8*theta*bi*(bi-1)*pow(DELTA,bi-2)*dtheta_dDelta;
        CoolPropDbl dDELTAbi3_dDelta2_dTau = bi*((bi-1)*pow(DELTA,bi-2)*dDELTA_dTau*dDELTA2_dDelta2 + pow(DELTA,bi-1)*dDELTA3_dDelta2_dTau+(bi-1)*((bi-2)*pow(DELTA,bi-3)*dDELTA_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau));
        
        derivs.alphar += ni*pow(DELTA, bi)*delta*PSI;
        
        derivs.dalphar_ddelta += ni*(pow(DELTA,bi)*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
        derivs.dalphar_dtau += ni*delta*(dDELTAbi_dTau*PSI+pow(DELTA,bi)*dPSI_dTau);
        
        derivs.d2alphar_ddelta2 += ni*(pow(DELTA,bi)*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
        derivs.d2alphar_ddelta_dtau += ni*(pow(DELTA,bi)*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
        derivs.d2alphar_dtau2 += ni*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
        
        derivs.d3alphar_ddelta3 += ni*(pow(DELTA,bi)*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
        CoolPropDbl Line1 = pow(DELTA,bi)*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
        CoolPropDbl Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
        CoolPropDbl Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        derivs.d3alphar_ddelta2_dtau += ni*(Line1+Line2+Line3);
        derivs.d3alphar_ddelta_dtau2 += ni*delta*(dDELTAbi2_dTau2*dPSI_dDelta+dDELTAbi3_dDelta_dTau2*PSI+2*dDELTAbi_dTau*dPSI2_dDelta_dTau+2.0*dDELTAbi2_dDelta_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI3_dDelta_dTau2+dDELTAbi_dDelta*dPSI2_dTau2)+ni*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
        derivs.d3alphar_dtau3 += ni*delta*(dDELTAbi3_dTau3*PSI+(3.0*dDELTAbi2_dTau2)*dPSI_dTau+(3*dDELTAbi_dTau)*dPSI2_dTau2+pow(DELTA,bi)*dPSI3_dTau3);
    }
}

CoolPropDbl ResidualHelmholtzNonAnalytic::base(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));

        s[i] = ni*pow(DELTA, bi)*delta*PSI;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        s[i] = ni*(pow(DELTA,bi)*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);

        s[i] = ni*delta*(dDELTAbi_dTau*PSI+pow(DELTA,bi)*dPSI_dTau);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    CoolPropDbl dDELTA2_dDelta2, dDELTAbi2_dDelta2;
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dDELTA_dDelta_over_delta_minus_1=(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;

        if (std::abs(delta-1) < 10*DBL_EPSILON){
            dDELTA2_dDelta2 = 0;
            dDELTAbi2_dDelta2 = 0;
        }
        else{
            dDELTA2_dDelta2 = dDELTA_dDelta_over_delta_minus_1+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
            dDELTAbi2_dDelta2 = bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        }

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta_dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        CoolPropDbl dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dTau2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        CoolPropDbl theta = (1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA = pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI = exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dTau = -2.0*Di*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi_dTau = -2.0*theta*bi*pow(DELTA,bi-1.0);
        CoolPropDbl dPSI2_dTau2 = (2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        CoolPropDbl dDELTAbi2_dTau2 = 2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        s[i] = ni*delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta3(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        CoolPropDbl dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        CoolPropDbl dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));

        CoolPropDbl dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        CoolPropDbl dPSI3_dDelta3=2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        CoolPropDbl PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
        CoolPropDbl dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
        CoolPropDbl dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
        CoolPropDbl dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta_dTau2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        CoolPropDbl dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;

        CoolPropDbl dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);

        CoolPropDbl dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        CoolPropDbl dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        CoolPropDbl dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        CoolPropDbl dPSI3_dDelta_dTau2 = 2*Di*(2*Di*pow(tau-1,2)-1)*dPSI_dDelta;
        CoolPropDbl dDELTAbi3_dDelta_dTau2 = 2*bi*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+4*pow(theta,2)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta+8*theta*bi*(bi-1)*pow(DELTA,bi-2)*dtheta_dDelta;

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*delta*(dDELTAbi2_dTau2*dPSI_dDelta+dDELTAbi3_dDelta_dTau2*PSI+2*dDELTAbi_dTau*dPSI2_dDelta_dTau+2.0*dDELTAbi2_dDelta_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI3_dDelta_dTau2+dDELTAbi_dDelta*dPSI2_dTau2)+ni*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

CoolPropDbl ResidualHelmholtzNonAnalytic::dDelta2_dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        CoolPropDbl dDELTAbi_dDelta;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        CoolPropDbl dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        CoolPropDbl dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        CoolPropDbl dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;

        CoolPropDbl dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        CoolPropDbl dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        CoolPropDbl dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));

        CoolPropDbl dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (std::abs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        //Following Terms added for this derivative
        CoolPropDbl dPSI3_dDelta2_dTau = (2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*dPSI_dTau;
        CoolPropDbl dDELTA_dTau = -2*theta;
        CoolPropDbl dDELTA2_dDelta_dTau = 2.0*Ai/(betai)*pow(pow(delta-1,2),1.0/(2.0*betai)-0.5);
        CoolPropDbl dDELTA3_dDelta2_dTau = 2.0*Ai*(betai-1)/(betai*betai)*pow(pow(delta-1,2),1/(2*betai)-1.0);

        CoolPropDbl dDELTAbim1_dTau = (bi-1)*pow(DELTA,bi-2)*dDELTA_dTau;
        CoolPropDbl dDELTAbim2_dTau = (bi-2)*pow(DELTA,bi-3)*dDELTA_dTau;
        CoolPropDbl Line11 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,bi-1)*dDELTA3_dDelta2_dTau;
        CoolPropDbl Line21 = (bi-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
        CoolPropDbl dDELTAbi3_dDelta2_dTau = bi*(Line11+Line21);

        CoolPropDbl Line1 = pow(DELTA,bi)*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
        CoolPropDbl Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
        CoolPropDbl Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        s[i] = ni*(Line1+Line2+Line3);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
CoolPropDbl ResidualHelmholtzNonAnalytic::dTau3(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        CoolPropDbl theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        CoolPropDbl DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        CoolPropDbl PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        CoolPropDbl dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        CoolPropDbl dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        CoolPropDbl dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        CoolPropDbl dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);
        CoolPropDbl dPSI3_dTau3=2.0*Di*PSI*(-4*Di*Di*pow(tau-1,3)+6*Di*(tau-1));
        CoolPropDbl dDELTAbi3_dTau3 = -12.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2)-8*pow(theta,3)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3.0);

        s[i] = ni*delta*(dDELTAbi3_dTau3*PSI+(3.0*dDELTAbi2_dTau2)*dPSI_dTau+(3*dDELTAbi_dTau )*dPSI2_dTau2+pow(DELTA,bi)*dPSI3_dTau3);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

void ResidualHelmholtzSAFTAssociating::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzSAFTAssociating",doc.GetAllocator());
    el.AddMember("a",a,doc.GetAllocator());
    el.AddMember("m",m,doc.GetAllocator());
    el.AddMember("epsilonbar",epsilonbar,doc.GetAllocator());
    el.AddMember("vbarn",vbarn,doc.GetAllocator());
    el.AddMember("kappabar",kappabar,doc.GetAllocator());
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::Deltabar(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->g(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dDeltabar_ddelta__consttau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->dg_deta(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta2__consttau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->d2g_deta2(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dDeltabar_dtau__constdelta(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*this->epsilonbar;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_dtau2__constdelta(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta_dtau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->dg_deta(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_dtau3__constdelta(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)3);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta_dtau2(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->dg_deta(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2)*this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta2_dtau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->d2g_deta2(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*pow(this->vbarn,(int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta3__consttau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    return this->d3g_deta3(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)3);
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::X(const CoolPropDbl &delta, const CoolPropDbl &Deltabar)
{
    return 2/(sqrt(1+4*Deltabar*delta)+1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_dDeltabar__constdelta(const CoolPropDbl &delta, const CoolPropDbl &Deltabar)
{
    CoolPropDbl X = this->X(delta,Deltabar);
    return -delta*X*X/(2*Deltabar*delta*X+1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_ddelta__constDeltabar(const CoolPropDbl &delta, const CoolPropDbl &Deltabar)
{
    CoolPropDbl X = this->X(delta,Deltabar);
    return -Deltabar*X*X/(2*Deltabar*delta*X+1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_dtau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    return this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_dtau__constdelta(tau, delta);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_ddelta(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
           + this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_ddelta__consttau(tau, delta));
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_dtau2(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl beta = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl d_dXdtau_dbeta = -delta*X*X/(2*Deltabar*delta*X+1);
    CoolPropDbl d_dXdtau_dDeltabar = 2*delta*delta*X*X*X/pow(2*Deltabar*delta*X+1,2)*beta;
    CoolPropDbl d_dXdtau_dX = -2*beta*delta*X*(Deltabar*delta*X+1)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXdtau_dX*dX_dDeltabar*beta+d_dXdtau_dDeltabar*beta+d_dXdtau_dbeta*dbeta_dtau;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_ddeltadtau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl beta = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    CoolPropDbl d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXddelta_dX*dX_dDeltabar*beta+d_dXddelta_dDeltabar*beta+d_dXddelta_dalpha*dalpha_dtau;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_ddelta2(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);

    CoolPropDbl dX_ddelta_constall = X*X*(2*Deltabar*Deltabar*X-alpha)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    CoolPropDbl d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);

    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Deltabar);

    CoolPropDbl val = (dX_ddelta_constall
            + d_dXddelta_dX*dX_ddelta
            + d_dXddelta_dX*dX_dDeltabar*alpha
            + d_dXddelta_dDeltabar*alpha
            + d_dXddelta_dalpha*dalpha_ddelta);
    return val;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_dtau3(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl Delta_ttt = this->d3Deltabar_dtau3__constdelta(tau, delta);
    CoolPropDbl dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    CoolPropDbl dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_dX*dX_dDelta*Delta_t+dXtt_dDelta*Delta_t + dXtt_dDelta_t*Delta_tt + dXtt_dDelta_tt*Delta_ttt;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddeltadtau2(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl Delta_dtt = this->d3Deltabar_ddelta_dtau2(tau, delta);
    CoolPropDbl dXtt_ddelta = pow(X, 2)*(-12*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 2*pow(Delta_t, 2)*X*delta*(-Delta*X*delta + 2)*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + 2*X*delta*(Delta*Delta_tt + 2*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    CoolPropDbl dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_ddelta + dXtt_dX*dX_ddelta + dXtt_dX*dX_dDelta*Delta_d + dXtt_dDelta*Delta_d + dXtt_dDelta_t*Delta_dt + dXtt_dDelta_tt*Delta_dtt;
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddelta2dtau(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    CoolPropDbl Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl Delta_ddt = this->d3Deltabar_ddelta2_dtau(tau, delta);
    CoolPropDbl dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    CoolPropDbl dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    CoolPropDbl dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_dX*dX_dDelta*Delta_t + dXdd_dDelta*Delta_t + dXdd_dDelta_d*Delta_dt + dXdd_dDelta_dd*Delta_ddt;
}

double Xdd(double X, double delta, double Delta, double Delta_d, double Delta_dd)
{
    return Delta*pow(X, 2)*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*delta*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*(2*Delta_d*X*pow(delta, 2) - 1)/pow(2*Delta*X*delta + 1, 2) - Delta_dd*pow(X, 2)*delta/(2*Delta*X*delta + 1) + pow(X, 2)*(2*pow(Delta, 2)*X - Delta_d)/pow(2*Delta*X*delta + 1, 2);
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddelta3(const CoolPropDbl &tau, const CoolPropDbl &delta)
{
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    CoolPropDbl Delta_ddd = this->d3Deltabar_ddelta3__consttau(tau, delta);

    CoolPropDbl dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    CoolPropDbl dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    CoolPropDbl dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    CoolPropDbl dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    CoolPropDbl dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_ddelta + dXdd_dX*(dX_ddelta + dX_dDelta*Delta_d) + dXdd_dDelta*Delta_d + dXdd_dDelta_d*Delta_dd + dXdd_dDelta_dd*Delta_ddd;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::g(const CoolPropDbl &eta)
{
    return 0.5*(2-eta)/pow(1-eta,(int)3);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dg_deta(const CoolPropDbl &eta)
{
    return 0.5*(5-2*eta)/pow(1-eta,(int)4);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2g_deta2(const CoolPropDbl &eta)
{
    return 3*(3-eta)/pow(1-eta,(int)5);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3g_deta3(const CoolPropDbl &eta)
{
    return 6*(7-2*eta)/pow(1-eta,(int)6);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::eta(const CoolPropDbl &delta){
    return this->vbarn*delta;
}

void ResidualHelmholtzSAFTAssociating::all(const CoolPropDbl &tau, const CoolPropDbl &delta, HelmholtzDerivatives &deriv) throw()
{
    if (disabled){return;}
    CoolPropDbl X = this->X(delta, this->Deltabar(tau, delta));
    CoolPropDbl X_t = this->dX_dtau(tau, delta);
    CoolPropDbl X_d = this->dX_ddelta(tau, delta);
    CoolPropDbl X_tt = this->d2X_dtau2(tau, delta);
    CoolPropDbl X_dd = this->d2X_ddelta2(tau, delta);
    CoolPropDbl X_dt = this->d2X_ddeltadtau(tau, delta);
    CoolPropDbl X_ttt = this->d3X_dtau3(tau, delta);
    CoolPropDbl X_dtt = this->d3X_ddeltadtau2(tau, delta);
    CoolPropDbl X_ddt = this->d3X_ddelta2dtau(tau, delta);
    CoolPropDbl X_ddd = this->d3X_ddelta3(tau, delta);
    
    deriv.alphar += this->m*this->a*((log(X)-X/2.0+0.5));
    deriv.dalphar_ddelta += this->m*this->a*(1/X-0.5)*this->dX_ddelta(tau, delta);
    deriv.dalphar_dtau += this->m*this->a*(1/X-0.5)*this->dX_dtau(tau, delta);
    deriv.d2alphar_dtau2 += this->m*this->a*((1/X-0.5)*X_tt-pow(X_t/X, 2));
    deriv.d2alphar_ddelta2 += this->m*this->a*((1/X-0.5)*X_dd-pow(X_d/X,2));
    deriv.d2alphar_ddelta_dtau += this->m*this->a*((-X_t/X/X)*X_d + X_dt*(1/X-0.5));
    deriv.d3alphar_dtau3 += this->m*this->a*((1/X-1.0/2.0)*X_ttt+(-X_t/pow(X,(int)2))*X_tt-2*(pow(X,(int)2)*(X_t*X_tt)-pow(X_t,(int)2)*(X*X_t))/pow(X,(int)4));
    deriv.d3alphar_ddelta_dtau2 += this->m*this->a*((1/X-1.0/2.0)*X_dtt-X_d/pow(X,(int)2)*X_tt-2*(pow(X,(int)2)*(X_t*X_dt)-pow(X_t,(int)2)*(X*X_d))/pow(X,(int)4));
    deriv.d3alphar_ddelta2_dtau += this->m*this->a*((1/X-1.0/2.0)*X_ddt-X_t/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dt)-pow(X_d,(int)2)*(X*X_t))/pow(X,(int)4));
    deriv.d3alphar_ddelta3 += this->m*this->a*((1/X-1.0/2.0)*X_ddd-X_d/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dd)-pow(X_d,(int)2)*(X*X_d))/pow(X,(int)4));
}

void IdealHelmholtzCP0PolyT::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","IdealGasCP0Poly", doc.GetAllocator());

    rapidjson::Value _c(rapidjson::kArrayType), _t(rapidjson::kArrayType);
    for (std::size_t i=0; i< N; ++i)
    {
        _c.PushBack(static_cast<double>(c[i]), doc.GetAllocator());
        _t.PushBack(static_cast<double>(t[i]), doc.GetAllocator());
    }
    el.AddMember("c",_c,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("Tc", static_cast<double>(Tc), doc.GetAllocator());
    el.AddMember("T0", static_cast<double>(T0), doc.GetAllocator());
}
CoolPropDbl IdealHelmholtzCP0PolyT::base(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (std::abs(t[i])<10*DBL_EPSILON)
        {
            sum += c[i]-c[i]*tau/tau0+c[i]*log(tau/tau0);
        }
        else if (std::abs(t[i]+1) < 10*DBL_EPSILON)
        {
            sum += c[i]*tau/Tc*log(tau0/tau)+c[i]/Tc*(tau-tau0);
        }
        else
        {
            sum += -c[i]*pow(Tc,t[i])*pow(tau,-t[i])/(t[i]*(t[i]+1))-c[i]*pow(T0,t[i]+1)*tau/(Tc*(t[i]+1))+c[i]*pow(T0,t[i])/t[i];
        }
    }
    return sum;
}
CoolPropDbl IdealHelmholtzCP0PolyT::dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (std::abs(t[i])<10*DBL_EPSILON)
        {
            sum += c[i]/tau-c[i]/tau0;
        }
        else if (std::abs(t[i]+1) < 10*DBL_EPSILON)
        {
            sum += c[i]/Tc*log(tau0/tau);
        }
        else
        {
            sum += c[i]*pow(Tc,t[i])*pow(tau,-t[i]-1)/(t[i]+1)-c[i]*pow(Tc,t[i])/(pow(tau0,t[i]+1)*(t[i]+1));
        }
    }
    return sum;
}
CoolPropDbl IdealHelmholtzCP0PolyT::dTau2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (std::abs(t[i])<10*DBL_EPSILON)
        {
            sum += -c[i]/(tau*tau);
        }
        else if (std::abs(t[i]+1) < 10*DBL_EPSILON)
        {
            sum += -c[i]/(tau*Tc);
        }
        else
        {
            sum += -c[i]*pow(Tc/tau,t[i])/(tau*tau);
        }
    }
    return sum;
}
CoolPropDbl IdealHelmholtzCP0PolyT::dTau3(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (std::abs(t[i])<10*DBL_EPSILON)
        {
            sum += 2*c[i]/(tau*tau*tau);
        }
        else if (std::abs(t[i]+1) < 10*DBL_EPSILON)
        {
            sum += c[i]/(tau*tau*Tc);
        }
        else
        {
            sum += c[i]*pow(Tc/tau,t[i])*(t[i]+2)/(tau*tau*tau);
        }
    }
    return sum;
}

void IdealHelmholtzCP0AlyLee::to_json(rapidjson::Value &el, rapidjson::Document &doc){
    el.AddMember("type","IdealGasHelmholtzCP0AlyLee",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType);
    for (std::size_t i=0; i<=4; ++i)
    {
        _n.PushBack(static_cast<double>(c[i]),doc.GetAllocator());
    }
    el.AddMember("c",_n,doc.GetAllocator());
    el.AddMember("Tc",static_cast<double>(Tc),doc.GetAllocator());
    el.AddMember("T0",static_cast<double>(T0),doc.GetAllocator());
}
CoolPropDbl IdealHelmholtzCP0AlyLee::base(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (!enabled){ return 0.0;}
    return -tau*(anti_deriv_cp0_tau2(tau)-anti_deriv_cp0_tau2(tau0))+(anti_deriv_cp0_tau(tau)-anti_deriv_cp0_tau(tau0));
}
CoolPropDbl IdealHelmholtzCP0AlyLee::dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (!enabled){ return 0.0;}
    return -(anti_deriv_cp0_tau2(tau) - anti_deriv_cp0_tau2(tau0));
}
CoolPropDbl IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau2(const CoolPropDbl &tau)
{
    return -c[0]/tau + 2*c[1]*c[2]/Tc/(exp(-2*c[2]*tau/Tc)-1) - 2*c[3]*c[4]/Tc*(exp(2*c[4]*tau/Tc)+1);
}
CoolPropDbl IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau(const CoolPropDbl &tau)
{
    CoolPropDbl term1 = c[0]*log(tau);
    CoolPropDbl term2 = 2*c[1]*c[2]*tau/(-Tc + Tc*exp(-2*c[2]*tau/Tc)) + c[1]*log(-1 + exp(-2*c[2]*tau/Tc)) + 2*c[1]*c[2]*tau/Tc;
    CoolPropDbl term3 = -c[3]*(Tc*exp(2*c[4]*tau/Tc)*log(exp(2*c[4]*tau/Tc) + 1) + Tc*log(exp(2*c[4]*tau/Tc) + 1) - 2*c[4]*tau*exp(2*c[4]*tau/Tc))/(Tc*(exp(2*c[4]*tau/Tc) + 1));
    return term1 + term2 + term3;
}
CoolPropDbl IdealHelmholtzCP0AlyLee::dTau2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (!enabled){ return 0.0;}
    return -c[0]/pow(tau,2) - c[1]*pow(c[2]/Tc/sinh(c[2]*tau/Tc),2) - c[3]*pow(c[4]/Tc/cosh(c[4]*tau/Tc),2);
}
CoolPropDbl IdealHelmholtzCP0AlyLee::dTau3(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
{
    if (!enabled){ return 0.0;}
    return 2*c[0]/pow(tau,3) + 2*c[1]*pow(c[2]/Tc,3)*cosh(c[2]*tau/Tc)/pow(sinh(c[2]*tau/Tc),3) + 2*c[3]*pow(c[4]/Tc,3)*sinh(c[4]*tau/Tc)/pow(cosh(c[4]*tau/Tc),3);
}

}; /* namespace CoolProp */


/*
IdealHelmholtzEnthalpyEntropyOffset EnthalpyEntropyOffset;
*/

#ifdef ENABLE_CATCH
#include <math.h>
#include "catch.hpp"
#include "crossplatform_shared_ptr.h"

class HelmholtzConsistencyFixture
{
public:
    CoolPropDbl numerical, analytic;

    shared_ptr<CoolProp::BaseHelmholtzTerm> PlanckEinstein, Lead, LogTau, IGPower, CP0Constant, CP0PolyT, SAFT, NonAnalytic;
    shared_ptr<CoolProp::ResidualHelmholtzGeneralizedExponential> Gaussian, Lemmon2005, Exponential, GERG2008, Power;

    HelmholtzConsistencyFixture(){
        Lead.reset(new CoolProp::IdealHelmholtzLead(1,3));
        LogTau.reset(new CoolProp::IdealHelmholtzLogTau(1.5));
        {
            std::vector<CoolPropDbl> n(4,0), t(4,1); n[0] = -0.1; n[2] = 0.1; t[1] = -1; t[2] = -2; t[3] = 2;
            IGPower.reset(new CoolProp::IdealHelmholtzPower(n,t));
        }
        {
            std::vector<CoolPropDbl> n(4,0), t(4,1), c(4,1), d(4,-1); n[0] = 0.1; n[2] = 0.5; t[0] = -1.5; t[1] = -1; t[2] = -2; t[3] = -2;
            PlanckEinstein.reset(new CoolProp::IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d));
        }
        {
            CoolPropDbl T0 = 273.15, Tc = 345.857, c = 1.0578, t = 0.33;
            CP0PolyT.reset(new CoolProp::IdealHelmholtzCP0PolyT(std::vector<CoolPropDbl>(1,c),
                                                                std::vector<CoolPropDbl>(1,t),
                                                                Tc, T0));
        }
        CP0Constant.reset(new CoolProp::IdealHelmholtzCP0Constant(4/8.314472,300,250));
        {
            CoolPropDbl beta[] = {1.24, 0.821, 15.45, 2.21, 437, 0.743},
                        d[] = {1, 1, 2, 2, 3, 3},
                        epsilon[] = {0.6734, 0.9239, 0.8636, 1.0507, 0.8482, 0.7522},
                        eta[] = {0.9667, 1.5154, 1.0591, 1.6642, 12.4856, 0.9662},
                        gamma[] = {1.2827, 0.4317, 1.1217, 1.1871, 1.1243, 0.4203},
                        n[] = {1.2198, -0.4883, -0.0033293, -0.0035387, -0.51172, -0.16882},
                        t[] = {1, 2.124, 0.4, 3.5, 0.5, 2.7};
            Gaussian.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Gaussian->add_Gaussian(std::vector<CoolPropDbl>(n,n+sizeof(n)/sizeof(n[0])),
                                   std::vector<CoolPropDbl>(d,d+sizeof(d)/sizeof(d[0])),
                                   std::vector<CoolPropDbl>(t,t+sizeof(t)/sizeof(t[0])),
                                   std::vector<CoolPropDbl>(eta,eta+sizeof(eta)/sizeof(eta[0])),
                                   std::vector<CoolPropDbl>(epsilon,epsilon+sizeof(epsilon)/sizeof(epsilon[0])),
                                   std::vector<CoolPropDbl>(beta,beta+sizeof(beta)/sizeof(beta[0])),
                                   std::vector<CoolPropDbl>(gamma,gamma+sizeof(gamma)/sizeof(gamma[0]))
                                   );
        }
        {
            CoolPropDbl d[] = {1, 1, 1, 2, 4, 1, 1, 2, 2, 3, 4, 5, 1, 5, 1, 2, 3, 5},
                        l[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 3, 3},
                        m[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.7, 7, 6},
                        n[] = {5.28076, -8.67658, 0.7501127, 0.7590023, 0.01451899, 4.777189, -3.330988, 3.775673, -2.290919, 0.8888268, -0.6234864, -0.04127263, -0.08455389, -0.1308752, 0.008344962, -1.532005, -0.05883649, 0.02296658},
                        t[]= {0.669, 1.05,2.75, 0.956, 1, 2, 2.75, 2.38, 3.37, 3.47, 2.63, 3.45, 0.72, 4.23, 0.2, 4.5, 29, 24};
            Lemmon2005.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Lemmon2005->add_Lemmon2005(std::vector<CoolPropDbl>(n, n+sizeof(n)/sizeof(n[0])),
                                       std::vector<CoolPropDbl>(d, d+sizeof(d)/sizeof(d[0])),
                                       std::vector<CoolPropDbl>(t, t+sizeof(t)/sizeof(t[0])),
                                       std::vector<CoolPropDbl>(l, l+sizeof(l)/sizeof(l[0])),
                                       std::vector<CoolPropDbl>(m, m+sizeof(m)/sizeof(m[0]))
                                       );
        }
        {
            CoolPropDbl d[] = {1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2},
                        l[] = {0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3},
                        n[] = {1.0038, -2.7662, 0.42921, 0.081363, 0.00024174, 0.48246, 0.75542, -0.00743, -0.4146, -0.016558, -0.10644, -0.021704},
                        t[] = {0.25, 1.25, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5};
            Power.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Power->add_Power(std::vector<CoolPropDbl>(n, n+sizeof(n)/sizeof(n[0])),
                             std::vector<CoolPropDbl>(d, d+sizeof(d)/sizeof(d[0])),
                             std::vector<CoolPropDbl>(t, t+sizeof(t)/sizeof(t[0])),
                             std::vector<CoolPropDbl>(l, l+sizeof(l)/sizeof(l[0]))
                             );
        }
        {

            CoolPropDbl a = 1, epsilonbar = 12.2735737, kappabar = 1.09117041e-05, m = 1.01871348, vbarn = 0.0444215309;
            SAFT.reset(new CoolProp::ResidualHelmholtzSAFTAssociating(a,m,epsilonbar,vbarn,kappabar));
        }
        {
            CoolPropDbl n[] = {-0.666422765408, 0.726086323499, 0.0550686686128},
                A[] = {0.7, 0.7, 0.7}, B[] = {0.3, 0.3, 1}, C[] = {10, 10, 12.5}, D[] = {275, 275, 275},
                a[] = {3.5, 3.5, 3}, b[] = {0.875, 0.925, 0.875}, beta[] = {0.3, 0.3, 0.3};
            NonAnalytic.reset(new CoolProp::ResidualHelmholtzNonAnalytic(std::vector<CoolPropDbl>(n, n+sizeof(n)/sizeof(n[0])),
                                                             std::vector<CoolPropDbl>(a, a+sizeof(a)/sizeof(a[0])),
                                                             std::vector<CoolPropDbl>(b, b+sizeof(b)/sizeof(b[0])),
                                                             std::vector<CoolPropDbl>(beta, beta+sizeof(beta)/sizeof(beta[0])),
                                                             std::vector<CoolPropDbl>(A, A+sizeof(A)/sizeof(A[0])),
                                                             std::vector<CoolPropDbl>(B, B+sizeof(B)/sizeof(B[0])),
                                                             std::vector<CoolPropDbl>(C, C+sizeof(C)/sizeof(C[0])),
                                                             std::vector<CoolPropDbl>(D, D+sizeof(D)/sizeof(D[0]))
                                                             ));
        }
        {
            CoolPropDbl d[] = {2, 2, 2, 0, 0, 0},
            g[] = {1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788},
            l[] = {2, 2, 2, 2, 2, 2},
            n[] = {-3.821884669859, 8.30345065618981, -4.4832307260286, -1.02590136933231, 2.20786016506394, -1.07889905203761},
            t[] = {3, 4, 5, 3, 4, 5};
            Exponential.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Exponential->add_Exponential(std::vector<CoolPropDbl>(n, n+sizeof(n)/sizeof(n[0])),
                                         std::vector<CoolPropDbl>(d, d+sizeof(d)/sizeof(n[0])),
                                         std::vector<CoolPropDbl>(t, t+sizeof(t)/sizeof(d[0])),
                                         std::vector<CoolPropDbl>(g, g+sizeof(g)/sizeof(t[0])),
                                         std::vector<CoolPropDbl>(l, l+sizeof(l)/sizeof(l[0]))
                                         );
        }
        {
            CoolPropDbl d[] = {1, 4, 1, 2, 2, 2, 2, 2, 3},
                t[] = {0.0, 1.85, 7.85, 5.4, 0.0, 0.75, 2.8, 4.45, 4.25},
                n[] = {-0.0098038985517335, 0.00042487270143005, -0.034800214576142, -0.13333813013896, -0.011993694974627, 0.069243379775168, -0.31022508148249, 0.24495491753226, 0.22369816716981},
                eta[] = {0.0, 0.0, 1.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0},
                epsilon[] = {0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
                beta[] = {0.0, 0.0, 1.0, 1.0, 2.5, 3.0, 3.0, 3.0, 3.0},
                gamma[] = {0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
            GERG2008.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            GERG2008->add_GERG2008Gaussian(std::vector<CoolPropDbl>(n, n+sizeof(n)/sizeof(n[0])),
                                           std::vector<CoolPropDbl>(d, d+sizeof(d)/sizeof(n[0])),
                                           std::vector<CoolPropDbl>(t, t+sizeof(t)/sizeof(d[0])),
                                           std::vector<CoolPropDbl>(eta, eta+sizeof(eta)/sizeof(eta[0])),
                                           std::vector<CoolPropDbl>(epsilon, epsilon+sizeof(epsilon)/sizeof(epsilon[0])),
                                           std::vector<CoolPropDbl>(beta, beta+sizeof(beta)/sizeof(beta[0])),
                                           std::vector<CoolPropDbl>(gamma, gamma+sizeof(gamma)/sizeof(gamma[0]))
                                           );
        }

    }
    void call(std::string d, shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta)
    {
              if (!d.compare("dTau")){return dTau(term,tau,delta,ddelta);}
        else if (!d.compare("dTau2")){return dTau2(term,tau,delta,ddelta);}
        else if (!d.compare("dTau3")){return dTau3(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta")){ return dDelta(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta2")){return dDelta2(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta3")){return dDelta3(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta_dTau")){ return dDelta_dTau(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta_dTau2")){ return dDelta_dTau2(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta2_dTau")){ return dDelta2_dTau(term,tau,delta,ddelta);}
        else{
            throw CoolProp::ValueError("don't understand deriv type");
        }
    }
    shared_ptr<CoolProp::BaseHelmholtzTerm> get(std::string t)
    {
        if (!t.compare("Lead")){return Lead;}
        else if (!t.compare("LogTau")){return LogTau;}
        else if (!t.compare("IGPower")){return IGPower;}
        else if (!t.compare("PlanckEinstein")){return PlanckEinstein;}
        else if (!t.compare("CP0Constant")){return CP0Constant;}
        else if (!t.compare("CP0PolyT")){return CP0PolyT;}

        else if (!t.compare("Gaussian")){return Gaussian;}
        else if (!t.compare("Lemmon2005")){return Lemmon2005;}
        else if (!t.compare("Power")){return Power;}
        else if (!t.compare("SAFT")){return SAFT;}
        else if (!t.compare("NonAnalytic")){return NonAnalytic;}
        else if (!t.compare("Exponential")){return Exponential;}
        else if (!t.compare("GERG2008")){return GERG2008;}
        else{
            throw CoolProp::ValueError(format("don't understand helmholtz type: %s",t.c_str()));
        }
    }
    void dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau){
        CoolPropDbl term_plus = term->base(tau + dtau, delta);
        CoolPropDbl term_minus = term->base(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau(tau, delta);
    };
    void dTau2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau){
        CoolPropDbl term_plus = term->dTau(tau + dtau, delta);
        CoolPropDbl term_minus = term->dTau(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau2(tau, delta);
    };
    void dTau3(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau){
        CoolPropDbl term_plus = term->dTau2(tau + dtau, delta);
        CoolPropDbl term_minus = term->dTau2(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau3(tau, delta);
    };
    void dDelta(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->base(tau, delta + ddelta);
        CoolPropDbl term_minus = term->base(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta(tau, delta);
    };
    void dDelta2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->dDelta(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta2(tau, delta);
    };
    void dDelta3(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->dDelta2(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta2(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta3(tau, delta);
    };
    void dDelta_dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->dTau(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dTau(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta_dTau(tau, delta);
    };
    void dDelta_dTau2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->dTau2(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dTau2(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta_dTau2(tau, delta);
    };
    void dDelta2_dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta){
        CoolPropDbl term_plus = term->dDelta_dTau(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta_dTau(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta2_dTau(tau, delta);
    };
    double err(double v1, double v2)
    {
        if (std::abs(v2) > 1e-15){
            return std::abs((v1-v2)/v2);
        }
        else{
            return std::abs(v1-v2);
        }
    }
};

std::string terms[] = {"Lead","LogTau","IGPower","PlanckEinstein","CP0Constant","CP0PolyT",
                       "Gaussian","Lemmon2005","Power","SAFT","NonAnalytic","Exponential",
                       "GERG2008"};
std::string derivs[] = {"dTau","dTau2","dTau3","dDelta","dDelta2","dDelta3","dDelta_dTau","dDelta_dTau2","dDelta2_dTau"};

TEST_CASE_METHOD(HelmholtzConsistencyFixture, "Helmholtz energy derivatives", "[helmholtz]")
{
    shared_ptr<CoolProp::BaseHelmholtzTerm> term;
    std::size_t n = sizeof(terms)/sizeof(terms[0]);
    for (std::size_t i = 0; i < n; ++i)
    {
        term = get(terms[i]);
        for (std::size_t j = 0; j < sizeof(derivs)/sizeof(derivs[0]); ++j)
        {
            call(derivs[j], term, 1.3, 0.7, 1e-6);
            CAPTURE(derivs[j]);
            CAPTURE(numerical);
            CAPTURE(analytic);
            CAPTURE(terms[i]);
            CHECK(err(analytic, numerical) < 1e-7);
        }
    }
}

#endif

