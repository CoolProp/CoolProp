#include <numeric>
#include "Helmholtz.h"

#ifdef __ANDROID__
#    undef _A
#    undef _B
#    undef _C
#    undef _D
#endif

namespace CoolProp {

CoolPropDbl kahanSum(const std::vector<CoolPropDbl>& x) {
    CoolPropDbl sum = x[0], y, t;
    CoolPropDbl c = 0.0;  //A running compensation for lost low-order bits.
    for (std::size_t i = 1; i < x.size(); ++i) {
        y = x[i] - c;       //So far, so good: c is zero.
        t = sum + y;        //Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y;  //(t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
        sum = t;            //Algebraically, c should always be zero. Beware eagerly optimising compilers!
    }
    return sum;
}
bool wayToSort(CoolPropDbl i, CoolPropDbl j) {
    return std::abs(i) > std::abs(j);
}

// define function to be applied coefficient-wise
double ramp(double x) {
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
void ResidualHelmholtzGeneralizedExponential::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    CoolPropDbl log_tau = log(tau), log_delta = log(delta), ndteu, one_over_delta = 1 / delta,
                one_over_tau = 1 / tau;  // division is much slower than multiplication, so do one division here

    // Maybe split the construction of u and other parts into two separate loops?
    // If both loops can get vectorized, could be worth it.
    const std::size_t N = elements.size();
    for (std::size_t i = 0; i < N; ++i) {
        ResidualHelmholtzGeneralizedExponentialElement& el = elements[i];
        CoolPropDbl ni = el.n, di = el.d, ti = el.t;

        // Set the u part of exp(u) to zero
        CoolPropDbl u = 0;
        CoolPropDbl du_ddelta = 0;
        CoolPropDbl du_dtau = 0;
        CoolPropDbl d2u_ddelta2 = 0;
        CoolPropDbl d2u_dtau2 = 0;
        CoolPropDbl d3u_ddelta3 = 0;
        CoolPropDbl d3u_dtau3 = 0;
        CoolPropDbl d4u_ddelta4 = 0;
        CoolPropDbl d4u_dtau4 = 0;

        if (delta_li_in_u) {
            CoolPropDbl ci = el.c, l_double = el.l_double;
            if (ValidNumber(l_double) && l_double > 0 && std::abs(ci) > DBL_EPSILON) {
                const CoolPropDbl u_increment = (el.l_is_int) ? -ci * powInt(delta, el.l_int) : -ci * pow(delta, l_double);
                const CoolPropDbl du_ddelta_increment = l_double * u_increment * one_over_delta;
                const CoolPropDbl d2u_ddelta2_increment = (l_double - 1) * du_ddelta_increment * one_over_delta;
                const CoolPropDbl d3u_ddelta3_increment = (l_double - 2) * d2u_ddelta2_increment * one_over_delta;
                const CoolPropDbl d4u_ddelta4_increment = (l_double - 3) * d3u_ddelta3_increment * one_over_delta;
                u += u_increment;
                du_ddelta += du_ddelta_increment;
                d2u_ddelta2 += d2u_ddelta2_increment;
                d3u_ddelta3 += d3u_ddelta3_increment;
                d4u_ddelta4 += d4u_ddelta4_increment;
            }
        }
        if (tau_mi_in_u) {
            CoolPropDbl omegai = el.omega, m_double = el.m_double;
            if (std::abs(m_double) > 0) {
                const CoolPropDbl u_increment = -omegai * pow(tau, m_double);
                const CoolPropDbl du_dtau_increment = m_double * u_increment * one_over_tau;
                const CoolPropDbl d2u_dtau2_increment = (m_double - 1) * du_dtau_increment * one_over_tau;
                const CoolPropDbl d3u_dtau3_increment = (m_double - 2) * d2u_dtau2_increment * one_over_tau;
                const CoolPropDbl d4u_dtau4_increment = (m_double - 3) * d3u_dtau3_increment * one_over_tau;
                u += u_increment;
                du_dtau += du_dtau_increment;
                d2u_dtau2 += d2u_dtau2_increment;
                d3u_dtau3 += d3u_dtau3_increment;
                d4u_dtau4 += d4u_dtau4_increment;
            }
        }
        if (eta1_in_u) {
            CoolPropDbl eta1 = el.eta1, epsilon1 = el.epsilon1;
            if (ValidNumber(eta1)) {
                u += -eta1 * (delta - epsilon1);
                du_ddelta += -eta1;
            }
        }
        if (eta2_in_u) {
            CoolPropDbl eta2 = el.eta2, epsilon2 = el.epsilon2;
            if (ValidNumber(eta2)) {
                u += -eta2 * POW2(delta - epsilon2);
                du_ddelta += -2 * eta2 * (delta - epsilon2);
                d2u_ddelta2 += -2 * eta2;
            }
        }
        if (beta1_in_u) {
            CoolPropDbl beta1 = el.beta1, gamma1 = el.gamma1;
            if (ValidNumber(beta1)) {
                u += -beta1 * (tau - gamma1);
                du_dtau += -beta1;
            }
        }
        if (beta2_in_u) {
            CoolPropDbl beta2 = el.beta2, gamma2 = el.gamma2;
            if (ValidNumber(beta2)) {
                u += -beta2 * POW2(tau - gamma2);
                du_dtau += -2 * beta2 * (tau - gamma2);
                d2u_dtau2 += -2 * beta2;
            }
        }

        ndteu = ni * exp(ti * log_tau + di * log_delta + u);

        const CoolPropDbl dB_delta_ddelta = delta * d2u_ddelta2 + du_ddelta;
        const CoolPropDbl d2B_delta_ddelta2 = delta * d3u_ddelta3 + 2 * d2u_ddelta2;
        const CoolPropDbl d3B_delta_ddelta3 = delta * d4u_ddelta4 + 3 * d3u_ddelta3;

        const CoolPropDbl B_delta = (delta * du_ddelta + di);
        const CoolPropDbl B_delta2 = delta * dB_delta_ddelta + (B_delta - 1) * B_delta;
        const CoolPropDbl dB_delta2_ddelta = delta * d2B_delta_ddelta2 + 2 * B_delta * dB_delta_ddelta;
        const CoolPropDbl B_delta3 = delta * dB_delta2_ddelta + (B_delta - 2) * B_delta2;
        const CoolPropDbl dB_delta3_ddelta = delta * delta * d3B_delta_ddelta3 + 3 * delta * B_delta * d2B_delta_ddelta2
                                             + 3 * delta * POW2(dB_delta_ddelta) + 3 * B_delta * (B_delta - 1) * dB_delta_ddelta;
        const CoolPropDbl B_delta4 = delta * dB_delta3_ddelta + (B_delta - 3) * B_delta3;

        const CoolPropDbl dB_tau_dtau = tau * d2u_dtau2 + du_dtau;
        const CoolPropDbl d2B_tau_dtau2 = tau * d3u_dtau3 + 2 * d2u_dtau2;
        const CoolPropDbl d3B_tau_dtau3 = tau * d4u_dtau4 + 3 * d3u_dtau3;

        const CoolPropDbl B_tau = (tau * du_dtau + ti);
        const CoolPropDbl B_tau2 = tau * dB_tau_dtau + (B_tau - 1) * B_tau;
        const CoolPropDbl dB_tau2_dtau = tau * d2B_tau_dtau2 + 2 * B_tau * dB_tau_dtau;
        const CoolPropDbl B_tau3 = tau * dB_tau2_dtau + (B_tau - 2) * B_tau2;
        const CoolPropDbl dB_tau3_dtau =
          tau * tau * d3B_tau_dtau3 + 3 * tau * B_tau * d2B_tau_dtau2 + 3 * tau * POW2(dB_tau_dtau) + 3 * B_tau * (B_tau - 1) * dB_tau_dtau;
        const CoolPropDbl B_tau4 = tau * dB_tau3_dtau + (B_tau - 3) * B_tau3;

        derivs.alphar += ndteu;

        derivs.dalphar_ddelta += ndteu * B_delta;
        derivs.dalphar_dtau += ndteu * B_tau;

        derivs.d2alphar_ddelta2 += ndteu * B_delta2;
        derivs.d2alphar_ddelta_dtau += ndteu * B_delta * B_tau;
        derivs.d2alphar_dtau2 += ndteu * B_tau2;

        derivs.d3alphar_ddelta3 += ndteu * B_delta3;
        derivs.d3alphar_ddelta2_dtau += ndteu * B_delta2 * B_tau;
        derivs.d3alphar_ddelta_dtau2 += ndteu * B_delta * B_tau2;
        derivs.d3alphar_dtau3 += ndteu * B_tau3;

        derivs.d4alphar_ddelta4 += ndteu * B_delta4;
        derivs.d4alphar_ddelta3_dtau += ndteu * B_delta3 * B_tau;
        derivs.d4alphar_ddelta2_dtau2 += ndteu * B_delta2 * B_tau2;
        derivs.d4alphar_ddelta_dtau3 += ndteu * B_delta * B_tau3;
        derivs.d4alphar_dtau4 += ndteu * B_tau4;
    }
    derivs.dalphar_ddelta *= one_over_delta;
    derivs.dalphar_dtau *= one_over_tau;
    derivs.d2alphar_ddelta2 *= POW2(one_over_delta);
    derivs.d2alphar_dtau2 *= POW2(one_over_tau);
    derivs.d2alphar_ddelta_dtau *= one_over_delta * one_over_tau;

    derivs.d3alphar_ddelta3 *= POW3(one_over_delta);
    derivs.d3alphar_dtau3 *= POW3(one_over_tau);
    derivs.d3alphar_ddelta2_dtau *= POW2(one_over_delta) * one_over_tau;
    derivs.d3alphar_ddelta_dtau2 *= one_over_delta * POW2(one_over_tau);

    derivs.d4alphar_ddelta4 *= POW4(one_over_delta);
    derivs.d4alphar_dtau4 *= POW4(one_over_tau);
    derivs.d4alphar_ddelta3_dtau *= POW3(one_over_delta) * one_over_tau;
    derivs.d4alphar_ddelta2_dtau2 *= POW2(one_over_delta) * POW2(one_over_tau);
    derivs.d4alphar_ddelta_dtau3 *= one_over_delta * POW3(one_over_tau);

    return;
};

void ResidualHelmholtzGeneralizedExponential::to_json(rapidjson::Value& el, rapidjson::Document& doc) {
    el.AddMember("type", "GeneralizedExponential", doc.GetAllocator());
    cpjson::set_double_array("n", n, el, doc);
    cpjson::set_double_array("t", t, el, doc);
    cpjson::set_double_array("d", d, el, doc);
    cpjson::set_double_array("eta1", eta1, el, doc);
    cpjson::set_double_array("eta2", eta2, el, doc);
    cpjson::set_double_array("beta1", beta1, el, doc);
    cpjson::set_double_array("beta2", beta2, el, doc);
    cpjson::set_double_array("gamma1", gamma1, el, doc);
    cpjson::set_double_array("gamma2", gamma2, el, doc);
    cpjson::set_double_array("epsilon1", epsilon1, el, doc);
    cpjson::set_double_array("epsilon2", epsilon2, el, doc);
    cpjson::set_double_array("l_double", l_double, el, doc);
    cpjson::set_int_array("l_int", l_int, el, doc);
}

void ResidualHelmholtzNonAnalytic::to_json(rapidjson::Value& el, rapidjson::Document& doc) {
    el.AddMember("type", "ResidualHelmholtzNonAnalytic", doc.GetAllocator());

    rapidjson::Value _n(rapidjson::kArrayType), _a(rapidjson::kArrayType), _b(rapidjson::kArrayType), _beta(rapidjson::kArrayType),
      _A(rapidjson::kArrayType), _B(rapidjson::kArrayType), _C(rapidjson::kArrayType), _D(rapidjson::kArrayType);
    for (unsigned int i = 0; i <= N; ++i) {
        ResidualHelmholtzNonAnalyticElement& elem = elements[i];
        _n.PushBack((double)elem.n, doc.GetAllocator());
        _a.PushBack((double)elem.a, doc.GetAllocator());
        _b.PushBack((double)elem.b, doc.GetAllocator());
        _beta.PushBack((double)elem.beta, doc.GetAllocator());
        _A.PushBack((double)elem.A, doc.GetAllocator());
        _B.PushBack((double)elem.B, doc.GetAllocator());
        _C.PushBack((double)elem.C, doc.GetAllocator());
        _D.PushBack((double)elem.D, doc.GetAllocator());
    }
    el.AddMember("n", _n, doc.GetAllocator());
    el.AddMember("a", _a, doc.GetAllocator());
    el.AddMember("b", _b, doc.GetAllocator());
    el.AddMember("beta", _beta, doc.GetAllocator());
    el.AddMember("A", _A, doc.GetAllocator());
    el.AddMember("B", _B, doc.GetAllocator());
    el.AddMember("C", _C, doc.GetAllocator());
    el.AddMember("D", _D, doc.GetAllocator());
}

void ResidualHelmholtzNonAnalytic::all(const CoolPropDbl& tau_in, const CoolPropDbl& delta_in, HelmholtzDerivatives& derivs) throw() {
    if (N == 0) {
        return;
    }

    // Here we want to hack this function just a tiny bit to avoid evaluation AT the critical point
    // If we are VERY close to the critical point, just offset us a tiny bit away
    CoolPropDbl tau = tau_in, delta = delta_in;
    if (std::abs(tau_in - 1) < 10 * DBL_EPSILON) {
        tau = 1.0 + 10 * DBL_EPSILON;
    }
    if (std::abs(delta_in - 1) < 10 * DBL_EPSILON) {
        delta = 1.0 + 10 * DBL_EPSILON;
    }

    for (unsigned int i = 0; i < N; ++i) {
        const ResidualHelmholtzNonAnalyticElement& el = elements[i];
        const CoolPropDbl ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        const CoolPropDbl Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        // Derivatives of theta (all others are zero) (OK - checked)
        // Do not factor because then when delta = 1 you are dividing by 0
        const CoolPropDbl theta = (1.0 - tau) + Ai * pow(POW2(delta - 1.0), 1.0 / (2.0 * betai));
        const CoolPropDbl dtheta_dTau = -1;
        const CoolPropDbl dtheta_dDelta = Ai / (betai)*pow(POW2(delta - 1), 1 / (2 * betai) - 1) * (delta - 1);

        const CoolPropDbl d2theta_dDelta2 = Ai / betai * (1 / betai - 1) * pow(POW2(delta - 1), 1 / (2 * betai) - 1);
        const CoolPropDbl d3theta_dDelta3 = Ai / betai * (2 - 3 / betai + 1 / POW2(betai)) * pow(POW2(delta - 1), 1 / (2 * betai)) / POW3(delta - 1);
        const CoolPropDbl d4theta_dDelta4 =
          Ai / betai * (-6 + 11 / betai - 6 / POW2(betai) + 1 / POW3(betai)) * pow(POW2(delta - 1), 1 / (2 * betai) - 2);

        // Derivatives of PSI (OK - checked)
        const CoolPropDbl PSI = exp(-Ci * POW2(delta - 1.0) - Di * POW2(tau - 1.0));
        const CoolPropDbl dPSI_dDelta_over_PSI = -2.0 * Ci * (delta - 1.0);
        const CoolPropDbl dPSI_dDelta = dPSI_dDelta_over_PSI * PSI;
        const CoolPropDbl dPSI_dTau_over_PSI = -2.0 * Di * (tau - 1.0);
        const CoolPropDbl dPSI_dTau = dPSI_dTau_over_PSI * PSI;
        const CoolPropDbl d2PSI_dDelta2_over_PSI = (2.0 * Ci * POW2(delta - 1.0) - 1.0) * 2.0 * Ci;
        const CoolPropDbl d2PSI_dDelta2 = d2PSI_dDelta2_over_PSI * PSI;
        const CoolPropDbl d3PSI_dDelta3 = 2 * Ci * PSI * (-4 * Ci * Ci * POW3(delta - 1) + 6 * Ci * (delta - 1));
        const CoolPropDbl d4PSI_dDelta4 = 4 * Ci * Ci * PSI * (4 * Ci * Ci * POW4(delta - 1) - 12 * Ci * POW2(delta - 1) + 3);
        const CoolPropDbl d2PSI_dTau2 = (2.0 * Di * POW2(tau - 1.0) - 1.0) * 2.0 * Di * PSI;
        const CoolPropDbl d3PSI_dTau3 = 2.0 * Di * PSI * (-4 * Di * Di * POW3(tau - 1) + 6 * Di * (tau - 1));
        const CoolPropDbl d4PSI_dTau4 = 4 * Di * Di * PSI * (4 * Di * Di * POW4(tau - 1) - 12 * Di * POW2(tau - 1) + 3);
        const CoolPropDbl d2PSI_dDelta_dTau = dPSI_dDelta * dPSI_dTau_over_PSI;
        const CoolPropDbl d3PSI_dDelta2_dTau = d2PSI_dDelta2 * dPSI_dTau_over_PSI;
        const CoolPropDbl d3PSI_dDelta_dTau2 = d2PSI_dTau2 * dPSI_dDelta_over_PSI;
        const CoolPropDbl d4PSI_dDelta_dTau3 = d3PSI_dTau3 * dPSI_dDelta_over_PSI;
        const CoolPropDbl d4PSI_dDelta2_dTau2 = d2PSI_dTau2 * d2PSI_dDelta2_over_PSI;
        const CoolPropDbl d4PSI_dDelta3_dTau = d3PSI_dDelta3 * dPSI_dTau_over_PSI;

        // Derivatives of DELTA (OK - Checked)
        const CoolPropDbl DELTA = POW2(theta) + Bi * pow(POW2(delta - 1.0), ai);
        const CoolPropDbl dDELTA_dTau = 2 * theta * dtheta_dTau;
        const CoolPropDbl dDELTA_dDelta = 2 * theta * dtheta_dDelta + 2 * Bi * ai * pow(POW2(delta - 1.0), ai - 1.0) * (delta - 1);
        const CoolPropDbl d2DELTA_dTau2 = 2;                                      // d2theta_dTau2 is zero and (dtheta_dtau)^2 = 1
        const CoolPropDbl d2DELTA_dDelta_dTau = 2 * dtheta_dTau * dtheta_dDelta;  // d2theta_dDelta2 is zero
        const CoolPropDbl d2DELTA_dDelta2 =
          2 * (theta * d2theta_dDelta2 + POW2(dtheta_dDelta) + Bi * (2 * ai * ai - ai) * pow(POW2(delta - 1.0), ai - 1.0));
        const CoolPropDbl d3DELTA_dTau3 = 0;
        const CoolPropDbl d3DELTA_dDelta_dTau2 = 0;
        const CoolPropDbl d3DELTA_dDelta2_dTau = 2 * dtheta_dTau * d2theta_dDelta2;
        const CoolPropDbl d3DELTA_dDelta3 = 2
                                            * (theta * d3theta_dDelta3 + 3 * dtheta_dDelta * d2theta_dDelta2
                                               + 2 * Bi * ai * (2 * ai * ai - 3 * ai + 1) * pow(POW2(delta - 1.0), ai - 1.0) / (delta - 1));

        const CoolPropDbl d4DELTA_dTau4 = 0;
        const CoolPropDbl d4DELTA_dDelta_dTau3 = 0;
        const CoolPropDbl d4DELTA_dDelta2_dTau2 = 0;
        const CoolPropDbl d4DELTA_dDelta3_dTau = 2 * dtheta_dTau * d3theta_dDelta3;
        const CoolPropDbl d4DELTA_dDelta4 = 2
                                            * (theta * d4theta_dDelta4 + 4 * dtheta_dDelta * d3theta_dDelta3 + 3 * POW2(d2theta_dDelta2)
                                               + 2 * Bi * ai * (4 * ai * ai * ai - 12 * ai * ai + 11 * ai - 3) * pow(POW2(delta - 1.0), ai - 2.0));

        const CoolPropDbl dDELTAbi_dDelta = bi * pow(DELTA, bi - 1.0) * dDELTA_dDelta;
        const CoolPropDbl dDELTAbi_dTau = -2.0 * theta * bi * pow(DELTA, bi - 1.0);
        const CoolPropDbl d2DELTAbi_dDelta2 = bi * (pow(DELTA, bi - 1) * d2DELTA_dDelta2 + (bi - 1.0) * pow(DELTA, bi - 2.0) * pow(dDELTA_dDelta, 2));
        const CoolPropDbl d3DELTAbi_dDelta3 =
          bi
          * (pow(DELTA, bi - 1) * d3DELTA_dDelta3 + d2DELTA_dDelta2 * (bi - 1) * pow(DELTA, bi - 2) * dDELTA_dDelta
             + (bi - 1)
                 * (pow(DELTA, bi - 2) * 2 * dDELTA_dDelta * d2DELTA_dDelta2
                    + pow(dDELTA_dDelta, 2) * (bi - 2) * pow(DELTA, bi - 3) * dDELTA_dDelta));
        const CoolPropDbl d2DELTAbi_dDelta_dTau =
          -Ai * bi * 2.0 / betai * pow(DELTA, bi - 1.0) * (delta - 1.0) * pow(pow(delta - 1.0, 2), 1.0 / (2.0 * betai) - 1.0)
          - 2.0 * theta * bi * (bi - 1.0) * pow(DELTA, bi - 2.0) * dDELTA_dDelta;
        const CoolPropDbl d2DELTAbi_dTau2 = 2.0 * bi * pow(DELTA, bi - 1.0) + 4.0 * pow(theta, 2) * bi * (bi - 1.0) * pow(DELTA, bi - 2.0);
        const CoolPropDbl d3DELTAbi_dTau3 =
          -12.0 * theta * bi * (bi - 1.0) * pow(DELTA, bi - 2) - 8 * pow(theta, 3) * bi * (bi - 1) * (bi - 2) * pow(DELTA, bi - 3);
        const CoolPropDbl d3DELTAbi_dDelta_dTau2 = 2 * bi * (bi - 1) * pow(DELTA, bi - 2) * dDELTA_dDelta
                                                   + 4 * pow(theta, 2) * bi * (bi - 1) * (bi - 2) * pow(DELTA, bi - 3) * dDELTA_dDelta
                                                   + 8 * theta * bi * (bi - 1) * pow(DELTA, bi - 2) * dtheta_dDelta;
        const CoolPropDbl d3DELTAbi_dDelta2_dTau =
          bi
          * ((bi - 1) * pow(DELTA, bi - 2) * dDELTA_dTau * d2DELTA_dDelta2 + pow(DELTA, bi - 1) * d3DELTA_dDelta2_dTau
             + (bi - 1)
                 * ((bi - 2) * pow(DELTA, bi - 3) * dDELTA_dTau * pow(dDELTA_dDelta, 2)
                    + pow(DELTA, bi - 2) * 2 * dDELTA_dDelta * d2DELTA_dDelta_dTau));

        // Fourth partials
        const CoolPropDbl DELTA_bi = pow(DELTA, bi);
        const CoolPropDbl d4DELTAbi_dTau4 =
          bi * DELTA_bi / DELTA
          * ((POW3(bi) - 6 * POW2(bi) + 11 * bi - 6) * POW4(dDELTA_dTau) / POW3(DELTA)
             + 6 * (bi * bi - 3 * bi + 2) * POW2(dDELTA_dTau / DELTA) * d2DELTA_dTau2 + 4 * (bi - 1) * dDELTA_dTau / DELTA * d3DELTA_dTau3
             + 3 * (bi - 1) * POW2(d2DELTA_dTau2) / DELTA + d4DELTA_dTau4);
        const CoolPropDbl d4DELTAbi_dDelta4 =
          bi * DELTA_bi / DELTA
          * ((POW3(bi) - 6 * POW2(bi) + 11 * bi - 6) * POW4(dDELTA_dDelta) / POW3(DELTA)
             + 6 * (bi * bi - 3 * bi + 2) * POW2(dDELTA_dDelta / DELTA) * d2DELTA_dDelta2 + 4 * (bi - 1) * dDELTA_dDelta / DELTA * d3DELTA_dDelta3
             + 3 * (bi - 1) * POW2(d2DELTA_dDelta2) / DELTA + d4DELTA_dDelta4);
        const CoolPropDbl d4DELTAbi_dDelta_dTau3 =
          bi * (bi - 1) * DELTA_bi / POW2(DELTA) * dDELTA_dDelta
            * ((bi - 1) * (bi - 2) * POW3(dDELTA_dTau) / POW2(DELTA) + 3 * (bi - 1) * dDELTA_dTau / DELTA * d2DELTA_dTau2 + d3DELTA_dTau3)
          + bi * DELTA_bi / DELTA
              * ((bi - 1) * (bi - 2)
                   * (3 * POW2(dDELTA_dTau) / POW2(DELTA) * d2DELTA_dDelta_dTau - 2 * POW3(dDELTA_dTau) / POW3(DELTA) * dDELTA_dDelta)
                 + 3 * (bi - 1)
                     * (dDELTA_dTau / DELTA * d3DELTA_dDelta_dTau2 + d2DELTA_dTau2 / DELTA * d2DELTA_dDelta_dTau
                        - dDELTA_dDelta / POW2(DELTA) * dDELTA_dTau * d2DELTA_dTau2)
                 + d4DELTA_dDelta_dTau3);
        const CoolPropDbl d4DELTAbi_dDelta3_dTau =
          bi * (bi - 1) * DELTA_bi / POW2(DELTA) * dDELTA_dTau
            * ((bi - 1) * (bi - 2) * POW3(dDELTA_dDelta) / POW2(DELTA) + 3 * (bi - 1) * dDELTA_dDelta / DELTA * d2DELTA_dDelta2 + d3DELTA_dDelta3)
          + bi * DELTA_bi / DELTA
              * ((bi - 1) * (bi - 2)
                   * (3 * POW2(dDELTA_dDelta) / POW2(DELTA) * d2DELTA_dDelta_dTau - 2 * POW3(dDELTA_dDelta) / POW3(DELTA) * dDELTA_dTau)
                 + 3 * (bi - 1)
                     * (dDELTA_dDelta / DELTA * d3DELTA_dDelta2_dTau + d2DELTA_dDelta2 / DELTA * d2DELTA_dDelta_dTau
                        - dDELTA_dTau / POW2(DELTA) * dDELTA_dDelta * d2DELTA_dDelta2)
                 + d4DELTA_dDelta3_dTau);
        const CoolPropDbl d4DELTAbi_dDelta2_dTau2 = bi * DELTA_bi / POW4(DELTA)
                                                    * ((POW3(bi) - 6 * bi * bi + 11 * bi - 6) * POW2(dDELTA_dDelta) * POW2(dDELTA_dTau)  // Yellow
                                                       + (bi - 1) * (bi - 2) * DELTA * POW2(dDELTA_dDelta) * d2DELTA_dTau2               // Orange
                                                       + (bi - 1) * (bi - 2) * DELTA * POW2(dDELTA_dTau) * d2DELTA_dDelta2               // Pink
                                                       + 4 * (bi - 1) * (bi - 2) * DELTA * dDELTA_dDelta * dDELTA_dTau * d2DELTA_dDelta_dTau  // Green
                                                       + 2 * (bi - 1) * POW2(DELTA * d2DELTA_dDelta_dTau)                   // Blue hi
                                                       + 2 * (bi - 1) * POW2(DELTA) * dDELTA_dTau * d3DELTA_dDelta2_dTau    // Blue sharp
                                                       + 2 * (bi - 1) * POW2(DELTA) * dDELTA_dDelta * d3DELTA_dDelta_dTau2  // Red sharp
                                                       + (bi - 1) * POW2(DELTA) * d2DELTA_dDelta2 * d2DELTA_dTau2           // black sharp
                                                       + POW3(DELTA) * d4DELTA_dDelta2_dTau2);

        derivs.alphar += delta * ni * DELTA_bi * PSI;

        // First partials
        derivs.dalphar_dtau += ni * delta * (DELTA_bi * dPSI_dTau + dDELTAbi_dTau * PSI);
        derivs.dalphar_ddelta += ni * (DELTA_bi * (PSI + delta * dPSI_dDelta) + dDELTAbi_dDelta * delta * PSI);

        // Second partials
        derivs.d2alphar_dtau2 += ni * delta * (d2DELTAbi_dTau2 * PSI + 2 * dDELTAbi_dTau * dPSI_dTau + DELTA_bi * d2PSI_dTau2);
        derivs.d2alphar_ddelta_dtau += ni
                                       * (DELTA_bi * (dPSI_dTau + delta * d2PSI_dDelta_dTau) + delta * dDELTAbi_dDelta * dPSI_dTau
                                          + dDELTAbi_dTau * (PSI + delta * dPSI_dDelta) + d2DELTAbi_dDelta_dTau * delta * PSI);
        derivs.d2alphar_ddelta2 += ni
                                   * (DELTA_bi * (2.0 * dPSI_dDelta + delta * d2PSI_dDelta2) + 2.0 * dDELTAbi_dDelta * (PSI + delta * dPSI_dDelta)
                                      + d2DELTAbi_dDelta2 * delta * PSI);

        // Third partials
        derivs.d3alphar_dtau3 +=
          ni * delta * (d3DELTAbi_dTau3 * PSI + 3 * d2DELTAbi_dTau2 * dPSI_dTau + 3 * dDELTAbi_dTau * d2PSI_dTau2 + DELTA_bi * d3PSI_dTau3);
        derivs.d3alphar_ddelta_dtau2 +=
          ni * delta
            * (d2DELTAbi_dTau2 * dPSI_dDelta + d3DELTAbi_dDelta_dTau2 * PSI + 2 * dDELTAbi_dTau * d2PSI_dDelta_dTau
               + 2.0 * d2DELTAbi_dDelta_dTau * dPSI_dTau + DELTA_bi * d3PSI_dDelta_dTau2 + dDELTAbi_dDelta * d2PSI_dTau2)
          + ni * (d2DELTAbi_dTau2 * PSI + 2.0 * dDELTAbi_dTau * dPSI_dTau + DELTA_bi * d2PSI_dTau2);
        derivs.d3alphar_ddelta3 +=
          ni
          * (DELTA_bi * (3 * d2PSI_dDelta2 + delta * d3PSI_dDelta3) + 3 * dDELTAbi_dDelta * (2 * dPSI_dDelta + delta * d2PSI_dDelta2)
             + 3 * d2DELTAbi_dDelta2 * (PSI + delta * dPSI_dDelta) + d3DELTAbi_dDelta3 * PSI * delta);
        CoolPropDbl Line1 =
          DELTA_bi * (2 * d2PSI_dDelta_dTau + delta * d3PSI_dDelta2_dTau) + dDELTAbi_dTau * (2 * dPSI_dDelta + delta * d2PSI_dDelta2);
        CoolPropDbl Line2 = 2 * dDELTAbi_dDelta * (dPSI_dTau + delta * d2PSI_dDelta_dTau) + 2 * d2DELTAbi_dDelta_dTau * (PSI + delta * dPSI_dDelta);
        CoolPropDbl Line3 = d2DELTAbi_dDelta2 * delta * dPSI_dTau + d3DELTAbi_dDelta2_dTau * delta * PSI;
        derivs.d3alphar_ddelta2_dtau += ni * (Line1 + Line2 + Line3);

        // Fourth partials
        derivs.d4alphar_dtau4 += ni * delta
                                 * (DELTA_bi * d4PSI_dTau4 + 4 * dDELTAbi_dTau * d3PSI_dTau3 + 6 * d2DELTAbi_dTau2 * d2PSI_dTau2
                                    + 4 * d3DELTAbi_dTau3 * dPSI_dTau + PSI * d4DELTAbi_dTau4);
        derivs.d4alphar_ddelta4 +=
          ni
          * (delta * DELTA_bi * d4PSI_dDelta4 + delta * PSI * d4DELTAbi_dDelta4 + 4 * delta * dDELTAbi_dDelta * d3PSI_dDelta3
             + 4 * delta * dPSI_dDelta * d3DELTAbi_dDelta3 + 6 * delta * d2DELTAbi_dDelta2 * d2PSI_dDelta2 + 4 * DELTA_bi * d3PSI_dDelta3
             + 4 * PSI * d3DELTAbi_dDelta3 + 12 * dDELTAbi_dDelta * d2PSI_dDelta2 + 12 * dPSI_dDelta * d2DELTAbi_dDelta2);

        derivs.d4alphar_ddelta_dtau3 +=
          ni
          * (delta * DELTA_bi * d4PSI_dDelta_dTau3 + delta * PSI * d4DELTAbi_dDelta_dTau3 + delta * dDELTAbi_dDelta * d3PSI_dTau3
             + 3 * delta * dDELTAbi_dTau * d3PSI_dDelta_dTau2 + delta * dPSI_dDelta * d3DELTAbi_dTau3 + 3 * delta * dPSI_dTau * d3DELTAbi_dDelta_dTau2
             + 3 * delta * d2DELTAbi_dDelta_dTau * d2PSI_dTau2 + 3 * delta * d2DELTAbi_dTau2 * d2PSI_dDelta_dTau + DELTA_bi * d3PSI_dTau3
             + PSI * d3DELTAbi_dTau3 + 3 * dDELTAbi_dTau * d2PSI_dTau2 + 3 * dPSI_dTau * d2DELTAbi_dTau2);
        derivs.d4alphar_ddelta3_dtau +=
          ni
          * (delta * DELTA_bi * d4PSI_dDelta3_dTau + delta * PSI * d4DELTAbi_dDelta3_dTau + 3 * delta * dDELTAbi_dDelta * d3PSI_dDelta2_dTau
             + delta * dDELTAbi_dTau * d3PSI_dDelta3 + 3 * delta * dPSI_dDelta * d3DELTAbi_dDelta2_dTau + delta * dPSI_dTau * d3DELTAbi_dDelta3
             + 3 * delta * d2DELTAbi_dDelta2 * d2PSI_dDelta_dTau + 3 * delta * d2DELTAbi_dDelta_dTau * d2PSI_dDelta2
             + 3 * DELTA_bi * d3PSI_dDelta2_dTau + 3 * PSI * d3DELTAbi_dDelta2_dTau + 6 * dDELTAbi_dDelta * d2PSI_dDelta_dTau
             + 3 * dDELTAbi_dTau * d2PSI_dDelta2 + 6 * dPSI_dDelta * d2DELTAbi_dDelta_dTau + 3 * dPSI_dTau * d2DELTAbi_dDelta2);
        derivs.d4alphar_ddelta2_dtau2 +=
          ni
          * (delta * DELTA_bi * d4PSI_dDelta2_dTau2 + delta * PSI * d4DELTAbi_dDelta2_dTau2 + 2 * delta * dDELTAbi_dDelta * d3PSI_dDelta_dTau2
             + 2 * delta * dDELTAbi_dTau * d3PSI_dDelta2_dTau + 2 * delta * dPSI_dDelta * d3DELTAbi_dDelta_dTau2
             + 2 * delta * dPSI_dTau * d3DELTAbi_dDelta2_dTau + delta * d2DELTAbi_dDelta2 * d2PSI_dTau2
             + 4 * delta * d2DELTAbi_dDelta_dTau * d2PSI_dDelta_dTau + delta * d2DELTAbi_dTau2 * d2PSI_dDelta2 + 2 * DELTA_bi * d3PSI_dDelta_dTau2
             + 2 * PSI * d3DELTAbi_dDelta_dTau2 + 2 * dDELTAbi_dDelta * d2PSI_dTau2 + 4 * dDELTAbi_dTau * d2PSI_dDelta_dTau
             + 2 * dPSI_dDelta * d2DELTAbi_dTau2 + 4 * dPSI_dTau * d2DELTAbi_dDelta_dTau);
    }
}

void ResidualHelmholtzGeneralizedCubic::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }

    derivs.alphar += m_abstractcubic->alphar(tau, delta, z, 0, 0);

    derivs.dalphar_ddelta += m_abstractcubic->alphar(tau, delta, z, 0, 1);
    derivs.dalphar_dtau += m_abstractcubic->alphar(tau, delta, z, 1, 0);

    derivs.d2alphar_ddelta2 += m_abstractcubic->alphar(tau, delta, z, 0, 2);
    derivs.d2alphar_ddelta_dtau += m_abstractcubic->alphar(tau, delta, z, 1, 1);
    derivs.d2alphar_dtau2 += m_abstractcubic->alphar(tau, delta, z, 2, 0);

    derivs.d3alphar_ddelta3 += m_abstractcubic->alphar(tau, delta, z, 0, 3);
    derivs.d3alphar_ddelta2_dtau += m_abstractcubic->alphar(tau, delta, z, 1, 2);
    derivs.d3alphar_ddelta_dtau2 += m_abstractcubic->alphar(tau, delta, z, 2, 1);
    derivs.d3alphar_dtau3 += m_abstractcubic->alphar(tau, delta, z, 3, 0);

    derivs.d4alphar_ddelta4 += m_abstractcubic->alphar(tau, delta, z, 0, 4);
    derivs.d4alphar_ddelta3_dtau += m_abstractcubic->alphar(tau, delta, z, 1, 3);
    derivs.d4alphar_ddelta2_dtau2 += m_abstractcubic->alphar(tau, delta, z, 2, 2);
    derivs.d4alphar_ddelta_dtau3 += m_abstractcubic->alphar(tau, delta, z, 3, 1);
    derivs.d4alphar_dtau4 += m_abstractcubic->alphar(tau, delta, z, 4, 0);
}

/**

Sympy code:

import sympy as sy
n,t,d,eta,beta,gamma,epsilon,b,tau,delta = sy.symbols('n,t,d,eta,beta,gamma,epsilon,b,tau,delta')
Fdelta = delta**d*sy.exp(eta*(delta-epsilon)**2)
Ftau = tau**t*sy.exp(1/(beta*(tau-gamma)**2+b))

alphar = n*Ftau*Fdelta
display(sy.ccode(Ftau))
display(sy.ccode(Fdelta))

display(sy.ccode(sy.simplify(sy.diff(Fdelta, delta, 1)*delta**1)))
display(sy.ccode(sy.simplify(sy.diff(Fdelta, delta, 2)*delta**2)))
display(sy.ccode(sy.simplify(sy.diff(Fdelta, delta, 3)*delta**3)))
display(sy.ccode(sy.simplify(sy.diff(Fdelta, delta, 4)*delta**4)))

display(sy.ccode(sy.simplify(sy.diff(Ftau, tau, 1)*tau**1)))
display(sy.ccode(sy.simplify(sy.diff(Ftau, tau, 2)*tau**2)))
display(sy.ccode(sy.simplify(sy.diff(Ftau, tau, 3)*tau**3)))
display(sy.ccode(sy.simplify(sy.diff(Ftau, tau, 4)*tau**4)))

*/
void ResidualHelmholtzGaoB::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }

    CoolPropDbl Ftau = 0, Fdelta = 0, taudFtaudtau = 0, tau2d2Ftaudtau2 = 0, tau3d3Ftaudtau3 = 0, tau4d4Ftaudtau4 = 0, deltadFdeltaddelta = 0,
                delta2d2Fdeltaddelta2 = 0, delta3d3Fdeltaddelta3 = 0, delta4d4Fdeltaddelta4 = 0;

    for (int i = 0; i < static_cast<int>(n.size()); ++i) {

        const CoolPropDbl n = this->n[i], t = this->t[i], d = this->d[i], eta = this->eta[i], beta = this->beta[i], gamma = this->gamma[i],
                          epsilon = this->epsilon[i], b = this->b[i];

        Ftau = pow(tau, t) * exp(1.0 / (b + beta * pow(-gamma + tau, 2)));
        Fdelta = pow(delta, d) * exp(eta * pow(delta - epsilon, 2));
        taudFtaudtau = (2 * beta * pow(tau, t + 1) * (gamma - tau) + t * pow(tau, t) * pow(b + beta * pow(gamma - tau, 2), 2))
                       * exp(1.0 / (b + beta * pow(gamma - tau, 2))) / pow(b + beta * pow(gamma - tau, 2), 2);
        tau2d2Ftaudtau2 = pow(tau, t)
                          * (4 * beta * t * tau * pow(b + beta * pow(gamma - tau, 2), 2) * (gamma - tau)
                             + 2 * beta * pow(tau, 2)
                                 * (4 * beta * (b + beta * pow(gamma - tau, 2)) * pow(gamma - tau, 2) + 2 * beta * pow(gamma - tau, 2)
                                    - pow(b + beta * pow(gamma - tau, 2), 2))
                             + t * pow(b + beta * pow(gamma - tau, 2), 4) * (t - 1))
                          * exp(1.0 / (b + beta * pow(gamma - tau, 2))) / pow(b + beta * pow(gamma - tau, 2), 4);
        tau3d3Ftaudtau3 =
          pow(tau, t)
          * (4 * pow(beta, 2) * pow(tau, 3) * (gamma - tau)
               * (12 * beta * (b + beta * pow(gamma - tau, 2)) * pow(gamma - tau, 2) + 2 * beta * pow(gamma - tau, 2)
                  - 6 * pow(b + beta * pow(gamma - tau, 2), 3) + pow(b + beta * pow(gamma - tau, 2), 2) * (12 * beta * pow(gamma - tau, 2) - 3))
             + 6 * beta * t * pow(tau, 2) * pow(b + beta * pow(gamma - tau, 2), 2)
                 * (4 * beta * (b + beta * pow(gamma - tau, 2)) * pow(gamma - tau, 2) + 2 * beta * pow(gamma - tau, 2)
                    - pow(b + beta * pow(gamma - tau, 2), 2))
             + 6 * beta * t * tau * pow(b + beta * pow(gamma - tau, 2), 4) * (gamma - tau) * (t - 1)
             + t * pow(b + beta * pow(gamma - tau, 2), 6) * (pow(t, 2) - 3 * t + 2))
          * exp(1.0 / (b + beta * pow(gamma - tau, 2))) / pow(b + beta * pow(gamma - tau, 2), 6);
        tau4d4Ftaudtau4 =
          pow(tau, t)
          * (16 * pow(beta, 2) * t * pow(tau, 3) * pow(b + beta * pow(gamma - tau, 2), 2) * (gamma - tau)
               * (12 * beta * (b + beta * pow(gamma - tau, 2)) * pow(gamma - tau, 2) + 2 * beta * pow(gamma - tau, 2)
                  - 6 * pow(b + beta * pow(gamma - tau, 2), 3) + pow(b + beta * pow(gamma - tau, 2), 2) * (12 * beta * pow(gamma - tau, 2) - 3))
             + pow(beta, 2) * pow(tau, 4)
                 * (pow(beta, 2) * (192 * b + 192 * beta * pow(gamma - tau, 2)) * pow(gamma - tau, 4) + 16 * pow(beta, 2) * pow(gamma - tau, 4)
                    + 96 * beta * pow(b + beta * pow(gamma - tau, 2), 3) * pow(gamma - tau, 2) * (4 * beta * pow(gamma - tau, 2) - 3)
                    + 48 * beta * pow(b + beta * pow(gamma - tau, 2), 2) * pow(gamma - tau, 2) * (12 * beta * pow(gamma - tau, 2) - 1)
                    + 24 * pow(b + beta * pow(gamma - tau, 2), 5) + pow(b + beta * pow(gamma - tau, 2), 4) * (-288 * beta * pow(gamma - tau, 2) + 12))
             + 12 * beta * t * pow(tau, 2) * pow(b + beta * pow(gamma - tau, 2), 4) * (t - 1)
                 * (4 * beta * (b + beta * pow(gamma - tau, 2)) * pow(gamma - tau, 2) + 2 * beta * pow(gamma - tau, 2)
                    - pow(b + beta * pow(gamma - tau, 2), 2))
             + 8 * beta * t * tau * pow(b + beta * pow(gamma - tau, 2), 6) * (gamma - tau) * (pow(t, 2) - 3 * t + 2)
             + t * pow(b + beta * pow(gamma - tau, 2), 8) * (pow(t, 3) - 6 * pow(t, 2) + 11 * t - 6))
          * exp(1.0 / (b + beta * pow(gamma - tau, 2))) / pow(b + beta * pow(gamma - tau, 2), 8);
        deltadFdeltaddelta = (d * pow(delta, d) + 2 * pow(delta, d + 1) * eta * (delta - epsilon)) * exp(eta * pow(delta - epsilon, 2));
        delta2d2Fdeltaddelta2 =
          pow(delta, d) * (4 * d * delta * eta * (delta - epsilon) + d * (d - 1) + 2 * pow(delta, 2) * eta * (2 * eta * pow(delta - epsilon, 2) + 1))
          * exp(eta * pow(delta - epsilon, 2));
        delta3d3Fdeltaddelta3 =
          pow(delta, d)
          * (6 * d * pow(delta, 2) * eta * (2 * eta * pow(delta - epsilon, 2) + 1) + 6 * d * delta * eta * (d - 1) * (delta - epsilon)
             + d * (pow(d, 2) - 3 * d + 2) + 4 * pow(delta, 3) * pow(eta, 2) * (delta - epsilon) * (2 * eta * pow(delta - epsilon, 2) + 3))
          * exp(eta * pow(delta - epsilon, 2));
        delta4d4Fdeltaddelta4 =
          pow(delta, d)
          * (16 * d * pow(delta, 3) * pow(eta, 2) * (delta - epsilon) * (2 * eta * pow(delta - epsilon, 2) + 3)
             + 12 * d * pow(delta, 2) * eta * (d - 1) * (2 * eta * pow(delta - epsilon, 2) + 1)
             + 8 * d * delta * eta * (delta - epsilon) * (pow(d, 2) - 3 * d + 2) + d * (pow(d, 3) - 6 * pow(d, 2) + 11 * d - 6)
             + pow(delta, 4) * pow(eta, 2) * (16 * pow(eta, 2) * pow(delta - epsilon, 4) + 48 * eta * pow(delta - epsilon, 2) + 12))
          * exp(eta * pow(delta - epsilon, 2));

        derivs.alphar += n * Ftau * Fdelta;

        derivs.dalphar_ddelta += n * Ftau * deltadFdeltaddelta / delta;
        derivs.dalphar_dtau += n * Fdelta * taudFtaudtau / tau;

        derivs.d2alphar_ddelta2 += n * Ftau * delta2d2Fdeltaddelta2 / POW2(delta);
        derivs.d2alphar_ddelta_dtau += n * taudFtaudtau * deltadFdeltaddelta / tau / delta;
        derivs.d2alphar_dtau2 += n * Fdelta * tau2d2Ftaudtau2 / POW2(tau);

        derivs.d3alphar_ddelta3 += n * Ftau * delta3d3Fdeltaddelta3 / POW3(delta);
        derivs.d3alphar_ddelta2_dtau += n * taudFtaudtau * delta2d2Fdeltaddelta2 / POW2(delta) / tau;
        derivs.d3alphar_ddelta_dtau2 += n * tau2d2Ftaudtau2 * deltadFdeltaddelta / POW2(tau) / delta;
        derivs.d3alphar_dtau3 += n * Fdelta * tau3d3Ftaudtau3 / POW3(tau);

        derivs.d4alphar_ddelta4 += n * Ftau * delta4d4Fdeltaddelta4 / POW4(delta);
        derivs.d4alphar_ddelta3_dtau += n * taudFtaudtau * delta3d3Fdeltaddelta3 / POW3(delta) / tau;
        derivs.d4alphar_ddelta2_dtau2 += n * tau2d2Ftaudtau2 * delta2d2Fdeltaddelta2 / POW2(delta) / POW2(tau);
        derivs.d4alphar_ddelta_dtau3 += n * tau3d3Ftaudtau3 * deltadFdeltaddelta / POW3(tau) / delta;
        derivs.d4alphar_dtau4 += n * Fdelta * tau4d4Ftaudtau4 / POW4(tau);
    }
}

ResidualHelmholtzXiangDeiters::ResidualHelmholtzXiangDeiters(const CoolPropDbl Tc, const CoolPropDbl pc, const CoolPropDbl rhomolarc,
                                                             const CoolPropDbl acentric, const CoolPropDbl R)
  : Tc(Tc), pc(pc), rhomolarc(rhomolarc), acentric(acentric), R(R) {
    double Zc = pc / (R * Tc * rhomolarc);
    theta = POW2(Zc - 0.29);

    // From Xiang & Deiters, doi:10.1016/j.ces.2007.11.029
    double _d[] = {1, 1, 1, 2, 3, 7, 1, 1, 2, 5, 1, 1, 4, 2};
    std::vector<CoolPropDbl> d(_d, _d + sizeof(_d) / sizeof(double));
    double _t[] = {0.25, 1.25, 1.5, 1.375, 0.25, 0.875, 0, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5};
    std::vector<CoolPropDbl> t(_t, _t + sizeof(_t) / sizeof(double));
    double _l[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3};
    std::vector<CoolPropDbl> l(_l, _l + sizeof(_l) / sizeof(double));
    double _g[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1};
    std::vector<CoolPropDbl> g(_g, _g + sizeof(_g) / sizeof(double));
    double _a0[] = {8.5740489E-01,  -3.2863233E+00, 1.6480939E+00,  -5.4524817E-02, 6.1623592E-02, 2.7389266E-04,  -6.0655087E-02,
                    -3.1811852E-02, -1.1550422E-01, -1.8610466E-02, -1.8348671E-01, 5.5071325E-03, -1.2268039E-02, -5.0433436E-03};
    std::vector<CoolPropDbl> a0(_a0, _a0 + sizeof(_a0) / sizeof(double));
    double _a1[] = {5.6200117E-01, 3.2439544E+00, -4.9628768E+00, -2.2132851E-01, 9.3388356E-02, 2.4940171E-05,  -1.7519204E-01,
                    8.9325660E-01, 2.9886613E+00, 1.0881387E-01,  -6.7166746E-01, 1.4477326E-01, -2.8716809E-01, -1.1478402E-01};
    std::vector<CoolPropDbl> a1(_a1, _a1 + sizeof(_a1) / sizeof(double));
    double _a2[] = {-8.1680511E+01, 4.6384732E+02,  -2.7970850E+02, 2.9317364E+01, -2.2324825E+01, -5.0932691E-02, -7.2836590E+00,
                    -2.2063100E+02, -3.0435126E+02, 5.8514719E+00,  1.7995451E+02, -1.0178400E+02, 4.0848053E+01,  1.2411984E+01};
    std::vector<CoolPropDbl> a2(_a2, _a2 + sizeof(_a2) / sizeof(double));

    phi0.add_Exponential(a0, d, t, g, l);
    phi1.add_Exponential(a1, d, t, g, l);
    phi2.add_Exponential(a2, d, t, g, l);

    enabled = true;
};

void ResidualHelmholtzXiangDeiters::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }

    HelmholtzDerivatives derivs0, derivs1, derivs2;

    // Calculate each of the derivative terms
    phi0.all(tau, delta, derivs0);
    phi1.all(tau, delta, derivs1);
    phi2.all(tau, delta, derivs2);

    // Add up the contributions
    derivs = derivs + derivs0 + derivs1 * acentric + derivs2 * theta;
}

void ResidualHelmholtzSAFTAssociating::to_json(rapidjson::Value& el, rapidjson::Document& doc) {
    el.AddMember("type", "ResidualHelmholtzSAFTAssociating", doc.GetAllocator());
    el.AddMember("a", a, doc.GetAllocator());
    el.AddMember("m", m, doc.GetAllocator());
    el.AddMember("epsilonbar", epsilonbar, doc.GetAllocator());
    el.AddMember("vbarn", vbarn, doc.GetAllocator());
    el.AddMember("kappabar", kappabar, doc.GetAllocator());
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::Deltabar(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->g(this->eta(delta)) * (exp(this->epsilonbar * tau) - 1) * this->kappabar;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dDeltabar_ddelta__consttau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->dg_deta(this->eta(delta)) * (exp(this->epsilonbar * tau) - 1) * this->kappabar * this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta2__consttau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->d2g_deta2(this->eta(delta)) * (exp(this->epsilonbar * tau) - 1) * this->kappabar * pow(this->vbarn, (int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dDeltabar_dtau__constdelta(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->g(this->eta(delta)) * this->kappabar * exp(this->epsilonbar * tau) * this->epsilonbar;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_dtau2__constdelta(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->g(this->eta(delta)) * this->kappabar * exp(this->epsilonbar * tau) * pow(this->epsilonbar, (int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta_dtau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->dg_deta(this->eta(delta)) * exp(this->epsilonbar * tau) * this->epsilonbar * this->kappabar * this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_dtau3__constdelta(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->g(this->eta(delta)) * this->kappabar * exp(this->epsilonbar * tau) * pow(this->epsilonbar, (int)3);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta_dtau2(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->dg_deta(this->eta(delta)) * this->kappabar * exp(this->epsilonbar * tau) * pow(this->epsilonbar, (int)2) * this->vbarn;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta2_dtau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->d2g_deta2(this->eta(delta)) * exp(this->epsilonbar * tau) * this->epsilonbar * this->kappabar * pow(this->vbarn, (int)2);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta3__consttau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    return this->d3g_deta3(this->eta(delta)) * (exp(this->epsilonbar * tau) - 1) * this->kappabar * pow(this->vbarn, (int)3);
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::X(const CoolPropDbl& delta, const CoolPropDbl& Deltabar) const {
    return 2 / (sqrt(1 + 4 * Deltabar * delta) + 1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_dDeltabar__constdelta(const CoolPropDbl& delta, const CoolPropDbl& Deltabar) const {
    CoolPropDbl X = this->X(delta, Deltabar);
    return -delta * X * X / (2 * Deltabar * delta * X + 1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_ddelta__constDeltabar(const CoolPropDbl& delta, const CoolPropDbl& Deltabar) const {
    CoolPropDbl X = this->X(delta, Deltabar);
    return -Deltabar * X * X / (2 * Deltabar * delta * X + 1);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_dtau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    return this->dX_dDeltabar__constdelta(delta, Deltabar) * this->dDeltabar_dtau__constdelta(tau, delta);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dX_ddelta(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
            + this->dX_dDeltabar__constdelta(delta, Deltabar) * this->dDeltabar_ddelta__consttau(tau, delta));
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_dtau2(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl beta = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl d_dXdtau_dbeta = -delta * X * X / (2 * Deltabar * delta * X + 1);
    CoolPropDbl d_dXdtau_dDeltabar = 2 * delta * delta * X * X * X / pow(2 * Deltabar * delta * X + 1, 2) * beta;
    CoolPropDbl d_dXdtau_dX = -2 * beta * delta * X * (Deltabar * delta * X + 1) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXdtau_dX * dX_dDeltabar * beta + d_dXdtau_dDeltabar * beta + d_dXdtau_dbeta * dbeta_dtau;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_ddeltadtau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl beta = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl d_dXddelta_dDeltabar = X * X * (2 * delta * delta * X * alpha - 1) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl d_dXddelta_dalpha = -delta * X * X / (2 * Deltabar * delta * X + 1);
    CoolPropDbl d_dXddelta_dX = -(Deltabar + delta * alpha) * 2 * (Deltabar * delta * X * X + X) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXddelta_dX * dX_dDeltabar * beta + d_dXddelta_dDeltabar * beta + d_dXddelta_dalpha * dalpha_dtau;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2X_ddelta2(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Deltabar = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Deltabar);
    CoolPropDbl alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);

    CoolPropDbl dX_ddelta_constall = X * X * (2 * Deltabar * Deltabar * X - alpha) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl d_dXddelta_dX = -(Deltabar + delta * alpha) * 2 * (Deltabar * delta * X * X + X) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl d_dXddelta_dDeltabar = X * X * (2 * delta * delta * X * alpha - 1) / pow(2 * Deltabar * delta * X + 1, 2);
    CoolPropDbl d_dXddelta_dalpha = -delta * X * X / (2 * Deltabar * delta * X + 1);

    CoolPropDbl dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Deltabar);

    CoolPropDbl val = (dX_ddelta_constall + d_dXddelta_dX * dX_ddelta + d_dXddelta_dX * dX_dDeltabar * alpha + d_dXddelta_dDeltabar * alpha
                       + d_dXddelta_dalpha * dalpha_ddelta);
    return val;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_dtau3(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl Delta_ttt = this->d3Deltabar_dtau3__constdelta(tau, delta);
    CoolPropDbl dXtt_dX = 2 * X * delta
                          * (-6 * Delta * pow(Delta_t, 2) * pow(X, 2) * pow(delta, 2) * (Delta * X * delta + 1)
                             + 3 * pow(Delta_t, 2) * X * delta * (2 * Delta * X * delta + 1) - Delta_tt * pow(2 * Delta * X * delta + 1, 3)
                             + X * delta * (Delta * Delta_tt + 3 * pow(Delta_t, 2)) * pow(2 * Delta * X * delta + 1, 2))
                          / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXtt_dDelta = 2 * pow(X, 3) * pow(delta, 2)
                              * (-6 * pow(Delta_t, 2) * X * delta * (Delta * X * delta + 1)
                                 - 3 * pow(Delta_t, 2) * X * delta * (2 * Delta * X * delta + 1) + Delta_tt * pow(2 * Delta * X * delta + 1, 2))
                              / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXtt_dDelta_t = 4 * Delta_t * pow(X, 3) * pow(delta, 2) * (3 * Delta * X * delta + 2) / pow(2 * Delta * X * delta + 1, 3);
    CoolPropDbl dXtt_dDelta_tt = -pow(X, 2) * delta / (2 * Delta * X * delta + 1);
    return dXtt_dX * dX_dDelta * Delta_t + dXtt_dDelta * Delta_t + dXtt_dDelta_t * Delta_tt + dXtt_dDelta_tt * Delta_ttt;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddeltadtau2(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    CoolPropDbl Delta_dtt = this->d3Deltabar_ddelta_dtau2(tau, delta);
    CoolPropDbl dXtt_ddelta =
      pow(X, 2)
      * (-12 * Delta * pow(Delta_t, 2) * pow(X, 2) * pow(delta, 2) * (Delta * X * delta + 1)
         + 2 * pow(Delta_t, 2) * X * delta * (-Delta * X * delta + 2) * (2 * Delta * X * delta + 1) - Delta_tt * pow(2 * Delta * X * delta + 1, 3)
         + 2 * X * delta * (Delta * Delta_tt + 2 * pow(Delta_t, 2)) * pow(2 * Delta * X * delta + 1, 2))
      / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXtt_dX = 2 * X * delta
                          * (-6 * Delta * pow(Delta_t, 2) * pow(X, 2) * pow(delta, 2) * (Delta * X * delta + 1)
                             + 3 * pow(Delta_t, 2) * X * delta * (2 * Delta * X * delta + 1) - Delta_tt * pow(2 * Delta * X * delta + 1, 3)
                             + X * delta * (Delta * Delta_tt + 3 * pow(Delta_t, 2)) * pow(2 * Delta * X * delta + 1, 2))
                          / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXtt_dDelta = 2 * pow(X, 3) * pow(delta, 2)
                              * (-6 * pow(Delta_t, 2) * X * delta * (Delta * X * delta + 1)
                                 - 3 * pow(Delta_t, 2) * X * delta * (2 * Delta * X * delta + 1) + Delta_tt * pow(2 * Delta * X * delta + 1, 2))
                              / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXtt_dDelta_t = 4 * Delta_t * pow(X, 3) * pow(delta, 2) * (3 * Delta * X * delta + 2) / pow(2 * Delta * X * delta + 1, 3);
    CoolPropDbl dXtt_dDelta_tt = -pow(X, 2) * delta / (2 * Delta * X * delta + 1);
    return dXtt_ddelta + dXtt_dX * dX_ddelta + dXtt_dX * dX_dDelta * Delta_d + dXtt_dDelta * Delta_d + dXtt_dDelta_t * Delta_dt
           + dXtt_dDelta_tt * Delta_dtt;
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddelta2dtau(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    CoolPropDbl Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    CoolPropDbl Delta_ddt = this->d3Deltabar_ddelta2_dtau(tau, delta);
    CoolPropDbl dXdd_dX =
      2 * X
      * (-6 * Delta * pow(X, 2) * delta * pow(Delta + Delta_d * delta, 2) * (Delta * X * delta + 1)
         - Delta_dd * delta * pow(2 * Delta * X * delta + 1, 3)
         + 2 * X * (2 * Delta * X * delta + 1)
             * (-Delta * Delta_d * delta * (2 * Delta_d * X * pow(delta, 2) - 1) - Delta * delta * (2 * pow(Delta, 2) * X - Delta_d)
                + Delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1) + Delta_d * delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1))
         + pow(2 * Delta * X * delta + 1, 2)
             * (3 * pow(Delta, 2) * X + Delta * Delta_dd * X * pow(delta, 2) + Delta * X * (Delta + Delta_d * delta)
                + pow(Delta_d, 2) * X * pow(delta, 2) + Delta_d * X * delta * (Delta + Delta_d * delta)
                + Delta_d * (2 * Delta_d * X * pow(delta, 2) - 1) - Delta_d))
      / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXdd_dDelta = pow(X, 3)
                              * (-8 * pow(Delta, 2) * Delta_d * pow(X, 2) * pow(delta, 3) + 8 * pow(Delta, 2) * Delta_dd * pow(X, 2) * pow(delta, 4)
                                 + 10 * pow(Delta, 2) * X * delta - 24 * Delta * pow(Delta_d, 2) * pow(X, 2) * pow(delta, 4)
                                 + 8 * Delta * Delta_d * X * pow(delta, 2) + 8 * Delta * Delta_dd * X * pow(delta, 3) + 8 * Delta
                                 - 18 * pow(Delta_d, 2) * X * pow(delta, 3) + 12 * Delta_d * delta + 2 * Delta_dd * pow(delta, 2))
                              / (16 * pow(Delta, 4) * pow(X, 4) * pow(delta, 4) + 32 * pow(Delta, 3) * pow(X, 3) * pow(delta, 3)
                                 + 24 * pow(Delta, 2) * pow(X, 2) * pow(delta, 2) + 8 * Delta * X * delta + 1);
    CoolPropDbl dXdd_dDelta_d =
      2 * pow(X, 2)
      * (2 * X * delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1) + (2 * Delta * X * delta + 1) * (2 * Delta_d * X * pow(delta, 2) - 1))
      / pow(2 * Delta * X * delta + 1, 3);
    CoolPropDbl dXdd_dDelta_dd = -pow(X, 2) * delta / (2 * Delta * X * delta + 1);

    return dXdd_dX * dX_dDelta * Delta_t + dXdd_dDelta * Delta_t + dXdd_dDelta_d * Delta_dt + dXdd_dDelta_dd * Delta_ddt;
}

double Xdd(double X, double delta, double Delta, double Delta_d, double Delta_dd) {
    return Delta * pow(X, 2) * (2 * Delta + 2 * Delta_d * delta) * (Delta * pow(X, 2) * delta + X) / pow(2 * Delta * X * delta + 1, 3)
           + Delta_d * pow(X, 2) * delta * (2 * Delta + 2 * Delta_d * delta) * (Delta * pow(X, 2) * delta + X) / pow(2 * Delta * X * delta + 1, 3)
           + Delta_d * pow(X, 2) * (2 * Delta_d * X * pow(delta, 2) - 1) / pow(2 * Delta * X * delta + 1, 2)
           - Delta_dd * pow(X, 2) * delta / (2 * Delta * X * delta + 1)
           + pow(X, 2) * (2 * pow(Delta, 2) * X - Delta_d) / pow(2 * Delta * X * delta + 1, 2);
}

CoolPropDbl ResidualHelmholtzSAFTAssociating::d3X_ddelta3(const CoolPropDbl& tau, const CoolPropDbl& delta) const {
    CoolPropDbl Delta = this->Deltabar(tau, delta);
    CoolPropDbl X = this->X(delta, Delta);
    CoolPropDbl dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    CoolPropDbl dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    CoolPropDbl Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    CoolPropDbl Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    CoolPropDbl Delta_ddd = this->d3Deltabar_ddelta3__consttau(tau, delta);

    CoolPropDbl dXdd_dX =
      2 * X
      * (-6 * Delta * pow(X, 2) * delta * pow(Delta + Delta_d * delta, 2) * (Delta * X * delta + 1)
         - Delta_dd * delta * pow(2 * Delta * X * delta + 1, 3)
         + 2 * X * (2 * Delta * X * delta + 1)
             * (-Delta * Delta_d * delta * (2 * Delta_d * X * pow(delta, 2) - 1) - Delta * delta * (2 * pow(Delta, 2) * X - Delta_d)
                + Delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1) + Delta_d * delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1))
         + pow(2 * Delta * X * delta + 1, 2)
             * (3 * pow(Delta, 2) * X + Delta * Delta_dd * X * pow(delta, 2) + Delta * X * (Delta + Delta_d * delta)
                + pow(Delta_d, 2) * X * pow(delta, 2) + Delta_d * X * delta * (Delta + Delta_d * delta)
                + Delta_d * (2 * Delta_d * X * pow(delta, 2) - 1) - Delta_d))
      / pow(2 * Delta * X * delta + 1, 4);
    CoolPropDbl dXdd_ddelta = pow(X, 2)
                              * (-24 * pow(Delta, 4) * pow(X, 3) * delta - 8 * pow(Delta, 3) * Delta_d * pow(X, 3) * pow(delta, 2)
                                 - 18 * pow(Delta, 3) * pow(X, 2) + 8 * pow(Delta, 2) * Delta_d * pow(X, 2) * delta
                                 - 4 * pow(Delta, 2) * Delta_dd * pow(X, 2) * pow(delta, 2) + 10 * Delta * pow(Delta_d, 2) * pow(X, 2) * pow(delta, 2)
                                 + 12 * Delta * Delta_d * X - 4 * Delta * Delta_dd * X * delta + 8 * pow(Delta_d, 2) * X * delta - Delta_dd)
                              / (16 * pow(Delta, 4) * pow(X, 4) * pow(delta, 4) + 32 * pow(Delta, 3) * pow(X, 3) * pow(delta, 3)
                                 + 24 * pow(Delta, 2) * pow(X, 2) * pow(delta, 2) + 8 * Delta * X * delta + 1);
    CoolPropDbl dXdd_dDelta = pow(X, 3)
                              * (-8 * pow(Delta, 2) * Delta_d * pow(X, 2) * pow(delta, 3) + 8 * pow(Delta, 2) * Delta_dd * pow(X, 2) * pow(delta, 4)
                                 + 10 * pow(Delta, 2) * X * delta - 24 * Delta * pow(Delta_d, 2) * pow(X, 2) * pow(delta, 4)
                                 + 8 * Delta * Delta_d * X * pow(delta, 2) + 8 * Delta * Delta_dd * X * pow(delta, 3) + 8 * Delta
                                 - 18 * pow(Delta_d, 2) * X * pow(delta, 3) + 12 * Delta_d * delta + 2 * Delta_dd * pow(delta, 2))
                              / (16 * pow(Delta, 4) * pow(X, 4) * pow(delta, 4) + 32 * pow(Delta, 3) * pow(X, 3) * pow(delta, 3)
                                 + 24 * pow(Delta, 2) * pow(X, 2) * pow(delta, 2) + 8 * Delta * X * delta + 1);
    CoolPropDbl dXdd_dDelta_d =
      2 * pow(X, 2)
      * (2 * X * delta * (Delta + Delta_d * delta) * (Delta * X * delta + 1) + (2 * Delta * X * delta + 1) * (2 * Delta_d * X * pow(delta, 2) - 1))
      / pow(2 * Delta * X * delta + 1, 3);
    CoolPropDbl dXdd_dDelta_dd = -pow(X, 2) * delta / (2 * Delta * X * delta + 1);

    return dXdd_ddelta + dXdd_dX * (dX_ddelta + dX_dDelta * Delta_d) + dXdd_dDelta * Delta_d + dXdd_dDelta_d * Delta_dd + dXdd_dDelta_dd * Delta_ddd;
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::g(const CoolPropDbl& eta) const {
    return 0.5 * (2 - eta) / pow(1 - eta, (int)3);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::dg_deta(const CoolPropDbl& eta) const {
    return 0.5 * (5 - 2 * eta) / pow(1 - eta, (int)4);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d2g_deta2(const CoolPropDbl& eta) const {
    return 3 * (3 - eta) / pow(1 - eta, (int)5);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::d3g_deta3(const CoolPropDbl& eta) const {
    return 6 * (7 - 2 * eta) / pow(1 - eta, (int)6);
}
CoolPropDbl ResidualHelmholtzSAFTAssociating::eta(const CoolPropDbl& delta) const {
    return this->vbarn * delta;
}

void ResidualHelmholtzSAFTAssociating::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& deriv) throw() {
    if (disabled) {
        return;
    }
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

    deriv.alphar += this->m * this->a * ((log(X) - X / 2.0 + 0.5));
    deriv.dalphar_ddelta += this->m * this->a * (1 / X - 0.5) * this->dX_ddelta(tau, delta);
    deriv.dalphar_dtau += this->m * this->a * (1 / X - 0.5) * this->dX_dtau(tau, delta);
    deriv.d2alphar_dtau2 += this->m * this->a * ((1 / X - 0.5) * X_tt - pow(X_t / X, 2));
    deriv.d2alphar_ddelta2 += this->m * this->a * ((1 / X - 0.5) * X_dd - pow(X_d / X, 2));
    deriv.d2alphar_ddelta_dtau += this->m * this->a * ((-X_t / X / X) * X_d + X_dt * (1 / X - 0.5));
    deriv.d3alphar_dtau3 += this->m * this->a
                            * ((1 / X - 1.0 / 2.0) * X_ttt + (-X_t / pow(X, (int)2)) * X_tt
                               - 2 * (pow(X, (int)2) * (X_t * X_tt) - pow(X_t, (int)2) * (X * X_t)) / pow(X, (int)4));
    deriv.d3alphar_ddelta_dtau2 += this->m * this->a
                                   * ((1 / X - 1.0 / 2.0) * X_dtt - X_d / pow(X, (int)2) * X_tt
                                      - 2 * (pow(X, (int)2) * (X_t * X_dt) - pow(X_t, (int)2) * (X * X_d)) / pow(X, (int)4));
    deriv.d3alphar_ddelta2_dtau += this->m * this->a
                                   * ((1 / X - 1.0 / 2.0) * X_ddt - X_t / pow(X, (int)2) * X_dd
                                      - 2 * (pow(X, (int)2) * (X_d * X_dt) - pow(X_d, (int)2) * (X * X_t)) / pow(X, (int)4));
    deriv.d3alphar_ddelta3 += this->m * this->a
                              * ((1 / X - 1.0 / 2.0) * X_ddd - X_d / pow(X, (int)2) * X_dd
                                 - 2 * (pow(X, (int)2) * (X_d * X_dd) - pow(X_d, (int)2) * (X * X_d)) / pow(X, (int)4));
}

void IdealHelmholtzCP0PolyT::to_json(rapidjson::Value& el, rapidjson::Document& doc) {
    el.AddMember("type", "IdealGasCP0Poly", doc.GetAllocator());

    rapidjson::Value _c(rapidjson::kArrayType), _t(rapidjson::kArrayType);
    for (std::size_t i = 0; i < N; ++i) {
        _c.PushBack(static_cast<double>(c[i]), doc.GetAllocator());
        _t.PushBack(static_cast<double>(t[i]), doc.GetAllocator());
    }
    el.AddMember("c", _c, doc.GetAllocator());
    el.AddMember("t", _t, doc.GetAllocator());
    el.AddMember("Tc", static_cast<double>(Tc), doc.GetAllocator());
    el.AddMember("T0", static_cast<double>(T0), doc.GetAllocator());
}

void IdealHelmholtzLead::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    derivs.alphar += log(delta) + a1 + a2 * tau;
    derivs.dalphar_ddelta += 1.0 / delta;
    derivs.dalphar_dtau += a2;
    derivs.d2alphar_ddelta2 += -1.0 / delta / delta;
    derivs.d3alphar_ddelta3 += 2 / delta / delta / delta;
    derivs.d4alphar_ddelta4 += -6 / POW4(delta);
}
void IdealHelmholtzEnthalpyEntropyOffset::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    derivs.alphar += a1 + a2 * tau;
    derivs.dalphar_dtau += a2;
}
void IdealHelmholtzLogTau::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {

    if (!enabled) {
        return;
    }
    derivs.alphar += a1 * log(tau);
    derivs.dalphar_dtau += a1 / tau;
    derivs.d2alphar_dtau2 += -a1 / tau / tau;
    derivs.d3alphar_dtau3 += 2 * a1 / tau / tau / tau;
    derivs.d4alphar_dtau4 += -6 * a1 / POW4(tau);
}
void IdealHelmholtzPower::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * pow(tau, t[i]);
        }
        derivs.alphar += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * t[i] * pow(tau, t[i] - 1);
        }
        derivs.dalphar_dtau += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * t[i] * (t[i] - 1) * pow(tau, t[i] - 2);
        }
        derivs.d2alphar_dtau2 += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * t[i] * (t[i] - 1) * (t[i] - 2) * pow(tau, t[i] - 3);
        }
        derivs.d3alphar_dtau3 += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * t[i] * (t[i] - 1) * (t[i] - 2) * (t[i] - 3) * pow(tau, t[i] - 4);
        }
        derivs.d4alphar_dtau4 += s;
    }
}
void IdealHelmholtzPlanckEinsteinGeneralized::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    // First pre-calculate exp(theta[i]*tau) for each contribution; used in each term
    std::vector<double> expthetatau(N);
    for (std::size_t i = 0; i < N; ++i) {
        expthetatau[i] = exp(theta[i] * tau);
    }

    if (!enabled) {
        return;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * log(c[i] + d[i] * expthetatau[i]);
        }
        derivs.alphar += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * theta[i] * d[i] * expthetatau[i] / (c[i] + d[i] * expthetatau[i]);
        }
        derivs.dalphar_dtau += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * POW2(theta[i]) * c[i] * d[i] * expthetatau[i] / pow(c[i] + d[i] * expthetatau[i], 2);
        }
        derivs.d2alphar_dtau2 += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            s += n[i] * POW3(theta[i]) * c[i] * d[i] * (c[i] - d[i] * expthetatau[i]) * expthetatau[i] / pow(c[i] + d[i] * expthetatau[i], 3);
        }
        derivs.d3alphar_dtau3 += s;
    }
    {
        CoolPropDbl s = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const CoolPropDbl para = c[i] + d[i] * expthetatau[i];
            const CoolPropDbl bracket = 6 * POW3(d[i]) * POW3(expthetatau[i]) - 12 * d[i] * d[i] * para * POW2(expthetatau[i])
                                        + 7 * d[i] * POW2(para) * expthetatau[i] - POW3(para);
            s += -n[i] * d[i] * pow(theta[i], 4) * bracket * expthetatau[i] / pow(c[i] + d[i] * expthetatau[i], 4);
        }
        derivs.d4alphar_dtau4 += s;
    }
}
void IdealHelmholtzCP0Constant::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    derivs.alphar += cp_over_R - cp_over_R * tau / tau0 + cp_over_R * log(tau / tau0);
    derivs.dalphar_dtau += cp_over_R / tau - cp_over_R / tau0;
    derivs.d2alphar_dtau2 += -cp_over_R / (tau * tau);
    derivs.d3alphar_dtau3 += 2 * cp_over_R / (tau * tau * tau);
    derivs.d4alphar_dtau4 += -6 * cp_over_R / POW4(tau);
}
void IdealHelmholtzCP0PolyT::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(t[i]) < 10 * DBL_EPSILON) {
                sum += c[i] - c[i] * tau / tau0 + c[i] * log(tau / tau0);
            } else if (std::abs(t[i] + 1) < 10 * DBL_EPSILON) {
                sum += c[i] * tau / Tc * log(tau0 / tau) + c[i] / Tc * (tau - tau0);
            } else {
                sum += -c[i] * pow(Tc, t[i]) * pow(tau, -t[i]) / (t[i] * (t[i] + 1)) - c[i] * pow(T0, t[i] + 1) * tau / (Tc * (t[i] + 1))
                       + c[i] * pow(T0, t[i]) / t[i];
            }
        }
        derivs.alphar += sum;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(t[i]) < 10 * DBL_EPSILON) {
                sum += c[i] / tau - c[i] / tau0;
            } else if (std::abs(t[i] + 1) < 10 * DBL_EPSILON) {
                sum += c[i] / Tc * log(tau0 / tau);
            } else {
                sum += c[i] * pow(Tc, t[i]) * pow(tau, -t[i] - 1) / (t[i] + 1) - c[i] * pow(Tc, t[i]) / (pow(tau0, t[i] + 1) * (t[i] + 1));
            }
        }
        derivs.dalphar_dtau += sum;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(t[i]) < 10 * DBL_EPSILON) {
                sum += -c[i] / (tau * tau);
            } else if (std::abs(t[i] + 1) < 10 * DBL_EPSILON) {
                sum += -c[i] / (tau * Tc);
            } else {
                sum += -c[i] * pow(Tc / tau, t[i]) / (tau * tau);
            }
        }
        derivs.d2alphar_dtau2 += sum;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(t[i]) < 10 * DBL_EPSILON) {
                sum += 2 * c[i] / (tau * tau * tau);
            } else if (std::abs(t[i] + 1) < 10 * DBL_EPSILON) {
                sum += c[i] / (tau * tau * Tc);
            } else {
                sum += c[i] * pow(Tc / tau, t[i]) * (t[i] + 2) / (tau * tau * tau);
            }
        }
        derivs.d3alphar_dtau3 += sum;
    }
    {
        double sum = 0;
        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(t[i]) < 10 * DBL_EPSILON) {
                sum += -6 * c[i] / POW4(tau);
            } else if (std::abs(t[i] + 1) < 10 * DBL_EPSILON) {
                sum += -3 * c[i] / (POW3(tau) * Tc);
            } else {
                sum += -c[i] * (t[i] + 2) * (t[i] + 3) * pow(Tc / tau, t[i]) / (tau * tau * tau * tau);
            }
        }
        derivs.d4alphar_dtau4 += sum;
    }
}

void IdealHelmholtzGERG2004Sinh::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    // Check that the reducing temperature value is provided
    CoolPropDbl T_red = _HUGE;
    if (ValidNumber(_Tr)) {
        T_red = _Tr;  // Primarily useful for testing
    } else if (ValidNumber(derivs.T_red)) {
        T_red = derivs.T_red;
    } else {
        throw ValueError("T_red needs to be stored somewhere for GERG2004Sinh");
    }
    CoolPropDbl Tci_over_Tr = Tc / T_red;

    double sum00 = 0, sum10 = 0, sum20 = 0, sum30 = 0, sum40 = 0;
    for (std::size_t i = 0; i < N; ++i) {
        CoolPropDbl t = theta[i] * Tci_over_Tr;
        sum00 += n[i] * log(std::abs(sinh(t * tau)));
        sum10 += n[i] * t / tanh(t * tau);
        sum20 += -n[i] * POW2(t) / POW2(sinh(t * tau));
        sum30 += -2 * n[i] * POW3(t) * (1 / tanh(t * tau) - 1 / POW3(tanh(t * tau)));
        sum40 += -2 * n[i] * POW4(t) * (1 - 4 / POW2(tanh(t * tau)) + 3 / POW4(tanh(t * tau)));
    }
    derivs.alphar += sum00;
    derivs.dalphar_dtau += sum10;
    derivs.d2alphar_dtau2 += sum20;
    derivs.d3alphar_dtau3 += sum30;
    derivs.d4alphar_dtau4 += sum40;
}

void IdealHelmholtzGERG2004Cosh::all(const CoolPropDbl& tau, const CoolPropDbl& delta, HelmholtzDerivatives& derivs) throw() {
    if (!enabled) {
        return;
    }
    // Check that the reducing temperature value is provided in the derivs structure
    CoolPropDbl T_red = _HUGE;
    if (ValidNumber(_Tr)) {
        T_red = _Tr;  // Primarily useful for testing
    } else if (ValidNumber(derivs.T_red)) {
        T_red = derivs.T_red;
    } else {
        throw ValueError("T_red needs to be stored somewhere for GERG2004Cosh");
    }
    CoolPropDbl Tci_over_Tr = Tc / T_red;

    double sum00 = 0, sum10 = 0, sum20 = 0, sum30 = 0, sum40 = 0;
    for (std::size_t i = 0; i < N; ++i) {
        CoolPropDbl t = theta[i] * Tci_over_Tr;
        sum00 += -n[i] * log(std::abs(cosh(t * tau)));
        sum10 += -n[i] * t * tanh(t * tau);
        sum20 += -n[i] * POW2(t) / POW2(cosh(t * tau));
        sum30 += -2 * n[i] * POW3(t) * (POW3(tanh(t * tau)) - tanh(t * tau));
        sum40 += 2 * n[i] * POW4(t) * (3 * POW4(tanh(t * tau)) - 4 * POW2(tanh(t * tau)) + 1);
    }
    derivs.alphar += sum00;
    derivs.dalphar_dtau += sum10;
    derivs.d2alphar_dtau2 += sum20;
    derivs.d3alphar_dtau3 += sum30;
    derivs.d4alphar_dtau4 += sum40;
}

//void IdealHelmholtzCP0AlyLee::to_json(rapidjson::Value &el, rapidjson::Document &doc){
//    el.AddMember("type","IdealGasHelmholtzCP0AlyLee",doc.GetAllocator());
//    rapidjson::Value _n(rapidjson::kArrayType);
//    for (std::size_t i=0; i<=4; ++i)
//    {
//        _n.PushBack(static_cast<double>(c[i]),doc.GetAllocator());
//    }
//    el.AddMember("c",_n,doc.GetAllocator());
//    el.AddMember("Tc",static_cast<double>(Tc),doc.GetAllocator());
//    el.AddMember("T0",static_cast<double>(T0),doc.GetAllocator());
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::base(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
//{
//    if (!enabled){ return 0.0;}
//    return -tau*(anti_deriv_cp0_tau2(tau)-anti_deriv_cp0_tau2(tau0))+(anti_deriv_cp0_tau(tau)-anti_deriv_cp0_tau(tau0));
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::dTau(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
//{
//    if (!enabled){ return 0.0;}
//    return -(anti_deriv_cp0_tau2(tau) - anti_deriv_cp0_tau2(tau0));
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau2(const CoolPropDbl &tau)
//{
//    return -c[0]/tau + 2*c[1]*c[2]/Tc/(exp(-2*c[2]*tau/Tc)-1) - 2*c[3]*c[4]/Tc*(exp(2*c[4]*tau/Tc)+1);
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau(const CoolPropDbl &tau)
//{
//    CoolPropDbl term1 = c[0]*log(tau);
//    CoolPropDbl term2 = 2*c[1]*c[2]*tau/(-Tc + Tc*exp(-2*c[2]*tau/Tc)) + c[1]*log(-1 + exp(-2*c[2]*tau/Tc)) + 2*c[1]*c[2]*tau/Tc;
//    CoolPropDbl term3 = -c[3]*(Tc*exp(2*c[4]*tau/Tc)*log(exp(2*c[4]*tau/Tc) + 1) + Tc*log(exp(2*c[4]*tau/Tc) + 1) - 2*c[4]*tau*exp(2*c[4]*tau/Tc))/(Tc*(exp(2*c[4]*tau/Tc) + 1));
//    return term1 + term2 + term3;
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::dTau2(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
//{
//    if (!enabled){ return 0.0;}
//    return -c[0]/pow(tau,2) - c[1]*pow(c[2]/Tc/sinh(c[2]*tau/Tc),2) - c[3]*pow(c[4]/Tc/cosh(c[4]*tau/Tc),2);
//}
//CoolPropDbl IdealHelmholtzCP0AlyLee::dTau3(const CoolPropDbl &tau, const CoolPropDbl &delta) throw()
//{
//    if (!enabled){ return 0.0;}
//    return 2*c[0]/pow(tau,3) + 2*c[1]*pow(c[2]/Tc,3)*cosh(c[2]*tau/Tc)/pow(sinh(c[2]*tau/Tc),3) + 2*c[3]*pow(c[4]/Tc,3)*sinh(c[4]*tau/Tc)/pow(cosh(c[4]*tau/Tc),3);
//}

}; /* namespace CoolProp */

/*
IdealHelmholtzEnthalpyEntropyOffset EnthalpyEntropyOffset;
*/

#ifdef ENABLE_CATCH
#    include <math.h>
#    include <catch2/catch_all.hpp>
#    include "crossplatform_shared_ptr.h"

class HelmholtzConsistencyFixture
{
   public:
    CoolPropDbl numerical, analytic;

    shared_ptr<CoolProp::BaseHelmholtzTerm> PlanckEinstein, Lead, LogTau, IGPower, CP0Constant, CP0PolyT, SAFT, NonAnalytic, Soave, PR, XiangDeiters,
      GaoB, GERG2004Cosh, GERG2004Sinh;
    shared_ptr<CoolProp::ResidualHelmholtzGeneralizedExponential> Gaussian, Lemmon2005, Exponential, GERG2008, Power;

    HelmholtzConsistencyFixture() {
        shared_ptr<AbstractCubic> _SRK(new SRK(300, 4e6, 0.3, 8.314461));
        _SRK->set_Tr(300);
        _SRK->set_rhor(4000);
        Soave.reset(new CoolProp::ResidualHelmholtzGeneralizedCubic(_SRK));

        shared_ptr<AbstractCubic> _PR(new PengRobinson(300, 4e6, 0.3, 8.314461));
        _PR->set_Tr(300);
        _PR->set_rhor(4000);
        PR.reset(new CoolProp::ResidualHelmholtzGeneralizedCubic(_PR));

        {
            CoolPropDbl beta[] = {0.3696, 0.2962}, epsilon[] = {0.4478, 0.44689}, eta[] = {-2.8452, -2.8342}, gamma[] = {1.108, 1.313},
                        n[] = {-1.6909858, 0.93739074}, t[] = {4.3315, 4.015}, d[] = {1, 1}, b[] = {1.244, 0.6826};
            GaoB.reset(new CoolProp::ResidualHelmholtzGaoB(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(t[0])),
              std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(d[0])), std::vector<CoolPropDbl>(eta, eta + sizeof(eta) / sizeof(eta[0])),
              std::vector<CoolPropDbl>(beta, beta + sizeof(beta) / sizeof(beta[0])),
              std::vector<CoolPropDbl>(gamma, gamma + sizeof(gamma) / sizeof(gamma[0])),
              std::vector<CoolPropDbl>(epsilon, epsilon + sizeof(epsilon) / sizeof(epsilon[0])),
              std::vector<CoolPropDbl>(b, b + sizeof(b) / sizeof(b[0]))));
        }

        XiangDeiters.reset(new CoolProp::ResidualHelmholtzXiangDeiters(300, 4e6, 4000, 0.3, 8.3144621));

        Lead.reset(new CoolProp::IdealHelmholtzLead(1, 3));
        LogTau.reset(new CoolProp::IdealHelmholtzLogTau(1.5));
        {
            std::vector<CoolPropDbl> n(4, 0), t(4, 1);
            n[0] = -0.1;
            n[2] = 0.1;
            t[1] = -1;
            t[2] = -2;
            t[3] = 2;
            IGPower.reset(new CoolProp::IdealHelmholtzPower(n, t));
        }
        {
            std::vector<CoolPropDbl> n(4, 0), t(4, 1), c(4, 1), d(4, -1);
            n[0] = 0.1;
            n[2] = 0.5;
            t[0] = -1.5;
            t[1] = -1;
            t[2] = -2;
            t[3] = -2;
            PlanckEinstein.reset(new CoolProp::IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d));
        }
        {
            std::vector<CoolPropDbl> c(3, 1), t(3, 0);
            t[1] = 1;
            t[2] = 2;
            c[1] = 2;
            c[2] = 3;
            CoolPropDbl T0 = 273.15, Tc = 345.857;
            CP0PolyT.reset(new CoolProp::IdealHelmholtzCP0PolyT(c, t, Tc, T0));
        }
        CP0Constant.reset(new CoolProp::IdealHelmholtzCP0Constant(4 / 8.314472, 300, 250));
        {
            // Nitrogen
            std::vector<CoolPropDbl> n(2, 0.0);
            n[0] = 0.137320000;
            n[1] = 0.900660000;
            std::vector<CoolPropDbl> theta(2, 0.0);
            theta[0] = 5.251822620;
            theta[1] = 13.788988208;
            CoolPropDbl rhomolar_crit = 11183.900000, T_crit = 126.192000000;
            GERG2004Cosh.reset(new CoolProp::IdealHelmholtzGERG2004Cosh(n, theta, T_crit));
            static_cast<CoolProp::IdealHelmholtzGERG2004Cosh*>(GERG2004Cosh.get())->set_Tred(T_crit * 1.3);
        }
        {
            // Nitrogen
            std::vector<CoolPropDbl> n(1, 0.0);
            n[0] = -0.146600000;
            std::vector<CoolPropDbl> theta(1, 0.0);
            theta[0] = -5.393067706;
            CoolPropDbl rhomolar_crit = 11183.900000, T_crit = 126.192000000;
            GERG2004Sinh.reset(new CoolProp::IdealHelmholtzGERG2004Sinh(n, theta, T_crit));
            static_cast<CoolProp::IdealHelmholtzGERG2004Sinh*>(GERG2004Sinh.get())->set_Tred(T_crit * 1.3);
        }

        {
            CoolPropDbl beta[] = {1.24, 0.821, 15.45, 2.21, 437, 0.743}, d[] = {1, 1, 2, 2, 3, 3},
                        epsilon[] = {0.6734, 0.9239, 0.8636, 1.0507, 0.8482, 0.7522}, eta[] = {0.9667, 1.5154, 1.0591, 1.6642, 12.4856, 0.9662},
                        gamma[] = {1.2827, 0.4317, 1.1217, 1.1871, 1.1243, 0.4203},
                        n[] = {1.2198, -0.4883, -0.0033293, -0.0035387, -0.51172, -0.16882}, t[] = {1, 2.124, 0.4, 3.5, 0.5, 2.7};
            Gaussian.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Gaussian->add_Gaussian(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(d[0])),
              std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(t[0])), std::vector<CoolPropDbl>(eta, eta + sizeof(eta) / sizeof(eta[0])),
              std::vector<CoolPropDbl>(epsilon, epsilon + sizeof(epsilon) / sizeof(epsilon[0])),
              std::vector<CoolPropDbl>(beta, beta + sizeof(beta) / sizeof(beta[0])),
              std::vector<CoolPropDbl>(gamma, gamma + sizeof(gamma) / sizeof(gamma[0])));
        }
        {
            CoolPropDbl d[] = {1, 1, 1, 2, 4, 1, 1, 2, 2, 3, 4, 5, 1, 5, 1, 2, 3, 5}, l[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 2, 3, 3},
                        m[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.7, 7, 6},
                        n[] = {5.28076,   -8.67658,   0.7501127,   0.7590023,   0.01451899, 4.777189,    -3.330988, 3.775673,    -2.290919,
                               0.8888268, -0.6234864, -0.04127263, -0.08455389, -0.1308752, 0.008344962, -1.532005, -0.05883649, 0.02296658},
                        t[] = {0.669, 1.05, 2.75, 0.956, 1, 2, 2.75, 2.38, 3.37, 3.47, 2.63, 3.45, 0.72, 4.23, 0.2, 4.5, 29, 24};
            Lemmon2005.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Lemmon2005->add_Lemmon2005(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(d[0])),
              std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(t[0])), std::vector<CoolPropDbl>(l, l + sizeof(l) / sizeof(l[0])),
              std::vector<CoolPropDbl>(m, m + sizeof(m) / sizeof(m[0])));
        }
        {
            CoolPropDbl d[] = {1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2}, l[] = {0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3},
                        n[] = {1.0038, -2.7662, 0.42921, 0.081363, 0.00024174, 0.48246, 0.75542, -0.00743, -0.4146, -0.016558, -0.10644, -0.021704},
                        t[] = {0.25, 1.25, 1.5, 0.25, 0.875, 2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5};
            Power.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Power->add_Power(std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(d[0])),
                             std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(t[0])), std::vector<CoolPropDbl>(l, l + sizeof(l) / sizeof(l[0])));
        }
        {

            CoolPropDbl a = 1, epsilonbar = 12.2735737, kappabar = 1.09117041e-05, m = 1.01871348, vbarn = 0.0444215309;
            SAFT.reset(new CoolProp::ResidualHelmholtzSAFTAssociating(a, m, epsilonbar, vbarn, kappabar));
        }
        {
            CoolPropDbl n[] = {-0.666422765408, 0.726086323499, 0.0550686686128}, A[] = {0.7, 0.7, 0.7}, B[] = {0.3, 0.3, 1}, C[] = {10, 10, 12.5},
                        D[] = {275, 275, 275}, a[] = {3.5, 3.5, 3}, b[] = {0.875, 0.925, 0.875}, beta[] = {0.3, 0.3, 0.3};
            NonAnalytic.reset(new CoolProp::ResidualHelmholtzNonAnalytic(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(a, a + sizeof(a) / sizeof(a[0])),
              std::vector<CoolPropDbl>(b, b + sizeof(b) / sizeof(b[0])), std::vector<CoolPropDbl>(beta, beta + sizeof(beta) / sizeof(beta[0])),
              std::vector<CoolPropDbl>(A, A + sizeof(A) / sizeof(A[0])), std::vector<CoolPropDbl>(B, B + sizeof(B) / sizeof(B[0])),
              std::vector<CoolPropDbl>(C, C + sizeof(C) / sizeof(C[0])), std::vector<CoolPropDbl>(D, D + sizeof(D) / sizeof(D[0]))));
        }
        {
            CoolPropDbl d[] = {2, 2, 2, 0, 0, 0}, g[] = {1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788, 1.65533788},
                        l[] = {2, 2, 2, 2, 2, 2},
                        n[] = {-3.821884669859, 8.30345065618981, -4.4832307260286, -1.02590136933231, 2.20786016506394, -1.07889905203761},
                        t[] = {3, 4, 5, 3, 4, 5};
            Exponential.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            Exponential->add_Exponential(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(n[0])),
              std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(d[0])), std::vector<CoolPropDbl>(g, g + sizeof(g) / sizeof(t[0])),
              std::vector<CoolPropDbl>(l, l + sizeof(l) / sizeof(l[0])));
        }
        {
            CoolPropDbl d[] = {1, 4, 1, 2, 2, 2, 2, 2, 3}, t[] = {0.0, 1.85, 7.85, 5.4, 0.0, 0.75, 2.8, 4.45, 4.25},
                        n[] = {-0.0098038985517335, 0.00042487270143005, -0.034800214576142, -0.13333813013896, -0.011993694974627,
                               0.069243379775168,   -0.31022508148249,   0.24495491753226,   0.22369816716981},
                        eta[] = {0.0, 0.0, 1.0, 1.0, 0.25, 0.0, 0.0, 0.0, 0.0}, epsilon[] = {0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
                        beta[] = {0.0, 0.0, 1.0, 1.0, 2.5, 3.0, 3.0, 3.0, 3.0}, gamma[] = {0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
            GERG2008.reset(new CoolProp::ResidualHelmholtzGeneralizedExponential());
            GERG2008->add_GERG2008Gaussian(
              std::vector<CoolPropDbl>(n, n + sizeof(n) / sizeof(n[0])), std::vector<CoolPropDbl>(d, d + sizeof(d) / sizeof(n[0])),
              std::vector<CoolPropDbl>(t, t + sizeof(t) / sizeof(d[0])), std::vector<CoolPropDbl>(eta, eta + sizeof(eta) / sizeof(eta[0])),
              std::vector<CoolPropDbl>(epsilon, epsilon + sizeof(epsilon) / sizeof(epsilon[0])),
              std::vector<CoolPropDbl>(beta, beta + sizeof(beta) / sizeof(beta[0])),
              std::vector<CoolPropDbl>(gamma, gamma + sizeof(gamma) / sizeof(gamma[0])));
        }
    }
    void call(std::string d, shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        if (!d.compare("dTau")) {
            return dTau(term, tau, delta, ddelta);
        } else if (!d.compare("dTau2")) {
            return dTau2(term, tau, delta, ddelta);
        } else if (!d.compare("dTau3")) {
            return dTau3(term, tau, delta, ddelta);
        } else if (!d.compare("dTau4")) {
            return dTau4(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta")) {
            return dDelta(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta2")) {
            return dDelta2(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta3")) {
            return dDelta3(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta4")) {
            return dDelta4(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta_dTau")) {
            return dDelta_dTau(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta_dTau2")) {
            return dDelta_dTau2(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta2_dTau")) {
            return dDelta2_dTau(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta3_dTau")) {
            return dDelta3_dTau(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta2_dTau2")) {
            return dDelta2_dTau2(term, tau, delta, ddelta);
        } else if (!d.compare("dDelta_dTau3")) {
            return dDelta_dTau3(term, tau, delta, ddelta);
        } else {
            throw CoolProp::ValueError("don't understand deriv type");
        }
    }
    shared_ptr<CoolProp::BaseHelmholtzTerm> get(std::string t) {
        if (!t.compare("Lead")) {
            return Lead;
        } else if (!t.compare("LogTau")) {
            return LogTau;
        } else if (!t.compare("IGPower")) {
            return IGPower;
        } else if (!t.compare("PlanckEinstein")) {
            return PlanckEinstein;
        } else if (!t.compare("CP0Constant")) {
            return CP0Constant;
        } else if (!t.compare("CP0PolyT")) {
            return CP0PolyT;
        } else if (!t.compare("GERG2004Cosh")) {
            return GERG2004Cosh;
        } else if (!t.compare("GERG2004Sinh")) {
            return GERG2004Sinh;
        } else if (!t.compare("SRK")) {
            return Soave;
        } else if (!t.compare("PengRobinson")) {
            return PR;
        } else if (!t.compare("XiangDeiters")) {
            return XiangDeiters;
        } else if (!t.compare("GaoB")) {
            return GaoB;
        }

        else if (!t.compare("Gaussian")) {
            return Gaussian;
        } else if (!t.compare("Lemmon2005")) {
            return Lemmon2005;
        } else if (!t.compare("Power")) {
            return Power;
        } else if (!t.compare("SAFT")) {
            return SAFT;
        } else if (!t.compare("NonAnalytic")) {
            return NonAnalytic;
        } else if (!t.compare("Exponential")) {
            return Exponential;
        } else if (!t.compare("GERG2008")) {
            return GERG2008;
        } else {
            throw CoolProp::ValueError(format("don't understand helmholtz type: %s", t.c_str()));
        }
    }
    void dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau) {
        CoolPropDbl term_plus = term->base(tau + dtau, delta);
        CoolPropDbl term_minus = term->base(tau - dtau, delta);
        numerical = (term_plus - term_minus) / (2 * dtau);
        analytic = term->dTau(tau, delta);
    };
    void dTau2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau) {
        CoolPropDbl term_plus = term->dTau(tau + dtau, delta);
        CoolPropDbl term_minus = term->dTau(tau - dtau, delta);
        numerical = (term_plus - term_minus) / (2 * dtau);
        analytic = term->dTau2(tau, delta);
    };
    void dTau3(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau) {
        CoolPropDbl term_plus = term->dTau2(tau + dtau, delta);
        CoolPropDbl term_minus = term->dTau2(tau - dtau, delta);
        numerical = (term_plus - term_minus) / (2 * dtau);
        analytic = term->dTau3(tau, delta);
    };
    void dTau4(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl dtau) {
        CoolPropDbl term_plus = term->dTau3(tau + dtau, delta);
        CoolPropDbl term_minus = term->dTau3(tau - dtau, delta);
        numerical = (term_plus - term_minus) / (2 * dtau);
        analytic = term->dTau4(tau, delta);
    };
    void dDelta(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->base(tau, delta + ddelta);
        CoolPropDbl term_minus = term->base(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta(tau, delta);
    };
    void dDelta2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta2(tau, delta);
    };
    void dDelta3(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta2(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta2(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta3(tau, delta);
    };
    void dDelta4(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta3(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta3(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta4(tau, delta);
    };
    void dDelta_dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dTau(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dTau(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta_dTau(tau, delta);
    };
    void dDelta_dTau2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dTau2(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dTau2(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta_dTau2(tau, delta);
    };
    void dDelta2_dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta_dTau(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta_dTau(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta2_dTau(tau, delta);
    };
    void dDelta3_dTau(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta2_dTau(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta2_dTau(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta3_dTau(tau, delta);
    };
    void dDelta2_dTau2(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dDelta_dTau2(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dDelta_dTau2(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta2_dTau2(tau, delta);
    };
    void dDelta_dTau3(shared_ptr<CoolProp::BaseHelmholtzTerm> term, CoolPropDbl tau, CoolPropDbl delta, CoolPropDbl ddelta) {
        CoolPropDbl term_plus = term->dTau3(tau, delta + ddelta);
        CoolPropDbl term_minus = term->dTau3(tau, delta - ddelta);
        numerical = (term_plus - term_minus) / (2 * ddelta);
        analytic = term->dDelta_dTau3(tau, delta);
    };
    double err(double v1, double v2) {
        if (std::abs(v2) > 1e-15) {
            return std::abs((v1 - v2) / v2);
        } else {
            return std::abs(v1 - v2);
        }
    }
};

std::string terms[] = {"Lead",         "LogTau",       "IGPower", "PlanckEinstein", "CP0Constant", "CP0PolyT", "Gaussian",
                       "Lemmon2005",   "Power",        "SAFT",    "NonAnalytic",    "Exponential", "GERG2008", "SRK",
                       "PengRobinson", "XiangDeiters", "GaoB",    "GERG2004Cosh",   "GERG2004Sinh"};
std::string derivs[] = {"dTau",         "dTau2",        "dTau3", "dDelta",       "dDelta2",       "dDelta3",      "dDelta_dTau",
                        "dDelta_dTau2", "dDelta2_dTau", "dTau4", "dDelta_dTau3", "dDelta2_dTau2", "dDelta3_dTau", "dDelta4"};

TEST_CASE_METHOD(HelmholtzConsistencyFixture, "Helmholtz energy derivatives", "[helmholtz]") {
    shared_ptr<CoolProp::BaseHelmholtzTerm> term;
    std::size_t n = sizeof(terms) / sizeof(terms[0]);
    for (std::size_t i = 0; i < n; ++i) {
        term = get(terms[i]);
        for (std::size_t j = 0; j < sizeof(derivs) / sizeof(derivs[0]); ++j) {
            if (terms[i] == "SAFT"
                && (derivs[j] == "dTau4" || derivs[j] == "dDelta_dTau3" || derivs[j] == "dDelta2_dTau2" || derivs[j] == "dDelta3_dTau"
                    || derivs[j] == "dDelta4")) {
                continue;
            }
            call(derivs[j], term, 1.3, 0.9, 1e-6);
            CAPTURE(derivs[j]);
            CAPTURE(numerical);
            CAPTURE(analytic);
            CAPTURE(terms[i]);
            CHECK(err(analytic, numerical) < 1e-8);
        }
    }
}

#endif
