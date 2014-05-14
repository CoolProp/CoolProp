#include <numeric>
#include "Helmholtz.h"

namespace CoolProp{

void ResidualHelmholtzPower::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzPower",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType);
    for (unsigned int i = 0; i <= N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        _n.PushBack((double)el.n,doc.GetAllocator());
        _d.PushBack((double)el.d,doc.GetAllocator());
        _t.PushBack((double)el.t,doc.GetAllocator());
        _l.PushBack((double)el.l,doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
}
long double ResidualHelmholtzPower::base(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li > 0){
            s[i] = ni*exp(ti*log_tau+di*log_delta-pow(delta,li));
        }
        else{
            s[i] = ni*exp(ti*log_tau+di*log_delta);
        }
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld, pow_delta_li;
        int li = el.l;
        if (li > 0){
            pow_delta_li = pow(delta, li);
            s[i] = ni*(di-lid*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else{
            s[i] = ni*di*exp(ti*log_tau+(di-1)*log_delta);
        }
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta-pow(delta, li));
        }
        else
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

long double ResidualHelmholtzPower::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*((di-li*pow_delta_li)*(di-1.0-li*pow_delta_li) - li*li*pow_delta_li)*exp(ti*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*exp(ti*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(di-li*pow_delta_li)*exp((ti-1)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*exp((ti-1)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            long double bracket = (di*(di-1)*(di-2)+pow_delta_li*(-2*li+6*di*li-3*di*di*li-3*di*li*li+3*li*li-li*li*li)+pow_delta_li*pow_delta_li*(3*di*li*li-3*li*li+3*li*li*li)-li*li*li*pow_delta_li*pow_delta_li*pow_delta_li);
            s[i] = ni*bracket*exp(ti*log_tau+(di-3)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*(di-2)*exp(ti*log_tau+(di-3)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(((di-li*pow_delta_li))*(di-1-li*pow_delta_li)-li*li*pow_delta_li)*exp((ti-1)*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*(di-1)*exp((ti-1)*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(ti-1)*(di-li*pow_delta_li)*exp((ti-2)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*ti*(ti-1)*di*exp((ti-2)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t;
        int li = el.l;
        if (li > 0)
        {
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta-pow(delta,li));
        }
        else
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

void ResidualHelmholtzExponential::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzExponential",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), _g(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        _n.PushBack((double)el.n, doc.GetAllocator());
        _d.PushBack((double)el.d, doc.GetAllocator());
        _t.PushBack((double)el.t, doc.GetAllocator());
        _l.PushBack((double)el.l, doc.GetAllocator());
        _g.PushBack((double)el.g, doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
    el.AddMember("g",_g,doc.GetAllocator());
}
long double ResidualHelmholtzExponential::base(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        s[i] = ni*exp(ti*log_tau+di*log_delta-gi*pow(delta, li));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dDelta(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*(di-gi*li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta-gi*pow(delta,li));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

long double ResidualHelmholtzExponential::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        // Typo in Span, 2000, re-derived from Sympy
        double bracket = di*di - 2*di*pow(delta,li)*gi*li - di + pow(delta,2*li)*gi*gi*li*li - pow(delta,li)*gi*li*li + pow(delta,li)*gi*li;
        s[i] = ni*bracket*exp(ti*log_tau+(di-2)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*ti*(di-gi*li*pow_delta_li)*exp((ti-1)*log_tau+(di-1)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta-gi*pow(delta,li));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        // >> n_i, tau, t_i, d_i, delta, g_i, l_i = symbols(' n_i tau t_i d_i delta g_i l_i')
        // >> phir = n_i*tau**t_i*delta**d_i*exp(-g_i*pow(delta,l_i))
        // >> simplify(diff(diff(diff(phir,delta),delta),delta))
        long double pow_delta_2li = pow(delta,2*li);
        long double pow_delta_3li = pow(delta,3*li);
        long double bracket = di*di*di - 3*di*di*pow_delta_li*gi*li - 3*di*di + 3*di*pow_delta_2li*gi*gi*li*li - 3*di*pow_delta_li*gi*li*li + 6*di*pow_delta_li*gi*li + 2*di - pow_delta_3li*gi*gi*gi*li*li*li + 3*pow_delta_2li*gi*gi*li*li*li - 3*pow_delta_2li*gi*gi*li*li - pow_delta_li*gi*li*li*li + 3*pow_delta_li*gi*li*li - 2*pow_delta_li*gi*li;
        s[i] = ni*bracket*exp(ti*log_tau+(di-3)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        // Typo in Span, 2000, re-derived from Sympy
        long double bracket = di*di - 2*di*pow(delta,li)*gi*li - di + pow(delta,2*li)*gi*gi*li*li - pow(delta, li)*gi*li*li + pow(delta,li)*gi*li;
        s[i] = ni*ti*bracket*exp((ti-1)*log_tau+(di-2)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*ti*(ti-1)*(di-gi*li*pow_delta_li)*exp((ti-2)*log_tau+(di-1)*log_delta-gi*pow_delta_li);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzExponential::dTau3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzExponentialElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, gi = el.g, li = el.l;
        pow_delta_li = pow(delta, li);
        s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta-gi*pow(delta,li));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};



void ResidualHelmholtzGaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzGaussian",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
        _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        _n.PushBack((double)el.n,doc.GetAllocator());
        _d.PushBack((double)el.d,doc.GetAllocator());
        _t.PushBack((double)el.t,doc.GetAllocator());
        _eta.PushBack((double)el.eta,doc.GetAllocator());
        _epsilon.PushBack((double)el.epsilon,doc.GetAllocator());
        _beta.PushBack((double)el.beta,doc.GetAllocator());
        _gamma.PushBack((double)el.gamma,doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("eta",_eta,doc.GetAllocator());
    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("gamma",_gamma,doc.GetAllocator());
}
long double ResidualHelmholtzGaussian::base(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.t/tau-2.0*el.beta*(tau-el.gamma));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(tau,el.t)*psi*(-2.0*el.eta*pow(delta,el.d)+4.0*pow(el.eta,2)*pow(delta,el.d)*pow(delta-el.epsilon,2)-4.0*el.d*el.eta*pow(delta,el.d-1)*(delta-el.epsilon)+el.d*(el.d-1.0)*pow(delta,el.d-2));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon))*(el.t/tau-2.0*el.beta*(tau-el.gamma));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        long double bracket = pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta;
        s[i] = el.n*pow(tau,el.t)*pow(delta,el.d-3)*psi*bracket;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(tau,el.t)*psi*(el.t/tau-2.0*el.beta*(tau-el.gamma))*(-2.0*el.eta*pow(delta,el.d)+4.0*pow(el.eta,2)*pow(delta,el.d)*pow(delta-el.epsilon,2)-4.0*el.d*el.eta*pow(delta,el.d-1)*(delta-el.epsilon)+el.d*(el.d-1.0)*pow(delta,el.d-2));
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon))*(pow(el.t-2.0*el.beta*tau*(tau-el.gamma),2)-el.t-2*el.beta*tau*tau)/tau/tau;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        long double dpsi_dTau = exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2))*(-2*el.beta*(tau-el.gamma));

        long double bracket = pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta;
        long double dbracket_dTau = 2*(el.t/tau-2.0*el.beta*(tau-el.gamma))*(-el.t/tau/tau-2*el.beta)+2*el.t/pow(tau,3);
        s[i] = el.n*pow(delta,el.d)*(el.t*pow(tau,el.t-1)*psi*bracket+pow(tau,el.t)*dpsi_dTau*bracket+pow(tau,el.t)*psi*dbracket_dTau);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

void ResidualHelmholtzGERG2008Gaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzGERG2008Gaussian",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
                     _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), 
                     _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        _n.PushBack((double)el.n, doc.GetAllocator());
        _d.PushBack((double)el.d, doc.GetAllocator());
        _t.PushBack((double)el.t, doc.GetAllocator());
        _eta.PushBack((double)el.eta, doc.GetAllocator());
        _epsilon.PushBack((double)el.epsilon, doc.GetAllocator());
        _beta.PushBack((double)el.beta, doc.GetAllocator());
        _gamma.PushBack((double)el.gamma, doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("eta",_eta,doc.GetAllocator());
    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("gamma",_gamma,doc.GetAllocator());
}
long double ResidualHelmholtzGERG2008Gaussian::base(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*pow(tau,el.t)*pow(delta,el.d)*psi;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dDelta(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*pow(tau,el.t)*pow(delta,el.d)*psi*(el.d/delta-2*el.eta*(delta-el.epsilon)-el.beta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dTau(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*pow(tau,el.t-1)*pow(delta,el.d)*psi;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dDelta2(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*pow(tau,el.t)*pow(delta,el.d)*psi*(pow(el.d/delta-2*el.eta*(delta-el.epsilon)-el.beta,2)-el.d/delta/delta-2*el.eta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*pow(tau,el.t-1)*pow(delta,el.d)*psi*(el.d/delta-2*el.eta*(delta-el.epsilon)-el.beta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dTau2(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*(el.t-1)*pow(tau,el.t-2)*pow(delta,el.d)*psi;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};


long double ResidualHelmholtzGERG2008Gaussian::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*pow(tau,el.t-1)*pow(delta,el.d)*psi*(pow(el.d/delta-2*el.eta*(delta-el.epsilon)-el.beta,2)-el.d/delta/delta-2*el.eta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*(el.t-1)*pow(tau,el.t-2)*pow(delta,el.d)*psi*(el.d/delta-2*el.eta*(delta-el.epsilon)-el.beta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGERG2008Gaussian::dTau3(const long double &tau, const long double &delta) throw()
{
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        long double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*(delta-el.gamma));
        s[i] = el.n*el.t*(el.t-1)*(el.t-2)*pow(tau,el.t-3)*pow(delta,el.d)*psi;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

void ResidualHelmholtzLemmon2005::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzLemmon2005",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), 
                     _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), 
                     _m(rapidjson::kArrayType);
    for (unsigned int i=0;i<=N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        _n.PushBack((double)el.n, doc.GetAllocator());
        _d.PushBack((double)el.d, doc.GetAllocator());
        _t.PushBack((double)el.t, doc.GetAllocator());
        _l.PushBack((double)el.l, doc.GetAllocator());
        _m.PushBack((double)el.m, doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
    el.AddMember("m",_m,doc.GetAllocator());
}

// Term and its derivatives
long double ResidualHelmholtzLemmon2005::base(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0)
            s[i] = ni*exp(ti*log_tau+di*log_delta-pow(delta, li)-pow(tau, mi));
        else if (li != 0 && mi == 0)
            s[i] = ni*exp(ti*log_tau+di*log_delta-pow(delta, li));
        else
            s[i] = ni*exp(ti*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta,li);
            long double pow_tau_mi = pow(tau,mi);
            s[i] = ni*(di-li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li>0 && mi == 0){
            long double pow_delta_li = pow(delta,li);
            s[i] = ni*(di-li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*exp(ti*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_tau_mi = pow(tau, mi);
            s[i] = ni*(ti-mi*pow_tau_mi)*exp((ti-1)*log_tau+di*log_delta-pow(delta,li)-pow_tau_mi);
        }
        else if (li != 0 && mi == 0)
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta-pow(delta,li));
        else
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){	
            long double pow_delta_li = pow(delta, li);
            long double pow_tau_mi = pow(tau, mi);
            s[i] = ni*((di-li*pow_delta_li)*(di-1.0-li*pow_delta_li) - li*li*pow_delta_li)*exp(ti*log_tau+(di-2)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0){
            long double pow_delta_li = pow(delta, li);
            s[i] = ni*((di-li*pow_delta_li)*(di-1.0-li*pow_delta_li) - li*li*pow_delta_li)*exp(ti*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*exp(ti*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta, li);
            long double pow_tau_mi = pow(tau, mi);
            s[i] = ni*(di-li*pow_delta_li)*(ti-mi*pow_tau_mi)*exp((ti-1)*log_tau+(di-1)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0){
            long double pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(di-li*pow_delta_li)*exp((ti-1)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*exp((ti-1)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_tau_mi = pow(tau, mi);
            long double bracket = (ti-mi*pow_tau_mi)*(ti-1-mi*pow_tau_mi)-mi*mi*pow_tau_mi;
            s[i] = ni*bracket*exp((ti-2)*log_tau+di*log_delta-pow(delta,li)-pow_tau_mi);
        }
        else if (li != 0 && mi == 0)
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta-pow(delta,li));
        else
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta, li);
            long double pow_tau_mi = pow(tau, mi);
            long double bracket = (di*(di-1)*(di-2)+pow_delta_li*(-2*li+6*di*li-3*di*di*li-3*di*li*li+3*li*li-li*li*li)+pow_delta_li*pow_delta_li*(3*di*li*li-3*li*li+3*li*li*li)-li*li*li*pow_delta_li*pow_delta_li*pow_delta_li);
            s[i] = ni*bracket*exp(ti*log_tau+(di-3)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0)
        {
            long double pow_delta_li = pow(delta,li);
            long double bracket = (di*(di-1)*(di-2)+pow_delta_li*(-2*li+6*di*li-3*di*di*li-3*di*li*li+3*li*li-li*li*li)+pow_delta_li*pow_delta_li*(3*di*li*li-3*li*li+3*li*li*li)-li*li*li*pow_delta_li*pow_delta_li*pow_delta_li);
            s[i] = ni*bracket*exp(ti*log_tau+(di-3)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*(di-2)*exp(ti*log_tau+(di-3)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta,li);
            long double pow_tau_mi = pow(tau,mi);
            long double bracket = (ti-mi*pow_tau_mi)*(((di-li*pow_delta_li))*(di-1-li*pow_delta_li)-li*li*pow_delta_li);
            s[i] = ni*bracket*exp((ti-1)*log_tau+(di-2)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0){
            long double pow_delta_li = pow(delta,li);
            long double bracket = ti*(((di-li*pow_delta_li))*(di-1-li*pow_delta_li)-li*li*pow_delta_li);
            s[i] = ni*bracket*exp((ti-1)*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*(di-1)*exp((ti-1)*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzLemmon2005::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta,li);
            long double pow_tau_mi = pow(tau,mi);
            // delta derivative of second tau derivative
            long double bracket = ((ti-mi*pow_tau_mi)*(ti-1-mi*pow_tau_mi)-mi*mi*pow_tau_mi)*(di-li*pow_delta_li);
            s[i] = ni*bracket*exp((ti-2)*log_tau+(di-1)*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0){
            long double pow_delta_li = pow(delta,li);
            s[i] = ni*ti*(ti-1)*(di-li*pow_delta_li)*exp((ti-2)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*ti*(ti-1)*di*exp((ti-2)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

/**
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}{\tau ^{{t_k} - 2}}\exp \left( { - {\delta ^{{l_k}}}} \right)\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}\exp \left( { - {\delta ^{{l_k}}}} \right){\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
Group all the terms that don't depend on \f$ \tau \f$
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = A{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] + \frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right]\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right] = ({t_k} - 2){\tau ^{{t_k} - 3}}\exp \left( { - {\tau ^{{m_k}}}} \right) + {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)( - {m_k}{\tau ^{{m_k} - 1}}) = \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] = \left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}} \right) + \left( { - m_k^2{\tau ^{{m_k} - 1}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^3{\tau ^{{m_k} - 1}} =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {{t_k} - {m_k}{\tau ^{{m_k}}} + {t_k} - 1 - {m_k}{\tau ^{{m_k}}} + {m_k}} \right] =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]} \right) + \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = \exp \left( { - {\tau ^{{m_k}}}} \right)\left[ { - {\tau ^{{t_k} - 2}}m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right] + \left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]} \right]
\f]
*/
long double ResidualHelmholtzLemmon2005::dTau3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    long double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i=0; i<N; ++i)
    {
        ResidualHelmholtzLemmon2005Element &el = elements[i];
        long double ni = el.n, ti = el.t, di = el.d;
        int li = el.l, mi = el.m;
        if (li != 0 && mi != 0){
            long double pow_delta_li = pow(delta,li);
            long double pow_tau_mi = pow(tau,mi);
            //long double bracket = -pow(tau,ti+mi-3)*mi*mi*(2*ti-2*mi*pow_tau_mi-1-mi)+((ti-2)*pow(tau,ti-3)-pow(tau,ti-2)*mi*pow(tau,mi-1))*((ti-mi*pow_tau_mi)*(ti-1-mi*pow_tau_mi)-mi*mi*pow_tau_mi);
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta-pow_delta_li-pow_tau_mi);
        }
        else if (li != 0 && mi == 0){
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta-pow(delta,li));
        }
        else
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

void ResidualHelmholtzNonAnalytic::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzNonAnalytic",doc.GetAllocator());

    rapidjson::Value _n(rapidjson::kArrayType), _a(rapidjson::kArrayType), _b(rapidjson::kArrayType), 
                     _beta(rapidjson::kArrayType), _A(rapidjson::kArrayType), _B(rapidjson::kArrayType), 
                     _C(rapidjson::kArrayType), _D(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        _n.PushBack((double)el.n, doc.GetAllocator());
        _a.PushBack((double)el.a, doc.GetAllocator());
        _b.PushBack((double)el.b, doc.GetAllocator());
        _beta.PushBack((double)el.beta, doc.GetAllocator());
        _A.PushBack((double)el.A, doc.GetAllocator());
        _B.PushBack((double)el.B, doc.GetAllocator());
        _C.PushBack((double)el.C, doc.GetAllocator());
        _D.PushBack((double)el.D, doc.GetAllocator());
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

long double ResidualHelmholtzNonAnalytic::base(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        
        s[i] = ni*pow(DELTA, bi)*delta*PSI;
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dDelta(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
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
long double ResidualHelmholtzNonAnalytic::dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        s[i] = ni*delta*(dDELTAbi_dTau*PSI+pow(DELTA,bi)*dPSI_dTau);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

long double ResidualHelmholtzNonAnalytic::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
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
long double ResidualHelmholtzNonAnalytic::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;
        
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
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
long double ResidualHelmholtzNonAnalytic::dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        s[i] = ni*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}

long double ResidualHelmholtzNonAnalytic::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        long double dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        long double dPSI3_dDelta3=2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        long double PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
        long double dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
        long double dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;        
        long double dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
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
long double ResidualHelmholtzNonAnalytic::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        
        long double dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        long double dPSI3_dDelta_dTau2 = 2*Di*(2*Di*pow(tau-1,2)-1)*dPSI_dDelta;
        long double dDELTAbi3_dDelta_dTau2 = 2*bi*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+4*pow(theta,2)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta+8*theta*bi*(bi-1)*pow(DELTA,bi-2)*dtheta_dDelta;

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
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

long double ResidualHelmholtzNonAnalytic::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        //Following Terms added for this derivative
        long double dPSI3_dDelta2_dTau = (2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*dPSI_dTau;
        long double dDELTA_dTau = -2*theta;
        long double dDELTA2_dDelta_dTau = 2.0*Ai/(betai)*pow(pow(delta-1,2),1.0/(2.0*betai)-0.5);
        long double dDELTA3_dDelta2_dTau = 2.0*Ai*(betai-1)/(betai*betai)*pow(pow(delta-1,2),1/(2*betai)-1.0);
        
        long double dDELTAbim1_dTau = (bi-1)*pow(DELTA,bi-2)*dDELTA_dTau;
        long double dDELTAbim2_dTau = (bi-2)*pow(DELTA,bi-3)*dDELTA_dTau;
        long double Line11 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,bi-1)*dDELTA3_dDelta2_dTau;
        long double Line21 = (bi-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
        long double dDELTAbi3_dDelta2_dTau = bi*(Line11+Line21);

        long double Line1 = pow(DELTA,bi)*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
        long double Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
        long double Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        s[i] = ni*(Line1+Line2+Line3);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dTau3(const long double &tau, const long double &delta) throw()
{
    if (N==0){return 0.0;}
    for (unsigned int i=0; i<N; ++i)
    {
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);
        long double dPSI3_dTau3=2.0*Di*PSI*(-4*Di*Di*pow(tau-1,3)+6*Di*(tau-1));
        long double dDELTAbi3_dTau3 = -12.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2)-8*pow(theta,3)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3.0);

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
long double ResidualHelmholtzSAFTAssociating::Deltabar(const long double &tau, const long double &delta)
{
    return this->g(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar;
}   
long double ResidualHelmholtzSAFTAssociating::dDeltabar_ddelta__consttau(const long double &tau, const long double &delta)
{
    return this->dg_deta(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*this->vbarn;
}
long double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta2__consttau(const long double &tau, const long double &delta)
{
    return this->d2g_deta2(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)2);
}
long double ResidualHelmholtzSAFTAssociating::dDeltabar_dtau__constdelta(const long double &tau, const long double &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*this->epsilonbar;
}
long double ResidualHelmholtzSAFTAssociating::d2Deltabar_dtau2__constdelta(const long double &tau, const long double &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2);
}
long double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta_dtau(const long double &tau, const long double &delta)
{
    return this->dg_deta(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*this->vbarn;
}
long double ResidualHelmholtzSAFTAssociating::d3Deltabar_dtau3__constdelta(const long double &tau, const long double &delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)3);
}
long double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta_dtau2(const long double &tau, const long double &delta)
{
    return this->dg_deta(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2)*this->vbarn;
}
long double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta2_dtau(const long double &tau, const long double &delta)
{
    return this->d2g_deta2(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*pow(this->vbarn,(int)2);
}
long double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta3__consttau(const long double &tau, const long double &delta)
{
    return this->d3g_deta3(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)3);
}

long double ResidualHelmholtzSAFTAssociating::X(const long double &delta, const long double &Deltabar)
{
    return 2/(sqrt(1+4*Deltabar*delta)+1);
}
long double ResidualHelmholtzSAFTAssociating::dX_dDeltabar__constdelta(const long double &delta, const long double &Deltabar)
{
    long double X = this->X(delta,Deltabar);
    return -delta*X*X/(2*Deltabar*delta*X+1);
}
long double ResidualHelmholtzSAFTAssociating::dX_ddelta__constDeltabar(const long double &delta, const long double &Deltabar)
{
    long double X = this->X(delta,Deltabar);
    return -Deltabar*X*X/(2*Deltabar*delta*X+1);
}
long double ResidualHelmholtzSAFTAssociating::dX_dtau(const long double &tau, const long double &delta)
{
    long double Deltabar = this->Deltabar(tau, delta);
    return this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_dtau__constdelta(tau, delta);
}
long double ResidualHelmholtzSAFTAssociating::dX_ddelta(const long double &tau, const long double &delta)
{
    long double Deltabar = this->Deltabar(tau, delta);
    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
           + this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_ddelta__consttau(tau, delta));
}
long double ResidualHelmholtzSAFTAssociating::d2X_dtau2(const long double &tau, const long double &delta)
{
    long double Deltabar = this->Deltabar(tau, delta);
    long double X = this->X(delta, Deltabar);
    long double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    long double d_dXdtau_dbeta = -delta*X*X/(2*Deltabar*delta*X+1);
    long double d_dXdtau_dDeltabar = 2*delta*delta*X*X*X/pow(2*Deltabar*delta*X+1,2)*beta;
    long double d_dXdtau_dX = -2*beta*delta*X*(Deltabar*delta*X+1)/pow(2*Deltabar*delta*X+1,2);
    long double dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
    long double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXdtau_dX*dX_dDeltabar*beta+d_dXdtau_dDeltabar*beta+d_dXdtau_dbeta*dbeta_dtau;
}
long double ResidualHelmholtzSAFTAssociating::d2X_ddeltadtau(const long double &tau, const long double &delta)
{
    long double Deltabar = this->Deltabar(tau, delta);
    long double X = this->X(delta, Deltabar);
    long double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    long double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    long double dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
    long double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    long double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    long double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    long double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXddelta_dX*dX_dDeltabar*beta+d_dXddelta_dDeltabar*beta+d_dXddelta_dalpha*dalpha_dtau;
}
long double ResidualHelmholtzSAFTAssociating::d2X_ddelta2(const long double &tau, const long double &delta)
{
    long double Deltabar = this->Deltabar(tau, delta);
    long double X = this->X(delta, Deltabar);
    long double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    long double dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);
    
    long double dX_ddelta_constall = X*X*(2*Deltabar*Deltabar*X-alpha)/pow(2*Deltabar*delta*X+1,2);
    long double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    long double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    long double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    
    long double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    long double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Deltabar);

    long double val = (dX_ddelta_constall
            + d_dXddelta_dX*dX_ddelta
            + d_dXddelta_dX*dX_dDeltabar*alpha
            + d_dXddelta_dDeltabar*alpha
            + d_dXddelta_dalpha*dalpha_ddelta);
    return val;
}   
long double ResidualHelmholtzSAFTAssociating::d3X_dtau3(const long double &tau, const long double &delta)
{
    long double Delta = this->Deltabar(tau, delta);
    long double X = this->X(delta, Delta);
    long double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    long double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    long double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    long double Delta_ttt = this->d3Deltabar_dtau3__constdelta(tau, delta);
    long double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    long double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    long double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    long double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_dX*dX_dDelta*Delta_t+dXtt_dDelta*Delta_t + dXtt_dDelta_t*Delta_tt + dXtt_dDelta_tt*Delta_ttt;
}
long double ResidualHelmholtzSAFTAssociating::d3X_ddeltadtau2(const long double &tau, const long double &delta)
{
    long double Delta = this->Deltabar(tau, delta);
    long double X = this->X(delta, Delta);
    long double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    long double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    long double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    long double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    long double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    long double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    long double Delta_dtt = this->d3Deltabar_ddelta_dtau2(tau, delta);
    long double dXtt_ddelta = pow(X, 2)*(-12*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 2*pow(Delta_t, 2)*X*delta*(-Delta*X*delta + 2)*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + 2*X*delta*(Delta*Delta_tt + 2*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    long double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    long double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    long double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    long double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_ddelta + dXtt_dX*dX_ddelta + dXtt_dX*dX_dDelta*Delta_d + dXtt_dDelta*Delta_d + dXtt_dDelta_t*Delta_dt + dXtt_dDelta_tt*Delta_dtt;
}

long double ResidualHelmholtzSAFTAssociating::d3X_ddelta2dtau(const long double &tau, const long double &delta)
{
    long double Delta = this->Deltabar(tau, delta);
    long double X = this->X(delta, Delta);
    long double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    long double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    long double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    long double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    long double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    long double Delta_ddt = this->d3Deltabar_ddelta2_dtau(tau, delta);
    long double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    long double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    long double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    long double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_dX*dX_dDelta*Delta_t + dXdd_dDelta*Delta_t + dXdd_dDelta_d*Delta_dt + dXdd_dDelta_dd*Delta_ddt;
}

double Xdd(double X, double delta, double Delta, double Delta_d, double Delta_dd)
{
    return Delta*pow(X, 2)*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*delta*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*(2*Delta_d*X*pow(delta, 2) - 1)/pow(2*Delta*X*delta + 1, 2) - Delta_dd*pow(X, 2)*delta/(2*Delta*X*delta + 1) + pow(X, 2)*(2*pow(Delta, 2)*X - Delta_d)/pow(2*Delta*X*delta + 1, 2);
}

long double ResidualHelmholtzSAFTAssociating::d3X_ddelta3(const long double &tau, const long double &delta)
{
    long double Delta = this->Deltabar(tau, delta);
    long double X = this->X(delta, Delta);
    long double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    long double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    long double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    long double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    long double Delta_ddd = this->d3Deltabar_ddelta3__consttau(tau, delta);

    long double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    long double dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    long double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    long double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    long double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_ddelta + dXdd_dX*(dX_ddelta + dX_dDelta*Delta_d) + dXdd_dDelta*Delta_d + dXdd_dDelta_d*Delta_dd + dXdd_dDelta_dd*Delta_ddd;
}
long double ResidualHelmholtzSAFTAssociating::g(const long double &eta)
{
    return 0.5*(2-eta)/pow(1-eta,(int)3);
}    
long double ResidualHelmholtzSAFTAssociating::dg_deta(const long double &eta)
{
    return 0.5*(5-2*eta)/pow(1-eta,(int)4);
}
long double ResidualHelmholtzSAFTAssociating::d2g_deta2(const long double &eta)
{
    return 3*(3-eta)/pow(1-eta,(int)5);
}   
long double ResidualHelmholtzSAFTAssociating::d3g_deta3(const long double &eta)
{
    return 6*(7-2*eta)/pow(1-eta,(int)6);
}   
long double ResidualHelmholtzSAFTAssociating::eta(const long double &delta){
    return this->vbarn*delta;
}
long double ResidualHelmholtzSAFTAssociating::base(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*((log(X)-X/2.0+0.5));
}
long double ResidualHelmholtzSAFTAssociating::dDelta(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*(1/X-0.5)*this->dX_ddelta(tau, delta);
}
long double ResidualHelmholtzSAFTAssociating::dTau(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*(1/X-0.5)*this->dX_dtau(tau, delta);
}
long double ResidualHelmholtzSAFTAssociating::dTau2(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_tau = this->dX_dtau(tau, delta);
    long double X_tautau = this->d2X_dtau2(tau, delta);
    return this->m*this->a*((1/X-0.5)*X_tautau-pow(X_tau/X, 2));
}
long double ResidualHelmholtzSAFTAssociating::dDelta2(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_delta = this->dX_ddelta(tau, delta);
    long double X_deltadelta = this->d2X_ddelta2(tau, delta);
    return this->m*this->a*((1/X-0.5)*X_deltadelta-pow(X_delta/X,2));
}
long double ResidualHelmholtzSAFTAssociating::dDelta_dTau(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_delta = this->dX_ddelta(tau, delta);
    long double X_tau = this->dX_dtau(tau, delta);
    long double X_deltatau = this->d2X_ddeltadtau(tau, delta);
    return this->m*this->a*((-X_tau/X/X)*X_delta+X_deltatau*(1/X-0.5));
}
long double ResidualHelmholtzSAFTAssociating::dTau3(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_t = this->dX_dtau(tau, delta);
    long double X_tt = this->d2X_dtau2(tau, delta);
    long double X_ttt = this->d3X_dtau3(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ttt+(-X_t/pow(X,(int)2))*X_tt-2*(pow(X,(int)2)*(X_t*X_tt)-pow(X_t,(int)2)*(X*X_t))/pow(X,(int)4));
}
long double ResidualHelmholtzSAFTAssociating::dDelta_dTau2(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_t = this->dX_dtau(tau, delta);
    long double X_d = this->dX_ddelta(tau, delta);
    long double X_tt = this->d2X_dtau2(tau, delta);
    long double X_dt = this->d2X_ddeltadtau(tau, delta);
    long double X_dtt = this->d3X_ddeltadtau2(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_dtt-X_d/pow(X,(int)2)*X_tt-2*(pow(X,(int)2)*(X_t*X_dt)-pow(X_t,(int)2)*(X*X_d))/pow(X,(int)4));
}
long double ResidualHelmholtzSAFTAssociating::dDelta2_dTau(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_t = this->dX_dtau(tau, delta);
    long double X_d = this->dX_ddelta(tau, delta);
    long double X_dd = this->d2X_ddelta2(tau, delta);
    long double X_dt = this->d2X_ddeltadtau(tau, delta);
    long double X_ddt = this->d3X_ddelta2dtau(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ddt-X_t/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dt)-pow(X_d,(int)2)*(X*X_t))/pow(X,(int)4));
}
long double ResidualHelmholtzSAFTAssociating::dDelta3(const long double &tau, const long double &delta) throw()
{
    if (disabled){return 0;}
    long double X = this->X(delta, this->Deltabar(tau, delta));
    long double X_d = this->dX_ddelta(tau, delta);
    long double X_dd = this->d2X_ddelta2(tau, delta);
    long double X_ddd = this->d3X_ddelta3(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ddd-X_d/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dd)-pow(X_d,(int)2)*(X*X_d))/pow(X,(int)4));
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
long double IdealHelmholtzCP0PolyT::base(const long double &tau, const long double &delta) throw()
{ 
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (fabs(t[i])<10*DBL_EPSILON)
        {
            sum += c[i]-c[i]*tau/tau0+c[i]*log(tau/tau0);
        }
        else if (fabs(t[i]+1) < 10*DBL_EPSILON)
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
long double IdealHelmholtzCP0PolyT::dTau(const long double &tau, const long double &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (fabs(t[i])<10*DBL_EPSILON)
        {
            sum += c[i]/tau-c[i]/tau0;
        }
        else if (fabs(t[i]+1) < 10*DBL_EPSILON)
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
long double IdealHelmholtzCP0PolyT::dTau2(const long double &tau, const long double &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (fabs(t[i])<10*DBL_EPSILON)
        {
            sum += -c[i]/(tau*tau);
        }
        else if (fabs(t[i]+1) < 10*DBL_EPSILON)
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
long double IdealHelmholtzCP0PolyT::dTau3(const long double &tau, const long double &delta) throw()
{
    double sum=0;
    for (std::size_t i = 0; i < N; ++i){
        if (fabs(t[i])<10*DBL_EPSILON)
        {
            sum += 2*c[i]/(tau*tau*tau);
        }
        else if (fabs(t[i]+1) < 10*DBL_EPSILON)
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
long double IdealHelmholtzCP0AlyLee::base(const long double &tau, const long double &delta) throw()
{	
    if (!enabled){ return 0.0;}
    return -tau*(anti_deriv_cp0_tau2(tau)-anti_deriv_cp0_tau2(tau0))+(anti_deriv_cp0_tau(tau)-anti_deriv_cp0_tau(tau0));
}
long double IdealHelmholtzCP0AlyLee::dTau(const long double &tau, const long double &delta) throw()
{
    if (!enabled){ return 0.0;}
    return -(anti_deriv_cp0_tau2(tau) - anti_deriv_cp0_tau2(tau0));
}
long double IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau2(const long double &tau)
{
    return -c[0]/tau  +  2*c[1]*c[2]/(Tc*(exp(-(2*c[2]*tau)/Tc)-1))  +  2*c[3]*c[4]/(Tc*(exp(-(2*c[4]*tau)/Tc)+1));
}
long double IdealHelmholtzCP0AlyLee::anti_deriv_cp0_tau(const long double &tau)
{
    long double term1 = c[0]*log(tau);
    long double term2 = 2*c[1]*c[2]*tau/(-Tc + Tc*exp(-2*c[2]*tau/Tc)) + c[1]*log(-1 + exp(-2*c[2]*tau/Tc)) + 2*c[1]*c[2]*tau/Tc;
    long double term3 = -c[3]*(Tc*exp(2*c[4]*tau/Tc)*log(exp(2*c[4]*tau/Tc) + 1) + Tc*log(exp(2*c[4]*tau/Tc) + 1) - 2*c[4]*tau*exp(2*c[4]*tau/Tc))/(Tc*(exp(2*c[4]*tau/Tc) + 1));
    return term1 + term2 + term3;
}
long double IdealHelmholtzCP0AlyLee::dTau2(const long double &tau, const long double &delta) throw()
{
    if (!enabled){ return 0.0;}
    return -c[0]/pow(tau,2) - c[1]*pow(c[2]/Tc/sinh(c[2]*tau/Tc),2) - c[3]*pow(c[4]/Tc/cosh(c[4]*tau/Tc),2);
}
long double IdealHelmholtzCP0AlyLee::dTau3(const long double &tau, const long double &delta) throw()
{
    if (!enabled){ return 0.0;}
    return 2*c[0]/pow(tau,3) + 2*c[1]*pow(c[2]/Tc,3)*cosh(c[2]*tau/Tc)/pow(sinh(c[2]*tau/Tc),3) + 2*c[3]*pow(c[4]/Tc,3)*sinh(c[4]*tau/Tc)/pow(cosh(c[4]*tau/Tc),3);
}

}; /* namespace CoolProp */


/*
IdealHelmholtzEnthalpyEntropyOffset EnthalpyEntropyOffset;
IdealHelmholtzPlanckEinstein PlanckEinstein;
IdealHelmholtzPlanckEinstein2 PlanckEinstein2;

IdealHelmholtzCP0Constant CP0Constant;
IdealHelmholtzCP0PolyT CP0PolyT;
IdealHelmholtzCP0AlyLee CP0AlyLee;
*/

#ifdef ENABLE_CATCH
#include <math.h>
#include "catch.hpp"

class HelmholtzConsistencyFixture
{
public:
    long double numerical, analytic;
    CoolProp::BaseHelmholtzTerm *Lead, *LogTau, *IGPower, *PlanckEinstein, *PlanckEinstein2;

    HelmholtzConsistencyFixture(){
        Lead = new CoolProp::IdealHelmholtzLead(1,3);
        LogTau = new CoolProp::IdealHelmholtzLogTau(1.5);
        {
            std::vector<long double> n(4,0), t(4,1); n[0] = -0.1; n[2] = 0.1; t[1] = -1; t[2] = -2; t[3] = 2;
            IGPower = new CoolProp::IdealHelmholtzPower(n,t);
        }
        {
            std::vector<long double> n(4,0), t(4,1); n[0] = -0.1; n[2] = 0.1; t[1] = -1; t[2] = -2; t[3] = 2;
            PlanckEinstein = new CoolProp::IdealHelmholtzPlanckEinstein(n, t);
        }
        {
            std::vector<long double> n(4,0), t(4,1), c(4,1); n[0] = -0.1; n[2] = 0.1; t[1] = -1; t[2] = -2; t[3] = 2;
            PlanckEinstein2 = new CoolProp::IdealHelmholtzPlanckEinstein2(n, t, c);
        }
    }
    void call(std::string d, CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double ddelta)
    {
              if (!d.compare("dTau")){return dTau(term,tau,delta,ddelta);}
        else if (!d.compare("dTau2")){return dTau2(term,tau,delta,ddelta);}
        else if (!d.compare("dTau3")){return dTau3(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta")){ return dDelta(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta2")){return dDelta2(term,tau,delta,ddelta);}
        else if (!d.compare("dDelta3")){return dDelta3(term,tau,delta,ddelta);}
        else{
            throw CoolProp::ValueError("don't understand deriv type");
        }
    }
    CoolProp::BaseHelmholtzTerm * get(std::string t)
    {
        if (!t.compare("Lead")){return Lead;}
        else if (!t.compare("LogTau")){return LogTau;}
        else if (!t.compare("IGPower")){return IGPower;}
        else if (!t.compare("PlanckEinstein")){return PlanckEinstein;}
        else if (!t.compare("PlanckEinstein2")){return PlanckEinstein2;}
        else{
            throw CoolProp::ValueError("don't understand helmholtz type");
        }
    }
    void dTau(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double dtau){
        long double term_plus = term->base(tau + dtau, delta);
        long double term_minus = term->base(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau(tau, delta);
    };
    void dTau2(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double dtau){
        long double term_plus = term->dTau(tau + dtau, delta);
        long double term_minus = term->dTau(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau2(tau, delta);
    };
    void dTau3(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double dtau){
        long double term_plus = term->dTau2(tau + dtau, delta);
        long double term_minus = term->dTau2(tau - dtau, delta);
        numerical = (term_plus - term_minus)/(2*dtau);
        analytic = term->dTau3(tau, delta);
    };
    void dDelta(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double ddelta){
        long double term_plus = term->base(tau, delta + ddelta);
        long double term_minus = term->base(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta(tau, delta);
    };
    void dDelta2(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double ddelta){
        long double term_plus = term->dDelta(tau, delta + ddelta);
        long double term_minus = term->dDelta(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta2(tau, delta);
    };
    void dDelta3(CoolProp::BaseHelmholtzTerm * term, long double tau, long double delta, long double ddelta){
        long double term_plus = term->dDelta2(tau, delta + ddelta);
        long double term_minus = term->dDelta2(tau, delta - ddelta);
        numerical = (term_plus - term_minus)/(2*ddelta);
        analytic = term->dDelta3(tau, delta);
    };
    double err(double v1, double v2)
    {
        if (fabs(v2) > 1e-15){
            return fabs((v1-v2)/v2);
        }
        else{
            return fabs(v1-v2);
        }
    }
};

std::string terms[] = {"Lead","LogTau","IGPower","PlanckEinstein","PlanckEinstein2"};
std::string derivs[] = {"dTau","dTau2","dTau3","dDelta","dDelta2","dDelta3"};

TEST_CASE_METHOD(HelmholtzConsistencyFixture, "Helmholtz energy derivatives", "[helmholtz]")
{
    CoolProp::BaseHelmholtzTerm * term;
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
            CHECK(err(analytic, numerical) < 1e-8);
        }
    }
}

#endif

