//
//  VTPRCubic.h
//  CoolProp
//
//  Created by Ian on 7/17/16.
//
//

#include "GeneralizedCubic.h"

#ifndef VTPRCubic_h
#define VTPRCubic_h

class VTPRCubic : public PengRobinson
{
private:
    UNIFAQ::UNIFAQMixture unifaq;
public:
    VTPRCubic(std::vector<double> Tc,
              std::vector<double> pc,
              std::vector<double> acentric,
              double R_u,
              const UNIFAQLibrary::UNIFAQParameterLibrary & lib
              )
    : PengRobinson(Tc, pc, acentric, R_u), unifaq(lib) {};
    
    VTPRCubic(double Tc,
              double pc,
              double acentric,
              double R_u,
              const UNIFAQLibrary::UNIFAQParameterLibrary & lib)
    : PengRobinson(std::vector<double>(1,Tc), std::vector<double>(1,pc), std::vector<double>(1,acentric), R_u), unifaq(lib) {};
    
    /// Get a reference to the managed UNIFAQ instance
    UNIFAQ::UNIFAQMixture &get_unifaq(){ return unifaq; }
    
    /// The attractive part in cubic EOS
    double a_alpha(double T, std::size_t i){
        return pow(1 + m_ii(i)*(1-sqrt(T/Tc[i])), 2);
    }
    /// Calculate the non-dimensionalized gE/RT term
    double gE_R_RT(){
        const std::vector<double> &z = unifaq.get_mole_fractions();
        double summer = 0;
        for (std::size_t i = 0; i < z.size(); ++i) {
            summer += z[i]*unifaq.ln_gamma_R(i);
        }
        return summer;
    }
    /// The co-volume for the i-th pure component
    double b_ii(std::size_t i){
        return 0.0778*R_u*Tc[i]/pc[i];
    }
    /// The attractive parameter for the i-th pure component
    double a_ii(std::size_t i){
        double a0 = 0.45724*pow(R_u*Tc[i], 2)/pc[i];
        double alpha = a_alpha(unifaq.get_temperature(), i);
        return a0*alpha;
    }
    double cm(){
        return 0;
    }
    double am_term(double tau, const std::vector<double> &x, std::size_t itau){
        if (itau == 0){
            set_temperature(T_r/tau);
            double _am,_bm; am_bm(_am, _bm);
            return _am;
        }
        else{
            throw CoolProp::NotImplementedError();
        }
    }
    double bm_term(const std::vector<double> &x){
        double _am,_bm; am_bm(_am, _bm);
        return _bm;
    }
    /// Calculate both am and bm because am and bm are dependent on each other
    void am_bm(double &am, double &bm){
        const std::vector<double> &z = unifaq.get_mole_fractions();
        double summeram = 0, summerbm = 0;
        for (std::size_t i = 0; i < z.size(); ++i){
            summeram += z[i]*a_ii(i)/b_ii(i);
            for (std::size_t j = 0; j < z.size(); ++j){
                summerbm += z[i]*z[j]*pow((pow(b_ii(i), 0.75) + pow(b_ii(j),0.75))/2.0, 4.0/3.0);
            }
        }
        bm = summerbm;
        am = bm*(summeram + gE_R_RT()/(-0.53087));
    };
    void set_temperature(const double T){ unifaq.set_temperature(T); }
    
    /// The residual non-dimensionalized Helmholtz energy \f$\alpha^r\f$
    double alphar(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta){
        unifaq.set_mole_fractions(x);
        unifaq.set_temperature(T_r/tau);
        double a_m,b_m; am_bm(a_m, b_m);
        double c_m = cm();
        if (itau == 0 && idelta == 0){
            return -log(1-delta*rho_r*(b_m-c_m)) - sqrt(2.0)*a_m*tau/(4*R_u*T_r*b_m)*log( (1+delta*rho_r*(b_m*(1+sqrt(2.0)+c_m))) / (1+delta*rho_r*(b_m*(1-sqrt(2.0)+c_m))) );
        }
        else if (itau == 1 && idelta == 0){
            double dtau = 0.01*tau;
            return (alphar(tau+dtau,delta,x,0,0) - alphar(tau-dtau,delta,x,0,0))/(2*dtau);
        }
        else if (itau == 2 && idelta == 0){
            double dtau = 0.01*tau;
            return (alphar(tau+dtau,delta,x,0,0) - 2*alphar(tau,delta,x,0,0) + alphar(tau-dtau,delta,x,0,0))/(dtau*dtau);
        }
        else if (itau == 1 && idelta == 1){
            double dtau = 0.01*tau, ddelta = 0.001*delta;
            return (alphar(tau+dtau,delta+ddelta,x,0,0) -alphar(tau-dtau,delta+ddelta,x,0,0) -alphar(tau+dtau,delta-ddelta,x,0,0) + alphar(tau-dtau,delta-ddelta,x,0,0))/(4*ddelta*dtau);
        }
        else if (itau == 0 && idelta == 1){
            double T= T_r/tau;
            double analytic = (1.0L/4.0L)*rho_r*(-4*R_u*T*(b_m - c_m)*(b_m*delta*rho_r*(c_m + 1 + sqrt(2.0)) + 1)*(b_m*delta*rho_r*(c_m - sqrt(2.0) + 1) + 1) + sqrt(2.0)*a_m*((b_m*delta*rho_r*(c_m + 1 + sqrt(2.0)) + 1)*(c_m - sqrt(2.0) + 1) - (b_m*delta*rho_r*(c_m - sqrt(2.0) + 1) + 1)*(c_m + 1 + sqrt(2.0)))*(delta*rho_r*(b_m - c_m) - 1))/(R_u*T*(delta*rho_r*(b_m - c_m) - 1)*(b_m*delta*rho_r*(c_m + 1 + sqrt(2.0)) + 1)*(b_m*delta*rho_r*(c_m - sqrt(2.0) + 1) + 1));
            double ddelta = 0.001*delta;
            double val = (alphar(tau,delta+ddelta,x,0,0) - alphar(tau,delta-ddelta,x,0,0))/(2*ddelta);
            return analytic;
        }
        else if (itau == 0 && idelta == 2){
            double ddelta = 0.001*delta;
            double val = (alphar(tau,delta+ddelta,x,0,0) - 2*alphar(tau,delta,x,0,0) + alphar(tau,delta-ddelta,x,0,0))/(ddelta*ddelta);
            return val;
        }
        else{
            return _HUGE;
            throw CoolProp::NotImplementedError();
        }
    };
    /// The first composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d_alphar_dxi(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, bool xN_independent){
        double dx = 1e-6;
        std::vector<double> xp = x, xm = x;
        xp[i] += dx; xm[i] -= dx;
        if (!xN_independent){ xp[xp.size()-1] -= dx; xm[xm.size()-1] += dx; }
        double val = (alphar(tau,delta,xp,itau,idelta) - alphar(tau,delta,xm,itau,idelta))/(2*dx);
        return val;
    };
    /// The second composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d2_alphar_dxidxj(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, bool xN_independent){
        double dx = 1e-6;
        std::vector<double> xp = x, xm = x;
        xp[i] += dx; xm[i] -= dx;
        if (!xN_independent){ xp[xp.size()-1] -= dx; xm[xm.size()-1] += dx; }
        double val = (d_alphar_dxi(tau,delta,xp,itau,idelta, j,xN_independent) - d_alphar_dxi(tau,delta,xm,itau,idelta,j,xN_independent))/(2*dx);
        return val;
    };
    /// The third composition derivative of \f$\alpha^r\f$ as well as derivatives with respect to \f$\tau\f$ and \f$\delta\f$
    double d3_alphar_dxidxjdxk(double tau, double delta, const std::vector<double> &x, std::size_t itau, std::size_t idelta, std::size_t i, std::size_t j, std::size_t k, bool xN_independent){
        throw CoolProp::NotImplementedError();
    };
};

#endif /* VTPRCubic_h */
