#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <memory>

#include "HumidAirProp.h"
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Solvers.h"
#include "CoolPropTools.h"
#include "Ice.h"
#include "CoolProp.h"
#include "crossplatform_shared_ptr.h"
#include "Exceptions.h"

#include <algorithm>    // std::next_permutation
#include <stdlib.h>
#include "math.h"
#include "time.h"
#include "stdio.h"
#include <string.h>
#include <iostream>
#include <list>

/// This is a stub overload to help with all the strcmp calls below and avoid needing to rewrite all of them
std::size_t strcmp(const std::string &s, const std::string e){
    return s.compare(e);
}

// This is a lazy stub function to avoid recoding all the strcpy calls below
void strcpy(std::string &s, const std::string e){
    s = e;
}

shared_ptr<CoolProp::HelmholtzEOSBackend> Water, Air;

namespace HumidAir
{


void check_fluid_instantiation()
{
    if (!Water.get()){
        Water.reset(new CoolProp::HelmholtzEOSBackend("Water"));
    }
    if (!Air.get()){
        Air.reset(new CoolProp::HelmholtzEOSBackend("Air"));
    }
};

enum givens{GIVEN_INVALID=0, GIVEN_TDP,GIVEN_PSIW, GIVEN_HUMRAT,GIVEN_VDA, GIVEN_VHA,GIVEN_TWB,GIVEN_RH,GIVEN_ENTHALPY,GIVEN_ENTHALPY_HA,GIVEN_ENTROPY,GIVEN_ENTROPY_HA, GIVEN_T,GIVEN_P,GIVEN_VISC,GIVEN_COND,GIVEN_CP,GIVEN_CPHA};

static double epsilon=0.621945,R_bar=8.314472;
static int FlagUseVirialCorrelations=0,FlagUseIsothermCompressCorrelation=0,FlagUseIdealGasEnthalpyCorrelations=0;
double f_factor(double T, double p);

// A couple of convenience functions that are needed quite a lot
static double MM_Air(void)
{
    check_fluid_instantiation();
    return Air->keyed_output(CoolProp::imolar_mass);
}
static double MM_Water(void)
{
    check_fluid_instantiation();
    return Water->keyed_output(CoolProp::imolar_mass);
}
static double B_Air(double T)
{
    check_fluid_instantiation();
    Air->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Air->keyed_output(CoolProp::iBvirial);
}
static double dBdT_Air(double T)
{
    check_fluid_instantiation();
    Air->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Air->keyed_output(CoolProp::idBvirial_dT);
}
static double B_Water(double T)
{
    check_fluid_instantiation();
    Water->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Water->keyed_output(CoolProp::iBvirial);
}
static double dBdT_Water(double T)
{
    check_fluid_instantiation();
    Water->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Water->keyed_output(CoolProp::idBvirial_dT);
}
static double C_Air(double T)
{
    check_fluid_instantiation();
    Air->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Air->keyed_output(CoolProp::iCvirial);
}
static double dCdT_Air(double T)
{
    check_fluid_instantiation();
    Air->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Air->keyed_output(CoolProp::idCvirial_dT);
}
static double C_Water(double T)
{
    check_fluid_instantiation();
    Water->update(CoolProp::DmolarT_INPUTS,1e-20,T);
    return Water->keyed_output(CoolProp::iCvirial);
}
static double dCdT_Water(double T)
{
    check_fluid_instantiation();
    Water->update(CoolProp::DmolarT_INPUTS,1e-12,T);
    return Water->keyed_output(CoolProp::idCvirial_dT);
}
void UseVirialCorrelations(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseVirialCorrelations=flag;
    }
    else
    {
        printf("UseVirialCorrelations takes an integer, either 0 (no) or 1 (yes)\n");
    }

}
void UseIsothermCompressCorrelation(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseIsothermCompressCorrelation=flag;
    }
    else
    {
        printf("UseIsothermCompressCorrelation takes an integer, either 0 (no) or 1 (yes)\n");
    }
}
void UseIdealGasEnthalpyCorrelations(int flag)
{
    if (flag==0 || flag==1)
    {
        FlagUseIdealGasEnthalpyCorrelations=flag;
    }
    else
    {
        printf("UseIdealGasEnthalpyCorrelations takes an integer, either 0 (no) or 1 (yes)\n");
    }
}
static double Brent_HAProps_T(const std::string &OutputName, const std::string &Input1Name, double Input1, const std::string &Input2Name, double Input2, double TargetVal, double T_min, double T_max)
{
    double T;
    class BrentSolverResids : public CoolProp::FuncWrapper1D
    {
    private:
        double Input1,Input2,TargetVal;
        std::string OutputName, Input1Name, Input2Name;
    public:
        BrentSolverResids(std::string OutputName, std::string Input1Name, double Input1, std::string Input2Name, double Input2, double TargetVal)
        {
            this->OutputName = OutputName;
            this->Input1Name = Input1Name;
            this->Input2Name = Input2Name;
            this->Input1 = Input1;
            this->Input2 = Input2;
            this->TargetVal = TargetVal;
        };
        ~BrentSolverResids(){};

        double call(double T){
            double val = HAPropsSI(OutputName, "T", T, Input1Name, Input1, Input2Name, Input2);
            return val - TargetVal;
        }
    };

    BrentSolverResids BSR = BrentSolverResids(OutputName, Input1Name, Input1, Input2Name, Input2, TargetVal);

    std::string errstr;
    // Now we need to check the bounds and make sure that they are ok (don't yield invalid output)
    // and actually bound the solution
    double r_min = BSR.call(T_min);
    bool T_min_valid = ValidNumber(r_min);
    double r_max = BSR.call(T_max);
    bool T_max_valid = ValidNumber(r_max);
    if (!T_min_valid && !T_max_valid){
        throw CoolProp::ValueError(format("Both T_min [%g] and T_max [%g] yield invalid output values in Brent_HAProps_T",T_min,T_max).c_str());
    }
    else if (T_min_valid && !T_max_valid){
        while (!T_max_valid){
            // Reduce T_max until it works
            T_max = 0.95*T_max + 0.05*T_min;
            r_max = BSR.call(T_max);
            T_max_valid = ValidNumber(r_max);
        }
    }
    else if (!T_min_valid && T_max_valid){
        while (!T_min_valid){
            // Increase T_min until it works
            T_min = 0.95*T_min + 0.05*T_max;
            r_min = BSR.call(T_min);
            T_min_valid = ValidNumber(r_min);
        }
    }
    // We will do a secant call if the values at T_min and T_max have the same sign
    if (r_min*r_max > 0){
        if (std::abs(r_min) < std::abs(r_max)){
            T = CoolProp::Secant(BSR, T_min, 0.01*T_min, 1e-7, 50, errstr);
        }
        else{
            T = CoolProp::Secant(BSR, T_max, -0.01*T_max, 1e-7, 50, errstr);
        }
    }
    else{
        T = CoolProp::Brent(BSR, T_min, T_max, 1e-7, 1e-4, 50, errstr);
    }
    
    return T;
}
static double Secant_Tdb_at_saturated_W(double psi_w, double p, double T_guess)
{
    double T;
    class BrentSolverResids : public CoolProp::FuncWrapper1D
    {
    private:
        double pp_water, psi_w, p, r;
    public:
        BrentSolverResids(double psi_w, double p) : psi_w(psi_w), p(p) { pp_water = psi_w*p; };
        ~BrentSolverResids(){};

        double call(double T){
            double p_ws;
            if (T>=273.16){
                // Saturation pressure [Pa]
                Water->update(CoolProp::QT_INPUTS, 0, T);
                p_ws= Water->keyed_output(CoolProp::iP);
            }
            else{
                // Sublimation pressure [Pa]
                p_ws=psub_Ice(T);
            }
            double f = f_factor(T, p);
            double pp_water_calc = f*p_ws;
            double psi_w_calc = pp_water_calc/p;
            r = (psi_w_calc - psi_w)/psi_w;
            return r;
        }
    };

    BrentSolverResids Resids(psi_w, p);

    std::string errstr; 
    T = CoolProp::Brent(Resids, 150, 350, 1e-16, 1e-7, 100, errstr);

    return T;
}

static double Secant_HAProps_T(const std::string &OutputName, const std::string &Input1Name, double Input1, const std::string &Input2Name, double Input2, double TargetVal, double T_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering T
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=5e-7,f=999,T=300,change;
    int iter=1;
    std::string sT = "T";

    while ((iter<=3 || (std::abs(f)>eps && std::abs(change)>1e-10)) && iter<100)
    {
        if (iter==1){x1=T_guess; T=x1;}
        if (iter==2){x2=T_guess+0.001; T=x2;}
        if (iter>2) {T=x2;}
            f=HAPropsSI(OutputName,sT,T,Input1Name,Input1,Input2Name,Input2)-TargetVal;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change = y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return T;
}

static double Secant_HAProps_W(const std::string &OutputName, const std::string &Input1Name, double Input1, const std::string &Input2Name, double Input2, double TargetVal, double W_guess)
{
    // Use a secant solve in order to yield a target output value for HAProps by altering humidity ratio
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-12,f=999,W=0.0001;
    int iter=1;
    
    if (OutputName == "B"){eps = 1e-7;}

    while ((iter<=3 || std::abs(f)>eps) && iter<100)
    {
        if (iter == 1){x1 = W_guess; W = x1;}
        if (iter == 2){x2 = W_guess*1.1; W = x2;}
        if (iter > 2) {W = x2;}
            f = HAPropsSI(OutputName,"W",W,Input1Name,Input1,Input2Name,Input2)-TargetVal;
        if (iter == 1){y1 = f;}
        if (iter > 1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return W;
}

// Mixed virial components
static double _B_aw(double T)
{
    check_fluid_instantiation();
    // Returns value in m^3/mol
    double a[]={0,0.665687e2,-0.238834e3,-0.176755e3};
    double b[]={0,-0.237,-1.048,-3.183};
    double rhobarstar=1000,Tstar=100;
    return 1/rhobarstar*(a[1]*pow(T/Tstar,b[1])+a[2]*pow(T/Tstar,b[2])+a[3]*pow(T/Tstar,b[3]))/1000; // Correlation has units of dm^3/mol, to convert to m^3/mol, divide by 1000
}

static double _dB_aw_dT(double T)
{
    check_fluid_instantiation();
    // Returns value in m^3/mol
    double a[]={0,0.665687e2,-0.238834e3,-0.176755e3};
    double b[]={0,-0.237,-1.048,-3.183};
    double rhobarstar=1000,Tstar=100;
    return 1/rhobarstar/Tstar*(a[1]*b[1]*pow(T/Tstar,b[1]-1)+a[2]*b[2]*pow(T/Tstar,b[2]-1)+a[3]*b[3]*pow(T/Tstar,b[3]-1))/1000; // Correlation has units of dm^3/mol/K, to convert to m^3/mol/K, divide by 1000
}

static double _C_aaw(double T)
{
    check_fluid_instantiation();
    // Function return has units of m^6/mol^2
    double c[]={0,0.482737e3,0.105678e6,-0.656394e8,0.294442e11,-0.319317e13};
    double rhobarstar=1000,Tstar=1,summer=0; int i;
    for (i=1;i<=5;i++)
    {
        summer+=c[i]*pow(T/Tstar,1-i);
    }
    return 1.0/rhobarstar/rhobarstar*summer/1e6; // Correlation has units of dm^6/mol^2, to convert to m^6/mol^2 divide by 1e6
}

static double _dC_aaw_dT(double T)
{
    check_fluid_instantiation();
    // Function return in units of m^6/mol^2/K
    double c[]={0,0.482737e3,0.105678e6,-0.656394e8,0.294442e11,-0.319317e13};
    double rhobarstar=1000,Tstar=1,summer=0; int i;
    for (i=2;i<=5;i++)
    {
        summer+=c[i]*(1-i)*pow(T/Tstar,-i);
    }
    return 1.0/rhobarstar/rhobarstar/Tstar*summer/1e6; // Correlation has units of dm^6/mol^2/K, to convert to m^6/mol^2/K divide by 1e6
}

static double _C_aww(double T)
{
    check_fluid_instantiation();
    // Function return has units of m^6/mol^2
    double d[]={0,-0.1072887e2,0.347804e4,-0.383383e6,0.334060e8};
    double rhobarstar=1,Tstar=1,summer=0; int i;
    for (i=1;i<=4;i++)
    {
        summer+=d[i]*pow(T/Tstar,1-i);
    }
    return -1.0/rhobarstar/rhobarstar*exp(summer)/1e6; // Correlation has units of dm^6/mol^2, to convert to m^6/mol^2 divide by 1e6
}

static double _dC_aww_dT(double T)
{
    check_fluid_instantiation();
     // Function return in units of m^6/mol^2/K
    double d[]={0,-0.1072887e2,0.347804e4,-0.383383e6,0.334060e8};
    double rhobarstar=1,Tstar=1,summer1=0,summer2=0; int i;
    for (i=1;i<=4;i++)
    {
        summer1+=d[i]*pow(T/Tstar,1-i);
    }
    for (i=2;i<=4;i++)
    {
        summer2+=d[i]*(1-i)*pow(T/Tstar,-i);
    }
    return -1.0/rhobarstar/rhobarstar/Tstar*exp(summer1)*summer2/1e6; // Correlation has units of dm^6/mol^2/K, to convert to m^6/mol^2/K divide by 1e6
}


static double B_m(double T, double psi_w)
{
    // Bm has units of m^3/mol
    double B_aa,B_ww,B_aw;
    if (FlagUseVirialCorrelations==1)
    {
        B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*pow(T,2)
        +4.104150642775e-10*pow(T,3) -1.192780880645e-12*pow(T,4) +2.134201312070e-15*pow(T,5)
        -2.157430412913e-18*pow(T,6) +9.453830907795e-22*pow(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*pow(T,2)
        +1.265864734412e-05*pow(T,3) -4.092175700300e-08*pow(T,4) +7.943925411344e-11*pow(T,5)
        -8.567808759123e-14*pow(T,6) +3.958203548563e-17*pow(T,7);
    }
    else
    {
        B_aa = B_Air(T); // [m^3/mol]
        B_ww = B_Water(T); // [m^3/mol]
    }

    B_aw=_B_aw(T); // [m^3/mol]
    return pow(1-psi_w,2)*B_aa+2*(1-psi_w)*psi_w*B_aw+psi_w*psi_w*B_ww;
}

static double dB_m_dT(double T, double psi_w)
{
    //dBm_dT has units of m^3/mol/K
    double dB_dT_aa,dB_dT_ww,dB_dT_aw;
    if (FlagUseVirialCorrelations)
    {
        dB_dT_aa=1.65159324353e-05 -3.026130954749e-07*T +2.558323847166e-09*pow(T,2) -1.250695660784e-11*pow(T,3) +3.759401946106e-14*pow(T,4) -6.889086380822e-17*pow(T,5) +7.089457032972e-20*pow(T,6) -3.149942145971e-23*pow(T,7);
        dB_dT_ww=0.65615868848 -1.487953162679e-02*T +1.450134660689e-04*pow(T,2) -7.863187630094e-07*pow(T,3) +2.559556607010e-09*pow(T,4) -4.997942221914e-12*pow(T,5) +5.417678681513e-15*pow(T,6) -2.513856275241e-18*pow(T,7);
    }
    else
    {
        dB_dT_aa=dBdT_Air(T); // [m^3/mol]
        dB_dT_ww=dBdT_Water(T); // [m^3/mol]
    }
        dB_dT_aw=_dB_aw_dT(T); // [m^3/mol]
    return pow(1-psi_w,2)*dB_dT_aa+2*(1-psi_w)*psi_w*dB_dT_aw+psi_w*psi_w*dB_dT_ww;
}

static double C_m(double T, double psi_w)
{
    // Cm has units of m^6/mol^2
    double C_aaa,C_www,C_aww,C_aaw;
    if (FlagUseVirialCorrelations)
    {
        C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*pow(T,2)
        -6.234878717893e-15*pow(T,3) +1.791668730770e-17*pow(T,4) -3.175283581294e-20*pow(T,5)
        +3.184306136120e-23*pow(T,6) -1.386043640106e-26*pow(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*pow(T,2)
        +7.687692259692e-07*pow(T,3) -2.571440816920e-09*pow(T,4) +5.147432221082e-12*pow(T,5)
        -5.708156494894e-15*pow(T,6) +2.704605721778e-18*pow(T,7);
    }
    else
    {
        C_aaa=C_Air(T); //[m^6/mol^2]
        C_www=C_Water(T); //[m^6/mol^2]
    }
    C_aaw=_C_aaw(T); //[m^6/mol^2]
    C_aww=_C_aww(T); //[m^6/mol^2]
    return pow(1-psi_w,3)*C_aaa+3*pow(1-psi_w,2)*psi_w*C_aaw+3*(1-psi_w)*psi_w*psi_w*C_aww+pow(psi_w,3)*C_www;
}

static double dC_m_dT(double T, double psi_w)
{
    // dCm_dT has units of m^6/mol^2/K

    double dC_dT_aaa,dC_dT_www,dC_dT_aww,dC_dT_aaw;
    // NDG for fluid EOS for virial terms
    if (FlagUseVirialCorrelations)
    {
        dC_dT_aaa=-2.46582342273e-10 +4.425401935447e-12*T -3.669987371644e-14*pow(T,2) +1.765891183964e-16*pow(T,3) -5.240097805744e-19*pow(T,4) +9.502177003614e-22*pow(T,5) -9.694252610339e-25*pow(T,6) +4.276261986741e-28*pow(T,7);
        dC_dT_www=0.0984601196142 -2.356713397262e-03*T +2.409113323685e-05*pow(T,2) -1.363083778715e-07*pow(T,3) +4.609623799524e-10*pow(T,4) -9.316416405390e-13*pow(T,5) +1.041909136255e-15*pow(T,6) -4.973918480607e-19*pow(T,7);
    }
    else
    {
        dC_dT_aaa=dCdT_Air(T); // [m^6/mol^2]
        dC_dT_www=dCdT_Water(T); // [m^6/mol^2]
    }
    dC_dT_aaw=_dC_aaw_dT(T); // [m^6/mol^2]
    dC_dT_aww=_dC_aww_dT(T); // [m^6/mol^2]
    return pow(1-psi_w,3)*dC_dT_aaa+3*pow(1-psi_w,2)*psi_w*dC_dT_aaw+3*(1-psi_w)*psi_w*psi_w*dC_dT_aww+pow(psi_w,3)*dC_dT_www;
}
double HumidityRatio(double psi_w)
{
    return psi_w*epsilon/(1-psi_w);
}

static double HenryConstant(double T)
{
    // Result has units of 1/Pa
    double p_ws,beta_N2,beta_O2,beta_Ar,beta_a,tau,Tr,Tc=647.096;
    Tr=T/Tc;
    tau=1-Tr;
    Water->update(CoolProp::QT_INPUTS, 1.0, T);
    p_ws = Water->keyed_output(CoolProp::iP); //[Pa]
    beta_N2=p_ws*exp(-9.67578/Tr+4.72162*pow(tau,0.355)/Tr+11.70585*pow(Tr,-0.41)*exp(tau));
    beta_O2=p_ws*exp(-9.44833/Tr+4.43822*pow(tau,0.355)/Tr+11.42005*pow(Tr,-0.41)*exp(tau));
    beta_Ar=p_ws*exp(-8.40954/Tr+4.29587*pow(tau,0.355)/Tr+10.52779*pow(Tr,-0.41)*exp(tau));
    beta_a=1/(0.7812/beta_N2+0.2095/beta_O2+0.0093/beta_Ar);
    return 1/(1.01325*beta_a);
}
double isothermal_compressibility(double T, double p)
{
    double k_T;

    if (T> 273.16)
    {
        if (FlagUseIsothermCompressCorrelation)
        {
            k_T = 1.6261876614E-22*pow(T,6) - 3.3016385196E-19*pow(T,5) + 2.7978984577E-16*pow(T,4)
                - 1.2672392901E-13*pow(T,3) + 3.2382864853E-11*pow(T,2) - 4.4318979503E-09*T + 2.5455947289E-07;
        }
        else
        {
            Water->update(CoolProp::PT_INPUTS, p, T);
            k_T = Water->keyed_output(CoolProp::iisothermal_compressibility);
        }
    }
    else
    {
        k_T = IsothermCompress_Ice(T,p); //[1/Pa]
    }
    return k_T;
}
double f_factor(double T, double p)
{
    double f=0,Rbar=8.314371,eps=1e-8;
    double x1=0,x2=0,x3,y1=0,y2,change=_HUGE;
    int iter=1;
    double p_ws,B_aa,B_aw,B_ww,C_aaa,C_aaw,C_aww,C_www,
        line1,line2,line3,line4,line5,line6,line7,line8,k_T,beta_H,LHS,RHS,psi_ws,
        vbar_ws;

    // Saturation pressure [Pa]
    if (T>273.16)
    {
        // It is liquid water
        Water->update(CoolProp::QT_INPUTS, 0, T);
        p_ws = Water->p();
        vbar_ws = 1.0/Water->keyed_output(CoolProp::iDmolar); //[m^3/mol]
        beta_H = HenryConstant(T); //[1/Pa]
    }
    else
    {
        // It is ice
        p_ws = psub_Ice(T); // [Pa]
        beta_H = 0;
        vbar_ws = dg_dp_Ice(T,p)*MM_Water(); //[m^3/mol]
    }

     k_T = isothermal_compressibility(T,p); //[1/Pa]

    // Hermann: In the iteration process of the enhancement factor in Eq. (3.25), k_T is set to zero for pw,s (T) > p.
    if (p_ws>p)
    {
        k_T=0;
        beta_H=0;
    }

    // NDG for fluid EOS for virial terms
    if (FlagUseVirialCorrelations)
    {
        B_aa=-0.000721183853646 +1.142682674467e-05*T -8.838228412173e-08*pow(T,2)
        +4.104150642775e-10*pow(T,3) -1.192780880645e-12*pow(T,4) +2.134201312070e-15*pow(T,5)
        -2.157430412913e-18*pow(T,6) +9.453830907795e-22*pow(T,7);
        B_ww=-10.8963128394 +2.439761625859e-01*T -2.353884845100e-03*pow(T,2)
        +1.265864734412e-05*pow(T,3) -4.092175700300e-08*pow(T,4) +7.943925411344e-11*pow(T,5)
        -8.567808759123e-14*pow(T,6) +3.958203548563e-17*pow(T,7);
        C_aaa=1.29192158975e-08 -1.776054020409e-10*T +1.359641176409e-12*pow(T,2)
        -6.234878717893e-15*pow(T,3) +1.791668730770e-17*pow(T,4) -3.175283581294e-20*pow(T,5)
        +3.184306136120e-23*pow(T,6) -1.386043640106e-26*pow(T,7);
        C_www=-0.580595811134 +1.365952762696e-02*T -1.375986293288e-04*pow(T,2)
        +7.687692259692e-07*pow(T,3) -2.571440816920e-09*pow(T,4) +5.147432221082e-12*pow(T,5)
        -5.708156494894e-15*pow(T,6) +2.704605721778e-18*pow(T,7);
    }
    else
    {
        B_aa = B_Air(T); // [m^3/mol]
        C_aaa = C_Air(T); // [m^6/mol^2]
        B_ww = B_Water(T); // [m^3/mol]
        C_www = C_Water(T); // [m^6/mol^2]
    }
    B_aw = _B_aw(T); //[m^3/mol]
    C_aaw = _C_aaw(T); //[m^6/mol^2]
    C_aww = _C_aww(T); //[m^6/mol^2]

    // Use a little secant loop to find f iteratively
    // Start out with a guess value of 1 for f
    while ((iter<=3 || change>eps) && iter<100)
    {
        if (iter==1){x1=1.00; f=x1;}
        if (iter==2){x2=1.00+0.000001; f=x2;}
        if (iter>2) {f=x2;}

            // Left-hand-side of Equation 3.25
            LHS=log(f);
            // Eqn 3.24
            psi_ws=f*p_ws/p;

            // All the terms forming the RHS of Eqn 3.25
            line1=((1+k_T*p_ws)*(p-p_ws)-k_T*(p*p-p_ws*p_ws)/2.0)/(Rbar*T)*vbar_ws+log(1-beta_H*(1-psi_ws)*p);
            line2=pow(1-psi_ws,2)*p/(Rbar*T)*B_aa-2*pow(1-psi_ws,2)*p/(Rbar*T)*B_aw-(p-p_ws-pow(1-psi_ws,2)*p)/(Rbar*T)*B_ww;
            line3=pow(1-psi_ws,3)*p*p/pow(Rbar*T,2)*C_aaa+(3*pow(1-psi_ws,2)*(1-2*(1-psi_ws))*p*p)/(2*pow(Rbar*T,2))*C_aaw;
            line4=-3*pow(1-psi_ws,2)*psi_ws*p*p/pow(Rbar*T,2)*C_aww-((3-2*psi_ws)*psi_ws*psi_ws*p*p-p_ws*p_ws)/(2*pow(Rbar*T,2))*C_www;
            line5=-(pow(1-psi_ws,2)*(-2+3*psi_ws)*psi_ws*p*p)/pow(Rbar*T,2)*B_aa*B_ww;
            line6=-(2*pow(1-psi_ws,3)*(-1+3*psi_ws)*p*p)/pow(Rbar*T,2)*B_aa*B_aw;
            line7=(6*pow(1-psi_ws,2)*psi_ws*psi_ws*p*p)/pow(Rbar*T,2)*B_ww*B_aw-(3*pow(1-psi_ws,4)*p*p)/(2*pow(Rbar*T,2))*B_aa*B_aa;
            line8=-(2*pow(1-psi_ws,2)*psi_ws*(-2+3*psi_ws)*p*p)/pow(Rbar*T,2)*B_aw*B_aw-(p_ws*p_ws-(4-3*psi_ws)*pow(psi_ws,3)*p*p)/(2*pow(Rbar*T,2))*B_ww*B_ww;
            RHS=line1+line2+line3+line4+line5+line6+line7+line8;

        if (iter==1){y1=LHS-RHS;}
        if (iter>1)
        {
            y2=LHS-RHS;
            x3=x2-y2/(y2-y1)*(x2-x1);
            change=std::abs(y2/(y2-y1)*(x2-x1));
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    if (f>=1.0)
        return f;
    else
        return 1.0;
}
void HAHelp(void)
{
    printf("Sorry, Need to update!");
}
int returnHumAirCode(const char * Code)
{
    if (!strcmp(Code,"GIVEN_TDP"))
        return GIVEN_TDP;
    else if (!strcmp(Code,"GIVEN_HUMRAT"))
        return GIVEN_HUMRAT;
    else if (!strcmp(Code,"GIVEN_TWB"))
        return GIVEN_TWB;
    else if (!strcmp(Code,"GIVEN_RH"))
        return GIVEN_RH;
    else if (!strcmp(Code,"GIVEN_ENTHALPY"))
        return GIVEN_ENTHALPY;
    else
    {
        fprintf(stderr,"Code to returnHumAirCode in HumAir.c [%s] not understood",Code);
        return -1;
    }
}
double Viscosity(double T, double p, double psi_w)
{
    /*
    Using the method of:

    P.T. Tsilingiris, 2009, Thermophysical and transport properties of humid air at temperature range between 0 and 100 oC, Energy Conversion and Management, 49, 1098-1010

    but using the detailed measurements for pure fluid from IAPWS formulations
    */
    double mu_a,mu_w,Phi_av,Phi_va,Ma,Mw;
    Mw=MM_Water();
    Ma=MM_Air();
    // Viscosity of dry air at dry-bulb temp and total pressure
    Air->update(CoolProp::PT_INPUTS,p,T);
    mu_a=Air->keyed_output(CoolProp::iviscosity);
    // Viscosity of pure saturated water at dry-bulb temperature
    Water->update(CoolProp::PQ_INPUTS,p,1);
    mu_w=Water->keyed_output(CoolProp::iviscosity);
    Phi_av=sqrt(2.0)/4.0*pow(1+Ma/Mw,-0.5)*pow(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2.0)/4.0*pow(1+Mw/Ma,-0.5)*pow(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
    return (1-psi_w)*mu_a/((1-psi_w)+psi_w*Phi_av)+psi_w*mu_w/(psi_w+(1-psi_w)*Phi_va);
}
double Conductivity(double T, double p, double psi_w)
{
    /*
    Using the method of:

    P.T. Tsilingiris, 2009, Thermophysical and transport properties of humid air at temperature range between 0 and 100 oC, Energy Conversion and Management, 49, 1098-1010

    but using the detailed measurements for pure fluid from IAPWS formulations
    */
    double mu_a,mu_w,k_a,k_w,Phi_av,Phi_va,Ma,Mw;
    Mw=MM_Water();
    Ma=MM_Air();

    // Viscosity of dry air at dry-bulb temp and total pressure
    Air->update(CoolProp::PT_INPUTS,p,T);
    mu_a=Air->keyed_output(CoolProp::iviscosity);
    k_a=Air->keyed_output(CoolProp::iconductivity);
    // Viscosity of pure saturated water at dry-bulb temperature
    Water->update(CoolProp::PQ_INPUTS,p,1);
    mu_w=Water->keyed_output(CoolProp::iviscosity);
    k_w=Water->keyed_output(CoolProp::iconductivity);
    Phi_av=sqrt(2.0)/4.0*pow(1+Ma/Mw,-0.5)*pow(1+sqrt(mu_a/mu_w)*pow(Mw/Ma,0.25),2); //[-]
    Phi_va=sqrt(2.0)/4.0*pow(1+Mw/Ma,-0.5)*pow(1+sqrt(mu_w/mu_a)*pow(Ma/Mw,0.25),2); //[-]
    return (1-psi_w)*k_a/((1-psi_w)+psi_w*Phi_av)+psi_w*k_w/(psi_w+(1-psi_w)*Phi_va);
}
double MolarVolume(double T, double p, double psi_w)
{
    // Output in m^3/mol
    int iter;
    double v_bar0, v_bar=0, R_bar=8.314472,x1=0,x2=0,x3,y1=0,y2,resid,eps,Bm,Cm;

    // -----------------------------
    // Iteratively find molar volume
    // -----------------------------

    // Start by assuming it is an ideal gas to get initial guess
    v_bar0=R_bar*T/p;

    //Bring outside the loop since not a function of v_bar
    Bm=B_m(T,psi_w);
    Cm=C_m(T,psi_w);

    iter=1; eps=1e-11; resid=999;
    while ((iter<=3 || std::abs(resid)>eps) && iter<100)
    {
        if (iter==1){x1=v_bar0; v_bar=x1;}
        if (iter==2){x2=v_bar0+0.000001; v_bar=x2;}
        if (iter>2) {v_bar=x2;}

            // want v_bar in m^3/mol and R_bar in J/mol-K
            resid = (p-(R_bar)*T/v_bar*(1+Bm/v_bar+Cm/(v_bar*v_bar)))/p;

        if (iter==1){y1=resid;}
        if (iter>1)
        {
            y2=resid;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return v_bar;
}
double IdealGasMolarEnthalpy_Water(double T, double vmolar)
{
    double hbar_w_0, tau, rhomolar, hbar_w;
    // Ideal-Gas contribution to enthalpy of water
    hbar_w_0 = -0.01102303806; //[J/mol]
    
    // Calculate the offset in the water enthalpy from a given state with a known (desired) enthalpy
    double Tref = 473.15, vmolarref = 0.038837428192186184, href = 51885.582451893446;
    Water->update(CoolProp::DmolarT_INPUTS,1/vmolarref,Tref);
    double tauref = Water->keyed_output(CoolProp::iT_reducing)/Tref; //[no units]
    double href_EOS = R_bar*Tref*(1+tauref*Water->keyed_output(CoolProp::idalpha0_dtau_constdelta));
    double hoffset = href - href_EOS;
    
    tau = Water->keyed_output(CoolProp::iT_reducing)/T;
    rhomolar = 1/vmolar; //[mol/m^3]
    Water->specify_phase(CoolProp::iphase_gas);
    Water->update_DmolarT_direct(rhomolar, T);
    Water->unspecify_phase();
    hbar_w = hbar_w_0 + hoffset + R_bar*T*(1+tau*Water->keyed_output(CoolProp::idalpha0_dtau_constdelta));
    return hbar_w;
}
double IdealGasMolarEntropy_Water(double T, double p)
{
    double sbar_w, tau, R_bar;
    R_bar = 8.314371; //[J/mol/K]
    
    // Calculate the offset in the water entropy from a given state with a known (desired) entropy
    double Tref = 473.15, pref = 101325, sref = 141.18297895840303;
    Water->update(CoolProp::DmolarT_INPUTS,pref/(R_bar*Tref),Tref);
    double tauref = Water->keyed_output(CoolProp::iT_reducing)/Tref; //[no units]
    double sref_EOS = R_bar*(tauref*Water->keyed_output(CoolProp::idalpha0_dtau_constdelta)-Water->keyed_output(CoolProp::ialpha0));
    double soffset = sref - sref_EOS;
    
    tau = Water->keyed_output(CoolProp::iT_reducing)/T;
    Water->specify_phase(CoolProp::iphase_gas);
    Water->update(CoolProp::DmolarT_INPUTS,p/(R_bar*T),T);
    Water->unspecify_phase();
    sbar_w = soffset + R_bar*(tau*Water->keyed_output(CoolProp::idalpha0_dtau_constdelta)-Water->keyed_output(CoolProp::ialpha0)); //[kJ/kmol/K]
    return sbar_w;
}
double IdealGasMolarEnthalpy_Air(double T, double vmolar)
{
    double hbar_a_0, tau, rhomolar, hbar_a, R_bar_Lemmon;
    // Ideal-Gas contribution to enthalpy of air
    hbar_a_0 = -7914.149298; //[J/mol]
    
    R_bar_Lemmon = 8.314510; //[J/mol/K]
    // Calculate the offset in the air enthalpy from a given state with a known (desired) enthalpy
    double Tref = 473.15, vmolarref = 0.038837428192186184, href = 13782.240592933371;
    Air->update(CoolProp::DmolarT_INPUTS, 1/vmolarref, Tref);
    double tauref = 132.6312/Tref; //[no units]
    double href_EOS = R_bar_Lemmon*Tref*(1+tauref*Air->keyed_output(CoolProp::idalpha0_dtau_constdelta)); 
    double hoffset = href - href_EOS;
    
    // Tj is given by 132.6312 K
    tau = 132.6312/T;
    rhomolar = 1/vmolar; //[mol/m^3]
    // Now calculate it based on the given inputs
    Air->specify_phase(CoolProp::iphase_gas);
    Air->update_DmolarT_direct(rhomolar, T);
    Air->unspecify_phase();
    hbar_a = hbar_a_0 + hoffset + R_bar_Lemmon*T*(1+tau*Air->keyed_output(CoolProp::idalpha0_dtau_constdelta)); //[J/mol]
    return hbar_a;
}
double IdealGasMolarEntropy_Air(double T, double vmolar_a)
{
    double sbar_0_Lem, tau, sbar_a, R_bar_Lemmon = 8.314510, T0=273.15, p0=101325, vmolar_a_0;

    // Ideal-Gas contribution to entropy of air
    sbar_0_Lem = -196.1375815; //[J/mol/K]

    vmolar_a_0 = R_bar_Lemmon*T0/p0; //[m^3/mol]

    // Calculate the offset in the air entropy from a given state with a known (desired) entropy
    double Tref = 473.15, vmolarref = 0.038837605637863169, sref = 212.22365283759311;
    Air->update(CoolProp::DmolarT_INPUTS, 1/vmolar_a_0, Tref);
    double tauref = 132.6312/Tref; //[no units]
    double sref_EOS = R_bar_Lemmon*(tauref*Air->keyed_output(CoolProp::idalpha0_dtau_constdelta)-Air->keyed_output(CoolProp::ialpha0))+R_bar_Lemmon*log(vmolarref/vmolar_a_0);
    double soffset = sref - sref_EOS;
    
    // Tj and rhoj are given by 132.6312 and 302.5507652 respectively
    tau = 132.6312/T; //[no units]
    
    Air->specify_phase(CoolProp::iphase_gas);
    Air->update_DmolarT_direct(1/vmolar_a_0,T);
    Air->unspecify_phase();
    sbar_a=sbar_0_Lem + soffset + R_bar_Lemmon*(tau*Air->keyed_output(CoolProp::idalpha0_dtau_constdelta)-Air->keyed_output(CoolProp::ialpha0))+R_bar_Lemmon*log(vmolar_a/vmolar_a_0); //[J/mol/K]

    return sbar_a; //[J/mol/K]
}

double MolarEnthalpy(double T, double p, double psi_w, double vmolar)
{
    // In units of kJ/kmol

    // vbar (molar volume) in m^3/kg

    double hbar_0, hbar_a, hbar_w, hbar, R_bar=8.314472;
    // ----------------------------------------
    //      Enthalpy
    // ----------------------------------------
    // Constant for enthalpy
    // Not clear why getting rid of this term yields the correct values in the table, but enthalpies are equal to an additive constant, so not a big deal
    hbar_0=0.0;//2.924425468; //[kJ/kmol]

    if (FlagUseIdealGasEnthalpyCorrelations){
        hbar_w = 2.7030251618E-03*T*T + 3.1994361015E+01*T + 3.6123174929E+04;
        hbar_a = 9.2486716590E-04*T*T + 2.8557221776E+01*T - 7.8616129429E+03;
    }
    else{
        hbar_w = IdealGasMolarEnthalpy_Water(T, vmolar);
        hbar_a = IdealGasMolarEnthalpy_Air(T, vmolar);
    }

    // If the user changes the reference state for water or Air, we need to ensure that the values returned from this 
    // function are always the same as the formulation expects.  Therefore we can use a state point for which we know what the 
    // enthalpy should be and then correct the calculated values for the enthalpy.

    hbar = hbar_0+(1-psi_w)*hbar_a+psi_w*hbar_w+R_bar*T*((B_m(T,psi_w)-T*dB_m_dT(T,psi_w))/vmolar+(C_m(T,psi_w)-T/2.0*dC_m_dT(T,psi_w))/(vmolar*vmolar));
    return hbar; //[J/mol]
}
double MassEnthalpy(double T, double p, double psi_w)
{
    double vmolar = MolarVolume(T, p, psi_w); //[m^3/mol_ha]
    double h_bar = MolarEnthalpy(T, p, psi_w, vmolar); //[J/mol_ha]
    double W = HumidityRatio(psi_w); //[kg_w/kg_da]
    double M_ha = MM_Water()*psi_w+(1-psi_w)*0.028966; // [kg_ha/mol_ha]
    // (1+W) is kg_ha/kg_da
    return h_bar*(1+W)/M_ha; //[J/kg_da]
}

double MolarEntropy(double T, double p, double psi_w, double v_bar)
{
    // In units of J/mol/K

    // Serious typo in RP-1485 - should use total pressure rather than
    // reference pressure in density calculation for water vapor molar entropy

    // vbar (molar volume) in m^3/mol
    double x1=0,x2=0,x3=0,y1=0,y2=0,eps=1e-8,f=999,R_bar_Lem=8.314510;
    int iter=1;
    double sbar_0,sbar_a=0,sbar_w=0,sbar,R_bar=8.314472,vbar_a_guess, Baa, Caaa,vbar_a=0;
    double B,dBdT,C,dCdT;
    // Constant for entropy
    sbar_0=0.02366427495;  //[J/mol/K]

    // Calculate vbar_a, the molar volume of dry air
    // B_m, C_m, etc. functions take care of the units
    Baa = B_m(T,0);
    B = B_m(T,psi_w);
    dBdT = dB_m_dT(T,psi_w);
    Caaa = C_m(T,0);
    C = C_m(T,psi_w);
    dCdT = dC_m_dT(T,psi_w);

    vbar_a_guess = R_bar_Lem*T/p; //[m^3/mol] since p in [Pa]

    while ((iter<=3 || std::abs(f)>eps) && iter<100)
    {
        if (iter==1){x1=vbar_a_guess; vbar_a=x1;}
        if (iter==2){x2=vbar_a_guess+0.001; vbar_a=x2;}
        if (iter>2) {vbar_a=x2;}
            f=R_bar_Lem*T/vbar_a*(1+Baa/vbar_a+Caaa/pow(vbar_a,2))-p;
        if (iter==1){y1=f;}
        if (iter>1)
        {
            y2=f;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
        if (iter>100){ return _HUGE; }
    }

    if (FlagUseIdealGasEnthalpyCorrelations){
        std::cout << "Not implemented" << std::endl;
    }
    else{
        sbar_w=IdealGasMolarEntropy_Water(T,p);
        sbar_a=IdealGasMolarEntropy_Air(T,vbar_a);
    }
    if (psi_w!=0){
        sbar = sbar_0+(1-psi_w)*sbar_a+psi_w*sbar_w-R_bar*( (B+T*dBdT)/v_bar+(C+T*dCdT)/(2*pow(v_bar,2))+(1-psi_w)*log(1-psi_w)+psi_w*log(psi_w));
    }
    else{
        sbar = sbar_0+sbar_a;
    }
    return sbar; //[kJ/kmol/K]
}

double DewpointTemperature(double T, double p, double psi_w)
{
    int iter;
    double p_w,eps,resid,Tdp=0,x1=0,x2=0,x3,y1=0,y2,T0;
    double p_ws_dp,f_dp;

    // Make sure it isn't dry air, return an impossible temperature otherwise
    if ((1-psi_w)<1e-16)
    {
        return -1;
    }
    // ------------------------------------------
    // Iteratively find the dewpoint temperature
    // ------------------------------------------

    // The highest dewpoint temperature possible is the dry-bulb temperature.
    // When they are equal, the air is saturated (R=1)

    p_w = psi_w*p;

    // 611.65... is the triple point pressure of water in Pa
    if (p_w > 611.6547241637944){
        Water->update(CoolProp::PQ_INPUTS, p, 1.0);
        T0 = Water->T()-1;
    }
    else{
        T0 = 268;
    }
    // A good guess for Tdp is that enhancement factor is unity, which yields
    // p_w_s = p_w, and get guess for T from saturation temperature

    iter=1; eps=1e-5; resid=999;
    while ((iter<=3 || std::abs(resid)>eps) && iter<100)
    {
        if (iter==1){x1 = T0; Tdp=x1;}
        if (iter==2){x2 = x1 + 0.1; Tdp=x2;}
        if (iter>2) {Tdp=x2;}

            if (Tdp >= 273.16)
            {
                // Saturation pressure at dewpoint [Pa]
                Water->update(CoolProp::QT_INPUTS, 0.0, Tdp);
                p_ws_dp = Water->keyed_output(CoolProp::iP);
            }
            else
            {
                // Sublimation pressure at icepoint [Pa]
                p_ws_dp=psub_Ice(Tdp);
            }
            // Enhancement Factor at dewpoint temperature [-]
            f_dp=f_factor(Tdp,p);
            // Error between target and actual pressure [Pa]
            resid=p_w-p_ws_dp*f_dp;

        if (iter==1){y1=resid;}
        if (iter>1)
        {
            y2=resid;
            x3=x2-y2/(y2-y1)*(x2-x1);
            y1=y2; x1=x2; x2=x3;
        }
        iter=iter+1;
    }
    return Tdp;
}

class WetBulbSolver : public CoolProp::FuncWrapper1D
{
private:
    double _T,_p,_W,LHS,RHS,v_bar_w,M_ha;
public:
    WetBulbSolver(double T, double p, double psi_w){
        _T = T;
        _p = p;
        _W = epsilon*psi_w/(1-psi_w);

        //These things are all not a function of Twb
        v_bar_w = MolarVolume(T,p,psi_w);
        M_ha = MM_Water()*psi_w+(1-psi_w)*0.028966;
        LHS = MolarEnthalpy(T,p,psi_w,v_bar_w)*(1+_W)/M_ha;
    };
    ~WetBulbSolver(){};
    double call(double Twb)
    {
        double epsilon=0.621945;
        double f_wb,p_ws_wb,p_s_wb,W_s_wb,h_w,M_ha_wb,psi_wb,v_bar_wb;

        // Enhancement Factor at wetbulb temperature [-]
        f_wb=f_factor(Twb,_p);
        if (Twb > 273.16)
        {
            // Saturation pressure at wetbulb temperature [Pa]
            Water->update(CoolProp::QT_INPUTS,0,Twb);
            p_ws_wb= Water->keyed_output(CoolProp::iP);
        }
        else
        {
            // Sublimation pressure at wetbulb temperature [kPa]
            p_ws_wb=psub_Ice(Twb);
        }

        // Vapor pressure
        p_s_wb = f_wb*p_ws_wb;
        // wetbulb humidity ratio
        W_s_wb = epsilon*p_s_wb/(_p-p_s_wb);
        // wetbulb water mole fraction
        psi_wb = W_s_wb/(epsilon+W_s_wb);
        if (Twb > 273.16)
        {
            // Enthalpy of water [J/kg_water]
            Water->update(CoolProp::PT_INPUTS, _p, Twb);
            h_w = Water->keyed_output(CoolProp::iHmass); //[J/kg_water]
        }
        else
        {
            // Enthalpy of ice [J/kg_water]
            h_w=h_Ice(Twb,_p);
        }
        // Mole masses of wetbulb and humid air

        M_ha_wb = MM_Water()*psi_wb+(1-psi_wb)*0.028966;
        v_bar_wb=MolarVolume(Twb,_p,psi_wb);
        RHS = (MolarEnthalpy(Twb,_p,psi_wb,v_bar_wb)*(1+W_s_wb)/M_ha_wb+(_W-W_s_wb)*h_w);
        if (!ValidNumber(LHS-RHS)){throw CoolProp::ValueError();}
        return LHS - RHS;
    }
};

class WetBulbTminSolver : public CoolProp::FuncWrapper1D
{
public:
    double p,hair_dry,r, RHS;
    WetBulbTminSolver(double p, double hair_dry){
        this->p = p;
        this->hair_dry = hair_dry;
    };
    ~WetBulbTminSolver(){};
    double call(double Ts)
    {
        RHS = HAPropsSI("H","T",Ts,"P",p,"R",1);
        if (!ValidNumber(RHS)){throw CoolProp::ValueError();}
        r = RHS - this->hair_dry;
        return r;
    }
};

double WetbulbTemperature(double T, double p, double psi_w)
{
    // ------------------------------------------
    // Iteratively find the wetbulb temperature
    // ------------------------------------------
    //
    // If the temperature is less than the saturation temperature of water
    // for the given atmospheric pressure, the highest wetbulb temperature that is possible is the dry bulb
    // temperature
    //
    // If the temperature is above the saturation temperature corresponding to the atmospheric pressure,
    // then the maximum value for the wetbulb temperature is the saturation temperature
    double Tmax = T;
    Water->update(CoolProp::PQ_INPUTS,p,1.0);
    double Tsat = Water->keyed_output(CoolProp::iT);
    if (T >= Tsat)
    {
        Tmax = Tsat;
    }

    // Instantiate the solver container class
    WetBulbSolver WBS(T, p, psi_w);

    std::string errstr;

    double return_val;
    try{
        return_val = Brent(WBS,Tmax+1,200, DBL_EPSILON, 1e-12, 50, errstr);

        // Solution obtained is out of range (T>Tmax)
        if (return_val > Tmax + 1) {throw CoolProp::ValueError();}
    }
    catch(std::exception &)
    {
        // The lowest wetbulb temperature that is possible for a given dry bulb temperature
        // is the saturated air temperature which yields the enthalpy of dry air at dry bulb temperature

        try{
            double hair_dry = MassEnthalpy(T,p,0);

            // Directly solve for the saturated temperature that yields the enthalpy desired
            WetBulbTminSolver WBTS(p,hair_dry);
            double Tmin = Brent(WBTS,210,Tsat-1,1e-12,1e-12,50,errstr);

            return_val = Brent(WBS,Tmin-30,Tmax-1,1e-12,1e-12,50,errstr);
        }
        catch(std::exception)
        {
            return_val = _HUGE;
        }
    }
    return return_val;
}
static givens Name2Type(const std::string &Name)
{
    if (!strcmp(Name,"Omega") || !strcmp(Name,"HumRat") || !strcmp(Name,"W"))
        return GIVEN_HUMRAT;
    else if (!strcmp(Name,"psi_w") || !strcmp(Name, "Y"))
        return GIVEN_PSIW;
    else if (!strcmp(Name,"Tdp") || !strcmp(Name,"T_dp") || !strcmp(Name,"DewPoint") || !strcmp(Name,"D"))
        return GIVEN_TDP;
    else if (!strcmp(Name,"Twb") || !strcmp(Name,"T_wb") || !strcmp(Name,"WetBulb") || !strcmp(Name,"B"))
        return GIVEN_TWB;
    else if (!strcmp(Name,"Enthalpy") || !strcmp(Name,"H") || !strcmp(Name,"Hda"))
        return GIVEN_ENTHALPY;
    else if (!strcmp(Name,"Hha"))
        return GIVEN_ENTHALPY_HA;
    else if (!strcmp(Name,"Entropy") || !strcmp(Name,"S"))
        return GIVEN_ENTROPY;        
    else if (!strcmp(Name,"RH") || !strcmp(Name,"RelHum") || !strcmp(Name,"R"))
        return GIVEN_RH;
    else if (!strcmp(Name,"Tdb") || !strcmp(Name,"T_db") || !strcmp(Name,"T"))
        return GIVEN_T;
    else if (!strcmp(Name,"P"))
        return GIVEN_P;
    else if (!strcmp(Name,"V") || !strcmp(Name,"Vda"))
        return GIVEN_VDA;
    else if (!strcmp(Name,"Vha"))
        return GIVEN_VHA;
    else if (!strcmp(Name,"mu") || !strcmp(Name,"Visc") || !strcmp(Name,"M"))
        return GIVEN_VISC;
    else if (!strcmp(Name,"k") || !strcmp(Name,"Conductivity") || !strcmp(Name,"K"))
        return GIVEN_COND;
    else if (!strcmp(Name,"C") || !strcmp(Name,"cp"))
        return GIVEN_CP;
    else if (!strcmp(Name,"Cha") || !strcmp(Name,"cp_ha"))
        return GIVEN_CPHA;
    else
        throw CoolProp::ValueError(format("Sorry, your input [%s] was not understood to Name2Type. Acceptable values are T,P,R,W,D,B,H,S,M,K and aliases thereof\n",Name.c_str()));
}
int TypeMatch(int TypeCode, const std::string &Input1Name, const std::string &Input2Name, const std::string &Input3Name)
{
    // Return the index of the input variable that matches the input, otherwise return -1 for failure
    if (TypeCode==Name2Type(Input1Name))
        return 1;
    if (TypeCode==Name2Type(Input2Name))
        return 2;
    if (TypeCode==Name2Type(Input3Name))
        return 3;
    else
        return -1;
}
double MoleFractionWater(double T, double p, int HumInput, double InVal)
{
    double p_ws,f,W,epsilon=0.621945,Tdp,p_ws_dp,f_dp,p_w_dp,p_s,RH;

    if (HumInput==GIVEN_HUMRAT) //(2)
    {
        W=InVal;
        return W/(epsilon+W);
    }
    else if (HumInput==GIVEN_RH)
    {
        if (T>=273.16)
        {
            // Saturation pressure [Pa]
            Water->update(CoolProp::QT_INPUTS,0,T);
            p_ws= Water->keyed_output(CoolProp::iP);;
        }
        else
        {
            // Sublimation pressure [Pa]
            p_ws=psub_Ice(T);
        }
        // Enhancement Factor [-]
        f=f_factor(T,p);

        // Saturation pressure [Pa]
        p_s=f*p_ws;
        RH=InVal;
        W=epsilon*RH*p_s/(p-RH*p_s);
        return W/(epsilon+W);
    }
    else if (HumInput==GIVEN_TDP)
    {
        Tdp=InVal;
        // Saturation pressure at dewpoint [Pa]
        if (Tdp>=273.16)
        {
            Water->update(CoolProp::QT_INPUTS,0,Tdp);
            p_ws_dp = Water->keyed_output(CoolProp::iP); //[Pa]
        }
        else{
            // Sublimation pressure [Pa]
            p_ws_dp=psub_Ice(Tdp);
        }

        // Enhancement Factor at dewpoint temperature [-]
        f_dp=f_factor(Tdp,p);
        // Water vapor pressure at dewpoint [Pa]
        p_w_dp=f_dp*p_ws_dp;
        // Water mole fraction [-]
        return p_w_dp/p;
    }
    else
    {
        return -1000000;
    }
}

double RelativeHumidity(double T, double p, double psi_w)
{
    double p_ws, f, p_s, W;
    if (T >= 273.16){
        // Saturation pressure [Pa]
        Water->update(CoolProp::QT_INPUTS, 0, T);
        p_ws = Water->keyed_output(CoolProp::iP); //[Pa]
    }
    else{
        // sublimation pressure [Pa]
        p_ws = psub_Ice(T);
    }
    // Enhancement Factor [-]
    f = f_factor(T,p);

    // Saturation pressure [Pa]
    p_s = f*p_ws;
    // Find humidity ratio
    W = HumidityRatio(psi_w);
    // Find relative humidity using W/e=phi*p_s/(p-phi*p_s)
    return W/epsilon*p/(p_s*(1+W/epsilon));
}

void convert_to_SI(const std::string &Name, double &val)
{
    switch(Name2Type(Name))
    {
        case GIVEN_COND:
        case GIVEN_ENTHALPY:
        case GIVEN_ENTHALPY_HA:
        case GIVEN_ENTROPY:
        case GIVEN_ENTROPY_HA:
        case GIVEN_CP:
        case GIVEN_CPHA:
        case GIVEN_P:
            val *= 1000; return;
        case GIVEN_T:
        case GIVEN_TDP:
        case GIVEN_TWB:
        case GIVEN_RH:
        case GIVEN_VDA:
        case GIVEN_VHA:
        case GIVEN_HUMRAT:
        case GIVEN_VISC:
        case GIVEN_PSIW:
            return;
        case GIVEN_INVALID:
            throw CoolProp::ValueError(format("invalid input to convert_to_SI"));
    }
}
void convert_from_SI(const std::string &Name, double &val)
{
    switch(Name2Type(Name))
    {
        case GIVEN_COND:
        case GIVEN_ENTHALPY:
        case GIVEN_ENTHALPY_HA:
        case GIVEN_ENTROPY:
        case GIVEN_ENTROPY_HA:
        case GIVEN_CP:
        case GIVEN_CPHA:
        case GIVEN_P:
            val /= 1000; return;
        case GIVEN_T:
        case GIVEN_TDP:
        case GIVEN_TWB:
        case GIVEN_RH:
        case GIVEN_VDA:
        case GIVEN_VHA:
        case GIVEN_HUMRAT:
        case GIVEN_VISC:
        case GIVEN_PSIW:
            return;
        case GIVEN_INVALID:
            throw CoolProp::ValueError(format("invalid input to convert_to_SI"));
    }
}
double HAProps(const std::string &OutputName, const std::string &Input1Name, double Input1, const std::string &Input2Name, double Input2, const std::string &Input3Name, double Input3)
{
    convert_to_SI(Input1Name, Input1);
    convert_to_SI(Input2Name, Input2);
    convert_to_SI(Input3Name, Input3);
    
    double out = HAPropsSI(OutputName, Input1Name, Input1, Input2Name, Input2, Input3Name, Input3);
    
    convert_from_SI(OutputName, out);
    
    return out;
}

double HAPropsSI(const std::string &OutputName, const std::string &Input1Name, double Input1, const std::string &Input2Name, double Input2, const std::string &Input3Name, double Input3)
{
    try
    {
        // Add a check to make sure that Air and Water fluid states have been properly instantiated
        check_fluid_instantiation();
        
        int iT, iW, iTdp, iRH, ip;
        givens In1Type, In2Type, In3Type, Type1, Type2, OutputType;
        double vals[3],p,T,RH,W,Tdp,psi_w,M_ha,v_bar,h_bar,s_bar,MainInputValue,SecondaryInputValue;
        double Value1,Value2,W_guess;
        std::string MainInputName, SecondaryInputName, Name1, Name2;

        vals[0]=Input1;
        vals[1]=Input2;
        vals[2]=Input3;
        
        OutputType=Name2Type(OutputName.c_str());

        // First figure out what kind of inputs you have, convert names to Macro expansions
        In1Type=Name2Type(Input1Name.c_str());
        In2Type=Name2Type(Input2Name.c_str());
        In3Type=Name2Type(Input3Name.c_str());
        
        if (OutputType == In1Type){return Input1;}
        if (OutputType == In2Type){return Input2;}
        if (OutputType == In3Type){return Input3;}

        // Pressure must be included
        ip=TypeMatch(GIVEN_P,Input1Name,Input2Name,Input3Name);
        if (ip>0)
            p=vals[ip-1];
        else
            return -1000;

        // -----------------------------------------------------------------------------------------------------
        // Check whether the remaining values give explicit solution for mole fraction of water - nice and fast
        // -----------------------------------------------------------------------------------------------------

        // Find the codes if they are there
        iT=   TypeMatch(GIVEN_T,Input1Name,Input2Name,Input3Name);
        iRH=  TypeMatch(GIVEN_RH,Input1Name,Input2Name,Input3Name);
        iW=   TypeMatch(GIVEN_HUMRAT,Input1Name,Input2Name,Input3Name);
        iTdp= TypeMatch(GIVEN_TDP,Input1Name,Input2Name,Input3Name);

        if (iT>0) // Found T (or alias) as an input
        {
            T=vals[iT-1];
            if (iRH>0) //Relative Humidity is provided
            {
                RH=vals[iRH-1];
                psi_w=MoleFractionWater(T,p,GIVEN_RH,RH);
            }
            else if (iW>0)
            {
                W=vals[iW-1];
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
            }
            else if (iTdp>0)
            {
                Tdp=vals[iTdp-1];
                psi_w=MoleFractionWater(T,p,GIVEN_TDP,Tdp);
            }
            else
            {
                // Temperature and pressure are known, figure out which variable holds the other value
                if (In1Type!=GIVEN_T && In1Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName, Input1Name);
                    SecondaryInputValue=Input1;
                }
                else if (In2Type!=GIVEN_T && In2Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName, Input2Name);
                    SecondaryInputValue=Input2;
                }
                else if (In3Type!=GIVEN_T && In3Type!=GIVEN_P)
                {
                    strcpy(SecondaryInputName, Input3Name);
                    SecondaryInputValue=Input3;
                }
                else{
                    return _HUGE;
                }
                // Find the value for W
                W_guess=0.0001;
                W=Secant_HAProps_W(SecondaryInputName,"P",p,"T",T,SecondaryInputValue,W_guess);
                // Mole fraction of water
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
                // And on to output...
            }
        }
        else
        {
            // Need to iterate to find dry bulb temperature since temperature is not provided

            // Pick one input, and alter T to match the other input

            // Get the variables and their values that are NOT pressure for simplicity
            // because you know you need pressure as an input and you already have
            // its value in variable p
            if (ip==1) // Pressure is in slot 1
            {
                strcpy(Name1,Input2Name);
                Value1=Input2;
                strcpy(Name2,Input3Name);
                Value2=Input3;
            }
            else if (ip==2) // Pressure is in slot 2
            {
                strcpy(Name1,Input1Name);
                Value1=Input1;
                strcpy(Name2,Input3Name);
                Value2=Input3;
            }
            else if (ip==3) // Pressure is in slot 3
            {
                strcpy(Name1,Input1Name);
                Value1=Input1;
                strcpy(Name2,Input2Name);
                Value2=Input2;
            }
            else{
            return _HUGE;
            }

            // Get the integer type codes
            Type1=Name2Type(Name1);
            Type2=Name2Type(Name2);
            
            // First see if humidity ratio is provided, will be the fastest option
            if (Type1 == GIVEN_HUMRAT)
            {
                // MainInput is the one that you are using in the call to HAProps
                MainInputValue=Value1; strcpy(MainInputName,Name1);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue=Value2; strcpy(SecondaryInputName,Name2);
            }
            else if (Type2 == GIVEN_HUMRAT)
            {
                // MainInput is the one that you are using in the call to HAPropsSI
                MainInputValue=Value2; strcpy(MainInputName, Name2);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue=Value1; strcpy(SecondaryInputName, Name1);
            }
            // Next, if one of the inputs is something that can potentially yield
            // an explicit solution at a given iteration of the solver, use it
            else if (Type1 == GIVEN_RH || Type1 == GIVEN_TDP)
            {
                // MainInput is the one that you are using in the call to HAProps
                MainInputValue = Value1; strcpy(MainInputName, Name1);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue = Value2; strcpy(SecondaryInputName, Name2);
            }
            else if (Type2 == GIVEN_RH || Type2 == GIVEN_TDP)
            {
                // MainInput is the one that you are using in the call to HAProps
                MainInputValue = Value2; strcpy(MainInputName, Name2);
                // SecondaryInput is the one that you are trying to match
                SecondaryInputValue=Value1; strcpy(SecondaryInputName, Name1);
            }
            else
            {
                CoolProp::ValueError("Sorry, but currently at least one of the variables as an input to HAPropsSI() must be temperature, relative humidity, humidity ratio, or dewpoint\n  Eventually will add a 2-D NR solver to find T and psi_w simultaneously, but not included now\n");
                return _HUGE;
            }

            double T_min = 200;
            double T_max = 450;
			
			if (Name2Type(MainInputName) == GIVEN_RH){
                if (MainInputValue < 1e-10){
                    T_max = 1000;
                }
				else{
                    T_max = CoolProp::PropsSI("T","P",p,"Q",0,"Water") - 1;
                }
			}
            // Minimum drybulb temperature is the drybulb temperature corresponding to saturated air for the humidity ratio
            // if the humidity ratio is provided
            else if (Name2Type(MainInputName) == GIVEN_HUMRAT){
                if (MainInputValue < 1e-10){
                    T_min = 135; // Around the critical point of dry air
                    T_max = 1000;
                }
                else{
                    // Calculate the saturated humid air water partial pressure;
                    double psi_w_sat = MoleFractionWater(T_min, p, GIVEN_HUMRAT, MainInputValue);
                    //double pp_water_sat = psi_w_sat*p; // partial pressure of water, which is equal to f*p_{w_s}
                    // Iteratively solve for temperature that will give desired pp_water_sat
                    T_min = Secant_Tdb_at_saturated_W(psi_w_sat, p, T_min);
                }
            }

			try{
                // Use the Brent's method solver to find T.  Slow but reliable
                T = Brent_HAProps_T(SecondaryInputName, "P", p, MainInputName, MainInputValue, SecondaryInputValue, T_min,T_max);
			}
			catch(std::exception &e){
				CoolProp::set_error_string(e.what());
				return _HUGE;
			}

            // If you want the temperature, return it
            if (Name2Type(OutputName)==GIVEN_T)
                return T;
            else
            {
                // Otherwise, find psi_w for further calculations in the following section
                W=HAPropsSI((char *)"W",(char *)"T",T,(char *)"P",p,MainInputName,MainInputValue);
                psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
            }
        }
        M_ha=(1-psi_w)*0.028966+MM_Water()*psi_w; //[kg_ha/mol_ha]

        // -----------------------------------------------------------------
        // Calculate and return the desired value for known set of T,p,psi_w
        // -----------------------------------------------------------------
        switch (OutputType){
            case GIVEN_VDA:
            {
                v_bar = MolarVolume(T,p,psi_w); //[m^3/mol_ha]
                W = HumidityRatio(psi_w); //[kg_w/kg_a]
                return v_bar*(1+W)/M_ha; //[m^3/kg_da]
            }
            case GIVEN_VHA:
            {
                v_bar = MolarVolume(T,p,psi_w); //[m^3/mol_ha]
                return v_bar/M_ha; //[m^3/kg_ha]
            }
            case GIVEN_PSIW:
            {
                return psi_w; //[mol_w/mol]
            }
            case GIVEN_ENTHALPY:
            {
                return MassEnthalpy(T,p,psi_w); //[J/kg_da]
            }
            case GIVEN_ENTHALPY_HA:
            {
                v_bar = MolarVolume(T, p, psi_w); //[m^3/mol_ha]
                h_bar = MolarEnthalpy(T, p, psi_w,v_bar); //[J/mol_ha]
                return h_bar/M_ha; //[J/kg_ha]
            }
            case GIVEN_ENTROPY_HA:
            {
                v_bar = MolarVolume(T, p, psi_w); //[m^3/mol_ha]
                s_bar = MolarEntropy(T, p, psi_w, v_bar); //[J/mol_da/K]
                W = HumidityRatio(psi_w); //[kg_w/kg_da]
                return s_bar*(1+W)/M_ha; //[kJ/kg_ha/K]
            }
            case GIVEN_TDP:{
                return DewpointTemperature(T,p,psi_w); //[K]
            }
            case GIVEN_TWB:{
                return WetbulbTemperature(T,p,psi_w); //[K]
            }
            case GIVEN_HUMRAT:{
                return HumidityRatio(psi_w);
            }
            case GIVEN_RH:{
                return RelativeHumidity(T,p,psi_w);
            }
            case GIVEN_VISC:{
                return Viscosity(T,p,psi_w);
            }
            case GIVEN_COND:{
                return Conductivity(T,p,psi_w);
            }
            case GIVEN_CPHA:{
                double v_bar1,v_bar2,h_bar1,h_bar2, cp_bar, dT = 1e-3;
                v_bar1=MolarVolume(T-dT,p,psi_w); //[m^3/mol_ha]
                h_bar1=MolarEnthalpy(T-dT,p,psi_w,v_bar1); //[kJ/kmol_ha]
                v_bar2=MolarVolume(T+dT,p,psi_w); //[m^3/mol_ha]
                h_bar2=MolarEnthalpy(T+dT,p,psi_w,v_bar2); //[kJ/kmol_ha]
                W=HumidityRatio(psi_w); //[kg_w/kg_da]
                cp_bar = (h_bar2-h_bar1)/(2*dT); //[J/mol_da/K]
                return cp_bar*(1+W)/M_ha; //[J/kg_ha/K]
            }
            case GIVEN_CP:{
                double v_bar1,v_bar2,h_bar1,h_bar2, cp_bar, dT = 1e-3;
                v_bar1=MolarVolume(T-dT,p,psi_w); //[m^3/mol_ha]
                h_bar1=MolarEnthalpy(T-dT,p,psi_w,v_bar1); //[kJ/kmol_ha]
                v_bar2=MolarVolume(T+dT,p,psi_w); //[m^3/mol_ha]
                h_bar2=MolarEnthalpy(T+dT,p,psi_w,v_bar2); //[kJ/kmol_ha]
                W=HumidityRatio(psi_w); //[kg_w/kg_da]
                cp_bar = (h_bar2-h_bar1)/(2*dT); //[J/mol_da/K]
                return cp_bar/M_ha; //[J/kg_da/K]
            }
            default:
                return _HUGE;
        }
    }
    catch (std::exception &e)
    {
        CoolProp::set_error_string(e.what());
        return _HUGE;
    }
    catch (...)
    {
        return _HUGE;
    }
}

double HAProps_Aux(const char* Name,double T, double p, double W, char *units)
{
    // This function provides some things that are not usually needed, but could be interesting for debug purposes.
    
    // Add a check to make sure that Air and Water fluid states have been properly instantiated
    check_fluid_instantiation();

    // Requires W since it is nice and fast and always defined.  Put a dummy value if you want something that doesn't use humidity

    // Takes temperature, pressure, and humidity ratio W as inputs;
    double psi_w,B_aa,C_aaa,B_ww,C_www,B_aw,C_aaw,C_aww,v_bar;

    try{
    if (!strcmp(Name,"Baa"))
    {
        B_aa=B_Air(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"Caaa"))
    {
        C_aaa=C_Air(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"Bww"))
    {
        B_ww=B_Water(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"Cwww"))
    {
        C_www=C_Water(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_www;
    }
    else if (!strcmp(Name,"dBaa"))
    {
        B_aa=dBdT_Air(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aa;
    }
    else if (!strcmp(Name,"dCaaa"))
    {
        C_aaa=dCdT_Air(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaa;
    }
    else if (!strcmp(Name,"dBww"))
    {
        B_ww=dBdT_Water(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_ww;
    }
    else if (!strcmp(Name,"dCwww"))
    {
        C_www=dCdT_Water(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_www;
    }
    else if (!strcmp(Name,"Baw"))
    {
        B_aw=_B_aw(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return B_aw;
    }
    else if (!strcmp(Name,"Caww"))
    {
        C_aww=_C_aww(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aww;
    }
    else if (!strcmp(Name,"Caaw"))
    {
        C_aaw=_C_aaw(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return C_aaw;
    }
    else if (!strcmp(Name,"dBaw"))
    {
        double dB_aw=_dB_aw_dT(T); // [m^3/mol]
        strcpy(units,"m^3/mol");
        return dB_aw;
    }
    else if (!strcmp(Name,"dCaww"))
    {
        double dC_aww=_dC_aww_dT(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return dC_aww;
    }
    else if (!strcmp(Name,"dCaaw"))
    {
        double dC_aaw=_dC_aaw_dT(T); // [m^6/mol^2]
        strcpy(units,"m^6/mol^2");
        return dC_aaw;
    }
    else if (!strcmp(Name,"beta_H"))
    {
        strcpy(units,"1/Pa");
        return HenryConstant(T);
    }
    else if (!strcmp(Name,"kT"))
    {
        strcpy(units,"1/Pa");
        if (T>273.16)
        {
            Water->update(CoolProp::PT_INPUTS, p, T);
            return Water->keyed_output(CoolProp::iisothermal_compressibility);
        }
        else
            return IsothermCompress_Ice(T,p); //[1/Pa]
    }
    else if (!strcmp(Name,"p_ws"))
    {
        strcpy(units,"Pa");
        if (T>273.16)
        {
            Water->update(CoolProp::QT_INPUTS, 0, T);
            return Water->keyed_output(CoolProp::iP);
        }
        else
            return psub_Ice(T);
    }
    else if (!strcmp(Name,"vbar_ws"))
    {
        strcpy(units,"m^3/mol");
        if (T>273.16)
        {
            Water->update(CoolProp::QT_INPUTS, 0, T);
            return 1.0/Water->keyed_output(CoolProp::iDmolar);
        }
        else
        {
            // It is ice
            return dg_dp_Ice(T,p)*MM_Water()/1000/1000; //[m^3/mol]
        }
    }
    else if (!strcmp(Name,"f"))
    {
        strcpy(units,"-");
        return f_factor(T,p);
    }
    // Get psi_w since everything else wants it
    psi_w=MoleFractionWater(T,p,GIVEN_HUMRAT,W);
    if (!strcmp(Name,"Bm"))
    {
        strcpy(units,"m^3/mol");
        return B_m(T,psi_w);
    }
    else if (!strcmp(Name,"Cm"))
    {
        strcpy(units,"m^6/mol^2");
        return C_m(T,psi_w);
    }
    else if (!strcmp(Name,"hvirial"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return 8.3145*T*((B_m(T,psi_w)-T*dB_m_dT(T,psi_w))/v_bar+(C_m(T,psi_w)-T/2.0*dC_m_dT(T,psi_w))/(v_bar*v_bar));
    }
    //else if (!strcmp(Name,"ha"))
    //{
    //    delta=1.1/322; tau=132/T;
    //    return 1+tau*DerivTerms("dphi0_dTau",tau,delta,"Water");
    //}
    //else if (!strcmp(Name,"hw"))
    //{
    //    //~ return Props('D','T',T,'P',p,"Water")/322; tau=647/T;
    //    delta=1000/322; tau=647/T;
    //    //~ delta=rho_Water(T,p,TYPE_TP);tau=647/T;
    //    return 1+tau*DerivTerms("dphi0_dTau",tau,delta,"Water");
    //}
    else if (!strcmp(Name,"hbaro_w"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return IdealGasMolarEnthalpy_Water(T,v_bar);
    }
    else if (!strcmp(Name,"hbaro_a"))
    {
        v_bar=MolarVolume(T,p,psi_w);
        return IdealGasMolarEnthalpy_Air(T,v_bar);
    }
    else
    {
        printf("Sorry I didn't understand your input [%s] to HAProps_Aux\n",Name);
        return -1;
    }
    }
    catch(std::exception &)
    {
        return _HUGE;
    }
    return _HUGE;
}
double cair_sat(double T)
{
    // Air saturation specific heat
    // Based on a correlation from EES, good from 250K to 300K.
    // No error bound checking is carried out
    // T: [K]
    // cair_s: [kJ/kg-K]
    return 2.14627073E+03-3.28917768E+01*T+1.89471075E-01*T*T-4.86290986E-04*T*T*T+4.69540143E-07*T*T*T*T;
}

double IceProps(const char* Name, double T, double p)
{
    if (!strcmp(Name,"s"))
    {
        return s_Ice(T,p*1000.0);
    }
    else if (!strcmp(Name,"rho"))
    {
        return rho_Ice(T,p*1000.0);
    }
    else if (!strcmp(Name,"h"))
    {
        return h_Ice(T,p*1000.0);
    }
    else
    {
        return 1e99;
    }
}

} /* namespace HumidAir */

#ifdef ENABLE_CATCH
#include <math.h>
#include "catch.hpp"

TEST_CASE("Check HA Virials from Table A.2.1","[RP1485]")
{
    SECTION("B_aa")
    {
        CHECK(std::abs(HumidAir::B_Air(-60+273.15)/(-33.065/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Air(0+273.15)/(-13.562/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Air(200+273.15)/(11.905/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Air(350+273.15)/(18.949/1e6)-1) < 1e-3);
    }
    SECTION("B_ww")
    {
        CHECK(std::abs(HumidAir::B_Water(-60+273.15)/(-11174/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Water(0+273.15)/(-2025.6/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Water(200+273.15)/(-200.52/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::B_Water(350+273.15)/(-89.888/1e6)-1) < 1e-3);
    }
    SECTION("B_aw")
    {
        CHECK(std::abs(HumidAir::_B_aw(-60+273.15)/(-68.306/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_B_aw(0+273.15)/(-38.074/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_B_aw(200+273.15)/(-2.0472/1e6)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_B_aw(350+273.15)/(7.5200/1e6)-1) < 1e-3);
    }

    SECTION("C_aaa")
    {
        CHECK(std::abs(HumidAir::C_Air(-60+273.15)/(2177.9/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::C_Air(0+273.15)/(1893.1/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::C_Air(200+273.15)/(1551.2/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::C_Air(350+273.15)/(1464.7/1e12)-1) < 1e-3);
    }
    SECTION("C_www")
    {
        CHECK(std::abs(HumidAir::C_Water(-60+273.15)/(-1.5162999202e-04)-1) < 1e-3); // Relaxed criterion for this parameter
        CHECK(std::abs(HumidAir::C_Water(0+273.15)/(-10981960/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::C_Water(200+273.15)/(-0.00000003713759442)-1) < 1e-3);
        CHECK(std::abs(HumidAir::C_Water(350+273.15)/(-0.000000001198914198)-1) < 1e-3);
    }
    SECTION("C_aaw")
    {
        CHECK(std::abs(HumidAir::_C_aaw(-60+273.15)/(1027.3/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aaw(0+273.15)/(861.02/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aaw(200+273.15)/(627.15/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aaw(350+273.15)/(583.79/1e12)-1) < 1e-3);
    }
    SECTION("C_aww")
    {
        CHECK(std::abs(HumidAir::_C_aww(-60+273.15)/(-1821432/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aww(0+273.15)/(-224234/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aww(200+273.15)/(-8436.5/1e12)-1) < 1e-3);
        CHECK(std::abs(HumidAir::_C_aww(350+273.15)/(-2486.9/1e12)-1) < 1e-3);
    }
}
TEST_CASE("Enhancement factor from Table A.3","[RP1485]")
{
    CHECK(std::abs(HumidAir::f_factor(-60+273.15,101325)/(1.00708)-1) < 1e-3);
    CHECK(std::abs(HumidAir::f_factor( 80+273.15,101325)/(1.00573)-1) < 1e-3);
    CHECK(std::abs(HumidAir::f_factor(-60+273.15,10000e3)/(2.23918)-1) < 1e-3);
    CHECK(std::abs(HumidAir::f_factor(300+273.15,10000e3)/(1.04804)-1) < 1e-3);
}
TEST_CASE("Isothermal compressibility from Table A.5","[RP1485]")
{
    CHECK(std::abs(HumidAir::isothermal_compressibility(-60+273.15,101325)/(0.10771e-9)-1) < 1e-3);
    CAPTURE(HumidAir::isothermal_compressibility( 80+273.15,101325));
    CHECK(std::abs(HumidAir::isothermal_compressibility( 80+273.15,101325)/(0.46009e-9)-1) < 1e-2);  // Relaxed criterion for this parameter
    CHECK(std::abs(HumidAir::isothermal_compressibility(-60+273.15,10000e3)/(0.10701e-9)-1) < 1e-3);
    CHECK(std::abs(HumidAir::isothermal_compressibility(300+273.15,10000e3)/(3.05896e-9)-1) < 1e-3);
}
TEST_CASE("Henry constant from Table A.6","[RP1485]")
{
    CHECK(std::abs(HumidAir::HenryConstant(0.010001+273.15)/(0.22600e-9)-1) < 1e-3);
    CHECK(std::abs(HumidAir::HenryConstant(300+273.15)/(0.58389e-9)-1) < 1e-3);
}

// A structure to hold the values for one call to HAProps
struct hel
{
public:
    std::string in1,in2,in3,out;
    double v1, v2, v3, expected;
    hel(std::string in1, double v1, std::string in2, double v2, std::string in3, double v3, std::string out, double expected)
    {
        this->in1 = in1; this->in2 = in2; this->in3 = in3;
        this->v1 = v1; this->v2 = v2; this->v3 = v3;
        this->expected = expected; this->out = out;
    };
};
hel table_A11[] ={hel("T",473.15,"W",0.00,"P",101325,"B",45.07+273.15),
                  hel("T",473.15,"W",0.00,"P",101325,"V",1.341),
                  hel("T",473.15,"W",0.00,"P",101325,"H",202520),
                  hel("T",473.15,"W",0.00,"P",101325,"S",555.8),
                  hel("T",473.15,"W",0.50,"P",101325,"B",81.12+273.15),
                  hel("T",473.15,"W",0.50,"P",101325,"V",2.416),
                  hel("T",473.15,"W",0.50,"P",101325,"H",1641400),
                  hel("T",473.15,"W",0.50,"P",101325,"S",4829.5),
                  hel("T",473.15,"W",1.00,"P",101325,"B",88.15+273.15),
                  hel("T",473.15,"W",1.00,"P",101325,"V",3.489),
                  hel("T",473.15,"W",1.00,"P",101325,"H",3079550),
                  hel("T",473.15,"W",1.00,"P",101325,"S",8889.0)};

hel table_A12[] ={hel("T",473.15,"W",0.00,"P",1e6,"B",90.47+273.15),
                  hel("T",473.15,"W",0.00,"P",1e6,"V",0.136),
                  hel("T",473.15,"W",0.00,"P",1e6,"H",201940),
                  hel("T",473.15,"W",0.00,"P",1e6,"S",-101.1), // Using CoolProp 4.2, this value seems incorrect from report
                  hel("T",473.15,"W",0.50,"P",1e6,"B",148.49+273.15),
                  hel("T",473.15,"W",0.50,"P",1e6,"V",0.243),
                  hel("T",473.15,"W",0.50,"P",1e6,"H",1630140),
                  hel("T",473.15,"W",0.50,"P",1e6,"S",3630.2),
                  hel("T",473.15,"W",1.00,"P",1e6,"B",159.92+273.15),
                  hel("T",473.15,"W",1.00,"P",1e6,"V",0.347),
                  hel("T",473.15,"W",1.00,"P",1e6,"H",3050210),
                  hel("T",473.15,"W",1.00,"P",1e6,"S",7141.3)};

hel table_A15[] ={hel("T",473.15,"W",0.10,"P",1e7,"B",188.92+273.15),
                  hel("T",473.15,"W",0.10,"P",1e7,"V",0.016),
                  hel("T",473.15,"W",0.10,"P",1e7,"H",473920),
                  hel("T",473.15,"W",0.10,"P",1e7,"S",-90.1),
                  hel("T",473.15,"W",0.10,"P",1e7,"R",0.734594),
                  };

class HAPropsConsistencyFixture
{
public:
    std::vector<hel> inputs;
    std::string in1,in2,in3,out;
    double v1, v2, v3, expected, actual;
    void set_table(hel h[], int nrow){
        int h1 = sizeof(h), h2 = sizeof(h[0]);
        inputs = std::vector<hel>(h, h + nrow);
    };
    void set_values(hel &h){
        this->in1 = h.in1; this->in2 = h.in2; this->in3 = h.in3;
        this->v1 = h.v1; this->v2 = h.v2; this->v3 = h.v3;
        this->expected = h.expected; this->out = h.out;
    };
    void call(){
        actual = HumidAir::HAPropsSI(out.c_str(), in1.c_str(), v1, in2.c_str(), v2, in3.c_str(), v3);
    }
};

TEST_CASE_METHOD(HAPropsConsistencyFixture, "ASHRAE RP1485 Tables", "[RP1485]")
{
    SECTION("Table A.15")
    {
        set_table(table_A15, 5);
        for (std::size_t i = 0; i < inputs.size(); ++i){
            set_values(inputs[i]);
            call();
            CAPTURE(inputs[i].in1);
            CAPTURE(inputs[i].v1);
            CAPTURE(inputs[i].in2);
            CAPTURE(inputs[i].v2);
            CAPTURE(inputs[i].in3);
            CAPTURE(inputs[i].v3);
            CAPTURE(out);
            CAPTURE(actual);
            CAPTURE(expected);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(std::abs(actual/expected-1) < 0.01);
        }
    }
    SECTION("Table A.11")
    {
        set_table(table_A11, 12);
        for (std::size_t i = 0; i < inputs.size(); ++i){
            set_values(inputs[i]);
            call();
            CAPTURE(inputs[i].in1);
            CAPTURE(inputs[i].v1);
            CAPTURE(inputs[i].in2);
            CAPTURE(inputs[i].v2);
            CAPTURE(inputs[i].in3);
            CAPTURE(inputs[i].v3);
            CAPTURE(out);
            CAPTURE(actual);
            CAPTURE(expected);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(std::abs(actual/expected-1) < 0.01);
        }
    }
    SECTION("Table A.12")
    {
        set_table(table_A12, 12);
        for (std::size_t i = 0; i < inputs.size(); ++i){
            set_values(inputs[i]);
            call();
            CAPTURE(inputs[i].in1);
            CAPTURE(inputs[i].v1);
            CAPTURE(inputs[i].in2);
            CAPTURE(inputs[i].v2);
            CAPTURE(inputs[i].in3);
            CAPTURE(inputs[i].v3);
            CAPTURE(out);
            CAPTURE(actual);
            CAPTURE(expected);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(std::abs(actual/expected-1) < 0.01);
        }
    }
}
TEST_CASE("Assorted tests","[HAPropsSI]")
{
    CHECK(ValidNumber(HumidAir::HAPropsSI("T", "H", 267769, "P", 104300, "W", 0.0)));
    CHECK(ValidNumber(HumidAir::HAPropsSI("T", "B", 252.84, "W", 5.097e-4, "P", 101325)));
    CHECK(ValidNumber(HumidAir::HAPropsSI("T", "B",290, "R", 1, "P", 101325)));
}
// a predicate implemented as a function:
bool is_not_a_pair (const std::set<std::size_t> &item) { return item.size() != 2; }

const int number_of_inputs = 6;
std::string inputs[number_of_inputs] = {"W","D","B","R","T","V"};//,"H","S"};

class ConsistencyTestData
{
public:
	bool is_built;
	std::vector<Dictionary> data;
	std::list<std::set<std::size_t> > inputs_list;
	ConsistencyTestData(){
		is_built = false;
	};
	void build(){
		if (is_built){return;}
		std::vector<std::size_t> indices(number_of_inputs);
		for (std::size_t i = 0; i < number_of_inputs; ++i){ indices[i] = i;}
		// Generate a powerset of all the permutations of all lengths of inputs
		std::set<std::size_t> indices_set(indices.begin(), indices.end());
		std::set<std::set<std::size_t> > inputs_powerset = powerset(indices_set);
		inputs_list = std::list<std::set<std::size_t> >(inputs_powerset.begin(), inputs_powerset.end());
		inputs_list.remove_if(is_not_a_pair);

		const int NT = 10, NW = 5;
		double p = 101325;
		for (double T = 210; T < 350; T += (350-210)/(NT-1))
		{
			double Wsat = HumidAir::HAPropsSI("W", "T", T, "P", p, "R", 1.0);
			for (double W = 1e-5; W < Wsat; W += (Wsat-1e-5)/(NW-1)){
				Dictionary vals;
				// Calculate all the values using T, W
				for (int i = 0; i < number_of_inputs; ++i){
					double v = HumidAir::HAPropsSI(inputs[i], "T", T, "P", p, "W", W);
					vals.add_number(inputs[i], v);
				}
				data.push_back(vals);
				std::cout << format("T %g W %g\n",T,W);
			}
		}
		is_built = true;
	};
} consistency_data;

/*
 * This test is incredibly slow, which is why it is currently commented out.  Many of the tests also fail
 * 
TEST_CASE("HAPropsSI", "[HAPropsSI]")
{
	consistency_data.build();
	double p = 101325;
	for (std::size_t i = 0; i < consistency_data.data.size(); ++i)
	{
		for (std::list<std::set<std::size_t> >::iterator iter = consistency_data.inputs_list.begin(); iter != consistency_data.inputs_list.end(); ++iter)
		{
			std::vector<std::size_t> pair(iter->begin(), iter->end());
			std::string i0 = inputs[pair[0]], i1 = inputs[pair[1]];
			double v0 = consistency_data.data[i].get_double(i0), v1 = consistency_data.data[i].get_double(i1);
            if ((i0 == "B" && i1 == "V") || (i1 == "B" && i0 == "V")){continue;}
			std::ostringstream ss2;
			ss2 << "Inputs: \"" << i0 << "\"," << v0 << ",\"" << i1 << "\"," << v1;
			SECTION(ss2.str(), ""){

				double T = consistency_data.data[i].get_double("T");
				double W = consistency_data.data[i].get_double("W");
				double Wcalc = HumidAir::HAPropsSI("W", i0, v0, i1, v1, "P", p);
				double Tcalc = HumidAir::HAPropsSI("T", i0, v0, i1, v1, "P", p);
				std::string err = CoolProp::get_global_param_string("errstring");
				CAPTURE(T);
				CAPTURE(W);
				CAPTURE(Tcalc);
				CAPTURE(Wcalc);
				CAPTURE(err);
				CHECK(std::abs(Tcalc - T) < 1e-1);
				CHECK(std::abs((Wcalc - W)/W) < 1e-3);
			}
		}
	}
}
 */

#endif /* CATCH_ENABLED */

