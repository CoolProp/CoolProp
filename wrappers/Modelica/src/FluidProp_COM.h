/*!
  Microsoft Visual C++ 2005 Client

  This is an example of a client application in Microsoft Visual C++ 2005 for FluidProp,
  a COM server module for the calculation of fluid properties. FluidProp is a common
  interface to GasMix, IF97, Refprop, StanMix, TPSI and is developed by Piero Colonna
  and Teus van der Stelt.

  The class defined in this file, IFluidProp_COM is the direct interface to the
  FluidProp COM server. This class should not be altered.

  Teus van der Stelt
  Energy Technology Section
  TU Delft

  July, 2004, for FluidProp 1
  January, 2006, for FluidProp 2
  April, 2007, for FluidProp 2.3
*/

#ifndef FluidProp_COM_h
#define FluidProp_COM_h

#ifndef comutil_h
#    include <comutil.h>
#endif

// An IFluidProp interface
//
interface IFluidProp_COM : public IDispatch {
   public:
    virtual void __stdcall CreateObject(BSTR ModelName, BSTR * ErrorMsg) = 0;
    virtual void __stdcall ReleaseObjects() = 0;

    virtual void __stdcall SetFluid(BSTR ModelName, long nComp, SAFEARRAY** sa_Comp, SAFEARRAY** sa_Conc, BSTR* ErrorMsg) = 0;
    virtual void __stdcall SetFluid_M(BSTR ModelName, long nComp, SAFEARRAY* sa_Comp, SAFEARRAY* sa_Conc, BSTR* ErrorMsg) = 0;

    virtual void __stdcall GetFluid(BSTR ModelName, long* nComp, SAFEARRAY** sa_Comp, SAFEARRAY** sa_Conc) = 0;
    virtual void __stdcall GetFluid_M(BSTR ModelName, long* nComp, SAFEARRAY** sa_Comp, SAFEARRAY** sa_Conc) = 0;

    virtual void __stdcall GetFluidNames(BSTR LongShort, BSTR ModelName, long* nComp, SAFEARRAY** sa_CompSet, BSTR* ErrorMsg) = 0;
    virtual void __stdcall GetFluidNames_M(BSTR LongShort, BSTR ModelName, long* nComp, SAFEARRAY** sa_CompSet, BSTR* ErrorMsg) = 0;

    virtual void __stdcall GetCompSet(BSTR ModelName, long* nComp, SAFEARRAY** sa_CompSet, BSTR* ErrorMsg) = 0;
    virtual void __stdcall GetCompSet_M(BSTR ModelName, long* nComp, SAFEARRAY** sa_CompSet, BSTR* ErrorMsg) = 0;

    virtual void __stdcall Pressure(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Temperature(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall SpecVolume(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Density(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Enthalpy(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Entropy(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall IntEnergy(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall VaporQual(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall LiquidCmp(BSTR InputSpec, double Input1, double Input2, SAFEARRAY** Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall LiquidCmp_M(BSTR InputSpec, double Input1, double Input2, SAFEARRAY** Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall VaporCmp(BSTR InputSpec, double Input1, double Input2, SAFEARRAY** Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall VaporCmp_M(BSTR InputSpec, double Input1, double Input2, SAFEARRAY** Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall HeatCapV(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall HeatCapP(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall SoundSpeed(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Alpha(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Beta(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Chi(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Fi(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Ksi(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Psi(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Zeta(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Theta(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Kappa(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Gamma(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall Viscosity(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall ThermCond(BSTR InputSpec, double Input1, double Input2, double* Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall AllProps(BSTR InputSpec, double Input1, double Input2, double* P, double* T, double* v, double* d, double* h, double* s,
                                    double* u, double* q, SAFEARRAY** x, SAFEARRAY** y, double* cv, double* cp, double* c, double* alpha,
                                    double* beta, double* chi, double* fi, double* ksi, double* psi, double* zeta, double* theta, double* kappa,
                                    double* gamma, double* eta, double* lambda, BSTR* ErrorMsg) = 0;
    virtual void __stdcall AllProps_M(BSTR InputSpec, double Input1, double Input2, double* P, double* T, double* v, double* d, double* h, double* s,
                                      double* u, double* q, SAFEARRAY** x, SAFEARRAY** y, double* cv, double* cp, double* c, double* alpha,
                                      double* beta, double* chi, double* fi, double* ksi, double* psi, double* zeta, double* theta, double* kappa,
                                      double* gamma, double* eta, double* lambda, BSTR* ErrorMsg) = 0;

    virtual void __stdcall AllPropsSat(BSTR InputSpec, double Input1, double Input2, double* P, double* T, double* v, double* d, double* h, double* s,
                                       double* u, double* q, SAFEARRAY** x, SAFEARRAY** y, double* cv, double* cp, double* c, double* alpha,
                                       double* beta, double* chi, double* fi, double* ksi, double* psi, double* zeta, double* theta, double* kappa,
                                       double* gamma, double* eta, double* lambda, double* d_liq, double* d_vap, double* h_liq, double* h_vap,
                                       double* T_sat, double* dd_liq_dP, double* dd_vap_dP, double* dh_liq_dP, double* dh_vap_dP, double* dT_sat_dP,
                                       BSTR* ErrorMsg) = 0;
    virtual void __stdcall AllPropsSat_M(BSTR InputSpec, double Input1, double Input2, double* P, double* T, double* v, double* d, double* h,
                                         double* s, double* u, double* q, SAFEARRAY** x, SAFEARRAY** y, double* cv, double* cp, double* c,
                                         double* alpha, double* beta, double* chi, double* fi, double* ksi, double* psi, double* zeta, double* theta,
                                         double* kappa, double* gamma, double* eta, double* lambda, double* d_liq, double* d_vap, double* h_liq,
                                         double* h_vap, double* T_sat, double* dd_liq_dP, double* dd_vap_dP, double* dh_liq_dP, double* dh_vap_dP,
                                         double* dT_sat_dP, BSTR* ErrorMsg) = 0;

    virtual void __stdcall Solve(BSTR FuncSpec, double FuncVal, BSTR InputSpec, long Target, double FixedVal, double MinVal, double MaxVal,
                                 double* Output, BSTR* ErrorMsg) = 0;

    virtual void __stdcall Mmol(double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Tcrit(double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Pcrit(double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Tmin(double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall Tmax(double* Output, BSTR* ErrorMsg) = 0;
    virtual void __stdcall AllInfo(double* M_mol, double* T_crit, double* P_crit, double* T_min, double* T_max, BSTR* ErrorMsg) = 0;

    virtual void __stdcall SetUnits(BSTR UnitSet, BSTR MassOrMole, BSTR Properties, BSTR Units, BSTR * ErrorMsg) = 0;
    virtual void __stdcall SetRefState(double T_ref, double P_ref, BSTR* ErrorMsg) = 0;
};

#endif
