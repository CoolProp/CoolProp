/*!
  Microsoft Visual C++ 2005

  This is an example of a client application in Microsoft Visual C++ 2005 for FluidProp,
  a COM server module for the calculation of fluid properties. FluidProp is a common
  interface to GasMix, IF97, Refprop, StanMix, TPSI and is developed by Piero Colonna
  and Teus van der Stelt.

  The class defined in this file, TFluidProp, is as a wrapper class for the base
  class IFluidProp_COM. TFluidProp hides specific COM server details like safe arrays
  (SAFEARRAY) and binary strings (BSTR) in IFluidProp_COM. In the TFluidProp class
  only standard C++ data types are used. This seriously facilitates working with the
  COM server.

  Teus van der Stelt
  Energy Technology Section
  TU Delft

  July, 2004, for FluidProp 1
  January, 2006, for FluidProp 2
  April, 2007, for FluidProp 2.3
*/

#ifndef FluidProp_IF_h
#define FluidProp_IF_h

#pragma comment(lib, "comsuppw.lib")

#include <string>
using namespace std;

#ifndef FluidProp_COM_h
#    include "FluidProp_COM.h"
#endif

// An TFluidProp interface
//
class TFluidProp
{
   public:
    TFluidProp();
    ~TFluidProp();

    void CreateObject(string ModelName, string* ErrorMsg);
    void ReleaseObjects();

    void SetFluid(string ModelName, int nComp, string* Comp, double* Conc, string* ErrorMsg);
    void GetFluid(string* ModelName, int* nComp, string* Comp, double* Conc, bool CompInfo = true);
    void GetFluidNames(string LongShort, string ModelName, int* nFluids, string* FluidNames, string* ErrorMsg);

    double Pressure(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Temperature(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double SpecVolume(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Density(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Enthalpy(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Entropy(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double IntEnergy(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double VaporQual(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double* LiquidCmp(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double* VaporCmp(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double HeatCapV(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double HeatCapP(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double SoundSpeed(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Alpha(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Beta(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Chi(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Fi(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Ksi(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Psi(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Zeta(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Theta(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Kappa(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Gamma(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double Viscosity(string InputSpec, double Input1, double Input2, string* ErrorMsg);
    double ThermCond(string InputSpec, double Input1, double Input2, string* ErrorMsg);

    void AllProps(string InputSpec, double Input1, double Input2, double& P, double& T, double& v, double& d, double& h, double& s, double& u,
                  double& q, double* x, double* y, double& cv, double& cp, double& c, double& alpha, double& beta, double& chi, double& fi,
                  double& ksi, double& psi, double& zeta, double& theta, double& kappa, double& gamma, double& eta, double& lambda, string* ErrorMsg);

    // Compute all the properties at once, including saturation properties
    void AllPropsSat(string InputSpec, double Input1, double Input2, double& P, double& T, double& v, double& d, double& h, double& s, double& u,
                     double& q, double* x, double* y, double& cv, double& cp, double& c, double& alpha, double& beta, double& chi, double& fi,
                     double& ksi, double& psi, double& zeta, double& theta, double& kappa, double& gamma, double& eta, double& lambda, double& d_liq,
                     double& d_vap, double& h_liq, double& h_vap, double& T_sat, double& dd_liq_dP, double& dd_vap_dP, double& dh_liq_dP,
                     double& dh_vap_dP, double& dT_sat_dP, string* ErrorMsg);

    double Solve(string FuncSpec, double FuncVal, string InputSpec, long Target, double FixedVal, double MinVal, double MaxVal, string* ErrorMsg);

    double Mmol(string* ErrorMsg);
    double Tcrit(string* ErrorMsg);
    double Pcrit(string* ErrorMsg);
    double Tmin(string* ErrorMsg);
    double Tmax(string* ErrorMsg);
    void AllInfo(double& Mmol, double& Tcrit, double& Pcrit, double& Tmin, double& Tmax, string* ErrorMsg);

    void SetUnits(string UnitSet, string MassOrMole, string Properties, string Units, string* ErrorMsg);
    void SetRefState(double T_ref, double P_ref, string* ErrorMsg);

   private:
    IClassFactory* ClassFactory;    // Pointer to class factory
    IFluidProp_COM* FluidProp_COM;  // Pointer to FluidProp interface
};

#endif
