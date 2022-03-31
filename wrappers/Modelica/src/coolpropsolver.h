#ifndef COOLPROPSOLVER_H_
#define COOLPROPSOLVER_H_

#include "basesolver.h"

//! CoolProp solver class
/*!
  This class defines a solver that calls out to the open-source CoolProp property database.

  libraryName = "CoolProp";

  Ian Bell (ian.h.bell@gmail.com)
  2012-2013
  University of Liege, Liege, Belgium
*/
class CoolPropSolver : public BaseSolver
{
   protected:
    class CoolPropStateClassSI* state;
    bool enable_TTSE, calc_transport, extend_twophase;
    int debug_level;
    double twophase_derivsmoothing_xend;
    double rho_smoothing_xend;
    long fluidType;

    virtual void preStateChange(void);
    virtual void postStateChange(ExternalThermodynamicState* const properties);

   public:
    CoolPropSolver(const std::string& mediumName, const std::string& libraryName, const std::string& substanceName);
    ~CoolPropSolver(){};
    virtual void setFluidConstants();

    virtual void setSat_p(double& p, ExternalSaturationProperties* const properties);
    virtual void setSat_T(double& T, ExternalSaturationProperties* const properties);

    virtual void setState_ph(double& p, double& h, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_pT(double& p, double& T, ExternalThermodynamicState* const properties);
    virtual void setState_dT(double& d, double& T, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_ps(double& p, double& s, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_hs(double& h, double& s, int& phase, ExternalThermodynamicState* const properties);

    virtual double Pr(ExternalThermodynamicState* const properties);
    virtual double T(ExternalThermodynamicState* const properties);
    virtual double a(ExternalThermodynamicState* const properties);
    virtual double beta(ExternalThermodynamicState* const properties);
    virtual double cp(ExternalThermodynamicState* const properties);
    virtual double cv(ExternalThermodynamicState* const properties);
    virtual double d(ExternalThermodynamicState* const properties);
    virtual double ddhp(ExternalThermodynamicState* const properties);
    virtual double ddph(ExternalThermodynamicState* const properties);
    virtual double eta(ExternalThermodynamicState* const properties);
    virtual double h(ExternalThermodynamicState* const properties);
    virtual double kappa(ExternalThermodynamicState* const properties);
    virtual double lambda(ExternalThermodynamicState* const properties);
    virtual double p(ExternalThermodynamicState* const properties);
    virtual int phase(ExternalThermodynamicState* const properties);
    virtual double s(ExternalThermodynamicState* const properties);
    virtual double d_der(ExternalThermodynamicState* const properties);
    virtual double isentropicEnthalpy(double& p, ExternalThermodynamicState* const properties);

    virtual double dTp(ExternalSaturationProperties* const properties);
    virtual double ddldp(ExternalSaturationProperties* const properties);
    virtual double ddvdp(ExternalSaturationProperties* const properties);
    virtual double dhldp(ExternalSaturationProperties* const properties);
    virtual double dhvdp(ExternalSaturationProperties* const properties);
    virtual double dl(ExternalSaturationProperties* const properties);
    virtual double dv(ExternalSaturationProperties* const properties);
    virtual double hl(ExternalSaturationProperties* const properties);
    virtual double hv(ExternalSaturationProperties* const properties);
    virtual double sigma(ExternalSaturationProperties* const properties);
    virtual double sl(ExternalSaturationProperties* const properties);
    virtual double sv(ExternalSaturationProperties* const properties);

    virtual double psat(ExternalSaturationProperties* const properties);
    virtual double Tsat(ExternalSaturationProperties* const properties);
};

#endif  // COOLPROPSOLVER_H_
