#ifndef BASESOLVER_H_
#define BASESOLVER_H_

#include "include.h"
#include "fluidconstants.h"
#include "CoolPropLib.h"
#include <string>

struct FluidConstants;

//! Base solver class.
/*!
  This is the base class for all external solver objects
  (e.g. TestSolver, FluidPropSolver). A solver object
  encapsulates the interface to external fluid property
  computation routines

  Francesco Casella, Christoph Richter, Roberto Bonifetto
  2006-2012
  Copyright Politecnico di Milano, TU Braunschweig, Politecnico
  di Torino
*/
class BaseSolver
{
   public:
    BaseSolver(const std::string& mediumName, const std::string& libraryName, const std::string& substanceName);
    virtual ~BaseSolver();

    double molarMass() const;
    double criticalTemperature() const;
    double criticalPressure() const;
    double criticalMolarVolume() const;
    double criticalDensity() const;
    double criticalEnthalpy() const;
    double criticalEntropy() const;

    virtual void setFluidConstants();

    virtual void setState_ph(double& p, double& h, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_pT(double& p, double& T, ExternalThermodynamicState* const properties);
    virtual void setState_dT(double& d, double& T, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_ps(double& p, double& s, int& phase, ExternalThermodynamicState* const properties);

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

    virtual void setSat_p(double& p, ExternalSaturationProperties* const properties);
    virtual void setSat_T(double& T, ExternalSaturationProperties* const properties);

    virtual void setBubbleState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const bubbleProperties);
    virtual void setDewState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const bubbleProperties);

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

    virtual bool computeDerivatives(ExternalThermodynamicState* const properties);

    virtual double psat(ExternalSaturationProperties* const properties);
    virtual double Tsat(ExternalSaturationProperties* const properties);

    //! Medium name
    std::string mediumName;
    //! Library name
    std::string libraryName;
    //! Substance name
    std::string substanceName;

   protected:
    //! Fluid constants
    FluidConstants _fluidConstants;
};

#endif  // BASESOLVER_H_
