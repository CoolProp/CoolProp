//! FluidProp solver interface class
/*!
 * This class defines a solver object encapsulating a FluidProp object
 *
 * The class will work if FluidProp is correctly installed, and if
 * the following files, defining the CFluidProp object, are included
 * in the C project:
 *   - FluidProp_IF.h
 *   - FluidProp_IF.cpp
 *   - FluidProp_COM.h
 *  These files are developed and maintained by TU Delft and distributed
 *  as a part of the FluidProp suite http://www.fluidprop.com
 *
 *  Compilation requires support of the COM libraries:
 *  http://en.wikipedia.org/wiki/Component_Object_Model
 *
 * To instantiate a specific FluidProp fluid, it is necessary to set
 * the libraryName and substanceNames package constants as in the
 * following example:
 *
 * libraryName = "FluidProp.RefProp";
 * substanceNames = {"H2O"};
 *
 * Instead of RefProp, it is possible to indicate TPSI, StanMix, etc.
 * Instead of H2O, it is possible to indicate any supported substance
 *
 * See also the solvermap.cpp code
 *
 * Francesco Casella, Christoph Richter, Roberto Bonifetto
 * 2006 - 2012
 * Copyright Politecnico di Milano, TU Braunschweig,
 * Politecnico di Torino
 ********************************************************************/

#ifndef FLUIDPROPSOLVER_H_
#define FLUIDPROPSOLVER_H_

#include "basesolver.h"

#if (FLUIDPROP == 1)

#    include "FluidProp_IF.h"

class FluidPropSolver : public BaseSolver
{
   public:
    FluidPropSolver(const string& mediumName, const string& libraryName, const string& substanceName);
    ~FluidPropSolver();
    virtual void setFluidConstants();

    virtual void setSat_p(double& p, ExternalSaturationProperties* const properties);
    virtual void setSat_T(double& T, ExternalSaturationProperties* const properties);

    virtual void setState_ph(double& p, double& h, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_pT(double& p, double& T, ExternalThermodynamicState* const properties);
    virtual void setState_dT(double& d, double& T, int& phase, ExternalThermodynamicState* const properties);
    virtual void setState_ps(double& p, double& s, int& phase, ExternalThermodynamicState* const properties);
    virtual void setBubbleState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const bubbleProperties);
    virtual void setDewState(ExternalSaturationProperties* const properties, int phase, ExternalThermodynamicState* const dewProperties);
    virtual double isentropicEnthalpy(double& p, ExternalThermodynamicState* const properties);

   protected:
    TFluidProp FluidProp;  // Instance of FluidProp wrapper object
    bool isError(string ErrorMsg);
};

#endif  // FLUIDPROP == 1

#endif /*FLUIDPROPSOLVER_H_*/
