#ifndef FLUIDCONSTANTS_H_
#define FLUIDCONSTANTS_H_

#include "include.h"

//! Fluid constants struct
/*!
  The fluid constants struct contains all the constant fluid properties
  that are returned by the external solver.

  Francesco Casella, Christoph Richter, Roberto Bonifetto
  2006-2012
  Copyright Politecnico di Milano, TU Braunschweig, Politecnico
  di Torino
*/

struct FluidConstants
{
    //! Molar mass
    double MM;
    //! Pressure at critical point
    double pc;
    //! Temperature at critical point
    double Tc;
    //! Density at critical point
    double dc;
    // The following two functions are currently only available internally
    // but do not have the required interface functions to be accessible from
    // Modelica.
    //! Specific enthalpy at critical point
    double hc;
    //! Specific entropy at critical point
    double sc;
};

#endif  // FLUIDCONSTANTS_H_
