/*!
  \file include.h
  \brief Main include file

  This is a main include file for the entire CoolProp
  project. It defines some important preprocessor variables that might
  have to be changed by the user.

  Uncomment the define directives as appropriate

  Ian Bell (ian.bell@ulg.ac.be)
  2012-2013
  University of Liege, Liege, Belgium

  Francesco Casella, Christoph Richter, Roberto Bonifetto
  2006-2012
  Copyright Politecnico di Milano, TU Braunschweig, Politecnico di Torino
*/
#ifndef INCLUDE_H_
#define INCLUDE_H_

/********************************************************************
 *                Start of user option selection
 ********************************************************************/

// Selection of Modelica compiler
//! Modelica compiler is Dymola
/*!
  Set this preprocessor variable to 1 if Dymola is the Modelica
  compiler that is going to be used with the compiled library.
  \sa OPEN_MODELICA
*/
#define DYMOLA 1

//! Modelica compiler is OpenModelica
/*!
  Set this preprocessor variable to 1 if OpenModelica is the Modelica
  compiler that is going to be used with the compiled library.
  \sa DYMOLA
*/
#define OPEN_MODELICA 0

// Selection of used external fluid property computation packages.
//! FluidProp solver
/*!
  Set this preprocessor variable to 1 to include the interface to the
  FluidProp solver developed and maintained by Francesco Casella.
*/
#if defined(WIN32) || defined(_WIN32)
#    define FLUIDPROP 1
#else
#    define FLUIDPROP 0
#endif

// Selection of used external fluid property computation packages.
//! CoolProp solver
/*!
  Set this preprocessor variable to 1 to include the interface to the
  CoolProp solver developed and maintained by Ian Bell (ian.h.bell@gmail.com).
*/
#define COOLPROP 1

// Selection of build type for this project
//! Build project into a DLL
/*!
  Set this preprocessor variable to 1 if the project is built into a
  dynamic link library. This setting influences the error reporting
  mechanism as well as the export statement.
*/
#define BUILD_DLL 0

//! Not a number
/*!
  This value is used as not a number value. It can be changed by
  the user if there is a more appropriate value.
*/
#define NAN 0xffffffff
#define ISNAN(x) (x == NAN)

/********************************************************************
 *                 End of user option selection
 *            Do not change anything below this line
 ********************************************************************/

// General purpose includes
#include <map>
using std::map;

#include <string>
using std::string;

// Include error handling
#include "errorhandling.h"

#endif /*INCLUDE_H_*/
