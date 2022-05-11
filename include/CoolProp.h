/*
This header file includes the high level API that is meant to be accessed via C++.  Functions may accept C++ types like std::vector

For the C-style wrapper, refer to CoolPropLib.h

\sa CoolPropLib.h
*/

/*! \mainpage CoolProp Core Code Documentation

Welcome to the home page for the C++ sources of CoolProp.  This information may be useful for developers or just the plain inquisitive

You might want to start by looking at CoolProp.h
*/

#ifndef CoolProp_H
#define CoolProp_H

#include <string>
#include <vector>
#include "DataStructures.h"

namespace CoolProp {

/// Return a value that does not depend on the thermodynamic state - this is a convenience function that does the call PropsSI(Output, "", 0, "", 0, FluidName)
/// @param FluidName The fluid name
/// @param Output The output parameter, one of "Tcrit","D","H",etc.
double Props1SI(std::string FluidName, std::string Output);
/**
     * @brief  Get a matrix of outputs that do not depend on the thermodynamic state - this is a convenience function that does the call PropsSImulti(Outputs, "", {0}, "", {0}, backend, fluids, fractions)
     * @param Outputs A vector of strings for the output parameters
     * @param backend The string representation of the backend (HEOS, REFPROP, INCOMP, etc.)
     * @param fluids The fluid name(s)
     * @param fractions The fractions (molar, mass, volume, etc.) of the components
     */     
std::vector<std::vector<double>> Props1SImulti(const std::vector<std::string>& Outputs, const std::string& backend, const std::vector<std::string>& fluids, const std::vector<double>& fractions);
/// Return a value that depends on the thermodynamic state
/// @param Output The output parameter, one of "T","D","H",etc.
/// @param Name1 The first state variable name, one of "T","D","H",etc.
/// @param Prop1 The first state variable value
/// @param Name2 The second state variable name, one of "T","D","H",etc.
/// @param Prop2 The second state variable value
/// @param FluidName The fluid name
double PropsSI(const std::string& Output, const std::string& Name1, double Prop1, const std::string& Name2, double Prop2,
               const std::string& FluidName);

/**
     * @brief Get a matrix of outputs for a given input.  Can handle both vector inputs as well as a vector of output strings
     * @param Outputs A vector of strings for the output parameters
     * @param Name1 The name of the first input variable
     * @param Prop1 A vector of the first input values
     * @param Name2 The name of the second input variable
     * @param Prop2 A vector of the second input values
     * @param backend The string representation of the backend (HEOS, REFPROP, INCOMP, etc.)
     * @param fluids The fluid name(s)
     * @param fractions The fractions (molar, mass, volume, etc.) of the components
     */
std::vector<std::vector<double>> PropsSImulti(const std::vector<std::string>& Outputs, const std::string& Name1, const std::vector<double>& Prop1,
                                              const std::string& Name2, const std::vector<double>& Prop2, const std::string& backend,
                                              const std::vector<std::string>& fluids, const std::vector<double>& fractions);

/// Get the debug level
/// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
int get_debug_level();
/// Set the debug level
/// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
void set_debug_level(int level);

/// Set the global error string
/// @param error The error string to use
void set_error_string(const std::string& error);
/// An internal function to set the global warning string
/// @param warning The string to set as the warning string
void set_warning_string(const std::string& warning);

/* \brief Extract a value from the saturation ancillary
     *
     * @param fluid_name The name of the fluid to be used - HelmholtzEOS backend only
     * @param output The desired output variable ("P" for instance for pressure)
     * @param Q The quality, 0 or 1
     * @param input The input variable ("T")
     * @param value The input value
     */
double saturation_ancillary(const std::string& fluid_name, const std::string& output, int Q, const std::string& input, double value);

/// Get a globally-defined string
/// @param ParamName A string, one of "version", "errstring", "warnstring", "gitrevision", "FluidsList", "fluids_list", "parameter_list","predefined_mixtures"
/// @returns str The string, or an error message if not valid input
std::string get_global_param_string(const std::string& ParamName);

/*/// Get a long that represents the fluid type
    /// @param FluidName The fluid name as a string
    /// @returns long element from global type enumeration
    long getFluidType(std::string FluidName);*/

/**
    \brief Get a string for a value from a fluid (numerical values for the fluid can be obtained from Props1SI function)

    @param FluidName The name of the fluid that is part of CoolProp, for instance "n-Propane"
    @param ParamName A string, can be in one of the terms described in the following table

    ParamName                    | Description
    --------------------------   | ----------------------------------------
    "aliases"                    | A comma separated list of aliases for the fluid
    "CAS", "CAS_number"          | The CAS number
    "ASHRAE34"                   | The ASHRAE standard 34 safety rating
    "REFPROPName","REFPROP_name" | The name of the fluid used in REFPROP
    "Bibtex-XXX"                 | A BibTeX key, where XXX is one of the bibtex keys used in get_BibTeXKey
    "pure"                       | "true" if the fluid is pure, "false" otherwise
    "formula"                    | The chemical formula of the fluid in LaTeX form if available, "" otherwise

    @returns The string, or an error message if not valid input
    */
std::string get_fluid_param_string(const std::string& FluidName, const std::string& ParamName);

/** \brief Check if the fluid name is valid
     *
     * @returns output Returns true if the fluid string is valid
     *
     * \note "gfreilgregre" -> false; "HEOS::Water" -> true; "Water" -> true
     *
     */
bool is_valid_fluid_string(const std::string& fluidstring);

/** \brief Add fluids as a JSON-formatted string
     *
     * @param backend The backend to which these should be added; e.g. "HEOS", "SRK", "PR"
     * @returns output Returns true if the fluids were able to be added
     *
     */
bool add_fluids_as_JSON(const std::string& backend, const std::string& fluidstring);

/**
    \brief Set the reference state based on a string representation

    @param FluidName The name of the fluid (Backend can be provided like "REFPROP::Water", or if no backend is provided, "HEOS" is the assumed backend)
    @param reference_state The reference state to use, one of

    Reference State | Description
    -------------   | -------------------
    "IIR"           | h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
    "ASHRAE"        | h = 0, s = 0 @ -40C saturated liquid
    "NBP"           | h = 0, s = 0 @ 1.0 bar saturated liquid
    "DEF"           | Reset to the default reference state for the fluid
    "RESET"         | Remove the offset

    The offset in the ideal gas Helmholtz energy can be obtained from
    \f[
    \displaystyle\frac{\Delta s}{R_u/M}+\frac{\Delta h}{(R_u/M)T}\tau
    \f]
    where \f$ \Delta s = s-s_{spec} \f$ and \f$ \Delta h = h-h_{spec} \f$
    */
void set_reference_stateS(const std::string& FluidName, const std::string& reference_state);

/// Set the reference state based on a thermodynamic state point specified by temperature and molar density
/// @param FluidName The name of the fluid
/// @param T Temperature at reference state [K]
/// @param rhomolar Molar density at reference state [mol/m^3]
/// @param hmolar0 Molar enthalpy at reference state [J/mol]
/// @param smolar0 Molar entropy at reference state [J/mol/K]
void set_reference_stateD(const std::string& FluidName, double T, double rhomolar, double hmolar0, double smolar0);

/// Return a string representation of the phase
/// @param Name1 The first state variable name, one of "T","D","H",etc.
/// @param Prop1 The first state variable value
/// @param Name2 The second state variable name, one of "T","D","H",etc.
/// @param Prop2 The second state variable value
/// @param FluidName The fluid name
/// \note Returns empty string if there was an error; use get_global_param_string("errstring") to retrieve the error
std::string PhaseSI(const std::string& Name1, double Prop1, const std::string& Name2, double Prop2, const std::string& FluidName);

/**
     * @brief Extract the backend from a string - something like "HEOS::Water" would split to "HEOS" and "Water".  If no backend is specified, the backend will be set to "?"
     * @param fluid_string The input string
     * @param backend The output backend, if none found, "?"
     * @param fluid The output fluid string (minus the backend string)
     */
void extract_backend(std::string fluid_string, std::string& backend, std::string& fluid);

/**
     * @brief Extract fractions (molar, mass, etc.) encoded in the string if any
     * @param fluid_string The input string
     * @param fractions The fractions
     * @return The fluids, as a '&' delimited string
     */
std::string extract_fractions(const std::string& fluid_string, std::vector<double>& fractions);

/// An internal function to extract the phase string, given the phase index;
/// Handy for printing the actual phase string in debug, warning, and error messages.
/// @param Phase The enumerated phase index to be looked up
std::string phase_lookup_string(phases Phase);

} /* namespace CoolProp */
#endif
