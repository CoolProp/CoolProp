/*
Add some pre-processor directives to this file so that it can either be built as 
usual, or if the COOLPROP_LIB macro is defined, it will export the functions in 
this file for building a static or dynamic library.  

The __stdcall calling convention is used by default.  By providing the macro CONVENTION, the 
calling convention can be changed at build time.

Any functions that are exported to DLL must also add the "EXPORT_CODE function_name CONVENTION ..." code 
to the CoolProp.cpp implementation.  See the Props function for instance
*/

/*! \mainpage CoolProp Core Code Documentation

Welcome to the home page for the C++ sources of CoolProp.  This information may be useful for developers or just the plain inquisitive

You might want to start by looking at CoolProp.h
*/

#ifndef CoolProp_H
#define CoolProp_H

    #include <string>
    #include <vector>

    namespace CoolProp {

    /// Return a value that depends on the thermodynamic state
	/// @param Output The output parameter, one of "T","D","H",etc.
	/// @param Name1 The first state variable name, one of "T","D","H",etc.
	/// @param Prop1 The first state variable value
	/// @param Name2 The second state variable name, one of "T","D","H",etc.
	/// @param Prop2 The second state variable value
	/// @param FluidName The fluid name
    double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName);
    /// Return a value that depends on the thermodynamic state
	/// @param Output The output parameter, one of "T","D","H",etc.
	/// @param Name1 The first state variable name, one of "T","D","H",etc.
	/// @param Prop1 The first state variable value
	/// @param Name2 The second state variable name, one of "T","D","H",etc.
	/// @param Prop2 The second state variable value
	/// @param FluidName The fluid name
    /// @param x The mole or mass fractions depending on the requirements of the backend
    double PropsSI(const std::string &Output, const std::string &Name1, double Prop1, const std::string &Name2, double Prop2, const std::string &FluidName, const std::vector<double> &x);
    /// Return a value that depends on the thermodynamic state
	/// @param Output The output parameter, one of "T","D","H",etc.
	/// @param Name1 The first state variable name, one of "T","D","H",etc.
	/// @param Prop1 The first state variable value
	/// @param Name2 The second state variable name, one of "T","D","H",etc.
	/// @param Prop2 The second state variable value
	/// @param FluidName The fluid name
    /// @param z The mole or mass fractions depending on the requirements of the backend
    std::vector<double> PropsSI(const std::string &Output, const std::string &Name1, const std::vector<double> &Prop1, const std::string &Name2, const std::vector<double> Prop2, const std::string &FluidName, const std::vector<double> &z);

    /**
    \overload 
    \sa PropsSI(std::string &Output, std::string &Name1, double Prop1, std::string &Name2, double Prop2, std::string &FluidName, const std::vector<double> &x);
    */
    double PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *FluidName, const std::vector<double> &x);

    /// Get the debug level
    /// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    int get_debug_level();
    /// Set the debug level
    /// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    void set_debug_level(int level);

    /// Set the global error string 
	/// @param error The error string to use
	void set_error_string(std::string error);
    /// An internal function to set the global warning string
	/// @param warning The string to set as the warning string
	void set_warning_string(std::string warning);

    /// Get a globally-defined string
	/// @param ParamName A string, one of "version", "errstring", "warnstring", "gitrevision", "FluidsList", "fluids_list"
	/// @returns str The string, or an error message if not valid input
	std::string get_global_param_string(std::string ParamName);

	/*/// Get a long that represents the fluid type
	/// @param FluidName The fluid name as a string
	/// @returns long element from global type enumeration
	long getFluidType(std::string FluidName);*/

	/// Get a string for a value from a fluid (numerical values can be obtained from Props1 function)
	/// @param FluidName The name of the fluid that is part of CoolProp, for instance "n-Propane"
	/// @param ParamName A string, one of "aliases", "CAS", "CAS_number", "ASHRAE34", "REFPROPName","REFPROP_name", 
	/// @returns str The string, or an error message if not valid input
	std::string get_fluid_param_string(std::string FluidName, std::string ParamName);

	/// Returns the BibTeX key from the bibtex library of CoolProp corresponding to the item requested
	/// @param FluidName The name of the fluid that is part of CoolProp, for instance "n-Propane"
	/// @param item The key that is desired, one of "EOS","CP0", "VISCOSITY", "CONDUCTIVITY", "ECS_LENNARD_JONES", "ECS_FITS", "SURFACE_TENSION"
	/// @returns key the BibTeX key
	std::string get_BibTeXKey(std::string FluidName, std::string item);
	
	/// Set the reference state based on a string representation of the desired reference state (consistent naming with REFPROP)
	/// @param FluidName The name of the fluid
	/// @param reference_state The reference state to use, one of "IIR" (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.) "ASHRAE" (h=0,s=0 @ -40C sat liq), "NBP" (h=0,s=0 @ 1.0 bar sat liq.)
	int set_reference_stateS(std::string FluidName, std::string reference_state);

	/// Set the reference state based on a thermodynamic state point
	/// @param FluidName The name of the fluid
	/// @param T Temperature at reference state [K]
	/// @param rho Density at reference state [mol/m^3]
	/// @param h0 Enthalpy at reference state [J/kg]
	/// @param s0 Entropy at references state [J/kg/K]
	int set_reference_stateD(std::string FluidName, double T, double rho, double h0, double s0);

    /*
    /// Return the phase of the given state point with temperature, pressure as inputs
	/// @param FluidName The name of the fluid
	/// @param T Temperature [K]
	/// @param p Pressure [kPa]
	/// @returns Phase as string, one of ""Two-Phase","Supercritical","Gas","Liquid"
	std::string Phase(std::string FluidName, double T, double p);
	/// Return the phase of the given state point with temperature, density as inputs
	/// @param FluidName The name of the fluid
	/// @param T Temperature [K]
	/// @param rho Density [kg/m^3]
	/// @returns Phase as string, one of ""Two-Phase","Supercritical","Gas","Liquid"
	std::string Phase_Trho(std::string FluidName, double T, double rho);
	/// Return the phase of the given state point with temperature, pressure as inputs
	/// @param FluidName The name of the fluid
	/// @param T Temperature [K]
	/// @param p Pressure [kPa]
	/// @returns Phase as string, one of ""Two-Phase","Supercritical","Gas","Liquid"
	std::string Phase_Tp(std::string FluidName, double T, double p);
    /// Return some low level derivative terms, see source for a complete list
	/// @param Term String, some options are "phir" (residual Helmholtz energy),"dphir_dDelta", "dphir_dTau", etc.
	/// @param T Temperature [K]
	/// @param rho Density [kg/m^3]
	/// @param FluidName String
	double DerivTerms(std::string Term, double T, double rho, std::string FluidName);
	/// Return some low level derivative terms, see source for a complete list
	/// @param iTerm long desired output
	/// @param T Temperature [K]
	/// @param rho Density [kg/m^3]
	/// @param pFluid Pointer to Fluid instance
	double DerivTerms(long iTerm, double T, double rho, Fluid * pFluid);*/

    } /* namespace CoolProp */
#endif

