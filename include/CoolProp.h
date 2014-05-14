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

    /*/// Return a fluid value that does not depend on the thermodynamic state
	/// @param FluidName The name of the fluid
	/// @param Output The name of the output parameter, some options are "Ttriple", "Tcrit", "pcrit", "Tmin", "molemass", "rhocrit", "accentric" (not all parameters are valid for all fluids)
	/// @returns val The value, or _HUGE if not valid
	double Props1SI(std::string FluidName,std::string Output);*/
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
	long getFluidType(std::string FluidName);
	/// Get a string for a value from a fluid (numerical values can be obtained from Props1 function)
	/// @param FluidName The name of the fluid
	/// @param ParamName A string, one of "aliases", "CAS", "CAS_number", "ASHRAE34", "REFPROPName","REFPROP_name", "TTSE_mode"		
	/// @returns str The string, or an error message if not valid input
	std::string get_fluid_param_string(std::string FluidName, std::string ParamName);
	/// Returns the BibTeX key from the bibtex library of CoolProp corresponding to the item requested
	/// @param FluidName The name of the fluid
	/// @param item String, one of "EOS","CP0", "VISCOSITY", "CONDUCTIVITY", "ECS_LENNARD_JONES", "ECS_FITS", "SURFACE_TENSION"
	/// @returns key the BibTeX key
	std::string get_BibTeXKey(std::string FluidName, std::string item);
	/// Add a REFPROP fluid to the fluid listing in CoolProp
	bool add_REFPROP_fluid(std::string FluidName);
	/// Get the parameter index for a given parameter
	/// @param param The string for an input or output parameter
	/// @returns index The parameter index (for use in IProps or elsewhere); -1 if not found
	long get_param_index(std::string param);
	/// Set the reference state for a pointer to a fluid (not exposed)
	/// @param pFluid A pointer to a Fluid instance
	/// @param reference_state The reference state to use, one of "IIR" (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.) "ASHRAE" (h=0,s=0 @ -40C sat liq), "NBP" (h=0,s=0 @ 1.0 bar sat liq.)
	//int set_reference_stateP(Fluid *pFluid, std::string reference_state);
	/// Set the reference state based on a string representation (consistent naming with RFPROP)
	/// @param FluidName The name of the fluid
	/// @param reference_state The reference state to use, one of "IIR" (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.) "ASHRAE" (h=0,s=0 @ -40C sat liq), "NBP" (h=0,s=0 @ 1.0 bar sat liq.)
	//int set_reference_stateS(std::string FluidName, std::string reference_state);*/
	/// Set the reference state based on a state point
	/// @param FluidName The name of the fluid
	/// @param T Temperature at reference state [K]
	/// @param rho Density at reference state [kg/m^3]
	/// @param h0 Enthalpy at reference state [kJ/kg]
	/// @param s0 Entropy at references state [kJ/kg/K]
	//int set_reference_stateD(std::string FluidName, double T, double rho, double h0, double s0);

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
    /*/// Return a fluid value that does not depend on the thermodynamic state
	/// @param FluidName The name of the fluid
	/// @param Output The name of the output parameter, some options are "Ttriple", "Tcrit", "pcrit", "Tmin", "molemass", "rhocrit", "accentric" (not all parameters are valid for all fluids)
	/// @returns val The value, or _HUGE if not valid
	double Props1(std::string FluidName,std::string Output);
    /// Return a value that depends on the thermodynamic state
	/// @param Output The output parameter, one of "T","D","H",etc.
	/// @param Name1 The first state variable name, one of "T","D","H",etc.
	/// @param Prop1 The first state variable value
	/// @param Name2 The second state variable name, one of "T","D","H",etc.
	/// @param Prop2 The second state variable value
	/// @param FluidName The fluid name
    double Props(std::string Output,std::string Name1, double Prop1, std::string Name2, double Prop2, std::string FluidName);
	/// Return a value that depends on the thermodynamic state
	/// @param Output The output parameter, one of "T","D","H",etc.
	/// @param Name1 The first state variable name, one of "T","D","H",etc.
	/// @param Prop1 The first state variable value
	/// @param Name2 The second state variable name, one of "T","D","H",etc.
	/// @param Prop2 The second state variable value
	/// @param FluidName The fluid name
	double Props(std::string Output,char Name1, double Prop1, char Name2, double Prop2, std::string FluidName);*/

	

	////    **************** DEPRECATION WARNING ***********************
	///// Nearly deprecated function
	//void set_phase(std::string Phase_str);

    } /* namespace CoolProp */
#endif

