/** \brief This file defines an interface for shared library (DLL) wrapping
 * 
 * In general the functions defined here take strings which are 0-terminated (C-style),  
 * vectors of doubles are passed as double* and length
 * These functions pass directly to equivalently named functions in CoolProp.h in the CoolProp namespace
 * that take std::string, vector<double> etc.
 * 
 * Functions with the call type like
 * EXPORT_CODE void CONVENTION AFunction(double, double);
 * will be exported to the DLL
 * 
 * The exact symbol that will be exported depends on the values of the preprocessor macros COOLPROP_LIB, EXPORT_CODE, CONVENTION, etc.
 * 
 * In order to have 100% control over the export macros, you can specify EXPORT_CODE and CONVENTION directly. Check out
 * CMakeLists.txt in the repo root to see some examples.
 * 
 */

#ifndef COOLPROPDLL_H
#define COOLPROPDLL_H

    #include "PlatformDetermination.h"

    #if defined(COOLPROP_LIB)
    #  ifndef EXPORT_CODE
    #    if defined(__ISWINDOWS__)
    #      define EXPORT_CODE extern "C" __declspec(dllexport)
    #    else
    #      define EXPORT_CODE extern "C"
    #    endif
    #  endif
    #  ifndef CONVENTION
    #    if defined(__ISWINDOWS__)
    #      define CONVENTION __stdcall
    #    else
    #      define CONVENTION
    #    endif
    #  endif
    #else
    #  ifndef EXPORT_CODE
    #    define EXPORT_CODE
    #  endif
    #  ifndef CONVENTION
    #    define CONVENTION
    #  endif
    #endif

    // Hack for PowerPC compilation to only use extern "C"
    #if defined(__powerpc__) || defined(EXTERNC)
    #  undef EXPORT_CODE
    #  define EXPORT_CODE extern "C"
    #endif
    
    #if defined(__powerpc__)
    // From https://rowley.zendesk.com/entries/46176--Undefined-reference-to-assert-error-message
    // The __assert function is an error handler function that is invoked when an assertion fails.
    // If you are writing a program that uses the assert macro then you must supply you own __assert error handler function. For example
    inline void __assert(const char *error)
    {
      while(1);
    }
    #endif

    /**
     * \overload
     * \sa \ref CoolProp::Props1SI(std::string, std::string)
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE double CONVENTION Props1SI(const char *FluidName, const char* Output);
    /**
     *\overload
     *\sa \ref CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE double CONVENTION PropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
    
    /**
     *\overload
     *\sa \ref CoolProp::PhaseSI(const std::string &, double, const std::string &, double, const std::string&)
     * 
     * \note This function returns the phase string in pre-allocated phase variable.  If buffer is not large enough, no copy is made
     */
    EXPORT_CODE long CONVENTION PhaseSI(const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref, char *phase, int n);
    
    /**
     *\overload
     *\sa \ref CoolProp::get_global_param_string
     * 
     * @returns error_code 1 = Ok 0 = error
     * 
     * \note This function returns the output string in pre-allocated char buffer.  If buffer is not large enough, no copy is made
     */
    EXPORT_CODE long CONVENTION get_global_param_string(const char *param, char *Output, int n);
    /**
     * \overload
     * \sa \ref CoolProp::get_parameter_information_string
     * \note This function returns the output string in pre-allocated char buffer.  If buffer is not large enough, no copy is made
     * 
     * @returns error_code 1 = Ok 0 = error
     */
    EXPORT_CODE long CONVENTION get_parameter_information_string(const char *key, char *Output, int n);
    /**
     * \overload
     * \sa \ref CoolProp::get_mixture_binary_pair_data
     */
    EXPORT_CODE long CONVENTION get_mixture_binary_pair_data(const char *CAS1, const char *CAS2, const char *key);
    /** 
     * \overload
     * \sa \ref CoolProp::get_fluid_param_string
     * 
     * @returns error_code 1 = Ok 0 = error
     */
    EXPORT_CODE long CONVENTION get_fluid_param_string(const char *fluid, const char *param, char *Output, int n);
    /**
     * \overload
     * \sa \ref CoolProp::set_reference_stateS
     * @returns error_code 1 = Ok 0 = error
     */
    EXPORT_CODE int CONVENTION set_reference_stateS(const char *Ref, const char *reference_state);
    /**
     * \overload
     * \sa \ref CoolProp::set_reference_stateD
     * @returns error_code 1 = Ok 0 = error
     */
    EXPORT_CODE int CONVENTION set_reference_stateD(const char *Ref, double T, double rho, double h0, double s0);
    /** \brief FORTRAN 77 style wrapper of the PropsSI function
     * \overload
     * \sa \ref CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE void CONVENTION propssi_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char * Ref, double *output);
    /** \brief An overload of \ref CoolProp::PropsSI with the mole fractions passed as an array of doubles and the number of 
     * \overload
     * \sa \ref CoolProp::PropsSIZ
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    //EXPORT_CODE double CONVENTION PropsSIZ(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *FluidName, const double *z, int n);

    /// Convert from degrees Fahrenheit to Kelvin (useful primarily for testing)
    EXPORT_CODE double CONVENTION F2K(double T_F);
    /// Convert from Kelvin to degrees Fahrenheit (useful primarily for testing)
    EXPORT_CODE double CONVENTION K2F(double T_K);
    /** \brief Get the index for a parameter "T", "P", etc.
     * 
     * @returns index The index as a long.  If input is invalid, returns -1
     */
    EXPORT_CODE long CONVENTION get_param_index(const char *param);
    /** \brief Redirect all output that would go to console (stdout) to a file
     */
    EXPORT_CODE long CONVENTION redirect_stdout(const char *file);

    // ---------------------------------
    // Getter and setter for debug level
    // ---------------------------------

    /// Get the debug level
    /// @returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    EXPORT_CODE int CONVENTION get_debug_level();
    /// Set the debug level
    /// @param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
    EXPORT_CODE void CONVENTION set_debug_level(int level);
    
    /* \brief Extract a value from the saturation ancillary
     * 
     * @param fluid_name The name of the fluid to be used - HelmholtzEOS backend only
     * @param output The desired output variable ("P" for instance for pressure)
     * @param Q The quality, 0 or 1
     * @param input The input variable ("T")
     * @param value The input value
     */
    EXPORT_CODE double CONVENTION saturation_ancillary(const char *fluid_name, const char *output, int Q, const char *input, double value);

    // ---------------------------------
    //        Humid Air Properties
    // ---------------------------------

    /** \brief DLL wrapper of the HAPropsSI function
     * \sa \ref HumidAir::HAPropsSI(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE double CONVENTION HAPropsSI(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Name3, double Prop3);

    /** \brief FORTRAN 77 style wrapper of the HAPropsSI function
     * \sa \ref HumidAir::HAPropsSI(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE void CONVENTION hapropssi_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char *Name3, const double *Prop3, double *output);
    
    /** \brief DLL wrapper of the HAProps function
     * 
     * \warning DEPRECATED!!
     * \sa \ref HumidAir::HAProps(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE double CONVENTION HAProps(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Name3, double Prop3);

    /** \brief FORTRAN 77 style wrapper of the HAProps function
     *
     * \warning DEPRECATED!!
     * \sa \ref HumidAir::HAProps(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     * 
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
    EXPORT_CODE void CONVENTION haprops_(const char *Output, const char *Name1, const double *Prop1, const char *Name2, const double *Prop2, const char *Name3, const double *Prop3, double *output);
    
    


    // *************************************************************************************
    // *************************************************************************************
    // *****************************  DEPRECATED *******************************************
    // *************************************************************************************
    // *************************************************************************************

    /**
    \overload
    \sa \ref Props(const char *Output, const char Name1, double Prop1, const char Name2, double Prop2, const char *Ref)
    */
    EXPORT_CODE double CONVENTION PropsS(const char *Output, const char *Name1, double Prop1, const char *Name2, double Prop2, const char *Ref);
    /**
    Works just like \ref CoolProp::PropsSI, but units are in KSI system.  This function is deprecated, no longer supported, and users should transition to using the PropsSI function
    */
    EXPORT_CODE double CONVENTION Props(const char *Output, const char Name1, double Prop1, const char Name2, double Prop2, const char *Ref);
    /**
    Works just like \ref CoolProp::Props1SI, but units are in KSI system.  This function is deprecated, no longer supported, and users should transition to using the Props1SI function
    */
    EXPORT_CODE double CONVENTION Props1(const char *FluidName, const char *Output);
    ///**
    //\overload
    // IsFluidType(std::string, std::string)
    //*/
    //EXPORT_CODE int CONVENTION IsFluidType(const char *Ref, const char *Type);

    // This version uses the indices in place of the strings for speed.  Get the parameter indices
    // from get_param_index('D') for instance and the Fluid index from get_Fluid_index('Air') for instance
    EXPORT_CODE double CONVENTION IPropsSI(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
    EXPORT_CODE double CONVENTION IProps(long iOutput, long iName1, double Prop1, long iName2, double Prop2, long iFluid);
#endif
