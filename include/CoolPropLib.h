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

    // See also http://stackoverflow.com/questions/5919996/how-to-detect-reliably-mac-os-x-ios-linux-windows-in-c-preprocessor
    // Copied verbatim from PlatformDetermination.h in order to have a single-include header
    #if _WIN64
    #  define __ISWINDOWS__
    #elif _WIN32
    #  define __ISWINDOWS__
    #elif __APPLE__
    #  define __ISAPPLE__
    #elif __linux || __unix || __posix
    #  define __ISLINUX__
    #elif __powerpc__
    #  define __ISPOWERPC__
    #else
    # pragma error
    #endif

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

    /// Convert from degrees Fahrenheit to Kelvin (useful primarily for testing)
    EXPORT_CODE double CONVENTION F2K(double T_F);
    /// Convert from Kelvin to degrees Fahrenheit (useful primarily for testing)
    EXPORT_CODE double CONVENTION K2F(double T_K);
    /** \brief Get the index for a parameter "T", "P", etc.
     * 
     * @returns index The index as a long.  If input is invalid, returns -1
     */
    EXPORT_CODE long CONVENTION get_param_index(const char *param);
    /** \brief Get the index for an input pair for AbstractState.update function
     * 
     * @returns index The index as a long.  If input is invalid, returns -1
     */
    EXPORT_CODE long CONVENTION get_input_pair_index(const char *param);
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

    /** \brief DLL wrapper of the cair_sat function
     * \sa \ref HumidAir::cair_sat(double);
     */
    EXPORT_CODE double CONVENTION cair_sat(double T);

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
    
    // ---------------------------------
    //        Low-level access
    // ---------------------------------
    
    /**
     * @brief Generate an AbstractState instance, return an integer handle to the state class generated to be used in the other low-level accessor functions
     * @param backend The backend you will use, "HEOS", "REFPROP", etc.
     * @param fluids '&' delimited list of fluids
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return A handle to the state class generated
     */
    EXPORT_CODE long CONVENTION AbstractState_factory(const char* backend, const char* fluids, long *errcode, char *message_buffer, const long buffer_length);
    /**
     * @brief Release a state class generated by the low-level interface wrapper
     * @param handle The integer handle for the state class stored in memory
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return 
     */
    EXPORT_CODE void CONVENTION AbstractState_free(const long handle, long *errcode, char *message_buffer, const long buffer_length);
    /**
     * @brief Set the fractions (mole, mass, volume) for the AbstractState
     * @param handle The integer handle for the state class stored in memory
     * @param fractions The array of fractions
     * @param N The length of the fractions array
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return 
     */
    EXPORT_CODE void CONVENTION AbstractState_set_fractions(const long handle, const double* fractions, const long N, long *errcode, char *message_buffer, const long buffer_length);
    /**
     * @brief Update the state of the AbstractState
     * @param handle The integer handle for the state class stored in memory
     * @param input_pair The integer value for the input pair obtained from XXXXXXXXXXXXXXXX
     * @param value1 The first input value
     * @param value2 The second input value
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return 
     */
    EXPORT_CODE void CONVENTION AbstractState_update(const long handle, const long input_pair, const double value1, const double value2, long *errcode, char *message_buffer, const long buffer_length);
    /**
     * @brief Get an output value from the AbstractState using an integer value for the desired output value
     * @param handle The integer handle for the state class stored in memory
     * @param param The integer value for the parameter you want
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return 
     */
    EXPORT_CODE double CONVENTION AbstractState_keyed_output(const long handle, const long param, long *errcode, char *message_buffer, const long buffer_length);

    /**
    * @brief Calculate a saturation derivative from the AbstractState using integer values for the desired parameters
    * @param handle The integer handle for the state class stored in memory
    * @param Of The parameter of which the derivative is being taken
    * @param Wrt The derivative with with respect to this parameter
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
    EXPORT_CODE double CONVENTION AbstractState_first_saturation_deriv(const long handle, const long Of, const long Wrt, long *errcode, char *message_buffer, const long buffer_length);
    
    /**
    * @brief Calculate the first partial derivative in homogeneous phases from the AbstractState using integer values for the desired parameters
    * @param handle The integer handle for the state class stored in memory
    * @param Of The parameter of which the derivative is being taken
    * @param Wrt The derivative with with respect to this parameter
    * @param Constant The parameter that is not affected by the derivative
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
    EXPORT_CODE double CONVENTION AbstractState_first_partial_deriv(const long handle, const long Of, const long Wrt, const long Constant, long *errcode, char *message_buffer, const long buffer_length);
    
    /**
    * @brief Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy)
    * @brief from the AbstractState using pointers as inputs and output to allow array computation.
    * @param handle The integer handle for the state class stored in memory
    * @param input_pair The integer value for the input pair obtained from get_input_pair_index
    * @param value1 The pointer to the array of the first input parameters
    * @param value2 The pointer to the array of the second input parameters
    * @param length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
    * @param T The pointer to the array of temperature
    * @param p The pointer to the array of pressure
    * @param rhomolar The pointer to the array of molar density
    * @param hmolar The pointer to the array of molar enthalpy
    * @param smolar The pointer to the array of molar entropy
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    *
    * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
    */
    EXPORT_CODE void CONVENTION AbstractState_update_and_common_out(const long handle, const long input_pair, const double* value1, const double* value2, const long length, double* T, double* p, double* rhomolar, double* hmolar, double* smolar, long *errcode, char *message_buffer, const long buffer_length);

    /**
    * @brief Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy)
    * @brief from the AbstractState using pointers as inputs and output to allow array computation.
    * @param handle The integer handle for the state class stored in memory
    * @param input_pair The integer value for the input pair obtained from get_input_pair_index
    * @param value1 The pointer to the array of the first input parameters
    * @param value2 The pointer to the array of the second input parameters
    * @param length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
    * @param outputs The 5-element vector of indices for the outputs desired
    * @param out1 The pointer to the array for the first output
    * @param out2 The pointer to the array for the second output
    * @param out3 The pointer to the array for the third output
    * @param out4 The pointer to the array for the fourth output
    * @param out5 The pointer to the array for the fifth output
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    *
    * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
    */
    EXPORT_CODE void CONVENTION AbstractState_update_and_5_out(const long handle, const long input_pair, const double* value1, const double* value2, const long length, long *outputs, double* out1, double* out2, double* out3, double* out4, double* out5, long *errcode, char *message_buffer, const long buffer_length);


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
