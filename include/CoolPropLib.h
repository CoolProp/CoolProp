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
#    define __ISWINDOWS__
#elif _WIN32
#    define __ISWINDOWS__
#elif __APPLE__
#    define __ISAPPLE__
#elif __linux || __unix || __posix
#    define __ISLINUX__
#elif __powerpc__
#    define __ISPOWERPC__
#else
#    pragma error
#endif

#if defined(COOLPROP_LIB)
#    ifndef EXPORT_CODE
#        if defined(__ISWINDOWS__)
#            define EXPORT_CODE extern "C" __declspec(dllexport)
#        else
#            define EXPORT_CODE extern "C"
#        endif
#    endif
#    ifndef CONVENTION
#        if defined(__ISWINDOWS__)
#            define CONVENTION __stdcall
#        else
#            define CONVENTION
#        endif
#    endif
#else
#    ifndef EXPORT_CODE
#        define EXPORT_CODE
#    endif
#    ifndef CONVENTION
#        define CONVENTION
#    endif
#endif

// Hack for PowerPC compilation to only use extern "C"
#if defined(__powerpc__) || defined(EXTERNC)
#    undef EXPORT_CODE
#    define EXPORT_CODE extern "C"
#endif

#if defined(__powerpc__)
// From https://rowley.zendesk.com/entries/46176--Undefined-reference-to-assert-error-message
// The __assert function is an error handler function that is invoked when an assertion fails.
// If you are writing a program that uses the assert macro then you must supply you own __assert error handler function. For example
inline void __assert(const char* error) {
    while (1)
        ;
}
#endif

/**
     * \overload
     * \sa \ref CoolProp::Props1SI(std::string, std::string)
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE double CONVENTION Props1SI(const char* FluidName, const char* Output);

/**
     *\overload
     *\sa \ref CoolProp::Props1SImulti(const std::vector<std::string>& Outputs, const std::string& backend, const std::vector<std::string>& fluids, const std::vector<double>& fractions)
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION Props1SImulti(const char* Outputs, char* backend, const char* FluidNames, const double* fractions,
                                          const long length_fractions, double* result, long* resdim1);
/**
     *\overload
     *\sa \ref CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE double CONVENTION PropsSI(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Ref);
/**
     *\overload
     *\sa \ref CoolProp::PropsSImulti(const std::vector<std::string>& Outputs, const std::string& Name1, const std::vector<double>& Prop1,
                                              const std::string& Name2, const std::vector<double>& Prop2, const std::string& backend,
                                              const std::vector<std::string>& fluids, const std::vector<double>& fractions)
     *
     * @param Outputs Delimited string separated by LIST_STRING_DELIMITER for the output parameters
     * @param Name1 The name of the first input variable
     * @param Prop1 A vector of the first input values 
     * @param size_Prop1 Size of Prop1 double*
     * @param Name2 The name of the second input variable
     * @param Prop2 A vector of the second input values 
     * @param size_Prop2 Size of Prop2 double*
     * @param backend 	The string representation of the backend (HEOS, REFPROP, INCOMP, etc.) 
     * @param FluidNames  Delimited string separated by LIST_STRING_DELIMITER for the fluid name(s)
     * @param fractions The fractions (molar, mass, volume, etc.) of the components
     * @param length_fractions Size of fractions double*
     * @param result Allocated memory for result vector
     * @param resdim1 result vector dimension 1 pointer, to check allocated space and return actual result size
     * @param resdim2 result vector dimension 2 pointer, to check allocated space and return actual result size
     * \note If there is an error, an empty vector will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION PropsSImulti(const char* Outputs, const char* Name1, double* Prop1, const long size_Prop1, const char* Name2,
                                         double* Prop2, const long size_Prop2, char* backend, const char* FluidNames, const double* fractions,
                                         const long length_fractions, double* result, long* resdim1, long* resdim2);
/**
     *\overload
     *\sa \ref CoolProp::PhaseSI(const std::string &, double, const std::string &, double, const std::string&)
     *
     * \note This function returns the phase string in pre-allocated phase variable.  If buffer is not large enough, no copy is made
     */
EXPORT_CODE long CONVENTION PhaseSI(const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Ref, char* phase, int n);

/**
     *\overload
     *\sa \ref CoolProp::get_global_param_string
     *
     * @returns error_code 1 = Ok 0 = error
     *
     * \note This function returns the output string in pre-allocated char buffer.  If buffer is not large enough, no copy is made
     */
EXPORT_CODE long CONVENTION get_global_param_string(const char* param, char* Output, int n);
/**
     * \overload
     * \sa \ref CoolProp::get_parameter_information_string
     * \note This function returns the output string in pre-allocated char buffer.  If buffer is not large enough, no copy is made
     *
     * @returns error_code 1 = Ok 0 = error
     */
EXPORT_CODE long CONVENTION get_parameter_information_string(const char* key, char* Output, int n);
/**
     * \overload
     * \sa \ref CoolProp::get_fluid_param_string
     *
     * @returns error_code 1 = Ok 0 = error
     */
EXPORT_CODE long CONVENTION get_fluid_param_string(const char* fluid, const char* param, char* Output, int n);
/** \brief Set configuration string
    * @param key The key to configure
    * @param val The value to set to the key
    * \note you can get the error message by doing something like get_global_param_string("errstring",output)
    */
EXPORT_CODE void CONVENTION set_config_string(const char* key, const char* val);
/** \brief Set configuration numerical value as double
    * @param key The key to configure
    * @param val The value to set to the key
    * \note you can get the error message by doing something like get_global_param_string("errstring",output)
    */
EXPORT_CODE void CONVENTION set_config_double(const char* key, const double val);
/** \brief Set configuration value as a boolean
     * @param key The key to configure
     * @param val The value to set to the key
     * \note you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION set_config_bool(const char* key, const bool val);
/**
     * @brief Set the departure functions in the departure function library from a string format
     * @param string_data The departure functions to be set, either provided as a JSON-formatted string
     *                    or as a string of the contents of a HMX.BNC file from REFPROP
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     *
     * @note By default, if a departure function already exists in the library, this is an error,
     *       unless the configuration variable OVERWRITE_DEPARTURE_FUNCTIONS is set to true
     */
EXPORT_CODE void CONVENTION set_departure_functions(const char* string_data, long* errcode, char* message_buffer, const long buffer_length);
/**
     * \overload
     * \sa \ref CoolProp::set_reference_stateS
     * @returns error_code 1 = Ok 0 = error
     */
EXPORT_CODE int CONVENTION set_reference_stateS(const char* Ref, const char* reference_state);
/**
     * \overload
     * \sa \ref CoolProp::set_reference_stateD
     * @returns error_code 1 = Ok 0 = error
     */
EXPORT_CODE int CONVENTION set_reference_stateD(const char* Ref, double T, double rhomolar, double hmolar0, double smolar0);
/** \brief FORTRAN 77 style wrapper of the PropsSI function
     * \overload
     * \sa \ref CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION propssi_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                     const char* Ref, double* output);

/// Convert from degrees Fahrenheit to Kelvin (useful primarily for testing)
EXPORT_CODE double CONVENTION F2K(double T_F);
/// Convert from Kelvin to degrees Fahrenheit (useful primarily for testing)
EXPORT_CODE double CONVENTION K2F(double T_K);
/** \brief Get the index for a parameter "T", "P", etc.
     *
     * @returns index The index as a long.  If input is invalid, returns -1
     */
EXPORT_CODE long CONVENTION get_param_index(const char* param);
/** \brief Get the index for an input pair for AbstractState.update function
     *
     * @returns index The index as a long.  If input is invalid, returns -1
     */
EXPORT_CODE long CONVENTION get_input_pair_index(const char* param);
/** \brief Redirect all output that would go to console (stdout) to a file
     */
EXPORT_CODE long CONVENTION redirect_stdout(const char* file);

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
EXPORT_CODE double CONVENTION saturation_ancillary(const char* fluid_name, const char* output, int Q, const char* input, double value);

// ---------------------------------
//        Humid Air Properties
// ---------------------------------

/** \brief DLL wrapper of the HAPropsSI function
     * \sa \ref HumidAir::HAPropsSI(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE double CONVENTION HAPropsSI(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Name3,
                                        double Prop3);

/** \brief Humid air saturation specific heat at 1 atmosphere, based on a correlation from EES.
     * \sa \ref HumidAir::cair_sat(double);
     *
     * @param T [K] good from 250K to 300K, no error bound checking is carried out.
     *
     * \note Equals partial derivative of enthalpy with respect to temperature at constant relative humidity of 100 percent and pressure of 1 atmosphere.
     */
EXPORT_CODE double CONVENTION cair_sat(double T);

/** \brief FORTRAN 77 style wrapper of the HAPropsSI function
     * \sa \ref HumidAir::HAPropsSI(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION hapropssi_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                       const char* Name3, const double* Prop3, double* output);

/** \brief DLL wrapper of the HAProps function
     *
     * \warning DEPRECATED!!
     * \sa \ref HumidAir::HAProps(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE double CONVENTION HAProps(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Name3,
                                      double Prop3);

/** \brief FORTRAN 77 style wrapper of the HAProps function
     *
     * \warning DEPRECATED!!
     * \sa \ref HumidAir::HAProps(const char *OutputName, const char *Input1Name, double Input1, const char *Input2Name, double Input2, const char *Input3Name, double Input3);
     *
     * \note If there is an error, a huge value will be returned, you can get the error message by doing something like get_global_param_string("errstring",output)
     */
EXPORT_CODE void CONVENTION haprops_(const char* Output, const char* Name1, const double* Prop1, const char* Name2, const double* Prop2,
                                     const char* Name3, const double* Prop3, double* output);

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
EXPORT_CODE long CONVENTION AbstractState_factory(const char* backend, const char* fluids, long* errcode, char* message_buffer,
                                                  const long buffer_length);
/**
     * @brief Get the fluid names for the AbstractState
     * @param handle The integer handle for the state class stored in memory
     * @param fluids LIST_STRING_DELIMETER (',') delimited list of fluids
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE void CONVENTION AbstractState_fluid_names(const long handle, char* fluids, long* errcode, char* message_buffer, const long buffer_length);
/**
     * @brief Release a state class generated by the low-level interface wrapper
     * @param handle The integer handle for the state class stored in memory
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE void CONVENTION AbstractState_free(const long handle, long* errcode, char* message_buffer, const long buffer_length);
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
EXPORT_CODE void CONVENTION AbstractState_set_fractions(const long handle, const double* fractions, const long N, long* errcode, char* message_buffer,
                                                        const long buffer_length);
/**
     * @brief Get the molar fractions for the AbstractState
     * @param handle The integer handle for the state class stored in memory
     * @param fractions The array of fractions
     * @param maxN The length of the buffer for the fractions
     * @param N number of fluids
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE void CONVENTION AbstractState_get_mole_fractions(const long handle, double* fractions, const long maxN, long* N, long* errcode,
                                                             char* message_buffer, const long buffer_length);
/**
     * @brief Get the molar fractions for the AbstractState and the desired saturated State
     * @param handle The integer handle for the state class stored in memory
     * @param saturated_state The string specifying the state (liquid or gas)
     * @param fractions The array of fractions
     * @param maxN The length of the buffer for the fractions
     * @param N number of fluids
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return 
     */
EXPORT_CODE void CONVENTION AbstractState_get_mole_fractions_satState(const long handle, const char* saturated_state, double* fractions,
                                                                      const long maxN, long* N, long* errcode, char* message_buffer,
                                                                      const long buffer_length);
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
EXPORT_CODE void CONVENTION AbstractState_update(const long handle, const long input_pair, const double value1, const double value2, long* errcode,
                                                 char* message_buffer, const long buffer_length);
/**
    * @brief Specify the phase to be used for all further calculations
    * @param handle The integer handle for the state class stored in memory
    * @param phase The string with the phase to use
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
EXPORT_CODE void CONVENTION AbstractState_specify_phase(const long handle, const char* phase, long* errcode, char* message_buffer,
                                                        const long buffer_length);
/**
    * @brief Unspecify the phase to be used for all further calculations
    * @param handle The integer handle for the state class stored in memory
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
EXPORT_CODE void CONVENTION AbstractState_unspecify_phase(const long handle, long* errcode, char* message_buffer, const long buffer_length);
/**
     * @brief Get an output value from the AbstractState using an integer value for the desired output value
     * @param handle The integer handle for the state class stored in memory
     * @param param The integer value for the parameter you want
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE double CONVENTION AbstractState_keyed_output(const long handle, const long param, long* errcode, char* message_buffer,
                                                         const long buffer_length);

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
EXPORT_CODE double CONVENTION AbstractState_first_saturation_deriv(const long handle, const long Of, const long Wrt, long* errcode,
                                                                   char* message_buffer, const long buffer_length);

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
EXPORT_CODE double CONVENTION AbstractState_first_partial_deriv(const long handle, const long Of, const long Wrt, const long Constant, long* errcode,
                                                                char* message_buffer, const long buffer_length);

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
EXPORT_CODE void CONVENTION AbstractState_update_and_common_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                                const long length, double* T, double* p, double* rhomolar, double* hmolar,
                                                                double* smolar, long* errcode, char* message_buffer, const long buffer_length);

/**
    * @brief Update the state of the AbstractState and get one output value (temperature, pressure, molar density, molar enthalpy and molar entropy)
    * @brief from the AbstractState using pointers as inputs and output to allow array computation.
    * @param handle The integer handle for the state class stored in memory
    * @param input_pair The integer value for the input pair obtained from get_input_pair_index
    * @param value1 The pointer to the array of the first input parameters
    * @param value2 The pointer to the array of the second input parameters
    * @param length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
    * @param output The indice for the output desired
    * @param out The pointer to the array for output
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    *
    * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
    */
EXPORT_CODE void CONVENTION AbstractState_update_and_1_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                           const long length, const long output, double* out, long* errcode, char* message_buffer,
                                                           const long buffer_length);

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
EXPORT_CODE void CONVENTION AbstractState_update_and_5_out(const long handle, const long input_pair, const double* value1, const double* value2,
                                                           const long length, long* outputs, double* out1, double* out2, double* out3, double* out4,
                                                           double* out5, long* errcode, char* message_buffer, const long buffer_length);

/**
    * @brief Set binary interraction parrameter for mixtures
    * @param handle The integer handle for the state class stored in memory
    * @param i indice of the first fluid of the binary pair
    * @param j indice of the second fluid of the binary pair
    * @param parameter string wit the name of the parameter
    * @param value the value of the binary interaction parameter
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
EXPORT_CODE void CONVENTION AbstractState_set_binary_interaction_double(const long handle, const long i, const long j, const char* parameter,
                                                                        const double value, long* errcode, char* message_buffer,
                                                                        const long buffer_length);

/**
    * @brief Set cubic's alpha function parameters
    * @param handle The integer handle for the state class stored in memory
    * @param i indice of the fluid the parramter should be applied too (for mixtures)
	* @param parameter the string specifying the alpha function to use, ex "TWU" for the TWU alpha function
    * @param c1 the first parameter for the alpha function
    * @param c2 the second parameter for the alpha function
    * @param c3 the third parameter for the alpha function
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
EXPORT_CODE void CONVENTION AbstractState_set_cubic_alpha_C(const long handle, const long i, const char* parameter, const double c1, const double c2,
                                                            const double c3, long* errcode, char* message_buffer, const long buffer_length);

/**
    * @brief Set some fluid parameter (ie volume translation for cubic)
    * @param handle The integer handle for the state class stored in memory
	* @param i indice of the fluid the parramter should be applied too (for mixtures)
	* @param parameter the string specifying the parameter to use, ex "cm" for volume translation
    * @param value the value of the parameter
    * @param errcode The errorcode that is returned (0 = no error, !0 = error)
    * @param message_buffer A buffer for the error code
    * @param buffer_length The length of the buffer for the error code
    * @return
    */
EXPORT_CODE void CONVENTION AbstractState_set_fluid_parameter_double(const long handle, const long i, const char* parameter, const double value,
                                                                     long* errcode, char* message_buffer, const long buffer_length);

/**
     * @brief Build the phase envelope
     * @param handle The integer handle for the state class stored in memory
     * @param level How much refining of the phase envelope ("none" to skip refining (recommended))
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     *
     * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
     */
EXPORT_CODE void CONVENTION AbstractState_build_phase_envelope(const long handle, const char* level, long* errcode, char* message_buffer,
                                                               const long buffer_length);

/**
     * @brief Get data from the phase envelope for the given mixture composition
     * @param handle The integer handle for the state class stored in memory
     * @param length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
     * @param T The pointer to the array of temperature (K)
     * @param p The pointer to the array of pressure (Pa)
     * @param rhomolar_vap The pointer to the array of molar density for vapor phase (m^3/mol)
     * @param rhomolar_liq The pointer to the array of molar density for liquid phase (m^3/mol)
     * @param x The compositions of the "liquid" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)
     * @param y The compositions of the "vapor" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     *
     * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
     */
EXPORT_CODE void CONVENTION AbstractState_get_phase_envelope_data(const long handle, const long length, double* T, double* p, double* rhomolar_vap,
                                                                  double* rhomolar_liq, double* x, double* y, long* errcode, char* message_buffer,
                                                                  const long buffer_length);

/**
     * @brief Get data from the phase envelope for the given mixture composition
     * @param handle The integer handle for the state class stored in memory
     * @param length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
     * @param maxComponents The number of fluid components for which memory is allocated
     * @param T The pointer to the array of temperature (K)
     * @param p The pointer to the array of pressure (Pa)
     * @param rhomolar_vap The pointer to the array of molar density for vapor phase (m^3/mol)
     * @param rhomolar_liq The pointer to the array of molar density for liquid phase (m^3/mol)
     * @param x The compositions of the "liquid" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)
     * @param y The compositions of the "vapor" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)
     * @param actual_length The number of elements actually stored in the arrays
     * @param actual_components The number of fluid components actually stored in the arrays
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     *
     * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
     */
EXPORT_CODE void CONVENTION AbstractState_get_phase_envelope_data_checkedMemory(const long handle, const long length, const long maxComponents, double* T,
                                                                  double* p, double* rhomolar_vap, double* rhomolar_liq, double* x, double* y,
                                                                  long* actual_length, long* actual_components, long* errcode, char* message_buffer,
                                                                  const long buffer_length);

/**
     * @brief Build the spinodal
     * @param handle The integer handle for the state class stored in memory
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE void CONVENTION AbstractState_build_spinodal(const long handle, long* errcode, char* message_buffer, const long buffer_length);

/**
     * @brief Get data for the spinodal curve
     * @param handle The integer handle for the state class stored in memory
     * @param length The number of elements stored in the arrays (all outputs MUST be the same length)
     * @param tau The pointer to the array of reciprocal reduced temperature
     * @param delta The pointer to the array of reduced density
     * @param M1 The pointer to the array of M1 values (when L1=M1=0, critical point)
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     *
     * @note If there is an error, no change in the output arrays will be made
     */
EXPORT_CODE void CONVENTION AbstractState_get_spinodal_data(const long handle, const long length, double* tau, double* delta, double* M1,
                                                            long* errcode, char* message_buffer, const long buffer_length);

/**
     * @brief Calculate all the critical points for a given composition
     * @param handle The integer handle for the state class stored in memory
     * @param length The length of the buffers passed to this function
     * @param T The pointer to the array of temperature (K)
     * @param p The pointer to the array of pressure (Pa)
     * @param rhomolar The pointer to the array of molar density (m^3/mol)
     * @param stable The pointer to the array of boolean flags for whether the critical point is stable (1) or unstable (0)
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     *
     * @note If there is an error in an update call for one of the inputs, no change in the output array will be made
     */
EXPORT_CODE void CONVENTION AbstractState_all_critical_points(const long handle, const long length, double* T, double* p, double* rhomolar,
                                                              long* stable, long* errcode, char* message_buffer, const long buffer_length);
/**
     * @brief Get an output value from the AbstractState using an integer value for the desired output value and desired saturated State
     * @param handle The integer handle for the state class stored in memory
     * @param saturated_state The string specifying the state (liquid or gas)
     * @param param The integer value for the parameter you want
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE double CONVENTION AbstractState_keyed_output_satState(const long handle, const char* saturated_state, const long param, long* errcode,
                                                                  char* message_buffer, const long buffer_length);
/**
     * @brief Return the name of the backend used in the AbstractState
     * @param handle The integer handle for the state class stored in memory
     * @param backend The char pointer the name is written to
     * @param errcode The errorcode that is returned (0 = no error, !0 = error)
     * @param message_buffer A buffer for the error code
     * @param buffer_length The length of the buffer for the error code
     * @return
     */
EXPORT_CODE void CONVENTION AbstractState_backend_name(const long handle, char* backend, long* errcode, char* message_buffer,
                                                       const long buffer_length);
/** 
     * \brief Add fluids as a JSON-formatted string
     * @param backend The backend to which these should be added; e.g. "HEOS", "SRK", "PR"
     * @param fluidstring The JSON-formatted string
     * @return
     *
     */
EXPORT_CODE void CONVENTION add_fluids_as_JSON(const char* backend, const char* fluidstring, long* errcode, char* message_buffer,
                                               const long buffer_length);

// *************************************************************************************
// *************************************************************************************
// *****************************  DEPRECATED *******************************************
// *************************************************************************************
// *************************************************************************************

/**
    \overload
    \sa \ref Props(const char *Output, const char Name1, double Prop1, const char Name2, double Prop2, const char *Ref)
    */
EXPORT_CODE double CONVENTION PropsS(const char* Output, const char* Name1, double Prop1, const char* Name2, double Prop2, const char* Ref);
/**
    Works just like \ref CoolProp::PropsSI, but units are in KSI system.  This function is deprecated, no longer supported, and users should transition to using the PropsSI function
    */
EXPORT_CODE double CONVENTION Props(const char* Output, const char Name1, double Prop1, const char Name2, double Prop2, const char* Ref);
/**
    Works just like \ref CoolProp::Props1SI, but units are in KSI system.  This function is deprecated, no longer supported, and users should transition to using the Props1SI function
    */
EXPORT_CODE double CONVENTION Props1(const char* FluidName, const char* Output);

#endif
