// CoolPropMathcad.cpp : Defines the exported functions for the DLL Add-in.
//

#include <string>
#include <cstring>

#ifndef NOMINMAX  // Kill windows' horrible min() and max() macros
#    define NOMINMAX
#endif
#include "mcadincl.h"
#undef NOMINMAX

enum
{
    MC_STRING = STRING
};             // substitute enumeration variable MC_STRING for STRING, use MC_STRING below
#undef STRING  // undefine STRING as it conflicts with STRING enum in fmtlib/format.h

#include "CoolProp.h"
#include "DataStructures.h"
#include "HumidAirProp.h"
#include <Backends/Helmholtz/MixtureParameters.h>

// Setup Dialog Window for debugging
HWND hwndDlg;  // Generic Dialog handle for pop-up message boxes (MessageBox) when needed

namespace CoolProp {
extern void apply_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule);
}

enum EC
{
    MUST_BE_REAL = 1,  // Mathcad Error Codes       v
    INSUFFICIENT_MEMORY,
    INTERRUPTED,
    ONLY_ONE_COLUMN,
    //---------------------------------------------------
    BAD_FLUID,        // CoolProp Error Codes from here v
    BAD_MIXTURE,
    BAD_IF97_FLUID,
    BAD_PARAMETER,
    BAD_PHASE,
    ONLY_ONE_PHASE_SPEC,
    BAD_REF,
    NON_TRIVIAL,
    NO_REFPROP,
    NOT_AVAIL,
    BAD_INPUT_PAIR,
    BAD_QUAL,
    TWO_PHASE,
    NON_TWO_PHASE,
    T_OUT_OF_RANGE,
    P_OUT_OF_RANGE,
    H_OUT_OF_RANGE,
    S_OUT_OF_RANGE,
    TP_SATURATION,
    HA_INPUTS,
    BAD_BINARY_PAIR,
    MISSING_BINARY_PAIR,
    BAD_RULE,
    PAIR_EXISTS,
    NO_SOLUTION,
    UNKNOWN,
    NUMBER_OF_ERRORS
};  // Dummy Code for Error Count

// table of error messages
// As of Mathcad Prime 10, these are now actually returned as Custom Error: messages
char* CPErrorMessageTable[NUMBER_OF_ERRORS] = {"Argument must be real",
                                               "Insufficient Memory",
                                               "Interrupted",
                                               "Only one column allowed in input array",
                                               "Invalid Fluid String",
                                               "Invalid predefined mixture",
                                               "IF97 Backend supports pure \"Water\" only",
                                               "Invalid Parameter String",
                                               "Invalid Phase String",
                                               "Only one input key phase specification allowed",
                                               "Cannot use this REF State with this fluid",
                                               "Input Parameter is Non-Trivial",
                                               "REFPROP not installed correctly",
                                               "This Output parameter is not available for this Fluid",
                                               "This Input Pair is not yet support for this Fluid",
                                               "Input vapor quality must be between 0 and 1",
                                               "Output variable not valid in two phase region",
                                               "Output variable only valid in two phase region",
                                               "Temperature out of range",
                                               "Pressure out of range",
                                               "Enthalpy out of range",
                                               "Entropy out of range",
                                               "Temperature-Pressure inputs in 2-phase region; use TQ or PQ",
                                               "At least one of the inputs must be [T], [R], [W], or [Tdp]",
                                               "Could not match binary pair",
                                               "Missing at least one set of binary interaction parameters.",
                                               "Mixing rule must be \"linear\" or \"Lorentz-Berthelot\".",
                                               "Specified binary pair already exists.",
                                               "No solution found for the given inputs and fluid.",
                                               "CoolProp Issue: Use get_global_param_string(\"errstring\") for more info.",
                                               "Error Count - Not Used"};

// Helper: allocate Mathcad string and copy contents
static char* AllocMathcadString(const std::string& s)
{
    // Must use MathcadAllocate(size) so Mathcad can track and release the memory properly.
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);
    if (c != nullptr) {
        // copy s into c, this process avoids the const-cast type which would result from instead
        // converting the string using s.c_str()
        // memcpy is fine for this because std::string in C++11+ is contiguous
        std::memcpy(c, s.data(), s.size());
        c[s.size()] = '\0';
    }
    return c;
}

// Helper: allocate Mathcad array and copy a vector-of-vectors into its real part
static LRESULT AllocateToMathcadArray(LPCOMPLEXARRAY dest, const std::vector<std::vector<double>>& Vec) {
    // Expect Vec to be non-empty (caller checks ValidNumber(Vec[0][0]) earlier)
    const size_t rows = Vec.size();
    const size_t cols = (rows > 0) ? Vec[0].size() : 0;

    // Make sure dest is an actual pointer before trying to allocate memory for it
    if (dest == nullptr) {
        return MAKELRESULT(INSUFFICIENT_MEMORY, 0);    // Return error if dest is not a valid pointer
    }

    // Allocate Mathcad array memory (real part only)
    if (!MathcadArrayAllocate(dest,
                              static_cast<int>(rows),  // rows = number of output variables
                              static_cast<int>(cols),  // cols = number of input values
                              TRUE,                    // allocate memory for the real part
                              FALSE))                  // do not allocate memory for the imaginary part
    {
        return MAKELRESULT(INSUFFICIENT_MEMORY, 1);
    }

    // Copy contents into the allocated real matrix
    for (size_t irow = 0; irow < rows; ++irow) {
        const auto& row = Vec[irow];
        // Assume each inner vector has 'cols' elements (consistent with PropsSImulti contract)
        for (size_t icol = 0; icol < cols; ++icol) {
            dest->hReal[icol][irow] = row[icol];  // Note the indexing order for Mathcad arrays: hReal[column][row]
        }
    }

    return 0;
}

// Helper: Get IEEE 754 double precision NaN value for returning in case of errors in array outputs
static double get_nan() {
    unsigned long long nan_pattern = 0xFFF8000000000000ULL;
    return *(double*)&nan_pattern;
}


// Helper: check that a complex scalar input is Real and return proper Mathcad error
static inline LRESULT CheckRealOrError(LPCCOMPLEXSCALAR val, int position)
{
    if (val->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, position);
    return 0;
}

// Helper: check that a complex array input is Real and return proper Mathcad error
static inline LRESULT CheckRealArrayOrError(LPCCOMPLEXARRAY val, int position) {
    if (val->cols != 1) return MAKELRESULT(ONLY_ONE_COLUMN, position);
    if (val->hImag != NULL) return MAKELRESULT(MUST_BE_REAL, position);
    return 0;
}

// Helper: Determine which of the valid delimiters were passed in a String.
static inline char FindDelimiter(std::string instr) {
    unsigned int nDelims = 0;  // Delimiter count
    static char dType =
      '\0';  // Initialize to 'null' character, which will be returned if no delimiters are found.
    static const char dels[] = {' ', ',', '&', ';', '|'};  // Valid delimiters list

    for (char d : dels) {                            // Loop through valid delimiters
        if (instr.find(d) != std::string::npos) {    //   If found
            nDelims++;                               //      increment 
            dType = d;                               //      and save the delimiter type
        }
    }

    if (nDelims > 1) dType = 'E';  // Multiple delimiters found, return 'null'E' character to trigger error

    return dType;  // Return the delimiter type found, or '\0' if none, or 'E' if multiple
   }


// This code executes the user function CP_get_global_param_string, which is a wrapper for
// the CoolProp.get_global_param_string() function, used to get a global string parameter from CoolProp
static LRESULT CP_get_global_param_string(LPMCSTRING ParamValue,  // output (value of parameter)
                                          LPCMCSTRING ParamName)  // name of parameter (string) to retrieve
{
    std::string s;
    // Invoke the std::string form of get_global_param_string() function, save result to a new string s
    try {
        s = CoolProp::get_global_param_string(ParamName->str);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("parameter") != std::string::npos)
            return MAKELRESULT(BAD_PARAMETER, 1);
        else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release, using Helper routine above
    char* c = AllocMathcadString(s);
    // assign the string to the function's output parameter
    ParamValue->str = c;

    // normal return
    return 0;
}

// This code executes the user function CP_get_fluid_param_string, which is a wrapper for
// the CoolProp.get_fluid_param_string() function, used to get a fluid string parameter from CoolProp
static LRESULT CP_get_fluid_param_string(LPMCSTRING ParamValue,  // output (value of parameter)
                                  LPCMCSTRING FluidName,  // name of fluid (string) to retrieve
                                  LPCMCSTRING ParamName)  // name of parameter (string) to retrieve
{
    std::string s;
    // Invoke the std::string form of get_fluid_param_string() function, save result to a new string s
    try {
        s = CoolProp::get_fluid_param_string(FluidName->str, ParamName->str);
    } catch (const CoolProp::ValueError& e) {
        // MSGBOX(e.what());
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("key") != std::string::npos)
            return MAKELRESULT(BAD_FLUID, 1);
        else if (emsg.find("parameter") != std::string::npos)
            return MAKELRESULT(BAD_PARAMETER, 2);
        else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release, using Helper routine above
    char* c = AllocMathcadString(s);
    // assign the string to the function's output parameter
    ParamValue->str = c;

    // normal return
    return 0;
}

// This code executes the user function CP_set_reference_state, which is a wrapper for
// the CoolProp.set_reference_stateS() function, used to set the H/S reference states
// based on a standard state string of "IIR", "ASHRAE", "NBP", or "DEF".
static LRESULT CP_set_reference_state(LPCOMPLEXSCALAR Conf,   // output (dummy value)
                               LPCMCSTRING FluidName,  // name of fluid (string) to retrieve
                               LPCMCSTRING StateStr)   // name of standard state (string) to set
{
    // Invoke the set_reference_stateS() function, no result from this void function.
    try {
        CoolProp::set_reference_stateS(FluidName->str, StateStr->str);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);

        // Normalize the error message to lowercase to simplify substring checks and
        // avoid repeating checks for different capitalizations.
        std::string emsg_l = emsg;
        std::transform(emsg_l.begin(), emsg_l.end(), emsg_l.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (emsg_l.find("key") != std::string::npos)
            return MAKELRESULT(BAD_FLUID, 1);
        else if ((emsg_l.find("cannot use") != std::string::npos) || (emsg_l.find("temperature to qt_flash") != std::string::npos))
            return MAKELRESULT(BAD_REF, 2);
        else if (emsg_l.find("reference state") != std::string::npos)
            return MAKELRESULT(BAD_PARAMETER, 2);
        else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // assign the dummy return value
    Conf->real = 0;

    // normal return
    return 0;
}

  // Helper: centralize text-matching logic for Props1SI error handling
  static LRESULT HandleProps1SIError(const std::string& emsg, const std::string& FluidString, const std::string& PropName)
  {
    auto contains = [&](const char* s) { return emsg.find(s) != std::string::npos; };
    std::string mixName = "";  // Temp variable to hold the name of the mixture string argument for error checking, if needed
    unsigned int errPos = 0;  // Temp variable to hold the position of the error argument for returning to Mathcad, if needed

    // Check for "valid fluid" error first, since it is the most common error and can have multiple causes
    // that require parsing the error message to determine the specific error type and position.
    if (contains("valid fluid")) {
        if (contains("Neither input")) {
            if (contains("REFPROP")) {  // REFPROP or REFPROP Fluid not found
                errPos = (FluidString.find("REFPROP") != std::string::npos) ? 1u : 2u;
                return MAKELRESULT(BAD_FLUID, errPos);
            }
            if (contains("IF97")) {  // Something other than "Water" was used to IF97
                errPos = (FluidString.find("IF97") != std::string::npos) ? 1u : 2u;
                return MAKELRESULT(BAD_IF97_FLUID, errPos);
            }
            if (lower(FluidString).find(".mix") != std::string::npos) {  // Mixture string error position 1
                errPos = 1u;
                mixName = FluidString;
            } else if (lower(PropName).find(".mix") != std::string::npos) {  // Mixture string error position 2
                errPos = 2u;
                mixName = PropName;
            }
            if (!mixName.empty()) {
                Dictionary dict;
                if (!CoolProp::is_predefined_mixture(mixName, dict))  // Mixture string was used, but not found in CoolProp's predefined mixtures
                    return MAKELRESULT(BAD_MIXTURE, errPos);
                else
                    return MAKELRESULT(MISSING_BINARY_PAIR, errPos); // Likely missing a binary interaction pair.
            }
            return MAKELRESULT(BAD_FLUID, 1);
        } else {  // "Both inputs"
            return MAKELRESULT(BAD_PARAMETER, 2);
        }
    }
    // Check for "invalid parameter" errors next...
    // The parameter name will be surrounded by brackets in the error message, so add brackets to the parameter name for searching
    std::string PropNameBrackets = "[" + PropName + "]";
    if (contains("Unable to use")) {
        errPos = (emsg.find(PropNameBrackets) != std::string::npos) ? 2u : 1u;
        if (contains("Output string is invalid")) return MAKELRESULT(BAD_PARAMETER, errPos);
        if (contains("non-trivial")) return MAKELRESULT(NON_TRIVIAL, errPos);
        if (contains("No outputs")) return MAKELRESULT(NOT_AVAIL, errPos);
        return MAKELRESULT(UNKNOWN, errPos);
    }
    // Otherwise, un-trapped error. Return generic error code and have user check get_global_param_string("errstring") for details
    return MAKELRESULT(UNKNOWN, 1);
}


// This code executes the user function CP_Props1SI, which is a wrapper for
// the CoolProp.PropsSI() function, used to simply extract a
// fluid-specific parameter that is not dependent on the state
static LRESULT CP_Props1SI(LPCOMPLEXSCALAR Prop,  // pointer to the result
                           LPCMCSTRING Fluid,     // string with a valid CoolProp fluid name
                           LPCMCSTRING PropName)  // a fluid property
{
    std::string PropStr(PropName->str);
    std::string FluidString = Fluid->str;

    // pass the arguments to the CoolProp.Props1() function
    Prop->real = CoolProp::Props1SI(Fluid->str, PropName->str);

    // Note: Props1SI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if Props1SI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        return HandleProps1SIError(emsg, FluidString, PropStr);
    }

    // normal return
    return 0;
}

// Helper: centralize text-matching logic for PropsSI and PhaseSI error handling
static LRESULT HandlePropsSIError(const std::string& emsg, const std::string& Prop1Name, const unsigned int o )
{
    // Parameter "o" is passed as a parameter # offset between PropsSI (6 parameters) and PhaseSI (5 parameters and no output string)
    // PropsSI calls this routine with o = 0u (zero), while PhaseSI calls this routine with o = 1u.
    auto contains = [&](const char* s) { return emsg.find(s) != std::string::npos; };
    unsigned int errPos = 0;
    CoolProp::parameters key1;

    // Check for invalid input pair specification first.
    if (contains("Input pair variable is invalid")) {
        errPos = !is_valid_parameter(Prop1Name, key1) ? 2u - o : 4u - o;
        return MAKELRESULT(BAD_PARAMETER, errPos);
    }
    if (contains("Input Name1")) return MAKELRESULT(BAD_PARAMETER, 2-o);
    if (contains("Input Name2")) return MAKELRESULT(BAD_PARAMETER, 4-o);
    if (contains("Phase can only be specified on one")) return MAKELRESULT(ONLY_ONE_PHASE_SPEC, 4-o);
    if (contains("valid phase")) {
        errPos = !is_valid_parameter(Prop1Name, key1) ? 2u-o : 4u-o;
        return MAKELRESULT(BAD_PHASE, errPos);
    }
    if (contains("This pair of inputs")) return MAKELRESULT(BAD_INPUT_PAIR, 2-o);
    if (contains("Input vapor quality")) return (Prop1Name == "Q") ? MAKELRESULT(BAD_QUAL, 3-o) : MAKELRESULT(BAD_QUAL, 5-o);
    if (contains("Output string is invalid")) return MAKELRESULT(BAD_PARAMETER, 1);
    if (contains("not valid in two phase region")) return MAKELRESULT(TWO_PHASE, 1);
    if (contains("only defined within the two-phase")) return MAKELRESULT(NON_TWO_PHASE, 1);
    if (contains("not implemented")) return MAKELRESULT(NOT_AVAIL, 1);

    // Check for invalid fluid errors next.
    if (contains("Predefined mixture") && contains("not found")) return MAKELRESULT(BAD_MIXTURE, 6u-o);
    if (contains("Initialize failed")) {
        errPos = 6u -o;
        if (contains("Could not match the binary pair")) return MAKELRESULT(MISSING_BINARY_PAIR, errPos);
        if (contains("REFPROP")) {
            if (contains("cannot use")) return MAKELRESULT(NO_REFPROP, errPos);
            return MAKELRESULT(BAD_FLUID, errPos);
        }
        if (contains("IF97")) return MAKELRESULT(BAD_IF97_FLUID, errPos);
        return MAKELRESULT(BAD_FLUID, errPos);
    }

    // Check for out of range errors next.
    if (contains("Temperature") || contains("below Tmelt(p)"))        // Cases where temperature is out of range.
        return (Prop1Name == "T") ? MAKELRESULT(T_OUT_OF_RANGE, 3u - o) : MAKELRESULT(T_OUT_OF_RANGE, 5u - o);
    if (contains("Saturation pressure"))                              // Cases at saturation for T-P
        return (Prop1Name == "P") ? MAKELRESULT(TP_SATURATION, 3u - o) : MAKELRESULT(TP_SATURATION, 5u - o);
    if (contains("Pressure") || contains("melting line T(p)"))        // Cases where pressure is out of range
        return (Prop1Name == "P") ? MAKELRESULT(P_OUT_OF_RANGE, 3u-o) : MAKELRESULT(P_OUT_OF_RANGE, 5u-o);
    if (contains("Enthalpy") || contains("solution because Hmolar"))  // Cases where enthalpy is out of range
        return ((Prop1Name == "H") || (Prop1Name == "Hmolar")) ? MAKELRESULT(H_OUT_OF_RANGE, 3u-o) : MAKELRESULT(H_OUT_OF_RANGE, 5u-o);
    if (contains("Entropy") || contains("solution because Smolar"))   // Cases where entropy is out of range
        return ((Prop1Name == "S") || (Prop1Name == "Smolar")) ? MAKELRESULT(S_OUT_OF_RANGE, 3u-o) : MAKELRESULT(S_OUT_OF_RANGE, 5u-o);

    // Un-trapped error, return generic error code and have user check get_global_param_string("errstring") for details
    return MAKELRESULT(UNKNOWN, 1);
}

// This code executes the user function CP_PropsSI, which is a wrapper for
// the CoolProp.PropsSI() function, used to extract a fluid-specific parameter that is dependent on the state
static LRESULT CP_PropsSI(LPCOMPLEXSCALAR Prop,         // pointer to the result
                   LPCMCSTRING OutputName,       // string with a valid CoolProp OutputName
                   LPCMCSTRING InputName1,       // CoolProp InputName1
                   LPCCOMPLEXSCALAR InputProp1,  // CoolProp InputProp1
                   LPCMCSTRING InputName2,       // CoolProp InputName2
                   LPCCOMPLEXSCALAR InputProp2,  // CoolProp InputProp2
                   LPCMCSTRING FluidName)        // CoolProp Fluid
{
    // unsigned int errPos = 0;
    std::string Prop1Name(InputName1->str);
    std::string Prop2Name(InputName2->str);
    std::string FluidString = FluidName->str;

    // check that the first scalar argument is real
    LRESULT r = CheckRealOrError(InputProp1, 3);
    if (r) return r;

    // check that the second scalar argument is real
    r = CheckRealOrError(InputProp2, 5);
    if (r) return r;

    // pass the arguments to the CoolProp.PropsSI() function
    Prop->real = CoolProp::PropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);

    // Note: PropsSI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if PropsSI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        return HandlePropsSIError(emsg, Prop1Name, 0u);
    }
    // normal return
    return 0;
}


// This code executes the user function CP_PropsSImulti, which is a wrapper for
// the CoolProp.PropsSImulti() function, used to extract a range of fluid-specific parameter dependent on the state ranges
static LRESULT CP_PropsSImulti(LPCOMPLEXARRAY Prop,          // pointer to the result matrix
                        LPCMCSTRING OutputName,       // string with delimited, valid CoolProp OutputName substrings for each output column
                        LPCMCSTRING InputName1,       // CoolProp InputName1
                        LPCCOMPLEXARRAY InputProp1,   // CoolProp InputProp1 (Array)
                        LPCMCSTRING InputName2,       // CoolProp InputName2
                        LPCCOMPLEXARRAY InputProp2,   // CoolProp InputProp2 (Array
                        LPCMCSTRING FluidName)        // CoolProp Fluid
{
    // unsigned int errPos = 0;
    std::string Prop1Name(InputName1->str);
    std::string Prop2Name(InputName2->str);
    const std::string OutStr(OutputName->str);
    const std::string FluidString(FluidName->str);

    // check that the first and second scalar arguments are real
    LRESULT r = CheckRealArrayOrError(InputProp1, 3);
    if (r) return r;

    r = CheckRealArrayOrError(InputProp2, 4);
    if (r) return r;

    // Convert the input arrays to vectors for passing to PropsSImulti.
    // The PropsSImulti function expects std::vector<double> for the input properties,
    std::vector<double> Prop1Vec(InputProp1->hReal[0], InputProp1->hReal[0] + InputProp1->rows);
    std::vector<double> Prop2Vec(InputProp2->hReal[0], InputProp2->hReal[0] + InputProp2->rows);

    //Parse the OutputName string into a vector of output names, splitting on the delimiter character
    const char del = FindDelimiter(OutStr);
    if (del == 'E') {  // Multiple delimiters found, return error
        CoolProp::set_error_string("Multiple delimiters found in OutputName string. Please use only one of the following delimiter types: space, comma, semicolon, ampersand, or pipe.");
        return MAKELRESULT(BAD_PARAMETER, 1);
    }
    const std::vector<std::string> OutNames = strsplit(OutStr, del);

    //Parse the fluid string to check for multiple fluids for multi-fluid support, splitting on the delimiter character
    std::string backend, fluid;
    CoolProp::extract_backend(FluidName->str, backend, fluid);
    std::vector<double> fractions(1, 1.0);
    // extract_fractions checks for has_fractions_in_string / has_solution_concentration; no need to double check
    const std::string fluid_string = CoolProp::extract_fractions(fluid, fractions);

    // With vectors obtained, pass the parameters to the CoolProp.PropsSI() function
    std::vector<std::vector<double>> IO =
      CoolProp::PropsSImulti(OutNames, InputName1->str, Prop1Vec, InputName2->str, Prop2Vec, backend, strsplit(fluid_string, '&'), fractions);

    // Note: PropsSImulti does not throw value exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE for all values that fail.
    // If there is only one input point and one output the return matrix with be empty and
    // we can handle the error with the same logic as PropsSI.
    if (IO.empty() || IO[0].empty()) {
        std::string emsg = CoolProp::get_global_param_string("errstring"); // Also clears the error string
        CoolProp::set_error_string(emsg);  // Reset error string so Mathcad can retrieve it
        // MessageBoxA(hwndDlg, emsg.c_str(), "CoolProp PropsSImulti Error", MB_OK | MB_ICONERROR);  // Pop up the error for debugging
        return HandlePropsSIError(emsg, Prop1Name, 0u);  // Show error without a parameter offset (0u).
    }

    // Return any _HUGE values as NaN to Mathcad to use NaN filtering in arrays.
    // PropsSImulti returns _HUGE for any output value that fails, and sets the error string accordingly,
    // so this should not interfere with valid outputs.  Use the get_nan() helper function above to ensure
    // we get a proper NaN value that is recognized as such by Mathcad.
    double NaN = get_nan();
    for (auto& row : IO) {
        for (auto& val : row) {
            if (!ValidNumber(val)) {
                val = NaN;
            }
        }
    }

    // Copy the results from the IO vector of vectors into the output complex array Prop.
    // Must use MathcadArrayAllocate() so Mathcad can track and release the memory properly, using Helper routine above
    LRESULT rc = AllocateToMathcadArray(Prop, IO);
    if (rc) return rc;
    
    // normal return
    return 0;
}


// This code executes the user function CP_PhaseSI, which is a wrapper for
// the CoolProp.PhaseSI() function, used to the fluid phase dependent on the state
static LRESULT CP_PhaseSI(LPMCSTRING PhaseStr,          // output (CoolProp phase string)
                   LPCMCSTRING InputName1,       // CoolProp InputName1
                   LPCCOMPLEXSCALAR InputProp1,  // CoolProp InputProp1
                   LPCMCSTRING InputName2,       // CoolProp InputName2
                   LPCCOMPLEXSCALAR InputProp2,  // CoolProp InputProp2
                   LPCMCSTRING FluidName)        // CoolProp Fluid
{
    std::string ph;
    std::string Prop1Name(InputName1->str);
    std::string Prop2Name(InputName2->str);
    std::string FluidString = FluidName->str;

    // check that the first scalar argument is real
    LRESULT r = CheckRealOrError(InputProp1, 2);
    if (r) return r;

    // check that the second scalar argument is real
    r = CheckRealOrError(InputProp2, 4);
    if (r) return r;

    // pass the arguments to the CoolProp.phaseSI() function
    ph = CoolProp::PhaseSI(InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);

    // Note: PhaseSI does not throw exceptions, but instead appends the global parameter "errstring"
    // to the returned phase string if there is an error and returns "unknown" for the phase string.
    //
    // Use string search to see if PhaseSI failed with an error message (appended ": <error>")...
    if (ph.find("unknown:") != std::string::npos) {
        std::string emsg = ph.substr(ph.find(":")+2,ph.length()-ph.find(":")+2);
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        // Use PropsSI error handling logic to determine the specific error message but with adjusted argument
        // positions for PhaseSI (which has no OutputName argument and thus shifts the positions of the InputName
        // and InputProp arguments by 1)
        return HandlePropsSIError(emsg, Prop1Name, 1u);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release, using Helper routine above
    char* c = AllocMathcadString(ph);
    // assign the string to the function's output parameter
    PhaseStr->str = c;

    // normal return
    return 0;
}


// This code executes the user function CP_HAPropsSI, which is a wrapper for
// the CoolProp.HAPropsSI() function, used to extract humid air properties in base-SI units
static LRESULT CP_HAPropsSI(LPCOMPLEXSCALAR Prop,         // pointer to the result
                     LPCMCSTRING OutputName,       // string with a valid CoolProp Output Name
                     LPCMCSTRING InputName1,       // CoolProp InputName1
                     LPCCOMPLEXSCALAR InputProp1,  // CoolProp InputProp1
                     LPCMCSTRING InputName2,       // CoolProp InputName2
                     LPCCOMPLEXSCALAR InputProp2,  // CoolProp InputProp2
                     LPCMCSTRING InputName3,       // CoolProp InputName3
                     LPCCOMPLEXSCALAR InputProp3)  // CoolProp InputProp3
{
    unsigned int errPos = 0;
    std::string OutName(OutputName->str);
    OutName = "[" + OutName + "]";
    std::string Prop1Name(InputName1->str);
    Prop1Name = "[" + Prop1Name + "]";
    std::string Prop2Name(InputName2->str);
    Prop2Name = "[" + Prop2Name + "]";
    std::string Prop3Name(InputName3->str);
    Prop3Name = "[" + Prop3Name + "]";

    // check that the 3 scalar arguments are real or throw error using inline helper function above
    LRESULT r = CheckRealOrError(InputProp1, 3);
    if (r) return r;

    r = CheckRealOrError(InputProp2, 5);
    if (r) return r;

    r = CheckRealOrError(InputProp3, 7);
    if (r) return r;

    // pass the arguments to the HumidAirProp.HAProps() function
    Prop->real =
      HumidAir::HAPropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, InputName3->str, InputProp3->real);

    // Note: HAPropsSI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if HAPropsSI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it

        unsigned int errPos = 0;
        if (emsg.find(OutName) != std::string::npos)
            errPos = 1;
        else if (emsg.find(Prop1Name) != std::string::npos)
            errPos = 2;
        else if (emsg.find(Prop2Name) != std::string::npos)
            errPos = 4;
        else if (emsg.find(Prop3Name) != std::string::npos)
            errPos = 6;

        if (errPos != 0)
            return MAKELRESULT(BAD_PARAMETER, errPos);
        else if (emsg.find("at least one of the variables") != std::string::npos)
            return MAKELRESULT(HA_INPUTS, 2);
        else
            // Only the Generic Error message supported at this time for HAPropsSI for any other errors
            return MAKELRESULT(UNKNOWN, 1);
    }

    // normal return
    return 0;
}

// ********************************************************************************************************
//   Pseudo-Low-Level Functions for mixture information & settings
// ********************************************************************************************************

// this code executes the user function CP_get_mixture_binary_pair_data, which is a wrapper for
// the CoolProp.get_mixture_binary_pair_data() function, used to get the requested binary pair
// interaction parameter (always returned as a string).
static LRESULT CP_get_mixture_binary_pair_data(LPMCSTRING Value,  // output string (string contains value of parameter)
                                        LPCMCSTRING CAS1,  // First component
                                        LPCMCSTRING CAS2,  // Second component
                                        LPCMCSTRING Key)   // name of the binary pair parameter (string) to retrieve
{
    std::string s;
    // Invoke the std::string form of get_global_param_string() function, save result to a new string s
    try {
        s = CoolProp::get_mixture_binary_pair_data(CAS1->str, CAS2->str, Key->str);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("parameter") != std::string::npos)
            return MAKELRESULT(BAD_PARAMETER, 3);
        else if (emsg.find("binary pair") != std::string::npos)
            return MAKELRESULT(BAD_BINARY_PAIR, 1);
        else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release, using Helper routine above
    char* c = AllocMathcadString(s);
    // assign the string to the function's output parameter
    Value->str = c;

    // normal return
    return 0;
}

// Helper: centralize check for predefined-mixture-not-found error messages
static inline bool PredefinedMixtureNotFound(const std::string& emsg) {
    return (emsg.find("Predefined mixture") != std::string::npos && emsg.find("not found") != std::string::npos);
}

// Return semicolon-delimited fluid names for a predefined mixture
static LRESULT CP_get_predefined_mixture_fluids(LPMCSTRING Comps,         // output string (semicolon-delimited fluid names)
                                                LPCMCSTRING MixtureName)  // name of predefined mixture (with or without .mix)
{
    std::string s;
    try {
        std::string key = MixtureName->str;
        if ((key.find(".mix") == std::string::npos) && (key.find(".MIX") == std::string::npos)) key += ".mix";

        Dictionary dict;
        if (!CoolProp::is_predefined_mixture(key, dict)) {
            throw CoolProp::ValueError(format("Predefined mixture [%s] not found", key.c_str()));
        }

        const std::vector<std::string>& fluids = dict.get_string_vector("fluids");
        s = strjoin(fluids, ";");
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (PredefinedMixtureNotFound(emsg))
            return MAKELRESULT(BAD_MIXTURE, 1);
        else
            return MAKELRESULT(UNKNOWN, 1);
    } catch (...) {
        CoolProp::set_error_string("Unknown error retrieving predefined mixture fluids");
        return MAKELRESULT(UNKNOWN, 1);
    }

    char* c = AllocMathcadString(s);
    Comps->str = c;
    return 0;
}

// Return a column vector of mole fractions for a predefined mixture
static LRESULT CP_get_predefined_mixture_mole_fractions(LPCOMPLEXARRAY Dest,      // output column vector of mole fractions
                                                        LPCMCSTRING MixtureName)  // name of predefined mixture (with or without .mix)
{
    std::vector<std::vector<double>> Vec;
    try {
        std::string key = MixtureName->str;
        if ((key.find(".mix") == std::string::npos) && (key.find(".MIX") == std::string::npos)) key += ".mix";

        Dictionary dict;
        if (!CoolProp::is_predefined_mixture(key, dict)) {
            throw CoolProp::ValueError(format("Predefined mixture [%s] not found", key.c_str()));
        }

        const std::vector<double>& mf = dict.get_double_vector("mole_fractions");
        // Convert to vector-of-vectors with one column
        Vec.reserve(mf.size());
        for (std::size_t i = 0; i < mf.size(); ++i) {
            Vec.push_back(std::vector<double>(1, mf[i]));
        }
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (PredefinedMixtureNotFound(emsg)) {
            return MAKELRESULT(BAD_MIXTURE, 1);
        }
        return MAKELRESULT(UNKNOWN, 1);
    } catch (...) {
        std::string emsg("Unknown error retrieving predefined mixture mole fractions");
        CoolProp::set_error_string(emsg);
        return MAKELRESULT(UNKNOWN, 1);
    }

    if (Vec.empty()) {
        // No mole fractions found -> return mixture error
        std::string emsg("No mole fractions found for specified predefined mixture");
        CoolProp::set_error_string(emsg);
        return MAKELRESULT(BAD_MIXTURE, 1);
    }

    LRESULT rc = AllocateToMathcadArray(Dest, Vec);
    if (rc) return rc;

    return 0;
}


// This code executes the user function CP_apply_simple_mixing_rule, which is a wrapper for
// the CoolProp.apply_simple_mixing_rule() function, used to set the mixing rule for a
// specific binary pair.
static LRESULT CP_apply_simple_mixing_rule(LPMCSTRING Msg,    // output string (verification message)
                                           LPCMCSTRING CAS1,  // First component
                                           LPCMCSTRING CAS2,  // Second component
                                           LPCMCSTRING Rule)  // Mixing rule, either 'linear' or 'Lorentz-Berthelot'
{
    std::string s = Rule->str;
    s.append(" mixing rule set.");
    // Invoke the std::string form of get_global_param_string() function, save result to a new string s
    try {
        CoolProp::apply_simple_mixing_rule(CAS1->str, CAS2->str, Rule->str);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("simple mixing rule") != std::string::npos) {
            return MAKELRESULT(BAD_RULE, 3);
        } else if (emsg.find("already in") != std::string::npos) {
            return MAKELRESULT(PAIR_EXISTS, 1);
        } else if (emsg.find("key") != std::string::npos) {
            if (emsg.find(CAS1->str) != std::string::npos) {
                return MAKELRESULT(BAD_FLUID, 1);
            } else if (emsg.find(CAS2->str) != std::string::npos) {
                return MAKELRESULT(BAD_FLUID, 2);
            } else
                return MAKELRESULT(UNKNOWN, 1);
        } else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release
    char* c = AllocMathcadString(s);
    // assign the string to the function's output parameter
    Msg->str = c;

    // normal return
    return 0;
}

// This code executes the user function CP_set_mixture_binary_pair_data, which is a wrapper for
// the CoolProp.set_mixture_binary_pair_data() function, used to set the mixing rule for a
// specific binary pair.
static LRESULT CP_set_mixture_binary_pair_data(LPMCSTRING Msg,          // output string (verification message)
                                        LPCMCSTRING CAS1,        // First component
                                        LPCMCSTRING CAS2,        // Second component
                                        LPCMCSTRING Param,       // Parameter Name String to set
                                        LPCCOMPLEXSCALAR Value)  // Parameter Value
{
    std::string s = Param->str;
    s.append(" parameter set.");

    // check that the first scalar argument is real
    LRESULT r = CheckRealOrError(Value, 4);
    if (r) return r;

    // Invoke the std::string form of get_global_param_string() function, save result to a new string s
    try {
        CoolProp::set_mixture_binary_pair_data(CAS1->str, CAS2->str, Param->str, Value->real);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("parameter") != std::string::npos) {
            return MAKELRESULT(BAD_PARAMETER, 3);
        } else if (emsg.find("key") != std::string::npos) {
            if (emsg.find(CAS1->str) != std::string::npos)
                return MAKELRESULT(BAD_FLUID, 1);
            else
                return MAKELRESULT(BAD_FLUID, 2);
        } else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // Must use MathcadAllocate(size) so Mathcad can track and release
    char* c = AllocMathcadString(s);
    // assign the string to the function's output parameter
    Msg->str = c;

    // normal return
    return 0;
}

// ********************************************************************************************************
// Fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
// ********************************************************************************************************

FUNCTIONINFO PropsParam = {
  "get_global_param_string",                                 // Name by which Mathcad will recognize the function
  "Name of the parameter to retrieve",                       // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_global_param_string,                   // Pointer to the function code.
  MC_STRING,                                                 // Returns a Mathcad string
  1,                                                         // Number of arguments
  {MC_STRING}                                                // Argument types
};

FUNCTIONINFO FluidParam = {
  "get_fluid_param_string",                                  // Name by which Mathcad will recognize the function
  "Fluid, Name of the parameter to retrieve",                // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_fluid_param_string,                    // Pointer to the function code.
  MC_STRING,                                                 // Returns a Mathcad string
  2,                                                         // Number of arguments
  {MC_STRING, MC_STRING}                                     // Argument types
};

FUNCTIONINFO RefState = {
  "set_reference_state",                                           // Name by which Mathcad will recognize the function
  "Fluid, Reference State String",                                 // Description of input parameters
  "Sets the reference state to either IIR, ASHRAE, NBP, or DEF.",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_set_reference_state,                             // Pointer to the function code.
  COMPLEX_SCALAR,                                                  // Returns a Mathcad complex scalar
  2,                                                               // Number of arguments
  {MC_STRING, MC_STRING}                                           // Argument types
};

FUNCTIONINFO Props1SI = {
  "Props1SI",                                                                                     // Name by which Mathcad will recognize the function
  "Fluid, Property Name",                                                                         // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is not dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_Props1SI,                                                                       // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                 // Returns a Mathcad complex scalar
  2,                                                                                              // Number of arguments
  {MC_STRING, MC_STRING}                                                                          // Argument types
};

FUNCTIONINFO PropsSI = {
  "PropsSI",                                                                                  // Name by which Mathcad will recognize the function
  "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name",  // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_PropsSI,                                                                    // Pointer to the function code.
  COMPLEX_SCALAR,                                                                             // Returns a Mathcad complex scalar
  6,                                                                                          // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING}                // Argument types
};

FUNCTIONINFO PropsSImulti = {
  "PropsSImulti",                                                                                  // Name by which Mathcad will recognize the function
  "Output Name, Input Name 1, Input Property 1 (Array), Input Name 2, Input Property 2 (Array), Fluid Name",  // Description of input parameters
  "Returns a range of fluid-specific parameters, where the parameters are dependent on the state ranges defined by the input property arrays",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_PropsSImulti,                                                                    // Pointer to the function code.
  COMPLEX_ARRAY,                                                                                   // Returns a Mathcad complex array
  6,                                                                                               // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_ARRAY, MC_STRING, COMPLEX_ARRAY, MC_STRING}                       // Argument types
};

FUNCTIONINFO PhaseSI = {
  "PhaseSI",                                                                                  // Name by which Mathcad will recognize the function
  "Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name",               // Description of input parameters
  "Returns the fluid phase, dependent on the fluid state",                                    // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_PhaseSI,                                                                    // Pointer to the function code.
  MC_STRING,                                                                                  // Returns a Mathcad String
  5,                                                                                          // Number of arguments
  {MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING}                           // Argument types
};

FUNCTIONINFO HAPropsSI = {
  "HAPropsSI",  // Name by which Mathcad will recognize the function
  "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3",  // Description of input parameters
  "Returns a parameter of humid air, where the parameter is dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_HAPropsSI,                                                                // Pointer to the function code.
  COMPLEX_SCALAR,                                                                           // Returns a Mathcad complex scalar
  7,                                                                                        // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR}  // Argument types
};

FUNCTIONINFO GetMixtureData = {
  "get_mixture_binary_pair_data",                            // Name by which Mathcad will recognize the function
  "CAS 1, CAS 2, Name of the parameter to retrieve",         // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_mixture_binary_pair_data,              // Pointer to the function code.
  MC_STRING,                                                 // Returns a Mathcad string
  3,                                                         // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING}                          // Argument types
};

FUNCTIONINFO ApplyMixingRule = {
  "apply_simple_mixing_rule",                   // Name by which Mathcad will recognize the function
  "CAS 1, CAS 2, Mixing Rule",                  // Description of input parameters
  "Sets a simple mixing rule for binary pair",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_apply_simple_mixing_rule,     // Pointer to the function code.
  MC_STRING,                                    // Returns a Mathcad string
  3,                                            // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING}             // Argument types
};

FUNCTIONINFO SetMixtureData = {
  "set_mixture_binary_pair_data",                             // Name by which Mathcad will recognize the function
  "CAS 1, CAS 2, Parameter Name, Parameter value",            // Description of input parameters
  "Sets the value of the specified binary mixing parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_set_mixture_binary_pair_data,               // Pointer to the function code.
  MC_STRING,                                                  // Returns a Mathcad string
  4,                                                          // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING, COMPLEX_SCALAR}           // Argument types
};

FUNCTIONINFO GetPredefFluids = {
  "get_predefined_mixture_fluids",                                      // Name by which Mathcad will recognize the function
  "MixtureName",                                                        // Description of input parameters
  "Returns semicolon-delimited fluid names in the predefined mixture",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_predefined_mixture_fluids,                        // Pointer to the function code.
  MC_STRING,                                                            // Returns a Mathcad string
  1,                                                                    // Number of arguments
  {MC_STRING}                                                           // Argument types
};

FUNCTIONINFO GetPredefMoleFracs = {
  "get_predefined_mixture_fractions",                                      // Name by which Mathcad will recognize the function
  "MixtureName",                                                           // Description of input parameters
  "Returns a column vector of mole fractions for the predefined mixture",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_predefined_mixture_mole_fractions,                   // Pointer to the function code.
  COMPLEX_ARRAY,                                                           // Returns a Mathcad complex array
  1,                                                                       // Number of arguments
  {MC_STRING}                                                              // Argument types
};


// ************************************************************************************
// DLL entry point code.
// ************************************************************************************
// Note 1: The _CRT_INIT function is needed if using Microsoft's WIN32 library and you want to use the C runtime library,
//         which we do for string handling and other utilities.  If you are using a different compiler or setup,
//         this may not be needed. It is only defined here if _WIN32 is defined.
// Note 2: The the _CRT_INIT function is not called directly in our code, but is instead called by the DllEntryPoint function
// Note 3: MS IntelliSense will warn "VCR001 Definition for _CRT_INIT not found", but this warning can be safely
//         ignored as a known false positive.  To disable, got to Tools > Options > Text Editor > C/C++ > IntelliSense and
//         set "Create declaration/definition suggestion level" to "Refactoring only".
#ifdef _WIN32
extern "C" BOOL WINAPI _CRT_INIT(HINSTANCE hinstDLL, DWORD dwReason, LPVOID lpReserved);
#endif

extern "C" BOOL WINAPI DllEntryPoint(HINSTANCE hDLL, DWORD dwReason, LPVOID lpReserved) {
    switch (dwReason) {
        case DLL_PROCESS_ATTACH:
            //
            // DLL is attaching to the address space of the current process.
            //
            if (!_CRT_INIT(hDLL, dwReason, lpReserved)) return FALSE;

            // Register the error message table
            if (!CreateUserErrorMessageTable(hDLL, NUMBER_OF_ERRORS, CPErrorMessageTable)) break;

            // ...and if the errors register OK, go ahead and register user function
            CreateUserFunction(hDLL, &PropsParam);
            CreateUserFunction(hDLL, &FluidParam);
            CreateUserFunction(hDLL, &RefState);
            CreateUserFunction(hDLL, &Props1SI);
            CreateUserFunction(hDLL, &PropsSI);
            CreateUserFunction(hDLL, &PropsSImulti);
            CreateUserFunction(hDLL, &PhaseSI);
            CreateUserFunction(hDLL, &HAPropsSI);
            CreateUserFunction(hDLL, &GetMixtureData);
            CreateUserFunction(hDLL, &SetMixtureData);
            CreateUserFunction(hDLL, &ApplyMixingRule);
            // Register the new helper functions for predefined mixtures
            CreateUserFunction(hDLL, &GetPredefFluids);
            CreateUserFunction(hDLL, &GetPredefMoleFracs);
            break;

        case DLL_THREAD_ATTACH:
        case DLL_THREAD_DETACH:
        case DLL_PROCESS_DETACH:

            if (!_CRT_INIT(hDLL, dwReason, lpReserved)) {
                Sleep(1000);  // Attempt to keep CRT_INIT from detaching before all threads are closed
                return FALSE;
            }
            break;
    }
    return TRUE;
}
