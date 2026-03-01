// CoolPropMathcad.cpp : Defines the exported functions for the DLL Add-in.
//

#include <string>
#include <cstring>
// #include <sstream>

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

namespace CoolProp {
extern void apply_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule);
}

enum EC
{
    MUST_BE_REAL = 1,  // Mathcad Error Codes       v
    INSUFFICIENT_MEMORY,
    INTERRUPTED,
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
    UNKNOWN,
    NUMBER_OF_ERRORS
};  // Dummy Code for Error Count

// table of error messages
// As of Mathcad Prime 10, these are now actually returned as Custom Error: messages
char* CPErrorMessageTable[NUMBER_OF_ERRORS] = {"Interrupted",
                                               "Insufficient Memory",
                                               "Argument must be real",
                                               "Invalid Fluid String",
                                               "Invalid predefined mixture or Binary Interaction Parameters Missing",
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
                                               "Missing at least one set of binary interaction parameters. Use get_global_param_string(\"errstring\") for more info.",
                                               "Mixing rule must be \"linear\" or \"Lorentz-Berthelot\".",
                                               "Specified binary pair already exists.",
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

// Helper: check that a complex scalar input is Real and return proper Mathcad error
static inline LRESULT CheckRealOrError(LPCCOMPLEXSCALAR val, int position)
{
    if (val->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, position);
    return 0;
}

// this code executes the user function CP_get_global_param_string, which is a wrapper for
// the CoolProp.get_global_param_string() function, used to get a global string parameter from CoolProp
LRESULT CP_get_global_param_string(LPMCSTRING ParamValue,  // output (value of parameter)
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

// this code executes the user function CP_get_fluid_param_string, which is a wrapper for
// the CoolProp.get_fluid_param_string() function, used to get a fluid string parameter from CoolProp
LRESULT CP_get_fluid_param_string(LPMCSTRING ParamValue,  // output (value of parameter)
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

// this code executes the user function CP_set_reference_state, which is a wrapper for
// the CoolProp.set_reference_stateS() function, used to set the H/S reference states
// based on a standard state string of "IIR", "ASHRAE", "NBP", or "DEF".
LRESULT CP_set_reference_state(LPCOMPLEXSCALAR Conf,   // output (dummy value)
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
  static LRESULT HandleProps1SIError(const std::string& emsg, const std::string& FluidString, const std::string& PropNameBrackets)
  {
    auto contains = [&](const char* s) { return emsg.find(s) != std::string::npos; };
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
            std::string fluid_l = lower(FluidString);
            if (fluid_l.find(".mix") != std::string::npos) return MAKELRESULT(BAD_MIXTURE, 1);  // Mixture string error position 1
            std::string prop_l = lower(PropNameBrackets);
            if (prop_l.find(".mix") != std::string::npos) return MAKELRESULT(BAD_MIXTURE, 2);  // Mixture string error position 2
            return MAKELRESULT(BAD_FLUID, 1);
        } else {  // "Both inputs"
            return MAKELRESULT(BAD_PARAMETER, 2);
        }
    }
    // Check for "invalid parameter" errors next...
    if (contains("Unable to use")) {
        errPos = (emsg.find(PropNameBrackets) != std::string::npos) ? 2u : 1u;
        if (contains("Output string is invalid")) return MAKELRESULT(BAD_PARAMETER, errPos);
        if (contains("non-trivial")) return MAKELRESULT(NON_TRIVIAL, errPos);
        if (contains("No outputs")) return MAKELRESULT(NOT_AVAIL, errPos);
        return MAKELRESULT(UNKNOWN, errPos);
    }
    //
    return MAKELRESULT(UNKNOWN, 1);
}


// this code executes the user function CP_Props1SI, which is a wrapper for
// the CoolProp.PropsSI() function, used to simply extract a
// fluid-specific parameter that is not dependent on the state
LRESULT CP_Props1SI(LPCOMPLEXSCALAR Prop,  // pointer to the result
                    LPCMCSTRING Fluid,     // string with a valid CoolProp fluid name
                    LPCMCSTRING PropName)  // a fluid property
{
    std::string PropNameBrackets(PropName->str);
    PropNameBrackets = "[" + PropNameBrackets + "]";
    std::string FluidString = Fluid->str;

    // pass the arguments to the CoolProp.Props1() function
    Prop->real = CoolProp::Props1SI(Fluid->str, PropName->str);

    // Note: Props1SI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if Props1SI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        return HandleProps1SIError(emsg, FluidString, PropNameBrackets);
    }

    // normal return
    return 0;
}

// Helper: centralize text-matching logic for PropsSI error handling
static LRESULT HandlePropsSIError(const std::string& emsg, const std::string& Prop1Name, CoolProp::parameters& key1)
{
    auto contains = [&](const char* s) { return emsg.find(s) != std::string::npos; };
    unsigned int errPos = 0;

    if (contains("Input pair variable is invalid")) {
        errPos = !is_valid_parameter(Prop1Name, key1) ? 2u : 4u;
        return MAKELRESULT(BAD_PARAMETER, errPos);
    }
    if (contains("Input Name1")) return MAKELRESULT(BAD_PARAMETER, 2);
    if (contains("Input Name2")) return MAKELRESULT(BAD_PARAMETER, 4);
    if (contains("Phase can only be specified on one")) return MAKELRESULT(ONLY_ONE_PHASE_SPEC, 4);
    if (contains("valid phase")) {
        errPos = !is_valid_parameter(Prop1Name, key1) ? 2u : 4u;
        return MAKELRESULT(BAD_PHASE, errPos);
    }
    if (contains("This pair of inputs")) return MAKELRESULT(BAD_INPUT_PAIR, 2);
    if (contains("Input vapor quality")) return (Prop1Name == "Q") ? MAKELRESULT(BAD_QUAL, 3) : MAKELRESULT(BAD_QUAL, 5);
    if (contains("Output string is invalid")) return MAKELRESULT(BAD_PARAMETER, 1);
    if (contains("not valid in two phase region")) return MAKELRESULT(TWO_PHASE, 1);
    if (contains("only defined within the two-phase")) return MAKELRESULT(NON_TWO_PHASE, 1);
    if (contains("not implemented")) return MAKELRESULT(NOT_AVAIL, 1);

    if (contains("Initialize failed")) {
        if (contains("Could not match the binary pair")) return MAKELRESULT(MISSING_BINARY_PAIR, 6);
        if (contains("REFPROP")) {
            if (contains("cannot use")) return MAKELRESULT(NO_REFPROP, 6);
            return MAKELRESULT(BAD_FLUID, 6);
        }
        if (contains("IF97")) return MAKELRESULT(BAD_IF97_FLUID, 6);
        return MAKELRESULT(BAD_FLUID, 6);
    }

    if (contains("Temperature")) return (Prop1Name == "T") ? MAKELRESULT(T_OUT_OF_RANGE, 3) : MAKELRESULT(T_OUT_OF_RANGE, 5);
    if (contains("Saturation pressure")) return (Prop1Name == "P") ? MAKELRESULT(TP_SATURATION, 3) : MAKELRESULT(TP_SATURATION, 5);
    if (contains("Pressure")) return (Prop1Name == "P") ? MAKELRESULT(P_OUT_OF_RANGE, 3) : MAKELRESULT(P_OUT_OF_RANGE, 5);
    if (contains("Enthalpy") || contains("solution because Hmolar"))
        return ((Prop1Name == "H") || (Prop1Name == "Hmolar")) ? MAKELRESULT(H_OUT_OF_RANGE, 3) : MAKELRESULT(H_OUT_OF_RANGE, 5);
    if (contains("Entropy") || contains("solution because Smolar"))
        return ((Prop1Name == "S") || (Prop1Name == "Smolar")) ? MAKELRESULT(S_OUT_OF_RANGE, 3) : MAKELRESULT(S_OUT_OF_RANGE, 5);

    return MAKELRESULT(UNKNOWN, 1);
}

// this code executes the user function CP_PropsSI, which is a wrapper for
// the CoolProp.PropsSI() function, used to extract a fluid-specific parameter that is dependent on the state
LRESULT CP_PropsSI(LPCOMPLEXSCALAR Prop,         // pointer to the result
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
    CoolProp::parameters key1;

    // check that the first scalar argument is real
    LRESULT r = CheckRealOrError(InputProp1, 3);
    if (r) return r;

    // check that the second scalar argument is real
    r = CheckRealOrError(InputProp2, 5);
    if (r) return r;

    // pass the arguments to the CoolProp.Props() function
    Prop->real = CoolProp::PropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);

    // Note: PropsSI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if PropsSI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        return HandlePropsSIError(emsg, Prop1Name, key1);
    }
    // normal return
    return 0;
}

// this code executes the user function CP_HAPropsSI, which is a wrapper for
// the CoolProp.HAPropsSI() function, used to extract humid air properties in base-SI units
LRESULT CP_HAPropsSI(LPCOMPLEXSCALAR Prop,         // pointer to the result
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

// this code executes the user function CP_get_mixture_binary_pair_data, which is a wrapper for
// the CoolProp.get_mixture_binary_pair_data() function, used to get the requested binary pair
// interaction parameter (always returned as a string).
LRESULT CP_get_mixture_binary_pair_data(LPMCSTRING Value,  // output string (string contains value of parameter)
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

// this code executes the user function CP_apply_simple_mixing_rule, which is a wrapper for
// the CoolProp.apply_simple_mixing_rule() function, used to set the mixing rule for a
// specific binary pair.
LRESULT CP_apply_simple_mixing_rule(LPMCSTRING Msg,    // output string (verification message)
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

// this code executes the user function CP_set_mixture_binary_pair_data, which is a wrapper for
// the CoolProp.set_mixture_binary_pair_data() function, used to set the mixing rule for a
// specific binary pair.
LRESULT CP_set_mixture_binary_pair_data(LPMCSTRING Msg,          // output string (verification message)
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

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO PropsParam = {
  "get_global_param_string",                                 // Name by which Mathcad will recognize the function
  "Name of the parameter to retrieve",                       // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_global_param_string,                   // Pointer to the function code.
  MC_STRING,                                                 // Returns a Mathcad string
  1,                                                         // Number of arguments
  {MC_STRING}                                                // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO FluidParam = {
  "get_fluid_param_string",                                  // Name by which Mathcad will recognize the function
  "Fluid, Name of the parameter to retrieve",                // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_fluid_param_string,                    // Pointer to the function code.
  MC_STRING,                                                 // Returns a Mathcad string
  2,                                                         // Number of arguments
  {MC_STRING, MC_STRING}                                     // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO RefState = {
  "set_reference_state",                                           // Name by which Mathcad will recognize the function
  "Fluid, Reference State String",                                 // Description of input parameters
  "Sets the reference state to either IIR, ASHRAE, NBP, or DEF.",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_set_reference_state,                             // Pointer to the function code.
  COMPLEX_SCALAR,                                                  // Returns a Mathcad complex scalar
  2,                                                               // Number of arguments
  {MC_STRING, MC_STRING}                                           // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO Props1SI = {
  "Props1SI",                                                                                     // Name by which Mathcad will recognize the function
  "Fluid, Property Name",                                                                         // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is not dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_Props1SI,                                                                       // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                 // Returns a Mathcad complex scalar
  2,                                                                                              // Number of arguments
  {MC_STRING, MC_STRING}                                                                          // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO PropsSI = {
  "PropsSI",                                                                                  // Name by which Mathcad will recognize the function
  "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name",  // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_PropsSI,                                                                    // Pointer to the function code.
  COMPLEX_SCALAR,                                                                             // Returns a Mathcad complex scalar
  6,                                                                                          // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING}                // Argument types
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

// ************************************************************************************
// DLL entry point code.
// ************************************************************************************
// The _CRT_INIT function is needed if you are using Microsoft's WIN32 compiler
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
            CreateUserFunction(hDLL, &HAPropsSI);
            CreateUserFunction(hDLL, &GetMixtureData);
            CreateUserFunction(hDLL, &SetMixtureData);
            CreateUserFunction(hDLL, &ApplyMixingRule);
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
