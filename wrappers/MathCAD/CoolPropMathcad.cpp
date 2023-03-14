// CoolPropMathcad.cpp : Defines the exported functions for the DLL Add-in.
//

#include <string>
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
    MUST_BE_REAL = 1,
    INSUFFICIENT_MEMORY,
    INTERRUPTED,  // Mathcad Error Codes
    BAD_FLUID,
    BAD_IF97_FLUID,
    BAD_PARAMETER,
    BAD_PHASE,  // CoolProp Error Codes
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
    BAD_RULE,
    PAIR_EXISTS,
    UNKNOWN,
    NUMBER_OF_ERRORS
};  // Dummy Code for Error Count

// table of error messages
// if user function never returns an
// error -- you do not need to create this
// table
char* CPErrorMessageTable[NUMBER_OF_ERRORS] = {"Interrupted",
                                               "Insufficient Memory",
                                               "Argument must be real",
                                               "Invalid Fluid String",
                                               "IF97 Backend supports pure water only",
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
                                               "Mixing rule must be \"linear\" or \"Lorentz-Berthelot\".",
                                               "Specified binary pair already exists.",
                                               "ERROR: Use get_global_param_string(\"errstring\") for more info",
                                               "Error Count - Not Used"};

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

    // Must use MathcadAllocate(size) so Mathcad can track and release
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);  // create a c-string (pointer) c with the same size as s
    // copy s into c, this process avoids the const-cast type which would result from instead
    // converting the string using s.c_str()
    std::copy(s.begin(), s.end(), c);
    c[s.size()] = '\0';
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

    // Must use MathcadAllocate(size) so Mathcad can track and release
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);  // create a c-string (pointer) c with the same size as s
    // copy s into c, this process avoids the const-cast type which would result from instead
    // converting the string using s.c_str()
    std::copy(s.begin(), s.end(), c);
    c[s.size()] = '\0';
    // assign the string to the function's output parameter
    ParamValue->str = c;

    // normal return
    return 0;
}

// this code executes the user function CP_set_reference_state, which is a wrapper for
// the CoolProp.set_reference_stateS() function, used to set the H/S reference states
// based on a standard state string of "IIR", "ASHRAE", "NBP", or "DEF".
LRESULT CP_set_reference_state(LPCOMPLEXSCALAR Conf,   // output (dummy value)
                               LPCMCSTRING FluidName,  // name of fluidr (string) to retrieve
                               LPCMCSTRING StateStr)   // name of standard state (string) to set
{
    // Invoke the set_reference_stateS() function, no result from this void function.
    try {
        CoolProp::set_reference_stateS(FluidName->str, StateStr->str);
    } catch (const CoolProp::ValueError& e) {
        std::string emsg(e.what());
        CoolProp::set_error_string(emsg);
        if (emsg.find("key") != std::string::npos)
            return MAKELRESULT(BAD_FLUID, 1);
        else if ((emsg.find("Cannot use") != std::string::npos) || (emsg.find("Temperature to QT_flash") != std::string::npos))
            return MAKELRESULT(BAD_REF, 2);
        else if ((emsg.find("reference state") != std::string::npos) || (emsg.find("Reference state") != std::string::npos)
                 || (emsg.find("Reference State") != std::string::npos))
            return MAKELRESULT(BAD_PARAMETER, 2);
        else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // assign the dummy return value
    Conf->real = 0;

    // normal return
    return 0;
}

// this code executes the user function CP_Props1SI, which is a wrapper for
// the CoolProp.PropsSI() function, used to simply extract a
// fluid-specific parameter that is not dependent on the state
LRESULT CP_Props1SI(LPCOMPLEXSCALAR Prop,  // pointer to the result
                    LPCMCSTRING Fluid,     // string with a valid CoolProp fluid name
                    LPCMCSTRING PropName)  // a fluid property
{
    unsigned int errPos;
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
        if (emsg.find("valid fluid") != std::string::npos) {
            if (emsg.find("Neither input") != std::string::npos) {
                if (emsg.find("REFPROP") != std::string::npos) {
                    // Fluid can be in either parameter location, find out which.
                    // It will be in brackets in the error message.
                    if (FluidString.find("REFPROP") != std::string::npos)
                        errPos = 1;  // [REFPROP::???] is in Fluid->str, i.e. position 1
                    else
                        errPos = 2;  // [REFPROP::???] is in PropName->str, i.e. position 2
                    return MAKELRESULT(BAD_FLUID, errPos);
                } else if (emsg.find("IF97") != std::string::npos) {
                    if (FluidString.find("IF97") != std::string::npos)
                        errPos = 1;  // [IF97::???] is in Fluid->str, i.e. position 1
                    else
                        errPos = 2;  // [IF97::???] is in PropName-str, i.e. position 2
                    return MAKELRESULT(BAD_IF97_FLUID, errPos);
                } else
                    return MAKELRESULT(BAD_FLUID, 1);
            } else  //  "Both inputs"
                return MAKELRESULT(BAD_PARAMETER, 2);
        } else if (emsg.find("Unable to use") != std::string::npos) {
            // PropName can be in either parameter location, find out which.
            // It will be in brackets in the error message.
            if (emsg.find(PropNameBrackets) != std::string::npos)
                errPos = 2;  // [PropName] is in error message, i.e. position 2
            else
                errPos = 1;  // [Fluid] variable is in error message, i.e. position 1
            // Now determine specific error type
            if (emsg.find("Output string is invalid") != std::string::npos)
                return MAKELRESULT(BAD_PARAMETER, errPos);
            else if (emsg.find("non-trivial") != std::string::npos)
                return MAKELRESULT(NON_TRIVIAL, errPos);
            else if (emsg.find("No outputs") != std::string::npos)
                return MAKELRESULT(NOT_AVAIL, errPos);
            else
                return MAKELRESULT(UNKNOWN, errPos);
        } else
            return MAKELRESULT(UNKNOWN, 1);
    }

    // normal return
    return 0;
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
    unsigned int errPos;
    std::string Prop1Name(InputName1->str);
    std::string Prop2Name(InputName2->str);
    std::string FluidString = FluidName->str;
    CoolProp::parameters key1;

    // check that the first scalar argument is real
    if (InputProp1->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 3);  // if not, display "Argument must be real" under scalar argument

    // check that the second scalar argument is real
    if (InputProp2->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 5);  // if not, display "Argument must be real" under scalar argument

    // pass the arguments to the CoolProp.Props() function
    Prop->real = CoolProp::PropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, FluidName->str);

    // Note: PropsSI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if PropsSI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it
        if (emsg.find("Input pair variable is invalid") != std::string::npos) {
            if (!is_valid_parameter(Prop1Name, key1))
                errPos = 2;  // First input parameter string; position 2.
            else             // must be the second input parameter that's bad
                errPos = 4;  // Second input parameter string; position 4.
            return MAKELRESULT(BAD_PARAMETER, errPos);
        } else if (emsg.find("Input Name1") != std::string::npos) {
            return MAKELRESULT(BAD_PARAMETER, 2);  // first position input parameter
        } else if (emsg.find("Input Name2") != std::string::npos) {
            return MAKELRESULT(BAD_PARAMETER, 4);  // second position input parameter
        } else if (emsg.find("Phase can only be specified on one") != std::string::npos) {
            return MAKELRESULT(ONLY_ONE_PHASE_SPEC, 4);  // second position parameter
        } else if (emsg.find("valid phase") != std::string::npos) {
            if (!is_valid_parameter(Prop1Name, key1))
                errPos = 2;                         // First input parameter string; position 2.
            else                                    // must be the second input parameter that's bad
                errPos = 4;                         // Second input parameter string; position 4.
            return MAKELRESULT(BAD_PHASE, errPos);  // second position parameter
        } else if (emsg.find("This pair of inputs") != std::string::npos) {
            return MAKELRESULT(BAD_INPUT_PAIR, 2);  // second position parameter
        } else if (emsg.find("Input vapor quality") != std::string::npos) {
            if (Prop1Name == "Q")
                return MAKELRESULT(BAD_QUAL, 3);  // First value position
            else
                return MAKELRESULT(BAD_QUAL, 5);  // Second value position
        } else if (emsg.find("Output string is invalid") != std::string::npos) {
            return MAKELRESULT(BAD_PARAMETER, 1);  // first position parameter
        } else if (emsg.find("not valid in two phase region") != std::string::npos) {
            return MAKELRESULT(TWO_PHASE, 1);  // first position parameter
        } else if (emsg.find("only defined within the two-phase") != std::string::npos) {
            return MAKELRESULT(NON_TWO_PHASE, 1);  // first position parameter
        } else if (emsg.find("not implemented") != std::string::npos) {
            return MAKELRESULT(NOT_AVAIL, 1);  // first position parameter
        } else if (emsg.find("Initialize failed") != std::string::npos) {
            if (emsg.find("REFPROP") != std::string::npos) {
                if (emsg.find("cannot use") != std::string::npos)
                    return MAKELRESULT(NO_REFPROP, 6);
                else
                    return MAKELRESULT(BAD_FLUID, 6);
            } else if (emsg.find("IF97") != std::string::npos) {
                return MAKELRESULT(BAD_IF97_FLUID, 6);
            } else
                return MAKELRESULT(BAD_FLUID, 6);
        } else if (emsg.find("Temperature") != std::string::npos) {
            if (Prop1Name == "T")
                return MAKELRESULT(T_OUT_OF_RANGE, 3);  // First value position
            else
                return MAKELRESULT(T_OUT_OF_RANGE, 5);  // Second value position
        } else if (emsg.find("Saturation pressure") != std::string::npos) {
            if (Prop1Name == "P")
                return MAKELRESULT(TP_SATURATION, 3);  // First value position
            else
                return MAKELRESULT(TP_SATURATION, 5);  // Second value position
        } else if (emsg.find("Pressure") != std::string::npos) {
            if (Prop1Name == "P")
                return MAKELRESULT(P_OUT_OF_RANGE, 3);  // First value position
            else
                return MAKELRESULT(P_OUT_OF_RANGE, 5);  // Second value position
        } else if ((emsg.find("Enthalpy") != std::string::npos) || (emsg.find("solution because Hmolar") != std::string::npos)) {
            if ((Prop1Name == "H") || (Prop1Name == "Hmolar"))
                return MAKELRESULT(H_OUT_OF_RANGE, 3);  // First value position
            else
                return MAKELRESULT(H_OUT_OF_RANGE, 5);  // Second value position
        } else if ((emsg.find("Entropy") != std::string::npos) || (emsg.find("solution because Smolar") != std::string::npos)) {
            if ((Prop1Name == "S") || (Prop1Name == "Smolar"))
                return MAKELRESULT(S_OUT_OF_RANGE, 3);  // First value position
            else
                return MAKELRESULT(S_OUT_OF_RANGE, 5);  // Second value position
        } else
            return MAKELRESULT(UNKNOWN, 1);
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
    unsigned int errPos;
    std::string OutName(OutputName->str);
    OutName = "[" + OutName + "]";
    std::string Prop1Name(InputName1->str);
    Prop1Name = "[" + Prop1Name + "]";
    std::string Prop2Name(InputName2->str);
    Prop2Name = "[" + Prop2Name + "]";
    std::string Prop3Name(InputName3->str);
    Prop3Name = "[" + Prop3Name + "]";

    // check that the first scalar argument is real
    if (InputProp1->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 3);  // if not, display "must be real" under scalar argument

    // check that the second scalar argument is real
    if (InputProp2->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 5);  // if not, display "must be real" under scalar argument

    // check that the third scalar argument is real
    if (InputProp3->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 7);  // if not, display "must be real" under scalar argument

    // pass the arguments to the HumidAirProp.HAProps() function
    Prop->real =
      HumidAir::HAPropsSI(OutputName->str, InputName1->str, InputProp1->real, InputName2->str, InputProp2->real, InputName3->str, InputProp3->real);

    // Note: HAPropsSI does not throw exceptions, but instead
    // sets global parameter "errstring" and returns _HUGE.
    // Use ValidNumber(val) to see if HAPropsSI failed with an error message...
    if (!ValidNumber(Prop->real)) {
        std::string emsg = CoolProp::get_global_param_string("errstring");
        CoolProp::set_error_string(emsg);  // reset error string so Mathcad can retrieve it

        errPos = 0;
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
                                        LPCMCSTRING CAS1,  // FIrst component
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

    // Must use MathcadAllocate(size) so Mathcad can track and release
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);  // create a c-string (pointer) c with the same size as s
    // copy s into c, this process avoids the const-cast type which would result from instead
    // converting the string using s.c_str()
    std::copy(s.begin(), s.end(), c);
    c[s.size()] = '\0';
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
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);  // create a c-string (pointer) c with the same size as s
    // copy s into c, this process avoids the const-cast type which would result from instead
    // converting the string using s.c_str()
    std::copy(s.begin(), s.end(), c);
    c[s.size()] = '\0';
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
    if (Value->imag != 0.0) return MAKELRESULT(MUST_BE_REAL, 4);  // if not, display "must be real" under scalar argument

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
    char* c = MathcadAllocate(static_cast<int>(s.size()) + 1);  // create a c-string (pointer) c with the same size as s
    // copy s into c, this process avoids the const-cast type which would result from instead
    // converting the string using s.c_str()
    std::copy(s.begin(), s.end(), c);
    c[s.size()] = '\0';
    // assign the string to the function's output parameter
    Msg->str = c;

    // normal return
    return 0;
}

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO PropsParam = {
  "get_global_param_string",                                 // Name by which MathCAD will recognize the function
  "Name of the parameter to retrieve",                       // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_global_param_string,                   // Pointer to the function code.
  MC_STRING,                                                 // Returns a MathCAD string
  1,                                                         // Number of arguments
  {MC_STRING}                                                // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO FluidParam = {
  "get_fluid_param_string",                                  // Name by which MathCAD will recognize the function
  "Fluid, Name of the parameter to retrieve",                // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_fluid_param_string,                    // Pointer to the function code.
  MC_STRING,                                                 // Returns a MathCAD string
  2,                                                         // Number of arguments
  {MC_STRING, MC_STRING}                                     // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO RefState = {
  "set_reference_state",                                           // Name by which MathCAD will recognize the function
  "Fluid, Reference State String",                                 // Description of input parameters
  "Sets the reference state to either IIR, ASHRAE, NBP, or DEF.",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_set_reference_state,                             // Pointer to the function code.
  COMPLEX_SCALAR,                                                  // Returns a MathCAD complex scalar
  2,                                                               // Number of arguments
  {MC_STRING, MC_STRING}                                           // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO Props1SI = {
  "Props1SI",                                                                                     // Name by which MathCAD will recognize the function
  "Fluid, Property Name",                                                                         // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is not dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_Props1SI,                                                                       // Pointer to the function code.
  COMPLEX_SCALAR,                                                                                 // Returns a MathCAD complex scalar
  2,                                                                                              // Number of arguments
  {MC_STRING, MC_STRING}                                                                          // Argument types
};

// fill out a FUNCTIONINFO structure with the information needed for registering the function with Mathcad
FUNCTIONINFO PropsSI = {
  "PropsSI",                                                                                  // Name by which MathCAD will recognize the function
  "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Fluid Name",  // Description of input parameters
  "Returns a fluid-specific parameter, where the parameter is dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_PropsSI,                                                                    // Pointer to the function code.
  COMPLEX_SCALAR,                                                                             // Returns a MathCAD complex scalar
  6,                                                                                          // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING}                // Argument types
};

FUNCTIONINFO HAPropsSI = {
  "HAPropsSI",  // Name by which MathCAD will recognize the function
  "Output Name, Input Name 1, Input Property 1, Input Name 2, Input Property 2, Input Name 3, Input Property 3",  // Description of input parameters
  "Returns a parameter of humid air, where the parameter is dependent on the fluid state",  // Description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_HAPropsSI,                                                                // Pointer to the function code.
  COMPLEX_SCALAR,                                                                           // Returns a MathCAD complex scalar
  7,                                                                                        // Number of arguments
  {MC_STRING, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR, MC_STRING, COMPLEX_SCALAR}  // Argument types
};

FUNCTIONINFO GetMixtureData = {
  "get_mixture_binary_pair_data",                            // Name by which MathCAD will recognize the function
  "CAS 1, CAS 2, Name of the parameter to retrieve",         // Description of input parameters
  "Returns the value of the requested CoolProps parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_get_mixture_binary_pair_data,              // Pointer to the function code.
  MC_STRING,                                                 // Returns a MathCAD string
  3,                                                         // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING}                          // Argument types
};

FUNCTIONINFO ApplyMixingRule = {
  "apply_simple_mixing_rule",                   // Name by which MathCAD will recognize the function
  "CAS 1, CAS 2, Mixing Rule",                  // Description of input parameters
  "Sets a simple mixing rule for binary pair",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_apply_simple_mixing_rule,     // Pointer to the function code.
  MC_STRING,                                    // Returns a MathCAD string
  3,                                            // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING}             // Argument types
};

FUNCTIONINFO SetMixtureData = {
  "set_mixture_binary_pair_data",                             // Name by which MathCAD will recognize the function
  "CAS 1, CAS 2, Parameter Name, Parameter value",            // Description of input parameters
  "Sets the value of the specified binary mixing parameter",  // description of the function for the Insert Function dialog box
  (LPCFUNCTION)CP_set_mixture_binary_pair_data,               // Pointer to the function code.
  MC_STRING,                                                  // Returns a MathCAD string
  4,                                                          // Number of arguments
  {MC_STRING, MC_STRING, MC_STRING, COMPLEX_SCALAR}           // Argument types
};

// ************************************************************************************
// DLL entry point code.
// ************************************************************************************
// The _CRT_INIT function is needed if you are using Microsoft's 32 bit compiler
#ifdef _WIN32
extern "C" BOOL WINAPI _CRT_INIT(HINSTANCE hinstDLL, DWORD dwReason, LPVOID lpReserved);
#endif

extern "C" BOOL WINAPI DllEntryPoint(HINSTANCE hDLL, DWORD dwReason, LPVOID lpReserved) {
    switch (dwReason) {
        case DLL_PROCESS_ATTACH:
            //
            // DLL is attaching to the address space of
            // the current process.
            //
            if (!_CRT_INIT(hDLL, dwReason, lpReserved)) return FALSE;

            // register the error message table
            // Note, that if your function never returns
            // an error -- you do not need to
            // register an error message table
            if (!CreateUserErrorMessageTable(hDLL, NUMBER_OF_ERRORS, CPErrorMessageTable)) break;

            // and if the errors register OK
            // go ahead and
            // register user function
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
