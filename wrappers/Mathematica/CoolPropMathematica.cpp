
/* Include required header */
#include "CoolProp.h"
#include "CoolPropLib.h"
#include "HumidAirProp.h"
#include "WolframLibrary.h"
#include "Backends/Cubics/CubicsLibrary.h"

/*
    outstring is used for all functions that have a char * output
    local string variables need to be disowned 
    but this cannot happen in the same function for the output function returned to Mathematica
    thus it has to be defined outside of the function
*/

long buffer_length = 20000;
char* outstring = (char*)malloc(sizeof(char) * buffer_length);
long* errcode = (long*)malloc(sizeof(long));
char* message_buffer = (char*)malloc(sizeof(char) * buffer_length);
long resdim1 = 20;
long resdim2 = 20;
long maxComponents = 100;

/* Return the version of Library Link */
extern "C" DLLEXPORT mint WolframLibrary_getVersion() {
    return WolframLibraryVersion;
}

/* Initialize Library */
extern "C" DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
    return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
extern "C" DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
    free(outstring);
    free(errcode);
    free(message_buffer);

    return;
}

/* Adds one to the input, returning the result */
extern "C" DLLEXPORT int plus_one(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;
    double x = MArgument_getReal(Args[0]);
    MArgument_setReal(Res, x + 1.0);
    return LIBRARY_NO_ERROR;
}

/* Change buffer size from Mathematica-side */
extern "C" DLLEXPORT int set_buffer_length(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    *errcode = 0;

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    char* buffer_name = MArgument_getUTF8String(Args[0]);
    mint new_buffer_length = MArgument_getInteger(Args[1]);

    if (new_buffer_length > 0) {
        if (strcmp("stringBuffer", buffer_name) == 0) {
            buffer_length = new_buffer_length;
            free(outstring);
            free(message_buffer);
            outstring = (char*)malloc(sizeof(char) * buffer_length);
            message_buffer = (char*)malloc(sizeof(char) * buffer_length);
        } else if (strcmp("resultMatrix", buffer_name) == 0) {
            resdim1 = new_buffer_length;
            resdim2 = new_buffer_length;
        } else if (strcmp("maxComponents", buffer_name) == 0) {
            maxComponents = new_buffer_length;
        } else {
            *errcode = 1;
            CoolProp::set_error_string(
              format("Buffer name [%s] does not exist. Options are: \"stringBuffer\", \"resultMatrix\" and \"maxComponents\".", buffer_name));
        }
    } else {
        *errcode = 1;
        CoolProp::set_error_string(format("New buffer length [%d] must be bigger than 0.", new_buffer_length));
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(buffer_name);

    return LIBRARY_NO_ERROR;
}

/* Try passing a Matrix back to Mathematica */
extern "C" DLLEXPORT int CreateMatrix(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    mint row, col, rows, cols, count;

    MTensor out;
    mint out_type = MType_Integer;
    mint out_rank = 2;
    mint out_dims[2];
    mint* out_data;
    int err;

    rows = MArgument_getInteger(Args[0]);
    cols = MArgument_getInteger(Args[1]);
    out_dims[0] = rows;
    out_dims[1] = cols;

    err = libData->MTensor_new(out_type, out_rank, out_dims, &out);
    out_data = libData->MTensor_getIntegerData(out);

    count = 1;
    for (row = 0; row < rows; row++) {
        for (col = 0; col < cols; col++) {
            out_data[col + cols * row] = count++;
        }
    }

    MArgument_setMTensor(Res, out);
    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int Props1SI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    char* Output = MArgument_getUTF8String(Args[0]);
    char* FluidName = MArgument_getUTF8String(Args[1]);

    double val = CoolProp::Props1SI(Output, FluidName);

    libData->UTF8String_disown(Output);
    libData->UTF8String_disown(FluidName);

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int Props1SImulti(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 4) return LIBRARY_FUNCTION_ERROR;

    char* Outputs = MArgument_getUTF8String(Args[0]);
    char* backend = MArgument_getUTF8String(Args[1]);
    char* FluidNames = MArgument_getUTF8String(Args[2]);
    MTensor mathematica_fractions = MArgument_getMTensor(Args[3]);

    MTensor out;
    mint out_type = MType_Real;
    mint out_rank = 1;
    mint out_dims[1];
    double* out_data;
    int err;

    double* fractions = libData->MTensor_getRealData(mathematica_fractions);
    long length_fractions = libData->MTensor_getFlattenedLength(mathematica_fractions);

    long* localdim1 = (long*)malloc(sizeof(long));
    *localdim1 = resdim1;

    double* result = (double*)malloc(sizeof(double) * *localdim1);

    Props1SImulti(Outputs, backend, FluidNames, fractions, length_fractions, result, localdim1);

    out_dims[0] = *localdim1;

    err = libData->MTensor_new(out_type, out_rank, out_dims, &out);
    out_data = libData->MTensor_getRealData(out);

    for (int i = 0; i < out_dims[0]; i++) {
        out_data[i] = result[i];
    }

    MArgument_setMTensor(Res, out);
    libData->UTF8String_disown(Outputs);
    libData->UTF8String_disown(FluidNames);
    libData->UTF8String_disown(backend);
    free(localdim1);

    free(result);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int PropsSI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 6) return LIBRARY_FUNCTION_ERROR;

    char* Output = MArgument_getUTF8String(Args[0]);
    char* Name1 = MArgument_getUTF8String(Args[1]);
    double Prop1 = MArgument_getReal(Args[2]);
    char* Name2 = MArgument_getUTF8String(Args[3]);
    double Prop2 = MArgument_getReal(Args[4]);
    char* FluidName = MArgument_getUTF8String(Args[5]);

    double val = CoolProp::PropsSI(Output, Name1, Prop1, Name2, Prop2, FluidName);

    libData->UTF8String_disown(Output);
    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(FluidName);

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int HAPropsSI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 7) return LIBRARY_FUNCTION_ERROR;

    char* Output = MArgument_getUTF8String(Args[0]);
    char* Name1 = MArgument_getUTF8String(Args[1]);
    double Prop1 = MArgument_getReal(Args[2]);
    char* Name2 = MArgument_getUTF8String(Args[3]);
    double Prop2 = MArgument_getReal(Args[4]);
    char* Name3 = MArgument_getUTF8String(Args[5]);
    double Prop3 = MArgument_getReal(Args[6]);

    double val = HumidAir::HAPropsSI(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3);

    libData->UTF8String_disown(Output);
    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(Name3);

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int PropsSImulti(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    if (Argc != 8) return LIBRARY_FUNCTION_ERROR;

    char* Outputs = MArgument_getUTF8String(Args[0]);
    char* Name1 = MArgument_getUTF8String(Args[1]);
    MTensor mathematica_Prop1 = MArgument_getMTensor(Args[2]);
    char* Name2 = MArgument_getUTF8String(Args[3]);
    MTensor mathematica_Prop2 = MArgument_getMTensor(Args[4]);
    char* backend = MArgument_getUTF8String(Args[5]);
    char* FluidNames = MArgument_getUTF8String(Args[6]);
    MTensor mathematica_fractions = MArgument_getMTensor(Args[7]);

    MTensor out;
    mint out_type = MType_Real;
    mint out_rank = 2;
    mint out_dims[2];
    double* out_data;
    int err;

    double* Prop1 = libData->MTensor_getRealData(mathematica_Prop1);
    long size_Prop1 = libData->MTensor_getFlattenedLength(mathematica_Prop1);
    double* Prop2 = libData->MTensor_getRealData(mathematica_Prop2);
    long size_Prop2 = libData->MTensor_getFlattenedLength(mathematica_Prop2);
    double* fractions = libData->MTensor_getRealData(mathematica_fractions);
    long length_fractions = libData->MTensor_getFlattenedLength(mathematica_fractions);

    long* localdim1 = (long*)malloc(sizeof(long));
    long* localdim2 = (long*)malloc(sizeof(long));
    *localdim1 = resdim1;
    *localdim2 = resdim2;

    double* result = (double*)malloc(sizeof(double) * *localdim1 * *localdim2);

    PropsSImulti(Outputs, Name1, Prop1, size_Prop1, Name2, Prop2, size_Prop2, backend, FluidNames, fractions, length_fractions, result, localdim1,
                 localdim2);

    out_dims[0] = *localdim1;
    out_dims[1] = *localdim2;

    err = libData->MTensor_new(out_type, out_rank, out_dims, &out);
    out_data = libData->MTensor_getRealData(out);

    for (int i = 0; i < out_dims[0]; i++) {
        for (int j = 0; j < out_dims[1]; j++) {
            out_data[j + out_dims[1] * i] = result[j + out_dims[1] * i];
        }
    }

    MArgument_setMTensor(Res, out);
    libData->UTF8String_disown(Outputs);
    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(FluidNames);
    libData->UTF8String_disown(backend);

    free(result);
    free(localdim1);
    free(localdim2);

    return LIBRARY_NO_ERROR;
}

extern "C" DLLEXPORT int PhaseSI(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(outstring, "");

    if (Argc != 5) return LIBRARY_FUNCTION_ERROR;

    char* Name1 = MArgument_getUTF8String(Args[0]);
    double Prop1 = MArgument_getReal(Args[1]);
    char* Name2 = MArgument_getUTF8String(Args[2]);
    double Prop2 = MArgument_getReal(Args[3]);
    char* FluidName = MArgument_getUTF8String(Args[4]);

    std::string val = CoolProp::PhaseSI(Name1, Prop1, Name2, Prop2, FluidName);

    strcpy(outstring, val.c_str());

    MArgument_setUTF8String(Res, outstring);

    libData->UTF8String_disown(Name1);
    libData->UTF8String_disown(Name2);
    libData->UTF8String_disown(FluidName);

    return LIBRARY_NO_ERROR;
}

// # ---------------------------------
// #        Information functions
// # ---------------------------------

/*
    get_global_param_string( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get a globally-defined string.

# Ref
ref get_global_param_string

# Arguments
* `key`: A string represents parameter name
"version"
"gitrevision"
"errstring"
"warnstring"
"FluidsList"
"incompressible_list_pure"
"incompressible_list_solution"
"mixture_binary_pairs_list"
"parameter_list"
"predefined_mixtures"
"HOME"
"cubic_fluids_schema"
*/

extern "C" DLLEXPORT int get_global_param_string(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    // to catch faulty input to get_global_param_string, the output has to be different from ""
    // as "" might also be the output for "errstring" or "warnstring"
    strcpy(outstring, "no string");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    char* param = MArgument_getUTF8String(Args[0]);

    get_global_param_string(param, outstring, buffer_length);

    MArgument_setUTF8String(Res, outstring);

    libData->UTF8String_disown(param);

    return LIBRARY_NO_ERROR;
}

/*
    get_input_pair_index( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get the index for an input pair "PT_INPUTS", "HmolarQ_INPUTS", etc, for `abstractstate_update()` function.

# Arguments
* `pair::AbstractString`: The name of an input pair, described in the following table

Input Pair          |Description
:-------------------|:------------------------------------------------
QT_INPUTS           |Molar quality, Temperature in K
QSmolar_INPUTS      |Molar quality, Entropy in J/mol/K
QSmass_INPUTS       |Molar quality, Entropy in J/kg/K
HmolarQ_INPUTS      |Enthalpy in J/mol, Molar quality
HmassQ_INPUTS       |Enthalpy in J/kg, Molar quality
DmassQ_INPUTS       |Molar density kg/m^3, Molar quality
DmolarQ_INPUTS      |Molar density in mol/m^3, Molar quality
PQ_INPUTS           |Pressure in Pa, Molar quality
PT_INPUTS           |Pressure in Pa, Temperature in K
DmassT_INPUTS       |Mass density in kg/m^3, Temperature in K
DmolarT_INPUTS      |Molar density in mol/m^3, Temperature in K
HmassT_INPUTS       |Enthalpy in J/kg, Temperature in K
HmolarT_INPUTS      |Enthalpy in J/mol, Temperature in K
SmassT_INPUTS       |Entropy in J/kg/K, Temperature in K
SmolarT_INPUTS      |Entropy in J/mol/K, Temperature in K
TUmass_INPUTS       |Temperature in K, Internal energy in J/kg
TUmolar_INPUTS      |Temperature in K, Internal energy in J/mol
DmassP_INPUTS       |Mass density in kg/m^3, Pressure in Pa
DmolarP_INPUTS      |Molar density in mol/m^3, Pressure in Pa
HmassP_INPUTS       |Enthalpy in J/kg, Pressure in Pa
HmolarP_INPUTS      |Enthalpy in J/mol, Pressure in Pa
PSmass_INPUTS       |Pressure in Pa, Entropy in J/kg/K
PSmolar_INPUTS      |Pressure in Pa, Entropy in J/mol/K
PUmass_INPUTS       |Pressure in Pa, Internal energy in J/kg
PUmolar_INPUTS      |Pressure in Pa, Internal energy in J/mol
DmassHmass_INPUTS   |Mass density in kg/m^3, Enthalpy in J/kg
DmolarHmolar_INPUTS |Molar density in mol/m^3, Enthalpy in J/mol
DmassSmass_INPUTS   |Mass density in kg/m^3, Entropy in J/kg/K
DmolarSmolar_INPUTS |Molar density in mol/m^3, Entropy in J/mol/K
DmassUmass_INPUTS   |Mass density in kg/m^3, Internal energy in J/kg
DmolarUmolar_INPUTS |Molar density in mol/m^3, Internal energy in J/mol
HmassSmass_INPUTS   |Enthalpy in J/kg, Entropy in J/kg/K
HmolarSmolar_INPUTS |Enthalpy in J/mol, Entropy in J/mol/K
SmassUmass_INPUTS   |Entropy in J/kg/K, Internal energy in J/kg
SmolarUmolar_INPUTS |Entropy in J/mol/K, Internal energy in J/mol
*/

extern "C" DLLEXPORT int get_input_pair_index(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    char* pair = MArgument_getUTF8String(Args[0]);

    long val = get_input_pair_index(pair);

    MArgument_setInteger(Res, val);

    libData->UTF8String_disown(pair);

    return LIBRARY_NO_ERROR;
}

/*
    get_input_pair_description

Get the information about an input pair for the given index.

# Arguments
* `index`: A integer index represents parameter
*/

extern "C" DLLEXPORT int get_input_pair_description(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(outstring, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long index = MArgument_getInteger(Args[0]);
    char* outtype = MArgument_getUTF8String(Args[1]);

    std::string val = "";

    try {
        // 35 is the number of elements in input_pairs, has to be updated, should input pairs be added
        // 0 is the default/invalid value
        if (index > 0 && index <= 35) {
            if (strcmp("short", outtype) == 0)
                val = CoolProp::get_input_pair_short_desc(static_cast<CoolProp::input_pairs>(index));
            else if (strcmp("long", outtype) == 0)
                val = CoolProp::get_input_pair_long_desc(static_cast<CoolProp::input_pairs>(index));
            else if (strcmp("both", outtype) == 0)
                val = CoolProp::get_input_pair_short_desc(static_cast<CoolProp::input_pairs>(index)) + ": "
                      + CoolProp::get_input_pair_long_desc(static_cast<CoolProp::input_pairs>(index));
            else
                throw CoolProp::ValueError(
                  format("Bad info string [%s] to get_parameter_information, options are \"short\", \"long\" or \"both\" (default)", outtype));
            strcpy(outstring, val.c_str());
        } else {
            throw CoolProp::ValueError(
              format("Unable to match the key [%d] in get_parameter_information for info [%s], range is 1 to 35.", index, outtype));
        }
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }

    MArgument_setUTF8String(Res, outstring);

    libData->UTF8String_disown(outtype);

    return LIBRARY_NO_ERROR;
}

/*
    get_param_index( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get the index as a long for a parameter "T", "P", etc, for `abstractstate_keyed_output()` function.

# Arguments
* `param`: A string represents parameter name, to see full list check "Table of string inputs to PropsSI function": http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table, or simply type `get_global_param_string("parameter_list")`
*/

extern "C" DLLEXPORT int get_param_index(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    char* param = MArgument_getUTF8String(Args[0]);

    long val = get_param_index(param);

    MArgument_setInteger(Res, val);

    libData->UTF8String_disown(param);

    return LIBRARY_NO_ERROR;
}

/*
    get_parameter_information

Get the information about a parameter for the given index.

# Arguments
* `index`: A integer index represents parameter
* `outtype="long"`: Output type, could be one of the `["IO", "short", "long", "units"]`, with a default value of "long"
*/

extern "C" DLLEXPORT int get_parameter_information(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(outstring, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long index = MArgument_getInteger(Args[0]);
    char* outtype = MArgument_getUTF8String(Args[1]);

    // std::string val = "empty string";

    try {
        std::string val = CoolProp::get_parameter_information(index, outtype);
        strcpy(outstring, val.c_str());
    } catch (std::exception& e) {
        CoolProp::set_error_string(e.what());
    } catch (...) {
        CoolProp::set_error_string("Undefined error");
    }

    MArgument_setUTF8String(Res, outstring);

    libData->UTF8String_disown(outtype);

    return LIBRARY_NO_ERROR;
}

/*
    add_fluids_as_JSON( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Add fluids as a JSON-formatted string

# Arguments
`backend` The backend to which these should be added; e.g. "HEOS", "SRK", "PR"
`fluidstring` The JSON-formatted string
*/

extern "C" DLLEXPORT int add_fluids_as_JSON(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    char* backend = MArgument_getUTF8String(Args[0]);
    char* fluidstring = MArgument_getUTF8String(Args[1]);

    add_fluids_as_JSON(backend, fluidstring, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(backend);
    libData->UTF8String_disown(fluidstring);

    return LIBRARY_NO_ERROR;
}

/*
    get_fluid_param_string( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get a string for a value from a fluid (numerical values for the fluid can be obtained from Props1SI function) 

# Arguments
`fluid` The name of the fluid that is part of CoolProp, for instance "n-Propane" 
`param`     	A string, can be in one of the terms described in the following table

ParamName 	Description
"aliases" 	A comma separated list of aliases for the fluid
"CAS", "CAS_number" 	The CAS number
"ASHRAE34" 	The ASHRAE standard 34 safety rating
"REFPROPName","REFPROP_name" 	The name of the fluid used in REFPROP
"Bibtex-XXX" 	A BibTeX key, where XXX is one of the bibtex keys used in get_BibTeXKey
"pure" 	"true" if the fluid is pure, "false" otherwise
"formula" 	The chemical formula of the fluid in LaTeX form if available, "" otherwise
*/

extern "C" DLLEXPORT int get_fluid_param_string(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {
    strcpy(outstring, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    char* fluid = MArgument_getUTF8String(Args[0]);
    char* param = MArgument_getUTF8String(Args[1]);

    long result = get_fluid_param_string(fluid, param, outstring, buffer_length);

    if (result == 0) {
        strcpy(outstring, "");
    }

    MArgument_setUTF8String(Res, outstring);

    libData->UTF8String_disown(fluid);
    libData->UTF8String_disown(param);

    return LIBRARY_NO_ERROR;
}

// # ---------------------------------
// #        AbstractState
// # ---------------------------------

/*
    AbstractState_factory( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Generate an AbstractState instance return an integer handle to the state class generated to be used in the other low-level accessor functions.

# Arguments
* `backend`: The backend you will use could be: `["HEOS", "REFPROP", "INCOMP", "IF97", "TREND", "HEOS&TTSE", "HEOS&BICUBIC", "SRK", "PR", "VTPR"]` etc.
* `fluids`: '&' delimited list of fluids. To get a list of possible values call `get_global_param_string(key)` with `key` one of the following: `["FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "predefined_mixtures"]`, also there is a list in CoolProp online documentation [List of Fluids](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids), or simply type `?CoolProp_fluids`
*/

extern "C" DLLEXPORT int AbstractState_factory(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    char* backend = MArgument_getUTF8String(Args[0]);
    char* fluids = MArgument_getUTF8String(Args[1]);

    long val = AbstractState_factory(backend, fluids, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, val);

    libData->UTF8String_disown(backend);
    libData->UTF8String_disown(fluids);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_free( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Release a state class generated by the low-level interface wrapper.

# Arguments
* `handle`: The integer handle for the state class stored in memory
*/

extern "C" DLLEXPORT int AbstractState_free(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);

    AbstractState_free(handle, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_update( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Update the state of the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `input_pair::AbstractString`: The name of an input pair
* `value1`: The first input value
* `value2`: The second input value
*/

extern "C" DLLEXPORT int AbstractState_update(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 4) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    long input_pair = MArgument_getInteger(Args[1]);
    double value1 = MArgument_getReal(Args[2]);
    double value2 = MArgument_getReal(Args[3]);

    AbstractState_update(handle, input_pair, value1, value2, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_keyed_output( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get an output value from the `AbstractState` using an integer value for the desired output value.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `param::Clong`: param The integer value for the parameter you want

*/

extern "C" DLLEXPORT int AbstractState_keyed_output(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    long param = MArgument_getInteger(Args[1]);

    double val = AbstractState_keyed_output(handle, param, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setReal(Res, val);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_set_fractions( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Set the fractions (mole, mass, volume) for the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `fractions`: The array of fractions
*/

extern "C" DLLEXPORT int AbstractState_set_fractions(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    MTensor mathematica_fractions = MArgument_getMTensor(Args[1]);

    double* fractions = libData->MTensor_getRealData(mathematica_fractions);
    int components_number = libData->MTensor_getFlattenedLength(mathematica_fractions);

    AbstractState_set_fractions(handle, fractions, components_number, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_set_binary_interaction_double( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Set binary interraction parameter for different mixtures model e.g.: "linear", "Lorentz-Berthelot"

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `componentI`: indice of the first fluid of the binary pair
* `componentJ`: indice of the second fluid of the binary pair
* `parameter`: string wit the name of the parameter, e.g.: "betaT", "gammaT", "betaV", "gammaV"
* `value`: the value of the binary interaction parameter
*/

extern "C" DLLEXPORT int AbstractState_set_binary_interaction_double(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 5) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    long componentI = MArgument_getInteger(Args[1]);
    long componentJ = MArgument_getInteger(Args[2]);
    char* parameter = MArgument_getUTF8String(Args[3]);
    double value = MArgument_getReal(Args[4]);

    AbstractState_set_binary_interaction_double(handle, componentI, componentJ, parameter, value, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(parameter);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_set_cubic_alpha_C( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Set cubic's alpha function parameters

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `componentI` indice of the fluid the parramter should be applied too (for mixtures)
* `param` parameter the string specifying the alpha function to use, ex "TWU" for the TWU alpha function
* `c` three parameters for the alpha function
*/

extern "C" DLLEXPORT int AbstractState_set_cubic_alpha_C(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 4) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    long componentI = MArgument_getInteger(Args[1]);
    char* param = MArgument_getUTF8String(Args[2]);
    MTensor mathematica_parameters = MArgument_getMTensor(Args[3]);

    double* parameters = libData->MTensor_getRealData(mathematica_parameters);

    AbstractState_set_cubic_alpha_C(handle, componentI, param, parameters[0], parameters[1], parameters[2], errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(param);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_get_mole_fractions( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get the molar fractions for the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
*/

extern "C" DLLEXPORT int AbstractState_get_mole_fractions(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);

    MTensor T0;
    mint i;
    int err = LIBRARY_NO_ERROR;

    double* fractions = (double*)malloc(sizeof(double) * maxComponents);
    long* components_number = (long*)malloc(sizeof(long));
    mint* dim = (mint*)malloc(sizeof(mint));

    *components_number = 0;

    AbstractState_get_mole_fractions(handle, fractions, maxComponents, components_number, errcode, message_buffer, buffer_length);

    if (*errcode == 0) {
        *dim = *components_number;
    } else {
        CoolProp::set_error_string(message_buffer);
        *dim = 0;
    }

    if (*dim == 0) {
        i = 1;
        *dim = 1;
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        err = libData->MTensor_setReal(T0, &i, 0);
    } else {
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        for (i = 1; i <= *dim && !err; i++) {
            err = libData->MTensor_setReal(T0, &i, fractions[i - 1]);
        }
    }

    MArgument_setMTensor(Res, T0);

    free(dim);
    free(components_number);
    free(fractions);

    return err;
}

// alternate version of AbstractState_get_mole_fractions
extern "C" DLLEXPORT int AbstractState_get_mole_fractions_ALTERNATE(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);

    MTensor T0;
    mint i;
    int err = LIBRARY_NO_ERROR;

    double* fractions = (double*)malloc(sizeof(double) * maxComponents);
    long* components_number = (long*)malloc(sizeof(long));
    mint* dim = (mint*)malloc(sizeof(mint));

    *components_number = 0;

    double* res;

    AbstractState_get_mole_fractions(handle, fractions, maxComponents, components_number, errcode, message_buffer, buffer_length);

    if (*components_number <= maxComponents) {
        *dim = *components_number;
    } else {
        *dim = 0;
    }

    if (*dim == 0) {
        i = 1;
        *dim = 1;
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        err = libData->MTensor_setReal(T0, &i, 0);
    } else {
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        res = libData->MTensor_getRealData(T0);
        for (i = 0; i < *dim && !err; i++) {
            res[i] = fractions[i];
        }
    }

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setMTensor(Res, T0);

    free(dim);
    free(components_number);
    free(fractions);

    return err;
}

/*
    AbstractState_fluid_names( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get the fluid names for the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
*/

extern "C" DLLEXPORT int AbstractState_fluid_names(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");
    strcpy(outstring, "");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);

    mint i;

    long* components_number = (long*)malloc(sizeof(long));

    *components_number = 0;

    AbstractState_fluid_names(handle, outstring, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setUTF8String(Res, outstring);

    free(components_number);

    return LIBRARY_NO_ERROR;
}

/*
    get_cubic_fluids_list( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Returns all fluids implemented in the CubicsLibrary.
*/

extern "C" DLLEXPORT int get_cubic_fluids_list(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    // to catch faulty input to get_global_param_string, the output has to be different from ""
    // as "" might also be the output for "errstring" or "warnstring"
    strcpy(outstring, "no string");

    if (Argc != 0) return LIBRARY_FUNCTION_ERROR;

    std::string val = CoolProp::CubicLibrary::get_cubic_fluids_list();

    strcpy(outstring, val.c_str());

    MArgument_setUTF8String(Res, outstring);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_specify_phase( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Specify the phase to be used for all further calculations.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `phase`: The string with the phase to use. Possible candidates are listed below:

Phase name                 |Condition
:--------------------------|-----------------------
phase_liquid               |
phase_gas                  |
phase_twophase             |
phase_supercritical        |
phase_supercritical_gas    |p < pc, T > Tc
phase_supercritical_liquid |p > pc, T < Tc
phase_critical_point       |p = pc, T = Tc
phase_unknown              |
phase_not_imposed          |
*/

extern "C" DLLEXPORT int AbstractState_specify_phase(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    char* phase = MArgument_getUTF8String(Args[1]);

    AbstractState_specify_phase(handle, phase, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(phase);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_build_phase_envelope( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Build the phase envelope.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `level`: How much refining of the phase envelope ("none" to skip refining (recommended) or "veryfine")

*/

extern "C" DLLEXPORT int AbstractState_build_phase_envelope(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    char* level = MArgument_getUTF8String(Args[1]);

    AbstractState_build_phase_envelope(handle, level, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setInteger(Res, *errcode);

    libData->UTF8String_disown(level);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_get_phase_envelope_data(handle::Clong, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar_vap::Array{Float64}, rhomolar_liq::Array{Float64}, x::Array{Float64}, y::Array{Float64})

Get data from the phase envelope for the given mixture composition.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `length`: The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
* `T`: The pointer to the array of temperature (K)
* `p`: The pointer to the array of pressure (Pa)
* `rhomolar_vap`: The pointer to the array of molar density for vapor phase (m^3/mol)
* `rhomolar_liq`: The pointer to the array of molar density for liquid phase (m^3/mol)
* `x`: The compositions of the "liquid" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)
* `y`: The compositions of the "vapor" phase (WARNING: buffer should be Ncomp*Npoints in length, at a minimum, but there is no way to check buffer length at runtime)

*/

extern "C" DLLEXPORT int AbstractState_get_phase_envelope_data(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 3) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    long length = MArgument_getInteger(Args[1]);
    long maxComponents = MArgument_getInteger(Args[2]);

    MTensor T0;
    mint i;
    mint pos[2];
    int err = LIBRARY_NO_ERROR;

    double* T = (double*)malloc(sizeof(double) * length);
    double* p = (double*)malloc(sizeof(double) * length);
    double* rhomolar_vap = (double*)malloc(sizeof(double) * length);
    double* rhomolar_liq = (double*)malloc(sizeof(double) * length);
    double* x = (double*)malloc(sizeof(double) * length * maxComponents);
    double* y = (double*)malloc(sizeof(double) * length * maxComponents);

    long* actual_length = (long*)malloc(sizeof(long));
    mint dims[2];
    double* res;
    long* actual_components = (long*)malloc(sizeof(long));

    *actual_length = 0;
    *actual_components = 0;

    AbstractState_get_phase_envelope_data(handle, length, maxComponents, T, p, rhomolar_vap, rhomolar_liq, x, y, actual_length, actual_components,
                                          errcode, message_buffer, buffer_length);

    if (*errcode == 0 && *actual_length > 0) {
        dims[0] = *actual_length;
        dims[1] = 4 + 2 * *actual_components;
    } else {
        CoolProp::set_error_string(message_buffer);
        dims[0] = 1;
        dims[1] = 0;
    }

    if (dims[1] == 0) {
        err = libData->MTensor_new(MType_Real, 1, &dims[0], &T0);
        // err = libData->MTensor_setReal( T0, pos, 0);
        res = libData->MTensor_getRealData(T0);
        res[0] = 0.;
    } else {
        err = libData->MTensor_new(MType_Real, 2, dims, &T0);
        res = libData->MTensor_getRealData(T0);

        for (i = 0; i < dims[0] && !err; i++) {
            res[0 + dims[1] * i] = T[i];
            res[1 + dims[1] * i] = p[i];
            res[2 + dims[1] * i] = rhomolar_vap[i];
            res[3 + dims[1] * i] = rhomolar_liq[i];

            for (int j = 0; j < *actual_components; j++) {
                res[4 + j + dims[1] * i] = x[i * *actual_components + j];
                res[4 + *actual_components + j + dims[1] * i] = y[i * *actual_components + j];
            }
        }
    }

    MArgument_setMTensor(Res, T0);

    free(T);
    free(p);
    free(rhomolar_vap);
    free(rhomolar_liq);
    free(x);
    free(y);
    free(actual_length);

    return err;
}

/*
    AbstractState_keyed_output_satState( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) 

Get an output value from the AbstractState using an integer value for the desired output value and desired saturated State

# Arguments
* `handle` The integer handle for the state class stored in memory
* `saturated_state` The string specifying the state (liquid or gas)
* `param` The integer value for the parameter you want
 */

extern "C" DLLEXPORT int AbstractState_keyed_output_satState(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 3) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    char* satState = MArgument_getUTF8String(Args[1]);
    long param = MArgument_getInteger(Args[2]);

    double val = AbstractState_keyed_output_satState(handle, satState, param, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setReal(Res, val);

    libData->UTF8String_disown(satState);

    return LIBRARY_NO_ERROR;
}

/*
    AbstractState_get_mole_fractions_satState( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Get the molar fractions for the AbstractState and the desired saturated State

# Arguments
* `handle` The integer handle for the state class stored in memory
* `saturated_state` The string specifying the state (liquid or gas)
* `fractions` The array of fractions
* `maxN` The length of the buffer for the fractions
* `N` number of fluids
 */
extern "C" DLLEXPORT int AbstractState_get_mole_fractions_satState(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");

    if (Argc != 2) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);
    char* satState = MArgument_getUTF8String(Args[1]);

    MTensor T0;
    mint i;
    int err = LIBRARY_NO_ERROR;

    double* fractions = (double*)malloc(sizeof(double) * maxComponents);
    long* components_number = (long*)malloc(sizeof(long));
    mint* dim = (mint*)malloc(sizeof(mint));

    *components_number = 0;

    double* res;

    AbstractState_get_mole_fractions_satState(handle, satState, fractions, maxComponents, components_number, errcode, message_buffer, buffer_length);

    if (*components_number <= maxComponents) {
        *dim = *components_number;
    } else {
        *dim = 0;
    }

    if (*dim == 0) {
        i = 1;
        *dim = 1;
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        err = libData->MTensor_setReal(T0, &i, 0);
    } else {
        err = libData->MTensor_new(MType_Real, 1, dim, &T0);
        res = libData->MTensor_getRealData(T0);
        for (i = 0; i < *dim && !err; i++) {
            res[i] = fractions[i];
        }
    }

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setMTensor(Res, T0);

    free(dim);
    free(components_number);
    free(fractions);
    libData->UTF8String_disown(satState);

    return err;
}

/*
    AbstractState_backend_name( WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)

Return the name of the backend used in the AbstractState

# Arguments
* `handle` The integer handle for the state class stored in memory
* `backend` The char pointer the name is written to
*/
extern "C" DLLEXPORT int AbstractState_backend_name(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res) {

    strcpy(message_buffer, "");
    strcpy(outstring, "");

    if (Argc != 1) return LIBRARY_FUNCTION_ERROR;

    long handle = MArgument_getInteger(Args[0]);

    AbstractState_backend_name(handle, outstring, errcode, message_buffer, buffer_length);

    if (*errcode != 0) {
        CoolProp::set_error_string(message_buffer);
    }

    MArgument_setUTF8String(Res, outstring);

    return LIBRARY_NO_ERROR;
}