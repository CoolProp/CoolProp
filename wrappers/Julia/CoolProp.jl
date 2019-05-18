#__precompile__()
module CoolProp
using Compat

errcode = Ref{Clong}(0)
const buffer_length = 20000
message_buffer = Array{UInt8}(undef, buffer_length)

const inputs_to_get_global_param_string = ["version", "gitrevision", "errstring", "warnstring", "FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "parameter_list", "predefined_mixtures", "HOME", "cubic_fluids_schema"]

# ---------------------------------
#        High-level functions
# ---------------------------------

"""
    PropsSI(fluid::AbstractString, output::AbstractString)

Return a value that does not depend on the thermodynamic state - this is a convenience function that does the call `PropsSI(output, "", 0, "", 0, fluid)`.

# Arguments
* `fluid::AbstractString`: The name of the fluid that is part of CoolProp, for instance "n-Propane", to get a list of possible values types call `get_global_param_string(key)` with `key` one of the following: `["FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "predefined_mixtures"]`, also there is a list in CoolProp online documentation [List of Fluids](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids), or simply type `?CoolProp_fluids`
* `output::AbstractString`: The name of parameter to evaluate. to see a list type `?CoolProp_parameters`

# Example
```julia
julia> PropsSI("n-Butane", "rhomolar_critical")
3922.769612987809
```

# Ref
CoolProp::Props1SI(std::string, std::string)
"""
function PropsSI(fluid::AbstractString, output::AbstractString)
    val = ccall( (:Props1SI, "CoolProp"), Cdouble, (Cstring, Cstring), fluid, output)
    if val == Inf
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

"""
    PropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

Return a value that depends on the thermodynamic state.
> For pure and pseudo-pure fluids, two state points are required to fix the state. The equations of state are based on T and ρ as state variables, so T, ρ will always be the fastest inputs. P, T will be a bit slower (3-10 times), and then comes inputs where neither T nor ρ are given, like p, h. They will be much slower. If speed is an issue, you can look into table-based interpolation methods using TTSE or bicubic interpolation.

# Arguments
* `fluid::AbstractString`: The name of the fluid that is part of CoolProp, for instance "n-Propane", to get a list of different passible fulid types call `get_global_param_string(key)` with `key` one of the following: `["FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "predefined_mixtures"]`, also there is a list in CoolProp online documentation [List of Fluids](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids)
* `output::AbstractString`: The name of parameter to evaluate. to see a list type `?CoolProp_parameters`
* `name1::AbstractString`: The name of parameter for first state point
* `value1::Real`: Value of the first state point
* `name2::AbstractString`: The name of parameter for second state point
* `value2::Real`: Value of the second state point

# Example
```julia
julia> PropsSI("D", "T", 300, "P", 101325, "n-Butane")
2.4325863624814326
julia> PropsSI("D", "T", 300, "P", 101325, "INCOMP::DEB") # incompressible pure
857.1454
julia> PropsSI("D", "T", 300, "P", 101325, "INCOMP::LiBr[0.23]") # incompressible mass-based binary mixture
1187.5438243617214
julia> PropsSI("D", "T", 300, "P", 101325, "INCOMP::ZM[0.23]") # incompressible volume-based binary mixtures
1028.7273860290911
julia> PropsSI("Dmass", "T", 300, "P", 101325, "R125[0.5]&R32[0.5]")
3.5413381483914512
julia> split(get_global_param_string("mixture_binary_pairs_list"), ', ')[1] # a random binary pair
"100-41-4&106-42-3"
julia> PropsSI("Dmass", "T", 300, "P", 101325, "100-41-4[0.5]&106-42-3[0.5]") # ethylbenzene[0.5]&p-Xylene[0.5]
857.7381127561846
```

# Ref
CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
"""
function PropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)
    val = ccall( (:PropsSI, "CoolProp"), Cdouble, (Cstring, Cstring, Cdouble, Cstring, Cdouble, Cstring), output, name1, value1, name2, value2, fluid)
    if val == Inf
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

"""
    PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)

Return a string representation of the phase. Valid states are: "liquid", "supercritical", "supercritical_gas", "supercritical_liquid", "critical_point", "gas", "twophase"

# Arguments
* `fluid::AbstractString`: The name of the fluid that is part of CoolProp, for instance "n-Propane", to get a list of different passible fulid types call `get_global_param_string(key)` with `key` one of the following: `["FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "predefined_mixtures"]`, also there is a list in CoolProp online documentation [List of Fluids](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids)
* `name1::AbstractString`: The name of parameter for first state point
* `value1::Real`: Value of the first state point
* `name2::AbstractString`: The name of parameter for second state point
* `value2::Real`: Value of the second state point

# Example
```julia
julia> PhaseSI("T", PropsSI("TCRIT", "Water"), "P", PropsSI("PCRIT", "Water"), "Water")
"critical_point"
julia> PhaseSI("T", 300, "Q", 1, "Water")
"twophase"
julia> PhaseSI("P", "T", 300, "Q", 1, "Water")
3536.806750422325
julia> PhaseSI("T", 300, "P", 3531, "Water")
"gas"
julia> PhaseSI("T", 300, "P", 3541, "Water")
"liquid"
```

# Ref
CoolProp::PhaseSI(const std::string &, double, const std::string &, double, const std::string&)
"""
function PhaseSI(name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, fluid::AbstractString)
    val = ccall( (:PhaseSI, "CoolProp"), Int32, (Cstring, Cdouble, Cstring, Cdouble, Cstring, Ptr{UInt8}, Int), name1, value1, name2, value2, fluid, message_buffer::Array{UInt8, 1}, buffer_length)
    val = unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
    if val == ""
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

#=No sample no test
"""
    set_departure_functions(string_data::AbstractString)

Set the departure functions in the departure function library from a string format

# Arguments
* `string_data`: The departure functions to be set, either provided as a JSON-formatted string or as a string of the contents of a HMX.BNC file from REFPROP

# Note
By default, if a departure function already exists in the library, this is an error, unless the configuration variable OVERWRITE_DEPARTURE_FUNCTIONS is set to true

# Ref
CoolProp::set_departure_functions(const char * string_data, long *errcode, char *message_buffer, const long buffer_length);
"""
function set_departure_functions(string_data::AbstractString)
  errcode = ref{Clong}(0);
  ccall( (:set_departure_functions, "CoolProp"), Nothing, (Cstring, Ptr{Clong}, Ptr{UInt8}, Int), string_data, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
  if errcode != 0
    error("CoolProp: ", unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1}))))
  end
end
=#

"""
    set_reference_state(ref::AbstractString, reference_state::AbstractString)

Set the reference state based on a string representation.

#Arguments
* `fluid::AbstractString`	The name of the fluid (Backend can be provided like "REFPROP::Water", or if no backend is provided, "HEOS" is the assumed backend)
* `reference_state::AbstractString`	The reference state to use, one of:

Reference State |	Description
:---------------|:-------------------------------------
"IIR"	          |h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
"ASHRAE"        |h = 0, s = 0 @ -40C saturated liquid
"NBP"	          |h = 0, s = 0 @ 1.0 bar saturated liquid
"DEF"	          |Reset to the default reference state for the fluid
"RESET"	        |Remove the offset

#Example
```julia
julia> h0=-15870000.0; # J/kg
julia> s0= 3887.0; #J/kg
julia> rho0=997.1;
julia> T0=298.15;
julia> M = PropsSI("molemass", "Water");
julia> set_reference_stateD("Water", T0, rho0/M, h0*M, s0*M);
julia> PropsSI("H", "T", T0, "P", 101325, "Water")
-1.5870107493843542e7
julia> set_reference_stateS("Water", "DEF");
julia> PropsSI("H", "T", T0, "P", 101325, "Water")
104920.1198093371
```

#Ref
set_reference_stateS(const std::string& FluidName, const std::string& reference_state)

#Note
The changing of the reference state should be part of the initialization of your program, and it is not recommended to change the reference state during the course of making calculations
"""
function set_reference_state(fluid::AbstractString, reference_state::AbstractString)
    val = ccall( (:set_reference_stateS, "CoolProp"), Cint, (Cstring, Cstring), fluid, reference_state)
    if val == 0
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

"""
    set_reference_state(fluid::AbstractString, T::Real, rhomolar::Real, hmolar0::Real, smolar0::Real)

Set the reference state based on a thermodynamic state point specified by temperature and molar density.

#Arguments
* `fluid::AbstractString`	The name of the fluid
* `T::Real`	Temperature at reference state [K]
* `rhomolar::Real`	Molar density at reference state [mol/m^3]
* `hmolar0::Real`	Molar enthalpy at reference state [J/mol]
* `smolar0::Real`	Molar entropy at reference state [J/mol/K]

#Ref
set_reference_stateD(const char* Ref, double T, double rhomolar, double hmolar0, double smolar0)
"""
function set_reference_state(fluid::AbstractString, T::Real, rhomolar::Real, hmolar0::Real, smolar0::Real)
    val = ccall( (:set_reference_stateD, "CoolProp"), Cint, (Cstring, Cdouble, Cdouble, Cdouble, Cdouble), fluid, T, rhomolar, hmolar0, smolar0)
    if val == 0
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

"""
    saturation_ancillary(fluid_name::AbstractString, output::AbstractString, quality::Integer, input::AbstractString, value::Real)

Extract a value from the saturation ancillary.

# Arguments
* `fluid_name`: The name of the fluid to be used - HelmholtzEOS backend only
* `output`: The desired output variable ("P" for instance for pressure)
* `quality`: The quality, 0 or 1
* `input`: The input variable ("T")
* `value`: The input value

#  Example
julia> saturation_ancillary("R410A","I",1,"T", 300)
0.004877519938463293

# Ref
double saturation_ancillary(const char* fluid_name, const char* output, int Q, const char* input, double value);
"""
function saturation_ancillary(fluid_name::AbstractString, output::AbstractString, quality::Integer, input::AbstractString, value::Real)
    val = ccall( (:saturation_ancillary, "CoolProp"), Cdouble, (Cstring, Cstring, Cint, Cstring, Cdouble), fluid_name, output, quality, input, value)
    if val == Inf
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

# ---------------------------------
#      Configuration functions
# ---------------------------------

"""
    set_config(key::AbstractString, val::AbstractString)

Set configuration string.

# Arguments
* `key::AbstractString`: The key to configure, following table shows possible `key` values, its default setting and its usage.
* `val::AbstractString`: The value to set to the key

Key                     |Default|Description
:-----------------------|:-----:|:---------------------------------------------------------------------------------
ALTERNATIVE_TABLES_DIRECTORY | "" | If provided, this path will be the root directory for the tabular data.  Otherwise, \${HOME}/.CoolProp/Tables is used
ALTERNATIVE_REFPROP_PATH | "" | An alternative path to be provided to the directory that contains REFPROP's fluids and mixtures directories.  If provided, the SETPATH function will be called with this directory prior to calling any REFPROP functions.
ALTERNATIVE_REFPROP_HMX_BNC_PATH | "" | An alternative path to the HMX.BNC file.  If provided, it will be passed into REFPROP's SETUP or SETMIX routines
VTPR_UNIFAC_PATH | "" | The path to the directory containing the UNIFAC JSON files.  Should be slash terminated
"""
function set_config(key::AbstractString, val::AbstractString)
    ccall( (:set_config_string, "CoolProp"), Nothing, (Cstring, Cstring), key, val)
    return get_global_param_string("errstring")
end

"""
    set_config(key::AbstractString, val::Real)

Set configuration numerical value as double.

# Arguments
* `key::AbstractString`: The key to configure, following table shows possible `key` values, its default setting and its usage.
* `val::Real`: The value to set to the key

Key                     |Default|Description
:-----------------------|:-----:|:---------------------------------------------------------------------------------
MAXIMUM_TABLE_DIRECTORY_SIZE_IN_GB |  1.0 | The maximum allowed size of the directory that is used to store tabular data
PHASE_ENVELOPE_STARTING_PRESSURE_PA |  100.0 | Starting pressure [Pa] for phase envelope construction
R_U_CODATA |  8.3144598 | The value for the ideal gas constant in J/mol/K according to CODATA 2014.  This value is used to harmonize all the ideal gas constants. This is especially important in the critical region.
SPINODAL_MINIMUM_DELTA |  0.5 | The minimal delta to be used in tracing out the spinodal; make sure that the EOS has a spinodal at this value of delta=rho/rho_r
"""
function set_config(key::AbstractString, val::Real)
    ccall( (:set_config_double, "CoolProp"), Nothing, (Cstring, Cdouble), key, val)
    return get_global_param_string("errstring")
end

"""
    set_config(key::AbstractString, val::Bool)

Set configuration value as a boolean.

# Arguments
* `key::AbstractString`: The key to configure, following table shows possible `key` values, its default setting and its usage.
* `val::Bool`: The value to set to the key

Key                     |Default|Description
:-----------------------|:-----:|:---------------------------------------------------------------------------------
NORMALIZE_GAS_CONSTANTS |  true | If true, for mixtures, the molar gas constant (R) will be set to the CODATA value
CRITICAL_WITHIN_1UK |  true | If true, any temperature within 1 uK of the critical temperature will be considered to be AT the critical point
CRITICAL_SPLINES_ENABLED |  true | If true, the critical splines will be used in the near-vicinity of the critical point
SAVE_RAW_TABLES |  false | If true, the raw, uncompressed tables will also be written to file
REFPROP_DONT_ESTIMATE_INTERACTION_PARAMETERS |  false | If true, if the binary interaction parameters in REFPROP are estimated, throw an error rather than silently continuing
REFPROP_IGNORE_ERROR_ESTIMATED_INTERACTION_PARAMETERS |  false | If true, if the binary interaction parameters in REFPROP are unable to be estimated, silently continue rather than failing
REFPROP_USE_GERG |  false | If true, rather than using the highly-accurate pure fluid equations of state, use the pure-fluid EOS from GERG-2008
REFPROP_USE_PENGROBINSON |  false | If true, rather than using the highly-accurate pure fluid equations of state, use the Peng-Robinson EOS
DONT_CHECK_PROPERTY_LIMITS |  false | If true, when possible, CoolProp will skip checking whether values are inside the property limits
HENRYS_LAW_TO_GENERATE_VLE_GUESSES |  false | If true, when doing water-based mixture dewpoint calculations, use Henry's Law to generate guesses for liquid-phase composition
OVERWRITE_FLUIDS |  false | If true, and a fluid is added to the fluids library that is already there, rather than not adding the fluid (and probably throwing an exception), overwrite it
OVERWRITE_DEPARTURE_FUNCTION |  false | If true, and a departure function to be added is already there, rather than not adding the departure function (and probably throwing an exception), overwrite it
OVERWRITE_BINARY_INTERACTION |  false | If true, and a pair of binary interaction pairs to be added is already there, rather than not adding the binary interaction pair (and probably throwing an exception), overwrite it
"""
function set_config(key::AbstractString, val::Bool)
    ccall( (:set_config_bool, "CoolProp"), Nothing, (Cstring, UInt8), key, val)
    return get_global_param_string("errstring")
end

export CoolProp_parameters, CoolProp_fluids;
"""
# CoolProp parameters table, to build run `CoolProp.buildparameters()`

$(isfile(abspath(@__FILE__, "..", "parameters.table")) ? readstring(abspath(@__FILE__, "..", "parameters.table")) : "")
"""
const CoolProp_parameters = "Type `?CoolProp_arameters` to get a list of all CoolProp parameters."
buildparameters() = begin
    logf = open("parameters.table", "w");
    println(logf, "Paramerer |Description |Unit |Comment ");
    println(logf, ":---------|:-----------|:----|:-------" );
    counter = 0;
    longunits = Set();
    for p in coolpropparameters
        longunit = get_parameter_information_string(p, "long") * " | " * get_parameter_information_string(p, "units");
        note = "";
        if (!in(longunit, longunits))
            push!(longunits, longunit);
        else
        note = " *Duplicated* "
        end
        for fluid in coolpropfluids
            try
                res = ("$(PropsSI(p, fluid))");
                note *= " **Constant Property** "
                break;
            catch err
            end
        end
        println(logf, "$p" * " | " * longunit * " | " * note);
    end
    close(logf);
end
"""
# CoolProp fluids table, to build run `CoolProp.buildfluids()`

$(isfile(abspath(@__FILE__, "..", "fluids.table")) ? readstring(abspath(@__FILE__, "..", "fluids.table")) : "")
"""
const CoolProp_fluids = "Type `?CoolProp_fluids` to get a list of all CoolProp fluids."
buildfluids() = begin
    logf = open("fluids.table", "w");
    println(logf, "ID |Name |Alias |CAS |Pure |Formula |BibTeX ");
    println(logf, ":--|:----|:-----|:---|:----|:-------|:------");
    id = 0;
    for fluid in coolpropfluids
        id+=1;
        print(logf, "$id | $fluid | $(get_fluid_param_string(fluid, "aliases"))");
        print(logf, " | $(get_fluid_param_string(fluid, "CAS"))");
        pure = get_fluid_param_string(fluid, "pure");
        print(logf, " | $pure");
        print(logf, " | $(get_fluid_param_string(fluid, "formula")) | ");
        for bi in ["BibTeX-CONDUCTIVITY", "BibTeX-EOS", "BibTeX-CP0", "BibTeX-SURFACE_TENSION","BibTeX-MELTING_LINE","BibTeX-VISCOSITY"]
            print(logf, " $bi:$(get_fluid_param_string(fluid, bi))");
        end
        print(logf, "\n");
    end
    close(logf);
end
# ---------------------------------
#       Information functions
# ---------------------------------

"""
    get_global_param_string(key::AbstractString)

Get a globally-defined string.

# Ref
ref CoolProp::get_global_param_string

# Arguments
* `key`: A string represents parameter name, could be one of $inputs_to_get_global_param_string
"""
function get_global_param_string(key::AbstractString)
    val = ccall( (:get_global_param_string, "CoolProp"), Clong, (Cstring, Ptr{UInt8}, Int), key, message_buffer::Array{UInt8, 1}, buffer_length)
    if val == 0
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
end

"""
    get_parameter_information_string(key::AbstractString, outtype::AbstractString)
    get_parameter_information_string(key::AbstractString)

Get information for a parameter.

# Arguments
* `key`: A string represents parameter name, to see full list check "Table of string inputs to PropsSI function": http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table, or simply type `get_global_param_string("parameter_list")`
* `outtype="long"`: Output type, could be one of the `["IO", "short", "long", "units"]`, with a default value of "long"

# Example
```julia
julia> get_parameter_information_string("HMOLAR")
"Molar specific enthalpy"
julia> get_parameter_information_string("HMOLAR", "units")
"J/mol"
```

# Note
A tabular output for this function is available with `?CoolProp_parameters`
"""
function get_parameter_information_string(key::AbstractString, outtype::AbstractString)
    message_buffer[1:length(outtype)+1] = [Vector{UInt8}(outtype); 0x00]
    val = ccall( (:get_parameter_information_string, "CoolProp"), Clong, (Cstring, Ptr{UInt8}, Int), key, message_buffer::Array{UInt8, 1}, buffer_length)
    if val == 0
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
end

function get_parameter_information_string(key::AbstractString)
    return get_parameter_information_string(key, "long")
end

"""
    get_fluid_param_string(fluid::AbstractString, param::AbstractString)

Get a string for a value from a fluid (numerical values for the fluid can be obtained from Props1SI function).

# Arguments
* `fluid`: The name of the fluid that is part of CoolProp, for instance "n-Propane"
* `param`: A string, can be in one of the terms described in the following table

ParamName                    | Description
:----------------------------|:----------------------------------------
"aliases"                    | A comma separated list of aliases for the fluid
"CAS", "CAS_number"          | The CAS number
"ASHRAE34"                   | The ASHRAE standard 34 safety rating
"REFPROPName", "REFPROP_name"| The name of the fluid used in REFPROP
"Bibtex-XXX"                 | A BibTeX key, where XXX is one of the bibtex keys used in get_BibTeXKey
"pure"                       | "true" if the fluid is pure, "false" otherwise
"formula"                    | The chemical formula of the fluid in LaTeX form if available, "" otherwise

# Note
A tabular output for this function is available with `?CoolProp_fluids`
"""
function get_fluid_param_string(fluid::AbstractString, param::AbstractString)
    val = ccall( (:get_fluid_param_string, "CoolProp"), Clong, (Cstring, Cstring, Ptr{UInt8}, Int), fluid, param, message_buffer::Array{UInt8, 1}, buffer_length)
    if val == 0
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
end

"""
    F2K(tf::Real)

Convert from degrees Fahrenheit to Kelvin (useful primarily for testing).
"""
function F2K(tf::Real)
    return ccall( (:F2K, "CoolProp"), Cdouble, (Cdouble,), tf)
end

"""
    K2F(tk::Real)

Convert from Kelvin to degrees Fahrenheit (useful primarily for testing).
"""
function K2F(tk::Real)
    return ccall( (:K2F, "CoolProp"), Cdouble, (Cdouble,), tk)
end

"""
    get_param_index(param::AbstractString)

Get the index as a long for a parameter "T", "P", etc, for `abstractstate_keyed_output()` function.

# Arguments
* `param`: A string represents parameter name, to see full list check "Table of string inputs to PropsSI function": http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table, or simply type `get_global_param_string("parameter_list")`
"""
function get_param_index(param::AbstractString)
    val = ccall( (:get_param_index, "CoolProp"), Clong, (Cstring,), param)
    if val == -1
        error("CoolProp: Unknown parameter: ", param)
    end
    return val
end

"""
    get_input_pair_index(pair::AbstractString)

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

# Example
```julia
julia> get_input_pair_index("PT_INPUTS")
9
```
"""
function get_input_pair_index(pair::AbstractString)
    val = ccall( (:get_input_pair_index, "CoolProp"), Clong, (Cstring,), pair)
    if val == -1
        error("CoolProp: Unknown input pair: ", pair)
    end
    return val
end
const coolpropparameters = map(Compat.String, split(get_global_param_string("parameter_list"),','));
const coolpropfluids = map(Compat.String, split(get_global_param_string("FluidsList"),','));

# ---------------------------------
# Getter and setter for debug level
# ---------------------------------

"""
    get_debug_level()

Get the debug level.

# Return value
Level The level of the verbosity for the debugging output (0-10) 0: no debgging output
"""
function get_debug_level()
    ccall( (:get_debug_level, "CoolProp"), Cint, () )
end

"""
    set_debug_level(level::Integer)

Set the debug level.

# Arguments
* `level::Integer`: The level of the verbosity for the debugging output (0-10) 0: no debgging output
"""
function set_debug_level(level::Integer) # change ::Int to ::Integer to make set_debug_level(get_debug_level()) works on different machine
    ccall( (:set_debug_level, "CoolProp"), Nothing, (Cint,), level)
end

# ---------------------------------
#        Humid Air Properties
# ---------------------------------

"""
    HAPropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, name3::AbstractString, value3::Real)

DLL wrapper of the HAPropsSI function.

# Arguments
* `output`: Output name for desired property, accepetd values are listed in the following table
* `name1`, `name2`, `name3`: Input name of  given state values, one must be "P"

# Note
Here, all outputs calculated as functions of these three: "Y"(Water mole fraction), "T" and "P", as "P" is mandatory so more performance achieved when "Y" and "T" is given (or at least one of them).

Parameter(s) name     |Description                                 |Unit      |Formula
:---------------------|:-------------------------------------------|:---------|:-------------------------------
Omega HumRat W        |Humidity ratio                              |          |
psi_w Y               |Water mole fraction                         |mol_w/mol |
Tdp T_dp DewPoint D   |Dew point temperature                       |K         |
Twb T_wb WetBulb B    |Wet bulb temperature                        |K         |
Enthalpy H Hda        |Enthalpy                                    |J/kg_da   |
Hha                   |Enthalpy per kg of humid air                |J/kg_ha   |
InternalEnergy U Uda  |Internal energy                             |J/kg_da   |
Uha                   |Internal energy per kg of humid air         |J/kg_ha   |
Entropy S Sda         |Entropy                                     |J/kg_da/K |
Sha                   |Entropy per kg of humid air                 |J/kg_ha/K |
RH RelHum R           |Relative humidity                           |          |
Tdb T_db T            |Temperature                                 |K         |
P                     |Pressure                                    |Pa        |
V Vda                 |Specific volume                             |m^3/kg_da |``MolarVolume*(1+HumidityRatio)/M_ha``
Vha                   |Specific volume per kg of humid air         |m^3/kg_ha |``MolarVolume/M_ha``
mu Visc M             |Viscosity                                   |          |
k Conductivity K      |Conductivity                                |          |
C cp                  |Heat cap. const. press.                     |J/kg_da/K |
Cha cp_ha             |Heat cap. const. press. per kg of humid air |J/kg_ha/K |
CV                    |Heat cap. const. vol.                       |J/kg_da/K |
CVha cv_ha            |Heat cap. const. vol. per kg of humid air   |J/kg_ha/K |
P_w                   |Partial pressure of water                   |          |
isentropic_exponent   |Isentropic exponent                         |          |
speed_of_sound        |Speed of sound                              |          |``sqrt(1/M_ha*cp/cv*dpdrho__constT)``
Z                     |Compressibility factor                      |          |``P*MolarVolume/(R*T)``

# Example
```julia
# Enthalpy (J per kg dry air) as a function of temperature, pressure,
# and relative humidity at STP
julia> h = HAPropsSI("H", "T", 298.15, "P", 101325, "R", 0.5)
50423.45039247604
# Temperature of saturated air at the previous enthalpy
julia> T = HAPropsSI("T", "P", 101325, "H", h, "R", 1.0)
290.9620891952412
# Temperature of saturated air - order of inputs doesn't matter
julia> T = HAPropsSI("T", "H", h, "R", 1.0, "P", 101325)
290.9620891952412
```

# Ref
HumidAir::HAPropsSI(const char* OutputName, const char* Input1Name, double Input1, const char* Input2Name, double Input2, const char* Input3Name, double Input3);
"""
function HAPropsSI(output::AbstractString, name1::AbstractString, value1::Real, name2::AbstractString, value2::Real, name3::AbstractString, value3::Real)
    val = ccall( (:HAPropsSI, "CoolProp"), Cdouble, (Cstring, Cstring, Cdouble, Cstring, Cdouble, Cstring, Cdouble), output, name1, value1, name2, value2, name3, value3)
    if val == Inf
        error("CoolProp: ", get_global_param_string("errstring"))
    end
    return val
end

"""
    cair_sat(t::Real)

Humid air saturation specific in [kJ/kg-K] heat at 1 atmosphere, based on a correlation from EES.

# Arguments
* `t`: T [K] good from 250K to 300K, no error bound checking is carried out.

# Ref
HumidAir::cair_sat(double);

# Note
Equals partial derivative of enthalpy with respect to temperature at constant relative humidity of 100 percent and pressure of 1 atmosphere.
"""
function cair_sat(t::Real)
    val = ccall( (:cair_sat, "CoolProp"), Cdouble, (Cdouble, ), t)
    return val;
end

function raise(errcode, message_buffer)
    if errcode[] != 0
        if errcode[] == 1
            error("CoolProp: ", unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer))))
        elseif errcode[] == 2
            error("CoolProp: message buffer too small")
        else # == 3
            error("CoolProp: unknown error")
        end
    end
end
# ---------------------------------
#        Low-level access
# ---------------------------------

"""
    AbstractState_factory(backend::AbstractString, fluids::AbstractString)

Generate an AbstractState instance return an integer handle to the state class generated to be used in the other low-level accessor functions.

# Arguments
* `backend`: The backend you will use could be: `["HEOS", "REFPROP", "INCOMP", "IF97", "TREND", "HEOS&TTSE", "HEOS&BICUBIC", "SRK", "PR", "VTPR"]` etc.
* `fluids`: '&' delimited list of fluids. To get a list of possible values call `get_global_param_string(key)` with `key` one of the following: `["FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "predefined_mixtures"]`, also there is a list in CoolProp online documentation [List of Fluids](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids), or simply type `?CoolProp_fluids`

# Example
```julia
julia> HEOS = AbstractState_factory("HEOS", "R245fa");
julia> TTSE = AbstractState_factory("HEOS&TTSE", "R245fa");
julia> BICU = AbstractState_factory("HEOS&BICUBIC", "R245fa");
julia> SRK = AbstractState_factory("SRK", "R245fa");
julia> PR = AbstractState_factory("PR", "R245fa");
```
"""
function AbstractState_factory(backend::AbstractString, fluids::AbstractString)
    AbstractState = ccall( (:AbstractState_factory, "CoolProp"), Clong, (Cstring, Cstring, Ref{Clong}, Ptr{UInt8}, Clong), backend, fluids, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return AbstractState
end

"""
    AbstractState_free(handle::Clong)

Release a state class generated by the low-level interface wrapper.

# Arguments
* `handle`: The integer handle for the state class stored in memory
"""
function AbstractState_free(handle::Clong)
    ccall( (:AbstractState_free, "CoolProp"), Nothing, (Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_set_fractions(handle::Clong, fractions::Array{Float64})

Set the fractions (mole, mass, volume) for the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `fractions`: The array of fractions

# Example
```julia
julia> handle = AbstractState_factory("HEOS", "Water&Ethanol");
julia> pq_inputs = get_input_pair_index("PQ_INPUTS");
julia> t = get_param_index("T");
julia> AbstractState_set_fractions(handle, [0.4, 0.6]);
julia> AbstractState_update(handle, pq_inputs, 101325, 0);
julia> AbstractState_keyed_output(handle, t)
352.3522212991724
julia> AbstractState_free(handle);
```
"""
function AbstractState_set_fractions(handle::Clong, fractions::Array{Float64})
    ccall( (:AbstractState_set_fractions, "CoolProp"), Nothing, (Clong, Ptr{Cdouble}, Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, fractions, length(fractions), errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_update(handle::Clong, input_pair::Clong, value1::Real, value2::Real)
    AbstractState_update(handle::Clong, input_pair::AbstractString, value1::Real, value2::Real)

Update the state of the AbstractState.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `input_pair::Clong`: The integer value for the input pair obtained from get_input_pair_index(param::AbstractString)
* `input_pair::AbstractString`: The name of an input pair
* `value1`: The first input value
* `value2`: The second input value

# Example
```julia
julia> handle = AbstractState_factory("HEOS", "Water&Ethanol");
julia> pq_inputs = get_input_pair_index("PQ_INPUTS");
julia> t = get_param_index("T");
julia> AbstractState_set_fractions(handle, [0.4, 0.6]);
julia> AbstractState_update(handle, pq_inputs, 101325, 0);
julia> AbstractState_keyed_output(handle, t)
352.3522212991724
julia> AbstractState_free(handle);
```
"""
function AbstractState_update(handle::Clong, input_pair::Clong, value1::Real, value2::Real)
    ccall( (:AbstractState_update, "CoolProp"), Nothing, (Clong, Clong, Cdouble, Cdouble, Ref{Clong}, Ptr{UInt8}, Clong), handle, input_pair, value1, value2, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

function AbstractState_update(handle::Clong, input_pair::AbstractString, value1::Real, value2::Real)
    AbstractState_update(handle::Clong, get_input_pair_index(input_pair), value1::Real, value2::Real)
    return nothing
end

"""
    AbstractState_keyed_output(handle::Clong, param::Clong)

Get an output value from the `AbstractState` using an integer value for the desired output value.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `param::Clong`: param The integer value for the parameter you want

# Note
See `AbstractState_output`
"""
function AbstractState_keyed_output(handle::Clong, param::Clong)
    output = ccall( (:AbstractState_keyed_output, "CoolProp"), Cdouble, (Clong, Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, param, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    if output == -Inf
        error("CoolProp: no correct state has been set with AbstractState_update")
    end
    return output
end

"""
    AbstractState_output(handle::Clong, param::AbstractString)

Get an output value from the `AbstractState` using an integer value for the desired output value. It is a convenience function that call `AbstractState_keyed_output`

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `param::AbstractString`: The name for the parameter you want
"""
function AbstractState_output(handle::Clong, param::AbstractString)
    return AbstractState_keyed_output(handle, get_param_index(param))
end


"""
    AbstractState_specify_phase(handle::Clong, phase::AbstractString)

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

# Example
```julia
julia> heos = AbstractState_factory("HEOS", "Water");
# Do a flash call that is a very low density state point, definitely vapor
julia> @time AbstractState_update(heos, "DmolarT_INPUTS", 1e-6, 300);
  0.025233 seconds (5.23 k allocations: 142.283 KB)
# Specify the phase - for some inputs (especially density-temperature), this will result in a
# more direct evaluation of the equation of state without checking the saturation boundary
julia> AbstractState_specify_phase(heos, "phase_gas");
# We try it again - a bit faster
julia> @time AbstractState_update(heos, "DmolarT_INPUTS", 1e-6, 300);
  0.000050 seconds (5 allocations: 156 bytes)
julia> AbstractState_free(heos);
```
"""
function AbstractState_specify_phase(handle::Clong, phase::AbstractString)
    ccall( (:AbstractState_specify_phase, "CoolProp"), Nothing, (Clong, Cstring, Ref{Clong}, Ptr{UInt8}, Clong), handle, phase, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_unspecify_phase(handle::Clong)

Unspecify the phase to be used for all further calculations.

# Arguments
* `handle`: The integer handle for the state class stored in memory
"""
function AbstractState_unspecify_phase(handle::Clong)
    ccall( (:AbstractState_unspecify_phase, "CoolProp"), Nothing, (Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_update_and_common_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, hmolar::Array{Float64}, smolar::Array{Float64})
    AbstractState_update_and_common_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, hmolar::Array{Float64}, smolar::Array{Float64})

Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `input_pair::Clong`: The integer value for the input pair obtained from get_input_pair_index
* `input_pair::AbstractString`:
* `value1`: The pointer to the array of the first input parameters
* `value2`: The pointer to the array of the second input parameters
* `length`: The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
* `T`: The pointer to the array of temperature
* `p`: The pointer to the array of pressure
* `rhomolar`: Array of molar density
* `hmolar`: The array of molar enthalpy
* `smolar`: Trray of molar entropy

# Example
```julia
julia> handle = AbstractState_factory("HEOS", "Water&Ethanol");
julia> pq_inputs = get_input_pair_index("PQ_INPUTS");
julia> AbstractState_set_fractions(handle, [0.4, 0.6]);
julia> T = [0.0]; p = [0.0]; rhomolar = [0.0]; hmolar = [0.0]; smolar = [0.0];
julia> AbstractState_update_and_common_out(handle, pq_inputs, [101325.0], [0.0], 1, T, p, rhomolar, hmolar, smolar);
julia> AbstractState_free(handle);
```
"""
function AbstractState_update_and_common_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, hmolar::Array{Float64}, smolar::Array{Float64})
    ccall( (:AbstractState_update_and_common_out, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Clong, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Clong}, Ptr{UInt8}, Clong), handle, input_pair, value1, value2, length, T, p, rhomolar, hmolar, smolar, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return T, p, rhomolar, hmolar, smolar
end

function AbstractState_update_and_common_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, hmolar::Array{Float64}, smolar::Array{Float64})
    return AbstractState_update_and_common_out(handle, get_input_pair_index(input_pair), value1, value2, length, T, p, rhomolar, hmolar, smolar)
end

function AbstractState_update_and_common_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer)
    T, p, rhomolar, hmolar, smolar = [fill(NaN,length) for i=1:5]
    return AbstractState_update_and_common_out(handle, input_pair, value1, value2, length, T, p, rhomolar, hmolar, smolar)
end

function AbstractState_update_and_common_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer)
    return AbstractState_update_and_common_out(handle, get_input_pair_index(input_pair), value1, value2, length)
end

"""
    AbstractState_update_and_1_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::Clong, out::Array{Float64})
    AbstractState_update_and_1_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::AbstractString, out::Array{Float64})

Update the state of the AbstractState and get one output value (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `input_pair::Clong`: The integer value for the input pair obtained from get_input_pair_index
* `input_pair::AbstractString`:
* `value1`: The pointer to the array of the first input parameters
* `value2`: The pointer to the array of the second input parameters
* `length`: The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
* `output`: The indice for the output desired
* `out`: The array for output
"""
function AbstractState_update_and_1_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::Clong, out::Array{Float64})
    ccall( (:AbstractState_update_and_1_out, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Clong, Clong, Ref{Cdouble}, Ref{Clong}, Ptr{UInt8}, Clong), handle, input_pair, value1, value2, length, output, out, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return out
end

function AbstractState_update_and_1_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::AbstractString, out::Array{Float64})
    return AbstractState_update_and_1_out(handle, get_input_pair_index(input_pair), value1, value2, length, get_param_index(output), out)
end

function AbstractState_update_and_1_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::Clong)
    out = fill(NaN,length)
    return AbstractState_update_and_1_out(handle, input_pair, value1, value2, length, output, out)
end

function AbstractState_update_and_1_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, output::AbstractString)
    return AbstractState_update_and_1_out(handle, get_input_pair_index(input_pair), value1, value2, length, get_param_index(output))
end

"""
    AbstractState_update_and_5_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{Clong}, out1::Array{Float64}, out2::Array{Float64}, out3::Array{Float64}, out4::Array{Float64}, out5::Array{Float64})
    AbstractState_update_and_5_out{S<:AbstractString}(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{S}, out1::Array{Float64}, out2::Array{Float64}, out3::Array{Float64}, out4::Array{Float64}, out5::Array{Float64})

Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `input_pair::Clong`: The integer value for the input pair obtained from get_input_pair_index
* `input_pair::AbstractString`:
* `value1`: The pointer to the array of the first input parameters
* `value2`: The pointer to the array of the second input parameters
* `length`: The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
* `outputs`: The 5-element vector of indices for the outputs desired
* `out1`: The array for the first output
* `out2`: The array for the second output
* `out3`: The array for the third output
* `out4`: The array for the fourth output
* `out5`: The array for the fifth output
"""
function AbstractState_update_and_5_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{Clong}, out1::Array{Float64}, out2::Array{Float64}, out3::Array{Float64}, out4::Array{Float64}, out5::Array{Float64})
    ccall( (:AbstractState_update_and_5_out, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Clong, Ref{Clong}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Clong}, Ptr{UInt8}, Clong), handle, input_pair, value1, value2, length, outputs, out1, out2, out3, out4, out5, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return out1, out2, out3, out4, out5
end

function AbstractState_update_and_5_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{S}, out1::Array{Float64}, out2::Array{Float64}, out3::Array{Float64}, out4::Array{Float64}, out5::Array{Float64}) where {S<:AbstractString}
    outputs_key = [get_param_index(outputs[k]) for k = 1:5]
    return AbstractState_update_and_5_out(handle, get_input_pair_index(input_pair), value1, value2, length, outputs_key, out1, out2, out3, out4, out5)
end

function AbstractState_update_and_5_out(handle::Clong, input_pair::Clong, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{Clong})
    out1, out2, out3, out4, out5 = [fill(NaN,length) for i=1:5]
    return AbstractState_update_and_5_out(handle, input_pair, value1, value2, length, outputs, out1, out2, out3, out4, out5)
end

function AbstractState_update_and_5_out(handle::Clong, input_pair::AbstractString, value1::Array{Float64}, value2::Array{Float64}, length::Integer, outputs::Array{S}) where {S<:AbstractString}
    outputs_key = [get_param_index(outputs[k]) for k = 1:5]
    return AbstractState_update_and_5_out(handle, get_input_pair_index(input_pair), value1, value2, length, outputs_key)
end

"""
    AbstractState_set_binary_interaction_double(handle::Clong, i::Int, j::Int, parameter::AbstractString, value::Float64)

Set binary interraction parameter for different mixtures model e.g.: "linear", "Lorentz-Berthelot"

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `i`: indice of the first fluid of the binary pair
* `j`: indice of the second fluid of the binary pair
* `parameter`: string wit the name of the parameter, e.g.: "betaT", "gammaT", "betaV", "gammaV"
* `value`: the value of the binary interaction parameter

# Example
```julia
julia> handle = AbstractState_factory("HEOS", "Water&Ethanol");
julia> AbstractState_set_binary_interaction_double(handle, 0, 1, "betaT", 0.987);
julia> pq_inputs = get_input_pair_index("PQ_INPUTS");
julia> t = get_param_index("T");
julia> AbstractState_set_fractions(handle, [0.4, 0.6]);
julia> AbstractState_update(handle, pq_inputs, 101325, 0);
julia> AbstractState_keyed_output(handle, t)
349.32634425309755
julia> AbstractState_free(handle);
```
"""
function AbstractState_set_binary_interaction_double(handle::Clong, i::Integer, j::Integer, parameter::AbstractString, value::Real)
    ccall( (:AbstractState_set_binary_interaction_double, "CoolProp"), Nothing, (Clong, Clong, Clong, Cstring, Cdouble, Ref{Clong}, Ptr{UInt8}, Clong), handle, i, j, parameter, value, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_set_cubic_alpha_C(handle::Clong, i::Integer, parameter::AbstractString, c1::Real, c2::Real, c3::Real)

Set cubic's alpha function parameters.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `i`: indice of the fluid the parameter should be applied too (for mixtures)
* `parameter`: the string specifying the alpha function to use, e.g. "TWU" for the Twu or "MC" for Mathias-Copeman alpha function.
* `c1`: the first parameter for the alpha function
* `c2`: the second parameter for the alpha function
* `c3`: the third parameter for the alpha function
"""
function AbstractState_set_cubic_alpha_C(handle::Clong, i::Integer, parameter::AbstractString, c1::Real, c2::Real, c3::Real)
    ccall( (:AbstractState_set_cubic_alpha_C, "CoolProp"), Nothing, (Clong, Clong, Cstring, Cdouble, Cdouble, Cdouble, Ref{Clong}, Ptr{UInt8}, Clong), handle, i, parameter, c1, c2, c3, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_set_fluid_parameter_double(handle::Clong, i::Integer, parameter::AbstractString, value::Real)

Set some fluid parameter (ie volume translation for cubic). Currently applied to the whole fluid not to components.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `i`: indice of the fluid the parameter should be applied to (for mixtures)
* `parameter`: string wit the name of the parameter, e.g. "c", "cm", "c_m" for volume translation parameter.
* `value`: the value of the parameter
"""
function AbstractState_set_fluid_parameter_double(handle::Clong, i::Integer, parameter::AbstractString, value::Real)
    ccall( (:AbstractState_set_fluid_parameter_double, "CoolProp"), Nothing, (Clong, Clong, Cstring, Cdouble, Ref{Clong}, Ptr{UInt8}, Clong), handle, i, parameter, value, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_first_saturation_deriv(handle::Clong, of::Clong, wrt::Clong)

Calculate a saturation derivative from the AbstractState using integer values for the desired parameters.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `of`: The parameter of which the derivative is being taken
* `wrt`: The derivative with with respect to this parameter

# Example
```julia
julia> as = AbstractState_factory("HEOS", "Water");
julia> AbstractState_update(as, "PQ_INPUTS", 15e5, 0);
julia> AbstractState_first_saturation_deriv(as, get_param_index("Hmolar"), get_param_index("P"))
0.0025636362140578207
```

# Ref
double CoolProp::AbstractState_first_saturation_deriv(const long handle, const long Of, const long Wrt, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_first_saturation_deriv(handle::Clong, of::Clong, wrt::Clong)
    output = ccall( (:AbstractState_first_saturation_deriv, "CoolProp"), Cdouble, (Clong, Clong, Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, of, wrt, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    if output == -Inf
        error("CoolProp: no correct state has been set with AbstractState_update")
    end
    return output
end

"""
    AbstractState_first_partial_deriv(handle::Clong, of::Clong, wrt::Clong, constant::Clong)

Calculate the first partial derivative in homogeneous phases from the AbstractState using integer values for the desired parameters.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `of`: The parameter of which the derivative is being taken
* `Wrt`: The derivative with with respect to this parameter
* `Constant`: The parameter that is not affected by the derivative

# Example
```julia
julia> as = AbstractState_factory("HEOS", "Water");
julia> AbstractState_update(as, "PQ_INPUTS", 15e5, 0);
julia> AbstractState_first_partial_deriv(as, get_param_index("Hmolar"), get_param_index("P"), get_param_index("S"))
2.07872526058326e-5
julia> AbstractState_first_partial_deriv(as, get_param_index("Hmolar"), get_param_index("P"), get_param_index("D"))
5.900781297636475e-5
```

# Ref
double CoolProp::AbstractState_first_partial_deriv(const long handle, const long Of, const long Wrt, const long Constant, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_first_partial_deriv(handle::Clong, of::Clong, wrt::Clong, constant::Clong)
    output = ccall( (:AbstractState_first_partial_deriv, "CoolProp"), Cdouble, (Clong, Clong, Clong, Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, of, wrt, constant, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    if output == -Inf
        error("CoolProp: no correct state has been set with AbstractState_update")
    end
    return output
end

"""
    AbstractState_build_phase_envelope(handle::Clong, level::AbstractString)

Build the phase envelope.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `level`: How much refining of the phase envelope ("none" to skip refining (recommended) or "veryfine")

# Note
If there is an error in an update call for one of the inputs, no change in the output array will be made

# Ref
CoolPRop::AbstractState_build_phase_envelope(const long handle, const char* level, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_build_phase_envelope(handle::Clong, level::AbstractString)
    ccall( (:AbstractState_build_phase_envelope, "CoolProp"), Nothing, (Clong, Cstring, Ref{Clong}, Ptr{UInt8}, Clong), handle, level, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
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

# Example
```julia
julia> HEOS=AbstractState_factory("HEOS","Methane&Ethane");
julia> length=200;
julia> t=zeros(length);p=zeros(length);x=zeros(2*length);y=zeros(2*length);rhomolar_vap=zeros(length);rhomolar_liq=zeros(length);
julia> AbstractState_set_fractions(HEOS, [0.2, 1 - 0.2])
julia> AbstractState_build_phase_envelope(HEOS, "none")
julia> AbstractState_get_phase_envelope_data(HEOS, length, t, p, rhomolar_vap, rhomolar_liq, x, y)
julia> AbstractState_free(HEOS)
```

# Note
If there is an error in an update call for one of the inputs, no change in the output array will be made

# Ref
CoolProp::AbstractState_get_phase_envelope_data(const long handle, const long length, double* T, double* p, double* rhomolar_vap, double* rhomolar_liq, double* x, double* y, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_get_phase_envelope_data(handle::Clong, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar_vap::Array{Float64}, rhomolar_liq::Array{Float64}, x::Array{Float64}, y::Array{Float64})
    ccall( (:AbstractState_get_phase_envelope_data, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Clong}, Ptr{UInt8}, Clong), handle, length, T, p, rhomolar_vap, rhomolar_liq, x, y, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return T, p, rhomolar_vap, rhomolar_liq, x, y
end

function AbstractState_get_phase_envelope_data(handle::Clong, length::Integer, ncomp::Integer)
    T, p, rhomolar_vap, rhomolar_liq = [fill(NaN,length) for i=1:5]
    x, y = [fill(NaN,length*ncomp) for i=1:2]
    return AbstractState_get_phase_envelope_data(handle, length, T, p, rhomolar_vap, rhomolar_liq, x, y)
end

"""
    AbstractState_build_spinodal(handle::Clong)

Build the spinodal.

# Arguments
* `handle`: The integer handle for the state class stored in memory

# Ref
CoolProp::AbstractState_build_spinodal(const long handle, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_build_spinodal(handle::Clong)
    ccall( (:AbstractState_build_spinodal, "CoolProp"), Nothing, (Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return nothing
end

"""
    AbstractState_get_spinodal_data(handle::Clong, length::Integer, tau::Array{Float64}, dalta::Array{Float64}, m1::Array{Float64})

Get data for the spinodal curve.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `length`: The number of elements stored in the arrays (all outputs MUST be the same length)
* `tau`: The pointer to the array of reciprocal reduced temperature
* `delta`: The pointer to the array of reduced density
* `m1`: The pointer to the array of M1 values (when L1=M1=0, critical point)

# Note
If there is an error, no change in the output arrays will be made

# Example
julia> HEOS=AbstractState_factory("HEOS","Methane&Ethane");
julia> AbstractState_set_fractions(HEOS, [0.1, 0.9]);
julia> AbstractState_build_spinodal(HEOS);
julia> tau, delta, m1 = AbstractState_get_spinodal_data(HEOS, 127);

# Ref
CoolProp::AbstractState_get_spinodal_data(const long handle, const long length, double* tau, double* delta, double* M1, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_get_spinodal_data(handle::Clong, length::Integer, tau::Array{Float64}, delta::Array{Float64}, m1::Array{Float64})
    ccall( (:AbstractState_get_spinodal_data, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Clong}, Ptr{UInt8}, Clong), handle, length, tau, delta, m1, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return tau, delta, m1;
end

function AbstractState_get_spinodal_data(handle::Clong, length::Integer)
    tau, delta, m1 = [fill(NaN,length) for i=1:3]
    return AbstractState_get_spinodal_data(handle, length, tau, delta, m1)
end

"""
    abstractState_all_critical_points(handle::Clong, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, stable::Array{Clong})

Calculate all the critical points for a given composition.

# Arguments
* `handle`: The integer handle for the state class stored in memory
* `length`: The length of the buffers passed to this function
* `T`: The pointer to the array of temperature (K)
* `p`: The pointer to the array of pressure (Pa)
* `rhomolar`: The pointer to the array of molar density (m^3/mol)
* `stable`: The pointer to the array of boolean flags for whether the critical point is stable (1) or unstable (0)

# Note
If there is an error in an update call for one of the inputs, no change in the output array will be made

# Ref
CoolProp::AbstractState_all_critical_points(const long handle, const long length, double* T, double* p, double* rhomolar, long* stable, long* errcode, char* message_buffer, const long buffer_length);
"""
function AbstractState_all_critical_points(handle::Clong, length::Integer, T::Array{Float64}, p::Array{Float64}, rhomolar::Array{Float64}, stable::Array{Clong})
    ccall( (:AbstractState_all_critical_points, "CoolProp"), Nothing, (Clong, Clong, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Clong}, Ref{Clong}, Ptr{UInt8}, Clong), handle, length, T, p, rhomolar, stable, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    raise(errcode, message_buffer)
    return T, p, rhomolar, stable
end

function AbstractState_all_critical_points(handle::Clong, length::Integer)
    T, p, rhomolar = [fill(NaN,length) for i=1:3]
    stable = zeros(Clong, length)
    return  AbstractState_all_critical_points(handle, length, T, p, rhomolar, stable)
end

for symorigin = [:PropsSI, :PhaseSI, :K2F, :F2K, :HAPropsSI, :AbstractState_factory, :AbstractState_free, :AbstractState_set_fractions, :AbstractState_update, :AbstractState_keyed_output, :AbstractState_output, :AbstractState_specify_phase, :AbstractState_unspecify_phase, :AbstractState_update_and_common_out, :AbstractState_update_and_1_out, :AbstractState_update_and_5_out, :AbstractState_set_binary_interaction_double, :AbstractState_set_cubic_alpha_C, :AbstractState_set_fluid_parameter_double, :AbstractState_first_saturation_deriv, :AbstractState_first_partial_deriv, :AbstractState_build_phase_envelope, :AbstractState_build_spinodal, :AbstractState_all_critical_points, :AbstractState_get_phase_envelope_data, :AbstractState_get_spinodal_data]
    sym = Symbol(lowercase(string(symorigin)))
    @eval const $sym = $symorigin
    @eval export $sym, $symorigin
end
const set_reference_stateS = set_reference_state
const set_reference_stateD = set_reference_state
const set_config_string = set_config
export set_reference_stateS, set_reference_stateD, set_reference_state
export get_global_param_string, get_parameter_information_string, get_fluid_param_string, get_param_index, get_input_pair_index, set_config
export saturation_ancillary, set_departure_functions, set_config_string, cair_sat
end #module
