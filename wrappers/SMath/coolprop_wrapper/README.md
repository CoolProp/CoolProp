# List of functions

Below there is the list of all the functions provided by this plugin inside [SMath Studio](http://en.smath.info).

****
```
CoolProp_get_fluid_param_string(FluidName, ParamName)
```
*Get a string for a value from a fluid.*
* `FluidName` The name of the fluid that is part of CoolProp, for instance `"n-Propane"`
* `ParamName` A string, can be in one of `"aliases"`, `"CAS"`, `"CAS_number"`, `"ASHRAE34"`, `"REFPROPName"`, `"REFPROP_name"`

****
```
CoolProp_get_global_param_string(ParamName)
```
*Get a globally-defined string.*
* `ParamName` A string, one of `"version"`, `"errstring"`, `"warnstring"`, `"gitrevision"`, `"FluidsList"`, `"fluids_list"`, `"parameter_list"`, `"predefined_mixtures"`

****
```
CoolProp_get_param_index(Name)
```
*Return the index of a parameter.*
* `Name`: The parameter name, one of `"Tcrit"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)

****
```
CoolProp_get_parameter_information_string(Key, Output)
```
*Get a parameter information string.*
* `Key` A string
* `Output` A string, one of `"IO"`, `"short"`, `"long"`, `"units"`

****
```
CoolProp_HAProps(Output, Name1, Prop1, Name2, Prop2, Name3, Prop3)
```
*Return a humid air property.*
* `Output` The output parameter, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Name1` The first state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop1` The first state variable value
* `Name2` The second state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop2` The second state variable value
* `Name3` The third state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop3` The third state variable value

****
```
CoolProp_Phase(Name1, Prop1, Name2, Prop2, FluidName)
```
*Return a string representation of the phase.*
* `Name1` The first state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop1` The first state variable value
* `Name2` The second state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop2` The second state variable value
* `FluidName` The fluid name

****
```
CoolProp_Props(Output, Name1, Prop1, Name2, Prop2, FluidName)
```
*Return a value that __depends__ on the thermodynamic state.*
* `Output` The output parameter, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Name1` The first state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop1` The first state variable value
* `Name2` The second state variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `Prop2` The second state variable value
* `FluidName` The fluid name

****
```
CoolProp_Props1(FluidName, Output)
```
*Return a value that does __not depends__ on the thermodynamic state.*
* `FluidName` The fluid name
* `Output` The output parameter, one of `"Tcrit"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)

****
```
CoolProp_saturation_ancillary(FluidName, output, Q, input, value)
```
*Extract a value from the saturation ancillary.*
* `FluidName` The name of the fluid to be used - HelmholtzEOS backend only
* `output` The desired output variable (\"P\" for instance for [pressure](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table))
* `Q` The mass vapor quality, 0 or 1
* `input` The input variable name, one of `"T"`, `"D"`, `"H"`, [etc...](http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table)
* `value` The second state variable value

****
```
CoolProp_set_reference_stateD(FluidName, T, RhoMolar, h0, s0)
```
*Set the reference state based on a thermodynamic state point specified by temperature and molar density.*
* `FluidName` The name of the fluid
* `T` Temperature at reference state [K]
* `RhoMolar` Density at reference state [mol/m^3]
* `h0` Enthalpy at reference state [J/mol]
* `s0` Entropy at references state [J/mol/K]

****
```
CoolProp_set_reference_stateS(FluidName, ReferenceState)
```
*Set the reference state based on a string representation*
* `FluidName` The name of the fluid
* `ReferenceState` The reference state to use, one of `"IIR"`, `"ASHRAE"`, `"NBP"`, `"DEF"`, `"RESET"`
