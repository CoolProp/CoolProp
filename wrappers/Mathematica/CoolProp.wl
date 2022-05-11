(* ::Package:: *)

(* ::Title:: *)
(*CoolProp *)


(* ::Subtitle:: *)
(*Wrapper for Mathematica*)


(* ::Text:: *)
(*This package contains the functions added to the CoolProp wrapper for Mathematica. *)


(* ::Subsection::Closed:: *)
(*General Remarks*)


(* ::Text:: *)
(*"Package Head" contains the general Package Building Blocks and usage messages for all implemented functions that are supposed to be called from outside the package. The setup for all these functions is the same:*)


(* ::Item:: *)
(*FunctionNamefunc: LibraryFunction, Namespace: CoolProp`Private`*)


(* ::Subitem:: *)
(*can be called by using complete function name CoolProp`Private`FunctionNamefunc, but this is not advised (no error handling)*)


(* ::Subitem:: *)
(*to move functions from Private to Public, simply add a usage message in the section "Package Head".*)


(* ::Item:: *)
(*FunctionName: Namespace: CoolProp`*)


(* ::Subitem:: *)
(*callable by just FunctionName when package is loaded*)


(* ::Subitem:: *)
(*uses pattern matching to assure you are only entering correct inputs to the LibraryFunction*)


(* ::Subitem:: *)
(*gets error messages generated inside CoolProp if the result is not as expected (empty string, Infinity, negative index number, error code (in case of void functions)) and returns $Failed*)


(* ::Subitem:: *)
(*overloaded function without pattern matching generates error message with false inputs and returns $Failed*)


(* ::Text:: *)
(*"Support Functions" contains Private functions that are used by other functions within the package.*)
(*"Direct Functions" contains functions that will take all information about the state you want to know about and return the value requested. PropsSI, PhaseSI and HAPropsSI.*)
(*"Information Functions" contains all informational functions currently implemented in the wrapper, i.e. index requests for input pairs or parameters.*)
(*"Low-Level Access: AbstractState" contains functions to interact with the low-level interface by using AbstractStates. Consult CoolProp Documentation for explanation. *)


(* ::Section::Closed:: *)
(*Package Head*)


BeginPackage["CoolProp`"]

SetBufferLength::usage = "SetBufferLength can change the buffer lengths defined in CoolPropMathematica.cpp. Options are: 
	\"stringBuffer\": string output and error messages, set to 20000
	\"resultMatrix\": dimensions of matrix output (e.g. PropsSImulti), both are set at the same time to the same value, set to 20
	\"maxComponents\":  maximum number of components supported by get/set fraction functions
SetBufferLength[bufferName, bufferLength]";
Props1SI::usage = "Return a value that does not depend on the thermodynamic state. 
PropsSI[param, fluid]
	param: name of return value
	fluid: name of fluid";
Props1SImulti::usage = "Get a matrix of outputs that do not depend on the thermodynamic state
Props1SImulti[outputs, backend, fluidNames, fractions]
	outputs: Delimited string separated by LIST_STRING_DELIMITER for the output parameters
	backend: The string representation of the backend (HEOS, REFPROP, INCOMP, etc.) 
	fluidNames: Delimited string separated by LIST_STRING_DELIMITER for the fluid name(s)
	frations: The fractions (molar, mass, volume, etc.) of the components";
PropsSI::usage = "Return a value that depends on the thermodynamic state. 
PropsSI[param, name1, value1, name2, value2, fluid]
	param: name of return value
	name1: first input variable name
	value1: value for name1
	name2: second input variable name
	value2: value for name2
	fluid: name of fluid";
PropsSImulti::usage = "Get a matrix of outputs for a given input.  Can handle both vector inputs as well as a vector of output strings
PropsSImulti[outputs, name1, values1, name2, values2, backend, fluidNames, fractions]
	outputs: Delimited string separated by LIST_STRING _DELIMITER for the output parameters
	name1: The name of the first input variable
	values1: A vector of the first input values
	name2: The name of the second input variable
	values2: A vector of the second input values 
	backend: The string representation of the backend (HEOS, REFPROP, INCOMP, etc.) 
	fluidNames: Delimited string separated by LIST_STRING _DELIMITER for the fluid name(s)
	frations: The fractions (molar, mass, volume, etc.) of the components";
PhaseSI::usage = "Return a string representation of the phase. 
PhaseSI[name1, value1, name2, value2, fluid]
	name1: first input variable name
	value1: value for name1
	name2: second input variable name
	value2: value for name2
	fluid: name of fluid";
HAPropsSI::usage = "Return a value that depends on the thermodynamic state for humid air. 
HAPropsSI[param, name1, value1, name2, value2, name3, value3]
	param: name of return value
	name1: first input variable name
	value1: value for name1
	name2: second input variable name
	value2: value for name2
	name3: third input variable name
	value3: value for name3";

GetGlobalParamString::usage = "Get a globally-defined string, one of \"version\", \"errstring\", \"warnstring\", \"gitrevision\", \"FluidsList\", \"fluids_list\", \"parameter_list\",\"predefined_mixtures\".
GetGlobalParamString[param]
	param: A string represents parameter name";
GetInputPairIndex::usage = "Get the index for an input pair for ASUpdate function.
GetInputPairIndex[inputPair]
	inputPair: The name of an input pair";
GetInputPairDescription::usage = "Get the information about an input pair for the given index.
GetInputPairDescription[index, outtype: \"long\"]
	index: index (result of GetInputPairIndex)
	outtype: Output type, can be \"long\", \"short\" or \"both\";"
GetParamIndex::usage = "Return the enum key corresponding to the parameter name (\"Dmolar\" for instance).
GetParamIndex[param]
	param:  A string represents parameter name";
GetParameterInformation::usage = "Get the information about a parameter for the given index.
GetParameterInformation[index, outtype: \"long\"]
	index: index (result of GetParamIndex
	outtype: Output type, could be one of the `[\"IO\", \"short\", \"long\", \"units\"]`, with a default value of \"long\"";
CubicFluidsList::usage = "Returns all fluids implemented in the CubicsLibrary.
CubicFluidsList[]";
AddFluidsAsJSON::usage = "Add fluids as a JSON-formatted string
AddFluidsAsJSON[backend, fluidstring]
	backend: The backend to which these should be added; e.g. \"HEOS\", \"SRK\", \"PR\"
	fluiddstring: The JSON-formatted string";
GetFluidsAsJSON::usage = "Return fluid as a JSON-formatted string
GetFluidsAsJSON[backend, fluidstring]
	backend: The backend from which these should be read; e.g. \"SRK\", \"PR\"
	fluiddstring: The fluid whose JSON is wanted";
FluidParamString::usage = "Get a string for a value from a fluid
FluidParamString[fluid, parameter]
	fluid: The name of the fluid that is part of CoolProp, for instance \"n-Propane\". Different backends accessible by prepending it e.g.: \"PR::Propane\"
	parameter: A string [aliases, CAS, CAS_number, ASHRAE34, REFPROPName, REFPROP_name, Bibtex-XXX, pure, formula";
	
ASFactory::usage = "Generate an AbstractState instance, return an integer handle to the state class generated to be used in the other low-level accessor functions.
ASFactory[backend, fluid]
	backend: which backend is used for the AbstractState
	fluid: which fluid is used for the AbstractState. Separate mixture components with \"&\"";
ASFree::usage = "Release a state class generated by the low-level interface wrapper.
ASFree[handle]
	handle: integer handle assigned by ASFactory";
ASUpdate::usage = "Update the state of the AbstractState.
ASUpdate[handle, inputPair, value1, value2]
	handle: integer handle assigned by ASFactory
	inputPair: The integer value for the input pair obtained from GetInputPairIndex
	value1: first input value
	value2: second input value";
ASKeyedOutput::usage = "Get an output value from the AbstractState using an integer value for the desired output value.
ASKeyedOutput[handle, param]
	handle: integer handle assigned by ASFactory
	param: integer value for the parameter to return, obtained from GetParamIndex";
ASKeyedOutputSatState::usage = "Get an output value from the AbstractState using an integer value for the desired output value and desired saturated state.
ASKeyedOutput[handle, param, satState]
	handle: integer handle assigned by ASFactory
	param: integer value for the parameter to return, obtained from GetParamIndex
	satSatte: string value specifying the state (\"liquid\" or \"gas\")";
ASSetFractions::usage = "Set the fractions (mole, mass, volume) for the AbstractState.
ASSetFractions[handle, fractions]
	handle: integer handle assigned by ASFactory
	fractions: list of fractions, use same order as when creating the AbstractState";		
ASGetMoleFractions::usage = "Get the mole fractions for the AbstractState.
ASSetFractions[handle]
	handle: integer handle assigned by ASFactory";
ASGetMoleFractionsSatState::usage = "Get the mole fractions for the AbstractState for the desired saturated state.
ASSetFractions[handle, satState]
	handle: integer handle assigned by ASFactory
	satSatte: string value specifying the state (\"liquid\" or \"gas\")";
ASSetBinaryInteraction::usage = "Set binary interaction parrameter for mixtures. 
ASSetBinaryInteraction[handle, componentI, componentJ, param, value]
	handle: integer handle assigned by ASFactory
	componentI: integer, position of first fluid in mixture (starts at 0)
	componentJ: integer, position of first second in mixture (starts at 0)
	param: string wit the name of the parameter, e.g.: \"betaT\", \"gammaT\", \"betaV\", \"gammaV\", \"kij\"";
ASSetCubicAlpha::usage = "Set cubic's alpha function parameters
	ASSetCubicAlpha[handle, componentI,  param, parameters]
		handle: integer handle assigned by AbstractStateFactory
		componentI: integer, indice of the fluid the parramter should be applied too (for mixtures)
		param: parameter the string specifying the alpha function to use, ex \"TWU\" for the TWU alpha function
		parameters: three parameters for the alpha function";
ASFluidNames::usage = "Get the fluid names for the AbstractState.
ASFluidNames[handle]
	handle: integer handle assigned by ASFactory";
ASSpecifyPhase::usage = "Impose a phase for the AbstractState.
ASSpecifyPhase[handle, phase]
	handle: integer handle assigned by ASFactory
	phase: string phase";
ASBuildPhaseEnvelope::usage = "Build the phase envelope for the AbstractState.
ASBuildPhaseEnvelope[handle, level]
	handle: integer handle assigned by ASFactory
	level: How much refining of the phase envelope (\"none\" to skip refining (recommended), \"veryfine\")";
ASReadPhaseEnvelopeData::usage = "Get data from the phase envelope for the given mixture composition.
ASReadPhaseEnvelopeData[handle, length:500, numberOfFluids:10]
	handle: integer handle assigned by ASFactory
	length: length of the buffer to read the phase envelope, can be set higher if necessary
	numberOfFluids: number of fluids used in the mixture, can be set higher if necessary";
ASBackendName::usage = "Returns the backend used for the AbstractState.
ASBackendName[handle]
	handle: integer handle assigned by ASFactory";

	
Begin["`Private`"]


(* ::Section::Closed:: *)
(*Support Functions*)


(* ::Subsection::Closed:: *)
(*CheckTypes*)


(* ::Text:: *)
(*CheckTypes is used by all Functions to check, which Arguments are the wrong Type.*)


CheckTypes::error = "Type \"`1`\" not implemented";
CheckTypes[function_, valueList_List, typeList_List]:=
	Module[
		{
			error = ""
		},
		If[
			Not@Switch[typeList[[#]], 
					"Real", NumericQ@valueList[[#]],
					"String", StringQ@valueList[[#]],
					"Integer", IntegerQ@valueList[[#]],
					"Boolean", BooleanQ@valueList[[#]],
					"ListOfReals", And[AllTrue[valueList[[#]], NumericQ], ListQ@valueList[[#]]],
					_, Message[CheckTypes::error, typeList[[#]]];False
					], 
			If[error == "", 
				error = "TypeError: ", 
				error = error <> ", "];
			If[And[ListQ@valueList[[#]], typeList[[#]] == "ListOfReals"],
				error = error <> "all elements of Argument " <> ToString@# <> " (" <>ToString@valueList[[#]] <> ")"<> " should be of type Number", 
				error = error <> "Argument " <> ToString@# <> " (value = " <>ToString@valueList[[#]] <> ", type = " <> ToString@Head@valueList[[#]] <> ")"<> " should be of type " <> typeList[[#]]]
			
		] &/@ Range@Length@typeList;
		
		function::typeError = error;
		Message[function::typeError];
		$Failed
	]


(* ::Subsection::Closed:: *)
(*SetBufferLength*)


(* ::Text:: *)
(*SetBufferLength can change the buffer lengths defined in CoolPropMathematica.cpp. Options are: *)


(* ::Item:: *)
(*"stringBuffer": string output and error messages, set to 20000*)


(* ::Item:: *)
(*"resultMatrix": dimensions of matrix output (e.g. PropsSImulti), both are set at the same time to the same value, set to 20*)


(* ::Item:: *)
(*"maxComponents":  maximum number of components supported by get/set fraction functions*)


SetBufferLengthfunc = LibraryFunctionLoad["libCoolProp", "set_buffer_length", {"UTF8String", Integer}, Integer];


SyntaxInformation[SetBufferLength] = {"ArgumentsPattern" -> {_,_}};


SetBufferLength[bufferName_String, bufferLength_Integer] :=
	Module[
		{
			result = SetBufferLengthfunc[bufferName, bufferLength]
		},
		If[result != 0,
			SetBufferLength::libraryFunctionError = GetGlobalParamString["errstring"]; Message[SetBufferLength::libraryFunctionError]; $Failed
		]
	];


legalTypeSetBufferLength = {"String", "Integer"};


SetBufferLength[bufferName_, bufferLength_] := CheckTypes[SetBufferLength, {bufferName, bufferLength}, legalTypeSetBufferLength]


(* ::Section:: *)
(*Direct Functions*)


(* ::Subsection::Closed:: *)
(*Props1SImulti*)


Props1SImultifunc = LibraryFunctionLoad["libCoolProp","Props1SImulti",{"UTF8String", "UTF8String", "UTF8String", {Real,_}}, {Real,_}];


SyntaxInformation[Props1SImulti] = {"ArgumentsPattern" -> {_,_,_,_}};


Props1SImulti[outputs_String, backend_String, fluidNames_String, fractions:{__?NumericQ}] :=
	Module[
		{
			result = Props1SImultifunc[outputs, backend, fluidNames, fractions]
		},
		If[result == {},
			Props1SImulti::libraryFunctionError = GetGlobalParamString["errstring"]; Message[Props1SImulti::libraryFunctionError]; $Failed,
			result
		]
	];


(* ::Text:: *)
(*legalTypeProps1SImulti contains the types, the arguments for Props1SImulti should have. "Real" is used, but they could be Integers and Reals. CheckType uses NumericQ to check if they are either one of those.*)


legalTypeProps1SImulti = {"String", "String", "String", "ListOfReals"};


Props1SImulti[outputs_, backend_, fluidNames_, fractions_] := CheckTypes[Props1SImulti, {outputs, backend, fluidNames, fractions}, legalTypeProps1SImulti]


(* ::Subsection::Closed:: *)
(*Props1SI*)


Props1SIfunc = LibraryFunctionLoad["libCoolProp","Props1SI",{"UTF8String","UTF8String"}, Real];


SyntaxInformation[Props1SI] = {"ArgumentsPattern" -> {_,_}};


Props1SI[param_String, fluid_String] :=
	Module[
		{
			result = Check[Props1SIfunc[param, fluid], $Failed]
		},
		If[Not@NumericQ@result,
			If[result === $Failed,
				Props1SI::libraryFunctionError = "Error evaluating LibraryFunction.",
				Props1SI::libraryFunctionError = GetGlobalParamString["errstring"]
			]; Message[Props1SI::libraryFunctionError]; $Failed,
			result]
	];


(* ::Text:: *)
(*legalTypeProps1SI contains the types, the arguments for Props1SI should have. "Real" is used, but they could be Integers and Reals. CheckType uses NumericQ to check if they are either one of those.*)


legalTypeProps1SI = {"String", "String"};


Props1SI[param_, fluid_] := CheckTypes[Props1SI, {param, fluid}, legalTypeProps1SI]


(* ::Subsection::Closed:: *)
(*PropsSI*)


PropsSIfunc = LibraryFunctionLoad["libCoolProp","PropsSI",{"UTF8String","UTF8String",Real,"UTF8String",Real,"UTF8String"}, Real];


SyntaxInformation[PropsSI] = {"ArgumentsPattern" -> {_,_,_,_,_,_}};


PropsSI[param_String, name1_String, value1_?NumericQ, name2_String, value2_?NumericQ, fluid_String] :=
	Module[
		{
			result = Check[PropsSIfunc[param, name1, value1, name2, value2, fluid], $Failed]
		},
		If[Not@NumericQ@result,
			If[result === $Failed,
				PropsSI::libraryFunctionError = "Error evaluating LibraryFunction.",
				PropsSI::libraryFunctionError = GetGlobalParamString["errstring"]
			]; Message[PropsSI::libraryFunctionError]; $Failed,
			result]
	];


(* ::Text:: *)
(*legalTypePropsSI contains the types, the arguments for PropsSI should have. "Real" is used, but they could be Integers and Reals. CheckType uses NumericQ to check if they are either one of those.*)


legalTypePropsSI = {"String", "String", "Real", "String", "Real", "String"};


PropsSI[param_, name1_, value1_, name2_, value2_, fluid_] := CheckTypes[PropsSI, {param, name1, value1, name2, value2, fluid}, legalTypePropsSI]


(* ::Subsection::Closed:: *)
(*PropsSImulti*)


PropsSImultifunc = LibraryFunctionLoad["libCoolProp","PropsSImulti",{"UTF8String", "UTF8String", {Real,_}, "UTF8String", {Real,_},  "UTF8String", "UTF8String", {Real,_}}, {Real,_}];


SyntaxInformation[PropsSImulti] = {"ArgumentsPattern" -> {_,_,_,_,_,_,_,_}};


PropsSImulti[outputs_String, name1_String, values1:{__?NumericQ}, name2_String, values2:{__?NumericQ}, backend_String, fluidNames_String, fractions:{__?NumericQ}] :=
	Module[
		{
			result = PropsSImultifunc[outputs, name1, values1, name2, values2, backend, fluidNames, fractions]
		},
		If[result == {},
			PropsSImulti::libraryFunctionError = GetGlobalParamString["errstring"]; Message[PropsSImulti::libraryFunctionError]; $Failed,
			result
		]
	];


(* ::Text:: *)
(*legalTypePropsSImulti contains the types, the arguments for PropsSImulti should have. "Real" is used, but they could be Integers and Reals. CheckType uses NumericQ to check if they are either one of those.*)


legalTypePropsSImulti = {"String", "String", "ListOfReals", "String", "ListOfReals", "String", "String", "ListOfReals"};


PropsSImulti[outputs_, name1_, values1_, name2_, values2_, backend_, fluidNames_, fractions_] := CheckTypes[PropsSImulti, {outputs, name1, values1, name2, values2, backend, fluidNames, fractions}, legalTypePropsSImulti]


(* ::Subsection::Closed:: *)
(*PhaseSI*)


validPhases = {"liquid", "supercritical","supercritical_gas", "supercritical_liquid", "critical_point", "gas", "twophase", "unknown", "not_imposed"};


PhaseSIfunc = LibraryFunctionLoad["libCoolProp","PhaseSI",{"UTF8String",Real,"UTF8String",Real,"UTF8String"}, "UTF8String"];


SyntaxInformation[PhaseSI] = {"ArgumentsPattern" -> {_,_,_,_,_}};


PhaseSI[name1_String, value1_?NumericQ, name2_String, value2_?NumericQ, fluid_String] :=
	Module[
		{
			result = PhaseSIfunc[name1, value1, name2, value2, fluid]
		},
		If[Not@MemberQ[validPhases, result],
			PhaseSI::libraryFunctionError = result; Message[PhaseSI::libraryFunctionError]; $Failed,
			result]
	];


legalTypePhaseSI = {"String", "Real", "String", "Real", "String"};


PhaseSI[name1_, value1_, name2_, value2_, fluid_] := CheckTypes[PhaseSI, {name1, value1, name2, value2, fluid}, legalTypePhaseSI]


(* ::Subsection::Closed:: *)
(*HAPropsSI*)


HAPropsSIfunc = LibraryFunctionLoad["libCoolProp","HAPropsSI",{"UTF8String","UTF8String",Real,"UTF8String",Real,"UTF8String",Real}, Real];


SyntaxInformation[HAPropsSI] = {"ArgumentsPattern" -> {_,_,_,_,_,_,_}};


HAPropsSI[param_String, name1_String, value1_?NumericQ, name2_String, value2_?NumericQ, name3_String, value3_?NumericQ] :=
	Module[
		{
			result = HAPropsSIfunc[param, name1, value1, name2, value2, name3, value3]
		},
		If[result == Infinity,
			HAPropsSI::libraryFunctionError = GetGlobalParamString["errstring"]; Message[HAPropsSI::libraryFunctionError]; $Failed,
			result]
	];


legalTypeHAPropsSI = {"String","String","Real","String","Real","String","Real"};


HAPropsSI[param_, name1_, value1_, name2_, value2_, name3_, value3_] := CheckTypes[HAPropsSI, {param, name1, value1, name2, value2, name3, value3}, legalTypeHAPropsSI]


(* ::Section:: *)
(*Information Functions*)


(* ::Subsection::Closed:: *)
(*GetGlobalParamString*)


(* ::Text:: *)
(*"version", "gitrevision", "errstring", "warnstring", "FluidsList", "incompressible_list_pure", "incompressible_list_solution", "mixture_binary_pairs_list", "parameter_list", "predefined_mixtures", "HOME", "cubic_fluids_schema"*)


GetGlobalParamStringfunc = LibraryFunctionLoad["libCoolProp","get_global_param_string",{"UTF8String"}, "UTF8String"];


SyntaxInformation[GetGlobalParamString] = {"ArgumentsPattern" -> {_}};


GetGlobalParamString[param_String] :=
	Module[
		{
			result = GetGlobalParamStringfunc[param]
		},
		If[result == "no string",
			GetGlobalParamString::libraryFunctionError = GetGlobalParamString["errstring"]; Message[GetGlobalParamString::libraryFunctionError]; $Failed,
			result]
	];


legalTypeGetGlobalParamString = {"String"};


GetGlobalParamString[param_] := CheckTypes[GetGlobalParamString, {param}, legalTypeGetGlobalParamString]


(* ::Subsection::Closed:: *)
(*GetInputPairIndex*)


(* ::Text:: *)
(*GetInputPairIndex returns the index for the given input pair. Used for the update function.*)


GetInputPairIndexfunc = LibraryFunctionLoad["libCoolProp","get_input_pair_index",{"UTF8String"}, Integer];


SyntaxInformation[GetInputPairIndex] = {"ArgumentsPattern" -> {_}};


GetInputPairIndex[inputPair_String] :=
	Module[
		{
			result = GetInputPairIndexfunc[inputPair]
		},
		If[result == -1,
			GetInputPairIndex::libraryFunctionError = GetGlobalParamString["errstring"]; Message[GetInputPairIndex::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeGetInputPairIndex = {"String"};


GetInputPairIndex[inputPair_] := CheckTypes[GetInputPairIndex, {inputPair}, legalTypeGetInputPairIndex]


(* ::Subsection::Closed:: *)
(*GetInputPairDescription*)


GetInputPairDescriptionfunc = LibraryFunctionLoad["libCoolProp","get_input_pair_description",{Integer, "UTF8String"}, "UTF8String"];


SyntaxInformation[GetInputPairDescription] = {"ArgumentsPattern" -> {_,_}};


GetInputPairDescription[index_Integer, outtype_String:"both"] :=
	Module[
		{
			result = GetInputPairDescriptionfunc[index, outtype]
		},
		If[result == "",
			GetInputPairDescription::libraryFunctionError = GetGlobalParamString["errstring"]; Message[GetInputPairDescription::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeGetInputPairDescription = {"Integer", "String"};


GetInputPairDescription[index_, outtype_] := CheckTypes[GetInputPairDescription, {index, outtype}, legalTypeGetInputPairDescription]


(* ::Subsection::Closed:: *)
(*GetParamIndex*)


(* ::Text:: *)
(*GetParamIndex returns the index for a given parameter.*)


GetParamIndexfunc = LibraryFunctionLoad["libCoolProp","get_param_index",{"UTF8String"}, Integer];


SyntaxInformation[GetParamIndex] = {"ArgumentsPattern" -> {_}};


GetParamIndex[param_String] :=
	Module[
		{
			result = GetParamIndexfunc[param]
		},
		If[result == -1,
			GetParamIndex::libraryFunctionError = GetGlobalParamString["errstring"]; Message[GetParamIndex::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeGetParamIndex = {"String"};


GetParamIndex[param_]:=CheckTypes[GetParamIndex,{param},legalTypeGetParamIndex]


(* ::Subsection::Closed:: *)
(*GetParameterInformation*)


GetParameterInformationfunc = LibraryFunctionLoad["libCoolProp","get_parameter_information",{Integer, "UTF8String"}, "UTF8String"];


SyntaxInformation[GetParameterInformation] = {"ArgumentsPattern" -> {_,_}};


GetParameterInformation[index_Integer, outtype_String:"long"] :=
	Module[
		{
			result = GetParameterInformationfunc[index, outtype]
		},
		If[result == "",
			GetParameterInformation::libraryFunctionError = GetGlobalParamString["errstring"]; Message[GetParameterInformation::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeGetParameterInformation = {"Integer", "String"};


GetParameterInformation[index_, outtype_] := CheckTypes[GetParameterInformation, {index, outtype}, legalTypeGetParameterInformation]


(* ::Subsection::Closed:: *)
(*CubicFluidsList*)


CubicFluidsList = LibraryFunctionLoad["libCoolProp","get_cubic_fluids_list",{},"UTF8String"];


(* ::Subsection::Closed:: *)
(*AddFluidsAsJSON*)


(* ::Text:: *)
(**)


AddFluidsAsJSONfunc = LibraryFunctionLoad["libCoolProp", "add_fluids_as_JSON", {"UTF8String", "UTF8String"}, Integer]


SyntaxInformation[AddFluidsAsJSON] = {"ArgumentsPattern" -> {_,_}};


AddFluidsAsJSON[backend_String, fluidstring_String] :=
	Module[
		{
			result = AddFluidsAsJSONfunc[backend, fluidstring]
		},
		If[result != 0,
			AddFluidsAsJSON::libraryFunctionError = GetGlobalParamString["errstring"]; Message[AddFluidsAsJSON::libraryFunctionError]; $Failed
		]
	];


legalTypeAddFluidsAsJSON = {"String", "String"};


AddFluidsAsJSON[backend_, fluidstring_] := CheckTypes[AddFluidsAsJSON, {backend, fluidstring}, legalTypeAddFluidsAsJSON]


(* ::Subsection::Closed:: *)
(*GetFluidsAsJSON*)


(* ::Text:: *)
(**)


SyntaxInformation[GetFluidsAsJSON] = {"ArgumentsPattern" -> {_,_}};


GetFluidsAsJSON[backend_String, fluidstring_String] := FluidParamString[backend<>"::"<>fluidstring, "JSON"]


legalTypeGetFluidsAsJSON = {"String", "String"};


GetFluidsAsJSON[backend_, fluidstring_] := CheckTypes[GetFluidsAsJSON, {backend, fluidstring}, legalTypeGetFluidsAsJSON]



(* ::Subsection::Closed:: *)
(*FluidParamString*)


(* ::Text:: *)
(**)


FluidParamStringfunc = LibraryFunctionLoad["libCoolProp", "get_fluid_param_string", {"UTF8String", "UTF8String"}, "UTF8String"];


SyntaxInformation[FluidParamString] = {"ArgumentsPattern" -> {_,_}};


FluidParamString[backend_String, fluidstring_String] :=
	Module[
		{
			result = FluidParamStringfunc[backend, fluidstring]
		},
		If[result == "",
			FluidParamString::libraryFunctionError = GetGlobalParamString["errstring"]; Message[FluidParamString::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeFluidParamString = {"String", "String"};


FluidParamString[backend_, fluidstring_] := CheckTypes[FluidParamString, {backend, fluidstring}, legalTypeFluidParamString]


(* ::Section:: *)
(*Low-Level Access: AbstractState*)


(* ::Subsection::Closed:: *)
(*ASFactory*)


(* ::Text:: *)
(*ASFactory was added to be able to use the low level interface per the high level interface. It returns the string handle for the Abstract State generated. This handle is later used to access the State.*)


ASFactoryfunc = LibraryFunctionLoad["libCoolProp","AbstractState_factory",{"UTF8String","UTF8String"}, Integer];


SyntaxInformation[ASFactory] = {"ArgumentsPattern" -> {_,_}};


ASFactory[backend_String, fluid_String] :=
	Module[
		{
			result = ASFactoryfunc[backend, fluid]
		},
		If[result == -1,
			ASFactory::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASFactory::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASFactory = {"String", "String"};


ASFactory[backend_, fluid_] := CheckTypes[ASFactory, {backend, fluid}, legalTypeASFactory]


(* ::Subsection::Closed:: *)
(*ASFree*)


(* ::Text:: *)
(*ASFree is the couterpart to ASFactory. It destroys the Abstract State specified by the handle. This frees the used memory.*)


ASFreefunc = LibraryFunctionLoad["libCoolProp","AbstractState_free",{Integer}, Integer];


SyntaxInformation[ASFree] = {"ArgumentsPattern" -> {_}};


ASFree[handle_Integer] :=
	Module[
		{
			result = ASFreefunc[handle]
		},
		If[result != 0,
			ASFree::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASFree::libraryFunctionError]; $Failed
		]
	];


legalTypeASFree = {"Integer"};


ASFree[handle_] := CheckTypes[ASFree, {handle}, legalTypeASFree]


(* ::Subsection::Closed:: *)
(*ASUpdate*)


(* ::Text:: *)
(*ASUpdate updates the Abstract State with index, input pair and both values.*)


ASUpdatefunc = LibraryFunctionLoad["libCoolProp","AbstractState_update",{Integer, Integer, Real, Real}, Integer];


SyntaxInformation[ASUpdate] = {"ArgumentsPattern" -> {_,_,_,_}};


ASUpdate[handle_Integer, inputPair_Integer, value1_?NumericQ, value2_?NumericQ] :=
	Module[
		{
			result = ASUpdatefunc[handle, inputPair, value1, value2]
		},
		If[result != 0,
			ASUpdate::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASUpdate::libraryFunctionError]; $Failed
		]
	];


legalTypeASUpdate = {"Integer", "Integer", "Real", "Real"};


ASUpdate[handle_, inputPair_, value1_, value2_] := CheckTypes[ASUpdate, {handle, inputPair, value1, value2}, legalTypeASUpdate]


(* ::Subsection::Closed:: *)
(*ASKeyedOutput*)


(* ::Text:: *)
(*ASKeyedOutput gives the value for the parameter for the set State.*)


ASKeyedOutputfunc = LibraryFunctionLoad["libCoolProp","AbstractState_keyed_output",{Integer, Integer}, Real];


SyntaxInformation[ASKeyedOutput] = {"ArgumentsPattern" -> {_,_}};


ASKeyedOutput[handle_Integer, param_Integer] :=
	Module[
		{
			result = ASKeyedOutputfunc[handle, param]
		},
		If[Or[result === Infinity, result === -Infinity, result === Indeterminate],
			ASKeyedOutput::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASKeyedOutput::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASKeyedOutput = {"Integer", "Integer"};


ASKeyedOutput[handle_, param_] := CheckTypes[ASKeyedOutput, {handle, param}, legalTypeASKeyedOutput]


(* ::Subsection::Closed:: *)
(*ASKeyedOutputSatState*)


(* ::Text:: *)
(*ASKeyedOutputSatState gives the value for the parameter for the set State.*)


ASKeyedOutputSatStatefunc = LibraryFunctionLoad["libCoolProp", "AbstractState_keyed_output_satState", {Integer, "UTF8String", Integer}, Real];


SyntaxInformation[ASKeyedOutputSatState] = {"ArgumentsPattern" -> {_,_,_}};


ASKeyedOutputSatState[handle_Integer, param_Integer, satState_String] :=
	Module[
		{
			result = ASKeyedOutputSatStatefunc[handle, satState, param]
		},
		If[Or[result === Infinity, result === -Infinity, result === Indeterminate],
			ASKeyedOutputSatState::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASKeyedOutputSatState::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASKeyedOutputSatState = {"Integer", "Integer", "String"};


ASKeyedOutputSatState[handle_, param_, satState_] := CheckTypes[ASKeyedOutputSatState, {handle, param, satState}, legalTypeASKeyedOutputSatState]


(* ::Subsection::Closed:: *)
(*ASSetFractions*)


(* ::Text:: *)
(*ASSetFractions sets the fractions in the Abstract State.*)


ASSetFractionsfunc = LibraryFunctionLoad["libCoolProp","AbstractState_set_fractions",{Integer, {Real,_}}, Integer];


SyntaxInformation[ASSetFractions] = {"ArgumentsPattern" -> {_,_}};


ASSetFractions[handle_Integer, fractions:{__?NumericQ}] :=
	Module[
		{
			result = ASSetFractionsfunc[handle, fractions]
		},
		If[result != 0,
			ASSetFractions::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASSetFractions::libraryFunctionError]; $Failed
		]
	];


legalTypeASSetFractions = {"Integer", "ListOfReals"};


ASSetFractions[handle_, fractions_] := CheckTypes[ASSetFractions, {handle, fractions}, legalTypeASSetFractions]


(* ::Subsection::Closed:: *)
(*ASGetMoleFractions*)


(* ::Text:: *)
(*ASGetMoleFractions gets the molar fractions in the Abstract State.*)


 ASGetMoleFractionsfunc = LibraryFunctionLoad["libCoolProp","AbstractState_get_mole_fractions",{Integer}, {Real,_}];


SyntaxInformation[ASGetMoleFractions] = {"ArgumentsPattern" -> {_}};


ASGetMoleFractions[handle_Integer] :=
	Module[
		{
			result = ASGetMoleFractionsfunc[handle]
		},
		If[result == {0.},
			ASGetMoleFractions::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASGetMoleFractions::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASGetMoleFractions = {"Integer"};


ASGetMoleFractions[handle_] := CheckTypes[ASGetMoleFractions, {handle}, legalTypeASGetMoleFractions]


(* ::Subsection::Closed:: *)
(*ASGetMoleFractionsSatState*)


(* ::Text:: *)
(*ASGetMoleFractions gets the molar fractions in the Abstract State.*)


 ASGetMoleFractionsSatStatefunc = LibraryFunctionLoad["libCoolProp","AbstractState_get_mole_fractions_satState",{Integer, "UTF8String"}, {Real,_}]


SyntaxInformation[ASGetMoleFractionsSatState] = {"ArgumentsPattern" -> {_,_}};


ASGetMoleFractionsSatState[handle_Integer, satState_String] :=
	Module[
		{
			result = ASGetMoleFractionsSatStatefunc[handle, satState]
		},
		If[result == {0.},
			ASGetMoleFractionsSatState::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASGetMoleFractionsSatState::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASGetMoleFractionsSatState = {"Integer", "String"};


ASGetMoleFractionsSatState[handle_, satState_] := CheckTypes[ASGetMoleFractionsSatState, {handle, satState}, legalTypeASGetMoleFractionsSatState]


(* ::Subsection::Closed:: *)
(*ASSetBinaryInteraction*)


(* ::Text:: *)
(*ASSetBinaryInteraction is used to set binary interaction parameters for different mixture models.*)


ASSetBinaryInteractionfunc = LibraryFunctionLoad["libCoolProp","AbstractState_set_binary_interaction_double",{Integer, Integer, Integer, "UTF8String", Real}, Integer]; 


SyntaxInformation[ASSetBinaryInteraction] = {"ArgumentsPattern" -> {_,_,_,_,_}};


ASSetBinaryInteraction[handle_Integer, componentI_Integer, componentJ_Integer, param_String, value_?NumericQ] :=
	Module[
		{
			result = ASSetBinaryInteractionfunc[handle, componentI, componentJ, param, value]
		},
		If[result != 0,
			ASSetBinaryInteraction::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASSetBinaryInteraction::libraryFunctionError]; $Failed
		]
	];


legalTypeASSetBinaryInteraction = {"Integer", "Integer", "Integer", "String", "Real"};


ASSetBinaryInteraction[handle_, componentI_, componentJ_, param_, value_] := CheckTypes[ASSetBinaryInteraction, {handle, componentI, componentJ, param, value}, legalTypeASSetBinaryInteraction]


(* ::Subsection::Closed:: *)
(*ASSetCubicAlpha*)


(* ::Text:: *)
(*ASSetCubicAlpha sets the fractions in the Abstract State.*)


ASSetCubicAlphafunc = LibraryFunctionLoad["libCoolProp","AbstractState_set_cubic_alpha_C",{Integer, Integer, "UTF8String", {Real,_}}, Integer];


SyntaxInformation[ASSetCubicAlpha] = {"ArgumentsPattern" -> {_,_,_,_}};


ASSetCubicAlpha[handle_Integer, componentI_Integer, param_String, parameters:{___?NumericQ}] :=
	Module[
		{
			result = ASSetCubicAlphafunc[handle, componentI, param, parameters]
		},
		If[result != 0,
			ASSetCubicAlpha::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASSetCubicAlpha::libraryFunctionError]; $Failed
		]
	];


legalTypeASSetCubicAlpha = {"Integer", "Integer", "String", "ListOfReals"};


ASSetCubicAlpha[handle_, componentI_, param_, parameters_] := CheckTypes[ASSetCubicAlpha, {handle, componentI, param, parameters}, legalTypeASSetCubicAlpha]


(* ::Subsection::Closed:: *)
(*ASFluidNames*)


ASFluidNamesfunc = LibraryFunctionLoad["libCoolProp","AbstractState_fluid_names",{Integer}, "UTF8String"];


SyntaxInformation[ASFluidNames] = {"ArgumentsPattern" -> {_}};


ASFluidNames[handle_Integer] :=
	Module[
		{
			result = ASFluidNamesfunc[handle]
		},
		If[result == "",
			ASFluidNames::libraryFunctionError = GetGlobalParamString["errstring"];Message[ASFluidNames::libraryFunctionError]; $Failed,
			StringSplit[result, ","]
		]
	];


legalTypeASFluidNames = {"Integer"};


ASFluidNames[handle_] := CheckTypes[ASFluidNames, {handle}, legalTypeASFluidNames]


(* ::Subsection::Closed:: *)
(*ASSpecifyPhase*)


(* ::Text:: *)
(*ASSpecifyPhase imposes the phase for Abstract State*)


ASSpecifyPhasefunc = LibraryFunctionLoad["libCoolProp","AbstractState_specify_phase",{Integer, "UTF8String"}, Integer];


SyntaxInformation[ASSpecifyPhase] = {"ArgumentsPattern" -> {_,_}};


ASSpecifyPhase[handle_Integer, phase_String] :=
	Module[
		{
			result = ASSpecifyPhasefunc[handle, phase]
		},
		If[result != 0,
			ASSpecifyPhase::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASSpecifyPhase::libraryFunctionError]; $Failed
		]
	];


legalTypeASSpecifyPhase = {"Integer", "String"};


ASSpecifyPhase[handle_, phase_] := CheckTypes[ASSpecifyPhase, {handle, phase}, legalTypeASSpecifyPhase]


(* ::Subsection::Closed:: *)
(*ASBuildPhaseEnvelope*)


(* ::Text:: *)
(*ASBuildPhaseEnvelope builds the phase envelope*)


ASBuildPhaseEnvelopefunc = LibraryFunctionLoad["libCoolProp", "AbstractState_build_phase_envelope", {Integer, "UTF8String"}, Integer];


SyntaxInformation[ASBuildPhaseEnvelope] = {"ArgumentsPattern" -> {_,_}};


ASBuildPhaseEnvelope[handle_Integer, level_String] :=
	Module[
		{
			result = ASBuildPhaseEnvelopefunc[handle, level]
		},
		If[result != 0,
			ASBuildPhaseEnvelope::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASBuildPhaseEnvelope::libraryFunctionError]; $Failed
		]
	];


legalTypeASBuildPhaseEnvelope = {"Integer", "String"};


ASBuildPhaseEnvelope[handle_, level_] := CheckTypes[ASBuildPhaseEnvelope, {handle, level}, legalTypeASBuildPhaseEnvelope]


(* ::Subsection::Closed:: *)
(*ASReadPhaseEnvelopeData*)


(* ::Text:: *)
(*ASReadPhaseEnvelopeData reads the phase envelope data (temperature, pressure, \[Rho]molar_vap and  \[Rho]molar_liq)*)


ASReadPhaseEnvelopeDatafunc = LibraryFunctionLoad["libCoolProp","AbstractState_get_phase_envelope_data",{Integer, Integer, Integer},{Real,_}];


SyntaxInformation[ASReadPhaseEnvelopeData] = {"ArgumentsPattern" -> {_}};


ASReadPhaseEnvelopeData[handle_Integer, length_Integer:500, numberOfFluids_Integer:10] :=
	Module[
		{
			result = ASReadPhaseEnvelopeDatafunc[handle, length, numberOfFluids]
		},
		If[result == {0},
			ASReadPhaseEnvelopeData::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASReadPhaseEnvelopeData::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASReadPhaseEnvelopeData = {"Integer", "Integer", "Integer"};


ASReadPhaseEnvelopeData[handle_, length_:500, numberOfFluids_:10] := CheckTypes[ASReadPhaseEnvelopeData, {handle, length, numberOfFluids}, legalTypeASReadPhaseEnvelopeData]


(* ::Subsection::Closed:: *)
(*ASBackendName*)


ASBackendNamefunc = LibraryFunctionLoad["libCoolProp","AbstractState_backend_name",{Integer}, "UTF8String"]


SyntaxInformation[ASBackendName] = {"ArgumentsPattern" -> {_}};


ASBackendName[handle_Integer] :=
	Module[
		{
			result = ASBackendNamefunc[handle]
		},
		If[result == "",
			ASBackendName::libraryFunctionError = GetGlobalParamString["errstring"]; Message[ASBackendName::libraryFunctionError]; $Failed,
			result
		]
	];


legalTypeASBackendName = {"Integer"};


ASBackendName[index_] := CheckTypes[ASBackendName, {index}, legalTypeASBackendName]


(* ::Section:: *)
(*Package End*)


End[];
EndPackage[];
