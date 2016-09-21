VERSION < v"0.4.0" && __precompile__()
module CoolProp

export PropsSI, PhaseSI, get_global_param_string, get_parameter_information_string,get_fluid_param_string,set_reference_stateS, get_param_index, get_input_pair_index, set_config_string, F2K, K2F, HAPropsSI, AbstractState_factory, AbstractState_free, AbstractState_set_fractions, AbstractState_update, AbstractState_specify_phase, AbstractState_unspecify_phase, AbstractState_keyed_output, AbstractState_output, AbstractState_update_and_common_out, AbstractState_update_and_1_out, AbstractState_update_and_5_out, AbstractState_set_binary_interaction_double, AbstractState_set_fluid_parameter_double

# Check the current Julia version to make this Julia 0.4 code compatible with older version
if VERSION <= VersionNumber(0,4)
  typealias UInt8 Uint8
  typealias Ref Ptr
  errcode = Array(Clong, 1)
  typealias AbstractString String
else
  errcode = Ref{Clong}(0)
end

const buffer_length = 10000
message_buffer = Array(UInt8, buffer_length)

# ---------------------------------
#        High-level functions
# ---------------------------------

# CoolProp::Props1SI(std::string, std::string)
function PropsSI(FluidName::AbstractString, Output::AbstractString)
  val = ccall( (:Props1SI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8}), FluidName,Output)
  if val == Inf
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return val
end

# CoolProp::PropsSI(const std::string &, const std::string &, double, const std::string &, double, const std::string&)
function PropsSI(Output::AbstractString, Name1::AbstractString, Value1::Real, Name2::AbstractString, Value2::Real, Fluid::AbstractString)
  val = ccall( (:PropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}), Output,Name1,Value1,Name2,Value2,Fluid)
  if val == Inf
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return val
end

# CoolProp::PhaseSI(const std::string &, double, const std::string &, double, const std::string&)
function PhaseSI(Name1::AbstractString, Value1::Real, Name2::AbstractString, Value2::Real, Fluid::AbstractString)
  val = ccall( (:PhaseSI, "CoolProp"), Int32, (Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}, Ptr{UInt8}, Int), Name1,Value1,Name2,Value2,Fluid,message_buffer::Array{UInt8,1},buffer_length)
  val = bytestring(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8,1})))
  if val == ""
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return val
end

# CoolProp::get_global_param_string
function get_global_param_string(Key::AbstractString)
  val = ccall( (:get_global_param_string, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Int), Key,message_buffer::Array{UInt8,1},buffer_length)
  if val == 0
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return bytestring(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8,1})))
end

# CoolProp::get_parameter_information_string
function get_parameter_information_string(Key::AbstractString,OutType::AbstractString)
  message_buffer[1:length(OutType)+1] = [OutType.data; 0x00]
  val = ccall( (:get_parameter_information_string, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Int), Key,message_buffer::Array{UInt8,1},buffer_length)
  if val == 0
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return bytestring(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8,1})))
end
function get_parameter_information_string(Key::AbstractString)
  return get_parameter_information_string(Key,"long")
end

# CoolProp::get_mixture_binary_pair_data
# Not valid yet
#function get_mixture_binary_pair_data(CAS1::AbstractString,CAS2::AbstractString,Key::AbstractString)
#  val = ccall( (:get_mixture_binary_pair_data, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Ptr{UInt8}), CAS1,CAS2,Key)
#  return val
#end

# CoolProp::get_fluid_param_string
function get_fluid_param_string(fluid::AbstractString,param::AbstractString)
  val = ccall( (:get_fluid_param_string, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Int), fluid,param,message_buffer::Array{UInt8,1},buffer_length)
  if val == 0
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return bytestring(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8,1})))
end

# CoolProp::set_reference_stateS
function set_reference_stateS(Ref::AbstractString, reference_state::AbstractString)
  val = ccall( (:set_reference_stateS, "CoolProp"), Cint, (Ptr{UInt8},Ptr{UInt8}), Ref,reference_state)
  if val == 0
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return val
end

# Get the index for a parameter "T", "P", etc.
# returns the index as a long.
function get_param_index(param::AbstractString)
  val = ccall( (:get_param_index, "CoolProp"), Clong, (Ptr{UInt8},), param)
  if val == -1
    error("CoolProp: Unknown parameter: ", param)
  end
  return val
end

# Get the index for an input pair for AbstractState.update function
# returns the index as a long.
function get_input_pair_index(param::AbstractString)
  val = ccall( (:get_input_pair_index, "CoolProp"), Clong, (Ptr{UInt8},), param)
  if val == -1
    error("CoolProp: Unknown input pair: ", param)
  end
  return val
end

# Set configuration string
function set_config_string(key::AbstractString,val::AbstractString)
  ccall( (:set_config_string, "CoolProp"), Void, (Ptr{UInt8},Ptr{UInt8}), key,val)
  return get_global_param_string("errstring")
end

# Get the debug level
# returns level The level of the verbosity for the debugging output (0-10) 0: no debgging output
function get_debug_level()
  ccall( (:get_debug_level, "CoolProp"), Cint, () )
end

# Set the debug level
# param level The level of the verbosity for the debugging output (0-10) 0: no debgging output
function set_debug_level(level::Int)
  ccall( (:set_debug_level, "CoolProp"), Void, (Cint,), level)
end

function F2K(TF::Real)
  return ccall( (:F2K, "CoolProp"), Cdouble, (Cdouble,), TF)
end

function K2F(TK::Real)
  return ccall( (:K2F, "CoolProp"), Cdouble, (Cdouble,), TK)
end

# ---------------------------------
#        Humid Air Properties
# ---------------------------------

function HAPropsSI(Output::AbstractString, Name1::AbstractString, Value1::Real, Name2::AbstractString, Value2::Real, Name3::AbstractString, Value3::Real)
  val = ccall( (:HAPropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8},Float64), Output,Name1,Value1,Name2,Value2,Name3,Value3)
  if val == Inf
    error("CoolProp: ", get_global_param_string("errstring"))
  end
  return val
end

# ---------------------------------
#        Low-level access
# ---------------------------------

# Generate an AbstractState instance, return an integer handle to the state class generated to be used in the other low-level accessor functions
# param backend The backend you will use, "HEOS", "REFPROP", etc.
# param fluids '&' delimited list of fluids
# return A handle to the state class generated
function AbstractState_factory(backend::AbstractString, fluids::AbstractString)
  AbstractState = ccall( (:AbstractState_factory, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Ref{Clong},Ptr{UInt8},Clong), backend,fluids,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return AbstractState
end

# Release a state class generated by the low-level interface wrapper
# param handle The integer handle for the state class stored in memory
function AbstractState_free(handle::Clong)
  ccall( (:AbstractState_free, "CoolProp"), Void, (Clong,Ref{Clong},Ptr{UInt8},Clong), handle,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

# Set the fractions (mole, mass, volume) for the AbstractState
# param handle The integer handle for the state class stored in memory
# param fractions The array of fractions
function AbstractState_set_fractions(handle::Clong,fractions::Array)
  ccall( (:AbstractState_set_fractions, "CoolProp"), Void, (Clong,Ptr{Cdouble},Clong,Ref{Clong},Ptr{UInt8},Clong), handle,fractions,length(fractions),errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

# Update the state of the AbstractState
# param handle The integer handle for the state class stored in memory
# param input_pair The integer value for the input pair obtained from get_input_pair_index(param::AbstractString)
# param value1 The first input value
# param value2 The second input value
function AbstractState_update(handle::Clong,input_pair::Clong,value1::Real,value2::Real)
  ccall( (:AbstractState_update, "CoolProp"), Void, (Clong,Clong,Cdouble,Cdouble,Ref{Clong},Ptr{UInt8},Clong), handle,input_pair,value1,value2,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end
function AbstractState_update(handle::Clong,input_pair::AbstractString,value1::Real,value2::Real)
  AbstractState_update(handle::Clong,get_input_pair_index(input_pair),value1::Real,value2::Real)
  return nothing
end

# Specify the phase to be used for all further calculations
# handle The integer handle for the state class stored in memory
# phase The string with the phase to use
# errcode The errorcode that is returned (0 = no error, !0 = error)
# message_buffer A buffer for the error code
# buffer_length The length of the buffer for the error code
function AbstractState_specify_phase(handle::Clong,phase::AbstractString)
  ccall( (:AbstractState_specify_phase, "CoolProp"), Void, (Clong,Ptr{UInt8},Ref{Clong},Ptr{UInt8},Clong), handle,phase,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

# Unspecify the phase to be used for all further calculations
# handle The integer handle for the state class stored in memory
# errcode The errorcode that is returned (0 = no error, !0 = error)
# message_buffer A buffer for the error code
# buffer_length The length of the buffer for the error code
function AbstractState_unspecify_phase(handle::Clong)
  ccall( (:AbstractState_unspecify_phase, "CoolProp"), Void, (Clong,Ref{Clong},Ptr{UInt8},Clong), handle,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

# Get an output value from the AbstractState using an integer value for the desired output value
# param handle The integer handle for the state class stored in memory
# param param The integer value for the parameter you want
function AbstractState_keyed_output(handle::Clong, param::Clong)
  output = ccall( (:AbstractState_keyed_output, "CoolProp"), Cdouble, (Clong,Clong,Ref{Clong},Ptr{UInt8},Clong), handle,param,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  elseif output == -Inf
    error("CoolProp: no correct state has been set with AbstractState_update")
  end
  return output
end
function AbstractState_output(handle::Clong, param::AbstractString)
  return AbstractState_keyed_output(handle, get_param_index(param))
end

# Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.
# handle The integer handle for the state class stored in memory
# input_pair The integer value for the input pair obtained from get_input_pair_index
# value1 The pointer to the array of the first input parameters
# value2 The pointer to the array of the second input parameters
# length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
# T The pointer to the array of temperature
# p The pointer to the array of pressure
# rhomolar The pointer to the array of molar density
# hmolar The pointer to the array of molar enthalpy
# smolar The pointer to the array of molar entropy
function AbstractState_update_and_common_out{R1<:Real,R2<:Real,F<:AbstractFloat}(handle::Clong, input_pair::Clong, value1::Array{R1}, value2::Array{R2}, length::Real, T::Array{F}, p::Array{F}, rhomolar::Array{F}, hmolar::Array{F}, smolar::Array{F})
  ccall( (:AbstractState_update_and_common_out, "CoolProp"), Void, (Clong,Clong,Ref{Cdouble},Ref{Cdouble},Clong,Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Clong},Ptr{UInt8},Clong), handle,input_pair,value1,value2,length,T,p,rhomolar,hmolar,smolar,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end
function AbstractState_update_and_common_out{R1<:Real,R2<:Real,F<:AbstractFloat}(handle::Clong, input_pair::AbstractString, value1::Array{R1}, value2::Array{R2}, length::Real, T::Array{F}, p::Array{F}, rhomolar::Array{F}, hmolar::Array{F}, smolar::Array{F})
  AbstractState_update_and_common_out(handle, get_input_pair_index(input_pair), value1, value2, length, T, p, rhomolar, hmolar, smolar)
  return nothing
end

# Update the state of the AbstractState and get one output value (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.
# handle The integer handle for the state class stored in memory
# input_pair The integer value for the input pair obtained from get_input_pair_index
# value1 The pointer to the array of the first input parameters
# value2 The pointer to the array of the second input parameters
# length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
# output The indice for the output desired
# out The pointer to the array for output
function AbstractState_update_and_1_out{R1<:Real,R2<:Real,F<:AbstractFloat}(handle::Clong, input_pair::Clong, value1::Array{R1}, value2::Array{R2}, length::Real, output::Clong, out::Array{F})
  ccall( (:AbstractState_update_and_1_out, "CoolProp"), Void, (Clong,Clong,Ref{Cdouble},Ref{Cdouble},Clong,Clong,Ref{Cdouble},Ref{Clong},Ptr{UInt8},Clong), handle,input_pair,value1,value2,length,output,out,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end
function AbstractState_update_and_1_out{R1<:Real,R2<:Real,F<:AbstractFloat}(handle::Clong, input_pair::AbstractString, value1::Array{R1}, value2::Array{R2}, length::Real, output::AbstractString, out::Array{F})
  AbstractState_update_and_1_out(handle, get_input_pair_index(input_pair), value1, value2, length, get_param_index(output), out)
  return nothing
end

# Update the state of the AbstractState and get an output value five common outputs (temperature, pressure, molar density, molar enthalpy and molar entropy) from the AbstractState using pointers as inputs and output to allow array computation.
# handle The integer handle for the state class stored in memory
# input_pair The integer value for the input pair obtained from get_input_pair_index
# value1 The pointer to the array of the first input parameters
# value2 The pointer to the array of the second input parameters
# length The number of elements stored in the arrays (both inputs and outputs MUST be the same length)
# outputs The 5-element vector of indices for the outputs desired
# out1 The pointer to the array for the first output
# out2 The pointer to the array for the second output
# out3 The pointer to the array for the third output
# out4 The pointer to the array for the fourth output
# out5 The pointer to the array for the fifth output
function AbstractState_update_and_5_out{R1<:Real,R2<:Real,F<:AbstractFloat}(handle::Clong, input_pair::Clong, value1::Array{R1}, value2::Array{R2}, length::Real, outputs::Array{Clong}, out1::Array{F}, out2::Array{F}, out3::Array{F}, out4::Array{F}, out5::Array{F})
  ccall( (:AbstractState_update_and_5_out, "CoolProp"), Void, (Clong,Clong,Ref{Cdouble},Ref{Cdouble},Clong,Ref{Clong},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Cdouble},Ref{Clong},Ptr{UInt8},Clong), handle,input_pair,value1,value2,length,outputs,out1,out2,out3,out4,out5,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end
function AbstractState_update_and_5_out{R1<:Real,R2<:Real,S<:AbstractString,F<:AbstractFloat}(handle::Clong, input_pair::AbstractString, value1::Array{R1}, value2::Array{R2}, length::Real, outputs::Array{S}, out1::Array{F}, out2::Array{F}, out3::Array{F}, out4::Array{F}, out5::Array{F})
  outputs_key = Array(Clong,5)
  for k = 1:5
    outputs_key[k] = get_param_index(outputs[k])
  end
  AbstractState_update_and_5_out(handle, get_input_pair_index(input_pair), value1, value2, length, outputs_key, out1, out2, out3, out4, out5)
  return nothing
end

# Set binary interraction parrameter for mixtures
# handle The integer handle for the state class stored in memory
# i indice of the first fluid of the binary pair
# j indice of the second fluid of the binary pair
# parameter string wit the name of the parameter
# value the value of the binary interaction parameter
function AbstractState_set_binary_interaction_double(handle::Clong,i::Int, j::Int, parameter::AbstractString, value::Cdouble)
  ccall( (:AbstractState_set_binary_interaction_double, "CoolProp"), Void, (Clong,Clong,Clong,Ptr{UInt8},Cdouble,Ref{Clong},Ptr{UInt8},Clong), handle,i,j,parameter,value,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

# Set some fluid parameter (ie volume translation for cubic)
# handle The integer handle for the state class stored in memory
# i indice of the fluid the parramter should be applied too (for mixtures)
# parameter string wit the name of the parameter
# value the value of the parameter
function AbstractState_set_fluid_parameter_double(handle::Clong, i::Int, parameter::AbstractString, value::Cdouble)
  ccall( (:AbstractState_set_fluid_parameter_double, "CoolProp"), Void, (Clong,Clong,Ptr{UInt8},Cdouble,Ref{Clong},Ptr{UInt8},Clong), handle,i,parameter,value,errcode,message_buffer::Array{UInt8,1},buffer_length)
  if errcode[] != 0
    if errcode[] == 1
      error("CoolProp: ", bytestring(convert(Ptr{UInt8}, pointer(message_buffer))))
    elseif errcode[] == 2
      error("CoolProp: message buffer too small")
    else # == 3
      error("CoolProp: unknown error")
    end
  end
  return nothing
end

end #module
