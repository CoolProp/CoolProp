module CoolProp

export F2K, K2F, HAPropsSI, PropsSI, PhaseSI, get_global_param_string, get_param_index

function F2K(TF::Number)
  return ccall( (:F2K, "CoolProp"), Cdouble, (Cdouble,), TF)
end

function K2F(TK::Number)
  return ccall( (:K2F, "CoolProp"), Cdouble, (Cdouble,), TK)
end

function HAPropsSI(Output::String, Name1::String, Value1::Number, Name2::String, Value2::Number, Name3::String, Value3::Number)
  val = ccall( (:HAPropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8},Float64), Output,Name1,Value1,Name2,Value2,Name3,Value3)
  if val == Inf
    error("CoolProp:", get_global_param_string("errstring"))
  end
  return val
end

function PropsSI(Output::String, Name1::String, Value1::Number, Name2::String, Value2::Number, Fluid::String)
  val = ccall( (:PropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}), Output,Name1,Value1,Name2,Value2,Fluid)
  if val == Inf
    error("CoolProp:", get_global_param_string("errstring"))
  end
  return val
end

function PhaseSI(Name1::String, Value1::Number, Name2::String, Value2::Number, Fluid::String)
  outstring = Array(UInt8, 255)
  val = ccall( (:PhaseSI, "CoolProp"), Int32, (Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}, Ptr{UInt8}, Int), Name1,Value1,Name2,Value2,Fluid,outstring,length(outstring))
  return bytestring(convert(Ptr{UInt8}, pointer(outstring)))
end

# This function returns the output string in pre-allocated char buffer.  If buffer is not large enough, no copy is made
function get_global_param_string(Key::String)
  Outstring = Array(UInt8, 255)
  val = ccall( (:get_global_param_string, "CoolProp"), Clong, (Ptr{UInt8},Ptr{UInt8},Int), Key,Outstring,length(Outstring))
  return bytestring(convert(Ptr{UInt8}, pointer(Outstring)))
end

# Get the index for a parameter "T", "P", etc.
# returns the index as a long.  If input is invalid, returns -1
function get_param_index(Param::String)
  return ccall( (:get_param_index, "CoolProp"), Clong, (Ptr{UInt8},), Param)
end

end #module
