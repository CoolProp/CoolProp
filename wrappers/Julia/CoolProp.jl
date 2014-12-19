module CoolProp

export F2K, PropsSI, PropsSI, PhaseSI

function F2K(TF::Float64)
  return ccall( (:F2K, "CoolProp"), Cdouble, (Cdouble,), TF)
end

function HAPropsSI(Output::String, Name1::String, Value1::Float64, Name2::String, Value2::Float64, Name3::String, Value3::Float64)
  return ccall( (:HAPropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8},Float64), Output,Name1,Value1,Name2,Value2,Name3,Value3)
end

function PropsSI(Output::String, Name1::String, Value1::Float64, Name2::String, Value2::Float64, Fluid::String)
  return ccall( (:PropsSI, "CoolProp"), Cdouble, (Ptr{UInt8},Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}), Output,Name1,Value1,Name2,Value2,Fluid)
end

function PhaseSI(Name1::String, Value1::Float64, Name2::String, Value2::Float64, Fluid::String)
  outstring = Array(UInt8, 255)
  val = ccall( (:PhaseSI, "CoolProp"), Int32, (Ptr{UInt8},Float64,Ptr{UInt8},Float64,Ptr{UInt8}, Ptr{UInt8}, Int), Name1,Value1,Name2,Value2,Fluid,outstring,length(outstring))
  return bytestring(convert(Ptr{UInt8}, outstring))
end

end #module