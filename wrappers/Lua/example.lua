local cp = require "coolprop"

print[[

Information Examples
====================
]]
print("Version", cp.version())
print("gitrevision", cp.gitrevision())
print("FluidsList", cp.FluidsList())
print("Debug Level", cp.get_debug_level())
print()
print[[

Usage Examples
==============
]]
print(cp.Props1SI("Water", "Phase"))
print(cp.PropsSI("C", "P", 101325, "T", 300, "Water"))
print(cp.PropsSI("d(Hmass)/d(T)|P", "P", 101325, "T", 300, "Water"))
print(cp.PropsSI("D", "P", 101325, "T", 300, "Air.mix"))
print(cp.PhaseSI("P", 101325, "Q", 0, "Water"))
print(cp.get_global_param_string("predefined_mixtures"))
print(cp.get_fluid_param_string("Water", "aliases"))
print(cp.get_param_index("T"))
print(cp.HAPropsSI('H','T',298.15,'P',101325,'R',0.5))
print(10, "Fahrenheits is", cp.F2K(10) , "Kelvins")
print(cp.F2K(10), "Kelvins is", cp.K2F(cp.F2K(10)), "Fahrenheits")
print()
print[[

Error Examples
==============
]]
print(cp.Props1SI("Error", "Phase"))
print(cp.PropsSI("Error", "T", 298.15, "P", 101325, "Nitrogen"))
print(cp.PhaseSI("Error", 101325, "Q", 0, "Water"))
print(cp.get_parameter_information_string("Error"))            -- What are the correct inputs for this?
print(cp.saturation_ancillary("Water", "C", 10, "Error", 300)) -- What are the correct inputs for this?
print(cp.HAPropsSI('Error','T',298.15,'P',101325,'R',0.5))
print(cp.get_param_index("Error"))
print(cp.get_fluid_param_string("Water", "wrong"))
print()