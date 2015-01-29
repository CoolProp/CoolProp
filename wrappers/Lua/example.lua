local c = require "coolprop"

print(c.PropsSI("XX", "T", 298.15, "P", 101325, "Nitrogen"))
print(c.PropsSI("C", "P", 101325, "T", 300, "Water"))
print(c.PropsSI("d(Hmass)/d(T)|P", "P", 101325, "T", 300, "Water"))
print(c.PropsSI("D", "P", 101325, "T", 300, "Air.mix"))
print(c.PhaseSI("PXX", 101325, "Q", 0, "Water"))
print(c.get_global_param_string("predefined_mixtures"))
print(c.get_param_index("T"))
print(c.HAPropsSI('H','T',298.15,'P',101325,'R',0.5))
print(c.HAPropsSI('HXXX','T',298.15,'P',101325,'R',0.5))
print(c.F2K(10))
print(c.K2F(c.F2K(10)))