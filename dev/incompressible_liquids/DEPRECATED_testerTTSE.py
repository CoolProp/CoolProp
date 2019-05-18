import CoolProp.CoolProp as CP
fluid = 'Propane'
print(CP.enable_TTSE_LUT(fluid))
print(CP.isenabled_TTSE_LUT(fluid))
print(CP.Props('H', 'P', 300, 'Q', 0, fluid))
print(CP.Props('H', 'P', 310, 'Q', 0, fluid))
print(CP.Props('H', 'P', 315, 'Q', 0, fluid))
#
fluid = 'TestSolution-0.3'
print(CP.enable_TTSE_LUT(fluid))
print(CP.isenabled_TTSE_LUT(fluid))
print(CP.Props('H', 'P', 3000, 'T', 280, fluid))
