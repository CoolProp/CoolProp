from __future__ import division

from CoolProp import CoolProp as CP

Fluid = 'CO2'
R = 8.314472 / CP.Props(Fluid, 'molemass')
Tc = CP.Props(Fluid, 'Tcrit')
pc = CP.Props(Fluid, 'pcrit')
w = CP.Props(Fluid, 'accentric')

a = 0.457235 * R**2 * Tc**2 / pc
b = 0.077796 * R * Tc / pc
kappa = 0.37464 + 1.54226 * w - 0.26992 * w**2

T = 298.15
rho = 1000
v = 1 / rho
Tr = T / Tc
alpha = (1 + kappa * (1 - Tr**0.5))**2

p = R * T / (v - b) - a * alpha / (v**2 + 2 * b * v - b**2)
print("%s %s" % (p, CP.Props('P', 'T', T, 'D', rho, Fluid)))
