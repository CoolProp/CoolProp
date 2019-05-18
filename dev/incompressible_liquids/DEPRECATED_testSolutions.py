#!/usr/bin/python
import sys
import CoolProp.CoolProp as cp
import numpy as np


def props(in1=None, in2=None, in3=None, in4=None, in5=None, in6=None):
    try:
        return cp.PropsU(in1=in1, in2=in2, in3=in3, in4=in4, in5=in5, in6=in6, in7="SI")
    except ValueError as ve:
        # print "Error in CoolProp, try adjusting the fluid or T and p:"
        print(ve)
        return -1.0 * np.NAN


# print "{0:14.8f}".format(CP.Props('V','D',13,'P',500,'n-Pentane'))
# print "{0:14.8f}".format(CP.Props('V','H',158,'P',1000,'TX22'))
#T = 300
T = float(sys.argv[1]) + 273.15
P = float(sys.argv[2]) * 1e5
print("Temperature: " + str(T - 273.15) + " C")
print("Pressure:    " + str(P / 1e5) + " bar")
print("")
Melinder = ["MEG", "MPG", "MEA", "MMA", "MGL", "MAM", "MKC", "MCA", "MMG", "MNA", "MKA", "MKF", "MLI"]
SecCool = ["ZiAC", "IceEA", "IcePG", "IceNA", "PK2000"]
Other = ["LiBr"]

fluids = []
# fluids.extend(Melinder)
# fluids.extend(SecCool)
fluids.extend(Other)

for fluid in fluids:
    print("Fluid: " + str(fluid))
    try:
        print("Density:    " + "{0:14.8f} kg/m3  ".format(props('D', 'T', T, 'P', P, fluid + '-20%')))
        print("Heat cap.:  " + "{0:14.8f} kJ/kg/K".format(props('C', 'T', T, 'P', P, fluid + '-20%') / 1e3))
        print("Th. cond.:  " + "{0:14.8f} W/m/K  ".format(props('L', 'T', T, 'P', P, fluid + '-20%')))
        print("Dyn. visc.: " + "{0:14.8f} mPas   ".format(props('V', 'T', T, 'P', P, fluid + '-20%') * 1e3))
        print("Enthalpy:   " + "{0:14.8f} kJ/kg  ".format(props('H', 'T', T, 'P', P, fluid + '-20%') / 1e3))
        print("In. energy: " + "{0:14.8f} kJ/kg  ".format(props('U', 'T', T, 'P', P, fluid + '-20%') / 1e3))
        print("Entropy:    " + "{0:14.8f} kJ/kg/K".format(props('S', 'T', T, 'P', P, fluid + '-20%') / 1e3))
        print("Saturation: " + "{0:14.8f} bar    ".format(props('Psat', 'T', T, 'P', P, fluid + '-20%') / 1e5))
        print("Freezing:   " + "{0:14.8f} C      ".format(props('Tfreeze', 'T', T, 'P', P, fluid + '-20%') - 273.15))
    except ValueError as ve:
        print("Error in CoolProp, try adjusting T and p:")
        print(ve)
    print("")
