import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import CoolProp
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from CoolProp.CoolProp import FluidsList

print("Testing TTSE for revision " + CoolProp.__gitrevision__)


def compareProperty(fluid="", p=0, what=""):
    global c_diff, c_unit, c_exce

    if p == 0:
        p = 0.75 * CP.Props(fluid, "pcrit")

    Delta_T = 50
    T_bub = CP.Props("T", "P", p, "Q", 0, fluid)
    T_dew = CP.Props("T", "P", p, "Q", 1, fluid)
    h_bub = CP.Props("H", "P", p, "Q", 0, fluid)
    h_dew = CP.Props("H", "P", p, "Q", 1, fluid)

    T_1 = T_bub - Delta_T
    if T_1 < CP.Props(fluid, "Tmin"):
        T_1 = CP.Props(fluid, "Tmin") + 0.5 * (T_bub - CP.Props(fluid, "Tmin"))
    T_2 = T_dew + Delta_T

    h_1 = CP.Props("H", "P", p, "T", T_1, fluid)
    h_2 = CP.Props("H", "P", p, "T", T_2, fluid)

    h_liq = numpy.linspace(h_1, h_bub, num=100)
    h_vap = numpy.linspace(h_dew, h_2, num=100)

    T_liq = CP.Props("T", "P", p, "H", h_liq, fluid) - T_bub
    T_vap = CP.Props("T", "P", p, "H", h_vap, fluid) - T_dew

    X_liq_STDV = CP.Props(what, "P", p, "H", h_liq, fluid)
    X_vap_STDV = CP.Props(what, "P", p, "H", h_vap, fluid)

    CP.enable_TTSE_LUT(fluid)
    X_liq_TTSE = CP.Props(what, "P", p, "H", h_liq, fluid)
    X_vap_TTSE = CP.Props(what, "P", p, "H", h_vap, fluid)

    if numpy.max([X_liq_STDV / X_liq_TTSE, X_vap_STDV / X_vap_TTSE]) > 1.25 or numpy.min([X_liq_STDV / X_liq_TTSE, X_vap_STDV / X_vap_TTSE]) < 0.75:
        c_diff += 1
        print("")
        print("There were problems with " + what + " for " + fluid)
        print("Relative difference liquid: " + str(numpy.mean((X_liq_STDV - X_liq_TTSE) / X_liq_STDV)))
        print("Relative difference vapour: " + str(numpy.mean((X_vap_STDV - X_vap_TTSE) / X_vap_STDV)))
        print("Average factor liquid: " + str(numpy.mean(X_liq_STDV / X_liq_TTSE)))
        print("Average factor vapour: " + str(numpy.mean(X_vap_STDV / X_vap_TTSE)))
        #plt.plot(numpy.append(T_liq,T_vap),numpy.append(X_liq_STDV,X_vap_STDV),label=what+", standard")
        #plt.plot(numpy.append(T_liq,T_vap),numpy.append(X_liq_TTSE,X_vap_TTSE),label=what+", TTSE")
        # plt.show(block=True)
        # plt.savefig("/home/jowr/tmp/viscosity/"+fluid+".png")
        # plt.clf()


#toTest  = ["L","V"]
toTest = "V"
#fluids = ["n-Pentane","R134a"]
fluids = FluidsList()
c_diff = 0
c_unit = 0
c_exce = 0

for fluid in fluids:
    try:
        compareProperty(fluid=fluid, what=toTest)
    except ValueError:
        c_exce += 1
        print("An exception occurred for " + toTest + " with " + fluid)

print("Finished testing TTSE:")
print("Errors occurred in " + str(c_exce) + " out of " + str(len(fluids)) + " fluids")
print("and differences occurred in " + str(c_diff) + " fluids.")
