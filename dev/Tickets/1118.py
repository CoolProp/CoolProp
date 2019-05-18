
from __future__ import print_function, division, absolute_import

from CoolProp import AbstractState

#fluid = "R407C"
#fluid = "R32[0.381109419953993]&R125[0.179558888662016]&R134A[0.439331691383991]"
#fluid = "R32[0.40]&R125[0.15]&R134A[0.45]"

fluids = ["R32", "R125", "R134a"]
moles = [0.381109419953993, 0.179558888662016, 0.439331691383991]

for backend in ["HEOS", "REFPROP"]:
    state = AbstractState(backend, '&'.join(fluids))
    state.set_mole_fractions(moles)
    print(state.get_mole_fractions())
    print(state.get_mass_fractions())
    try:
        print("Normal routines")
        print(state.T_critical())
        print(state.p_critical())
        print(True)
    except:
        try:
            print("Routine for all points")
            states = state.all_critical_points()
            for crit_state in states:
                print(crit_state.T)
                print(crit_state.p)
                print(crit_state.stable)
        except:
            print("All failed")

print("\n-----------------\n")

# import numpy as np
# import CoolProp
# from CoolProp.CoolProp import PyCriticalState
# def get_critical_point(state):
#     crit_state = PyCriticalState()
#     crit_state.T = np.nan
#     crit_state.p = np.nan
#     crit_state.rhomolar = np.nan
#     crit_state.stable = False
#     try:
#         crit_state.T = state.T_critical()
#         crit_state.p = state.p_critical()
#         crit_state.rhomolar = state.rhomolar_critical()
#         crit_state.stable = True
#     except:
#         try:
#             for crit_state_tmp in state.all_critical_points():
#                 if crit_state_tmp.stable and (crit_state_tmp.T > crit_state.T or not np.isfinite(crit_state.T)):
#                     crit_state.T = crit_state_tmp.T
#                     crit_state.p = crit_state_tmp.p
#                     crit_state.rhomolar = crit_state_tmp.rhomolar
#                     crit_state.stable = crit_state_tmp.stable
#         except:
#             raise ValueError("Could not calculate the critical point data.")
#     new_state = AbstractState(state.backend_name(), '&'.join(state.fluid_names()))
#     masses = state.get_mass_fractions()
#     if len(masses)>1:
#         new_state.set_mass_fractions(masses) # Uses mass fraction to work with incompressibles
#     if np.isfinite(crit_state.p) and np.isfinite(crit_state.T):
#         new_state.specify_phase(CoolProp.iphase_critical_point)
#         new_state.update(CoolProp.PT_INPUTS, crit_state.p, crit_state.T)
#         #new_state.update(CoolProp.DmolarT_INPUTS, crit_state.rhomolar, crit_state.T)
#         return new_state
#     raise ValueError("Could not calculate the critical point data.")

from CoolProp.Plots.Common import process_fluid_state, get_critical_point
state = process_fluid_state("HEOS::R32[0.381109419953993]&R125[0.179558888662016]&R134A[0.439331691383991]")
cstate = get_critical_point(state)
print(cstate.T())


import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot

fluid = "R407C"
fluid = "R32[0.381109419953993]&R125[0.179558888662016]&R134A[0.439331691383991]"
#fluid = "R32[0.40]&R125[0.15]&R134A[0.45]"

plot_RP = PropertyPlot("REFPROP::" + fluid, 'PH', unit_system='EUR', tp_limits='ACHP')
plot_RP.calc_isolines()
plot_RP.draw()

plot_CP = PropertyPlot("HEOS::" + fluid, 'PH', unit_system='EUR', tp_limits='ACHP')
plot_CP.calc_isolines()
plot_CP.draw()

plt.show(block=True)
