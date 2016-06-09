
from __future__ import print_function,division,absolute_import

from CoolProp import AbstractState

#fluid = "R407C"
#fluid = "R32[0.381109419953993]&R125[0.179558888662016]&R134A[0.439331691383991]"
#fluid = "R32[0.40]&R125[0.15]&R134A[0.45]"

fluids = ["R32","R125","R134a"]
moles = [0.381109419953993,0.179558888662016,0.439331691383991]

for backend in ["HEOS","REFPROP"]:
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


import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot

fluid = "R407C"
#fluid = "R32[0.381109419953993]&R125[0.179558888662016]&R134A[0.439331691383991]"
#fluid = "R32[0.40]&R125[0.15]&R134A[0.45]"

plot_RP = PropertyPlot("REFPROP::"+fluid, 'PH', unit_system='EUR', tp_limits='ACHP')
plot_RP.calc_isolines()
plot_RP.draw()

plot_CP = PropertyPlot("HEOS::"+fluid, 'PH', unit_system='EUR', tp_limits='ACHP')
plot_CP.calc_isolines()
plot_CP.draw()

plt.show(block=False)