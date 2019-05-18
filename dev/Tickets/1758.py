import numpy as np
import random
import CoolProp.CoolProp as CP
import time

random.seed("coolprop_test")
p = 101325  # 1 atmosphere
T = np.random.uniform(120, 400, 10000) + 273.15  # Random points from 120 to 400 deg C, gas phase only

# Make sure the objects exist and create tables if needed
normal_state = CP.AbstractState("HEOS", "H2O")
tabular_state = CP.AbstractState("BICUBIC&HEOS", "H2O")

# Measure execution speed
results = {}

tmp = time.time()
for Ti in T:
    rho = CP.PropsSI("Dmass", "P", p, "T", Ti, "H2O")
results["1. PropsSI"] = time.time() - tmp

tmp = time.time()
for Ti in T:
    normal_state.update(CP.PT_INPUTS, p, Ti)
    rho = normal_state.keyed_output(CP.iDmass)
results["2. HEOS"] = time.time() - tmp

tmp = time.time()
for Ti in T:
    tabular_state.update(CP.PT_INPUTS, p, Ti)
    rho = tabular_state.keyed_output(CP.iDmass)
results["3. Tables"] = time.time() - tmp

# for k in sorted(results):
#    print("{0} : {1} ms".format(k, results[k]*1e3))
#print("\nDo NOT do this!")
tmp = time.time()
for Ti in T:
    normal_state = CP.AbstractState("HEOS", "H2O")
    normal_state.update(CP.PT_INPUTS, p, Ti)
    rho = normal_state.keyed_output(CP.iDmass)
results["4. HEOS (create state)"] = time.time() - tmp

tmp = time.time()
for Ti in T:
    tabular_state = CP.AbstractState("BICUBIC&HEOS", "H2O")
    tabular_state.update(CP.PT_INPUTS, p, Ti)
    rho = tabular_state.keyed_output(CP.iDmass)
results["5. Tables (create state)"] = time.time() - tmp

for k in sorted(results):
    print("{0} : {1} ms".format(k, results[k] * 1e3))
