import json
import time

import CoolProp.CoolProp as CP
import numpy as np

if __name__ == "__main__":
    n = int(1e6)

    FLUID = "MEG"

    fluid_json = json.load(open(f"dev/incompressible_liquids/json/{FLUID}.json"))  # noqa: SIM115
    fluid_json["name"] = fluid_json["name"] + "_mod"
    fluid_json.pop("density_cheb", None)
    fluid_json.pop("specific_heat_cheb", None)
    CP.add_fluids_as_JSON("INCOMP", json.dumps([fluid_json]))

    fluid = CP.AbstractState("INCOMP", f"{FLUID}_mod")
    fluid.set_mass_fractions([0.35])

    print("\nStarting CoolProp polynomial tests")
    print("-" * 40)

    T = np.linspace(300, 350, n)
    fluid.update(CP.PT_INPUTS, 1e5, T[0])
    h1, s1 = fluid.hmass(), fluid.smass()
    fluid.update(CP.PT_INPUTS, 1e5, T[-1])
    h2, s2 = fluid.hmass(), fluid.smass()
    h = np.linspace(h1, h2, n)
    s = np.linspace(s1, s2, n)

    Tfreeze = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 0.0, _T)
        Tfreeze.append(fluid.keyed_output(CP.iT_freeze))
    end_time = time.time()
    print(f"Tfreeze(T,p) elapsed time: {end_time - start_time} seconds")

    mu = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        mu.append(fluid.viscosity())
    end_time = time.time()
    print(f"mu(T,p) elapsed time: {end_time - start_time} seconds")

    k = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        k.append(fluid.conductivity())
    end_time = time.time()
    print(f"k(T,p) elapsed time: {end_time - start_time} seconds")

    d = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        d.append(fluid.rhomass())
    end_time = time.time()
    print(f"d(T,p) elapsed time: {end_time - start_time} seconds")

    cp = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        cp.append(fluid.cpmass())
    end_time = time.time()
    print(f"cp(T,p) elapsed time: {end_time - start_time} seconds")

    Touth = []
    start_time = time.time()
    for _h in h:
        fluid.update(CP.HmassP_INPUTS, _h, 1e5)
        Touth.append(fluid.T())
    end_time = time.time()
    print(f"T(p,h) elapsed time: {end_time - start_time} seconds")

    Touts = []
    start_time = time.time()
    for _s in s:
        fluid.update(CP.PSmass_INPUTS, 1e5, _s)
        Touts.append(fluid.T())
    end_time = time.time()
    print(f"T(p,s) elapsed time: {end_time - start_time} seconds")

    ###########################################################################################################

    print("-" * 40)
    print("\nStarting CoolProp chebychev tests")
    print("-" * 40)

    fluid = CP.AbstractState("INCOMP", f"{FLUID}")
    fluid.set_mass_fractions([0.35])

    T = np.linspace(300, 350, n)
    fluid.update(CP.PT_INPUTS, 1e5, T[0])
    h1, s1 = fluid.hmass(), fluid.smass()
    fluid.update(CP.PT_INPUTS, 1e5, T[-1])
    h2, s2 = fluid.hmass(), fluid.smass()
    h = np.linspace(h1, h2, n)
    s = np.linspace(s1, s2, n)

    # psat = []
    # start_time = time.time()
    # for _T in T:
    #     fluid.update(CP.QT_INPUTS, 0.0, _T)
    #     psat.append(fluid.p())
    # end_time = time.time()
    # print(f"psat(Q,T) elapsed time: {end_time - start_time} seconds")

    Tfreeze = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 0.0, _T)
        Tfreeze.append(fluid.keyed_output(CP.iT_freeze))
    end_time = time.time()
    print(f"Tfreeze(T,p) elapsed time: {end_time - start_time} seconds")

    mu = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        mu.append(fluid.viscosity())
    end_time = time.time()
    print(f"mu(T,p) elapsed time: {end_time - start_time} seconds")

    k = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        k.append(fluid.conductivity())
    end_time = time.time()
    print(f"k(T,p) elapsed time: {end_time - start_time} seconds")

    d = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        d.append(fluid.rhomass())
    end_time = time.time()
    print(f"d(T,p) elapsed time: {end_time - start_time} seconds")

    cp = []
    start_time = time.time()
    for _T in T:
        fluid.update(CP.PT_INPUTS, 1e5, _T)
        cp.append(fluid.cpmass())
    end_time = time.time()
    print(f"cp(T,p) elapsed time: {end_time - start_time} seconds")

    Touth = []
    start_time = time.time()
    for _h in h:
        fluid.update(CP.HmassP_INPUTS, _h, 1e5)
        Touth.append(fluid.T())
    end_time = time.time()
    print(f"T(p,h) elapsed time: {end_time - start_time} seconds")

    Touts = []
    start_time = time.time()
    for _s in s:
        fluid.update(CP.PSmass_INPUTS, 1e5, _s)
        Touts.append(fluid.T())
    end_time = time.time()
    print(f"T(p,s) elapsed time: {end_time - start_time} seconds")

    print("-" * 40)
    print("\n")
