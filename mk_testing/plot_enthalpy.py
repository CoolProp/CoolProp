import json

import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

T_C_orig = np.linspace(10, 100, 10)
pressures = {
    # "0 bar": 0,
    "1 bar": 101325,
    "6 bar": 6e5,
    "16 bar": 1.6e6,
    "100 bar": 1e7,
}

water = json.load(open("dev/incompressible_liquids/json/Water.json"))  # noqa: SIM115
CP.add_fluids_as_JSON("INCOMP", json.dumps([water]))

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
ax = axes[0]
cp = [CP.PropsSI("C", "T", T + 273.15, "P", 101325, "INCOMP::Water") for T in T_C_orig]
print(cp)

water800 = json.load(open("dev/incompressible_liquids/json/Water800.json"))  # noqa: SIM115
water800["name"] += "800"

CP.add_fluids_as_JSON("INCOMP", json.dumps([water800]))
cp = [
    CP.PropsSI("C", "T", T + 273.15, "P", 101325, "INCOMP::Water800") for T in T_C_orig
]
print(cp)


ax.plot(T_C_orig, cp, ".", label="cp")
for label, p in pressures.items():
    dhdT = [
        CP.PropsSI("H", "T", T + 273.15, "P", p, "INCOMP::Water")
        - CP.PropsSI("H", "T", T - 1.0 + 273.15, "P", p, "INCOMP::Water")
        for T in T_C_orig
    ]
    ax.plot(T_C_orig, dhdT, label="dhdT " + label)
ax.set_xlabel("Temperature [\u00b0C]")
ax.set_ylabel("Enthalpy [J/kg]")
ax.set_title("h(P)-h(P0) vs. temperature")
ax.legend()
ax.grid(True)

ax = axes[1]
for label, p in pressures.items():
    h_water = [
        CP.PropsSI("H", "T", T + 273.15, "P", p, "INCOMP::Water") for T in T_C_orig
    ]
    ax.plot(T_C_orig, h_water, label=label)
ax.set_xlabel("Temperature [\u00b0C]")
ax.set_ylabel("Enthalpy [J/kg]")
ax.set_title("Water enthalpy vs. temperature")
ax.legend()
ax.grid(True)

plt.show(block=False)

print("test")
