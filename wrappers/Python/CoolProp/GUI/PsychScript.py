
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.HumidAirProp import HAPropsSI

# All humid-air calls are SI (the v8 interface): pressure in Pa, enthalpy in
# J/kg of dry air, temperatures in K.  (The non-SI ``HAProps`` was removed in
# CoolProp v8; bd CoolProp-r9sq.27.)
p = 101325.0  # Pa

Tdb = np.linspace(-10, 55, 100) + 273.15

# Make the figure and the axes
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes((0.15, 0.15, 0.8, 0.8))
ax.set_xlim(Tdb[0] - 273.15, Tdb[-1] - 273.15)
ax.set_ylim(0, 0.03)
ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")

# Saturation line
w = [HAPropsSI('W', 'T', T, 'P', p, 'R', 1.0) for T in Tdb]
ax.plot(Tdb - 273.15, w, lw=2)

# Constant relative-humidity lines
RHValues = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for RH in RHValues:
    w = [HAPropsSI('W', 'T', T, 'P', p, 'R', RH) for T in Tdb]
    ax.plot(Tdb - 273.15, w, 'r', lw=1)

# Constant-enthalpy lines (enthalpy in J/kg of dry air)
for H in [-20e3, -10e3, 0, 10e3, 20e3, 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 90e3]:
    # Line goes from saturation to zero humidity ratio for this enthalpy
    T1 = HAPropsSI('T', 'H', H, 'P', p, 'R', 1.0) - 273.15
    T0 = HAPropsSI('T', 'H', H, 'P', p, 'R', 0.0) - 273.15
    w1 = HAPropsSI('W', 'H', H, 'P', p, 'R', 1.0)
    w0 = HAPropsSI('W', 'H', H, 'P', p, 'R', 0.0)
    ax.plot(np.r_[T1, T0], np.r_[w1, w0], 'r', lw=1)

plt.show()
