
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.HumidAirProp import HAProps

Tdb = np.linspace(-10, 55, 100) + 273.15

# Make the figure and the axes
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes((0.15, 0.15, 0.8, 0.8))
ax.set_xlim(Tdb[0] - 273.15, Tdb[-1] - 273.15)
ax.set_ylim(0, 0.03)
ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")

# Saturation line
w = [HAProps('W', 'T', T, 'P', 101.325, 'R', 1.0) for T in Tdb]
ax.plot(Tdb - 273.15, w, lw=2)

# Humidity lines
RHValues = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
for RH in RHValues:
    w = [HAProps('W', 'T', T, 'P', 101.325, 'R', RH) for T in Tdb]
    ax.plot(Tdb - 273.15, w, 'r', lw=1)

# Humidity lines
for H in [-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
    # Line goes from saturation to zero humidity ratio for this enthalpy
    T1 = HAProps('T', 'H', H, 'P', 101.325, 'R', 1.0) - 273.15
    T0 = HAProps('T', 'H', H, 'P', 101.325, 'R', 0.0) - 273.15
    w1 = HAProps('W', 'H', H, 'P', 101.325, 'R', 1.0)
    w0 = HAProps('W', 'H', H, 'P', 101.325, 'R', 0.0)
    ax.plot(np.r_[T1, T0], np.r_[w1, w0], 'r', lw=1)

plt.show()
