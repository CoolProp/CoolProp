import numpy as np
import CoolProp
import matplotlib.pyplot as plt

AS = CoolProp.AbstractState('HEOS','Water')

# Saturated liquid
AS.update(CoolProp.PQ_INPUTS, 101325, 0)
Ts = AS.T()
h_fg = AS.saturated_vapor_keyed_output(CoolProp.iHmass) - AS.saturated_liquid_keyed_output(CoolProp.iHmass)
cl = AS.cpmass()

# Subcooled liquid
x, y = [], []
for T in np.linspace(Ts - 30, Ts - 0.1, 1000):
    AS.update(CoolProp.PT_INPUTS, 101325, T)
    x.append(-cl*(Ts-T)/h_fg)
    y.append(AS.first_partial_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
plt.plot(x, y, label = 'Subcooled')

# Two-phase derivatives (splined)
x, y1 = [], []
for Q in np.linspace(0, 0.3, 1000):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    y1.append(AS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3))
plt.plot(x, y1, label = 'Two-phase (splined)')

# Two-phase derivatives (normal)
x, y1 = [], []
for Q in np.linspace(0.0, 0.6, 1000):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    y1.append(AS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
plt.plot(x, y1, label = 'Two-phase')

AS = CoolProp.AbstractState('TTSE&HEOS','Water')
# Two-phase derivatives (normal, tabular)
x, y1 = [], []
for Q in np.linspace(0.0, 0.6, 1000):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    y1.append(AS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
plt.plot(x, y1, label = 'Two-phase (tabular)', lw = 4)

plt.title(r'$d\rho/dh|p$')
plt.xlabel('vapor quality (-)')
plt.ylabel(r'$d\rho/dh|p$')
plt.ylim(-0.005, 0.005)
plt.legend(loc='best')
plt.savefig('1238_test_two_phase_tabular.png')
plt.show()