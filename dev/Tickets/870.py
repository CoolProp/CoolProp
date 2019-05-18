import numpy as np
import CoolProp
import matplotlib.pyplot as plt

AS = CoolProp.AbstractState('HEOS', 'Water')

steps = 100

# Saturated liquid
AS.update(CoolProp.PQ_INPUTS, 101325, 0)
Ts = AS.T()
h_fg = AS.saturated_vapor_keyed_output(CoolProp.iHmass) - AS.saturated_liquid_keyed_output(CoolProp.iHmass)
cl = AS.cpmass()

# Subcooled liquid
x, y = [], []
for T in np.linspace(Ts - 30, Ts - 0.1, steps):
    AS.update(CoolProp.PT_INPUTS, 101325, T)
    x.append(-cl * (Ts - T) / h_fg)
    y.append(AS.first_partial_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
plt.plot(x, y, label='Subcooled', color='gray')

# Two-phase derivatives (splined)
x, y1 = [], []
for Q in np.linspace(0, 0.3, steps):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    y1.append(AS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3))
plt.plot(x, y1, label='Two-phase (splined)')

# Two-phase derivatives (normal)
x, y1 = [], []
for Q in np.linspace(0.0, 0.6, steps):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    y1.append(AS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
plt.plot(x, y1, label='Two-phase (normal)')

AS = CoolProp.AbstractState('TTSE&HEOS', 'Water')

# Two-phase derivatives (splined, tabular)
x, y1 = [], []
for Q in np.linspace(0, 0.3, steps):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    try:
        y1.append(AS.first_two_phase_deriv_splined(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP, 0.3))
    except Exception as e:
        print(e)
        y1.append(np.NAN)
        break
plt.plot(x, y1, label='Two-phase (splined, tabular)', ls='--', lw=3)

# Two-phase derivatives (normal, tabular)
x, y1 = [], []
for Q in np.linspace(0.0, 0.6, steps):
    AS.update(CoolProp.PQ_INPUTS, 101325, Q)
    x.append(AS.Q())
    try:
        y1.append(AS.first_two_phase_deriv(CoolProp.iDmass, CoolProp.iHmass, CoolProp.iP))
    except Exception as e:
        print(e)
        y1.append(np.NAN)
        break
plt.plot(x, y1, label='Two-phase (normal, tabular)', ls='--', lw=3)

plt.title(r'$d\rho/dh|p$')
plt.xlabel('vapor quality (-)')
plt.ylabel(r'$d\rho/dh|p$')
plt.ylim(-0.005, 0.005)
plt.legend(loc='best')
plt.savefig('870_test_two_phase_tabular.png')
plt.show()
