import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

HEOS = CP.AbstractState('HEOS','Methane&Ethane')

for x0 in [0.02, 0.2, 0.4, 0.6, 0.8, 0.98]:
    HEOS.set_mole_fractions([x0, 1 - x0])
    try:
        HEOS.build_phase_envelope("dummy")
        PE = HEOS.get_phase_envelope_data()
        plt.plot(PE.T, PE.p, '-')
    except ValueError as VE:
        print(VE)

plt.xlabel('Temperature [K]')
plt.ylabel('Pressure [Pa]')
plt.yscale('log')
plt.savefig('methane-ethane.pdf')
plt.savefig('methane-ethane.png')
plt.close()