import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import itertools, numpy as np

variables = ['T', 'P', 'H', 'S', 'U', 'D']

for fluid in ['Water', 'CO2', 'n-Propane', 'MDM']:
    T = np.linspace(CP.PropsSI(fluid, 'Tmin'), CP.PropsSI(fluid, 'T_critical') - 1e-5, 1000)

    fig = plt.figure(1, figsize=(10, 10), dpi=100)
    for i, types in enumerate(itertools.combinations(variables, 2)):
        ax = plt.subplot(5, 3, i + 1)

        types = list(types)

        if types[0] in ['T', 'P'] and types != ['T', 'P']:
            types[0], types[1] = types[1], types[0]

        xL = CP.PropsSI(types[0], 'T', T, 'Q', 0, fluid)
        yL = CP.PropsSI(types[1], 'T', T, 'Q', 0, fluid)
        xV = CP.PropsSI(types[0], 'T', T, 'Q', 1, fluid)
        yV = CP.PropsSI(types[1], 'T', T, 'Q', 1, fluid)
        Tc = CP.PropsSI(fluid, 'T_critical')
        xc = CP.PropsSI(types[0], 'T', Tc - 1e-6, 'Q', 1, fluid)
        yc = CP.PropsSI(types[1], 'T', Tc - 1e-6, 'Q', 1, fluid)

        ax.plot(xL, yL, 'k')
        ax.plot(xV, yV, 'k')
        ax.plot(xc, yc, 'o')

        if types[0] in ['P', 'D']:
            ax.set_xscale('log')
        if types[1] in ['P', 'D']:
            ax.set_yscale('log')

        plt.title(' '.join(types))
        plt.xlabel(types[0])
        plt.ylabel(types[1])

    plt.tight_layout()
    plt.savefig('AllFlash' + fluid + '.pdf')
    plt.show()
