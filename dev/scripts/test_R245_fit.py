from CoolProp import CoolProp as CP
from PDSim.misc.datatypes import Collector
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fluid = 'R245fa'
Ttriple = CP.Props(fluid, 'Ttriple')
Tcrit = CP.Props(fluid, 'Tcrit')

RHO, TTT, RHO0, TTT0, ERR = Collector(), Collector(), Collector(), Collector(), Collector()

rhomax = CP.Props('D', 'T', Ttriple, 'Q', 0, 'R245fa')
# Build a database of "experimental" data
for T in np.linspace(Ttriple, Tcrit + 50, 80):
    for rho in np.linspace(1e-10, rhomax, 80):
        if (T > Tcrit or rho > CP.rhosatL_anc(fluid, T) or rho < CP.rhosatV_anc(fluid, T)):
            h = CP.Props('H', 'T', T, 'D', rho, 'R245fa')
            p = CP.Props('P', 'T', T, 'D', rho, 'R245fa')

            RHO << rho
            TTT << T
            ERR << h
            TTT0 << p

fig = plt.figure()
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(np.array(RHO.vec), np.array(TTT.vec), ERR.vec)
ax2 = fig.add_subplot(122, projection='3d')
ax2.scatter(np.array(RHO.vec), np.array(TTT.vec), TTT0.vec)
plt.show()
