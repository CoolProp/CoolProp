import matplotlib
matplotlib.use('WXAgg')
from matplotlib import cm
import matplotlib.pyplot as plt

import numpy as np
import CoolProp
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(2, 2))
ax = fig.add_subplot(111, projection='3d')

NT = 1000
NR = 1000
rho, t = np.logspace(np.log10(2e-3), np.log10(1100), NR), np.linspace(275.15, 700, NT)
RHO, T = np.meshgrid(rho, t)

P = CoolProp.CoolProp.PropsSI('P', 'D', RHO.reshape((NR * NT, 1)), 'T', T.reshape((NR * NT, 1)), 'REFPROP-Water').reshape(NT, NR)

Tsat = np.linspace(273.17, 647.0, 100)
psat = CoolProp.CoolProp.PropsSI('P', 'Q', 0, 'T', Tsat, 'Water')
rhoL = CoolProp.CoolProp.PropsSI('D', 'Q', 0, 'T', Tsat, 'Water')
rhoV = CoolProp.CoolProp.PropsSI('D', 'Q', 1, 'T', Tsat, 'Water')

ax.plot_surface(np.log(RHO), T, np.log(P), cmap=cm.jet, edgecolor='none')
ax.plot(np.log(rhoL), Tsat, np.log(psat), color='k', lw=2)
ax.plot(np.log(rhoV), Tsat, np.log(psat), color='k', lw=2)

ax.text(0.3, 800, 22, "CoolProp", size=12)
ax.set_frame_on(False)
ax.set_axis_off()
ax.view_init(22, -136)
ax.set_xlabel(r'$\ln\rho$ ')
ax.set_ylabel('$T$')
ax.set_zlabel('$p$')
plt.tight_layout()
plt.savefig('_static/PVTCP.png', transparent=True)
plt.savefig('_static/PVTCP.pdf', transparent=True)
plt.close()
