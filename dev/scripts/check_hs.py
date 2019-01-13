from CoolProp.Plots.Plots import hs
import CoolProp
from CoolProp.CoolProp import Props
import matplotlib.pyplot as plt
import numpy as np

Fluid = 'Nitrogen'

fig = plt.figure()
ax = fig.add_subplot(111)

for Fluid in CoolProp.__fluids__:

    if Fluid == 'SES36':
        continue

    ax.cla()

    h_crit = Props('H', 'T', Props(Fluid, "Tcrit"), 'D', Props(Fluid, "rhocrit"), Fluid)
    s_crit = Props('S', 'T', Props(Fluid, "Tcrit"), 'D', Props(Fluid, "rhocrit"), Fluid)

    hL_Tmin = Props('H', 'T', Props(Fluid, "Tmin"), 'Q', 0, Fluid)
    hV_Tmin = Props('H', 'T', Props(Fluid, "Tmin"), 'Q', 1, Fluid)
    sL_Tmin = Props('S', 'T', Props(Fluid, "Tmin"), 'Q', 0, Fluid)
    sV_Tmin = Props('S', 'T', Props(Fluid, "Tmin"), 'Q', 1, Fluid)

    hs(Fluid, axis=ax)
    plt.plot(s_crit, h_crit, 'rd')
    plt.plot([sL_Tmin, sV_Tmin], [hL_Tmin, hV_Tmin], '--')
    plt.gca().axhline(h_crit)
    plt.gca().axhline(hV_Tmin)

    # Two-Phase
    for T in np.linspace(Props(Fluid, "Tmin") + 0.1, Props(Fluid, "Tcrit") - 1e-3, 30):
        for Q in np.linspace(0, 1, 30):
            try:
                h = Props("H", 'Q', Q, 'T', T, Fluid)
                s = Props("S", 'Q', Q, 'T', T, Fluid)
                T = Props("T", 'S', s, 'H', h, Fluid)
                # ax.plot(s,h,'o',mfc='none')
            except ValueError as VE:
                print(T, Q, '|||', '"T","S",', s, ',"H",', h, ',"' + Fluid + '"', '|||', VE)
                ax.plot(s, h, 'o', mfc='none')

    for h in np.linspace(hL_Tmin, hV_Tmin + 1500, 100):
        for s in np.linspace(sL_Tmin + 0.01, sV_Tmin, 100):
            try:
                h_pmax = Props('H', 'S', s, 'P', 6 * Props(Fluid, 'pcrit'), Fluid)
            except ValueError:
                h_pmax = 0
            htriple_s = (hV_Tmin - hL_Tmin) / (sV_Tmin - sL_Tmin) * (s - sL_Tmin) + hL_Tmin
            if h < htriple_s or h > h_pmax: continue
            try:
                T = Props("T", 'S', s, 'H', h, Fluid)
                # ax.plot(s,h,'o',mfc='none',ms=6)
            except ValueError:
                ax.plot(s, h, 's', mfc='none')


# if Fluid =='Propane':
##         ps = Props("P",'T',300,'Q',0,Fluid);
##         hL = Props("H",'Q',0,'T',300,Fluid);
##         sL = Props("S",'Q',0,'T',300,Fluid);
##         h = Props("H",'P',ps,'T',299.5,Fluid);
##         s = Props("S",'P',ps,'T',299.5,Fluid);
# print s,h,sL,hL
# plt.plot(s,h,'o')
# plt.plot(sL,hL,'d')
# plt.gca().axvline(s)
# plt.gca().axhline(268.75968916316691)

    fig.savefig('figs/' + Fluid + '.png', dpi=200)
    fig.savefig('figs/' + Fluid + '.pdf')
