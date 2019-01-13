from CoolProp.Plots import Ph
import CoolProp
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker
import numpy as np
import random

fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_axes((0.08, 0.1, 0.32, 0.83))
ax2 = fig.add_axes((0.50, 0.1, 0.32, 0.83))

Ref = 'R245fa'

BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS', Ref)
TTSE = CoolProp.AbstractState('TTSE&HEOS', Ref)
EOS = CoolProp.AbstractState('HEOS', Ref)
MM = EOS.molar_mass()
print(MM)

T = np.linspace(CP.PropsSI(Ref, 'Tmin') + 0.1, CP.PropsSI(Ref, 'Tcrit') - 0.01, 300)
pV = CP.PropsSI('P', 'T', T, 'Q', 1, Ref)
hL = CP.PropsSI('Hmolar', 'T', T, 'Q', 0, Ref)
hV = CP.PropsSI('Hmolar', 'T', T, 'Q', 1, Ref)

HHH1, PPP1, EEE1 = [], [], []
HHH2, PPP2, EEE2 = [], [], []

cNorm = colors.LogNorm(vmin=1e-12, vmax=10)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet'))

for a_useless_counter in range(40000):

    h = random.uniform(150000 * MM, 590000 * MM)
    p = 10**random.uniform(np.log10(100000), np.log10(7000000))
    CP.set_debug_level(0)
    try:

        EOS.update(CoolProp.HmolarP_INPUTS, h, p)
        rhoEOS = EOS.rhomolar(); TEOS = EOS.T()

        TTSE.update(CoolProp.HmolarP_INPUTS, h, p)
        rhoTTSE = TTSE.rhomolar(); TTTSE = TTSE.T()

        BICUBIC.update(CoolProp.HmolarP_INPUTS, h, p)
        rhoBICUBIC = BICUBIC.rhomolar(); TBICUBIC = BICUBIC.T()

        errorTTSE = abs(rhoTTSE / rhoEOS - 1) * 100
        errorBICUBIC = abs(rhoBICUBIC / rhoEOS - 1) * 100
        if errorTTSE > 100 or errorTTSE < 1e-12:
            print("%s %s %s" % (h, p, errorTTSE))

        HHH1.append(h)
        PPP1.append(p)
        EEE1.append(errorTTSE)

        HHH2.append(h)
        PPP2.append(p)
        EEE2.append(errorBICUBIC)

    except ValueError as VE:
        print('ERROR %s' % VE)
        pass

print('done')
SC1 = ax1.scatter(HHH1, PPP1, s=8, c=EEE1, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)
SC2 = ax2.scatter(HHH2, PPP2, s=8, c=EEE2, edgecolors='none', cmap=plt.get_cmap('jet'), norm=cNorm)

ax1.set_title('Error in Density from TTSE')
ax2.set_title('Error in Density from Bicubic')

for ax in [ax1, ax2]:

    ax.set_xlim(250000 * MM, 550000 * MM)
    ax.set_ylim(100000, 7000000)

    ax.set_yscale('log')

    ticks = [100000, 200000, 400000, 600000, 800000, 1000000, 2000000, 4000000, 6000000]
    labels = [str(tick) for tick in ticks]
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ticks = [150000 * MM, 250000 * MM, 350000 * MM, 450000 * MM, 550000 * MM]
    labels = [str(tick) for tick in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    ax.tick_params(axis='y', which='minor', left='off')

    ax.set_xticklabels(ax.get_xticks() / MM / 1e3)
    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_yticklabels(ax.get_yticks() / 10**3)
    ax.set_ylabel('Pressure [kPa]')

    ax.plot(hL, pV, 'k', lw=4)
    ax.plot(hV, pV, 'k', lw=4)

cbar_ax = fig.add_axes([0.85, 0.15, 0.06, 0.7])
CB = fig.colorbar(SC1, cax=cbar_ax)
CB.set_label(r'$(\rho/\rho_{EOS}-1)\times 100$ [%]')

plt.savefig('TTSE_BICUBIC.png', dpi=300, transparent=True)
plt.savefig('TTSE_BICUBIC.pdf')
plt.close()
