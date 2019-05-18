import matplotlib
matplotlib.use('WXAgg')
import CoolProp
from CoolProp.Plots import Ph
from CoolProp.Plots.Plots import Trho, Ps, PT, Prho
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import random
import numpy as np
from math import log, exp
random.seed()


def check_Trho(N=5000, param='P', fluid='R245fa'):
    values = []
    CP.enable_TTSE_LUT(fluid)
    try:
        CP.Props('D', 'P', CP.Props(fluid, 'ptriple') + 1, 'Q', 1, fluid)
    except:
        return []
    # CP.set_TTSESinglePhase_LUT_size(fluid,500,500)
    hmin, hmax, pmin, pmax = CP.get_TTSESinglePhase_LUT_range(fluid)
    for i in range(N):
        x1 = random.random()
        h = x1 * hmin + (1 - x1) * hmax
        x2 = random.random()
        logp = x2 * log(pmin) + (1 - x2) * log(pmax)
        p = exp(logp)
        try:
            try:
                # Get the T,rho from the EOS directly without the LUT
                CP.disable_TTSE_LUT(fluid)
                s = CP.Props('S', 'P', p, 'H', h, fluid)
                T = CP.Props('T', 'P', p, 'H', h, fluid)
                rho = CP.Props('D', 'P', p, 'H', h, fluid)
            except:
                print('EOS failed: %s %s' % (p, h))
                raise
            # Now get p,h from the T,rho
            CP.enable_TTSE_LUT(fluid)
            val = CP.Props(param, 'T', T, 'D', rho, fluid)
            CP.disable_TTSE_LUT(fluid)
            valREFPROP = CP.Props(param, 'T', T, 'D', rho, fluid)
            # print T,rho,val,valREFPROP,(val/valREFPROP-1)*100
            if abs(val - valREFPROP) > 0.00001:
                raise ValueError
        except ValueError:
            print('TTSE failed: %s %s' % (T, rho))
            values.append((T, rho, 0, 0))
            pass
    return values


fluid = 'R245fa'
values = check_Trho(param='P', fluid=fluid)
if len(values) == 0:
    print('good')

T, rho, values_withTTSE, values_noTTSE = zip(*values)

CP.disable_TTSE_LUT(fluid)
Trho(fluid)

plt.plot(rho, T, '.', mfc='none')
plt.savefig('Trho_TTSE_Validation.png', dpi=300)
plt.gca().set_xscale('log')
plt.show()

quit()


def check_Pother(N=5000, param='T', other='S', fluid='R245fa'):
    values = []
    CP.enable_TTSE_LUT(fluid)
    try:
        CP.Props('D', 'P', CP.Props(fluid, 'ptriple') + 1, 'Q', 1, fluid)
    except:
        return []
    # CP.set_TTSESinglePhase_LUT_size(fluid,500,500)
    hmin, hmax, pmin, pmax = CP.get_TTSESinglePhase_LUT_range(fluid)
    for i in range(N):
        x1 = random.random()
        h = x1 * hmin + (1 - x1) * hmax
        x2 = random.random()
        logp = x2 * log(pmin) + (1 - x2) * log(pmax)
        p = exp(logp)
        try:
            try:
                # Get the T,rho from the EOS directly without the LUT
                CP.disable_TTSE_LUT(fluid)
                s = CP.Props('S', 'P', p, 'H', h, fluid)
                T = CP.Props('T', 'P', p, 'H', h, fluid)
                rho = CP.Props('D', 'P', p, 'H', h, fluid)
            except:
                print('EOS failed: %s %s' % (p, h))
                raise
            # Now get p,h from the T,rho
            CP.enable_TTSE_LUT(fluid)
            if other == 'S':
                other_val = s
            elif other == 'T':
                other_val = T
            elif other == 'D':
                other_val = rho
            else:
                raise ValueError
            val = CP.Props(param, 'P', p, other, other_val, fluid)
        except ValueError:
            print('TTSE failed: %s %s' % (p, other_val))
            values.append((p, other_val, 0, 0))
            pass
    return values


fluid = 'R245fa'
for other in ['D', 'T', 'S']:
    if other == 'D':
        param = 'T'
    elif other == 'T':
        param = 'D'
    elif other == 'S':
        param = 'T'
    values = check_Pother(other=other)
    if len(values) == 0:
        print('good: %s' % other)
        continue

    p, othervals, values_withTTSE, values_noTTSE = zip(*values)

    CP.disable_TTSE_LUT(fluid)
    if other == 'D':
        Prho(fluid)
    elif other == 'S':
        Ps(fluid)
    elif other == 'T':
        PT(fluid)

    plt.plot(othervals, p, '.', mfc='none')
    plt.gca().set_yscale('log')
    plt.savefig('P' + other + '_TTSE_Validation.png', dpi=300)
    plt.show()

quit()


def check(N=5000, param='D', fluid='R245fa'):
    values = []
    CP.enable_TTSE_LUT(fluid)

    try:
        CP.Props('D', 'P', CP.Props(fluid, 'ptriple') + 1, 'Q', 1, fluid)
    except:
        return []
    # CP.set_TTSESinglePhase_LUT_size(fluid,500,500)
    hmin, hmax, pmin, pmax = CP.get_TTSESinglePhase_LUT_range(fluid)
    for i in range(N):
        x1 = random.random()
        h = x1 * hmin + (1 - x1) * hmax
        x2 = random.random()
        logp = x2 * log(pmin) + (1 - x2) * log(pmax)
        p = exp(logp)
        try:
            CP.enable_TTSE_LUT(fluid)
            value_withTTSE = CP.Props(param, 'P', p, 'H', h, fluid)
            CP.disable_TTSE_LUT(fluid)
            value_noTTSE = CP.Props(param, 'P', p, 'H', h, fluid)
            values.append((h, p, value_withTTSE, value_noTTSE))
        except ValueError:
            pass

    return values


Ncols = 10
Nrows = 10
for parameter in ['D', 'T', 'S', 'C']:
    fig = plt.figure(figsize=(40, 40))
    for Index, fluid in enumerate(sorted(CoolProp.__fluids__)):
        print(fluid)
        ax = fig.add_subplot(Ncols, Nrows, Index + 1)
        ax.set_title(fluid)
        values = check(fluid=fluid, param=parameter)
        if len(values) == 0:
            continue
        h, p, values_withTTSE, values_noTTSE = zip(*values)
        ax = fig.add_subplot(Ncols, Nrows, Index + 1)
        CP.disable_TTSE_LUT(fluid)
        Ph(fluid)
        CP.enable_TTSE_LUT(fluid)
        error = (np.array(values_withTTSE) / np.array(values_noTTSE) - 1) * 100
        plt.scatter(h, p, s=8, c=np.abs(error), norm=LogNorm(), edgecolor='none', vmin=1e-16, vmax=10)
        plt.gca().set_yscale('log')
        plt.colorbar()

    plt.savefig(parameter + '_TTSE.png', dpi=200)
    plt.tight_layout()
    plt.close()
