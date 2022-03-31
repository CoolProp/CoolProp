import CoolProp
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import random

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_axes((0.08,0.1,0.32,0.83))
ax2 = fig.add_axes((0.50,0.1,0.32,0.83))

Ref = 'helium'

BICUBIC = CoolProp.AbstractState('BICUBIC&HEOS',Ref)
TTSE = CoolProp.AbstractState('TTSE&HEOS',Ref)
EOS = CoolProp.AbstractState('HEOS',Ref)

T = np.linspace(CP.PropsSI(Ref,'Tmin')+0.1, CP.PropsSI(Ref,'Tcrit')-0.01, 300)
DL = CP.PropsSI('Dmolar','T',T,'Q',0,Ref)
DV = CP.PropsSI('Dmolar','T',T,'Q',1,Ref)
uL = CP.PropsSI('Umolar','T',T,'Q',0,Ref)
uV = CP.PropsSI('Umolar','T',T,'Q',1,Ref)

UUU, DDD, EEE1, EEE2  = [], [], [], []

cNorm  = colors.LogNorm(vmin=1e-12, vmax=10)
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = plt.get_cmap('jet'))

for a_useless_counter in range(40000):

    umolar = random.uniform(-50,100)#random.uniform(-50,50000)
    rhomolar = 10**random.uniform(2,5.5) #10**random.uniform(-10,5.5)
    CP.set_debug_level(0)
    try:

        EOS.update(CoolProp.DmolarUmolar_INPUTS, rhomolar, umolar)
        pEOS = EOS.p(); TEOS = EOS.T()

        TTSE.update(CoolProp.DmolarUmolar_INPUTS, rhomolar, umolar)
        pTTSE = TTSE.p(); TTTSE = TTSE.T()

        BICUBIC.update(CoolProp.DmolarUmolar_INPUTS, rhomolar, umolar)
        pBICUBIC = BICUBIC.p(); TBICUBIC = BICUBIC.T()

        errorTTSE = abs(pTTSE/pEOS-1)*100
        errorBICUBIC = abs(pBICUBIC/pEOS-1)*100
        if errorTTSE > 100 or errorTTSE < 1e-12:
            print(umolar, rhomolar, errorTTSE)

        UUU.append(umolar)
        DDD.append(rhomolar)
        EEE1.append(errorTTSE)

        EEE2.append(errorBICUBIC)

    except ValueError as VE:
        print('ERROR', VE)
        
    except:
        print("Unkwnown error")
        raise
    

SC1 = ax1.scatter(DDD, UUU, s = 8, c = EEE1, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)
SC2 = ax2.scatter(DDD, UUU, s = 8, c = EEE2, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)

ax1.set_title('Error in Pressure from TTSE')
ax2.set_title('Error in Pressure from Bicubic')

for ax in [ax1, ax2]:

    ax.set_xscale('log')
    ax.set_ylabel('Specific internal energy [J/mol]')
    ax.set_xlabel('Molar density [m^3/mol^3]')

    ax.plot(DL,uL,'k',lw = 4)
    ax.plot(DV,uV,'k',lw = 4)

cbar_ax = fig.add_axes([0.85, 0.15, 0.06, 0.7])
CB = fig.colorbar(SC1, cax=cbar_ax)
CB.set_label(r'$(p/p_{EOS}-1)\times 100$ [%]')

plt.show()