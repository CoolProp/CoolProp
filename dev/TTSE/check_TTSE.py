from CoolProp.Plots import Ph
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker
import numpy as np
import random

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_axes((0.08,0.1,0.32,0.83))
ax2 = fig.add_axes((0.50,0.1,0.32,0.83))

Ref = 'R245fa'

T = np.linspace(CP.Props(Ref,'Tmin')+0.1,CP.Props(Ref,'Tcrit')-0.01,300)
pV = CP.Props('P','T',T,'Q',1,Ref)
hL = CP.Props('H','T',T,'Q',0,Ref)
hV = CP.Props('H','T',T,'Q',1,Ref)

HHH1, PPP1, EEE1 = [], [], []
HHH2, PPP2, EEE2 = [], [], []

cNorm  = colors.LogNorm(vmin=1e-12, vmax=10)
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = plt.get_cmap('jet'))

for a_useless_counter in range(40000):
        
    h = random.uniform(100,590)
    p = 10**random.uniform(np.log10(100),np.log10(7000))
    
    try:
        # Using the EOS
        CP.disable_TTSE_LUT(Ref)
        rhoEOS = CP.Props('D','P',p,'H',h,Ref)
        TEOS = CP.Props('T','P',p,'H',h,Ref)
##         cpEOS = CP.Props('C','P',p,'H',h,Ref)
        
        # Using the TTSE method
        CP.enable_TTSE_LUT(Ref)
        CP.set_TTSE_mode(Ref,"TTSE")
        rhoTTSE = CP.Props('D','P',p,'H',h,Ref)
        TTTSE = CP.Props('T','P',p,'H',h,Ref)
##         cpTTSE = CP.Props('C','P',p,'H',h,Ref)
        
        # Using the Bicubic method
        CP.enable_TTSE_LUT(Ref)
        CP.set_TTSE_mode(Ref,"BICUBIC")
        rhoBICUBIC = CP.Props('D','P',p,'H',h,Ref)
        TBICUBIC = CP.Props('T','P',p,'H',h,Ref)
##         cpBICUBIC = CP.Props('C','P',p,'H',h,Ref)
        
##         errorTTSE = abs(TTTSE/TEOS-1)*100
##         errorBICUBIC = abs(TBICUBIC/TEOS-1)*100
        errorTTSE = abs(rhoTTSE/rhoEOS-1)*100
        errorBICUBIC = abs(rhoBICUBIC/rhoEOS-1)*100
##         errorTTSE = abs(cpTTSE/cpEOS-1)*100
##         errorBICUBIC = abs(cpBICUBIC/cpEOS-1)*100
        HHH1.append(h)
        PPP1.append(p)
        EEE1.append(errorTTSE)
        
        HHH2.append(h)
        PPP2.append(p)
        EEE2.append(errorBICUBIC)
        
    except ValueError as VE:
        print VE
        pass
    
SC1 = ax1.scatter(HHH1, PPP1, s = 8, c = EEE1, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)
SC2 = ax2.scatter(HHH2, PPP2, s = 8, c = EEE2, edgecolors = 'none', cmap = plt.get_cmap('jet'), norm = cNorm)

ax1.set_title('Error in Density from TTSE')
ax2.set_title('Error in Density from Bicubic')

for ax in [ax1, ax2]:
    
    ax.set_xlim(150, 550)
    ax.set_ylim(100, 7000)

    ax.set_yscale('log')
    
    ticks = [100,200,400,600,800,1000,2000, 4000, 6000]
    labels = [str(tick) for tick in ticks]
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    
    ticks = [150,250,350,450,550]
    labels = [str(tick) for tick in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    ax.tick_params(axis='y',which='minor', left='off')

    ax.set_xlabel('Enthalpy [kJ/kg]')
    ax.set_ylabel('Pressure [kPa]')

    ax.plot(hL,pV,'k',lw = 4)
    ax.plot(hV,pV,'k',lw = 4)

cbar_ax = fig.add_axes([0.85, 0.15, 0.06, 0.7])
CB = fig.colorbar(SC1, cax=cbar_ax)
CB.set_label(r'$(\rho/\rho_{EOS}-1)\times 100$ [%]')
#CB.set_label(r'$(T/T_{EOS}-1)\times 100$ [%]')

plt.savefig('TTSE_BICUBIC.png', dpi = 300, transparent = True)
plt.savefig('TTSE_BICUBIC.pdf')
plt.close()