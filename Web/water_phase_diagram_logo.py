import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt

Water = CP.AbstractState("HEOS", "Water")
pc = Water.keyed_output(CP.iP_critical)
Tc = Water.keyed_output(CP.iT_critical)
Tmin = 200
Tmax = 1000
pmax = Water.keyed_output(CP.iP_max)
pt = 611.657
Tt = 273.16
fillcolor = 'g'

fig = plt.figure(figsize = (3,3))
ax = fig.add_axes((0.02,0.02,0.98,0.98))
lw = 4
melt_args = dict(lw = lw)
TT = []
PP = list(np.logspace(np.log10(pt), np.log10(pmax),1000))
for p in PP:
    TT.append(Water.melting_line(CP.iT, CP.iP, p))

#Zone VI
for T in np.linspace(max(TT), 355):
    TT.append(T)
    theta = T/273.31
    pi = 1-1.07476*(1-theta**4.6)
    p = pi*632.4e6
    PP.append(p)

#Zone VII
for T in np.linspace(355, 715):
    TT.append(T)
    theta = T/355
    lnpi = 0.173683e1*(1-1/theta)-0.544606e-1*(1-theta**5)+0.806106e-7*(1-theta**22)
    p = np.exp(lnpi)*2216e6
    PP.append(p)
plt.plot(TT,PP,'g',**melt_args)

# Supercritical liquid
#plt.fill_betweenx(PP,TT,Tc,color=fillcolor,alpha = 0.9, where =  np.logical_and(np.array(PP)>pc, np.array(TT)<Tc) )

Ts = np.linspace(273.16, Tc, 1000)
ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',[0]*len(Ts),'Water',[1])

#~ TT = []
#~ PP = np.logspace(np.log10(pt), np.log10(pc), 1000)
#~ for p in PP:
    #~ TT.append(Water.melting_line(CP.iT, CP.iP, p))
#~ plt.plot(TT,PP,'g',**melt_args)

# Liquid
#plt.fill(np.r_[TT,Ts[::-1]],np.r_[PP,ps[::-1]],fillcolor, alpha = 0.7)
# Gas
#plt.fill(np.r_[Ts,Tc,Ts[0]],np.r_[ps,pt,ps[0]],fillcolor, alpha = 0.1)
# Supercritical gas
#plt.fill(np.r_[Tc,Tmax, Tmax, Tc, Tc],np.r_[pc, pc, pt, pt, pc],fillcolor, alpha = 0.3)
# Supercritical
pe = 1.3e10
#plt.fill(np.r_[Tc, Tc, Tmax, Tmax],np.r_[pc, pe, pe, pc],fillcolor, alpha = 0.5)

plt.plot(Ts,ps,'k',lw = lw)
plt.plot(Ts,ps,'w',lw = lw/2)

plt.plot(Tc,pc,'ro')

TD,DD,PD = [], [], []
for T in np.linspace(310,1000,200):
    for p in np.logspace(np.log10(1000),np.log10(1e9),200):
        D = CP.CoolProp.PropsSI('D','T',T,'P',p,'Water')
        TD.append(T)
        DD.append(np.log10(D))
        PD.append(p)
        
nm = matplotlib.colors.Normalize(min(DD), max(DD))
plt.scatter(TD, PD, c = DD, edgecolor = 'none', s = 6, cmap=plt.cm.get_cmap('Blues'), norm = nm)

plt.ylim(611,max(PP)*1.2)
plt.gca().set_yscale('log')    
#plt.gca().set_xlim(Tmin*0.99, Tmax)
plt.gca().set_xlim(240, 1000)
plt.gca().axis('off')
plt.savefig('WaterPhaseDiagram.pdf')
plt.savefig('WaterPhaseDiagram.png', transparent = True)
plt.close()