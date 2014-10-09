import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import scipy.interpolate

# Prepare the constants 
Water = CP.AbstractState("HEOS", "Water")
pc = Water.keyed_output(CP.iP_critical)
Tc = Water.keyed_output(CP.iT_critical)
T_min = 200
T_max = 1000
p_max = Water.keyed_output(CP.iP_max)
p_triple = 611.657
T_triple = 273.16

# Prepare the data for the melting line
steps = 2000

TT = []
PP = list(np.logspace(np.log10(p_triple), np.log10(p_max),steps))
for p in PP:
    TT.append(Water.melting_line(CP.iT, CP.iP, p))

#Zone VI
for T in np.linspace(max(TT), 355, int(steps/10)):
    TT.append(T)
    theta = T/273.31
    pi = 1-1.07476*(1-theta**4.6)
    p = pi*632.4e6
    PP.append(p)

#Zone VII
for T in np.linspace(355, 715, int(steps/10)):
    TT.append(T)
    theta = T/355
    lnpi = 0.173683e1*(1-1/theta)-0.544606e-1*(1-theta**5)+0.806106e-7*(1-theta**22)
    p = np.exp(lnpi)*2216e6
    PP.append(p)

# Changes number of points
steps = int(steps/10.0)

p_melt   = np.logspace(np.log10(np.min(PP)), np.log10(np.max(PP)),steps)
T_melt_f = scipy.interpolate.interp1d(np.log10(PP),TT)
#T_melt_f = scipy.interpolate.spline(np.log10(PP),TT,np.log10(p_melt))
T_melt = T_melt_f(np.log10(p_melt))
#T_melt = np.array(TT)
#p_melt = np.array(PP)

#
# Prepare the data for the saturation line
T_sat = np.linspace(273.16, Tc, len(T_melt))
p_sat = CP.CoolProp.PropsSI('P','T',T_sat,'Q',[0]*len(T_sat),'Water',[1])

#
# Prepare density data
TT,DD,PP = [], [], []
for T in np.linspace(T_min, T_max, steps):
    for p in np.logspace(np.log10(p_triple), np.log10(p_max), steps):
        Tm = scipy.interpolate.interp1d(p_melt, T_melt)(p)
        if T < Tm: continue
        D = CP.CoolProp.PropsSI('D','T',T,'P',p,'Water')
        TT.append(T)
        DD.append(np.log10(D))
        PP.append(p)

tt = np.linspace(T_min, T_max, steps)
pp = np.logspace(np.log10(p_triple), np.log10(p_max), steps)
tt, pp = np.meshgrid(tt, pp)
dd = np.empty(tt.shape)
dd[:][:] = np.NAN

nr,nc = tt.shape

for i in range(nr):
    for j in range(nc):
        Tm = T_melt_f(np.log10(pp[i][j]))
        if tt[i][j] < Tm: continue
        D = CP.CoolProp.PropsSI('D','T',tt[i][j],'P',pp[i][j],'Water')
        dd[i][j] = np.log10(D)

#
# Define colours etc
lw = 3
melt_args = dict(color = 'orange', lw = lw, solid_capstyle = 'round')
sat_args  = melt_args.copy()

nm = matplotlib.colors.Normalize(min(DD), max(DD))
rho_args = dict(cmap=plt.cm.get_cmap('Blues'), norm = nm)

fig = plt.figure(figsize = (3,3))
ax = fig.add_axes((0.0,0.0,1.0,1.0))

plt.plot(T_melt, p_melt, **melt_args)
plt.plot(T_sat,  p_sat,  **sat_args )
plt.scatter(TT, PP, c=DD, edgecolor = 'none', s = 6, **rho_args )
#plt.contourf(tt, pp, dd, steps, **rho_args )

#ax.set_ylim(p_triple,p_max)
ax.set_yscale('log') 
ax.set_xlim(T_min, T_max)
ax.axis('off')

plt.savefig('WaterPhaseDiagram.pdf')
plt.savefig('WaterPhaseDiagram.png', transparent = True)
plt.close()