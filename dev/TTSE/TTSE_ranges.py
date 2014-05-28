from CoolProp.Plots import Ph
import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker
import numpy as np
import random

#fig = plt.figure(figsize=(10,5))
#ax1 = fig.add_axes((0.08,0.1,0.32,0.83))
#ax2 = fig.add_axes((0.50,0.1,0.32,0.83))
#
#cNorm  = colors.LogNorm(vmin=1e-12, vmax=10)
#scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = plt.get_cmap('jet'))
#
#
#
#plt.savefig('TTSE_BICUBIC.png', dpi = 300, transparent = True)
#plt.savefig('TTSE_BICUBIC.pdf')
#plt.close()

fluid = 'water'

PRINT = True

if PRINT:
    points  = 500
    rstride = 1
    cstride = 1
    antiAli = True
    toFile  = True
else:
    points  = 200
    rstride = 5
    cstride = 5
    antiAli = False 
    toFile  = False
    

Tmin = CP.PropsSI('Tmin','T',0,'P',0,fluid)
Tcri = CP.PropsSI('Tcrit','T',0,'P',0,fluid)
Tmax = CP.PropsSI('Tcrit','T',0,'P',0,fluid) * 2.0 

## Get phase boundary
Y_TP = np.linspace(Tmin, Tcri-0.5, 0.5*points)

X_TP = CP.PropsSI('D','T',Y_TP,'Q',0,fluid)
Z_TP = CP.PropsSI('P','T',Y_TP,'Q',0,fluid)

X_TP = np.append(X_TP, CP.PropsSI('D','T',Y_TP[::-1],'Q',1,fluid))
Z_TP = np.append(Z_TP, CP.PropsSI('P','T',Y_TP[::-1],'Q',1,fluid))

X_TP = np.log10(1.0/X_TP)
Z_TP = np.log10(1.0*Z_TP)

Y_TP = np.append(Y_TP, Y_TP[::-1])


pmin = CP.PropsSI('ptriple','T',0,'P',0,fluid) + 1

vmin = np.log10(1.0/CP.PropsSI('D','T',Tmin,'P',pmin,fluid))
vmax = np.max(X_TP)



X = np.linspace(vmin,vmax,points) # Volume
Y = np.linspace(Tmin, Tmax,points) # Temperature
X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)





def fillZ(X,Y):
    Z = np.empty_like(X)
    for i in range(len(X)):
        Z[i] = np.log10(CP.PropsSI('P','T',Y[i],'D',1.0/np.power(10.0,X[i]),fluid))
    return Z



Z    = fillZ(X, Y)


from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')


#mpl.rcParams['legend.fontsize'] = 10
#
#fig = plt.figure()
##ax = fig.gca(projection='3d')
#ax = Axes3D(fig)
#theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
#z = np.linspace(-2, 2, 100)
#r = z**2 + 1
#x = r * np.sin(theta)
#y = r * np.cos(theta)
#ax.plot(x, y, z, label='parametric curve')
#ax.legend()


from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)
#
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
#fig.colorbar(surf, shrink=0.5, aspect=5)




fig  = plt.figure()
ax   = fig.gca(projection='3d')
#cset = ax.contour(X, Y, Z, cmap=cm.coolwarm)
#ax.clabel(cset, fontsize=9, inline=1)

ax.plot(X_TP,Y_TP,Z_TP,color='black')

cmap = plt.get_cmap('jet')

ax.plot_surface(
  X, Y, Z, 
  rstride=rstride, 
  cstride=cstride, 
  alpha=0.5, 
  cmap=cmap,
  linewidth=0, 
  antialiased=antiAli
  )

#xlim = ax.get_xlim()
#ylim = ax.get_ylim()
#zlim = ax.get_zlim()
#
##cset = ax.contour(X, Y, Z, zdir='z', offset=zlim[0], cmap=cmap)
##cset = ax.contour(X, Y, Z, zdir='x', offset=xlim[0], cmap=cmap)
##cset = ax.contour(X, Y, Z, zdir='y', offset=ylim[1], cmap=cmap)
#
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)
#ax.set_zlim(zlim)



ax.set_xlabel(r'log $v$')
ax.set_ylabel(r'$T$')
ax.set_zlabel(r'log $p$')

#fig  = plt.figure()
#ax   = fig.gca(projection='3d')
#surf = ax.plot_surface(
#         X, Y, Z, 
#         #rstride=1, 
#         #cstride=1, 
#         #cmap=plt.get_cmap('jet'),
#         linewidth=0, 
#         antialiased=False
#         )
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)
if toFile:
    plt.savefig('TTSE_RANGES.png', dpi = 300, transparent = True)
    #plt.savefig('TTSE_RANGES.pdf')
    plt.close()
else:
    plt.show()
