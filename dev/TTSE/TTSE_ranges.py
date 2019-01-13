from CoolProp.Plots import Ph
import CoolProp.CoolProp as CP
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker
import numpy as np
import random
import scipy.ndimage
import scipy.interpolate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib._pylab_helpers


def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = scipy.interpolate.interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B


def get_Z(X_in, Y_in, fluid):
    '''
    Just a wrapper to call CoolProp
    '''
    return False


def fill_Z(X, Y):
    '''
    Use 2D arrays X and Y to calculate a Z value
    '''
    Z = np.empty_like(X)
    for i in range(len(X)):
        for j in range(len(X[i])):
            Z[i, j] = get_Z(X[i, j], Y[i, j], fluid)
        Z[i] = fill_nan(Z[i])
    tr = np.transpose(Z)
    for i in range(len(tr)):
        tr[i] = fill_nan(tr[i])
    return np.transpose(tr)


def plotLineAndProjection(ax, X, Y, Z, draw='CXYZ', color='black', xlim=None, ylim=None, zlim=None):
    alpha = 0.5

    if 'C' in draw:
        ax.plot(X, Y, zs=Z, color=color)
    # else:
    #    ax.plot(X,Y,Z,color=color,alpha=0)

    if xlim == None: xlim = ax.get_xlim()
    if ylim == None: ylim = ax.get_ylim()
    if zlim == None: zlim = ax.get_zlim()

    if 'X' in draw:
        constArr = np.ones_like(X) * xlim[0]
        ax.plot(constArr, Y, zs=Z, color=color, alpha=alpha)
        ax.set_xlim(xlim)
    if 'Y' in draw:
        constArr = np.ones_like(Y) * ylim[1]
        ax.plot(X, constArr, zs=Z, color=color, alpha=alpha)
        ax.set_ylim(ylim)
    if 'Z' in draw:
        constArr = np.ones_like(Z) * zlim[0]
        ax.plot(X, Y, zs=constArr, color=color, alpha=alpha)
        ax.set_zlim(zlim)


def make3Dlpot(X, Y, Z=None, ax=None, invert='', draw='CXYZ', color='blue', xlim=None, ylim=None, zlim=None):
    '''
    A simple wrapper around the plotting routines

    invert: which axes to invert could be 'X' or 'YZ' or ...
    draw: which lines to draw could be 'C' or 'YZ' or ...
          C for curve and the rest are projections in the different dimensions
    '''

    if Z == None: Z = fill_Z(X, Y)

    # Now we have all data and need to resample it for a smoother plot
    resFactor = 10
    Xr = scipy.ndimage.zoom(X, resFactor, order=1)
    Yr = scipy.ndimage.zoom(Y, resFactor, order=1)
    Zr = scipy.ndimage.zoom(Z, resFactor, order=1)

    # Make the plot
    if ax == None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    else:
        fig = ''

    if 'X' in invert: ax.invert_xaxis()
    if 'Y' in invert: ax.invert_yaxis()
    if 'Z' in invert: ax.invert_zaxis()

    # Reduce data again and call the plotting wrapper
    stride = np.round(len(Xr) / 3.0)

    for i in range(len(Xr)):
        if np.mod(i, stride) == 0:
            plotLineAndProjection(ax, Xr[i], Yr[i], Zr[i], draw=draw, color=color, xlim=xlim, ylim=ylim, zlim=zlim)
    plotLineAndProjection(ax, Xr[-1], Yr[-1], Zr[-1], draw=draw, color=color, xlim=xlim, ylim=ylim, zlim=zlim)

    for i in range(len(Xr[0])):
        if np.mod(i, stride) == 0:
            Xi = [row[i] for row in Xr]
            Yi = [row[i] for row in Yr]
            Zi = [row[i] for row in Zr]
            plotLineAndProjection(ax, Xi, Yi, Zi, draw=draw, color=color, xlim=xlim, ylim=ylim, zlim=zlim)
    Xi = [row[-1] for row in Xr]
    Yi = [row[-1] for row in Yr]
    Zi = [row[-1] for row in Zr]
    plotLineAndProjection(ax, Xi, Yi, Zi, draw=draw, color=color, xlim=xlim, ylim=ylim, zlim=zlim)

    #cset = ax.contour(X, Y, Z, cmap=cm.coolwarm)
    #ax.clabel(cset, fontsize=9, inline=1)

    # In case we need a surface plot
    #from matplotlib import cm
    #from matplotlib.ticker import LinearLocator, FormatStrFormatter
    #
    # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    #        linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    #
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    #
    #fig.colorbar(surf, shrink=0.5, aspect=5)

    # Alternative plotting solution
    #cmap = plt.get_cmap('jet')
    #
    # if PRINT:
    #    ax.plot_surface(
    #      X, Y, Z,
    #      rstride=100,
    #      cstride=100,
    #      alpha=0.5,
    #      cmap=cmap,
    #      linewidth=0,
    #      antialiased=True
    #      )
    # else:
    #    #ax.plot_surface(
    #    ax.plot_wireframe(
    #      X, Y, Z,
    #      rstride=100,
    #      cstride=100,
    #      #alpha=0.5,
    #      #cmap=cmap,
    #      #linewidth=0,
    #      #antialiased=False
    #      )

    # Make the plot
    #fig  = plt.figure()
    #ax   = fig.gca(projection='3d')
    ##cset = ax.contour(X, Y, Z, cmap=cm.coolwarm)
    ##ax.clabel(cset, fontsize=9, inline=1)
    #
    # ax.plot_wireframe(
    #  X, Y, Z,
    #  rstride=np.round(resFactor*10),
    #  cstride=np.round(resFactor*10),
    #  #alpha=0.5,
    #  #cmap=cmap,
    #  #linewidth=0,
    #  #antialiased=False
    #  )

    # Plot the two-phase dome and its projections
    # ax.plot(X_TP,Y_TP,Z_TP,color='black')

    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()
    #zlim = ax.get_zlim()
    #
    #constArr = np.ones_like(X_TP) * xlim[0]
    # ax.plot(constArr,Y_TP,Z_TP,color='black')
    #constArr = np.ones_like(Y_TP) * ylim[1]
    # ax.plot(X_TP,constArr,Z_TP,color='black')
    #constArr = np.ones_like(Z_TP) * zlim[0]
    # ax.plot(X_TP,Y_TP,constArr,color='black')

    # Or do we want contour lines?
    ##cset = ax.contour(X, Y, Z, zdir='z', offset=zlim[0], cmap=cmap)
    ##cset = ax.contour(X, Y, Z, zdir='x', offset=xlim[0], cmap=cmap)
    ##cset = ax.contour(X, Y, Z, zdir='y', offset=ylim[1], cmap=cmap)
    #
    # ax.set_xlim(xlim)
    # ax.set_ylim(ylim)
    # ax.set_zlim(zlim)

    #majorFormatter = matplotlib.ticker.EngFormatter("J/kg")
    # majorFormatter = matplotlib.ticker.
    majorFormatter = matplotlib.ticker.ScalarFormatter()
    # majorFormatter.set_scientific(True)
    majorFormatter.set_powerlimits((-4, 4))
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.zaxis.set_major_formatter(majorFormatter)
    return fig, ax, Z


def getlim(key, dicts, fac=1):
    min = np.min([dict[key] / fac for dict in dicts])
    max = np.max([dict[key] / fac for dict in dicts])
    return [np.floor(min) * fac, np.ceil(max) * fac]


#########################################################################
# Here starts the real script
#########################################################################


fluid = 'water'
CP.enable_TTSE_LUT(fluid)

PRINT = False

if PRINT:
    points = 200
    dpi = 300
    toFile = True
else:
    points = 100
    dpi = 75
    toFile = False


Tmin = CP.PropsSI('Tmin', 'T', 0, 'P', 0, fluid) + 5.0
Tcri = CP.PropsSI('Tcrit', 'T', 0, 'P', 0, fluid)
Tmax = CP.PropsSI('Tcrit', 'T', 0, 'P', 0, fluid) * 2.0

# Get phase boundary
T_TP = np.linspace(Tmin, Tcri - 0.5, 0.5 * points)
D_TP = CP.PropsSI('D', 'T', T_TP, 'Q', 0, fluid)
P_TP = CP.PropsSI('P', 'T', T_TP, 'Q', 0, fluid)
H_TP = CP.PropsSI('H', 'T', T_TP, 'Q', 0, fluid)
S_TP = CP.PropsSI('S', 'T', T_TP, 'Q', 0, fluid)

D_TP = np.append(D_TP, CP.PropsSI('D', 'T', T_TP[::-1], 'Q', 1, fluid))
P_TP = np.append(P_TP, CP.PropsSI('P', 'T', T_TP[::-1], 'Q', 1, fluid))
H_TP = np.append(H_TP, CP.PropsSI('H', 'T', T_TP[::-1], 'Q', 1, fluid))
S_TP = np.append(S_TP, CP.PropsSI('S', 'T', T_TP[::-1], 'Q', 1, fluid))

logV_TP = np.log10(1.0 / D_TP)
logP_TP = np.log10(1.0 * P_TP)
T_TP = np.append(T_TP, T_TP[::-1])


# These are the default TTSE ranges
D_L = CP.PropsSI('D', 'T', Tmin, 'Q', 0, fluid)
D_V = CP.PropsSI('D', 'T', Tmin, 'Q', 1, fluid)
H_L = CP.PropsSI('H', 'T', Tmin, 'Q', 0, fluid)
H_V = CP.PropsSI('H', 'T', Tmin, 'Q', 1, fluid)
P_L = CP.PropsSI('P', 'T', Tmin, 'Q', 0, fluid)
P_V = CP.PropsSI('P', 'T', Tmin, 'Q', 1, fluid)

Hmin = H_L
Hmax = H_L + (H_V - H_L) * 2.0

#Pmin = CP.PropsSI('ptriple','T',0,'P',0,fluid) + 1
Pmin = np.min([P_L, P_V])
Pmax = CP.PropsSI('pcrit', 'T', 0, 'P', 0, fluid) * 2.0  # should be p_reduce

Dmin = CP.PropsSI('D', 'H', Hmin, 'P', Pmax, fluid)
Dmax = CP.PropsSI('D', 'H', Hmax, 'P', Pmin, fluid)


#########################################################################
# Start with the first diagram, hps
#########################################################################
# Set the ranges for the plot
X_TP = H_TP
Y_TP = logP_TP
Z_TP = S_TP

Xmin = Hmin
Xmax = Hmax
Ymin = np.log10(Pmin)
Ymax = np.log10(Pmax)

X = np.linspace(Xmin, Xmax, points)  # Enthalpy
Y = np.linspace(Ymin, Ymax, points)  # Pressure
X, Y = np.meshgrid(X, Y)


def get_Z(X_in, Y_in, fluid, out='S'):
    '''
    Just a wrapper to call CoolProp
    '''
    X = X_in
    Y = np.power(10.0, Y_in)
    Z = np.NAN
    try:
        Z = CP.PropsSI(out, 'H', X, 'P', Y, fluid)
    except(ValueError):
        print("CoolProp failed, returning NAN")
        Z = np.NAN
    return Z


# Now we also need the variables
Z = fill_Z(X, Y)
get_Z_old = get_Z


def get_Z(X_in, Y_in, fluid):
    return get_Z_old(X_in, Y_in, fluid, out='D')


logVdata = np.log10(1.0 / fill_Z(X, Y))


def get_Z(X_in, Y_in, fluid):
    return get_Z_old(X_in, Y_in, fluid, out='T')


Tdata = fill_Z(X, Y)

HPSdict = {'H': X, 'P': Y, 'S': Z, 'V': logVdata, 'T': Tdata, 'H_TP': X_TP, 'P_TP': Y_TP, 'S_TP': Z_TP}


#########################################################################
# Start with the next diagram, vTp
#########################################################################

# Set the ranges for the plot
X_TP = logV_TP
Y_TP = T_TP
Z_TP = logP_TP

#Xmin = np.log10(1.0/Dmax)
#Xmax = np.log10(1.0/Dmin)
Xmin = np.min(X_TP)
Xmax = np.max(X_TP)
Ymin = np.min(Y_TP)
Ymax = Tmax

X = np.linspace(Xmin, Xmax, points)  # Volume
Y = np.linspace(Ymin, Ymax, points)  # Temperature
X, Y = np.meshgrid(X, Y)


def get_Z(X_in, Y_in, fluid, out='P'):
    '''
    Just a wrapper to call CoolProp
    '''
    X = 1.0 / np.power(10.0, X_in)
    Y = Y_in
    Z = np.NAN
    try:
        Z = np.log10(CP.PropsSI(out, 'D', X, 'T', Y, fluid))
    except(ValueError):
        print("CoolProp failed, returning NAN")
        Z = np.NAN
    return Z


# Now we also need the variables
Z = fill_Z(X, Y)
get_Z_old = get_Z


def get_Z(X_in, Y_in, fluid):
    return np.power(10.0, get_Z_old(X_in, Y_in, fluid, out='H'))


Hdata = fill_Z(X, Y)


def get_Z(X_in, Y_in, fluid):
    return np.power(10.0, get_Z_old(X_in, Y_in, fluid, out='S'))


Sdata = fill_Z(X, Y)

VTPdict = {'V': X, 'T': Y, 'P': Z, 'H': Hdata, 'S': Sdata, 'V_TP': X_TP, 'T_TP': Y_TP, 'P_TP': Z_TP}


#########################################################################
# Now we have all the data and can start mixing the
# different definitions
#########################################################################

dicts = [HPSdict, VTPdict]

Hlim = getlim('H', dicts, fac=1e6)
Plim = getlim('P', dicts)
Slim = getlim('S', dicts, fac=1e3)
Vlim = getlim('V', dicts)
Tlim = getlim('T', dicts, fac=1e2)

figHPS, axHPS, Z = make3Dlpot(HPSdict['H'], HPSdict['P'], Z=HPSdict['S'], invert='Z', draw='XYZ', color='blue', xlim=Hlim, ylim=Plim, zlim=Slim[::-1])
# Plot the two-phase dome and its projections
plotLineAndProjection(axHPS, HPSdict['H_TP'], HPSdict['P_TP'], HPSdict['S_TP'])
# .. and add the other plot
make3Dlpot(VTPdict['H'], VTPdict['P'], Z=VTPdict['S'], ax=axHPS, draw='XYZ', color='red')

axHPS.set_xlabel(r'$h$')
axHPS.set_ylabel(r'log $p$')
axHPS.set_zlabel(r'$s$')


figVTP, axVTP, Z = make3Dlpot(VTPdict['V'], VTPdict['T'], Z=VTPdict['P'], draw='XYZ', color='red', xlim=Vlim, ylim=Tlim, zlim=Plim)
# Plot the two-phase dome and its projections
plotLineAndProjection(axVTP, VTPdict['V_TP'], VTPdict['T_TP'], VTPdict['P_TP'])
# .. and add the other plot
make3Dlpot(HPSdict['V'], HPSdict['T'], Z=HPSdict['P'], ax=axVTP, draw='XYZ', color='blue')

axVTP.set_xlabel(r'log $v$')
axVTP.set_ylabel(r'$T$')
axVTP.set_zlabel(r'log $p$')


if toFile:
    figHPS.savefig('phs-' + fluid + '.png', dpi=dpi, transparent=True)
    figHPS.savefig('phs-' + fluid + '.pdf', transparent=True)
    # figHPS.close()
    figVTP.savefig('pvT-' + fluid + '.png', dpi=dpi, transparent=True)
    figVTP.savefig('pvT-' + fluid + '.pdf', transparent=True)
    # figVTP.close()
    # plt.savefig('TTSE_RANGES.pdf')
else:
    plt.show()
    # figHPS.show()
    # figVTP.show()
