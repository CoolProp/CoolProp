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

def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = scipy.interpolate.interp1d(inds[good], A[good],bounds_error=False)
    B = np.where(np.isfinite(A),A,f(inds))
    return B


def get_Z(X_in,Y_in,fluid):
    '''
    Just a wrapper to call CoolProp
    '''
    return False 


def fill_Z(X,Y):
    '''
    Use 2D arrays X and Y to calculate a Z value
    '''
    Z = np.empty_like(X)
    for i in range(len(X)):
        for j in range(len(X[i])):
            Z[i,j] = get_Z(X[i,j],Y[i,j],fluid)
        Z[i] = fill_nan(Z[i])
    tr = np.transpose(Z) 
    for i in range(len(tr)):
        tr[i] = fill_nan(tr[i])
    return np.transpose(tr)


def plotLineAndProjection(ax,X,Y,Z,draw='CXYZ',color='black'):
    alpha = 0.5
    
    if 'C' in draw:
        ax.plot(X,Y,zs=Z,color=color)
    #else:
    #    ax.plot(X,Y,Z,color=color,alpha=0)
        
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()
    
    if 'X' in draw:
        constArr = np.ones_like(X) * xlim[0]  
        ax.plot(constArr,Y,zs=Z,color=color,alpha=alpha)
    if 'Y' in draw:
        constArr = np.ones_like(Y) * ylim[1]
        ax.plot(X,constArr,zs=Z,color=color,alpha=alpha)
    if 'Z' in draw:
        constArr = np.ones_like(Z) * zlim[0]  
        ax.plot(X,Y,zs=constArr,color=color,alpha=alpha)


def make3Dlpot(X,Y,Z=None,ax=None,invert='',draw='CXYZ',color='blue'):
    '''
    A simple wrapper around the plotting routines
    
    invert: which axes to invert could be 'X' or 'YZ' or ...
    draw: which lines to draw could be 'C' or 'YZ' or ... 
          C for curve and the rest are projections in the different dimensions 
    '''
    
    if Z == None: Z = fill_Z(X, Y)
    
    ## Now we have all data and need to resample it for a smoother plot
    resFactor = 10
    Xr = scipy.ndimage.zoom(X, resFactor, order=1)
    Yr = scipy.ndimage.zoom(Y, resFactor, order=1)
    Zr = scipy.ndimage.zoom(Z, resFactor, order=1)
    
    
    ## Make the plot
    if ax==None: 
        fig  = plt.figure()
        ax   = fig.gca(projection='3d')
    else:
        fig = ''
    
    if 'X' in invert: ax.invert_xaxis()
    if 'Y' in invert: ax.invert_yaxis()
    if 'Z' in invert: ax.invert_zaxis()
    
    ## Reduce data again and call the plotting wrapper
    stride=np.round(len(Xr)/3.0)
    
    for i in range(len(Xr)):
        if np.mod(i,stride)==0:
            plotLineAndProjection(ax,Xr[i],Yr[i],Zr[i],draw=draw,color=color)
    plotLineAndProjection(ax,Xr[-1],Yr[-1],Zr[-1],draw=draw,color=color)
            
    for i in range(len(Xr[0])):
        if np.mod(i,stride)==0:
            Xi = [row[i] for row in Xr]
            Yi = [row[i] for row in Yr]
            Zi = [row[i] for row in Zr]
            plotLineAndProjection(ax,Xi,Yi,Zi,draw=draw,color=color)
    Xi = [row[-1] for row in Xr]
    Yi = [row[-1] for row in Yr]
    Zi = [row[-1] for row in Zr]
    plotLineAndProjection(ax,Xi,Yi,Zi,draw=draw,color=color)
    
    #cset = ax.contour(X, Y, Z, cmap=cm.coolwarm)
    #ax.clabel(cset, fontsize=9, inline=1)
    
#    cmap = plt.get_cmap('jet')
#    ax.plot_surface(
#      X, Y, Z, 
#      rstride=np.round(resFactor*points/dpi), 
#      cstride=np.round(resFactor*points/dpi), 
#      alpha=0.5, 
#      cmap=cmap,
#      linewidth=0, 
#      #antialiased=False
#      )

#
#    stride=np.round(len(X)/7)
#    
#    if 'C' in draw:
#        ax.plot_wireframe(
#          X, Y, Z, 
#          rstride=stride, 
#          cstride=stride, 
#          alpha=0.5, 
#          color = color,
#          #cmap=cmap,
#          #linewidth=0, 
#          #antialiased=False
#          )
#    else:
#        ax.plot_wireframe(
#          X, Y, Z, 
#          #rstride=stride, 
#          #cstride=stride, 
#          #alpha=0.5, 
#          #cmap=cmap,
#          linewidth=0, 
#          #antialiased=False
#          )    
    ## In case we need a surface plot
    #from matplotlib import cm
    #from matplotlib.ticker import LinearLocator, FormatStrFormatter
    #
    #surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
    #        linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    #
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    #
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    
    
    ## Alternative plotting solution
    #cmap = plt.get_cmap('jet')
    #
    #if PRINT:
    #    ax.plot_surface(
    #      X, Y, Z, 
    #      rstride=100, 
    #      cstride=100, 
    #      alpha=0.5, 
    #      cmap=cmap,
    #      linewidth=0, 
    #      antialiased=True
    #      )
    #else:
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
    
    
    ### Make the plot
    #fig  = plt.figure()
    #ax   = fig.gca(projection='3d')
    ##cset = ax.contour(X, Y, Z, cmap=cm.coolwarm)
    ##ax.clabel(cset, fontsize=9, inline=1)    
    #    
    #ax.plot_wireframe(
    #  X, Y, Z, 
    #  rstride=np.round(resFactor*10), 
    #  cstride=np.round(resFactor*10), 
    #  #alpha=0.5, 
    #  #cmap=cmap,
    #  #linewidth=0, 
    #  #antialiased=False
    #  )
    
    
    ## Plot the two-phase dome and its projections 
    #ax.plot(X_TP,Y_TP,Z_TP,color='black')
    
    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()
    #zlim = ax.get_zlim()
    #
    #constArr = np.ones_like(X_TP) * xlim[0]  
    #ax.plot(constArr,Y_TP,Z_TP,color='black')
    #constArr = np.ones_like(Y_TP) * ylim[1]  
    #ax.plot(X_TP,constArr,Z_TP,color='black')
    #constArr = np.ones_like(Z_TP) * zlim[0]  
    #ax.plot(X_TP,Y_TP,constArr,color='black')
    
    ## Or do we want contour lines?
    ##cset = ax.contour(X, Y, Z, zdir='z', offset=zlim[0], cmap=cmap)
    ##cset = ax.contour(X, Y, Z, zdir='x', offset=xlim[0], cmap=cmap)
    ##cset = ax.contour(X, Y, Z, zdir='y', offset=ylim[1], cmap=cmap)
    #
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    #ax.set_zlim(zlim)
    
    #majorFormatter = matplotlib.ticker.EngFormatter("J/kg")
    #majorFormatter = matplotlib.ticker.
    majorFormatter = matplotlib.ticker.ScalarFormatter()
    majorFormatter.set_scientific(True)
    majorFormatter.set_powerlimits((2,2))
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.zaxis.set_major_formatter(majorFormatter)
    
    return fig, ax, Z






fluid = 'water'
CP.enable_TTSE_LUT(fluid)

PRINT = False

if PRINT:
    points  = 200
    dpi     = 300
    toFile  = True
else:
    points  = 100
    dpi     = 75
    toFile  = False
    

Tmin = CP.PropsSI('Tmin' ,'T',0,'P',0,fluid) + 5.0
Tcri = CP.PropsSI('Tcrit','T',0,'P',0,fluid)
Tmax = CP.PropsSI('Tcrit','T',0,'P',0,fluid) * 2.0 

## Get phase boundary
T_TP = np.linspace(Tmin, Tcri-0.5, 0.5*points)
D_TP = CP.PropsSI('D','T',T_TP,'Q',0,fluid)
P_TP = CP.PropsSI('P','T',T_TP,'Q',0,fluid)
H_TP = CP.PropsSI('H','T',T_TP,'Q',0,fluid)
S_TP = CP.PropsSI('S','T',T_TP,'Q',0,fluid)

D_TP = np.append(D_TP, CP.PropsSI('D','T',T_TP[::-1],'Q',1,fluid))
P_TP = np.append(P_TP, CP.PropsSI('P','T',T_TP[::-1],'Q',1,fluid))
H_TP = np.append(H_TP, CP.PropsSI('H','T',T_TP[::-1],'Q',1,fluid))
S_TP = np.append(S_TP, CP.PropsSI('S','T',T_TP[::-1],'Q',1,fluid))

logV_TP = np.log10(1.0/D_TP)
logP_TP = np.log10(1.0*P_TP)
T_TP = np.append(T_TP, T_TP[::-1])


## These are the default TTSE ranges
D_L = CP.PropsSI('D','T',Tmin,'Q',0,fluid)
D_V = CP.PropsSI('D','T',Tmin,'Q',1,fluid)
H_L = CP.PropsSI('H','T',Tmin,'Q',0,fluid)
H_V = CP.PropsSI('H','T',Tmin,'Q',1,fluid)
P_L = CP.PropsSI('P','T',Tmin,'Q',0,fluid)
P_V = CP.PropsSI('P','T',Tmin,'Q',1,fluid)

Hmin = H_L
Hmax = H_L + (H_V-H_L)*2.0     

#Pmin = CP.PropsSI('ptriple','T',0,'P',0,fluid) + 1
Pmin = np.min([P_L,P_V])
Pmax = CP.PropsSI('pcrit','T',0,'P',0,fluid) * 2.0 # should be p_reduce

Dmin = CP.PropsSI('D','H',Hmin,'P',Pmax,fluid)
Dmax = CP.PropsSI('D','H',Hmax,'P',Pmin,fluid) 


## Set the ranges for the plot 
X_TP = H_TP
Y_TP = logP_TP
Z_TP = S_TP

Xmin = Hmin
Xmax = Hmax
Ymin = np.log10(Pmin)
Ymax = np.log10(Pmax)

X = np.linspace(Xmin, Xmax, points) # Enthalpy
Y = np.linspace(Ymin, Ymax, points) # Pressure
X, Y = np.meshgrid(X, Y)

def get_Z(X_in,Y_in,fluid,out='S'):
    '''
    Just a wrapper to call CoolProp
    '''
    X = X_in
    Y = np.power(10.0,Y_in)
    Z = np.NAN
    try:
        Z = CP.PropsSI(out,'H',X,'P',Y,fluid)
    except(ValueError):
        print "CoolProp failed, returning NAN"
        Z = np.NAN
    return Z

figHPS, axHPS, Z = make3Dlpot(X,Y,invert='Z',draw='XZ')
## Plot the two-phase dome and its projections 
plotLineAndProjection(axHPS,X_TP,Y_TP,Z_TP)

make3Dlpot(X,Y,Z=Z,ax=axHPS,draw='Y')

axHPS.set_xlabel(r'$h$')
axHPS.set_ylabel(r'log $p$')
axHPS.set_zlabel(r'$s$')

## Now we also need thet other variables
get_Z_old = get_Z

def get_Z(X_in,Y_in,fluid):
    return get_Z_old(X_in,Y_in,fluid,out='D')
logVdata = np.log10(1.0/fill_Z(X,Y))

def get_Z(X_in,Y_in,fluid):
    return get_Z_old(X_in,Y_in,fluid,out='T')
Tdata = fill_Z(X,Y)

HPSdict = {'H': X, 'P': Y, 'S': Z, 'V': logVdata, 'T': Tdata}



#########################################################################
# Start with the next diagram
#
#
#
#########################################################################

## Set the ranges for the plot 
X_TP = logV_TP
Y_TP = T_TP
Z_TP = logP_TP

#Xmin = np.log10(1.0/Dmax)
#Xmax = np.log10(1.0/Dmin)
Xmin = np.min(X_TP)
Xmax = np.max(X_TP)
Ymin = np.min(Y_TP)
Ymax = Tmax

X = np.linspace(Xmin, Xmax, points) # Volume
Y = np.linspace(Ymin, Ymax, points) # Temperature
X, Y = np.meshgrid(X, Y)

def get_Z(X_in,Y_in,fluid,out='P'):
    '''
    Just a wrapper to call CoolProp
    '''
    X = 1.0/np.power(10.0,X_in)
    Y = Y_in
    Z = np.NAN
    try:
        Z = np.log10(CP.PropsSI(out,'D',X,'T',Y,fluid))
    except(ValueError):
        print "CoolProp failed, returning NAN"
        Z = np.NAN
    return Z

figVTP, axVTP, Z = make3Dlpot(X,Y,draw='XZ')
## Plot the two-phase dome and its projections 
plotLineAndProjection(axVTP,X_TP,Y_TP,Z_TP)

make3Dlpot(X,Y,Z=Z,ax=axVTP,draw='Y')

axVTP.set_xlabel(r'log $v$')
axVTP.set_ylabel(r'$T$')
axVTP.set_zlabel(r'log $p$')

## Now we also need thet other variables
get_Z_old = get_Z

def get_Z(X_in,Y_in,fluid):
    return get_Z_old(X_in,Y_in,fluid,out='H')

Hdata = fill_Z(X,Y)

def get_Z(X_in,Y_in,fluid):
    return get_Z_old(X_in,Y_in,fluid,out='S')

Sdata = fill_Z(X,Y)

VTPdict = {'V': X, 'T': Y, 'P': Z, 'H': Hdata, 'S': Sdata}



#########################################################################
# Now we have all the data and can start mixing the 
# different definitions
#
#
#########################################################################

#make3Dlpot(fig=figHPS,HPSdict,Y,invert='Z',draw='XZ')

make3Dlpot(HPSdict['V'],HPSdict['T'],Z=HPSdict['P'],ax=axVTP,draw='XYZ',color='red')

make3Dlpot(VTPdict['H'],VTPdict['P'],Z=VTPdict['S'],ax=axHPS,draw='C',color='red')

#make3Dlpot(X,Y,Z=Z,ax=axVTP,draw='C')

#plotLineAndProjection(axVTP,HPSdict['V'][0],HPSdict['T'][0],HPSdict['P'][0])

#
#H_data    = HPSgrid[0]
#logP_data = HPSgrid[1]
#logV_data = np.empty_like(H_data)
#T_data    = np.empty_like(H_data)
#
#   
#for i in range(len(H_data)):
#    for j in range(len(H_data[0])):
#        logV_data[i,j] = np.log10(1.0/get_Z_old(H_data[i,j],logP_data[i,j],fluid,out='D'))
#        T_data[i,j]    = get_Z_old(H_data[i,j],logP_data[i,j],fluid,out='T')
#
#    
#rstride=np.round(len(H_data)/3.0)
#for i in range(len(H_data)):
#    if np.mod(i,rstride)==0:
#        #ax.plot(logV_data[i],T_data[i],logP_data[i],color='red')
#        plotLineAndProjection(ax,logV_data[i],T_data[i],logP_data[i],color='red')
#plotLineAndProjection(ax,logV_data[-1],T_data[-1],logP_data[-1],color='red')
#        
#for i in range(len(H_data[0])):
#    if np.mod(i,rstride)==0:
#        X = [row[i] for row in logV_data]
#        Y = [row[i] for row in T_data]
#        Z = [row[i] for row in logP_data]
#        plotLineAndProjection(ax,X,Y,Z,color='red')
#X = [row[-1] for row in logV_data]
#Y = [row[-1] for row in T_data]
#Z = [row[-1] for row in logP_data]
#plotLineAndProjection(ax,X,Y,Z,color='red')    



#figures=[manager.canvas.figure
#         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
#
#for i, figure in enumerate(figures):
#    figure.savefig('figure%d.png' % i)


#fig.colorbar(surf, shrink=0.5, aspect=5)
if toFile:
    figHPS.savefig('phs.png', dpi = dpi, transparent = True)
    #figHPS.close()
    figVTP.savefig('pvT.png', dpi = dpi, transparent = True)
    #figVTP.close()
    #plt.savefig('TTSE_RANGES.pdf')
else:
    plt.show()
    #figHPS.show()
    #figVTP.show()
