"""
This file implements a psychrometric chart for air at 1 atm
"""

from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.Plots.Plots import InlineLabel 
import matplotlib, numpy, textwrap

import_template=(
"""
import numpy, matplotlib
from CoolProp.HumidAirProp import HAPropsSI
from CoolProp.Plots.Plots import InlineLabel 

p = 101325
Tdb = numpy.linspace(-10,60,100)+273.15

#Make the figure and the axes
fig=matplotlib.pyplot.figure(figsize=(10,8))
ax=fig.add_axes((0.1,0.1,0.85,0.85))
"""
)

closure_template=(
"""
matplotlib.pyplot.show()
"""                  
)

Tdb = numpy.linspace(-10,60,100)+273.15
p = 101325

class PlotFormatting(object):
    
    def plot(self,ax):
        ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
        ax.set_ylim(0,0.03)
        ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
        ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")
        
    def __str__(self):
        return textwrap.dedent("""
            ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
            ax.set_ylim(0,0.03)
            ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
            ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")
            """)
        
class SaturationLine(object):
    
    def plot(self,ax):
        w = [HAPropsSI('W','T',T,'P',p,'R',1.0) for T in Tdb]
        ax.plot(Tdb-273.15,w,lw=2)
        
    def __str__(self):
        return textwrap.dedent("""
               # Saturation line
               w = [HAPropsSI('W','T',T,'P',p,'R',1.0) for T in Tdb]
               ax.plot(Tdb-273.15,w,lw=2)
               """
               )
        
class HumidityLabels(object):
    def __init__(self,RH_values,h):
        self.RH_values = RH_values
        self.h = h
    
    def plot(self,ax):
        xv = Tdb #[K]
        for RH in self.RH_values:
            yv = [HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb]
            y = HAPropsSI('W','P',p,'H',self.h,'R',RH)
            T_K,w,rot = InlineLabel(xv, yv, y=y, axis = ax)
            string = r'$\phi$='+'{s:0.0f}'.format(s=RH*100)+'%'
            #Make a temporary label to get its bounding box
            bbox_opts = dict(boxstyle='square,pad=0.0',fc='white',ec='None',alpha = 0.5)
            ax.text(T_K-273.15,w,string,rotation = rot,ha ='center',va='center',bbox=bbox_opts)
    
    def __str__(self):
        return textwrap.dedent("""
                xv = Tdb #[K]
                for RH in {RHValues:s}:
                    yv = [HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb]
                    y = HAPropsSI('W','P',p,'H',{h:f},'R',RH)
                    T_K,w,rot = InlineLabel(xv, yv, y=y, axis = ax)
                    string = r'$\phi$='+{s:s}+'%'
                    bbox_opts = dict(boxstyle='square,pad=0.0',fc='white',ec='None',alpha = 0.5)
                    ax.text(T_K-273.15,w,string,rotation = rot,ha ='center',va='center',bbox=bbox_opts)
                """.format(h=self.h, RHValues=str(self.RH_values), s="'{s:0.0f}'.format(s=RH*100)")
                )
        
class HumidityLines(object):
    
    def __init__(self,RH_values):
        self.RH_values = RH_values
        
    def plot(self,ax):
        for RH in self.RH_values:
            w = [HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb]
            ax.plot(Tdb-273.15,w,'r',lw=1)
        
    def __str__(self):
        return textwrap.dedent("""
               # Humidity lines
               RHValues = {RHValues:s}
               for RH in RHValues:
                   w = [HAPropsSI('W','T',T,'P',p,'R',RH) for T in Tdb]
                   ax.plot(Tdb-273.15,w,'r',lw=1)
               """.format(RHValues=str(self.RH_values))
               )

class EnthalpyLines(object):
    
    def __init__(self,H_values):
        self.H_values = H_values
        
    def plot(self,ax):
        for H in self.H_values:
            #Line goes from saturation to zero humidity ratio for this enthalpy
            T1 = HAPropsSI('T','H',H,'P',p,'R',1.0)-273.15
            T0 = HAPropsSI('T','H',H,'P',p,'R',0.0)-273.15
            w1 = HAPropsSI('W','H',H,'P',p,'R',1.0)
            w0 = HAPropsSI('W','H',H,'P',p,'R',0.0)
            ax.plot(numpy.r_[T1,T0],numpy.r_[w1,w0],'r',lw=1)
        
    def __str__(self):
        return textwrap.dedent("""
               # Humidity lines
               for H in {HValues:s}:
                   #Line goes from saturation to zero humidity ratio for this enthalpy
                   T1 = HAPropsSI('T','H',H,'P',p,'R',1.0)-273.15
                   T0 = HAPropsSI('T','H',H,'P',p,'R',0.0)-273.15
                   w1 = HAPropsSI('W','H',H,'P',p,'R',1.0)
                   w0 = HAPropsSI('W','H',H,'P',p,'R',0.0)
                   ax.plot(numpy.r_[T1,T0],numpy.r_[w1,w0],'r',lw=1)
               """.format(HValues=str(self.H_values))
               )
    
if __name__=='__main__':
    fig=matplotlib.pyplot.figure(figsize=(10,8))
    ax=fig.add_axes((0.1,0.1,0.85,0.85))
    ax.set_xlim(Tdb[0]-273.15,Tdb[-1]-273.15)
    ax.set_ylim(0,0.03)
    ax.set_xlabel(r"Dry bulb temperature [$^{\circ}$C]")
    ax.set_ylabel(r"Humidity ratio ($m_{water}/m_{dry\ air}$) [-]")

    SL = SaturationLine()
    SL.plot(ax)

    RHL = HumidityLines([0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    RHL.plot(ax)
    
    RHLabels = HumidityLabels([0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], h=65000)
    RHLabels.plot(ax)
     
    HL = EnthalpyLines(range(-20000,100000,10000))
    HL.plot(ax)
    
    PF = PlotFormatting()
    PF.plot(ax)

    matplotlib.pyplot.show()

    fp = open('PsychScript.py','w')
    for chunk in [import_template,SL,RHL,HL,PF,RHLabels,closure_template]:
        fp.write(str(chunk))
    fp.close()
    execfile('PsychScript.py')
