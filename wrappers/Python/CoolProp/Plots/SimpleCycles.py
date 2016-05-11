# -*- coding: utf-8 -*-

from __future__ import print_function, division
from six import with_metaclass

import matplotlib,numpy

import numpy as np

import CoolProp
from CoolProp.CoolProp import PropsSI
from .Common import BasePlot, process_fluid_state, PropertyDict, SIunits
import warnings
from abc import ABCMeta


def SimpleCycle(Ref,Te,Tc,DTsh,DTsc,eta_a,Ts_Ph='Ph',skipPlot=False,axis=None):
    """
    This function plots a simple four-component cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : A string for the refrigerant
    * Te : Evap Temperature in K
    * Tc : Condensing Temperature in K
    * DTsh : Evaporator outlet superheat in K
    * DTsc : Condenser outlet subcooling in K
    * eta_a : Adiabatic efficiency of compressor (no units) in range [0,1]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """
    
    warnings.warn("This function has been deprecated. Please consider converting it to an object inheriting from \"BaseCycle\".",DeprecationWarning)
    from scipy.optimize import newton
    
    T=numpy.zeros((6))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    T2s=newton(lambda T: PropsSI('S','T',T,'P',pc,Ref)-s[1],T[1]+30)
    h2s=PropsSI('H','T',T2s,'P',pc,Ref)
    h[2]=(h2s-h[1])/eta_a+h[1]
    T[2]=PropsSI('T','H',h[2],'P',pc,Ref)
    s[2]=PropsSI('S','T',T[2],'P',pc,Ref)

    sbubble_c=PropsSI('S','P',pc,'Q',0,Ref)
    sdew_c=PropsSI('S','P',pc,'Q',1,Ref)
    sbubble_e=PropsSI('S','P',pe,'Q',0,Ref)
    sdew_e=PropsSI('S','P',pe,'Q',1,Ref)
    T[3]=Tc-DTsc
    h[3]=PropsSI('H','T',T[3],'P',pc,Ref)
    s[3]=PropsSI('S','T',T[3],'P',pc,Ref)
    h[4]=h[3]
    h[5]=h[1]
    s[5]=s[1]
    T[5]=T[1]
    p=[numpy.nan,pe,pc,pc,pe,pe]
    COP=(h[1]-h[4])/(h[2]-h[1])
    COPH=(h[2]-h[3])/(h[2]-h[1])

    hsatL=PropsSI('H','T',Te,'Q',0,Ref)
    hsatV=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL=PropsSI('S','T',Te,'Q',0,Ref)
    ssatV=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL=1/PropsSI('D','T',Te,'Q',0,Ref)
    vsatV=1/PropsSI('D','T',Te,'Q',1,Ref)
    x=(h[4]-hsatL)/(hsatV-hsatL)
    s[4]=x*ssatV+(1-x)*ssatL
    T[4]=x*Te+(1-x)*Te

    print(COP,COPH)
    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
        elif Ts_Ph in ['Ts','ts']:
            s=list(s)
            T=list(T)
            s.insert(5,sdew_e)
            T.insert(5,Te)
            s.insert(3,sbubble_c)
            T.insert(3,Tc)
            s.insert(3,sdew_c)
            T.insert(3,Tc)
            ax.plot(s[1::],T[1::],'b')
        else:
            raise TypeError('Type of Ts_Ph invalid')

def TwoStage(Ref,Q,Te,Tc,DTsh,DTsc,eta_oi,f_p,Tsat_ic,DTsh_ic,Ts_Ph='Ph',prints=False,skipPlot=False,axis=None,**kwargs):
    """
    This function plots a two-stage cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : Refrigerant [string]
    * Q : Cooling capacity [W]
    * Te : Evap Temperature [K]
    * Tc : Condensing Temperature [K]
    * DTsh : Evaporator outlet superheat [K]
    * DTsc : Condenser outlet subcooling [K]
    * eta_oi : Adiabatic efficiency of compressor (no units) in range [0,1]
    * f_p : fraction of compressor power lost as ambient heat transfer in range [0,1]
    * Tsat_ic : Saturation temperature corresponding to intermediate pressure [K]
    * DTsh_ic : Superheating at outlet of intermediate stage [K]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * prints : True to print out some values
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """
    
    warnings.warn("This function has been deprecated. PLease consider converting it to an object inheriting from \"BaseCycle\".",DeprecationWarning)

    T=numpy.zeros((8))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    rho=numpy.zeros_like(T)
    T[0]=numpy.NAN
    s[0]=numpy.NAN
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    pic=PropsSI('P','T',Tsat_ic,'Q',1.0,Ref)
    Tbubble_c=PropsSI('T','P',pc,'Q',0,Ref)
    Tbubble_e=PropsSI('T','P',pe,'Q',0,Ref)

    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    rho[1]=PropsSI('D','T',T[1],'P',pe,Ref)
    T[5]=Tbubble_c-DTsc
    h[5]=PropsSI('H','T',T[5],'P',pc,Ref)
    s[5]=PropsSI('S','T',T[5],'P',pc,Ref)
    rho[5]=PropsSI('D','T',T[5],'P',pc,Ref)
    mdot=Q/(h[1]-h[5])

    rho1=PropsSI('D','T',T[1],'P',pe,Ref)
    h2s=PropsSI('H','S',s[1],'P',pic,Ref)
    Wdot1=mdot*(h2s-h[1])/eta_oi
    h[2]=h[1]+(1-f_p)*Wdot1/mdot
    T[2]=PropsSI('T','H',h[2],'P',pic,Ref)
    s[2]=PropsSI('S','T',T[2],'P',pic,Ref)
    rho[2]=PropsSI('D','T',T[2],'P',pic,Ref)
    T[3]=288
    p[3]=pic
    h[3]=PropsSI('H','T',T[3],'P',pic,Ref)
    s[3]=PropsSI('S','T',T[3],'P',pic,Ref)
    rho[3]=PropsSI('D','T',T[3],'P',pic,Ref)
    rho3=PropsSI('D','T',T[3],'P',pic,Ref)
    h4s=PropsSI('H','T',s[3],'P',pc,Ref)
    Wdot2=mdot*(h4s-h[3])/eta_oi
    h[4]=h[3]+(1-f_p)*Wdot2/mdot
    T[4]=PropsSI('T','H',h[4],'P',pc,Ref)
    s[4]=PropsSI('S','T',T[4],'P',pc,Ref)
    rho[4]=PropsSI('D','T',T[4],'P',pc,Ref)

    sbubble_e=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    sbubble_c=PropsSI('S','T',Tbubble_c,'Q',0,Ref)
    sdew_e=PropsSI('S','T',Te,'Q',1,Ref)
    sdew_c=PropsSI('S','T',Tc,'Q',1,Ref)

    hsatL=PropsSI('H','T',Tbubble_e,'Q',0,Ref)
    hsatV=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    ssatV=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL=1/PropsSI('D','T',Tbubble_e,'Q',0,Ref)
    vsatV=1/PropsSI('D','T',Te,'Q',1,Ref)
    x=(h[5]-hsatL)/(hsatV-hsatL)
    s[6]=x*ssatV+(1-x)*ssatL
    T[6]=x*Te+(1-x)*Tbubble_e
    rho[6]=1.0/(x*vsatV+(1-x)*vsatL)

    h[6]=h[5]
    h[7]=h[1]
    s[7]=s[1]
    T[7]=T[1]
    p=[numpy.nan,pe,pic,pic,pc,pc,pe,pe]
    COP=Q/(Wdot1+Wdot2)
    RE=h[1]-h[6]

    if prints==True:
        print('x5:',x)
        print('COP:', COP)
        print('COPH', (Q+Wdot1+Wdot2)/(Wdot1+Wdot2))
        print(T[2]-273.15,T[4]-273.15,p[2]/p[1],p[4]/p[3])
        print(mdot,mdot*(h[4]-h[5]),pic)
        print('Vdot1',mdot/rho1,'Vdisp',mdot/rho1/(3500/60.)*1e6/0.7)
        print('Vdot2',mdot/rho3,'Vdisp',mdot/rho3/(3500/60.)*1e6/0.7)
        print(mdot*(h[4]-h[5]),Tc-273.15)
        for i in range(1,len(T)-1):
            print('%d & %g & %g & %g & %g & %g \\\\' %(i,T[i]-273.15,p[i],h[i],s[i],rho[i]))
    else:
        print(Tsat_ic,COP)

    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        else:
            ax=axis
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
        elif Ts_Ph in ['Ts','ts']:
            s_copy=s.copy()
            T_copy=T.copy()
            for i in range(1,len(s)-1):
                ax.plot(s[i],T[i],'bo',mfc='b',mec='b')
                dT=[0,-5,5,-20,5,5,5]
                ds=[0,0.05,0,0,0,0,0]
                ax.text(s[i]+ds[i],T[i]+dT[i],str(i))

            s=list(s)
            T=list(T)
            s.insert(7,sdew_e)
            T.insert(7,Te)
            s.insert(5,sbubble_c)
            T.insert(5,Tbubble_c)
            s.insert(5,sdew_c)
            T.insert(5,Tc)

            ax.plot(s,T)
            s=s_copy
            T=T_copy
        else:
            raise TypeError('Type of Ts_Ph invalid')
    return COP

def EconomizedCycle(Ref,Qin,Te,Tc,DTsh,DTsc,eta_oi,f_p,Ti,Ts_Ph='Ts',skipPlot=False,axis=None,**kwargs):
    """
    This function plots an economized cycle, on the current axis, or that given by the optional parameter *axis*

    Required parameters:

    * Ref : Refrigerant [string]
    * Qin : Cooling capacity [W]
    * Te : Evap Temperature [K]
    * Tc : Condensing Temperature [K]
    * DTsh : Evaporator outlet superheat [K]
    * DTsc : Condenser outlet subcooling [K]
    * eta_oi : Adiabatic efficiency of compressor (no units) in range [0,1]
    * f_p : fraction of compressor power lost as ambient heat transfer in range [0,1]
    * Ti : Saturation temperature corresponding to intermediate pressure [K]

    Optional parameters:

    * Ts_Ph : 'Ts' for a Temperature-Entropy plot, 'Ph' for a Pressure-Enthalpy
    * axis : An axis to use instead of the active axis
    * skipPlot : If True, won't actually plot anything, just print COP

    """
    
    warnings.warn("This function has been deprecated. Please consider converting it to an object inheriting from \"BaseCycle\".",DeprecationWarning)
    from scipy.optimize import newton

    m=1

    T=numpy.zeros((11))
    h=numpy.zeros_like(T)
    p=numpy.zeros_like(T)
    s=numpy.zeros_like(T)
    rho=numpy.zeros_like(T)

    T[0]=numpy.NAN
    s[0]=numpy.NAN
    T[1]=Te+DTsh
    pe=PropsSI('P','T',Te,'Q',1.0,Ref)
    pc=PropsSI('P','T',Tc,'Q',1.0,Ref)
    pi=PropsSI('P','T',Ti,'Q',1.0,Ref)
    p[1]=pe
    h[1]=PropsSI('H','T',T[1],'P',pe,Ref)
    s[1]=PropsSI('S','T',T[1],'P',pe,Ref)
    rho[1]=PropsSI('D','T',T[1],'P',pe,Ref)
    h2s=PropsSI('H','S',s[1],'P',pi,Ref)
    wdot1=(h2s-h[1])/eta_oi
    h[2]=h[1]+(1-f_p[0])*wdot1
    p[2]=pi
    #T[2]=T_hp(Ref,h[2],pi,T2s)
    T[2]=PropsSI('T','H',h[2],'P',pi,Ref)
    
    
    s[2]=PropsSI('S','T',T[2],'P',pi,Ref)
    rho[2]=PropsSI('D','T',T[2],'P',pi,Ref)

    T[5]=Tc-DTsc
    h[5]=PropsSI('H','T',T[5],'P',pc,Ref)
    s[5]=PropsSI('S','T',T[5],'P',pc,Ref)
    rho[5]=PropsSI('D','T',T[5],'P',pc,Ref)

    p[5]=pc
    p[6]=pi
    h[6]=h[5]

    p[7]=pi
    p[8]=pi
    p[6]=pi
    T[7]=Ti
    h[7]=PropsSI('H','T',Ti,'Q',1,Ref)
    s[7]=PropsSI('S','T',Ti,'Q',1,Ref)
    rho[7]=PropsSI('D','T',Ti,'Q',1,Ref)
    T[8]=Ti
    h[8]=PropsSI('H','T',Ti,'Q',0,Ref)
    s[8]=PropsSI('S','T',Ti,'Q',0,Ref)
    rho[8]=PropsSI('D','T',Ti,'Q',0,Ref)
    x6=(h[6]-h[8])/(h[7]-h[8]) #Vapor Quality
    s[6]=s[7]*x6+s[8]*(1-x6)
    rho[6]=1.0/(x6/rho[7]+(1-x6)/rho[8])
    T[6]=Ti

    #Injection mass flow rate
    x=m*(h[6]-h[8])/(h[7]-h[6])


    p[3]=pi
    h[3]=(m*h[2]+x*h[7])/(m+x)
    #T[3]=T_hp(Ref,h[3],pi,T[2])
    T[3]=PropsSI('T','H',h[3],'P',pi,Ref)
    s[3]=PropsSI('S','T',T[3],'P',pi,Ref)
    rho[3]=PropsSI('D','T',T[3],'P',pi,Ref)
    T4s=newton(lambda T: PropsSI('S','T',T,'P',pc,Ref)-s[3],T[2]+30)
    h4s=PropsSI('H','T',T4s,'P',pc,Ref)
    p[4]=pc
    wdot2=(h4s-h[3])/eta_oi
    h[4]=h[3]+(1-f_p[1])*wdot2
    #T[4]=T_hp(Ref,h[4],pc,T4s)
    T[4]=PropsSI('T','H',h[4],'P',pc,Ref)
    s[4]=PropsSI('S','T',T[4],'P',pc,Ref)
    rho[4]=PropsSI('D','T',T[4],'P',pc,Ref)

    p[9]=pe
    h[9]=h[8]
    T[9]=Te
    hsatL_e=PropsSI('H','T',Te,'Q',0,Ref)
    hsatV_e=PropsSI('H','T',Te,'Q',1,Ref)
    ssatL_e=PropsSI('S','T',Te,'Q',0,Ref)
    ssatV_e=PropsSI('S','T',Te,'Q',1,Ref)
    vsatL_e=1/PropsSI('D','T',Te,'Q',0,Ref)
    vsatV_e=1/PropsSI('D','T',Te,'Q',1,Ref)
    x9=(h[9]-hsatL_e)/(hsatV_e-hsatL_e) #Vapor Quality
    s[9]=ssatV_e*x9+ssatL_e*(1-x9)
    rho[9]=1.0/(x9*vsatV_e+(1-x9)*vsatL_e)

    s[10]=s[1]
    T[10]=T[1]
    h[10]=h[1]
    p[10]=p[1]

    Tbubble_e=Te
    Tbubble_c=Tc
    sbubble_e=PropsSI('S','T',Tbubble_e,'Q',0,Ref)
    sbubble_c=PropsSI('S','T',Tbubble_c,'Q',0,Ref)
    sdew_e=PropsSI('S','T',Te,'Q',1,Ref)
    sdew_c=PropsSI('S','T',Tc,'Q',1,Ref)

    Wdot1=m*wdot1
    Wdot2=(m+x)*wdot2
    if skipPlot==False:
        if axis==None:
            ax=matplotlib.pyplot.gca()
        else:
            ax=axis
        if Ts_Ph in ['ph','Ph']:
            ax.plot(h,p)
            ax.set_yscale('log')
        elif Ts_Ph in ['Ts','ts']:
            ax.plot(numpy.r_[s[7],s[3]],numpy.r_[T[7],T[3]],'b')
            s_copy=s.copy()
            T_copy=T.copy()
            dT=[0,-5,5,-12,5,12,-12,0,0,0]
            ds=[0,0.05,0.05,0,0.05,0,0.0,0.05,-0.05,-0.05]
            for i in range(1,len(s)-1):
                ax.plot(s[i],T[i],'bo',mfc='b',mec='b')
                ax.text(s[i]+ds[i],T[i]+dT[i],str(i),ha='center',va='center')

            s=list(s)
            T=list(T)
            s.insert(10,sdew_e)
            T.insert(10,Te)
            s.insert(5,sbubble_c)
            T.insert(5,Tbubble_c)
            s.insert(5,sdew_c)
            T.insert(5,Tc)
            ax.plot(s,T,'b')

            s=s_copy
            T=T_copy
        else:
            raise TypeError('Type of Ts_Ph invalid')

    COP=m*(h[1]-h[9])/(m*(h[2]-h[1])+(m+x)*(h[4]-h[3]))
    for i in range(1,len(T)-1):
            print('%d & %g & %g & %g & %g & %g \\\\' %(i,T[i]-273.15,p[i],h[i],s[i],rho[i]))
    print(x,m*(h[1]-h[9]),(m*(h[2]-h[1])+(m+x)*(h[4]-h[3])),COP)
    mdot=Qin/(h[1]-h[9])
    mdot_inj=x*mdot
    print('x9',x9,)
    print('Qcond',(mdot+mdot_inj)*(h[4]-h[5]),'T4',T[4]-273.15)
    print(mdot,mdot+mdot_inj)
    f=3500/60.
    eta_v=0.7
    print('Vdisp1: ',mdot/(rho[1]*f*eta_v)*1e6,'cm^3')
    print('Vdisp2: ',(mdot+mdot_inj)/(rho[1]*f*eta_v)*1e6,'cm^3')
    return COP


    #class SimpleCycle(object):
    #    """A class that calculates a simple thermodynamic cycle"""
    #    def __init__(self, *args, **kwargs):
    #        object.__init__(self, *args, **kwargs)
    # (states, steps, fluid):
    
#             Parameters
#         ----------        
#         x_type : int, str
#             Either a letter or an integer that specifies the property type for the x-axis
#         y_type : int, str
#             Either a letter or an integer that specifies the property type for the y-axis
#         states : list 
#             A collection of state points that follows a fixed scheme defined 
#             in the implementing subclass.
#         fluid_ref : str, CoolProp.AbstractState
#             The fluid property provider, either a subclass of CoolProp.AbstractState
#             or a string that can be used to generate a CoolProp.AbstractState instance
#             via :func:`Common.process_fluid_state`.
#         steps : int
#             The number of steps used for going from one state to another
#         
#         for more properties, see :class:`CoolProp.Plots.Common.Base2DObject`.  

# # See http://stackoverflow.com/questions/1061283/lt-instead-of-cmp
# class ComparableMixin:
#     """A mixin class that implements all comparing mathods except for __lt__"""
#     def __eq__(self, other):
#         return not self<other and not other<self
#     def __ne__(self, other):
#         return self<other or other<self
#     def __gt__(self, other):
#         return other<self
#     def __ge__(self, other):
#         return not self<other
#     def __le__(self, other):
#         return not other<self

class StatePoint(PropertyDict):
    """A simple fixed dimension dict represented by an object with attributes"""
    
    # Significant digits in SI units
    ROUND_DECIMALS = {
      CoolProp.iDmass : 5,
      CoolProp.iHmass : 5,
      CoolProp.iP     : 2,
      CoolProp.iSmass : 5,
      CoolProp.iT     : 5,
      CoolProp.iUmass : 5,
      CoolProp.iQ     : 5
    }
    
    def __iter__(self):
        """Make sure we always iterate in the same order"""
        keys = [CoolProp.iDmass,CoolProp.iHmass,CoolProp.iP,CoolProp.iSmass,CoolProp.iT]
        for key in sorted(keys):
            yield key
    
    def __str__(self):        
        return str(self.__dict__)
    
    def __prop_compare(self,other,typ):
        # TODO
        if   self[typ] is None and other[typ] is None: return 0
        elif self[typ] is None and other[typ] is not None: return -1
        elif self[typ] is not None and other[typ] is None: return 1
        else:
            A = np.round(self[typ] , self.ROUND_DECIMALS[typ])
            B = np.round(other[typ], self.ROUND_DECIMALS[typ])
            if A>B: return 1
            elif A<B: return -1
            elif A==B: return 0
            else: raise ValueError("Comparison failed.")

    def __eq__(self, other):
        for i in self:
            if not self.__prop_compare(other,i) == 0:
                return False 
        return True
    
    def __hash__(self):
        return hash(repr(self))
    

class StateContainer(object):
    """A collection of values for the main properties, built to mixin with :class:`CoolProp.Plots.Common.PropertyDict`
    
    Examples
    --------
    This container has overloaded accessor methods. Just pick your own flavour 
    or mix the styles as you like:
    
    >>> from __future__ import print_function
    >>> import CoolProp
    >>> from CoolProp.PLots.SimpleCycles import StateContainer
    >>> T0 = 300.000; p0 = 200000.000; h0 = 112745.749; s0 = 393.035
    >>> cycle_states = StateContainer()
    >>> cycle_states[0,'H'] = h0
    >>> cycle_states[0]['S'] = s0
    >>> cycle_states[0][CoolProp.iP] = p0
    >>> cycle_states[0,CoolProp.iT] = T0
    >>> cycle_states[1,"T"] = 300.064
    >>> print(cycle_states)
    Stored State Points:
    state        T (K)       p (Pa)    ρ (kg/m³)     h (J/kg)   s (J/kg/K)
        0      300.000   200000.000      996.601   112745.749      393.035
        1      300.064            -            -            -            -
    
    """
    
    def __init__(self,unit_system=SIunits()):
        self._points = {}
        self._units  = unit_system
    
    @property
    def points(self): return self._points
    @points.setter
    def points(self, value): self._points = value
    
    @property
    def units(self): return self._units
    @units.setter
    def units(self, value): self._units = value
    
    def get_point(self, index, SI=True):
        if SI: 
            state = self[index]
        else:
            state = self[index]
            for i in state: 
                state[i] = self.units[i].from_SI(state[i])
        return state 
    
    def set_point(self, index, value, SI=True):
        if SI: 
            self._points[index] = value
        else:
            for i in value: 
                self._points[index][i] = self.units[i].to_SI(value[i])

    def _list_like(self, value):
        """Try to detect a list-like structure excluding strings"""
        return (not hasattr(value, "strip") and
            (hasattr(value, "__getitem__") or
            hasattr(value, "__iter__")))
        # return is_sequence(value) # use from pandas.core.common import is_sequence
    
    def __len__(self):
        """Some cheating to get the correct behaviour"""
        return len(self._points)
    
    def __iter__(self):
        """Make sure we iterate in the righ order"""
        for key in sorted(self._points):
            yield key
    
    def __getitem__(self, index):
        """Another tweak that changes the default access path"""
        if self._list_like(index):
            if len(index)==0: raise IndexError("Received empty index.")
            elif len(index)==1: return self._points[index[0]]
            elif len(index)==2: return self._points[index[0]][index[1]]
            else: raise IndexError("Received too long index.")
        return self._points[index]
    
    def __setitem__(self, index, value):
        """Another tweak that changes the default access path"""
        if self._list_like(index):
            if len(index)==0: raise IndexError("Received empty index.")
            elif len(index)==1: self._points[index[0]]           = value
            elif len(index)==2:
                # safeguard against empty entries
                if index[0] not in self._points:
                    self._points[index[0]] = StatePoint() 
                self._points[index[0]][index[1]] = value
            else: raise IndexError("Received too long index.")
        else:
            self._points[index] = value 
        
    def __str__(self):
        out = "Stored State Points:\n"
        keys = True 
        for i in self._points:
            if keys:
                row = ["{0:>5s}".format("state")]
                for j in self._points[i]:
                    label = u"{0:s} ({1:s})".format(self.units[j].symbol,self.units[j].unit)
                    row.append(u"{0:>11s}".format(label))
                out = out + "  ".join(row) + "\n"
                keys = False
            row = ["{0:>5s}".format(str(i))]
            for j in self._points[i]:
                try:
                    row.append(u"{0:11.3f}".format(self.units[j].from_SI(self._points[i][j])))
                except:
                    row.append(u"{0:>11s}".format("-"))
            out = out + "  ".join(row) + "\n"
        return out.encode('utf8', 'replace')
    
    def append(self,new):
        i = 0 + self.__len__()
        for j in new:
            self[i,j] = new[j]
        return self
    
    def extend(self,new):
        i = 0 + self.__len__()
        for j in new:
            for k in new[j]:
                self[i,k] = new[j][k]
            i = i +1
        return self 
            
    @property
    def D(self): return np.array([self._points[k].D for k in self])
    @property
    def H(self): return np.array([self._points[k].H for k in self])
    @property
    def P(self): return np.array([self._points[k].P for k in self])
    @property
    def S(self): return np.array([self._points[k].S for k in self])
    @property
    def T(self): return np.array([self._points[k].T for k in self])
    @property
    def U(self): return np.array([self._points[k].U for k in self])
    @property
    def Q(self): return np.array([self._points[k].Q for k in self])
        
    

class BaseCycle(BasePlot):
    """A simple thermodynamic cycle, should not be used on its own."""
    
    # Define the iteration keys
    PROPERTIES = {
      CoolProp.iDmass : 'density', 
      CoolProp.iHmass : 'specific enthalpy', 
      CoolProp.iP     : 'pressure', 
      CoolProp.iSmass : 'specific entropy', 
      CoolProp.iT     : 'temperature'
    }
    
    STATECOUNT=0
    """A list of accepted numbers of states"""
    
    STATECHANGE=None
    """A list of lists of tuples that defines how the state transitions 
    behave for the corresponding entry in BaseCycle.STATECOUNT"""
    
    def __init__(self, fluid_ref, graph_type, unit_system='EUR', **kwargs):
        """Initialises a simple cycle calculator
        
        Parameters
        ----------        
        fluid_ref : str, CoolProp.AbstractState
            The fluid property provider, either a subclass of CoolProp.AbstractState
            or a string that can be used to generate a CoolProp.AbstractState instance
            via :func:`Common.process_fluid_state`.
        graph_type : string
            The graph type to be plotted, like \"PH\" or \"TS\"
        unit_system : string, ['EUR','KSI','SI']
            Select the units used for the plotting.  'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
                
        for more properties, see :class:`CoolProp.Plots.Common.BasePlot`.        
        """
        self._cycle_states = StateContainer()
        self._steps = 2
        BasePlot.__init__(self, fluid_ref, graph_type, unit_system, **kwargs)
        
    
    @property
    def cycle_states(self): return self._cycle_states
    @cycle_states.setter
    def cycle_states(self, value):
        if len(value) != self.STATECOUNT:
            raise ValueError("Your number of states ({0:d}) is not in the list of allowed state counts: {1:s}.".format(len(value),str(self.STATECOUNT)))
        self._cycle_states = value
    
    @property
    def steps(self): return self._steps
    @steps.setter
    def steps(self, value): self._steps = int(max([value,2]))

    @BasePlot.system.setter
    def system(self, value): 
        if value in self.UNIT_SYSTEMS:
            self._system = self.UNIT_SYSTEMS[value]
        elif isinstance(value, PropertyDict):
            self._system = value
        else:
            raise ValueError("Invalid unit_system input \"{0:s}\", expected a string from {1:s}".format(str(value),str(self.UNIT_SYSTEMS.keys())))
        self._cycle_states._units = self._system
    
    
    def valid_states(self):
        """Check the formats of BaseCycle.STATECOUNT and BaseCycle.STATECHANGE"""
        if len(self.STATECHANGE) != self.STATECOUNT: 
            raise ValueError("Invalid number of states and or state change operations")
        return True 
    
    def fill_states(self,objs=None):
        """Try to populate all fields in the state objects"""
        
        if objs is None: 
            objs = self._cycle_states
            local = True
        else: 
            local = False
        
        for i in objs:
            full = True
            for j in objs[i]:
                if objs[i][j] is None:
                    full = False
            if full: continue
            if (objs[i][CoolProp.iDmass] is not None and 
              objs[i][CoolProp.iT] is not None):
                self._state.update(CoolProp.DmassT_INPUTS, objs[i][CoolProp.iDmass], objs[i][CoolProp.iT])
            elif (objs[i][CoolProp.iP] is not None and 
              objs[i][CoolProp.iHmass] is not None):
                self._state.update(CoolProp.HmassP_INPUTS, objs[i][CoolProp.iHmass], objs[i][CoolProp.iP])
            elif (objs[i][CoolProp.iP] is not None and 
              objs[i][CoolProp.iSmass] is not None):
                self._state.update(CoolProp.PSmass_INPUTS, objs[i][CoolProp.iP], objs[i][CoolProp.iSmass])
            else:
                warnings.warn("Please fill the state[{0:s}] manually.".format(str(i)))
                continue
            for j in objs[i]:
                if objs[i][j] is None:
                    objs[i][j] = self._state.keyed_output(j)
        
        if local: self._cycle_states = objs
        return objs

    
    
    def state_change(self,in1,in2,start,ty1='lin',ty2='lin'):
        """Calculates a state change defined by the properties in1 and in2
        
        Uses self.states[start] and self.states[start+1] (or self.states[0]) to define 
        the process and interpolates between the values. 
        
        Parameters
        ----------
        in1 : int 
            The index of the first defined property.
        in2 : int 
            The index of the second defined property.
        start : int 
            The index of the start state. 
        ty1 : str
            The key that defines the type of state change for in1, lin or log.
        ty2 : str 
            The key that defines the type of state change for in2, lin or log.
            
        Returns
        -------
        scalar or array_like 
            a list of the length of self.steps+1 that describes the process. It includes start and end state.
        """
        self.fill_states()
        end = start + 1
        if end >= len(self.cycle_states): end -= len(self.cycle_states) 
        start = self.cycle_states[start]
        end = self.cycle_states[end]
        #
        val = []
        inv = [in1,in2]
        typ = [ty1,ty2]
        for i,v in enumerate(inv):
            if typ[i] == 'lin':
                val.append(np.linspace(start[v], end[v], self.steps))
            elif typ[i] == 'log':
                val.append(np.logspace(np.log10(start[v]), np.log10(end[v]), self.steps))
            else:
                raise ValueError("Unknow range generator {0:s}".format(str(typ[i])))
        
        sc = StateContainer(self._system)
        for i,_ in enumerate(val[0]):
            sc[i,inv[0]] = val[0][i]
            sc[i,inv[1]] = val[1][i]
            
        return self.fill_states(sc)     
        
    def get_state_change(self, index):
        return self.STATECHANGE[index](self)
    
    def get_state_changes(self):
        sc = self.get_state_change(0)
        for i in range(1,self.STATECOUNT):
            sc.extend(self.get_state_change(i))
        return sc 
                
        
class BasePowerCycle(BaseCycle):
    """A thermodynamic cycle for power producing processes.
    
    Defines the basic properties and methods to unify access to 
    power cycle-related quantities. 
    """
    
    def __init__(self, fluid_ref='HEOS::Water', graph_type='TS', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BaseCycle` for details."""
        BaseCycle.__init__(self, fluid_ref, graph_type, **kwargs)
    
    def eta_carnot(self):
        """Carnot efficiency
         
        Calculates the Carnot efficiency for the specified process, :math:`\eta_c = 1 - \frac{T_c}{T_h}`.
         
        Returns
        -------
        float
        """
        Tvector = self._cycle_states.T
        return 1. - np.min(Tvector) / np.max(Tvector)
    
    def eta_thermal(self):
        """Thermal efficiency
         
        The thermal efficiency for the specified process(es), :math:`\eta_{th} = \frac{\dot{W}_{exp} - \dot{W}_{pum}}{\dot{Q}_{in}}`.
         
        Returns
        -------
        float
        """
        raise NotImplementedError("Implement it in the subclass.")


class SimpleRankineCycle(BasePowerCycle):
    """A simple Rankine cycle *without* regeneration"""
    STATECOUNT=4
    STATECHANGE=[
      lambda inp: BaseCycle.state_change(inp,'S','P',0,ty1='log',ty2='log'), # Pumping process
      lambda inp: BaseCycle.state_change(inp,'H','P',1,ty1='lin',ty2='lin'), # Heat addition
      lambda inp: BaseCycle.state_change(inp,'H','P',2,ty1='log',ty2='log'), # Expansion
      lambda inp: BaseCycle.state_change(inp,'H','P',3,ty1='lin',ty2='lin')  # Heat removal
      ]
    
    def __init__(self, fluid_ref='HEOS::Water', graph_type='TS', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BasePowerCycle` for details."""
        BasePowerCycle.__init__(self, fluid_ref, graph_type, **kwargs)
    
    def simple_solve(self, T0, p0, T2, p2, eta_exp, eta_pum, fluid=None, SI=True):
        """" 
        A simple Rankine cycle calculation
        
        Parameters
        ----------
        T0 : float
            The coldest point, before the pump
        p0 : float
            The lowest pressure, before the pump
        T2 : float
            The hottest point, before the expander
        p2 : float
            The highest pressure, before the expander
        eta_exp : float
            Isentropic expander efficiency 
        eta_pum : float
            Isentropic pump efficiency
        
        Examples
        --------
        >>> import CoolProp
        >>> from CoolProp.Plots import PropertyPlot
        >>> pp = PropertyPlot('HEOS::Water', 'TS', unit_system='EUR')
        >>> cycle = SimpleRankineCycle('HEOS::Water', 'TS', unit_system='EUR')
        >>> T0 = 300
        >>> pp.state.update(CoolProp.QT_INPUTS,0.0,T0+15)
        >>> p0 = pp.state.keyed_output(CoolProp.iP)
        >>> T2 = 700
        >>> pp.state.update(CoolProp.QT_INPUTS,1.0,T2-150)
        >>> p2 = pp.state.keyed_output(CoolProp.iP)
        >>> cycle.simple_solve(T0, p0, T2, p2, 0.7, 0.8, SI=True)
        >>> cycle.steps = 50
        >>> sc = cycle.get_state_changes()
        >>> pp.draw_process(sc)
        
        """
        if fluid is not None: self.state = process_fluid_state(fluid)
        if self._state is None: 
            raise ValueError("You have specify a fluid before you calculate.")
        
        cycle_states = StateContainer(unit_system=self._system)
        
        if not SI:
            Tc = self._system[CoolProp.iT].to_SI
            pc = self._system[CoolProp.iP].to_SI
            T0 = Tc(T0)
            p0 = pc(p0)
            T2 = Tc(T2)
            p2 = pc(p2)
        
        # Subcooled liquid
        self.state.update(CoolProp.PT_INPUTS,p0,T0)
        h0 = self.state.hmass()
        s0 = self.state.smass()
        # Just a showcase for the different accessor methods
        cycle_states[0,'H'] = h0
        cycle_states[0]['S'] = s0
        cycle_states[0][CoolProp.iP] = p0
        cycle_states[0,CoolProp.iT] = T0
        
        # Pressurised liquid
        p1 = p2
        self.state.update(CoolProp.PSmass_INPUTS,p1,s0)
        h1 = h0 + (self.state.hmass() - h0) / eta_pum
        self.state.update(CoolProp.HmassP_INPUTS,h1,p1)
        s1 = self.state.smass()
        T1 = self.state.T()
        cycle_states[1,'H'] = h1
        cycle_states[1,'S'] = s1
        cycle_states[1,'P'] = p1
        cycle_states[1,'T'] = T1
        
        # Evaporated vapour    
        self.state.update(CoolProp.PT_INPUTS,p2,T2)
        h2 = self.state.hmass()
        s2 = self.state.smass()
        cycle_states[2,'H'] = h2
        cycle_states[2,'S'] = s2
        cycle_states[2,'P'] = p2
        cycle_states[2,'T'] = T2
        
        # Expanded gas
        p3 = p0
        self.state.update(CoolProp.PSmass_INPUTS,p3,s2)
        h3 = h2 - eta_exp * (h2 - self.state.hmass())
        self.state.update(CoolProp.HmassP_INPUTS,h3,p3)
        s3 = self.state.smass()
        T3 = self.state.T()
        cycle_states[3,'H'] = h3
        cycle_states[3,'S'] = s3
        cycle_states[3,'P'] = p3
        cycle_states[3,'T'] = T3
    
        w_net = h2 - h3
        q_boiler = h2 - h1
        eta_th = w_net / q_boiler
        
        self.cycle_states = cycle_states
        self.fill_states()
        
        
    def eta_thermal(self):
        """Thermal efficiency
         
        The thermal efficiency for the specified process(es), :math:`\eta_{th} = \frac{\dot{W}_{exp} - \dot{W}_{pum}}{\dot{Q}_{in}}`.
         
        Returns
        -------
        float
        """    
        w_net = self.cycle_states[2].H - self.cycle_states[3].H - (self.cycle_states[1].H - self.cycle_states[0].H) 
        q_boiler = self.cycle_states[2].H - self.cycle_states[1].H
        return w_net / q_boiler
        

