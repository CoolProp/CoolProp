from __future__ import print_function, division
from six import with_metaclass

import matplotlib,numpy

import numpy as np

import CoolProp
from CoolProp.CoolProp import PropsSI
from scipy.optimize import newton
from .Common import BasePlot, _process_fluid_state, UnitSystem, SIunits

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
    T[2]=T_hp(Ref,h[2],pi,T2s)
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
    T[3]=T_hp(Ref,h[3],pi,T[2])
    s[3]=PropsSI('S','T',T[3],'P',pi,Ref)
    rho[3]=PropsSI('D','T',T[3],'P',pi,Ref)
    T4s=newton(lambda T: PropsSI('S','T',T,'P',pc,Ref)-s[3],T[2]+30)
    h4s=PropsSI('H','T',T4s,'P',pc,Ref)
    p[4]=pc
    wdot2=(h4s-h[3])/eta_oi
    h[4]=h[3]+(1-f_p[1])*wdot2
    T[4]=T_hp(Ref,h[4],pc,T4s)
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
#             via :func:`Common._process_fluid_state`.
#         steps : int
#             The number of steps used for going from one state to another
#         
#         for more properties, see :class:`CoolProp.Plots.Common.Base2DObject`.  

class StateContainer(object):
    """A collection of values for the main properties, built to mixin with :class:`CoolProp.Plots.Common.UnitSystem`"""
    
    def __init__(self,unit_system=SIunits()):
        self.values = {}
        self._system = unit_system
        
    @property
    def system(self): return self._system
    @system.setter
    def system(self, value): self._system = value
        
    @property
    def dimensions(self):
        return self._system.dimensions
    
    @property
    def d(self): return self._d
    @d.setter
    def d(self, value): self._d = value
    @property
    def h(self): return self._h
    @h.setter
    def h(self, value): self._h = value
    @property
    def p(self): return self._p
    @p.setter
    def p(self, value): self._p = value
    @property
    def s(self): return self._s
    @s.setter
    def s(self, value): self._s = value
    @property
    def t(self): return self._t
    @t.setter
    def t(self, value): self._t = value
    
    @property
    def values(self): 
        return {
      CoolProp.iDmass : self._d,
      CoolProp.iHmass : self._h,
      CoolProp.iP     : self._p,
      CoolProp.iSmass : self._s,
      CoolProp.iT     : self._t,
    }
    
    @values.setter
    def values(self,values): 
        self.d = None 
        self.h = None 
        self.p = None 
        self.s = None 
        self.t = None 
        for i in values:
            if i == CoolProp.iDmass : self.d = values[i] 
            if i == CoolProp.iHmass : self.h = values[i] 
            if i == CoolProp.iP     : self.p = values[i] 
            if i == CoolProp.iSmass : self.s = values[i] 
            if i == CoolProp.iT     : self.t = values[i] 
        
    def get_si_states(self):
        return self.values
    
    def get_local_states(self):
        states = {}
        for n in self.values:
            states[n] = self.dimensions[n].from_SI(self.values[n])
        return states
    
    def set_si_states(self,values):
        self.values = values
    
    def set_local_states(self,values):
        states = {}
        for n in values:
            states[n] = self.dimensions[n].to_SI(np.asarray(values[n]))
        self.values = states
        
    def __len__(self):
        """Some cheating to get the correct behaviour"""
        return np.min([len(np.asarray(c)) for c in self.values.values()])
    
    def __getitem__(self, index):
        """Another tweak that changes the default access path"""
        state = {}
        for n in self.values:
            state[n] = self.dimensions[n].from_SI(self.values[n][index])
        return state
    
    def __setitem__(self, index, values):
        """Another tweak that changes the default access path"""
        for n in values:
            if n == CoolProp.iDmass : self.d[index] = self.dimensions[n].to_SI(values[n])
            if n == CoolProp.iHmass : self.h[index] = self.dimensions[n].to_SI(values[n])
            if n == CoolProp.iP     : self.p[index] = self.dimensions[n].to_SI(values[n])
            if n == CoolProp.iSmass : self.s[index] = self.dimensions[n].to_SI(values[n])
            if n == CoolProp.iT     : self.t[index] = self.dimensions[n].to_SI(values[n])
            

class BaseCycle(BasePlot):
    """A simple thermodynamic cycle, should not be used on its own."""
    
    # Define the iteration keys
    PROPERTIES = {
      CoolProp.iDmass:'density', 
      CoolProp.iHmass:'specific enthalpy', 
      CoolProp.iP:'pressure', 
      CoolProp.iSmass:'specific entropy', 
      CoolProp.iT:'temperature'
    }
    
    STATECOUNTS=[0]
    """A list of accepted numbers of states"""
    
    STATECHANGE=[None]
    """A list of lists of tuples that defines how the state transitions 
    behave for the corresponding entry in BaseCycle.STATECOUNTS"""
    
    def __init__(self, fluid_ref, graph_type, unit_system='EUR', **kwargs):
        """Initialises a simple cycle calculator
        
        Parameters
        ----------        
        fluid_ref : str, CoolProp.AbstractState
            The fluid property provider, either a subclass of CoolProp.AbstractState
            or a string that can be used to generate a CoolProp.AbstractState instance
            via :func:`Common._process_fluid_state`.
        graph_type : string
            The graph type to be plotted, like \"PH\" or \"TS\"
        unit_system : string, ['EUR','KSI','SI']
            Select the units used for the plotting.  'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
                
        for more properties, see :class:`CoolProp.Plots.Common.BasePlot`.        
        """
        self._cycle_states = StateContainer(SIunits())
        BasePlot.__init__(self, fluid_ref, graph_type, unit_system, **kwargs)
        
    
    @property
    def cycle_states(self): return self._cycle_states
    @cycle_states.setter
    def cycle_states(self, value):
        if len(value) not in self.STATECOUNTS:
            raise ValueError("Your number of states ({0:d}) is not in the list of allowed state counts: {1:s}.".format(len(value),str(self.STATECOUNTS)))
        self._cycle_states = value
    
    @property
    def steps(self): return self._steps
    @steps.setter
    def steps(self, value): self._steps = int(max([value,1.0]))

    @BasePlot.system.setter
    def system(self, value): 
        if value in self.UNIT_SYSTEMS:
            self._system = self.UNIT_SYSTEMS[value]
        elif isinstance(value, UnitSystem):
            self._system = value
        else:
            raise ValueError("Invalid unit_system input \"{0:s}\", expected a string from {1:s}".format(str(value),str(self.UNIT_SYSTEMS.keys())))
        self._cycle_states.system = self._system       
    
    
    def valid_states(self):
        """Check the formats of BaseCycle.STATECOUNTS and BaseCycle.STATECHANGE"""
        for i,sn in enumerate(self.STATECOUNTS):
            if len(self.STATECHANGE[i]) != sn: 
                raise ValueError("Invalid number of states and or state change operations")
        return True 
    
    
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
        raise NotImplementedError()
    

class SimpleRankineCycle(BaseCycle):
    """A simple Rankine cycle *without* regeneration"""
    STATECOUNTS=[4]
    STATECHANGE=[[
      lambda: BaseCycle.state_change('S','P',0,ty1='lin',ty2='log'), # Pumping process
      lambda: BaseCycle.state_change('H','P',1,ty1='lin',ty2='lin'), # Heat addition
      lambda: BaseCycle.state_change('S','P',2,ty1='lin',ty2='log'), # Expansion
      lambda: BaseCycle.state_change('H','P',3,ty1='lin',ty2='lin')  # Heat removal
      ]]
    
    def __init__(self, fluid_ref='HEOS::Water', graph_type='TS', **kwargs):
        """see :class:`CoolProp.Plots.SimpleCycles.BaseCycle` for details.      
        """
        BaseCycle.__init__(self, fluid_ref, graph_type, **kwargs)
        
    
    
    def simple_solve(self, T0, p0, T2, p2, eta_exp, eta_pum, fluid=None):
        if fluid is not None: self.state = _process_fluid_state(fluid)
        if self._state is None: 
            raise ValueError("You have specify a fluid before you calculate.")
        
        self.state.update(CoolProp.PT_INPUTS,p0,T0) 
        h0 = self.state.hmass()
        s0 = self.state.smass()
    
        p1 = p2
        self.state.update(CoolProp.PSmass_INPUTS,p1,s0)
        h1 = h0 + (self.state.hmass() - h0) / eta_pum
        self.state.update(CoolProp.HmassP_INPUTS,h1,p1)
        s1 = self.state.smass()
        T1 = self.state.T()
    
        self.state.update(CoolProp.PT_INPUTS,p2,T2)
        h2 = self.state.hmass()
        s2 = self.state.smass()
    
        p3 = p2
        self.state.update(CoolProp.PSmass_INPUTS,p3,s2)
        h3 = h2 - eta_exp * (h2 - self.state.hmass())
        self.state.update(CoolProp.HmassP_INPUTS,h3,p3)
        s3 = self.state.smass()
        T3 = self.state.T()
    
        w_net = h2 - h3
        q_boiler = h2 - h1
        eta_th = w_net / q_boiler
        
        print(eta_th)
    
        



# def SimpleRankineCycle(states, steps, fluid):
#     """A simple Rankine cycle, with optional regeneration
#     
#     Calculates a simple Rankine cycle including some basic cycle-related quantities. 
#     
#     Parameters
#     ----------
#         states : list, dict 
#             A collection of states that follows a fixed scheme:
#                 0) Subcooled liquid
#                 1) Pressurised liquid
#                 2) Preheated liquid    if regenrative cycle, else superheated state
#                 3) Superheated state   if regenrative cycle, else expanded state
#                 4) Expanded state      if regenrative cycle, else NULL
#                 5) Desuperheated state if regenrative cycle, else NULL  
#     
#     
#     """
#     state = _process_fluid_state(fluid)
# 
# 
#     state.update(CoolProp.PT_INPUTS,p1,T1)
#     h1 = state.hmass()
#     s1 = state.smass()
# 
#     p2 = p3
#     state.update(CoolProp.PSmass_INPUTS,p2,s1)
#     h2 = h1 + (state.hmass() - h1) / epsilon_p
#     state.update(CoolProp.HmassP_INPUTS,h2,p2)
#     s2 = state.smass()
#     T2 = state.T()
# 
#     state.update(CoolProp.PT_INPUTS,p3,T3)
#     h3 = state.hmass()
#     s3 = state.smass()
# 
#     p4 = p1
#     state.update(CoolProp.PSmass_INPUTS,p4,s3)
#     h4 = h3 - epsilon_e * (h3 - state.hmass())
#     state.update(CoolProp.HmassP_INPUTS,h4,p4)
#     s4 = state.smass()
#     T4 = state.T()
# 
#     w_net = h3 - h4
#     q_boiler = h3 - h2
#     eta_c = w_net / q_boiler
# 
#     #Ts = PropsPlot(fluid, 'Ts')
#     #Ts.draw_isolines('P', [p1, p3], num=10)
#     #Ts.set_axis_limits([0., 12., 200., 900.])
# 
#     #axObj.plot(s_tp/1e3,T_tp-273.15 , color=plotterObj._black, ls='-', alpha=1.0)
# 
#     isoObj  = IsoLines(fluid, "Ts", "Q")
#     isoqual = isoObj.get_isolines([0.0,1.0], num=2)
# 
#     x = np.append(isoqual[ 0]['x'],isoqual[-1]['x'][::-1])/1e3
#     y = np.append(isoqual[ 0]['y'],isoqual[-1]['y'][::-1])-273.15
#     axObj.plot(x,y, color=plotterObj._black, ls='-', alpha=1.0)
# 
#     isoObj  = IsoLines(fluid, "Ts", "P")
#     prange  = [p1,2e5,5e5,10e5,p3]
#     isobars = isoObj.get_isolines(prange, num=len(prange))
# 
#     p = -1
#     for c,i in enumerate(isobars):
#         x = i['x']/1e3
#         y = i['y']-273.15
#         dp = prange[c]/1e5 - p
#         p = prange[c]/1e5
#         s = PropsSI('S','P',p*1e5,'Q',0.5,fluid)/1e3
#         #print "Delta p: {0}".format(dp)
#         if abs(dp)>0.8: #c%2==0 :
#             axObj.plot(  x,    y, color=plotterObj._black, ls='-', alpha=0.50)
#             if label:
#                 putXLabel(xv=x, yv=y, x=s, text="{0:3.1f} bar".format(p), axis=axObj)
# 
#     #for i in range(len(x)):
#     #    axObj.plot(  x[i]/1e3,   y[i]-273.15, color=plotterObj._black, ls='-', alpha=0.5)
#     #    putXLabel(xv=x[i]/1e3,yv=y[i]-273.15, x=0, text="", axis=axObj)
# 
# 
# 
#     # Create the process lines
#     A = []
#     A.append({'H':h1,'P':p1,'S':s1,'T':T1})
#     A.append({'H':h2,'P':p2,'S':s2,'T':T2})
#     A.append({'H':h3,'P':p3,'S':s3,'T':T3})
#     A.append({'H':h4,'P':p4,'S':s4,'T':T4})
# 
#     A.append(A[0].copy())
# 
#     processes = []
# 
#     for i in range(len(A)-1):
#         s = np.linspace(         A[i]['S'] ,          A[i+1]['S'] ,num=points)
#         p = np.logspace(np.log10(A[i]['P']), np.log10(A[i+1]['P']),num=points)
#         dic = {}
#         dic['P'] = p
#         dic['S'] = s
#         dic['T'] = PropsSI('T','P',p,'S',s,fluid)
#         processes.append(dic)
# 
#     x = []
#     y = []
#     for lin in processes:
#         #axObj.plot(lin['S']/1e3,lin['T']-273.15,color=plotterObj._black, linestyle='--')
#         x.extend(lin['S']/1e3)
#         y.extend(lin['T']-273.15)
# 
#     plotterObj.plotData([x],[y],ax=axObj,legend=False)
# 
#     x = np.array([s1,s2,s3,s4])
#     y = np.array([T1,T2,T3,T4])
# 
#     #print x
#     #print y
#     #print " "
# 
#     plotterObj.plotData([x/1e3],[y-273.15],ax=axObj,legend=False)
# 
#     #axObj.plot(x/1e3,y-273.15,'o',color=plotterObj._black)
# 
#     #plotterObj.drawLegend(ax=axObj,loc=0) # the relative size of legend markers vs. original
#     axObj.set_xlabel(ur"Specific entropy $s$ / \si{\kilo\joule\per\kilo\gram\per\kelvin}")
#     axObj.set_ylabel(ur"Temperature $T$ / \si{\celsius}")
#     axObj.set_xlim([-0.25,1.60])
#     axObj.set_ylim([-25,325])
# 
#     #plotterObj.plotData([x], [y], ax=axObj)
# 
#     #ax = Ts.axis
#     #ax.text(s1/1000., T1,' 1', fontsize=10, rotation=0, color='r')
#     #ax.text(s2/1000., T2,' 2', fontsize=10, rotation=0, color='r')
#     #ax.text(s3/1000., T3,' 3', fontsize=10, rotation=0, color='r')
#     #ax.text(s4/1000., T4,' 4', fontsize=10, rotation=0, color='r')
#     #ax.text(8., 850., "Efficiency: %.1f%%" %(eta_c*100.))
#     #ax.text(8., 800., "Net work: %d kJ/kg" %(w_net/1000))
#     #ax.text(8., 750., "Heat input: %d kJ/kg" %(q_boiler/1000))
# 
# simPlotterObj = BasePlotter()
# figPV = simPlotterObj.getFigure(**sixupProps)
# simPlotterObj.ccycle = simPlotterObj.multiplyCycle(simPlotterObj.getColorCycle(length=3),doubles=2)
# simPlotterObj.scycle = cycle(['-'])
# simPlotterObj.mcycle = cycle(['None'])





# if __name__=='__main__':
#     
#     cycle = SimpleRankineCycle()
#     
#     cycle.simple_solve(300, 2e5, 550, 10e5, 0.7, 0.8)
    
    
#     from CoolProp.Plots import Ph,Ts
# 
#     Ref='R290'
#     fig=matplotlib.pyplot.figure(figsize=(4,3))
#     ax=fig.add_axes((0.15,0.15,0.8,0.8))
#     Ph(Ref,Tmin=273.15-30,hbounds=[0,600],axis=ax)
#     COP=TwoStage('Propane',10000,273.15-5,273.15+43.3,5,7,0.7,0.3,15+273.15,3,prints = True)
#     matplotlib.pyplot.show()
# 
#     Ref='R290'
#     fig=matplotlib.pyplot.figure(figsize=(4,3))
#     ax=fig.add_axes((0.15,0.15,0.8,0.8))
#     Ph(Ref,Tmin=273.15-30,hbounds=[0,600],axis=ax)
#     COP=SimpleCycle(Ref,273.15-5,273.15+45,5,7,0.7,Ts_Ph='Ph')
#     matplotlib.pyplot.show()
# 
#     Ref='R410A'
#     fig=matplotlib.pyplot.figure(figsize=(4,3))
#     ax=fig.add_axes((0.15,0.15,0.8,0.8))
#     Ts(Ref,Tmin=273.15-100,sbounds=[0,600],axis=ax)
#     COP=SimpleCycle(Ref,273.15-5,273.15+45,5,7,0.7,Ts_Ph='Ts')
#     matplotlib.pyplot.show()




##     for x in numpy.linspace(0,1):
##         Ref='REFPROP-MIX:R152A[%g]&R32[%g]' %(x,1-x)
##         COP=SimpleCycle(273.15+8,273.15+44,5,7,0.7,skipPlot=True,Ts_Ph='Ph')
##     matplotlib.pyplot.show()
