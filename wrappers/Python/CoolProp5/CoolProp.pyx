#cython: embedsignature = True, c_string_type=str, c_string_encoding=ascii
from __future__ import division
#
# This file provides wrapper functions of all the CoolProp functions
#
# Each of the functions from the CoolProp header are renamed in cython code to
# an underscored name so that the same name can be used in the exposed functions below

import cython
cimport cython

import math
import warnings

from libcpp.string cimport string
from libcpp.vector cimport vector

from constants import *
from constants_header cimport *

cpdef bint iterable(object a):
    if _numpy_supported:
        return isinstance(a,(list,tuple, np.ndarray))
    else:
        return isinstance(a,(list,tuple))

cpdef ndarray_or_iterable(object input):
    if _numpy_supported:
        return np.array(input)
    else:
        return input

include "HumidAirProp.pyx"
include "AbstractState.pyx"
    
# def set_reference_state(string_like FluidName, *args):
#     """
#     Accepts one of two signatures:
#     
#     Type #1:
#     
#     set_reference_state(FluidName,reference_state)
#     
#     FluidName The name of the fluid
#     param reference_state The reference state to use, one of 
#     
#     ==========   ===========================================
#     ``IIR``      (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.)
#     ``ASHRAE``   (h=0,s=0 @ -40C sat liq)
#     ``NBP``      (h=0,s=0 @ 1.0 bar sat liq.)
#     ==========   ===========================================
#     
#     Type #2:
#     
#     set_reference_state(FluidName,T0,rho0,h0,s0)
#     
#     ``FluidName`` The name of the fluid
#     
#     ``T0`` The temperature at the reference point [K]
#     
#     ``rho0`` The density at the reference point [kg/m^3]
#     
#     ``h0`` The enthalpy at the reference point [J/kg]
#     
#     ``s0`` The entropy at the reference point [J/kg]
#     """
#     
#     cdef bytes _param
#     cdef int retval
#     
#     if len(args) == 1:
#         _param = args[0].encode('ascii')
#         retval = _set_reference_stateS(FluidName, _param)
#     elif len(args) == 4:
#         retval = _set_reference_stateD(FluidName, args[0], args[1], args[2], args[3])
#     else:
#         raise ValueError('Invalid number of inputs')
#     
#     if retval < 0:
#         raise ValueError('Unable to set reference state')
#         
#     
# cpdef long get_Fluid_index(string_like Fluid):
#     """
#     Gets the integer index of the given CoolProp fluid (primarily for use in ``IProps`` function)
#     """
#     return _get_Fluid_index(Fluid)
#     
# cpdef double IProps(long iOutput, long iInput1, double Input1, long iInput2, double Input2, long iFluid) except *:
#     """
#     This is a more computationally efficient version of the Props() function as it uses integer keys for the input and output codes as well as the fluid index for the fluid.  It can only be used with CoolProp fluids.  An example of how it should be used::
#     
#         #  These should be run once in the header of your file
#         from CoolProp.CoolProp import IProps, get_Fluid_index
#         from CoolProp import param_constants
#         iPropane = get_Fluid_index('Propane')
#         
#         #  This should be run using the cached values - much faster !
#         IProps(param_constants.iP,param_constants.iT,0.8*Tc,param_constants.iQ,1,iPropane)
#         
#     The reason that this function is significantly faster than Props is that it skips all the string comparisons which slows down the Props function quite a lot.  At the C++ level, IProps doesn't use any strings and operates on integers and floating point values
#     """
#     cdef double val = _IProps(iOutput, iInput1, Input1, iInput2, Input2, iFluid)
#     
#     if math.isinf(val) or math.isnan(val):
#         err_string = _get_global_param_string('errstring')
#         if not len(err_string) == 0:
#             raise ValueError("{err:s} :: inputs were :{iin1:d},{in1:g},{iin2:d},{in2:g},{iFluid:d}".format(err= err_string,iin1=iInput1,in1=Input1,iin2=iInput2,in2=Input2,iFluid = iFluid))
#         else:
#             raise ValueError("IProps failed ungracefully with inputs:\"{in1:g}\",\"{in2:g}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(in1=Input1,in2=Input2))
#     else:
#         return val

cpdef get_global_param_string(string param):
    return _get_global_param_string(param)
     
cpdef get_fluid_param_string(string_like fluid, string_like param):
    return _get_fluid_param_string(fluid, param)
#     
# cpdef __Props_err1(in1,in2,errstr):
#     if not len(errstr) == 0:
#         raise ValueError("{err:s} :: inputs were :\"{in1:s}\",\"{in2:s}\"".format(err= errstr,in1=in1,in2=in2))
#     else:
#         raise ValueError("Props failed ungracefully with inputs:\"{in1:s}\",\"{in2:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(in1=in1,in2=in2))
#         
cpdef __Props_err2(in1, in2, in3, in4, in5, in6):
    errstr = _get_global_param_string('errstring')
    if not len(errstr) == 0:
        raise ValueError("{err:s} :: inputs were:\"{in1:s}\",\"{in2:s}\",{in3:0.16e},\"{in4:s}\",{in5:0.16e},\"{in6:s}\"".format(err=errstr,in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))
    else:
        raise ValueError("Props failed ungracefully :: inputs were:\"{in1:s}\",\"{in2:s}\",{in3:0.16e},\"{in4:s}\",{in5:0.16e},\"{in6:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))

# cpdef Props(in1, in2, in3, in4, in5, in6, in7 = None):
#     """
#     $$Props$$
#     """ 
#     if in7 is None:
#         return _Props(in1, in2, in3, in4, in5, in6) 
#     else:
#         return _Props(in1, in2, in3, in4, in5, in6, in7)
cpdef PropsSI(in1, in2, in3, in4, in5, in6, in7 = None):
    """
    $$PropsSI$$
    """ 
    cdef double val
    if in7 is None:
        val = _PropsSI(in1, in2, in3, in4, in5, in6)
        if not _ValidNumber(val): 
            __Props_err2(in1, in2, in3, in4, in5, in6)
        else:
            return val
    else:
        return _PropsSI(in1, in2, in3, in4, in5, in6, in7)

cpdef list FluidsList():
    """
    Return a list of strings of all fluid names
    
    Returns
    -------
    FluidsList : list of strings of fluid names
        All the fluids that are included in CoolProp
    
    Notes
    -----
    
    Here is an example::
        
       In [0]: from CoolProp.CoolProp import FluidsList
    
       In [1]: FluidsList()
       
    """ 
    return _get_global_param_string("FluidsList").split(',')

cpdef get_aliases(str Fluid):
    """
    Return a comma separated string of aliases for the given fluid
    """
    cdef bytes _Fluid = Fluid.encode('ascii')
    return [F.encode('ascii') for F in (_get_fluid_param_string(_Fluid,'aliases').encode('ascii')).decode('ascii').split(',')]

cpdef string get_REFPROPname(str Fluid):
    """
    Return the REFPROP compatible name for the fluid
    
    Some fluids do not use the REFPROP name.  For instance, 
    ammonia is R717, and propane is R290.  You can still can still call CoolProp
    using the name ammonia or R717, but REFPROP requires that you use a limited
    subset of names.  Therefore, this function that returns the REFPROP compatible
    name.  To then use this to call REFPROP, you would do something like::
    
       In [0]: from CoolProp.CoolProp import get_REFPROPname, PropsSI
    
       In [1]: get_REFPROPname('R290')
       
       In [2]: PropsSI('D', 'T', 300, 'P', 300, Fluid)
    """
    cdef bytes _Fluid = Fluid.encode('ascii')
    return _get_fluid_param_string(_Fluid,'REFPROP_name')

# cpdef string get_BibTeXKey(str Fluid, str key):
#     """
#     Return the BibTeX key for the given fluid.
#     
#     The possible keys are
#     
#     * ``EOS``
#     * ``CP0``
#     * ``VISCOSITY``
#     * ``CONDUCTIVITY``
#     * ``ECS_LENNARD_JONES``
#     * ``ECS_FITS``
#     * ``SURFACE_TENSION``
#     
#     BibTeX keys refer to the BibTeX file in the trunk/CoolProp folder
#     
#     Returns
#     -------
#     key, string
#          empty string if Fluid not in CoolProp, "Bad key" if key is invalid
#     """
#     return _get_BibTeXKey(Fluid, key)
cpdef string get_errstr():
    """
    Return the current error string
    """
    return _get_global_param_string("errstring")
    
cpdef set_debug_level(int level):
    """
    Set the current debug level as integer in the range [0,10]
    
    Parameters
    ----------
    level : int
        If level is 0, no output will be written to screen, if >0, 
        some output will be written to screen.  The larger level is, 
        the more verbose the output will be
    """
    _set_debug_level(level)

cpdef int get_debug_level():
    """
    Return the current debug level as integer
    
    Returns
    -------
    level : int
        If level is 0, no output will be written to screen, if >0, 
        some output will be written to screen.  The larger level is, 
        the more verbose the output will be
    """
    return _get_debug_level()

# cpdef bint IsFluidType(string Fluid, string Type):
#     """
#     Check if a fluid is of a given type
#     
#     Valid types are:
#     
#     * ``Brine``
#     * ``PseudoPure`` (or equivalently ``PseudoPureFluid``)
#     * ``PureFluid``
#     """
#     cdef bytes _Fluid = Fluid.encode('ascii')
#     cdef bytes _Type = Type.encode('ascii')
#     if _IsFluidType(_Fluid,_Type):
#         return True
#     else:
#         return False
#     

# from math import pow as pow_

#A dictionary mapping parameter index to string for use with non-CoolProp fluids
cdef dict paras = {iDmass : 'D',
                   iQ : 'Q',
                   imolar_mass : 'M',
                   iT : 'T',
                   iHmass : 'H',
                   iP : 'P',
                   iCpmass : 'C',
#                    iC0 : 'C0',
                   iCvmass : 'O',
                   iviscosity : 'V',
                   iconductivity : 'L',
                   iSmass : 'S',
                   iUmass : 'U',
#                    iDpdT : 'dpdT'
}

# cdef dict paras_inverse = {v:k for k,v in paras.iteritems()}

cdef class State:
    """
    A class that contains all the code that represents a thermodynamic state
    
    The motivation for this class is that it is useful to be able to define the
    state once using whatever state inputs you like and then be able to calculate
    other thermodynamic properties with the minimum of computational work.
    
    Let's suppose that you have inputs of pressure and temperature and you want
    to calculate the enthalpy and pressure.  Since the Equations of State are
    all explicit in temperature and density, each time you call something like::
    
        h = PropsSI('H','T',T','P',P,Fluid)
        s = PropsSI('S','T',T','P',P,Fluid)
        
    the solver is used to carry out the T-P flash calculation. And if you wanted
    entropy as well you could either intermediately calculate ``T``, ``rho`` and then use
    ``T``, ``rho`` in the EOS in a manner like::
    
        rho = PropsSI('D','T',T','P',P,Fluid)
        h = PropsSI('H','T',T','D',rho,Fluid)
        s = PropsSI('S','T',T','D',rho,Fluid)
        
    Instead in this class all that is handled internally. So the call to update
    sets the internal variables in the most computationally efficient way possible
    """
        
    def __init__(self, object Fluid, dict StateDict, object phase = None, backend = None):
        """
        Parameters
        ----------
        Fluid : string
        StateDict : dictionary
            The state of the fluid - passed to the update function
        phase : string
            DEPRECATED : this input is ignored
        backend : string
            The CoolProp backend that should be used, one of "HEOS" (default), "REFPROP", "INCOMP", "BRINE", etc.
        """
        cdef string _Fluid = Fluid
        
        if _Fluid == <string>'none':
            return
        else:
            if backend is None:
                backend = '?'
            self.set_Fluid(Fluid, backend)
#         
#         self.phase = _phase
#         #Parse the inputs provided
#         self.update(StateDict)
#         #Set the phase flag
#         if self.phase == <string>'Gas' or self.phase == <string>'Liquid' or self.phase == <string>'Supercritical':
#             if self.is_CPFluid and (self.phase == <string>'Gas' or self.phase == <string>'Liquid'):
#                 self.CPS.flag_SinglePhase = True
#             elif not self.is_CPFluid and phase is not None:
#                 _set_phase(self.phase)
            
#     def __reduce__(self):
#         d={}
#         d['Fluid']=self.Fluid
#         d['T']=self.T_
#         d['rho']=self.rho_
#         d['phase'] = self.phase
#         return rebuildState,(d,)
        
    cpdef set_Fluid(self, string Fluid, string backend):
        self.pAS = AbstractState(backend, Fluid)
        
#     cpdef update_ph(self, double p, double h):
#         """
#         Use the pressure and enthalpy directly
#         
#         Parameters
#         ----------
#         p: float
#             Pressure (absolute) [kPa]
#         h: float
#             Enthalpy [kJ/kg]
#         
#         """
#         p = _toSIints(iP, p, _get_standard_unit_system());
#         h = _toSIints(iH, h, _get_standard_unit_system());
#         self.p_ = p
#         cdef double T
#         
#         if self.is_CPFluid:
#             self.CPS.update(iP, p, iH, h)
#             self.T_ = self.CPS.T()
#             self.rho_ = self.CPS.rho()
#         else:
#             T = _PropsSI('T','P',p,'H',h,self.Fluid)
#             if abs(T)<1e90:
#                 self.T_=T
#             else:
#                 errstr = _get_global_param_string('errstring')
#                 raise ValueError(errstr)
#             self.rho_ = _Props('D','P',p,'H',h,self.Fluid)
            
    cpdef update_Trho(self, double T, double rho):
        """
        Just use the temperature and density directly for speed
        
        Parameters
        ----------
        T: float
            Temperature [K]
        rho: float
            Density [kg/m^3]
            
        """
        self.T_ = T
        self.rho_ = rho
        
        self.pAS.update(DmassT_INPUTS, rho, T)
        
#     cpdef update(self, dict params):
#         """
#         Parameters
#         params, dictionary 
#             A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
#             for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
#         """
#         
#         # Convert to integer_pair input
#             
#         cdef double p, val1, val2
#         cdef long iInput1, iInput2
#         cdef bytes errstr
#             
#         if self.is_CPFluid:
#             items = list(params.items())
#             iInput1 = paras_inverse[items[0][0]]
#             iInput2 = paras_inverse[items[1][0]]
#             # Convert to SI units
#             val1 = _toSIints(iInput1, items[0][1], _get_standard_unit_system());
#             val2 = _toSIints(iInput2, items[1][1], _get_standard_unit_system());
#             try: 
#                 self.CPS.update(iInput1, val1, iInput2, val2)
#             except:
#                 raise
#             self.T_ = self.CPS.T()
#             self.p_ =  _fromSIints(iP, self.CPS.p(), _get_standard_unit_system());
#             self.rho_ = self.CPS.rho()
#             
#             if not _ValidNumber(self.T_) or not _ValidNumber(self.p_) or not _ValidNumber(self.rho_):
#                 raise ValueError(str(params))
#             return
#         
#         #Get the density if T,P provided, or pressure if T,rho provided
#         if 'P' in params:
#             self.p_=params['P']
#             rho = _Props('D','T',self.T_,'P',self.p_,self.Fluid)
#             
#             if abs(rho) < 1e90:
#                 self.rho_=rho
#             else:
#                 errstr = _get_global_param_string('errstring')
#                 raise ValueError(errstr)
#         elif 'D' in params:
#             self.rho_=params['D']
#             p = _Props('P','T',self.T_,'D',self.rho_,self.Fluid)
#             
#             if abs(p)<1e90:
#                 self.p_=p
#             else:
#                 errstr = _get_global_param_string('errstring')
#                 raise ValueError(errstr+str(params))
#         elif 'Q' in params:
#             p = _Props('P','T',self.T_,'Q',params['Q'],self.Fluid)
#             self.rho_ = _Props('D','T',self.T_,'Q',params['Q'],self.Fluid)
#             
#             if abs(self.rho_)<1e90:
#                 pass
#             else:
#                 errstr = _get_global_param_string('errstring')
#                 raise ValueError(errstr+str(params))
#         else:
#             raise KeyError("Dictionary must contain the key 'T' and one of 'P' or 'D'")
        
#     cpdef long Phase(self) except *:
#         """
#         Returns an integer flag for the phase of the fluid, where the flag value
#         is one of iLiquid, iSupercritical, iGas, iTwoPhase
#         
#         These constants are defined in the phase_constants module, and are imported
#         into this module
#         """
#         
#         if self.is_CPFluid:
#             return self.CPS.phase()
#         else:
#             raise NotImplementedError("Phase not defined for fluids other than CoolProp fluids")
        
#     cpdef double Props(self, long iOutput) except *: 
#         if iOutput<0:
#             raise ValueError('Your output is invalid') 
#         
#         if self.is_CPFluid:
#             val = self.CPS.keyed_output(iOutput)
#             return _fromSIints(iOutput,val,_get_standard_unit_system());
#         else:
#             return _Props(paras[iOutput],'T',self.T_,'D',self.rho_,self.Fluid)
#             
#     cpdef double get_Q(self) except *:
#         """ Get the quality [-] """
#         return self.Props(iQ)
#     property Q:
#         """ The quality [-] """
#         def __get__(self):
#             return self.get_Q()
#     
#     cpdef double get_MM(self) except *:
#         """ Get the mole mass [kg/kmol] or [g/mol] """
#         return self.Props(iMM)
#     
#     cpdef double get_rho(self) except *:
#         """ Get the density [kg/m^3] """ 
#         return self.Props(iD)
#     property rho:
#         """ The density [kg/m^3] """
#         def __get__(self):
#             return self.Props(iD)
#             
#     cpdef double get_p(self) except *:
#         """ Get the pressure [kPa] """ 
#         return self.Props(iP)
#     property p:
#         """ The pressure [kPa] """
#         def __get__(self):
#             return self.get_p()
#     
#     cpdef double get_T(self) except *: 
#         """ Get the temperature [K] """
#         return self.Props(iT)
#     property T:
#         """ The temperature [K] """
#         def __get__(self):
#             return self.get_T()
#     
#     cpdef double get_h(self) except *: 
#         """ Get the specific enthalpy [kJ/kg] """
#         return self.Props(iH)
#     property h:
#         """ The specific enthalpy [kJ/kg] """
#         def __get__(self):
#             return self.get_h()
#           
#     cpdef double get_u(self) except *: 
#         """ Get the specific internal energy [kJ/kg] """
#         return self.Props(iU)
#     property u:
#         """ The internal energy [kJ/kg] """
#         def __get__(self):
#             return self.get_u()
#             
#     cpdef double get_s(self) except *: 
#         """ Get the specific enthalpy [kJ/kg/K] """
#         return self.Props(iS)
#     property s:
#         """ The specific enthalpy [kJ/kg/K] """
#         def __get__(self):
#             return self.get_s()
#     
#     cpdef double get_cp0(self) except *:
#         """ Get the specific heat at constant pressure for the ideal gas [kJ/kg/K] """
#         return self.Props(iC0)
#     
#     cpdef double get_cp(self) except *: 
#         """ Get the specific heat at constant pressure  [kJ/kg/K] """
#         return self.Props(iC)
#     property cp:
#         """ The specific heat at constant pressure  [kJ/kg/K] """
#         def __get__(self):
#             return self.get_cp()
#             
#     cpdef double get_cv(self) except *: 
#         """ Get the specific heat at constant volume  [kJ/kg/K] """
#         return self.Props(iO)
#     property cv:
#         """ The specific heat at constant volume  [kJ/kg/K] """
#         def __get__(self):
#             return self.get_cv()
#         
#     cpdef double get_speed_sound(self) except *: 
#         """ Get the speed of sound  [m/s] """
#         return self.Props(iA)
#             
#     cpdef double get_visc(self) except *:
#         """ Get the viscosity, in [Pa-s]"""
#         return self.Props(iV)
#     property visc:
#         """ The viscosity, in [Pa-s]"""
#         def __get__(self):
#             return self.get_visc()

#     cpdef double get_cond(self) except *:
#         """ Get the thermal conductivity, in [kW/m/K]"""
#         return self.Props(iL)
#     property k:
#         """ The thermal conductivity, in [kW/m/K]"""
#         def __get__(self):
#             return self.get_cond()
#         
#     cpdef get_Tsat(self, double Q = 1):
#         """ 
#         Get the saturation temperature, in [K]
#         
#         Returns ``None`` if pressure is not within the two-phase pressure range 
#         """
#         if self.p_ > _Props1(self.Fluid,'pcrit') or self.p_ < _Props1(self.Fluid,'ptriple'):
#             return None 
#         else:
#             return _Props('T', 'P', self.p_, 'Q', Q, self.Fluid)
#     property Tsat:
#         """ The saturation temperature (dew) for the given pressure, in [K]"""
#         def __get__(self):
#             return self.get_Tsat(1.0)
#         
#     cpdef get_superheat(self):
#         """ 
#         Get the amount of superheat above the saturation temperature corresponding to the pressure, in [K]
#         
#         Returns ``None`` if pressure is not within the two-phase pressure range 
#         """
#         
#         Tsat = self.get_Tsat(1) #dewpoint temp
#         
#         if Tsat is not None:
#             return self.T_-Tsat
#         else:
#             return None
#     property superheat:
#         """ 
#         The amount of superheat above the saturation temperature corresponding to the pressure, in [K]
#         
#         Returns ``None`` if pressure is not within the two-phase pressure range 
#         """
#         def __get__(self):    
#             return self.get_superheat()
#         
#     cpdef get_subcooling(self):
#         """ 
#         Get the amount of subcooling below the saturation temperature corresponding to the pressure, in [K]
#         
#         Returns ``None`` if pressure is not within the two-phase pressure range 
#         """
#         
#         Tsat = self.get_Tsat(0) #bubblepoint temp
#         
#         if Tsat is not None:
#             return Tsat - self.T_
#         else:
#             return None
#     property subcooling:
#         """ 
#         The amount of subcooling below the saturation temperature corresponding to the pressure, in [K]
#         
#         Returns ``None`` if pressure is not within the two-phase pressure range 
#         """
#         def __get__(self):    
#             return self.get_subcooling()
#             
#     property Prandtl:
#         """ The Prandtl number (cp*mu/k) [-] """
#         def __get__(self):
#             return self.cp * self.visc / self.k
#             
#     cpdef double get_dpdT(self) except *:
#         if self.is_CPFluid:
#             return _fromSIints(iDERdp_dT__rho, self.CPS.dpdT_constrho(), _get_standard_unit_system());
#         else:
#             raise ValueError("get_dpdT not supported for fluids that are not in CoolProp")
#     property dpdT:
#         def __get__(self):
#             return self.get_dpdT()
#         
    cpdef speed_test(self, int N):
        from time import clock
#         cdef int i
#         cdef char * k
#         cdef long ikey
#         cdef string Fluid = self.Fluid
#         cdef long IT = 'T'
#         cdef long ID = 'D'
#         import CoolProp as CP
#         
#         print 'Call to the Python call layer (CoolProp.CoolProp.Props)'
#         print "'M' involves basically no computational effort and is a good measure of the function call overhead"
#         keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
#         for key in keys:
#             t1=clock()
#             for i in range(N):
#                 CP.Props(key,'T',self.T_,'D',self.rho_,Fluid)
#             t2=clock()
#             print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
#             
#         print 'Direct c++ call to CoolProp without the Python call layer (_Props function)'
#         print "'M' involves basically no computational effort and is a good measure of the function call overhead"
#         keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
#         for key in keys:
#             t1=clock()
#             for i in range(N):
#                 _Props(key,'T',self.T_,'D',self.rho_,Fluid)
#             t2=clock()
#             print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
#         
#         print 'Call to the c++ layer through IProps'
#         keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
#         for key in keys:
#             t1=clock()
#             for i in range(N):
#                 _IProps(key,iT,self.T_,iD,self.rho_,self.iFluid)
#             t2=clock()
#             print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
#             
        print 'Call to the c++ layer using integers'
        keys = [iHmass, iP,iSmass,iUmass]
        for key in keys:
            t1=clock()
            for i in range(N):
                self.pAS.update(DmassT_INPUTS,self.rho_,self.T_)
                self.pAS.keyed_output(key)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
#         
#         keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
#         isenabled = _isenabled_TTSE_LUT(<bytes>Fluid)
#         _enable_TTSE_LUT(<bytes>Fluid)
#         _IProps(iH,iT,self.T_,iD,self.rho_,self.iFluid)
#         
#         print 'Call using TTSE with T,rho'
#         print "'M' involves basically no computational effort and is a good measure of the function call overhead"
#         for ikey in keys:
#             t1=clock()
#             self.CPS.update(iT,self.T_,iD,self.rho_)
#             for i in range(N):
#                 self.CPS.keyed_output(ikey)
#             t2=clock()
#             print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[ikey],(t2-t1)/N*1e6)
#             
#         print 'Call using TTSE with p,h'
#         print "'M' involves basically no computational effort and is a good measure of the function call overhead"
#         cdef double hh = self.h
#         for ikey in keys:
#             t1=clock()
#             self.CPS.update(iP,self.p_,iH,hh)
#             for i in range(N):
#                 self.CPS.keyed_output(ikey)
#             t2=clock()
#             print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[ikey],(t2-t1)/N*1e6)
#         
#         print 'Using CoolPropStateClass with T,rho with LUT'
#         keys = [iH,iP,iC,iO,iDpdT]
#         t1=clock()
#         for i in range(N):
#             self.CPS.update(iT,self.T_,iD,self.rho_)
#             for ikey in keys:
#                 self.CPS.keyed_output(ikey)
#         t2=clock()
#         print 'Elapsed time for {0:d} calls of iH,iP,iC,iO,iDpdT takes {1:g} us/call'.format(N,(t2-t1)/N*1e6)
#         
#         if not isenabled:
#             _disable_TTSE_LUT(<bytes>Fluid)
#     
#     def __str__(self):
#         """
#         Return a string representation of the state
#         """
#         units={'T': 'K', 
#                'p': 'kPa', 
#                'rho': 'kg/m^3',
#                'Q':'kg/kg',
#                'h':'kJ/kg',
#                'u':'kJ/kg',
#                's':'kJ/kg/K',
#                'visc':'Pa-s',
#                'k':'kW/m/K',
#                'cp':'kJ/kg/K',
#                'cv':'kJ/kg/K',
#                'dpdT':'kPa/K',
#                'Tsat':'K',
#                'superheat':'K',
#                'subcooling':'K',
#         }
#         s='phase = '+self.phase+'\n'
#         for k in ['T','p','rho','Q','h','u','s','visc','k','cp','cv','dpdT','Prandtl','superheat','subcooling']:
#             if k in units:
#                 s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
#             else:
#                 s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
#         return s.rstrip()
#         
#     cpdef State copy(self):
#         """
#         Make a copy of this State class
#         """
#         cdef State S = State(self.Fluid,dict(T=self.T_,D=self.rho_))
#         S.phase = self.phase
#         return S
#     
# def rebuildState(d):
#     S=State(d['Fluid'],{'T':d['T'],'D':d['rho']},phase=d['phase'])
#     return S
#     
