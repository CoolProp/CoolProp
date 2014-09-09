#cython: embedsignature = True, c_string_type = unicode, c_string_encoding = ascii
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

try:
    import numpy as np
    _numpy_supported = True
except ImportError:
    _numpy_supported = False

from libcpp.string cimport string
from libcpp.vector cimport vector

from constants import *
from constants_header cimport *
cimport constants_header

cdef bint iterable(object a):
    """
    If numpy is supported, this function retuns true if the argument is a 
    numpy array or another iterable, otherwise just checks if list or tuple
    """
    if _numpy_supported:
        return isinstance(a,(list, tuple, np.ndarray))
    else:
        return isinstance(a,(list, tuple))

cdef ndarray_or_iterable(object input):
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

cpdef int get_parameter_index(string key):
    return _get_parameter_index(key)
    
cpdef string get_parameter_information(int key, string info):
    return _get_parameter_information(key, info)

cpdef get_global_param_string(string param):
    return _get_global_param_string(param)
     
cpdef get_fluid_param_string(string fluid, string param):
    return _get_fluid_param_string(fluid, param)
     
cpdef __Props_err1(fcn, in1,in2):
    errstr = _get_global_param_string('errstring')
    if not len(errstr) == 0:
        raise ValueError("{err:s} :: inputs were :\"{in1:s}\",\"{in2:s}\"".format(err= errstr,in1=in1,in2=in2))
    else:
        raise ValueError("{fcn:s} failed ungracefully with inputs:\"{in1:s}\",\"{in2:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(fcn=fcn, in1=in1,in2=in2))
        
cpdef __Props_err2(fcn, in1, in2, in3, in4, in5, in6):
    errstr = _get_global_param_string('errstring')
    if not len(errstr) == 0:
        raise ValueError("{err:s}".format(err=errstr))
    else:
        raise ValueError("{fcn:s} failed ungracefully :: inputs were:\"{in1:s}\",\"{in2:s}\",{in3:0.16e},\"{in4:s}\",{in5:0.16e},\"{in6:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(fcn = fcn, in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))

cpdef Props(in1, in2, in3 = None, in4 = None, in5 = None, in6 = None):
    """
    ${CoolProp::Props}
    """ 
    import warnings
    dep_warning = "Props() function is deprecated; Use the PropsSI() function"
    warnings.warn_explicit(dep_warning, category=UserWarning, filename='CoolProp.pyx', lineno = -1)
    if len(in2) != 1: 
        raise ValueError('Length of input name #1 must be 1 character')
    if len(in4) != 1: 
        raise ValueError('Length of input name #2 must be 1 character')
    cdef char* c1 = (<bytes>in2)
    cdef char* c2 = (<bytes>in4)
    val = _Props(in1, c1[0], in3, c2[0], in5, in6)
    if not _ValidNumber(val):
        __Props_err2("Props", in1, in2, in3, in4, in5, in6)
    else:
        return val

cpdef PropsSI(in1, in2, in3 = None, in4 = None, in5 = None, in6 = None, in7 = None):
    """
    $$PropsSI$$
    """ 
    cdef double val
    
    # Two parameter inputs
    if in3 is None and in4 is None and in5 is None and in6 is None and in7 is None:
        val = _Props1SI(in1, in2)
        if not _ValidNumber(val):
            __Props_err1("PropsSI", in1, in2)
        else:
            return val
    # Six parameter inputs
    elif in7 is None:
        if iterable(in3) and iterable(in5):
            # This version takes iterables
            return _PropsSII(in1, in2, in3, in4, in5, in6)
        elif iterable(in3) and not(iterable(in5)):
            i5 = [in5]*len(in3)
            # This version takes iterables
            return _PropsSII(in1, in2, in3, in4, i5, in6)
        elif iterable(in5) and not(iterable(in3)):
            i3 = [in3]*len(in5)
            # This version takes iterables
            return _PropsSII(in1, in2, i3, in4, in5, in6)
        else:
            # This version takes doubles
            val = _PropsSI(in1, in2, in3, in4, in5, in6)
            if not _ValidNumber(val):
                __Props_err2("PropsSI", in1, in2, in3, in4, in5, in6)
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

cpdef get_aliases(string Fluid):
    """
    Return a comma separated string of aliases for the given fluid
    """
    cdef bytes _Fluid = Fluid.encode('ascii')
    return [F for F in _get_fluid_param_string(_Fluid, 'aliases').split(',')]

cpdef string get_REFPROPname(string Fluid):
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
    return _get_fluid_param_string(Fluid,'REFPROP_name')

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

cdef toSI(constants_header.parameters key, double val):
    """
    Convert a value in kSI system to SI system (supports a limited subset of variables)
    """
    if key in [iT, iDmass]:
        return val
    elif key in [iP, iHmass, iSmass, iUmass]:
        return val*1000
    else:
        raise KeyError('key is invalid to toSI')

#A dictionary mapping parameter index to string for use with non-CoolProp fluids
cdef dict paras = {iDmass : 'D',
                   iQ : 'Q',
                   imolar_mass : 'M',
                   iT : 'T',
                   iHmass : 'H',
                   iP : 'P',
                   iCpmass : 'C',
                   iCp0mass : 'C0',
                   iCvmass : 'O',
                   iviscosity : 'V',
                   iconductivity : 'L',
                   ispeed_sound: 'A',
                   iSmass : 'S',
                   iUmass : 'U'
}

cdef dict paras_inverse = {v:k for k,v in paras.iteritems()}

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
                backend = u'?'
            self.set_Fluid(Fluid, backend)
        self.Fluid = _Fluid
        self.phase = phase
        if phase is None:
            self.phase = u'??'.encode('ascii')
        # Parse the inputs provided
        self.update(StateDict)
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
        
    cpdef update_ph(self, double p, double h):
        """
        Use the pressure and enthalpy directly
        
        Parameters
        ----------
        p: float
            Pressure (absolute) [kPa]
        h: float
            Enthalpy [kJ/kg]
        """
        self.pAS.update(HmassP_INPUTS, p*1000, h*1000)
        self.T_ = self.pAS.T()
        self.rho_ = self.pAS.rhomass()
            
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
        
    cpdef update(self, dict params):
        """
        Parameters
        params, dictionary 
            A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
            for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
        """
        
        # Convert to integer_pair input
            
        cdef double p, val1, val2, o1 = 0, o2 = 0
        cdef long iInput1, iInput2
        cdef bytes errstr
        cdef constants_header.input_pairs input_pair
        
        # Convert inputs to input pair
        items = list(params.items())
        key1 = paras_inverse[items[0][0]]
        key2 = paras_inverse[items[1][0]]
        # Convert to SI units
        val1 = toSI(key1, items[0][1])
        val2 = toSI(key2, items[1][1])
        
        input_pair = _generate_update_pair(key1, val1, key2, val2, o1, o2)
        self.pAS.update(input_pair, o1, o2);
        
        self.T_ = self.pAS.T()
        self.p_ =  self.pAS.p()/1000;
        self.rho_ = self.pAS.rhomass()
        
    cpdef long Phase(self) except *:
        """
        Returns an integer flag for the phase of the fluid, where the flag value
        is one of iLiquid, iSupercritical, iGas, iTwoPhase
        
        These constants are defined in the phase_constants module, and are imported
        into this module
        """
        
        if self.is_CPFluid:
            return self.pAS.phase()
        else:
            raise NotImplementedError("Phase not defined for fluids other than CoolProp fluids")
        
    cpdef double Props(self, constants_header.parameters iOutput) except *: 
        if iOutput<0:
            raise ValueError('Your output is invalid') 
        return self.pAS.keyed_output(iOutput)
            
    cpdef double get_Q(self) except *:
        """ Get the quality [-] """
        return self.Props(iQ)
    property Q:
        """ The quality [-] """
        def __get__(self):
            return self.get_Q()
    
    cpdef double get_MM(self) except *:
        """ Get the mole mass [kg/kmol] or [g/mol] """
        return self.Props(imolar_mass)*1000
    property MM:
        """ The molar mass [kg/kmol] or [g/mol] """
        def __get__(self):
            return self.get_MM()
            
    cpdef double get_rho(self) except *:
        """ Get the density [kg/m^3] """ 
        return self.Props(iDmass)
    property rho:
        """ The density [kg/m^3] """
        def __get__(self):
            return self.Props(iDmass)
            
    cpdef double get_p(self) except *:
        """ Get the pressure [kPa] """ 
        return self.Props(iP)/1000
    property p:
        """ The pressure [kPa] """
        def __get__(self):
            return self.get_p()
    
    cpdef double get_T(self) except *: 
        """ Get the temperature [K] """
        return self.Props(iT)
    property T:
        """ The temperature [K] """
        def __get__(self):
            return self.get_T()
    
    cpdef double get_h(self) except *: 
        """ Get the specific enthalpy [kJ/kg] """
        return self.Props(iHmass)/1000
    property h:
        """ The specific enthalpy [kJ/kg] """
        def __get__(self):
            return self.get_h()
          
    cpdef double get_u(self) except *: 
        """ Get the specific internal energy [kJ/kg] """
        return self.Props(iUmass)/1000
    property u:
        """ The internal energy [kJ/kg] """
        def __get__(self):
            return self.get_u()
            
    cpdef double get_s(self) except *: 
        """ Get the specific enthalpy [kJ/kg/K] """
        return self.Props(iSmass)/1000
    property s:
        """ The specific enthalpy [kJ/kg/K] """
        def __get__(self):
            return self.get_s()
    
    cpdef double get_cp0(self) except *:
        """ Get the specific heat at constant pressure for the ideal gas [kJ/kg/K] """
        return self.Props(iCp0mass)/1000
    property cp0:
        """ The ideal-gas specific heat at constant pressure  [kJ/kg/K] """
        def __get__(self):
            return self.get_cp0()
    
    cpdef double get_cp(self) except *: 
        """ Get the specific heat at constant pressure  [kJ/kg/K] """
        return self.Props(iCpmass)/1000
    property cp:
        """ The specific heat at constant pressure  [kJ/kg/K] """
        def __get__(self):
            return self.get_cp()
            
    cpdef double get_cv(self) except *: 
        """ Get the specific heat at constant volume  [kJ/kg/K] """
        return self.Props(iCvmass)/1000
    property cv:
        """ The specific heat at constant volume  [kJ/kg/K] """
        def __get__(self):
            return self.get_cv()
        
    cpdef double get_speed_sound(self) except *: 
        """ Get the speed of sound  [m/s] """
        return self.Props(ispeed_sound)
            
    cpdef double get_visc(self) except *:
        """ Get the viscosity, in [Pa-s]"""
        return self.Props(iviscosity)
    property visc:
        """ The viscosity, in [Pa-s]"""
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self) except *:
        """ Get the thermal conductivity, in [kW/m/K]"""
        return self.Props(iconductivity)/1000
    property k:
        """ The thermal conductivity, in [kW/m/K]"""
        def __get__(self):
            return self.get_cond()
        
    cpdef get_Tsat(self, double Q = 1):
        """ 
        Get the saturation temperature, in [K]
        
        Returns ``None`` if pressure is not within the two-phase pressure range 
        """
        if self.p_ > _Props('pcrit','T',0,'P',0,self.Fluid) or self.p_ < _Props('ptriple','T',0,'P',0, self.Fluid):
            return None 
        else:
            return _Props('T', 'P', self.p_, 'Q', Q, self.Fluid)
    property Tsat:
        """ The saturation temperature (dew) for the given pressure, in [K]"""
        def __get__(self):
            return self.get_Tsat(1.0)
        
    cpdef get_superheat(self):
        """ 
        Get the amount of superheat above the saturation temperature corresponding to the pressure, in [K]
        
        Returns ``None`` if pressure is not within the two-phase pressure range 
        """
        
        Tsat = self.get_Tsat(1) #dewpoint temp
        
        if Tsat is not None:
            return self.T_-Tsat
        else:
            return None
    property superheat:
        """ 
        The amount of superheat above the saturation temperature corresponding to the pressure, in [K]
        
        Returns ``None`` if pressure is not within the two-phase pressure range 
        """
        def __get__(self):    
            return self.get_superheat()
        
    cpdef get_subcooling(self):
        """ 
        Get the amount of subcooling below the saturation temperature corresponding to the pressure, in [K]
        
        Returns ``None`` if pressure is not within the two-phase pressure range 
        """
        
        Tsat = self.get_Tsat(0) #bubblepoint temp
        
        if Tsat is not None:
            return Tsat - self.T_
        else:
            return None
    property subcooling:
        """ 
        The amount of subcooling below the saturation temperature corresponding to the pressure, in [K]
        
        Returns ``None`` if pressure is not within the two-phase pressure range 
        """
        def __get__(self):    
            return self.get_subcooling()
            
    property Prandtl:
        """ The Prandtl number (cp*mu/k) [-] """
        def __get__(self):
            return self.cp * self.visc / self.k
            
    cpdef double get_dpdT(self) except *:
        return self.pAS.first_partial_deriv(iP, iT, iDmolar)/1000;
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
        
    cpdef speed_test(self, int N):
        from time import clock
        cdef int i
        cdef char * k
        cdef long ikey
        cdef bytes Fluid = self.Fluid
        cdef long IT = 'T'
        cdef long ID = 'D'
        import CoolProp as CP
        
        print 'Call to the Python call layer (CoolProp.CoolProp.Props)'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
        for key in keys:
            t1=clock()
            for i in range(N):
                CP.Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
            
        print 'Direct c++ call to CoolProp without the Python call layer (_Props function)'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0','dpdT']
        for key in keys:
            t1=clock()
            for i in range(N):
                _Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
        
        print 'Call to the c++ layer through IProps'
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
        for key in keys:
            t1=clock()
            for i in range(N):
                _IProps(key,iT,self.T_,iD,self.rho_,self.iFluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
            
        print 'Call to the c++ layer using integers'
        keys = [iHmass, iP,iSmass,iUmass]
        for key in keys:
            t1=clock()
            for i in range(N):
                self.pAS.update(DmassT_INPUTS,self.rho_,self.T_)
                self.pAS.keyed_output(key)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
        
        #~ keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0,iDpdT]
        #~ isenabled = _isenabled_TTSE_LUT(<bytes>Fluid)
        #~ _enable_TTSE_LUT(<bytes>Fluid)
        #~ _IProps(iH,iT,self.T_,iD,self.rho_,self.iFluid)
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
    def __str__(self):
        """
        Return a string representation of the state
        """
        units={'T': 'K', 
               'p': 'kPa', 
               'rho': 'kg/m^3',
               'Q':'kg/kg',
               'h':'kJ/kg',
               'u':'kJ/kg',
               's':'kJ/kg/K',
               'visc':'Pa-s',
               'k':'kW/m/K',
               'cp':'kJ/kg/K',
               'cp0':'kJ/kg/K',
               'cv':'kJ/kg/K',
               'dpdT':'kPa/K',
               'Tsat':'K',
               'superheat':'K',
               'subcooling':'K',
               'MM':'kg/kmol'
        }
        s='phase = '+self.phase+'\n'
        for k in ['T','p','rho','Q','h','u','s','visc','k','cp','cp0','cv','dpdT','Prandtl','superheat','subcooling','MM']:
            if k in units:
                s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
        return s.rstrip()
        
    cpdef State copy(self):
        """
        Make a copy of this State class
        """
        cdef State S = State(self.Fluid,dict(T=self.T_,D=self.rho_))
        S.phase = self.phase
        return S
    
def rebuildState(d):
    S=State(d['Fluid'],{'T':d['T'],'D':d['rho']},phase=d['phase'])
    return S
    
