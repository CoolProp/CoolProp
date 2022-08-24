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

#
# Ensure that all deprecation warnings triggered in __main__ are displayed by default
# Note: Python deprecation warning behavior varies by version
#       3.7+       - Show all deprecation warnings
#       3.2 to 3.6 - Hide all deprecation warnings
#       before 3.2 - Show all deprecation warnings
#       The filter below will ensure that all dep. warnings show
#       in all versions of Python
#
warnings.filterwarnings('default', category=DeprecationWarning, module='__main__')

from .typedefs cimport CoolPropDbl

try:
    import numpy as np
    _numpy_supported = True
except ImportError:
    _numpy_supported = False

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Python.h":
    char* __FILE__

cdef extern from "Python.h":
    int __LINE__

cdef extern from "CoolPropTools.h":
    double get_HUGE()

cdef extern from "CoolPropTools.h" namespace "CoolProp":
    bint _ValidNumber "ValidNumber"(double)

cdef extern from "Configuration.h" namespace "CoolProp":
    string _get_config_as_json_string "CoolProp::get_config_as_json_string"() except +
    void _set_config_as_json_string "CoolProp::set_config_as_json_string"(string) except +
    string _config_key_description "CoolProp::config_key_description"(string) except +

    void _set_config_string "CoolProp::set_config_string"(constants_header.configuration_keys,string) except +
    void _set_config_double "CoolProp::set_config_double"(constants_header.configuration_keys,double) except +
    void _set_config_bool "CoolProp::set_config_bool"(constants_header.configuration_keys,bint) except +
    void _set_config_int "CoolProp::set_config_int"(constants_header.configuration_keys,int) except +

    string _get_config_string "CoolProp::get_config_string"(constants_header.configuration_keys) except +
    double _get_config_double "CoolProp::get_config_double"(constants_header.configuration_keys) except +
    bint _get_config_bool "CoolProp::get_config_bool"(constants_header.configuration_keys) except +
    int _get_config_int "CoolProp::get_config_int"(constants_header.configuration_keys) except +

cdef extern from "DataStructures.h" namespace "CoolProp":
    string _get_parameter_information "CoolProp::get_parameter_information"(int, string) except +
    int _get_parameter_index "CoolProp::get_parameter_index"(string) except +
    int _get_phase_index "CoolProp::get_phase_index"(string) except +
    bint _is_trivial_parameter "CoolProp::is_trivial_parameter"(int) except +
    constants_header.input_pairs _generate_update_pair "CoolProp::generate_update_pair"(constants_header.parameters key1, double value1, constants_header.parameters key2, double value2, double &out1, double &out2) except +

cdef extern from "CoolPropLib.h":
    double _Props "Props"(const char* Output, const char Name1, double Prop1, const char Name2, double Prop2, const char* Ref)

cdef extern from "CoolProp.h" namespace "CoolProp":
    double _Props1SI "CoolProp::Props1SI"(string Ref, string Output)
    double _PropsSI "CoolProp::PropsSI"(string Output, string Name1, double Prop1, string Name2, double Prop2, string FluidName)
    string _PhaseSI "CoolProp::PhaseSI"(string Name1, double Prop1, string Name2, double Prop2, string FluidName)
    vector[vector[double]] _PropsSImulti "CoolProp::PropsSImulti"(vector[string] Outputs, string Name1, vector[double] Prop1, string Name2, vector[double] Prop2, string backend, vector[string] FluidName, vector[double] fractions)
    string _get_global_param_string "CoolProp::get_global_param_string"(string ParamName) except +
    int _get_debug_level "CoolProp::get_debug_level"() except +
    void _set_debug_level "CoolProp::set_debug_level"(int level) except +
    string _get_fluid_param_string "CoolProp::get_fluid_param_string"(string ParamName, string FluidName) except +
    void _extract_backend "CoolProp::extract_backend"(string input, string backend, string fluids) except +
    string _extract_fractions "CoolProp::extract_fractions"(string input, vector[double] fractions) except +
    void _set_reference_stateS "CoolProp::set_reference_stateS"(string, string) except +
    void _set_reference_stateD "CoolProp::set_reference_stateD"(string, double, double, double, double) except +
    double _saturation_ancillary "CoolProp::saturation_ancillary"(string, string, int, string, double) except +
    bint _add_fluids_as_JSON "CoolProp::add_fluids_as_JSON"(const string backend, const string JSON) except +

cdef extern from "HumidAirProp.h" namespace "HumidAir":
    double _HAPropsSI "HumidAir::HAPropsSI"(string OutputName, string Input1Name, double Input1, string Input2Name, double Input2, string Input3Name, double Input3)
    double _HAProps "HumidAir::HAProps"(string OutputName, string Input1Name, double Input1, string Input2Name, double Input2, string Input3Name, double Input3)
    double _HAProps_Aux "HumidAir::HAProps_Aux"(const char* Name,double T, double p, double W, char *units)
    double _cair_sat "HumidAir::cair_sat"(double T)

cdef extern from "Backends/Helmholtz/MixtureParameters.h" namespace "CoolProp":
    string _get_mixture_binary_pair_data "CoolProp::get_mixture_binary_pair_data"(const string CAS1, const string CAS2, const string key) except +
    void _set_mixture_binary_pair_data "CoolProp::set_mixture_binary_pair_data"(const string CAS1, const string CAS2, const string key, const double val) except +
    void _apply_simple_mixing_rule "CoolProp::apply_simple_mixing_rule"(const string &CAS1, const string &CAS2, const string &rule) except +
    void _set_departure_functions "CoolProp::set_departure_functions"(const string &functions) except +
    void _set_interaction_parameters "CoolProp::set_interaction_parameters"(const string &data) except +

cdef extern from "Backends/PCSAFT/PCSAFTLibrary.h" namespace "CoolProp":
    string _get_mixture_binary_pair_pcsaft "CoolProp::get_mixture_binary_pair_pcsaft"(const string CAS1, const string CAS2, const string key) except +
    void _set_mixture_binary_pair_pcsaft "CoolProp::set_mixture_binary_pair_pcsaft"(const string CAS1, const string CAS2, const string key, const double val) except +

from .constants import *
from .constants_header cimport *
from . cimport constants_header

cdef bint iterable(object a):
    """
    If numpy is supported, this function returns true if the argument is a
    numpy array or another iterable, otherwise just checks if list or tuple
    """
    if _numpy_supported:
        return isinstance(a,(list, tuple, np.ndarray))
    else:
        return isinstance(a,(list, tuple))

cdef ndarray_or_iterable(object input):
    if _numpy_supported:
        return np.squeeze(np.array(input))
    else:
        return input

include "HumidAirProp.pyx"
include "AbstractState.pyx"

def set_reference_state(string FluidName, *args):
    """
    Accepts one of two signatures:

    Type #1 (A Python wrapper of :cpapi:`CoolProp::set_reference_stateS`):

    set_reference_state(FluidName,reference_state)

    FluidName The name of the fluid
    param reference_state The reference state to use, one of

    ==========   ===========================================
    ``IIR``      (h=200 kJ/kg, s=1 kJ/kg/K at 0C sat. liq.)
    ``ASHRAE``   (h=0,s=0 @ -40C sat liq)
    ``NBP``      (h=0,s=0 @ 1.0 bar sat liq.)
    ==========   ===========================================

    Type #2 (A Python wrapper of :cpapi:`CoolProp::set_reference_stateD`):

    set_reference_state(FluidName,T0,rhomolar,hmolar0,smolar0)

    .. note::

        Only supported for internal backend currently

    ``FluidName`` The name of the fluid

    ``T0`` The temperature at the reference point [K]

    ``rhomolar`` The density at the reference point [mol/m^3]

    ``hmolar0`` The enthalpy at the reference point [J/mol]

    ``smolar0`` The entropy at the reference point [J/mol/K]
    """

    cdef bytes _param
    cdef int retval

    if len(args) == 1:
        _set_reference_stateS(FluidName, args[0])
    elif len(args) == 4:
        _set_reference_stateD(FluidName, args[0], args[1], args[2], args[3])
    else:
        raise ValueError(b'Invalid number of inputs')

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

cpdef tuple generate_update_pair(constants_header.parameters key1, double value1, constants_header.parameters key2, double value2):
    """
    This function will generate an input pair to the update() function given the key, value pairs for both inputs
    """
    cdef constants_header.input_pairs pair
    cdef double out1 = -1000000000, out2 = -100000000000
    pair = _generate_update_pair(key1, value1, key2, value2, out1, out2)
    return pair, out1, out2

cpdef string get_config_as_json_string():
    """
    Obtain a json formulation of the internal configuration in CoolProp

    Values can be set by passing a modified json library (converted to string) to set_config_as_json_string
    """
    return _get_config_as_json_string()

cpdef string config_key_description(string key):
    """
    Obtain the string description for a configuration key.  Python wrapper of C++ function :cpapi:`CoolProp::config_key_description`
    """
    return _config_key_description(key)

cpdef set_config_as_json_string(string s):
    """
    Set the internal configuration in CoolProp from a json data string

    Current state can be obtained by calling get_config_as_json_string
    """
    _set_config_as_json_string(s)

cpdef set_config_double(constants_header.configuration_keys key, double value):
    """ Set configuration key that is a double-precision float;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_double` """
    _set_config_double(key, value)

cpdef set_config_string(constants_header.configuration_keys key, string value):
    """ Set a configuration key that is a string;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_string` """
    _set_config_string(key, value)

cpdef set_config_bool(constants_header.configuration_keys key, bint value):
    """ Set a configuration key that is a boolean;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_bool` """
    _set_config_bool(key, value)

cpdef set_config_int(constants_header.configuration_keys key, int value):
    """ Set a configuration key that is an integer;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_int` """
    _set_config_int(key, value)

cpdef double get_config_double(constants_header.configuration_keys key):
    """ Get a configuration key that is a double-precision float;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_double` """
    return _get_config_double(key)

cpdef string get_config_string(constants_header.configuration_keys key):
    """ Get a configuration key that is a string;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_string` """
    return _get_config_string(key)

cpdef bint get_config_bool(constants_header.configuration_keys key):
    """ Get a configuration key that is a boolean;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_bool` """
    return _get_config_bool(key)

cpdef int get_config_int(constants_header.configuration_keys key):
    """ Get a configuration key that is an integer;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_int` """
    return _get_config_int(key)

cpdef int get_parameter_index(string key):
    return _get_parameter_index(key)

cpdef int get_phase_index(string key):
    return _get_phase_index(key)

cpdef string get_parameter_information(int key, string info):
    return _get_parameter_information(key, info)

cpdef string get_mixture_binary_pair_data(CAS1, CAS2, key) except *:
    """
    Obtain mixture interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::get_mixture_binary_pair_data`
    """
    return _get_mixture_binary_pair_data(CAS1, CAS2, key)

cpdef set_mixture_binary_pair_data(CAS1, CAS2, key, val):
    """
    Set mixture interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::set_mixture_binary_pair_data`
    """
    _set_mixture_binary_pair_data(CAS1, CAS2, key, val)

cpdef string get_mixture_binary_pair_pcsaft(CAS1, CAS2, key) except *:
    """
    Obtain mixture PC-SAFT interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::get_mixture_binary_pair_pcsaft`
    """
    _get_mixture_binary_pair_pcsaft(CAS1, CAS2, key)

cpdef set_mixture_binary_pair_pcsaft(CAS1, CAS2, key, val):
    """
    Set mixture PC-SAFT interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::set_mixture_binary_pair_pcsaft`
    """
    _set_mixture_binary_pair_pcsaft(CAS1, CAS2, key, val)

cpdef add_fluids_as_JSON(backend, JSONstring):
    """
    Add fluids in a JSON-formatted string format. Python wrapper of C++ function :cpapi:`CoolProp::add_fluids_as_JSON`
    """
    _add_fluids_as_JSON(backend, JSONstring)

cpdef get_global_param_string(string param):
    return _get_global_param_string(param)

cpdef is_trivial_parameter(int key):
    return _is_trivial_parameter(key)

cpdef get_fluid_param_string(string fluid, string param):
    return _get_fluid_param_string(fluid, param)

cpdef apply_simple_mixing_rule(CAS1, CAS2, rule):
    """
    Apply simple mixing rule.  Currently linear or Lorentz-Berthelot.  Python wrapper of C++ function :cpapi:`CoolProp::apply_simple_mixing_rule`
    """
    _apply_simple_mixing_rule(CAS1, CAS2, rule)

cpdef set_departure_functions(functions):
    """
    Specify the departure terms as JSON. Python wrapper of C++ function :cpapi:`CoolProp::set_departure_functions`
    """
    _set_departure_functions(functions)

cpdef set_interaction_parameters(data):
    """
    Specify the binary interaction terms as JSON. Python wrapper of C++ function :cpapi:`CoolProp::set_interaction_parameters`
    """
    _set_interaction_parameters(data)    

cpdef double saturation_ancillary(string name, string output, int Q, string input, double value):
    """
    Return a value from the saturation ancillary equations; python wrapper of :cpapi:`CoolProp::saturation_ancillary`
    """
    return _saturation_ancillary(name, output, Q, input, value)

cpdef __Props_err1(fcn, in1,in2):
    errstr = _get_global_param_string(b'errstring')
    if not len(errstr) == 0:
        raise ValueError("{err:s} :: inputs were :\"{in1:s}\",\"{in2:s}\"".format(err= errstr,in1=in1,in2=in2))
    else:
        raise ValueError("{fcn:s} failed ungracefully with inputs:\"{in1:s}\",\"{in2:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(fcn=fcn, in1=in1,in2=in2))

cpdef __Props_err2(fcn, in1, in2, in3, in4, in5, in6):
    errstr = _get_global_param_string(b'errstring')
    if not len(errstr) == 0:
        raise ValueError("{err:s}".format(err=errstr))
    else:
        raise ValueError("{fcn:s} failed ungracefully :: inputs were:\"{in1:s}\",\"{in2:s}\",{in3:0.16e},\"{in4:s}\",{in5:0.16e},\"{in6:s}\"; please file a ticket at https://github.com/CoolProp/CoolProp/issues".format(fcn = fcn, in1=in1,in2=in2,in3=in3,in4=in4,in5=in5,in6=in6))

cpdef Props(in1, in2, in3 = None, in4 = None, in5 = None, in6 = None):
    """
    This function is deprecated, use PropsSI instead
    """
    import warnings
    # Issue deprecation warning....
    dep_warning = "Props() function is deprecated; Use the PropsSI() function"
    warnings.warn(dep_warning, DeprecationWarning)
    # ...but process Props() function anyway...
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

cpdef PhaseSI(in1, in2, in3, in4, in5):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::PhaseSI`

    Does not support vectorization of the inputs like PropsSI
    """
    return _PhaseSI(in1, in2, in3, in4, in5)

cpdef PropsSI(in1, in2, in3 = None, in4 = None, in5 = None, in6 = None, in7 = None):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::PropsSI` .
    """
    cdef vector[string] vin1
    cdef vector[double] fractions, vval1, vval2
    cdef double val
    cdef string backend, fluid, delimitedfluids
    cdef bool is_iterable1, is_iterable3, is_iterable5

    # Two parameter inputs
    if in3 is None and in4 is None and in5 is None and in6 is None and in7 is None:
        val = _Props1SI(in1, in2)
        if not _ValidNumber(val):
            __Props_err1("PropsSI", in1, in2)
        else:
            return val
    # Six parameter inputs
    elif in7 is None:
        is_iterable1 = iterable(in1)
        is_iterable3 = iterable(in3)
        is_iterable5 = iterable(in5)

        if _numpy_supported and is_iterable3 and isinstance(in3, np.ndarray) and (np.prod(in3.shape) != max(in3.shape)):
            raise ValueError("Input 3 is not one-dimensional")
        if _numpy_supported and is_iterable5 and isinstance(in5, np.ndarray) and (np.prod(in5.shape) != max(in5.shape)):
            raise ValueError("Input 5 is not one-dimensional")

        if is_iterable1 or is_iterable3 or is_iterable5:
            # Prepare the output datatype
            if not is_iterable1:
                vin1.push_back(in1)
            else:
                vin1 = in1

            # Resize state variable inputs
            if is_iterable3 and is_iterable5:
                if len(in3) != len(in5):
                    raise TypeError("Sizes of Prop1 {n1:d} and Prop2 {n2:d} to PropsSI are not the same".format(n1 = len(in3), n2 = len(in5)))
                else:
                    vval1 = in3
                    vval2 = in5
            elif is_iterable3 and not is_iterable5:
                vval1 = in3
                vval2.resize(len(in3))
                templist = [in5]*len(in3)
                vval2 = templist
            elif is_iterable5 and not is_iterable3:
                vval1.resize(len(in5))
                templist = [in3]*len(in5)
                vval1 = templist
                vval2 = in5
            else:
                vval1.resize(1)
                vval1[0] = in3
                vval2.resize(1)
                vval2[0] = in5

            # Extract the backend and the fluid from the input string
            _extract_backend(in6, backend, fluid)

            # Extract the fractions
            fractions.push_back(1.0)
            delimitedfluids = _extract_fractions(fluid, fractions)

            # Extract the fluids
            fluids = delimitedfluids.split('&')

            # Call the function - this version takes iterables
            outmat = _PropsSImulti(vin1, in2, vval1, in4, vval2, backend, fluids, fractions)

            # Check that we got some output
            if outmat.empty():
                raise ValueError(_get_global_param_string(b'errstring'))

            return ndarray_or_iterable(outmat)
        else:
            # This version takes doubles
            val = _PropsSI(in1, in2, in3, in4, in5, in6)
            if not _ValidNumber(val):
                __Props_err2("PropsSI", in1, in2, in3, in4, in5, in6)
            else:
                return val
    else:
        raise ValueError('input #7 cannot be provided')

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
    return _get_global_param_string(b"FluidsList").split(',')

cpdef get_aliases(string Fluid):
    """
    Return a comma separated string of aliases for the given fluid
    """
    cdef bytes _Fluid = Fluid.encode('ascii')
    return [F for F in _get_fluid_param_string(_Fluid, b'aliases').split(',')]

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
    return _get_fluid_param_string(Fluid,b'REFPROP_name')

cpdef string get_BibTeXKey(string Fluid, string key) except *:
    """
    Return the BibTeX key for the given fluid.

    The possible keys are

    * ``EOS``
    * ``CP0``
    * ``VISCOSITY``
    * ``CONDUCTIVITY``
    * ``ECS_LENNARD_JONES``
    * ``ECS_FITS``
    * ``SURFACE_TENSION``
    * ``MELTING_LINE``

    BibTeX keys refer to the BibTeX file in the trunk/CoolProp folder

    Returns
    -------
    key, string
         empty string if Fluid not in CoolProp, "Bad key" if key is invalid
    """
    return _get_fluid_param_string(Fluid, b"BibTeX-"+key.encode('ascii'))

cpdef string get_errstr():
    """
    Return the current error string
    """
    return _get_global_param_string(b"errstring")

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


cpdef extract_backend(string in_str):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::extract_backend` .
    """
    cdef string bck, fld
    # Extract the backend and the fluid from the input string
    _extract_backend(in_str, bck, fld)
    return bck, fld


cpdef extract_fractions(string flds):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::extract_fractions` .
    """
    cdef vector[double] frcs
    cdef string del_flds
    # Extract the fractions
    #frcs.clear()
    frcs.push_back(1.0)
    del_flds = _extract_fractions(flds, frcs)
    # Extract the fluids
    fluids = del_flds.split('&')
    return fluids,frcs


cdef toSI(constants_header.parameters key, double val):
    """
    Convert a value in kSI system to SI system (supports a limited subset of variables)
    """
    if key in [iT, iDmass, iQ]:
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

    .. warning::

        This class is deprecated.  You should use :py:class:`CoolProp.AbstractState` instead

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

    def __init__(self, object _Fluid, dict StateDict, object phase = None, backend = None):
        """
        Parameters
        ----------
        Fluid : string
        StateDict : dictionary
            The state of the fluid - passed to the update function; if None, does not do a state update
        phase : string
            DEPRECATED : this input is ignored
        backend : string
            The CoolProp backend that should be used, one of "HEOS" (default), "REFPROP", "INCOMP", "BRINE", etc.
        """
        cdef string Fluid = _Fluid


        if Fluid == b'none':
            return
        else:
            if b'::' in <bytes>Fluid:
                backend, _Fluid = (<bytes>Fluid).split(b'::')
            elif backend is None:
                backend = u'?'

            self.set_Fluid(_Fluid, backend)
        self.Fluid = Fluid

        # Parse the inputs provided
        if StateDict is not None:
            self.update(StateDict)

        if phase is None:
            self.phase = b'??'
        else:
            self.phase = phase.encode('ascii')

        # Set the phase flag
        if self.phase.lower() == 'gas':
            self.pAS.specify_phase(constants_header.iphase_gas)
        elif self.phase.lower() == 'liquid':
            self.pAS.specify_phase(constants_header.iphase_liquid)

#     def __reduce__(self):
#         d={}
#         d['Fluid']=self.Fluid
#         d['T']=self.T_
#         d['rho']=self.rho_
#         d['phase'] = self.phase
#         return rebuildState,(d,)

    cpdef set_Fluid(self, string Fluid, string backend):

        cdef object _Fluid = Fluid
        cdef object _backend = backend
        cdef bint set_fractions = False
        new_fluid = []
        fracs = []
        if '[' in _Fluid and ']' in _Fluid:
            pairs = _Fluid.split('&')
            for pair in pairs:
                fluid, frac = pair.split('[')
                new_fluid.append(fluid)
                fracs.append(float(frac.strip(']')))
            _Fluid = '&'.join(new_fluid)
            set_fractions = True

        self.pAS = AbstractState(_backend, _Fluid)
        if set_fractions:
            self.pAS.set_mole_fractions(fracs)

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
        self.pAS.update(HmassP_INPUTS, h*1000, p*1000)
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
        return self.pAS.Q()
    property Q:
        """ The quality [-] """
        def __get__(self):
            return self.get_Q()

    cpdef double get_MM(self) except *:
        """ Get the mole mass [kg/kmol] or [g/mol] """
        return self.pAS.molar_mass()*1000
    property MM:
        """ The molar mass [kg/kmol] or [g/mol] """
        def __get__(self):
            return self.get_MM()

    cpdef double get_rho(self) except *:
        """ Get the density [kg/m^3] """
        return self.pAS.rhomass()
    property rho:
        """ The density [kg/m^3] """
        def __get__(self):
            return self.Props(iDmass)

    cpdef double get_p(self) except *:
        """ Get the pressure [kPa] """
        return self.pAS.p()/1000
    property p:
        """ The pressure [kPa] """
        def __get__(self):
            return self.get_p()

    cpdef double get_T(self) except *:
        """ Get the temperature [K] """
        return self.pAS.T()
    property T:
        """ The temperature [K] """
        def __get__(self):
            return self.get_T()

    cpdef double get_h(self) except *:
        """ Get the specific enthalpy [kJ/kg] """
        return self.pAS.hmass()/1000
    property h:
        """ The specific enthalpy [kJ/kg] """
        def __get__(self):
            return self.get_h()

    cpdef double get_u(self) except *:
        """ Get the specific internal energy [kJ/kg] """
        return self.pAS.umass()/1000
    property u:
        """ The internal energy [kJ/kg] """
        def __get__(self):
            return self.get_u()

    cpdef double get_s(self) except *:
        """ Get the specific enthalpy [kJ/kg/K] """
        return self.pAS.smass()/1000
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
        return self.pAS.cpmass()/1000
    property cp:
        """ The specific heat at constant pressure  [kJ/kg/K] """
        def __get__(self):
            return self.get_cp()

    cpdef double get_cv(self) except *:
        """ Get the specific heat at constant volume  [kJ/kg/K] """
        return self.pAS.cvmass()/1000
    property cv:
        """ The specific heat at constant volume  [kJ/kg/K] """
        def __get__(self):
            return self.get_cv()

    cpdef double get_speed_sound(self) except *:
        """ Get the speed of sound  [m/s] """
        return self.Props(ispeed_sound)

    cpdef double get_visc(self) except *:
        """ Get the viscosity, in [Pa-s]"""
        return self.pAS.viscosity()
    property visc:
        """ The viscosity, in [Pa-s]"""
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self) except *:
        """ Get the thermal conductivity, in [kW/m/K]"""
        return self.pAS.conductivity()/1000
    property k:
        """ The thermal conductivity, in [kW/m/K]"""
        def __get__(self):
            return self.get_cond()

    cpdef get_Tsat(self, double Q = 1):
        """
        Get the saturation temperature, in [K]

        Returns ``None`` if pressure is not within the two-phase pressure range
        """
        cdef State state = State(self.Fluid, None)
        cdef double pc = state.Props(iP_critical)
        cdef double pt
        try:
            pt = state.Props(iP_triple)
        except ValueError:
            pt = -1
        if _ValidNumber(pc) and _ValidNumber(pt):
            if self.p_ > 0.001*pc or self.p_ < 0.001*pt:
                return None
            else:
                state.update(dict(P=self.p_,Q=Q))
                return state.T
        else:
            state.update(dict(P=self.p_,Q=Q))
            return state.T
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
        cdef long IT = b'T'
        cdef long ID = b'D'
        import CoolProp as CP

        print('Call to the Python call layer (CoolProp.CoolProp.Props)')
        print("'M' involves basically no computational effort and is a good measure of the function call overhead")
        keys = ['H','P','S','U','C','O','V','L','M','d(P)/d(T)|Dmolar']
        for key in keys:
            t1=clock()
            for i in range(N):
                CP.PropsSI(key,b'T',self.T_,b'D',self.rho_,Fluid)
            t2=clock()
            print('Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6))

        print('Direct c++ call to CoolProp without the Python call layer (_Props function)')
        print("'M' involves basically no computational effort and is a good measure of the function call overhead")
        keys = ['H','P','S','U','C','O','V','L','M','C0','d(P)/d(T)|Dmolar']
        for key in keys:
            t1=clock()
            for i in range(N):
                _PropsSI(key,b'T',self.T_,b'D',self.rho_,Fluid)
            t2=clock()
            print('Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6))

        print('Call to the c++ layer using integers')
        keys = [iHmass, iP,iSmass,iUmass]
        for key in keys:
            t1=clock()
            for i in range(N):
                self.pAS.update(DmassT_INPUTS,self.rho_,self.T_)
                self.pAS.keyed_output(key)
            t2=clock()
            print('Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6))

        print('Call to the AbstractState for molar mass (fast)')
        t1=clock()
        for i in range(N):
            self.pAS.keyed_output(imolar_mass)
        t2=clock()
        print('Elapsed time for {0:d} calls at {1:g} us/call'.format(N, (t2-t1)/N*1e6))

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
        units = {
        'T': 'K',
        'p': 'kPa',
        'rho': 'kg/m^3',
        'Q': 'kg/kg',
        'h': 'kJ/kg',
        'u': 'kJ/kg',
        's': 'kJ/kg/K',
        'visc': 'Pa-s',
        'k': 'kW/m/K',
        'cp': 'kJ/kg/K',
        'cp0': 'kJ/kg/K',
        'cv': 'kJ/kg/K',
        'dpdT': 'kPa/K',
        'Tsat': 'K',
        'superheat': 'K',
        'subcooling': 'K',
        'MM': 'kg/kmol'
        }
        s = 'phase = '+self.phase.decode('ascii')+'\n'
        for k in ['T','p','rho','Q','h','u','s','visc','k','cp','cp0','cv','dpdT','Prandtl','superheat','subcooling','MM']:
            if k in units:
                s += k + ' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s += k + ' = '+str(getattr(self,k))+' NO UNITS'+'\n'
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
