from libcpp.string cimport string
import cython
cimport cython

from libcpp.vector cimport vector

include "AbstractState.pxd"
    
cdef extern from "Python.h":
    char* __FILE__

cdef extern from "Python.h":
    int __LINE__
    
cdef extern from "CoolPropTools.h" namespace "CoolProp":
    bint _ValidNumber "ValidNumber"(double)
    
cdef extern from "Configuration.h" namespace "CoolProp":    
    string _get_config_as_json_string "CoolProp::get_config_as_json_string"() except +
    void _set_config_as_json_string "CoolProp::set_config_as_json_string"(string) except +
    
cdef extern from "DataStructures.h" namespace "CoolProp":    
    string _get_mixture_binary_pair_data "CoolProp::get_mixture_binary_pair_data"(const string CAS1, const string CAS2, const string key) except +
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
    
    # Convenience functions from v4
#     long _get_parameter_index "CoolProp::get_parameter_index" (string param)
#     int _IsFluidType "IsFluidType"(char* Ref, char* Type)
#     string _get_BibTeXKey "CoolProp::get_BibTeXKey"(string Ref, string key)
#     long _get_Fluid_index "CoolProp::get_Fluid_index" (string Fluid)
#     double _IProps "CoolProp::IProps"(long Output, long Name1, double Prop1, long Name2, double Prop2, long Ref)
 
cdef extern from "HumidAirProp.h" namespace "HumidAir":
    double _HAPropsSI "HumidAir::HAPropsSI"(string OutputName, string Input1Name, double Input1, string Input2Name, double Input2, string Input3Name, double Input3)
    double _HAProps "HumidAir::HAProps"(string OutputName, string Input1Name, double Input1, string Input2Name, double Input2, string Input3Name, double Input3)
    double _HAProps_Aux "HumidAir::HAProps_Aux"(const char* Name,double T, double p, double W, char *units)
    double _cair_sat "HumidAir::cair_sat"(double T)
       
cdef class State:
    cdef AbstractState pAS
    cdef readonly bytes Fluid, phase
    cdef int iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_
    
    cpdef set_Fluid(self, string Fluid, string backend)
    cpdef speed_test(self, int N)
    cpdef update(self, dict params)
    cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
    cpdef State copy(self)
    cpdef double Props(self, constants_header.parameters iOutput) except *
    cpdef long Phase(self) except *
    cpdef double get_Q(self) except *
    cpdef double get_T(self) except *
    cpdef double get_p(self) except *
    cpdef double get_h(self) except *
    cpdef double get_rho(self) except *
    cpdef double get_s(self) except *
    cpdef double get_u(self) except *
    cpdef double get_visc(self) except *
    cpdef double get_cond(self) except *
    cpdef double get_cp(self) except *
    cpdef double get_cp0(self) except *
    cpdef double get_cv(self) except *
    cpdef double get_MM(self) except *
    cpdef double get_dpdT(self) except *
    cpdef double get_speed_sound(self) except *
    cpdef get_Tsat(self, double Q = *)
    cpdef get_subcooling(self)
    cpdef get_superheat(self)
