from libcpp.string cimport string
import cython
cimport cython

from libcpp.vector cimport vector

include "AbstractState.pxd"
    
cdef extern from "CoolPropTools.h" namespace "CoolProp":
    bint _ValidNumber "ValidNumber"(double)
    
cdef extern from "DataStructures.h" namespace "CoolProp":    
    string _get_parameter_information "CoolProp::get_parameter_information"(int, string) except +
    int _get_parameter_index "CoolProp::get_parameter_index"(string) except +
    
cdef extern from "CoolProp.h" namespace "CoolProp":
    double _Props1SI "CoolProp::Props1SI"(string Ref, string Output)
    double _PropsSI "CoolProp::PropsSI"(string Output, string Name1, double Prop1, string Name2, double Prop2, string FluidName) 
    vector[double] _PropsSI "CoolProp::PropsSI"(string Output, string Name1, vector[double] Prop1, string Name2, vector[double] Prop2, string FluidName, vector[double] fractions)
    vector[double] _PropsSII "CoolProp::PropsSI"(string Output, string Name1, vector[double] Prop1, string Name2, vector[double] Prop2, string FluidName)
    string _get_global_param_string "CoolProp::get_global_param_string"(string ParamName) except +
    
#     double _Props "CoolProp::Props"(string Output, string Name1, double Prop1, string Name2, double Prop2, string Ref)
#     double _Props1 "CoolProp::Props1"(string Ref, string Output)
#     string _get_fluid_param_string "CoolProp::get_fluid_param_string"(string ParamName, string FluidName)
    
    #long _get_parameter_index "CoolProp::get_parameter_index" (string param)
    int _get_debug_level "CoolProp::get_debug_level"()
    void _set_debug_level "CoolProp::set_debug_level"(int level)
    
    # Convenience functions
#     int _IsFluidType "IsFluidType"(char* Ref, char* Type)
#     string _get_BibTeXKey "CoolProp::get_BibTeXKey"(string Ref, string key)
#     long _get_Fluid_index "CoolProp::get_Fluid_index" (string Fluid)
#     double _IProps "CoolProp::IProps"(long Output, long Name1, double Prop1, long Name2, double Prop2, long Ref)
 
cdef extern from "HumidAirProp.h":
    double _HAPropsSI "HumidAir::HAPropsSI"(string OutputName, string Input1Name, double Input1, string Input2Name, double Input2, string Input3Name, double Input3)
#     double _HAProps_Aux "HAProps_Aux"(const char* Name,double T, double p, double W, char *units)
#     double _cair_sat "cair_sat"(double T)
       
cdef class State:
    cdef AbstractState pAS
    cdef readonly string Fluid, phase
    cdef int iFluid,iParam1,iParam2,iOutput
    cdef double T_, rho_, p_
    cdef readonly bint is_CPFluid
    
    cpdef set_Fluid(self, string Fluid, string backend)
#     cpdef speed_test(self, int N)
#     cpdef update(self, dict params)
#     cpdef update_ph(self, double p, double h)
    cpdef update_Trho(self, double T, double rho)
#     cpdef State copy(self)
#     cpdef double Props(self, long iOutput) except *
#     cpdef long Phase(self) except *
    
#     cpdef double get_Q(self) except *
#     cpdef double get_T(self) except *
#     cpdef double get_p(self) except *
#     cpdef double get_h(self) except *
#     cpdef double get_rho(self) except *
#     cpdef double get_s(self) except *
#     cpdef double get_u(self) except *
#     cpdef double get_visc(self) except *
#     cpdef double get_cond(self) except *
#     cpdef double get_cp(self) except *
#     cpdef double get_cp0(self) except *
#     cpdef double get_cv(self) except *
#     cpdef double get_MM(self) except *
#     cpdef double get_dpdT(self) except *
#     cpdef double get_speed_sound(self) except *
#     cpdef get_Tsat(self, double Q = *)
#     cpdef get_subcooling(self)
#     cpdef get_superheat(self)cdef class State:
#     cdef CoolPropStateClassSI CPS
#     cdef readonly string Fluid, phase
#     cdef int iFluid,iParam1,iParam2,iOutput
#     cdef double T_, rho_, p_
#     cdef readonly bint is_CPFluid
#     
#     cpdef set_Fluid(self, string_like Fluid)
    cpdef speed_test(self, int N)
#     cpdef update(self, dict params)
#     cpdef update_ph(self, double p, double h)
#     cpdef update_Trho(self, double T, double rho)
#     cpdef State copy(self)
#     cpdef double Props(self, long iOutput) except *
#     cpdef long Phase(self) except *
#     
#     cpdef double get_Q(self) except *
#     cpdef double get_T(self) except *
#     cpdef double get_p(self) except *
#     cpdef double get_h(self) except *
#     cpdef double get_rho(self) except *
#     cpdef double get_s(self) except *
#     cpdef double get_u(self) except *
#     cpdef double get_visc(self) except *
#     cpdef double get_cond(self) except *
#     cpdef double get_cp(self) except *
#     cpdef double get_cp0(self) except *
#     cpdef double get_cv(self) except *
#     cpdef double get_MM(self) except *
#     cpdef double get_dpdT(self) except *
#     cpdef double get_speed_sound(self) except *
#     cpdef get_Tsat(self, double Q = *)
#     cpdef get_subcooling(self)
#     cpdef get_superheat(self)
