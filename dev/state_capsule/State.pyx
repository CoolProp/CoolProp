# distutils: language = c++
# cython: language_level=3
"""
Frozen `State` compat shim (full surface — increment 2).

PDSim cimports this instead of CoolProp's old Cython `State`.  Every call is
forwarded through the `CoolProp._capi` PyCapsule the (nanobind) core exports;
there is no Cython AbstractState anywhere.  Link-free: the function table is
fetched from the already-imported core at import time.

Unit convention is the legacy one pinned by dev/pdsim_cimport_contract:
  - State getters return kPa / kJ (get_p = p()/1000, get_h/u/s/cp/cp0/cv = …/1000,
    get_dpdT = first_partial_deriv(iP,iT,iDmolar)/1000); get_rho / get_speed_sound
    / get_T are SI; cached p_ is kPa, rho_/T_ are SI.
  - Props() and everything via .pAS are SI.
"""
from cpython.pycapsule cimport PyCapsule_GetPointer

cdef extern from "CoolProp/detail/state_capi.h":
    ctypedef struct CoolProp_StateCAPI:
        void* (*make)(const char*, const char*)
        void (*destroy)(void*)
        void (*update)(void*, long, double, double)
        double (*keyed_output)(void*, long)
        double (*first_partial_deriv)(void*, long, long, long)
        const char* (*last_error)()

import CoolProp as _CP  # the nanobind core: exports `_capi` + the int-convertible enums
if not hasattr(_CP, "_capi"):
    # package layout: the core is the CoolProp.CoolProp submodule, re-exported
    # by CoolProp/__init__.py (which `from .CoolProp import *` skips _capi).
    from CoolProp import CoolProp as _CP
cdef CoolProp_StateCAPI* _C = <CoolProp_StateCAPI*>PyCapsule_GetPointer(_CP._capi, "CoolProp._capi")

# Canonical key / input-pair indices, read once from the core's enums.
cdef long _iT = int(_CP.iT)
cdef long _iP = int(_CP.iP)
cdef long _iHmass = int(_CP.iHmass)
cdef long _iSmass = int(_CP.iSmass)
cdef long _iUmass = int(_CP.iUmass)
cdef long _iDmass = int(_CP.iDmass)
cdef long _iDmolar = int(_CP.iDmolar)
cdef long _iCpmass = int(_CP.iCpmass)
cdef long _iCp0mass = int(_CP.iCp0mass)
cdef long _iCvmass = int(_CP.iCvmass)
cdef long _imolar_mass = int(_CP.imolar_mass)
cdef long _iviscosity = int(_CP.iviscosity)
cdef long _iconductivity = int(_CP.iconductivity)
cdef long _ispeed_sound = int(_CP.ispeed_sound)
cdef long _iQ = int(_CP.iQ)
cdef long _DmassT = int(_CP.DmassT_INPUTS)
cdef long _QT = int(_CP.QT_INPUTS)
cdef long _PT = int(_CP.PT_INPUTS)
cdef long _PQ = int(_CP.PQ_INPUTS)
cdef long _HmassP = int(_CP.HmassP_INPUTS)


cdef inline void _raise_if_error() except *:
    # CoolProp threw a C++ exception across the capsule; re-raise it as Python
    # (PDSim's ODE solver catches these and rejects the step -- like the old State).
    cdef const char* e = _C.last_error()
    if e is not NULL:
        raise ValueError("CoolProp: " + e.decode("utf-8", "replace"))


cdef inline double _ko(void* h, long key) except *:
    cdef double v = _C.keyed_output(h, key)
    _raise_if_error()
    return v


cdef inline void _upd(void* h, long pair, double v1, double v2) except *:
    _C.update(h, pair, v1, v2)
    _raise_if_error()


cdef inline double _fpd(void* h, long Of, long Wrt, long Constant) except *:
    cdef double v = _C.first_partial_deriv(h, Of, Wrt, Constant)
    _raise_if_error()
    return v


cdef tuple _resolve_pair(dict params):
    # Map a legacy State dict (kPa/kJ values for P/H/S) to an SI update triple.
    cdef set keys = set(params.keys())
    if keys == {'T', 'D'}:
        return (_DmassT, params['D'], params['T'])
    elif keys == {'T', 'Q'}:
        return (_QT, params['Q'], params['T'])
    elif keys == {'T', 'P'}:
        return (_PT, params['P'] * 1000.0, params['T'])
    elif keys == {'P', 'Q'}:
        return (_PQ, params['P'] * 1000.0, params['Q'])
    elif keys == {'P', 'H'}:
        return (_HmassP, params['H'] * 1000.0, params['P'] * 1000.0)
    raise ValueError("State shim: unsupported input pair " + repr(sorted(keys)))


cdef class _AbstractStateView:
    # Non-owning view of the same handle: the AbstractState methods PDSim
    # reaches through `State.pAS` (SI units, like the real AbstractState).
    cpdef double keyed_output(self, long key) except *:
        return _ko(self.handle,key)
    cpdef double rhomass(self) except *:
        return _ko(self.handle,_iDmass)
    cpdef double cpmass(self) except *:
        return _ko(self.handle,_iCpmass)
    cpdef double cvmass(self) except *:
        return _ko(self.handle,_iCvmass)
    cpdef double T(self) except *:
        return _ko(self.handle,_iT)
    cpdef double p(self) except *:
        return _ko(self.handle,_iP)
    cpdef update(self, long input_pair, double value1, double value2):
        _upd(self.handle,input_pair, value1, value2)
    cpdef double first_partial_deriv(self, long Of, long Wrt, long Constant) except *:
        return _fpd(self.handle,Of, Wrt, Constant)
    cpdef fluid_names(self):
        return self.fluids.split('&')


cdef class State:
    def __cinit__(self, object Fluid, dict params):
        self._fluids = Fluid if isinstance(Fluid, str) else Fluid.decode()
        self.handle = _C.make(b"HEOS", self._fluids.encode())
        if self.handle == NULL:
            _raise_if_error()
            raise ValueError("State: could not create AbstractState for " + self._fluids)
        self.Fluid = self._fluids.encode()
        self.phase = b""
        self.pAS = _AbstractStateView.__new__(_AbstractStateView)
        self.pAS.handle = self.handle
        self.pAS.fluids = self._fluids
        self.update(params)

    def __dealloc__(self):
        if self.handle != NULL:
            _C.destroy(self.handle)

    cdef _refresh(self):
        self.T_ = _ko(self.handle,_iT)             # K
        self.p_ = _ko(self.handle,_iP) / 1000.0    # kPa
        self.rho_ = _ko(self.handle,_iDmass)       # kg/m^3 (SI)

    cpdef update(self, dict params):
        cdef long pair
        cdef double v1, v2
        pair, v1, v2 = _resolve_pair(params)
        _upd(self.handle,pair, v1, v2)
        self._refresh()

    cpdef update_Trho(self, double T, double rho):
        _upd(self.handle, _DmassT, rho, T)
        self._refresh()

    cpdef update_ph(self, double p, double h):  # p [kPa], h [kJ/kg]
        _upd(self.handle,_HmassP, h * 1000.0, p * 1000.0)
        self._refresh()

    cpdef State copy(self):
        return State(self._fluids, {'T': self.T_, 'D': self.rho_})

    cpdef double Props(self, long key) except *:
        return _ko(self.handle,key)                # SI

    cpdef double get_T(self) except *:
        return _ko(self.handle,_iT)                # K
    cpdef double get_p(self) except *:
        return _ko(self.handle,_iP) / 1000.0       # kPa
    cpdef double get_h(self) except *:
        return _ko(self.handle,_iHmass) / 1000.0   # kJ/kg
    cpdef double get_rho(self) except *:
        return _ko(self.handle,_iDmass)            # kg/m^3 (SI)
    cpdef double get_s(self) except *:
        return _ko(self.handle,_iSmass) / 1000.0
    cpdef double get_u(self) except *:
        return _ko(self.handle,_iUmass) / 1000.0
    cpdef double get_cp(self) except *:
        return _ko(self.handle,_iCpmass) / 1000.0
    cpdef double get_cp0(self) except *:
        return _ko(self.handle,_iCp0mass) / 1000.0
    cpdef double get_cv(self) except *:
        return _ko(self.handle,_iCvmass) / 1000.0
    cpdef double get_MM(self) except *:
        return _ko(self.handle,_imolar_mass) * 1000.0   # g/mol (legacy State convention; core is kg/mol SI)
    cpdef double get_dpdT(self) except *:
        return _fpd(self.handle,_iP, _iT, _iDmolar) / 1000.0
    cpdef double get_visc(self) except *:
        return _ko(self.handle,_iviscosity)
    cpdef double get_cond(self) except *:
        return _ko(self.handle,_iconductivity) / 1000.0   # kW/m/K (legacy State convention; core is W/m/K SI)
    cpdef double get_speed_sound(self) except *:
        return _ko(self.handle,_ispeed_sound)      # m/s (SI)
    cpdef double get_Q(self) except *:
        return _ko(self.handle,_iQ)

    # Python-level property accessors (the legacy State exposes these in
    # addition to the get_* methods; PDSim's .py code uses state.p, state.T, ...)
    property Q:
        def __get__(self): return self.get_Q()
    property MM:
        def __get__(self): return self.get_MM()
    property rho:
        def __get__(self): return self.get_rho()
    property p:
        def __get__(self): return self.get_p()
    property T:
        def __get__(self): return self.get_T()
    property h:
        def __get__(self): return self.get_h()
    property u:
        def __get__(self): return self.get_u()
    property s:
        def __get__(self): return self.get_s()
    property cp0:
        def __get__(self): return self.get_cp0()
    property cp:
        def __get__(self): return self.get_cp()
    property cv:
        def __get__(self): return self.get_cv()
    property visc:
        def __get__(self): return self.get_visc()
    property k:
        def __get__(self): return self.get_cond()
    property dpdT:
        def __get__(self): return self.get_dpdT()
    property speed_sound:
        def __get__(self): return self.get_speed_sound()
