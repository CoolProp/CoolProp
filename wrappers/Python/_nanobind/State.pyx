# distutils: language = c++
# cython: language_level=3
"""
Frozen ``State`` compatibility shim for the nanobind-based CoolProp.

PDSim (and other downstream code) ``cimport``/``import`` this module exactly as
they did the legacy Cython ``State``.  Every call is forwarded through the
``CoolProp._capi`` PyCapsule that the nanobind core exports; there is no Cython
``AbstractState`` anywhere.  The shim is link-free -- the C-ABI function table is
fetched from the already-imported core at import time -- so it contains no
CoolProp C++ and does not duplicate the core in the wheel.

Unit convention (the legacy one, pinned by ``dev/state_capsule`` parity test):
  - the kPa/kJ getters: ``get_p`` = p/1000, ``get_h``/``u``/``s``/``cp``/``cp0``/
    ``cv`` = .../1000, ``get_dpdT`` = first_partial_deriv(iP,iT,iDmolar)/1000;
  - ``get_MM`` returns g/mol and ``get_cond`` returns kW/m/K (both 1000x the SI
    value the core reports);
  - ``get_rho``/``get_T``/``get_visc``/``get_speed_sound`` are SI; cached ``p_``
    is kPa, ``T_``/``rho_`` are SI;
  - ``Props()`` and everything reached via ``.pAS`` are SI.
"""
from cpython.pycapsule cimport PyCapsule_GetPointer
from libcpp.vector cimport vector

cdef extern from "CoolProp/detail/state_capi.h":
    ctypedef struct CoolProp_StateCAPI:
        void* (*make)(const char*, const char*)
        void (*destroy)(void*)
        void (*update)(void*, long, double, double)
        double (*keyed_output)(void*, long)
        double (*first_partial_deriv)(void*, long, long, long)
        const char* (*last_error)()
        void (*set_mole_fractions)(void*, const double*, long)

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
cdef long _iPhase = int(_CP.iPhase)
cdef long _iP_critical = int(_CP.iP_critical)
cdef long _iP_triple = int(_CP.iP_triple)
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


cdef inline double _keyed_output(void* handle, long key) except *:
    # One keyed output (SI), with the C++ error re-raised as Python.
    cdef double v = _C.keyed_output(handle, key)
    _raise_if_error()
    return v


cdef inline void _update(void* handle, long input_pair, double value1, double value2) except *:
    # Update the handle's state from an (input_pair, value1, value2) triple.
    _C.update(handle, input_pair, value1, value2)
    _raise_if_error()


cdef inline double _first_partial_deriv(void* handle, long Of, long Wrt, long Constant) except *:
    # First partial derivative d(Of)/d(Wrt) at constant `Constant` (SI).
    cdef double v = _C.first_partial_deriv(handle, Of, Wrt, Constant)
    _raise_if_error()
    return v


cdef inline void _set_mole_fractions(void* handle, list fracs) except *:
    # Set the mixture mole fractions on the handle (used for bracketed fluids).
    cdef vector[double] z
    cdef double f
    for f in fracs:
        z.push_back(f)
    _C.set_mole_fractions(handle, z.data(), <long>z.size())
    _raise_if_error()


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
    """Non-owning view of the same ``AbstractState`` handle as its owning ``State``.

    Reached as ``State.pAS``; mirrors the subset of the real ``AbstractState``
    that downstream code uses, and (unlike the ``State`` getters) returns plain
    SI units.  It does not own the handle, so it has no destructor.
    """
    cpdef double keyed_output(self, long key) except *:
        """Value of CoolProp parameter ``key`` for the current state [SI]."""
        return _keyed_output(self.handle, key)
    cpdef double rhomass(self) except *:
        """Mass density [kg/m^3, SI]."""
        return _keyed_output(self.handle, _iDmass)
    cpdef double cpmass(self) except *:
        """Mass constant-pressure specific heat [J/kg/K, SI]."""
        return _keyed_output(self.handle, _iCpmass)
    cpdef double cvmass(self) except *:
        """Mass constant-volume specific heat [J/kg/K, SI]."""
        return _keyed_output(self.handle, _iCvmass)
    cpdef double T(self) except *:
        """Temperature [K, SI]."""
        return _keyed_output(self.handle, _iT)
    cpdef double p(self) except *:
        """Pressure [Pa, SI]."""
        return _keyed_output(self.handle, _iP)
    cpdef update(self, long input_pair, double value1, double value2):
        """Update the underlying state from an SI ``input_pair`` and two values."""
        _update(self.handle, input_pair, value1, value2)
    cpdef double first_partial_deriv(self, long Of, long Wrt, long Constant) except *:
        """First partial derivative d(Of)/d(Wrt) at constant ``Constant`` [SI]."""
        return _first_partial_deriv(self.handle, Of, Wrt, Constant)
    cpdef fluid_names(self):
        """List of the constituent fluid names."""
        return self.fluids.split('&')


cdef class State:
    """Drop-in replacement for the legacy Cython ``CoolProp.State.State``.

    Holds an ``AbstractState`` (created for backend ``"HEOS"``) via the C-ABI
    capsule and exposes the legacy kPa/kJ getter surface plus the ``.pAS`` view
    (SI).  Constructed from a fluid name and an input dict, e.g.
    ``State("Water", {"T": 300, "D": 1000})`` or ``State("R134a", {"T": 300, "P": 500})``
    (pressures in kPa).  See the module docstring for the full unit convention.

    Cached attributes (refreshed on every update): ``T_`` [K], ``p_`` [kPa],
    ``rho_`` [kg/m^3]; ``Fluid`` (bytes), ``phase`` (bytes); ``pAS`` -- the
    :class:`_AbstractStateView` over the same handle.
    """
    def __cinit__(self, object Fluid, params=None, object phase=None, object backend=None):
        # Widened to the legacy State.__init__ surface (bd CoolProp-r9sq.26):
        # optional ``params`` (None -> no state update, used by get_Tsat/copy),
        # an explicit ``backend``, a ``BACKEND::Fluid`` fluid string, and the
        # (legacy-deprecated) ``phase`` flag -- recorded on ``.phase`` but, since
        # the capsule exposes no specify_phase, not imposed on the state.
        cdef str _fl = Fluid if isinstance(Fluid, str) else Fluid.decode()
        cdef str _backend
        self.handle = NULL
        self.pAS = None
        self.phase = b"" if phase is None else (phase.encode() if isinstance(phase, str) else phase)
        if _fl == 'none':  # legacy no-op construction (e.g. the bare State get_Tsat builds)
            self._fluids = 'none'
            self.Fluid = b'none'
            return
        if '::' in _fl:
            _backend, _fl = _fl.split('::')
        elif backend is None:
            _backend = "HEOS"
        elif isinstance(backend, str):
            _backend = backend
        else:
            _backend = backend.decode()
        self.set_Fluid(_fl, _backend)
        if params is not None:
            self.update(params)

    def __dealloc__(self):
        if self.handle != NULL:
            _C.destroy(self.handle)

    cpdef set_Fluid(self, Fluid, backend):
        """(Re)create the underlying ``AbstractState`` for ``Fluid`` on ``backend``.

        Plain names, ``&``-joined mixtures, and bracketed mole-fraction mixtures
        (e.g. ``R32[0.5]&R134a[0.5]``) all work -- the bracket fractions are parsed
        out and set on the handle, exactly like the legacy ``State.set_Fluid``.
        """
        cdef str _fl = Fluid if isinstance(Fluid, str) else Fluid.decode()
        cdef str _backend
        if backend is None:
            _backend = "HEOS"
        elif isinstance(backend, str):
            _backend = backend
        else:
            _backend = backend.decode()
        # Parse bracketed mole fractions, e.g. "R32[0.5]&R134a[0.5]" -> names
        # "R32&R134a" + fracs [0.5, 0.5] (mirrors the legacy set_Fluid).
        cdef list names = []
        cdef list fracs = []
        cdef bint set_fractions = False
        if '[' in _fl and ']' in _fl:
            for pair in _fl.split('&'):
                name, frac = pair.split('[')
                names.append(name)
                fracs.append(float(frac.strip(']')))
            _fl = '&'.join(names)
            set_fractions = True
        if self.handle != NULL:
            _C.destroy(self.handle)
            self.handle = NULL
        self.handle = _C.make(_backend.encode(), _fl.encode())
        if self.handle == NULL:
            _raise_if_error()
            raise ValueError("State: could not create AbstractState for " + _fl)
        if set_fractions:
            _set_mole_fractions(self.handle, fracs)
        self._fluids = _fl
        self.Fluid = _fl.encode()
        if self.pAS is None:
            self.pAS = _AbstractStateView.__new__(_AbstractStateView)
        self.pAS.handle = self.handle
        self.pAS.fluids = self._fluids

    cdef _refresh(self):
        # Refresh the cached T_/p_/rho_ after a state change.
        self.T_ = _keyed_output(self.handle, _iT)             # K
        self.p_ = _keyed_output(self.handle, _iP) / 1000.0    # kPa
        self.rho_ = _keyed_output(self.handle, _iDmass)       # kg/m^3 (SI)

    cpdef update(self, dict params):
        """Update the state from a legacy input dict.

        ``params`` holds exactly two of T/D/P/Q/H, e.g. ``{'T': K, 'D': kg/m^3}``,
        ``{'T': K, 'P': kPa}``, ``{'T': K, 'Q': -}``, ``{'P': kPa, 'Q': -}`` or
        ``{'P': kPa, 'H': kJ/kg}`` (P/H given in the legacy kPa/kJ units).
        """
        cdef long pair
        cdef double v1, v2
        pair, v1, v2 = _resolve_pair(params)
        _update(self.handle, pair, v1, v2)
        self._refresh()

    cpdef update_Trho(self, double T, double rho):
        """Update from temperature ``T`` [K] and mass density ``rho`` [kg/m^3]."""
        _update(self.handle, _DmassT, rho, T)
        self._refresh()

    cpdef update_ph(self, double p, double h):
        """Update from pressure ``p`` [kPa] and mass enthalpy ``h`` [kJ/kg]."""
        _update(self.handle, _HmassP, h * 1000.0, p * 1000.0)
        self._refresh()

    cpdef State copy(self):
        """Return an independent ``State`` for the same fluid at this (T, rho)."""
        return State(self._fluids, {'T': self.T_, 'D': self.rho_})

    cpdef double Props(self, long key) except *:
        """Raw keyed output for CoolProp parameter ``key`` [SI]."""
        if key < 0:  # legacy parity: ValueError on an invalid (negative) key
            raise ValueError('Your output is invalid')
        return _keyed_output(self.handle, key)

    cpdef long Phase(self) except *:
        """Integer phase flag for the current state (an ``iphase_*`` constant)."""
        return <long>_keyed_output(self.handle, _iPhase)

    cpdef get_Tsat(self, double Q=1):
        """Saturation temperature [K] at the current pressure (``Q=1`` dew, ``Q=0``
        bubble), or ``None`` if the pressure is outside the two-phase range."""
        cdef State state = State(self._fluids, None)
        cdef double pc, pt
        try:
            pc = state.Props(_iP_critical)   # SI Pa
        except ValueError:
            pc = -1
        try:
            pt = state.Props(_iP_triple)     # SI Pa
        except ValueError:
            pt = -1
        # self.p_ is kPa; pc/pt are Pa -> compare against 0.001*p [kPa].
        if pc > 0:
            if self.p_ > 0.001 * pc or (pt > 0 and self.p_ < 0.001 * pt):
                return None
        state.update({'P': self.p_, 'Q': Q})
        return state.T

    cpdef get_superheat(self):
        """Superheat above the dew temperature [K], or ``None`` outside two-phase."""
        Tsat = self.get_Tsat(1)
        return None if Tsat is None else self.T_ - Tsat

    cpdef get_subcooling(self):
        """Subcooling below the bubble temperature [K], or ``None`` outside two-phase."""
        Tsat = self.get_Tsat(0)
        return None if Tsat is None else Tsat - self.T_

    cpdef double get_T(self) except *:
        """Temperature [K]."""
        return _keyed_output(self.handle, _iT)
    cpdef double get_p(self) except *:
        """Pressure [kPa]."""
        return _keyed_output(self.handle, _iP) / 1000.0
    cpdef double get_h(self) except *:
        """Mass enthalpy [kJ/kg]."""
        return _keyed_output(self.handle, _iHmass) / 1000.0
    cpdef double get_rho(self) except *:
        """Mass density [kg/m^3]."""
        return _keyed_output(self.handle, _iDmass)
    cpdef double get_s(self) except *:
        """Mass entropy [kJ/kg/K]."""
        return _keyed_output(self.handle, _iSmass) / 1000.0
    cpdef double get_u(self) except *:
        """Mass internal energy [kJ/kg]."""
        return _keyed_output(self.handle, _iUmass) / 1000.0
    cpdef double get_cp(self) except *:
        """Mass constant-pressure specific heat [kJ/kg/K]."""
        return _keyed_output(self.handle, _iCpmass) / 1000.0
    cpdef double get_cp0(self) except *:
        """Mass ideal-gas specific heat [kJ/kg/K]."""
        return _keyed_output(self.handle, _iCp0mass) / 1000.0
    cpdef double get_cv(self) except *:
        """Mass constant-volume specific heat [kJ/kg/K]."""
        return _keyed_output(self.handle, _iCvmass) / 1000.0
    cpdef double get_MM(self) except *:
        """Molar mass [g/mol] (legacy convention; the core reports kg/mol)."""
        return _keyed_output(self.handle, _imolar_mass) * 1000.0
    cpdef double get_dpdT(self) except *:
        """(dp/dT) at constant density [kPa/K]."""
        return _first_partial_deriv(self.handle, _iP, _iT, _iDmolar) / 1000.0
    cpdef double get_visc(self) except *:
        """Viscosity [Pa.s]."""
        return _keyed_output(self.handle, _iviscosity)
    cpdef double get_cond(self) except *:
        """Thermal conductivity [kW/m/K] (legacy convention; the core reports W/m/K)."""
        return _keyed_output(self.handle, _iconductivity) / 1000.0
    cpdef double get_speed_sound(self) except *:
        """Speed of sound [m/s]."""
        return _keyed_output(self.handle, _ispeed_sound)
    cpdef double get_Q(self) except *:
        """Vapor quality [-]."""
        return _keyed_output(self.handle, _iQ)

    # Python-level property accessors (the legacy State exposes these alongside
    # the get_* methods; PDSim's .py code uses state.p, state.T, ...).
    property Q:
        """Vapor quality [-]."""
        def __get__(self): return self.get_Q()
    property MM:
        """Molar mass [g/mol]."""
        def __get__(self): return self.get_MM()
    property rho:
        """Mass density [kg/m^3]."""
        def __get__(self): return self.get_rho()
    property p:
        """Pressure [kPa]."""
        def __get__(self): return self.get_p()
    property T:
        """Temperature [K]."""
        def __get__(self): return self.get_T()
    property h:
        """Mass enthalpy [kJ/kg]."""
        def __get__(self): return self.get_h()
    property u:
        """Mass internal energy [kJ/kg]."""
        def __get__(self): return self.get_u()
    property s:
        """Mass entropy [kJ/kg/K]."""
        def __get__(self): return self.get_s()
    property cp0:
        """Mass ideal-gas specific heat [kJ/kg/K]."""
        def __get__(self): return self.get_cp0()
    property cp:
        """Mass constant-pressure specific heat [kJ/kg/K]."""
        def __get__(self): return self.get_cp()
    property cv:
        """Mass constant-volume specific heat [kJ/kg/K]."""
        def __get__(self): return self.get_cv()
    property visc:
        """Viscosity [Pa.s]."""
        def __get__(self): return self.get_visc()
    property k:
        """Thermal conductivity [kW/m/K]."""
        def __get__(self): return self.get_cond()
    property dpdT:
        """(dp/dT) at constant density [kPa/K]."""
        def __get__(self): return self.get_dpdT()
    property speed_sound:
        """Speed of sound [m/s]."""
        def __get__(self): return self.get_speed_sound()
    property Tsat:
        """Saturation (dew) temperature at the current pressure [K]."""
        def __get__(self): return self.get_Tsat(1.0)
    property superheat:
        """Superheat above the dew temperature [K] (None outside two-phase)."""
        def __get__(self): return self.get_superheat()
    property subcooling:
        """Subcooling below the bubble temperature [K] (None outside two-phase)."""
        def __get__(self): return self.get_subcooling()
    property Prandtl:
        """Prandtl number cp*mu/k [-]."""
        def __get__(self): return self.cp * self.visc / self.k
