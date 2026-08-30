# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import

import matplotlib.pyplot as plt
import numpy as np
from abc import ABCMeta
import warnings

import CoolProp
from CoolProp import AbstractState
from CoolProp import CoolProp as CP
from CoolProp.CoolProp import PropsSI, extract_backend, extract_fractions, PyCriticalState


def get_critical_point(state):
    crit_state = PyCriticalState()
    crit_state.T = np.nan
    crit_state.p = np.nan
    crit_state.rhomolar = np.nan
    crit_state.rhomolar = np.nan
    crit_state.stable = False
    try:
        crit_state.T = state.T_critical()
        crit_state.p = state.p_critical()
        crit_state.rhomolar = state.rhomolar_critical()
        crit_state.stable = True
    except:
        try:
            for crit_state_tmp in state.all_critical_points():
                if crit_state_tmp.stable and (crit_state_tmp.T > crit_state.T or not np.isfinite(crit_state.T)):
                    crit_state.T = crit_state_tmp.T
                    crit_state.p = crit_state_tmp.p
                    crit_state.rhomolar = crit_state_tmp.rhomolar
                    crit_state.stable = crit_state_tmp.stable
        except:
            raise ValueError("Could not calculate the critical point data.")
    new_state = AbstractState(state.backend_name(), '&'.join(state.fluid_names()))
    masses = state.get_mass_fractions()
    if len(masses) > 1:
        new_state.set_mass_fractions(masses)  # Uses mass fraction to work with incompressibles
        # try: new_state.build_phase_envelope("dummy")
        # except: pass
    msg = ""
    if np.isfinite(crit_state.p) and np.isfinite(crit_state.T):
        try:
            new_state.specify_phase(CoolProp.iphase_critical_point)
            new_state.update(CoolProp.PT_INPUTS, crit_state.p, crit_state.T)
            return new_state
        except Exception as e:
            msg += str(e) + " - "
            pass
        try:
            new_state.update(CoolProp.PT_INPUTS, crit_state.p, crit_state.T)
            return new_state
        except Exception as e:
            msg += str(e) + " - "
            pass
    if np.isfinite(crit_state.rhomolar) and np.isfinite(crit_state.T):
        try:
            new_state.specify_phase(CoolProp.iphase_critical_point)
            new_state.update(CoolProp.DmolarT_INPUTS, crit_state.rhomolar, crit_state.T)
            return new_state
        except Exception as e:
            msg += str(e) + " - "
            pass
        try:
            new_state.update(CoolProp.DmolarT_INPUTS, crit_state.rhomolar, crit_state.T)
            return new_state
        except Exception as e:
            msg += str(e) + " - "
            pass
    raise ValueError("Could not calculate the critical point data. " + msg)


def interpolate_values_1d(x, y, x_points=None, kind='linear'):
    try:
        from scipy.interpolate.interpolate import interp1d
        if x_points is None:
            return interp1d(x, y, kind=kind)(x[np.isfinite(x)])
        else:
            return interp1d(x, y, kind=kind)(x_points)
    except ImportError:
        if kind != 'linear':
            warnings.warn(
              "You requested a non-linear interpolation, but SciPy is not available. Falling back to linear interpolation.",
              UserWarning)
        if x_points is None:
            return np.interp((x[np.isfinite(x)]), x, y)
        else:
            return np.interp(x_points, x, y)


def is_string(in_obj):
    return isinstance(in_obj, str)
    # except:
    #    return False


def process_fluid_state(fluid_ref, fractions='mole'):
    """Check input for state object or fluid string

    Parameters
    ----------
        fluid_ref : str, CoolProp.AbstractState
        fractions : str, switch to set mass, volu or mole fractions

    Returns
    -------
        CoolProp.AbstractState
    """
    # Process the fluid and set self._state
    if is_string(fluid_ref):
        backend, fluids = extract_backend(fluid_ref)
        fluids, fractions = extract_fractions(fluids)
        state = AbstractState(backend, '&'.join(fluids))
        if len(fluids) > 1 and len(fluids) == len(fractions):
            if fractions == 'mass': state.set_mass_fractions(fractions)
            elif fractions == 'volu': state.set_volu_fractions(fractions)
            else: state.set_mole_fractions(fractions)
        return state
    elif isinstance(fluid_ref, AbstractState):
        return fluid_ref
    raise TypeError("Invalid fluid_ref input, expected a string or an abstract state instance.")


def _get_index(prop):
    if is_string(prop):
        return CP.get_parameter_index(prop)
    elif isinstance(prop, int):
        return prop
    else:
        raise ValueError("Invalid input, expected a string or an int, not {0:s}.".format(str(prop)))


class BaseQuantity(object):
    """A very basic property that can convert an input to and from a
    given unit system, note that the conversion from SI units starts
    with a multiplication. If you need to remove an offset, use the
    off_SI property.
    Examples with temperature:
    celsius = BaseQuantity(add_SI=-273.15)
    fahrenheit = BaseQuantity(add_SI=32.0, mul_SI=1.8, off_SI=-273.15)
    Examples with pressure:
    bar = BaseQuantity(mul_SI=1e-5)
    psi = BaseQuantity(mul_SI=0.000145037738)
    """

    def __init__(self, add_SI=0.0, mul_SI=1.0, off_SI=0.0):
        self._add_SI = add_SI
        self._mul_SI = mul_SI
        self._off_SI = off_SI

    @property
    def add_SI(self): return self._add_SI

    @add_SI.setter
    def add_SI(self, value): self._add_SI = value

    @property
    def mul_SI(self): return self._mul_SI

    @mul_SI.setter
    def mul_SI(self, value): self._mul_SI = value

    @property
    def off_SI(self): return self._off_SI

    @off_SI.setter
    def off_SI(self, value): self._off_SI = value

    def from_SI(self, value): return ((value + self.off_SI) * self.mul_SI) + self.add_SI

    def to_SI(self, value): return (value - self.add_SI) / self.mul_SI - self.off_SI


class BaseDimension(BaseQuantity):
    """A dimension is a class that extends the BaseQuantity and adds a label, a symbol and a unit label"""

    def __init__(self, add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='', symbol='', unit=''):
        self._label = label
        self._symbol = symbol
        self._unit = unit
        super(BaseDimension, self).__init__(add_SI=add_SI, mul_SI=mul_SI, off_SI=off_SI)

    @property
    def label(self): return self._label

    @label.setter
    def label(self, value): self._label = value

    @property
    def symbol(self): return self._symbol

    @symbol.setter
    def symbol(self, value): self._symbol = value

    @property
    def unit(self): return self._unit

    @unit.setter
    def unit(self, value): self._unit = value


class PropertyDict(object, metaclass=ABCMeta):
    """A collection of dimensions for all the required quantities"""

    def __init__(self):
        self._D = None
        self._H = None
        self._P = None
        self._S = None
        self._T = None
        self._U = None
        self._Q = None

    @property
    def D(self): return self._D

    @D.setter
    def D(self, value): self._D = value

    @property
    def H(self): return self._H

    @H.setter
    def H(self, value): self._H = value

    @property
    def P(self): return self._P

    @P.setter
    def P(self, value): self._P = value

    @property
    def S(self): return self._S

    @S.setter
    def S(self, value): self._S = value

    @property
    def T(self): return self._T

    @T.setter
    def T(self, value): self._T = value

    @property
    def U(self): return self._U

    @U.setter
    def U(self, value): self._U = value

    @property
    def Q(self): return self._Q

    @Q.setter
    def Q(self, value): self._Q = value

    @property
    def dimensions(self):
        return {
      CoolProp.iDmass: self._D,
      CoolProp.iHmass: self._H,
      CoolProp.iP: self._P,
      CoolProp.iSmass: self._S,
      CoolProp.iT: self._T,
      CoolProp.iUmass: self._U,
      CoolProp.iQ: self._Q
    }

    def __getitem__(self, index):
        """Allow for property access via square brackets"""
        idx = _get_index(index)
        if idx == CoolProp.iDmass: return self.D
        elif idx == CoolProp.iHmass: return self.H
        elif idx == CoolProp.iP: return self.P
        elif idx == CoolProp.iSmass: return self.S
        elif idx == CoolProp.iT: return self.T
        elif idx == CoolProp.iUmass: return self.U
        elif idx == CoolProp.iQ: return self.Q
        else: raise IndexError("Unknown index \"{0:s}\".".format(str(index)))

    def __setitem__(self, index, value):
        """Allow for property access via square brackets"""
        idx = _get_index(index)
        if idx == CoolProp.iDmass: self.D = value
        elif idx == CoolProp.iHmass: self.H = value
        elif idx == CoolProp.iP: self.P = value
        elif idx == CoolProp.iSmass: self.S = value
        elif idx == CoolProp.iT: self.T = value
        elif idx == CoolProp.iUmass: self.U = value
        elif idx == CoolProp.iQ: self.Q = value
        else: raise IndexError("Unknown index \"{0:s}\".".format(str(index)))


class SIunits(PropertyDict):
    def __init__(self):
        super(SIunits, self).__init__()
        self._D = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Density', symbol=u'd', unit=u'kg/m3')
        self._H = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Enthalpy', symbol=u'h', unit=u'J/kg')
        self._P = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Pressure', symbol=u'p', unit=u'Pa')
        self._S = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Entropy', symbol=u's', unit=u'J/kg/K')
        self._T = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Temperature', symbol=u'T', unit=u'K')
        self._U = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Specific Internal Energy', symbol=u'u', unit=u'J/kg')
        self._Q = BaseDimension(add_SI=0.0, mul_SI=1.0, off_SI=0.0, label='Vapour Quality', symbol=u'x', unit=u'')


class KSIunits(SIunits):
    def __init__(self):
        super(KSIunits, self).__init__()
        self.H.mul_SI = 1e-3
        self.H.unit = u'kJ/kg'
        self.P.mul_SI = 1e-3
        self.P.unit = u'kPa'
        self.S.mul_SI = 1e-3
        self.S.unit = u'kJ/kg/K'
        self.U.mul_SI = 1e-3
        self.U.unit = u'kJ/kg'


class EURunits(KSIunits):
    def __init__(self):
        super(EURunits, self).__init__()
        self.P.mul_SI = 1e-5
        self.P.unit = u'bar'
        self.T.add_SI = -273.15
        self.T.unit = u'deg C'


class Base2DObject(object, metaclass=ABCMeta):
    """A container for shared settings and constants for the
    isolines and the property plots."""

    # A list of supported plot
    TS = CoolProp.iT * 10 + CoolProp.iSmass
    PH = CoolProp.iP * 10 + CoolProp.iHmass
    HS = CoolProp.iHmass * 10 + CoolProp.iSmass
    PS = CoolProp.iP * 10 + CoolProp.iSmass
    PD = CoolProp.iP * 10 + CoolProp.iDmass
    TD = CoolProp.iT * 10 + CoolProp.iDmass
    PT = CoolProp.iP * 10 + CoolProp.iT
    PU = CoolProp.iP * 10 + CoolProp.iUmass

    PLOTS = {
      'TS': TS,
      'PH': PH,
      'HS': HS,
      'PS': PS,
      'PD': PD,
      'TD': TD,
      'PT': PT,
    }

    PLOTS_INV = {v: k for k, v in PLOTS.items()}

#     # A list of supported plot
#     @property
#     def TS(self): return type(self).TS
#     @property
#     def PH(self): return CoolProp.iP*10     + CoolProp.iHmass
#     @property
#     def HS(self): return CoolProp.iHmass*10 + CoolProp.iSmass
#     @property
#     def PS(self): return CoolProp.iP*10     + CoolProp.iSmass
#     @property
#     def PD(self): return CoolProp.iP*10     + CoolProp.iDmass
#     @property
#     def TD(self): return CoolProp.iT*10     + CoolProp.iDmass
#     @property
#     def PT(self): return CoolProp.iP*10     + CoolProp.iT
#     @property
#     def PU(self): return CoolProp.iP*10     + CoolProp.iUmass

    def __init__(self, x_type, y_type, state=None, small=None):
        self._x_index = _get_index(x_type)
        self._y_index = _get_index(y_type)
        self._critical_state = None
        if small is not None: self._small = small
        else: self._small = 1e-7
        if state is not None: self.state = state
        else: self._state = None

    # A list of supported plot
    @property
    def x_index(self): return self._x_index

    @property
    def y_index(self): return self._y_index

    @property
    def critical_state(self):
        if self._critical_state is None and self._state is not None:
            self._critical_state = get_critical_point(self._state)
        return self._critical_state

    @property
    def state(self): return self._state

    @state.setter
    def state(self, value):
        self._state = process_fluid_state(value)
        # try: self._state.build_phase_envelope("dummy")
        # except: pass
        self._critical_state = None
        #self._T_small = self._state.trivial_keyed_output(CoolProp.iT_critical)*self._small
        #self._P_small = self._state.trivial_keyed_output(CoolProp.iP_critical)*self._small
        self._T_small = self.critical_state.keyed_output(CoolProp.iT) * self._small
        self._P_small = self.critical_state.keyed_output(CoolProp.iP) * self._small

    def _get_sat_bounds(self, kind, smin=None, smax=None):
        """Generates limits for the saturation line in either T or p determined
        by 'kind'. If smin or smax are provided, values will be checked
        against the allowable range for the EOS and a warning might be
        generated. Returns a tuple containing (xmin, xmax)"""

        # TODO: REFPROP backend does not have ptriple.
        T_triple = self._state.trivial_keyed_output(CoolProp.iT_triple)
        try:
            T_min = self._state.trivial_keyed_output(CoolProp.iT_min)
        except:
            T_min = T_triple
        self._state.update(CoolProp.QT_INPUTS, 0, max([T_triple, T_min]) + self._T_small)
        kind = _get_index(kind)
        if kind == CoolProp.iP:
            fluid_min = self._state.keyed_output(CoolProp.iP) + self._P_small
            fluid_max = self.critical_state.keyed_output(CoolProp.iP) - self._P_small
        elif kind == CoolProp.iT:
            fluid_min = self._state.keyed_output(CoolProp.iT) + self._T_small
            fluid_max = self.critical_state.keyed_output(CoolProp.iT) - self._T_small
        else:
            raise ValueError("Saturation boundaries have to be defined in T or P, but not in {0:s}".format(str(kind)))

        if smin is not None:
            if fluid_min < smin < fluid_max:
                sat_min = smin
            else:
                warnings.warn(
                  "Your minimum {0:s} has been ignored, {1:f} is not between {2:f} and {3:f}".format(self.PROPERTIES[kind], smin, fluid_min, fluid_max),
                  UserWarning)
                sat_min = fluid_min
        else:
            sat_min = fluid_min

        if smax is not None:
            if fluid_min < smax < fluid_max:
                sat_max = smax
            else:
                warnings.warn(
                  "Your maximum {0:s} has been ignored, {1:f} is not between {2:f} and {3:f}".format(self.PROPERTIES[kind], smax, fluid_min, fluid_max),
                  UserWarning)
                sat_max = fluid_max
        else:
            sat_max = fluid_max

        return sat_min, sat_max


#: Binary interaction parameters that have to be carried onto a cloned state.
#: Rebuilding from the fluid names alone gets the defaults, so a state the
#: caller customised with set_binary_interaction_double or
#: apply_simple_mixing_rule would otherwise be traced with a different fluid
#: model than it is flashed with.
_BINARY_INTERACTION_KEYS = ('betaT', 'gammaT', 'betaV', 'gammaV', 'Fij', 'kij')


#: A two-phase state whose two phases differ in density by less than this,
#: relatively, is the trivial solution: the VLE solver collapsed both phases
#: onto the same root rather than converging on a saturation state.  Measured
#: over Q sweeps of eight fluids, genuine states bottom out at 3.4e-3 (Air,
#: right at its critical point) while collapsed ones come in at ~1e-6, so this
#: sits about thirty times clear of both.
_TRIVIAL_SOLUTION_RTOL = 1e-4


def _phases_are_distinct(state):
    """Whether a two-phase state really has two phases

    Near the critical point the mixture QT and PQ solvers can return the
    trivial solution, reporting the requested quality and a state whose
    saturated liquid and vapour are the same root.  The result is a plausible
    looking point that is nowhere near the one asked for -- on an R513A Q=0.5
    line it puts the last point 93 J/kg/K along in entropy, which is the jag
    at the top of the dome.
    """
    try:
        rho_liq = state.saturated_liquid_keyed_output(CoolProp.iDmolar)
        rho_vap = state.saturated_vapor_keyed_output(CoolProp.iDmolar)
    except Exception:
        return True         # nothing to check against; the caller's other tests apply
    if not (np.isfinite(rho_liq) and np.isfinite(rho_vap)) or rho_liq <= 0.0 or rho_vap <= 0.0:
        return False
    return (rho_liq - rho_vap) > _TRIVIAL_SOLUTION_RTOL * max(rho_liq, rho_vap)


def _clone_state(state):
    """Return a fresh AbstractState with the same backend, fluids and model"""
    fluids = state.fluid_names()
    new_state = AbstractState(state.backend_name(), '&'.join(fluids))
    if len(fluids) > 1:
        new_state.set_mole_fractions(state.get_mole_fractions())
        for i in range(len(fluids)):
            for j in range(i + 1, len(fluids)):
                for key in _BINARY_INTERACTION_KEYS:
                    try:
                        value = state.get_binary_interaction_double(i, j, key)
                    except Exception:
                        continue
                    new_state.set_binary_interaction_double(i, j, key, value)
    return new_state


class IsoLineTracer(object):
    """Walk along an isoline, warm-starting every point from the previous one.

    :func:`IsoLine.calc_range` evaluates an isoline at consecutive points, so
    the solution at one point is an excellent initial guess for the next one.
    That turns each point into a much cheaper problem than the cold
    two-dimensional flash it replaces:

    * A single-phase point is solved as a two-by-two Newton iteration in
      ``(T, rhomolar)``.  Every property is an explicit function of
      ``(T, rhomolar)`` for a Helmholtz-energy EOS, so the iteration needs no
      flash at all, only ``DmolarT_INPUTS`` updates on a state whose phase has
      been imposed.  Imposing the phase skips the saturation call that
      dominates a mixture update and, just as importantly, keeps the iterates
      on the single-phase root instead of snapping to a two-phase solution.
    * A two-phase point is bracketed by the two saturation states at the
      constant temperature or pressure of the input pair, and the vapour
      quality between them is found by a bracketed one-dimensional solve.

    The tracer needs temperature or pressure among the two inputs -- see
    :func:`supports` -- because that is what lets it bracket the saturation
    curve.  Every point it cannot bracket, on either side, raises so that the
    caller falls back to a plain :func:`AbstractState.update` call: a
    phase-imposed Newton iteration will happily return the metastable
    extension of the EOS inside the dome, and only the bracket rules that out.
    """

    #: The Newton iteration has settled once the relative step in both T and
    #: rhomolar is this small.  That alone is not convergence -- see
    #: RESIDUAL_RTOL -- but it is what stops the iteration.
    NEWTON_RTOL = 1e-10
    #: Residual accepted at that point, relative to how much the property
    #: varies over the iteration variables.
    RESIDUAL_RTOL = 1e-9
    NEWTON_ITMAX = 30
    #: A single Newton step may not move T or rhomolar by more than this
    #: fraction of their current value.
    NEWTON_TRUST = 0.25
    #: Bracketed solve for the vapour quality of a two-phase point.
    QUALITY_XTOL = 1e-12
    QUALITY_RTOL = 1e-8
    QUALITY_ITMAX = 60
    #: A saturation bracket narrower than this (relative) cannot resolve a
    #: quality: a pure fluid at fixed T brackets p by psat on both sides, and
    #: an azeotropic mixture is barely better.  Such a point is single-phase
    #: as far as the tracer is concerned.
    BRACKET_RTOL = 1e-8
    #: How far outside the saturation densities a converged single-phase root
    #: may sit before it is rejected as the wrong branch.
    BRANCH_RTOL = 1e-6

    #: The reported extent of the two-phase region is a traced curve, so its
    #: end is not exactly the cricondenbar (or cricondentherm).  A point is
    #: only called supercritical this far beyond it.
    DOME_MARGIN = 1.05
    #: A phase envelope that ends with its two densities still this far apart
    #: did not run to the critical point, so where its dome ends is unknown
    #: and has to be found the expensive way.
    ENVELOPE_CLOSED_RTOL = 0.35

    LIQUID = -1
    VAPOUR = 1

    @classmethod
    def supports(cls, index1, index2):
        """Whether an isoline with these two inputs can be traced

        The tracer needs temperature or pressure among the inputs, because
        that is what lets it put a saturation bracket around a point and so
        tell a two-phase state from a single-phase one.  Without it the Newton
        iteration would run on a phase-imposed state with nothing to stop it
        returning the metastable extension of the EOS inside the dome -- a
        plausible-looking answer that is simply wrong.  ``HmassSmass_INPUTS``
        and ``DmassSmass_INPUTS`` are the pairs this rules out.
        """
        return CoolProp.iT in (index1, index2) or CoolProp.iP in (index1, index2)

    def __init__(self, state, index1, index2, iso_index):
        """
        Parameters
        ----------
        state : CoolProp.AbstractState
            The state the isoline belongs to; it is cloned, never modified.
        index1, index2 : int
            Parameter indices of the two inputs, in the order the input pair
            expects them.
        iso_index : int
            Parameter index the isoline holds constant.
        """
        if index1 == index2:
            raise ValueError("The two inputs of an isoline cannot be the same property.")
        if not self.supports(index1, index2):
            raise ValueError("An isoline whose inputs include neither temperature "
                             "nor pressure cannot be traced.")
        self._i1 = index1
        self._i2 = index2
        self._sat = _clone_state(state)
        self._sp = _clone_state(state)
        self._check_clone(state)
        # Not every backend honours an imposed phase; the tracer is still
        # correct without it, just slower and less robust, so fail soft.
        try:
            self._sp.specify_phase(CoolProp.iphase_gas)
        except Exception:
            pass  # backend without an imposed phase; correctness does not depend on it
        # Using the isoline's own constant as the saturation variable makes the
        # bracket constant along the whole line, so it is computed only once.
        if iso_index in (CoolProp.iT, CoolProp.iP) and iso_index in (index1, index2):
            self._sat_index = iso_index
        elif CoolProp.iT in (index1, index2):
            self._sat_index = CoolProp.iT
        else:
            self._sat_index = CoolProp.iP
        self._eos_bounds = self._find_eos_bounds()
        envelope_bounds, self._saturation_bounds = self._find_saturation_bounds()
        self._dome_limit = envelope_bounds.get(self._sat_index)
        self._dome_T_limit, self._dome_p_limit = self._find_dome_limits()
        self._sat_value = None
        self._bracket = None
        self._bracket_failure = None
        self._T = np.nan
        self._rhomolar = np.nan
        self._side = None
        self._out = None

    def _check_clone(self, state):
        """Refuse to trace if the clone is not the caller's fluid model

        The numeric binary interaction parameters are copied across, but a
        departure function installed with ``set_binary_interaction_string``
        cannot be read back from Python at all, so the clone would quietly
        fall back on the default one and the isoline would mix two models --
        traced points from the defaults, flashed points from the caller's.
        Pressure at a given ``(T, rhomolar)`` is a direct function of the
        model, so one comparison settles it.

        The one case this cannot cover: a caller whose state is two-phase or
        has never been updated gives nothing to compare against, and a
        departure function is invisible either way, so a mixture customised
        that way and handed over in that state would be traced with the
        default model while its fallback points were flashed with the
        caller's.  Everything reachable through ``set_binary_interaction_double``
        and ``apply_simple_mixing_rule`` is copied outright and unaffected,
        and ``PropertyPlot`` leaves its state single-phase, so the probe does
        run on the path that matters.
        """
        try:
            if state.Q() >= 0.0:
                return
            T, rhomolar, expected = state.T(), state.rhomolar(), state.p()
        except Exception:
            return  # nothing to compare the clone against; see the docstring
        if not (np.isfinite(T) and np.isfinite(rhomolar) and np.isfinite(expected)):
            return
        self._sp.update(CoolProp.DmolarT_INPUTS, rhomolar, T)
        got = self._sp.p()
        if abs(got - expected) > 1e-9 * max(abs(expected), abs(got)):
            raise ValueError("The traced copy of this state does not reproduce it; "
                             "the fluid model cannot be cloned.")

    @property
    def state(self):
        """The state holding the most recently traced point

        Note that its ``phase()`` reports the phase the tracer imposed on it,
        not one it determined; ``Q()`` is meaningful, ``phase()`` is not.
        """
        if self._out is None:
            raise ValueError("No point has been traced yet.")
        return self._out

    def keyed_output(self, key):
        """Read a property of the most recently traced point"""
        return self.state.keyed_output(key)

    def set_guess(self, state):
        """Adopt a state the caller solved itself as the warm start

        A two-phase state says nothing about the single-phase root the Newton
        iteration follows, so it is ignored.  Both ``phase()`` and ``Q()`` are
        consulted: the caller's state is shared and may have had a phase
        imposed on it, in which case ``phase()`` reports that rather than what
        was found.
        """
        try:
            if state.phase() == CoolProp.iphase_twophase:
                return
            quality = state.Q()
            # phase() reports whatever the caller imposed on the shared state,
            # so do not rely on it alone.
            if np.isfinite(quality) and 0.0 <= quality <= 1.0:
                return
            T, rhomolar = state.T(), state.rhomolar()
        except Exception:
            return  # unusable state; keep whatever warm start we already had
        if np.isfinite(T) and np.isfinite(rhomolar) and T > 0.0 and rhomolar > 0.0:
            self._T = T
            self._rhomolar = rhomolar
            self._side = None

    def update(self, value1, value2):
        """Trace the state with ``index1 == value1`` and ``index2 == value2``

        Raises if the point could not be traced.  The warm start survives,
        because it holds the last point that *did* converge and that is still
        the best guess available for the next one -- dropping it on every
        failure lets a single awkward point strand the whole rest of the line
        in a region where there is no saturation state to restart from.
        """
        try:
            self._out = None
            targets = {self._i1: value1, self._i2: value2}
            bracket = self._get_bracket(targets)
            if bracket is not None and bracket['two_phase']:
                self._solve_quality(bracket)
            else:
                self._solve_single_phase(targets, bracket)
        except Exception:
            self._out = None
            raise

    # -------------------------------------------------------------------
    # Saturation bracket
    # -------------------------------------------------------------------
    def _find_eos_bounds(self):
        """The range the equation of state is actually valid over

        Cheap -- these are trivial parameters -- and used twice: to keep a
        diverged phase envelope from setting a meaningless bound, and to stop
        the Newton iteration reporting a root far outside the fitted range.
        With the phase imposed there is nothing else to stop it: an isochore
        on R504 will happily extrapolate to tens of thousands of kelvin and
        report finite numbers there.
        """
        bounds = {}
        for index, key in ((CoolProp.iT, CoolProp.iT_max), (CoolProp.iP, CoolProp.iP_max)):
            try:
                value = self._sp.trivial_keyed_output(key)
            except Exception:
                continue
            if np.isfinite(value) and value > 0.0:
                bounds[index] = value
        try:
            value = self._sp.trivial_keyed_output(CoolProp.iT_min)
            if np.isfinite(value) and value > 0.0:
                bounds['T_min'] = value
        except Exception:
            pass  # no lower limit published; the runaway guard uses the upper one
        return bounds

    def _find_saturation_bounds(self):
        """Two bounds on how far the two-phase region reaches, in T and in p

        They answer different questions and must not be conflated.

        The *inner* one is the phase envelope's own extent.  It decides what a
        failed saturation call means: beyond it there is probably no dome, so
        the point is traced without a bracket -- and then checked against the
        saturation densities at its own temperature by
        :func:`_check_outside_the_dome`, which is what actually makes that
        safe.  Being merely probable is therefore fine, and keeping it tight
        is what keeps the fallback rate down.

        The *outer* one also takes in the largest critical point among the
        components, and decides whether a saturation state the solver returned
        can be believed.  Beyond the cricondenbar the mixture PQ and QT
        solvers do not fail: they return a state that reports the requested
        pressure and quality exactly and sits hundreds of kelvin from the
        dome.  This one has to be generous, because the envelope trace stops
        short of the real cricondentherm for several predefined blends --
        R504 by about 20 K -- and those genuine dome states must not be thrown
        away as solver noise.

        Building the envelope is also what makes the saturation calls on
        ``self._sat`` converge in the first place -- the mixture PQ and QT
        flashes take their guesses from it -- so this is deliberately done on
        that object and deliberately not cached between tracers.  Caching it
        made the first isoline of a fluid behave differently from the rest.
        """
        envelope = {}
        try:
            self._sat.build_phase_envelope("none")
            data = self._sat.get_phase_envelope_data()
            envelope = {CoolProp.iT: max(data.T), CoolProp.iP: max(data.p)}
        except Exception:
            pass  # no envelope: every point must be bracketed or handed back
        outer = dict(envelope)
        try:
            for index, name in ((CoolProp.iT, 'Tcrit'), (CoolProp.iP, 'pcrit')):
                critical = max(PropsSI(name, fluid) for fluid in self._sat.fluid_names())
                outer[index] = max(outer.get(index, 0.0), critical)
        except Exception:
            pass  # no component critical data; the envelope alone has to do
        def clean(candidate):
            """Drop unusable bounds, clamp to the EOS range, apply the margin"""
            bounds = {}
            for index, value in candidate.items():
                if not (np.isfinite(value) and value > 0.0):
                    continue
                # A phase envelope can diverge rather than fail -- CO2/nitrogen
                # traces one out to 18000 K -- and an unclamped bound from it
                # makes every plausibility test below a tautology.
                limit = self._eos_bounds.get(index)
                bounds[index] = min(value, limit) if limit is not None else value
            return {k: v * self.DOME_MARGIN for k, v in bounds.items()}
        return clean(envelope), clean(outer)

    def _find_dome_limits(self):
        """The temperature, and pressure, above which there is certainly no dome

        This one has to be right rather than generous: it is what decides
        whether a saturation call that could not answer counts as "there is no
        dome here" or as "this point cannot be checked, hand it back".  Too
        high and every supercritical point pays for a flash; too low and roots
        inside the dome are accepted unexamined.

        A phase envelope that ran to the critical point ends with its liquid
        and vapour densities converging, and then its own maximum temperature
        is the answer for free.  One that stopped early -- R504's ends with
        them still a factor of six apart, 37 K below the real critical point
        -- says nothing, and the critical point has to be located directly,
        which is slow enough to be worth avoiding for every other fluid.
        """
        try:
            data = self._sat.get_phase_envelope_data()
            temperatures = np.asarray(data.T, dtype=float)
            hot = int(np.argmax(temperatures))
            rho_liq = float(np.asarray(data.rhomolar_liq, dtype=float)[hot])
            rho_vap = float(np.asarray(data.rhomolar_vap, dtype=float)[hot])
            closed = abs(rho_liq - rho_vap) <= self.ENVELOPE_CLOSED_RTOL * max(rho_liq, rho_vap)
            if closed:
                # The trace covered the whole dome, so its own extremes are the
                # cricondentherm and the cricondenbar.
                return (temperatures[hot] * self.DOME_MARGIN,
                        float(np.max(np.asarray(data.p, dtype=float))) * self.DOME_MARGIN)
        except Exception:
            pass  # envelope unusable; fall through to locating the critical point
        # An envelope that stopped early bounds neither.  The critical point can
        # be located directly, which is slow enough to be worth avoiding when
        # the envelope already answered, and gives the temperature only: the
        # cricondenbar can lie above the critical pressure, so no pressure
        # bound is claimed here and pressure simply stops short-circuiting.
        try:
            stable = [point.T for point in self._sat.all_critical_points() if point.stable]
            if stable:
                return max(stable) * self.DOME_MARGIN, None
        except Exception:
            pass  # critical point not locatable; fall back on the outer bound
        return self._saturation_bounds.get(CoolProp.iT), None

    def _is_plausible_saturation(self, T, p):
        """Whether a state the saturation solver returned can be on the dome"""
        for index, value in ((CoolProp.iT, T), (CoolProp.iP, p)):
            bound = self._saturation_bounds.get(index)
            if bound is not None and not value <= bound:
                return False
        return True

    def _saturation_endpoint(self, sat_value, quality):
        """Move ``self._sat`` onto the saturation curve

        These deliberately do not go through ``update_with_guesses``.  Seeding
        the saturation solver from the previous point is about twice as fast,
        but near the critical point it converges on a state that is not the
        bubble or dew point and reports no error while doing so, and every
        branch decision below is built on this bracket.  The saturation calls
        are a small part of the cost of a traced isoline anyway.
        """
        if self._sat_index == CoolProp.iT:
            self._sat.update(CoolProp.QT_INPUTS, quality, sat_value)
        else:
            self._sat.update(CoolProp.PQ_INPUTS, sat_value, quality)

    def _get_bracket(self, targets):
        """Bracket the target property between its saturated liquid and vapour values

        Returns a dict describing the bracket and whether the point falls
        inside it, or ``None`` if the point is beyond the dome entirely, where
        an unbracketed single-phase root is safe because there is no second
        root to confuse it with.  Raises if the dome should be there and could
        not be found, so that the caller flashes the point instead.
        """
        sat_value = targets[self._sat_index]
        target_index = self._i2 if self._sat_index == self._i1 else self._i1
        target = targets[target_index]

        if self._sat_value is None or sat_value != self._sat_value:
            # Record the outcome and the value it belongs to together.  Storing
            # the value first and letting the build raise past it would leave
            # the previous outcome -- None, meaning "no dome, nothing to
            # check" -- cached against the new value, and since an isoline in T
            # or p has the same saturation value at every point, one failed
            # build would disable the branch check for the whole line.
            bracket, failure = None, None
            try:
                bracket = self._build_bracket(sat_value, target_index)
            except Exception as e:
                failure = str(e)
            self._sat_value = sat_value
            self._bracket = bracket
            self._bracket_failure = failure
        if self._bracket_failure is not None:
            raise ValueError(self._bracket_failure)
        if self._bracket is None:
            # Above the dome; a single-phase root here needs no checking
            # because there is no other root to be confused with.
            return None

        bracket = dict(self._bracket)
        v_liq, v_vap = bracket['liq'][0], bracket['vap'][0]
        bracket['target'] = target
        bracket['two_phase'] = (not bracket['degenerate']
                                and min(v_liq, v_vap) < target < max(v_liq, v_vap))
        # Which single-phase branch the point sits on.  Following the bracket
        # rather than the previous point is what keeps the iteration on the
        # correct root where an isoline crosses the saturation curve.  The
        # midpoint and the slope are used rather than the two bracket values
        # because for an azeotropic mixture the bubble and dew values can
        # agree to within round-off, and then their order is noise: R513A's
        # bubble and dew pressures come within 1.4e-9 relative of each other
        # near 305 K, below BRACKET_RTOL.
        slope = bracket['slope']
        if bracket['two_phase']:
            bracket['side'] = None
        elif slope is not None:
            bracket['side'] = self.LIQUID if (target - 0.5 * (v_liq + v_vap)) * slope > 0 \
                else self.VAPOUR
        elif not bracket['degenerate']:
            # No slope, but a bracket this wide still separates the branches
            bracket['side'] = self.LIQUID if (target - v_liq) * (v_liq - v_vap) > 0 \
                else self.VAPOUR
        else:
            raise ValueError("Cannot tell which side of the saturation curve "
                             "this point is on.")
        return bracket

    def _build_bracket(self, sat_value, target_index):
        """Saturation states at ``sat_value``, or None if there is no dome there

        Raises rather than returning None when the dome should exist but could
        not be found: an unbracketed point is only safe where there is nothing
        to bracket, and everywhere else the caller has to flash it instead.
        """
        # The bracket is attempted even beyond the reported extent of the dome,
        # because that extent is only a lower bound and a saturation state
        # that really does exist is better evidence than the limit.  A call
        # that fails, or that returns a state the envelope says cannot be on
        # the dome, then falls back on the limit to decide whether it came up
        # empty for want of a dome or for want of a solver.
        ends = {}
        try:
            for name, quality in (('liq', 0.0), ('vap', 1.0)):
                self._saturation_endpoint(sat_value, quality)
                ends[name] = (self._sat.keyed_output(target_index),
                              self._sat.T(), self._sat.rhomolar())
                if not self._is_plausible_saturation(self._sat.T(), self._sat.p()):
                    raise ValueError(
                        "The saturation solver returned a state at T={0:g} K, "
                        "p={1:g} Pa, beyond any possible dome.".format(
                            self._sat.T(), self._sat.p()))
        except Exception:
            if self._dome_limit is not None and sat_value > self._dome_limit:
                return None
            raise
        v_liq, v_vap = ends['liq'][0], ends['vap'][0]
        if not (np.isfinite(v_liq) and np.isfinite(v_vap)):
            raise ValueError("The saturation states are not finite.")
        # A saturation solver that converged on the wrong branch, or onto the
        # trivial solution with both phases at one density, would make every
        # branch decision below meaningless.
        if not ends['liq'][2] > ends['vap'][2] * (1.0 + _TRIVIAL_SOLUTION_RTOL):
            raise ValueError("The saturation states are not on their own branches.")
        degenerate = abs(v_vap - v_liq) <= self.BRACKET_RTOL * max(abs(v_liq), abs(v_vap), 1e-30)
        return {'target_index': target_index, 'liq': ends['liq'], 'vap': ends['vap'],
                'degenerate': degenerate,
                'slope': self._branch_slope(target_index, ends['liq'])}

    def _branch_slope(self, target_index, liquid_end):
        """Sign of the target property along density at the saturated liquid

        A point whose target property lies on the same side of the saturation
        curve as this slope points sits on the liquid branch.  Returns None if
        the derivative is unavailable, in which case the branch is left to the
        previous point.
        """
        try:
            self._sp.update(CoolProp.DmolarT_INPUTS, liquid_end[2], liquid_end[1])
            slope = self._sp.first_partial_deriv(target_index, CoolProp.iDmolar, self._sat_index)
        except Exception:
            return None
        if not np.isfinite(slope) or slope == 0.0:
            return None
        return slope

    # -------------------------------------------------------------------
    # Two-phase point: solve for the vapour quality
    # -------------------------------------------------------------------
    def _quality_residual(self, quality, bracket):
        """How far the two-phase state at this quality is from the target"""
        if self._sat_index == CoolProp.iT:
            self._sat.update(CoolProp.QT_INPUTS, quality, self._sat_value)
        else:
            self._sat.update(CoolProp.PQ_INPUTS, self._sat_value, quality)
        return self._sat.keyed_output(bracket['target_index']) - bracket['target']

    def _solve_quality(self, bracket):
        """Regula falsi (Illinois) on the vapour quality between the saturation states"""
        target = bracket['target']
        a, fa = 0.0, bracket['liq'][0] - target
        b, fb = 1.0, bracket['vap'][0] - target
        width = abs(fb - fa)
        if fa * fb > 0.0:
            raise ValueError("The saturation bracket does not contain the target.")
        for _ in range(self.QUALITY_ITMAX):
            quality = b - fb * (b - a) / (fb - fa)
            # The endpoints are already known and a regula falsi can creep onto
            # them, so keep every trial strictly inside the bracket.
            lo, hi = min(a, b), max(a, b)
            if not np.isfinite(quality) or not lo < quality < hi:
                quality = 0.5 * (lo + hi)
            fq = self._quality_residual(quality, bracket)
            if fq * fb < 0.0:
                a, fa = b, fb
            else:
                fa *= 0.5  # Illinois: stops the retained endpoint from stalling
            b, fb = quality, fq
            if fq == 0.0 or abs(b - a) < self.QUALITY_XTOL:
                break
        else:
            raise ValueError("The quality solver did not converge.")
        # A bracket built on a bad saturation state can be wide enough to
        # appear to contain a target the two-phase region never reaches, and
        # then the iteration runs out of interval with the residual intact.
        if abs(fb) > self.QUALITY_RTOL * width:
            raise ValueError("The quality solver ran out of bracket with the "
                             "residual unresolved.")
        self._out = self._sat
        # A two-phase point says nothing about the single-phase root.
        self._T = np.nan
        self._rhomolar = np.nan
        self._side = None

    # -------------------------------------------------------------------
    # Single-phase point: Newton iteration in (T, rhomolar)
    # -------------------------------------------------------------------
    def _seed(self, bracket):
        """Pick the initial (T, rhomolar), preferring the previous point"""
        side = None if bracket is None else bracket['side']
        if np.isfinite(self._T) and np.isfinite(self._rhomolar) \
                and (side is None or self._side is None or side == self._side):
            return self._T, self._rhomolar
        if bracket is not None:
            end = bracket['liq'] if side == self.LIQUID else bracket['vap']
            return end[1], end[2]
        raise ValueError("No initial guess is available for the Newton iteration.")

    def _solve_single_phase(self, targets, bracket):
        """Newton iteration in (T, rhomolar) onto both target properties

        Raises unless it settles on a root that is checked, in range, and on
        the side of the saturation curve the bracket says it should be.
        """
        T, rhomolar = self._seed(bracket)
        i1, i2, state = self._i1, self._i2, self._sp
        v1, v2 = targets[i1], targets[i2]
        for _ in range(self.NEWTON_ITMAX):
            state.update(CoolProp.DmolarT_INPUTS, rhomolar, T)
            f1 = state.keyed_output(i1) - v1
            f2 = state.keyed_output(i2) - v2
            J11 = state.first_partial_deriv(i1, CoolProp.iT, CoolProp.iDmolar)
            J12 = state.first_partial_deriv(i1, CoolProp.iDmolar, CoolProp.iT)
            J21 = state.first_partial_deriv(i2, CoolProp.iT, CoolProp.iDmolar)
            J22 = state.first_partial_deriv(i2, CoolProp.iDmolar, CoolProp.iT)
            det = J11 * J22 - J12 * J21
            if not np.isfinite(det) or det == 0.0:
                raise ValueError("The Newton iteration hit a singular Jacobian.")
            dT = (J12 * f2 - J22 * f1) / det
            dD = (J21 * f1 - J11 * f2) / det
            if not (np.isfinite(dT) and np.isfinite(dD)):
                raise ValueError("The Newton iteration produced a non-finite step.")
            if abs(dT) <= self.NEWTON_RTOL * T and abs(dD) <= self.NEWTON_RTOL * rhomolar:
                # A small step is not by itself a small residual: where the
                # Jacobian is large -- cv near the critical point, say -- the
                # iteration can stall with the residual intact.  Scale each
                # residual by how much its property moves over the iteration
                # variables, which is finite and non-zero even for an enthalpy
                # that happens to sit near its reference value.
                scale1 = abs(J11) * T + abs(J12) * rhomolar
                scale2 = abs(J21) * T + abs(J22) * rhomolar
                if not (np.isfinite(scale1) and np.isfinite(scale2)):
                    raise ValueError("The residual cannot be scaled at this point.")
                if abs(f1) > self.RESIDUAL_RTOL * scale1 \
                        or abs(f2) > self.RESIDUAL_RTOL * scale2:
                    raise ValueError("The Newton iteration stalled with the "
                                     "residual unresolved.")
                self._check_in_range(T)
                side = None if bracket is None else bracket['side']
                self._check_branch(bracket, side, T, rhomolar)
                self._T, self._rhomolar = T, rhomolar
                self._side = side
                self._out = state
                return
            damp = 1.0
            if abs(dT) > self.NEWTON_TRUST * T:
                damp = min(damp, self.NEWTON_TRUST * T / abs(dT))
            if abs(dD) > self.NEWTON_TRUST * rhomolar:
                damp = min(damp, self.NEWTON_TRUST * rhomolar / abs(dD))
            T += damp * dT
            rhomolar += damp * dD
            if not (T > 0.0 and rhomolar > 0.0):
                raise ValueError("The Newton iteration left the physical domain.")
        raise ValueError("The Newton iteration did not converge.")

    #: How far past the equation of state's stated temperature range a root
    #: may sit.  Deliberately loose: this is not a validity check.  CoolProp
    #: answers a PT flash well outside the fitted range and so would the
    #: fallback this hands back to, so a tighter bound would only make traced
    #: and flashed points on the same isoline obey different rules.  What it
    #: catches is an iteration that has run away entirely -- an R504 isochore
    #: reaching four hundred times the maximum temperature.
    RANGE_FACTOR = 2.0

    def _check_in_range(self, T):
        """Reject a root the Newton iteration has run away with"""
        upper = self._eos_bounds.get(CoolProp.iT)
        lower = self._eos_bounds.get('T_min')
        if (upper is not None and T > upper * self.RANGE_FACTOR) \
                or (lower is not None and T < lower / self.RANGE_FACTOR):
            raise ValueError("The Newton iteration converged at T={0:g} K, far "
                             "outside the range of the equation of state.".format(T))

    def _check_branch(self, bracket, side, T, rhomolar):
        """Reject a root that ended up inside the two-phase dome

        With a bracket this is a comparison against the two saturation
        densities at the constant temperature (or pressure) it was built at: a
        liquid root is denser than the saturated liquid, a vapour root thinner
        than the saturated vapour, and a root that is neither is metastable or
        unstable and goes back to the caller to flash.

        Without one the point was taken to be beyond the dome, and that is
        checked here rather than assumed.  The reported extent of the phase
        envelope is not reliable enough to assume it: for several predefined
        blends the envelope trace stops well short of the real cricondentherm
        -- R504 by about 20 K -- and every point in between would otherwise be
        accepted unexamined, which is exactly where the dome is.
        """
        if bracket is None:
            self._check_outside_the_dome(T, rhomolar)
            return
        if side is None:
            return
        if side == self.LIQUID:
            limit = bracket['liq'][2]
            wrong = rhomolar < limit * (1.0 - self.BRANCH_RTOL)
        else:
            limit = bracket['vap'][2]
            wrong = rhomolar > limit * (1.0 + self.BRANCH_RTOL)
        if wrong:
            raise ValueError("The Newton iteration converged on the wrong side "
                             "of the saturation curve.")

    def _check_outside_the_dome(self, T, rhomolar):
        """Confirm there is no two-phase region around this root

        A saturation call that succeeds here gives the two densities to
        compare against directly.  One that fails is not the reassurance it
        looks like: when the isoline is indexed by temperature this is the
        same call, at the same temperature, that :func:`_build_bracket`
        already gave up on, so it can add nothing.  Treating that as "no dome"
        left a band of about ten kelvin on R504 -- between where the
        saturation solver stops converging and the real critical point at
        335.5 K, which its phase envelope stops 37 K short of -- in which
        roots inside the dome were accepted unexamined.  So a failure only
        counts as an answer above the temperature where a dome could exist at
        all; below it the point goes back to the caller to flash.
        """
        if self._dome_T_limit is not None and T > self._dome_T_limit:
            return          # above the critical point; there is no dome to be in
        if self._sat_index == CoolProp.iP and self._dome_p_limit is not None \
                and self._sat_value is not None and self._sat_value > self._dome_p_limit:
            return          # above the cricondenbar; this isobar has no dome
        undecided = None
        try:
            self._sat.update(CoolProp.QT_INPUTS, 0.0, T)
            rho_liq, p_liq = self._sat.rhomolar(), self._sat.p()
            self._sat.update(CoolProp.QT_INPUTS, 1.0, T)
            rho_vap, p_vap = self._sat.rhomolar(), self._sat.p()
        except Exception:
            undecided = "the saturation solver did not converge"
        if undecided is None and not self._is_plausible_saturation(T, max(p_liq, p_vap)):
            undecided = "the saturation solver left the dome behind"
        if undecided is None and not rho_liq > rho_vap:
            undecided = "the saturation states are not on their own branches"
        if undecided is not None:
            bound = self._dome_T_limit
            if bound is None or T <= bound:
                raise ValueError("This root cannot be shown to be outside the "
                                 "two-phase dome: {0}.".format(undecided))
            return          # above any possible dome, so there is none to be in
        if rho_vap < rhomolar < rho_liq:
            raise ValueError("The Newton iteration converged inside the "
                             "two-phase dome.")


class IsoLine(Base2DObject):
    """An object that holds the functions to calculate a line of
    a constant property in the dimensions of a property plot. This
    class only uses SI units."""

    # Normally we calculate a sweep in x-dimensions, but
    # sometimes a sweep in y-dimensions is better.
    XY_SWITCH = {
      CoolProp.iDmass: {Base2DObject.TS: True, Base2DObject.PH: True, Base2DObject.HS: False, Base2DObject.PS: True, Base2DObject.PD: None, Base2DObject.TD: None, Base2DObject.PT: False},
      CoolProp.iHmass: {Base2DObject.TS: False, Base2DObject.PH: None, Base2DObject.HS: None, Base2DObject.PS: True, Base2DObject.PD: True, Base2DObject.TD: False, Base2DObject.PT: False},
      CoolProp.iP: {Base2DObject.TS: False, Base2DObject.PH: None, Base2DObject.HS: False, Base2DObject.PS: None, Base2DObject.PD: None, Base2DObject.TD: False, Base2DObject.PT: None},
      CoolProp.iSmass: {Base2DObject.TS: None, Base2DObject.PH: True, Base2DObject.HS: None, Base2DObject.PS: None, Base2DObject.PD: True, Base2DObject.TD: False, Base2DObject.PT: True},
      CoolProp.iT: {Base2DObject.TS: None, Base2DObject.PH: True, Base2DObject.HS: False, Base2DObject.PS: False, Base2DObject.PD: False, Base2DObject.TD: None, Base2DObject.PT: None},
      CoolProp.iQ: {Base2DObject.TS: True, Base2DObject.PH: True, Base2DObject.HS: True, Base2DObject.PS: True, Base2DObject.PD: True, Base2DObject.TD: True, Base2DObject.PT: False}
    }

    # Abort interpolation if there are not enough
    # valid entries.
    VALID_REQ = 5.0 / 100.0

    def __init__(self, i_index, x_index, y_index, value=0.0, state=None, tracing=True):
        super().__init__(x_index, y_index, state)
        self._i_index = _get_index(i_index)
        if value is not None: self.value = value
        else: self._value = None
        self._x = None
        self._y = None
        self._tracing = bool(tracing)

    @property
    def i_index(self): return self._i_index

    @property
    def tracing(self):
        """Whether :class:`IsoLineTracer` is used to walk along the isoline

        Set this to False to calculate every point with an independent flash,
        which is much slower but does not depend on the neighbouring points.
        """
        return self._tracing

    @tracing.setter
    def tracing(self, value):
        """Turn isoline tracing on or off for this line"""
        self._tracing = bool(value)

    @property
    def value(self): return self._value

    @value.setter
    def value(self, value): self._value = float(value)

    @property
    def x(self): return self._x

    @x.setter
    def x(self, value): self._x = np.array(value)

    @property
    def y(self): return self._y

    @y.setter
    def y(self, value): self._y = np.array(value)

    def get_update_pair(self):
        """Processes the values for the isoproperty and the graph dimensions
        to figure which should be used as inputs to the state update. Returns
        a tuple with the indices for the update call and the property constant.
        For an isobar in a Ts-diagram it returns the default order and the
        correct constant for the update pair:
        get_update_pair(CoolProp.iP,CoolProp.iSmass,CoolProp.iT) -> (0,1,2,CoolProp.PSmass_INPUTS)
        other values require switching and swapping.
        """
        # Figure out if x or y-dimension should be used
        switch = self.XY_SWITCH[self.i_index][self.y_index * 10 + self.x_index]

        if switch is None:
            raise ValueError("This isoline cannot be calculated!")
        elif switch is False:
            pair, out1, _ = CP.generate_update_pair(self.i_index, 0.0, self.x_index, 1.0)
        elif switch is True:
            pair, out1, _ = CP.generate_update_pair(self.i_index, 0.0, self.y_index, 1.0)
        else:
            raise ValueError("Unknown error!")

        if out1 == 0.0:  # Correct order
            swap = False
        else:  # Wrong order
            swap = True

        if not switch and not swap:
            return 0, 1, 2, pair
        elif switch and not swap:
            return 0, 2, 1, pair
        elif not switch and swap:
            return 1, 0, 2, pair
        elif switch and swap:
            return 1, 2, 0, pair
        else:
            raise ValueError("Check the code, this should not happen!")

    def calc_sat_range(self, Trange=None, Prange=None, num=200):
        if Trange is not None:
            two = np.array(Trange)
            one = np.resize(np.array(self.value), two.shape)
            pair = CoolProp.QT_INPUTS
        elif Prange is not None:
            one = np.array(Prange)
            two = np.resize(np.array(self.value), one.shape)
            pair = CoolProp.PQ_INPUTS
        else:
            T_lo, T_hi = self._get_sat_bounds(CoolProp.iT)
            two = np.linspace(T_lo, T_hi, num)
            one = np.resize(np.array(self.value), two.shape)
            pair = CoolProp.QT_INPUTS

        Tcrit = self.critical_state.keyed_output(CoolProp.iT)
        Pcrit = self.critical_state.keyed_output(CoolProp.iP)
        Dcrit = self.critical_state.keyed_output(CoolProp.iDmass)
        try:
            #self.state.update(CoolProp.DmassT_INPUTS, Dcrit, Tcrit)
            #xcrit = self.state.keyed_output(self._x_index)
            #ycrit = self.state.keyed_output(self._y_index)
            xcrit = self.critical_state.keyed_output(self._x_index)
            ycrit = self.critical_state.keyed_output(self._y_index)
        except:
            warnings.warn(
              "An error occurred for the critical inputs, skipping it.",
              UserWarning)
            xcrit = np.nan
            ycrit = np.nan

        X = np.empty_like(one)
        Y = np.empty_like(one)

        err = False
        for index, _ in np.ndenumerate(one):
            try:
                self.state.update(pair, one[index], two[index])
                if not _phases_are_distinct(self.state):
                    raise ValueError(
                      "The saturation solver returned the trivial solution, with both "
                      "phases at the same density.")
                X[index] = self.state.keyed_output(self._x_index)
                Y[index] = self.state.keyed_output(self._y_index)
            except Exception as e:
                if (pair == CoolProp.QT_INPUTS and abs(two[index] - Tcrit) < 1e0) or \
                   (pair == CoolProp.PQ_INPUTS and abs(one[index] - Pcrit) < 1e2):
                    X[index] = xcrit
                    Y[index] = ycrit
                    warnings.warn(
                  "An error occurred for near critical inputs {0:f}, {1:f} with index {2:s}: {3:s}".format(one[index], two[index], str(index), str(e)),
                  UserWarning)
                    pass

                warnings.warn(
                  "An error occurred for inputs {0:f}, {1:f} with index {2:s}: {3:s}".format(one[index], two[index], str(index), str(e)),
                  UserWarning)
                X[index] = np.nan
                Y[index] = np.nan
                err = True
        self.x = X; self.y = Y
        return

    def calc_range(self, xvals=None, yvals=None):

        if self.i_index == CoolProp.iQ:
            warnings.warn(
                "Please use \"calc_sat_range\" to calculate saturation and isoquality lines. Input ranges are discarded.",
                UserWarning)
            if xvals is not None: self.calc_sat_range(num=xvals.size)
            elif yvals is not None: self.calc_sat_range(num=yvals.size)
            else: self.calc_sat_range()
            return

        ipos, xpos, ypos, pair = self.get_update_pair()

        order = [ipos, xpos, ypos]
        idxs = [v for (_, v) in sorted(zip(order, [self.i_index, self.x_index, self.y_index]))]
        vals = [v for (_, v) in sorted(zip(order, [np.array(self.value), xvals, yvals]))]
        if vals[0] is None or vals[1] is None:
            raise ValueError("One required input is missing, make sure to supply the correct xvals ({0:s}) or yvals ({1:s}).".format(str(xvals), str(yvals)))

        if vals[0].size > vals[1].size:
            vals[1] = np.resize(vals[1], vals[0].shape)
        elif vals[0].size < vals[1].size:
            vals[0] = np.resize(vals[0], vals[1].shape)

        vals[2] = np.empty_like(vals[0])
        err = False
        # Consecutive points of an isoline are close together, so each one can
        # be warm-started from the last instead of being flashed cold.  For a
        # mixture that is orders of magnitude faster and, because the iteration
        # never has to find the region on its own, far less likely to leave a
        # gap in the line.  Anything the tracer will not take is flashed.
        tracer = None
        if self.tracing and IsoLineTracer.supports(idxs[0], idxs[1]):
            try:
                tracer = IsoLineTracer(self.state, idxs[0], idxs[1], self.i_index)
            except Exception as e:
                warnings.warn(
                  "Could not set up the isoline tracer, falling back to flash calculations: {0:s}".format(str(e)),
                  UserWarning)
        for index, _ in np.ndenumerate(vals[0]):
            traced = False
            if tracer is not None:
                try:
                    tracer.update(vals[0][index], vals[1][index])
                    vals[2][index] = tracer.keyed_output(idxs[2])
                    traced = np.isfinite(vals[2][index])
                except Exception:
                    traced = False
            if not traced:
                try:
                    self.state.update(pair, vals[0][index], vals[1][index])
                    vals[2][index] = self.state.keyed_output(idxs[2])
                    if tracer is not None:
                        tracer.set_guess(self.state)
                except Exception as e:
                    warnings.warn(
                      "An error occurred for inputs {0:f}, {1:f} with index {2:s}: {3:s}".format(vals[0][index], vals[1][index], str(index), str(e)),
                      UserWarning)
                    vals[2][index] = np.nan
                    err = True

        for i, v in enumerate(idxs):
            if v == self.x_index: self.x = vals[i]
            if v == self.y_index: self.y = vals[i]

    def sanitize_data(self):
        """Fill the series via interpolation"""
        validx = None; validy = None
        countx = None; county = None
        if self.x is not None:
            validx = np.isfinite(self.x)
            countx = float(self.x.size)
        else:
            raise ValueError("The x-axis is not populated, calculate values before you interpolate.")
        if self.y is not None:
            validy = np.isfinite(self.y)
            county = float(self.y.size)
        else:
            raise ValueError("The y-axis is not populated, calculate values before you interpolate.")

        if min([np.sum(validx) / countx, np.sum(validy) / county]) < self.VALID_REQ:
            warnings.warn(
              "Poor data quality, there are not enough valid entries for x ({0:f}/{1:f}) or y ({2:f}/{3:f}).".format(np.sum(validx), countx, np.sum(validy), county),
              UserWarning)
        # TODO: use filter and cubic splines!
        #filter = np.logical_and(np.isfinite(self.x),np.isfinite(self.y))
        if np.sum(validy) > np.sum(validx):
            self.x = interpolate_values_1d(self.y, self.x, x_points=self.y[validy])
            self.y = self.y[validy]
        else:
            self.y = interpolate_values_1d(self.x, self.y, x_points=self.x[validx])
            self.x = self.x[validx]


class BasePlot(Base2DObject):
    """The base class for all plots. It can be instantiated itself, but provides many
    general facilities to be used in the different plots. """

    # Define the iteration keys
    PROPERTIES = {
      CoolProp.iDmass: 'density',
      CoolProp.iHmass: 'specific enthalpy',
      CoolProp.iP: 'pressure',
      CoolProp.iSmass: 'specific entropy',
      CoolProp.iT: 'temperature',
      CoolProp.iUmass: 'specific internal energy'
    }

    # Define the unit systems
    UNIT_SYSTEMS = {
      'SI': SIunits(),
      'KSI': KSIunits(),
      'EUR': EURunits()
    }

    LINE_PROPS = {
      CoolProp.iT: dict(color='Darkred', lw=0.25),
      CoolProp.iP: dict(color='DarkCyan', lw=0.25),
      CoolProp.iHmass: dict(color='DarkGreen', lw=0.25),
      CoolProp.iDmass: dict(color='DarkBlue', lw=0.25),
      CoolProp.iSmass: dict(color='DarkOrange', lw=0.25),
      CoolProp.iQ: dict(color='black', lw=0.25)
    }

    ID_FACTOR = 10.0  # Values below this number are interpreted as factors
    HI_FACTOR = 2.25  # Upper default limits: HI_FACTOR*T_crit and HI_FACTOR*p_crit
    LO_FACTOR = 1.01  # Lower default limits: LO_FACTOR*T_triple and LO_FACTOR*p_triple

    TP_LIMITS = {
      'NONE': [None, None, None, None],
      'DEF': [LO_FACTOR, HI_FACTOR, LO_FACTOR, HI_FACTOR],
      'ACHP': [173.15, 493.15, 0.25e5, HI_FACTOR],
      'ORC': [273.15, 673.15, 0.25e5, HI_FACTOR]
    }

    def __init__(self, fluid_ref, graph_type, unit_system='KSI', tp_limits='DEF', **kwargs):

        # Process the graph_type and set self._x_type and self._y_type
        graph_type = graph_type.upper()
        graph_type = graph_type.replace(r'RHO', r'D')
        if graph_type not in Base2DObject.PLOTS:
            raise ValueError("Invalid graph_type input, expected a string from {0:s}".format(str(self.PLOTS)))

        # Process the unit_system and set self._system
        self.system = unit_system
        # Process the plotting range based on T and p
        self.limits = tp_limits
        # Other properties
        self.figure = kwargs.pop('figure', plt.figure(tight_layout=True))
        self.axis = kwargs.pop('axis', self.figure.add_subplot(111))
        self.props = kwargs.pop('props', None)

        # call the base class
        state = process_fluid_state(fluid_ref)
        Base2DObject.__init__(self, graph_type[1], graph_type[0], state, **kwargs)

    @property
    def system(self): return self._system

    @system.setter
    def system(self, value):
        value = value.upper()
        if value in self.UNIT_SYSTEMS: self._system = self.UNIT_SYSTEMS[value]
        else: raise ValueError("Invalid input, expected a string from {0:s}".format(str(self.UNIT_SYSTEMS.keys())))

    @property
    def limits(self):
        """Returns [Tmin,Tmax,pmin,pmax] as value or factors"""
        return self._limits

    @limits.setter
    def limits(self, value):
        if is_string(value):
            value = value.upper()
        if value in self.TP_LIMITS:
            self._limits = self.TP_LIMITS[value]
        elif len(value) == 4:
            self._limits = value
        else:
            raise ValueError("Invalid input, expected a list with 4 items or a string from {0:s}".format(str(self.TP_LIMITS.keys())))

    @property
    def figure(self): return self._figure

    @figure.setter
    def figure(self, value): self._figure = value

    @property
    def axis(self): return self._axis

    @axis.setter
    def axis(self, value): self._axis = value

    @property
    def props(self): return self._props

    @props.setter
    def props(self, value):
        self._props = self.LINE_PROPS.copy()
        if value is not None:
            self._props.update(value)

    def __sat_bounds(self, kind, smin=None, smax=None):
        warnings.warn(
          "You called the deprecated function \"__sat_bounds\", \
consider replacing it with \"_get_sat_bounds\".",
          DeprecationWarning)
        return self._get_sat_bounds(kind, smin, smax)

    def _get_iso_label(self, isoline, unit=True):
        if self._system is not None:
            dim = self._system[isoline.i_index]
            return str(r"$" + dim.symbol + "=" + str(dim.from_SI(isoline.value)) + "$ " + dim.unit if unit else "$").strip()
        return str(isoline.value).strip()

    # def _get_phase_envelope(self):
    #
    #HEOS = CoolProp.AbstractState("HEOS", fluid)
    # HEOS.build_phase_envelope("")
    #PED = HEOS.get_phase_envelope_data()
    #plt.plot(PED.T, np.log(PED.p))
    # plt.show()

    def _plot_default_annotations(self):
#         def filter_fluid_ref(fluid_ref):
#             fluid_ref_string = fluid_ref
#             if fluid_ref.startswith('REFPROP-MIX'):
#                 end = 0
#                 fluid_ref_string = ''
#                 while fluid_ref.find('[', end + 1) != -1:
#                     start = fluid_ref.find('&', end + 1)
#                     if end == 0:
#                         start = fluid_ref.find(':', end + 1)
#                     end = fluid_ref.find('[', end + 1)
#                     fluid_ref_string = ' '.join([fluid_ref_string,
#                                                 fluid_ref[start+1:end], '+'])
#                 fluid_ref_string = fluid_ref_string[0:len(fluid_ref_string)-2]
#             return fluid_ref_string
#
#         if len(self.graph_type) == 2:
#             y_axis_id = self.graph_type[0]
#             x_axis_id = self.graph_type[1]
#         else:
#             y_axis_id = self.graph_type[0]
#             x_axis_id = self.graph_type[1:len(self.graph_type)]
#
#         tl_str = "%s - %s Graph for %s"
#         if not self.axis.get_title():
#             self.axis.set_title(tl_str % (self.AXIS_LABELS[self.unit_system][y_axis_id][0],
#                                           self.AXIS_LABELS[self.unit_system][x_axis_id][0],
#                                           filter_fluid_ref(self.fluid_ref)))
        if self._x_index in [CoolProp.iDmass, CoolProp.iP]:
            self.axis.set_xscale('log')
        if self._y_index in [CoolProp.iDmass, CoolProp.iP]:
            self.axis.set_yscale('log')

        if not self.axis.get_xlabel():
            dim = self._system[self._x_index]
            self.xlabel((dim.label + u" $" + dim.symbol + u"$ / " + dim.unit).strip())
        if not self.axis.get_ylabel():
            dim = self._system[self._y_index]
            self.ylabel((dim.label + u" $" + dim.symbol + u"$ / " + dim.unit).strip())

    def title(self, title):
        self.axis.set_title(title)

    def xlabel(self, xlabel):
        self.axis.set_xlabel(xlabel)

    def ylabel(self, ylabel):
        self.axis.set_ylabel(ylabel)

    def grid(self, b=None, **kwargs):
        g_map = {'on': True, 'off': False}
        if b is not None:
            b = g_map[b.lower()]
        if not kwargs:  # len=0
            self.axis.grid(b)
        else:
            self.axis.grid(kwargs)

    def set_Tp_limits(self, limits):
        """Set the limits for the graphs in temperature and pressure, based on
        the active units: [Tmin, Tmax, pmin, pmax]"""
        dim = self._system[CoolProp.iT]
        limits[0] = dim.to_SI(limits[0])
        limits[1] = dim.to_SI(limits[1])
        dim = self._system[CoolProp.iP]
        limits[2] = dim.to_SI(limits[2])
        limits[3] = dim.to_SI(limits[3])
        self.limits = tuple(limits)

    def get_Tp_limits(self):
        """Get the limits for the graphs in temperature and pressure, based on
        the active units: [Tmin, Tmax, pmin, pmax]"""
        limits = self._get_Tp_limits()
        dim = self._system[CoolProp.iT]
        limits[0] = dim.from_SI(limits[0])
        limits[1] = dim.from_SI(limits[1])
        dim = self._system[CoolProp.iP]
        limits[2] = dim.from_SI(limits[2])
        limits[3] = dim.from_SI(limits[3])
        return limits

    def _get_Tp_limits(self):
        """Get the limits for the graphs in temperature and pressure, based on
        SI units: [Tmin, Tmax, pmin, pmax]"""
        T_lo, T_hi, P_lo, P_hi = self.limits
        Ts_lo, Ts_hi = self._get_sat_bounds(CoolProp.iT)
        Ps_lo, Ps_hi = self._get_sat_bounds(CoolProp.iP)

        if T_lo is None: T_lo = 0.0
        elif T_lo < self.ID_FACTOR: T_lo *= Ts_lo
        if T_hi is None: T_hi = 1e6
        elif T_hi < self.ID_FACTOR: T_hi *= Ts_hi
        if P_lo is None: P_lo = 0.0
        elif P_lo < self.ID_FACTOR: P_lo *= Ps_lo
        if P_hi is None: P_hi = 1e10
        elif P_hi < self.ID_FACTOR: P_hi *= Ps_hi

        try: T_lo = np.nanmax([T_lo, self._state.trivial_keyed_output(CoolProp.iT_min)])
        except: pass
        try: T_hi = np.nanmin([T_hi, self._state.trivial_keyed_output(CoolProp.iT_max)])
        except: pass
        try: P_lo = np.nanmax([P_lo, self._state.trivial_keyed_output(CoolProp.iP_min)])
        except: pass
        try: P_hi = np.nanmin([P_hi, self._state.trivial_keyed_output(CoolProp.iP_max)])
        except: pass

        return [T_lo, T_hi, P_lo, P_hi]

    def set_axis_limits(self, limits):
        """Set the limits of the internal axis object based on the active units,
        takes [xmin, xmax, ymin, ymax]"""
        self.axis.set_xlim([limits[0], limits[1]])
        self.axis.set_ylim([limits[2], limits[3]])

    def _set_axis_limits(self, limits):
        """Set the limits of the internal axis object based on SI units,
        takes [xmin, xmax, ymin, ymax]"""
        dim = self._system[self._x_index]
        self.axis.set_xlim([dim.from_SI(limits[0]), dim.from_SI(limits[1])])
        dim = self._system[self._y_index]
        self.axis.set_ylim([dim.from_SI(limits[2]), dim.from_SI(limits[3])])

    def get_axis_limits(self, x_index=None, y_index=None):
        """Returns the previously set limits or generates them and
        converts the default values to the selected unit system.
        Returns a list containing [xmin, xmax, ymin, ymax]"""
        if x_index is None: x_index = self._x_index
        if y_index is None: y_index = self._y_index

        if x_index != self.x_index or y_index != self.y_index or \
          self.axis.get_autoscalex_on() or self.axis.get_autoscaley_on():
            # One of them is not set or we work on a different set of axes
            T_lo, T_hi, P_lo, P_hi = self._get_Tp_limits()

            X = [0.0] * 4; Y = [0.0] * 4
            i = -1
            for T in [T_lo, T_hi]:
                for P in [P_lo, P_hi]:
                    i += 1
                    try:
                        self._state.update(CoolProp.PT_INPUTS, P, T)
                        # TODO: include a check for P and T?
                        X[i] = self._state.keyed_output(x_index)
                        Y[i] = self._state.keyed_output(y_index)
                    except:
                        X[i] = np.nan; Y[i] = np.nan
            # Figure out what to update
            dim = self._system[x_index]
            x_lim = [dim.from_SI(np.nanmin(X)), dim.from_SI(np.nanmax(X))]
            dim = self._system[y_index]
            y_lim = [dim.from_SI(np.nanmin(Y)), dim.from_SI(np.nanmax(Y))]
            # Either update the axes limits or get them
            if x_index == self._x_index:
                if self.axis.get_autoscalex_on():
                    self.axis.set_xlim(x_lim)
                else:
                    x_lim = self.axis.get_xlim()
            if y_index == self._y_index:
                if self.axis.get_autoscaley_on():
                    self.axis.set_ylim(y_lim)
                else:
                    y_lim = self.axis.get_ylim()
        else:  # We only asked for the real axes limits and they are set already
            x_lim = self.axis.get_xlim()
            y_lim = self.axis.get_ylim()

        return [x_lim[0], x_lim[1], y_lim[0], y_lim[1]]

    def _get_axis_limits(self, x_index=None, y_index=None):
        """Get the limits of the internal axis object in SI units
        Returns a list containing [xmin, xmax, ymin, ymax]"""
        if x_index is None: x_index = self._x_index
        if y_index is None: y_index = self._y_index
        limits = self.get_axis_limits(x_index, y_index)
        dim = self._system[x_index]
        limits[0] = dim.to_SI(limits[0])
        limits[1] = dim.to_SI(limits[1])
        dim = self._system[y_index]
        limits[2] = dim.to_SI(limits[2])
        limits[3] = dim.to_SI(limits[3])
        return limits

    @staticmethod
    def generate_ranges(itype, imin, imax, num):
        """Generate a range for a certain property"""
        if itype in [CoolProp.iP, CoolProp.iDmass]:
            return np.logspace(np.log2(imin), np.log2(imax), num=num, base=2.)
        return np.linspace(imin, imax, num=num)

    def _get_conversion_data(self):
        [Axmin, Axmax, Aymin, Aymax] = self._get_axis_limits()
        DELTAX_axis = Axmax - Axmin
        DELTAY_axis = Aymax - Aymin
        width = self.figure.get_figwidth()
        height = self.figure.get_figheight()
        pos = self.axis.get_position().get_points()
        [[Fxmin, Fymin], [Fxmax, Fymax]] = pos
        DELTAX_fig = width * (Fxmax - Fxmin)
        DELTAY_fig = height * (Fymax - Fymin)
        return [[Axmin, Axmax, Aymin, Aymax, Fxmin, Fxmax, Fymin, Fymax], [DELTAX_axis, DELTAY_axis, DELTAX_fig, DELTAY_fig]]

    def _to_pixel_coords(self, xv, yv):
        [[Axmin, Axmax, Aymin, Aymax, Fxmin, Fxmax, Fymin, Fymax], [DELTAX_axis, DELTAY_axis, DELTAX_fig, DELTAY_fig]] = self._get_conversion_data()
        # Convert coords to pixels
        x = (xv - Axmin) / DELTAX_axis * DELTAX_fig + Fxmin
        y = (yv - Aymin) / DELTAY_axis * DELTAY_fig + Fymin
        return x, y

    def _to_data_coords(self, xv, yv):
        [[Axmin, Axmax, Aymin, Aymax, Fxmin, Fxmax, Fymin, Fymax], [DELTAX_axis, DELTAY_axis, DELTAX_fig, DELTAY_fig]] = self._get_conversion_data()
        # Convert back to measurements
        x = (xv - Fxmin) / DELTAX_fig * DELTAX_axis + Axmin
        y = (yv - Fymin) / DELTAY_fig * DELTAY_axis + Aymin
        return x, y

    @staticmethod
    def get_x_y_dydx(xv, yv, x):
        """Get x and y coordinates and the linear interpolation derivative"""
        # Old implementation:
        # Get the rotation angle
        #f = interp1d(xv, yv)
        #y = f(x)
        #h = 0.00001*x
        #dy_dx = (f(x+h)-f(x-h))/(2*h)
        # return x,y,dy_dx
        if len(xv) == len(yv) and len(yv) > 1:  # assure same length
            if len(xv) == len(yv) and len(yv) == 2:  # only two points
                if np.min(xv) < x < np.max(xv):
                    dx = xv[1] - xv[0]
                    dy = yv[1] - yv[0]
                    dydx = dy / dx
                    y = yv[0] + dydx * (x - xv[0])
                    return x, y, dydx
                else:
                    raise ValueError("Your coordinate has to be between the input values.")
            else:
                limit = 1e-10                    # avoid hitting a point directly
                diff = np.array(xv) - x        # get differences
                index = np.argmin(diff * diff)  # nearest neighbour
                if (xv[index] < x < xv[index + 1]      # nearest below, positive inclination
                  or xv[index] > x > xv[index + 1]):   # nearest above, negative inclination
                    if diff[index] < limit:
                        index = [index - 1, index + 1]
                    else:
                        index = [index, index + 1]
                elif (xv[index - 1] < x < xv[index]    # nearest above, positive inclination
                  or xv[index - 1] > x > xv[index]):   # nearest below, negative inclination
                    if diff[index] < limit:
                        index = [index - 1, index + 1]
                    else:
                        index = [index - 1, index]
                xvnew = xv[index]
                yvnew = yv[index]
                return BasePlot.get_x_y_dydx(xvnew, yvnew, x)  # Allow for a single recursion
        else:
            raise ValueError("You have to provide the same amount of x- and y-pairs with at least two entries each.")

    def _inline_label(self, xv, yv, x=None, y=None):
        """
        This will give the coordinates and rotation required to align a label with
        a line on a plot in SI units.

        Exactly one of ``x`` or ``y`` must be supplied; the other is found by
        interpolation along ``(xv, yv)``.
        """
        if (x is None) == (y is None):
            raise ValueError("Provide exactly one of `x` or `y`; the other is interpolated.")

        if y is None:
            trash = 0
            (xv, yv) = self._to_pixel_coords(xv, yv)
            # x is provided but y isn't
            (x, trash) = self._to_pixel_coords(x, trash)

            # Get the rotation angle and y-value
            x, y, dy_dx = BasePlot.get_x_y_dydx(xv, yv, x)
            rot = np.arctan(dy_dx) / np.pi * 180.

        else:
            # y is provided, but x isn't
            _xv = xv[::-1]
            _yv = yv[::-1]
            # Find x by interpolation
            x = interpolate_values_1d(yv, xv, x_points=y)
            trash = 0
            (xv, yv) = self._to_pixel_coords(xv, yv)
            (x, trash) = self._to_pixel_coords(x, trash)

            # Get the rotation angle and y-value
            x, y, dy_dx = BasePlot.get_x_y_dydx(xv, yv, x)
            rot = np.arctan(dy_dx) / np.pi * 180.
        (x, y) = self._to_data_coords(x, y)
        return (x, y, rot)

    def inline_label(self, xv, yv, x=None, y=None):
        """
        This will give the coordinates and rotation required to align a label with
        a line on a plot in axis units.
        """
        dimx = self._system[self._x_index]
        xv = dimx.to_SI(xv)
        if x is not None: x = dimx.to_SI(x)
        dimy = self._system[self._y_index]
        yv = dimy.to_SI(yv)
        if y is not None: y = dimy.to_SI(y)
        (x, y, rot) = self._inline_label(xv, yv, x, y)
        x = dimx.from_SI(x)
        y = dimy.from_SI(y)
        return (x, y, rot)

    def show(self):
        plt.show()

    def savefig(self, *args, **kwargs):
        self.figure.savefig(*args, **kwargs)


if __name__ == "__main__":
    for sys in [SIunits(), KSIunits(), EURunits()]:
        print(sys.H.label)
        print(sys.H.to_SI(20))
        print(sys.P.label)
        print(sys.P.to_SI(20))

    # i_index, x_index, y_index, value=None, state=None)
    iso = IsoLine('T', 'H', 'P')
    print(iso.get_update_pair())

    state = AbstractState("HEOS", "water")
    iso = IsoLine('T', 'H', 'P', 300.0, state)
    hr = PropsSI("H", "T", [290, 310], "P", [1e5, 1e5], "water")
    pr = np.linspace(0.9e5, 1.1e5, 3)
    iso.calc_range(hr, pr)
    print(iso.x, iso.y)

    iso = IsoLine('Q', 'H', 'P', 0.0, state)
    iso.calc_range(hr, pr); print(iso.x, iso.y)
    iso = IsoLine('Q', 'H', 'P', 1.0, state)
    iso.calc_range(hr, pr); print(iso.x, iso.y)

    # bp = BasePlot(fluid_ref, graph_type, unit_system = 'KSI', **kwargs):
    bp = BasePlot('n-Pentane', 'PH', unit_system='EUR')
    # print(bp._get_sat_bounds('P'))
    # print(bp._get_iso_label(iso))
    print(bp.get_axis_limits())

        # get_update_pair(CoolProp.iP,CoolProp.iSmass,CoolProp.iT) -> (0,1,2,CoolProp.PSmass_INPUTS)
        # other values require switching and swapping
        # get_update_pair(CoolProp.iSmass,CoolProp.iP,CoolProp.iHmass) -> (1,0,2,CoolProp.PSmass_INPUTS)
