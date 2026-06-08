# AUTO-GENERATED — DO NOT EDIT BY HAND.
# Regenerate with: dev/stubs/gen_stubs.sh
# Source: wrappers/Python/CoolProp/CoolProp.pyx (which includes HumidAirProp.pyx
# and AbstractState.pyx) via stubgen-pyx + dev/stubs/postprocess.py.
#
# Enum-like arguments (parameters / phases / input_pairs / configuration_keys)
# are plain ``int`` at the Python boundary: pass the module-level constants from
# CoolProp.constants (e.g. PT_INPUTS, iT, iphase_liquid), which are integers.
from __future__ import annotations

from typing import Any, overload

from numpy import ndarray


class ChebyshevExpansion:

    def __init__(self, xmin: float, xmax: float, coef: list[float]):
        ...

    def xmin(self):
        ...

    def xmax(self):
        ...

    def coeff(self):
        ...

    def eval_many(self, x, y):
        ...

    def solve_for_x(self, y: float, a: float, b: float, bits: int, max_iter: int, boundstytol: float):
        ...

    def solve_for_x_many(self, y, a: float, b: float, bits: int, max_iter: int, boundstytol: float, x, counts):
        ...

class ChebyshevApproximation1D:

    def __init__(self, expansions):
        ...

    def xmin(self):
        ...

    def xmax(self):
        ...

    def is_monotonic(self):
        ...

    def eval_many(self, x, y):
        ...

    def get_x_for_y(self, y: float, bits: int, max_iter: int, boundstytol: float):
        ...

    def count_x_for_y_many(self, y, bits: int, max_iter: int, boundstytol: float, counts):
        ...

    def monotonic_intervals(self):
        ...

class SuperAncillary:

    def __init__(self, json_as_string: str):
        ...

    def eval_sat(self, T: float, prop: str, Q: int):
        ...

    def eval_sat_many(self, T, prop: str, Q: int, y):
        ...

class PyPhaseEnvelopeData:
    TypeI: bool
    iTsat_max: int
    ipsat_max: int
    icrit: int
    T: list
    p: list
    lnT: list
    lnp: list
    rhomolar_liq: list
    rhomolar_vap: list
    lnrhomolar_liq: list
    lnrhomolar_vap: list
    hmolar_liq: list
    hmolar_vap: list
    smolar_liq: list
    smolar_vap: list
    Q: list
    x: list
    y: list
    K: list

class PyCriticalState:
    T: float
    p: float
    rhomolar: float
    hmolar: float
    smolar: float
    stable: bool

class PyGuessesStructure:
    T: float
    p: float
    rhomolar: float
    hmolar: float
    smolar: float
    rhomolar_liq: float
    rhomolar_vap: float
    x: list
    y: list

    def __init__(self):
        ...

class PySpinodalData:
    tau: list[float]
    delta: list[float]
    M1: list[float]

class AbstractState:
    """
    This class is a one-to-one python wrapper of the :cpapi:`AbstractState` class
    """

    def __init__(self, backend: str, fluid: str):
        ...

    def fluid_param_string(self, key: str):
        """ Get a fluid parameter string - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fluid_param_string` """

    def name(self):
        """ Get the fluid name - wrapper of c++ function :cpapi:`CoolProp::AbstractState::name` """

    def backend_name(self):
        """ Get the backend name - wrapper of c++ function :cpapi:`CoolProp::AbstractState::backend_name` """

    def build_options_json(self):
        """ Canonical-JSON string of the options this instance was built with.

        Returns an empty string for backends that don't accept factory-string
        options.  For options-aware backends (e.g. SVDSBTL), the return value
        round-trips through :py:func:`AbstractState.factory`: feeding it back
        as the ``?<options>`` suffix on the factory string reproduces the
        construction byte-for-byte.

        Wrapper of c++ function :cpapi:`CoolProp::AbstractState::build_options_json` .
        """

    def fluid_names(self):
        """ Get the list of fluid names - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fluid_names` """

    def phase(self) -> int:
        """ Get the phase as key value- wrapper of c++ function :cpapi:`CoolProp::AbstractState::phase` """

    def specify_phase(self, phase: int):
        """ Specify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::specify_phase` """

    def unspecify_phase(self):
        """ Unspecify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::unspecify_phase` """

    def change_EOS(self, i: int, EOS_name: str):
        """ Change the EOS for one component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::change_EOS` """

    def apply_simple_mixing_rule(self, i: int, j: int, model: str):
        """ Apply a simple mixing rule - wrapper of c++ function :cpapi:`CoolProp::AbstractState::apply_simple_mixing_rule` """

    def set_binary_interaction_double(self, CAS1: str | int, CAS2: str | int, parameter: str, val: float):
        """ Set a double precision interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_binary_interaction_double` """

    def get_binary_interaction_double(self, CAS1: str | int, CAS2: str | int, parameter: str) -> float:
        """ Get a double precision interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_binary_interaction_double` """

    def set_binary_interaction_string(self, CAS1: str | int, CAS2: str | int, parameter: str, val: str):
        """ Set a string interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_binary_interaction_string` """

    def get_binary_interaction_string(self, CAS1: str, CAS2: str, parameter: str) -> str:
        """ Get a string interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_binary_interaction_string` """

    def set_cubic_alpha_C(self, i: int, parameter: str, c1: float, c2: float, c3: float):
        """ Set alpha function (MC or Twu) and fluid specific parameters - wrapper of c++ function :cpapi:`CoolProp::AbstractCubicBackend::set_cubic_alpha_C` """

    def set_fluid_parameter_double(self, i: int, parameter: str, val: float):
        """ Set a fluid parameter that is a double-precision number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_fluid_parameter_double` """

    def get_fluid_parameter_double(self, i: int, parameter: str) -> float:
        """ Get a fluid parameter that is a double-precision number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_fluid_parameter_double` """

    def update(self, ipair: int, Value1: float, Value2: float):
        """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """

    def update_QT_pure_superanc(self, Q: float, T: float):
        """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """

    def update_with_guesses(self, ipair: int, Value1: float, Value2: float, guesses: PyGuessesStructure):
        """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """

    def fast_evaluate(self, input_pair: int, val1, val2, outputs, out, status, imposed_phase: int=...):
        """
        Vectorized cache-bypassing batch evaluation - wrapper of c++ function
        :cpapi:`CoolProp::AbstractState::fast_evaluate`.

        This is a zero-allocation fast path supported on tabular backends
        (BICUBIC, TTSE) and the IF97 backend. All buffers are caller-allocated
        and validated for shape before the call:

        :param input_pair: One of ``HmolarP_INPUTS``, ``PT_INPUTS`` (tabular)
            or ``PT_INPUTS``/``HmolarP_INPUTS``/``HmassP_INPUTS`` (IF97).
        :param val1: 1D contiguous array of first inputs, length ``N_inputs``.
        :param val2: 1D contiguous array of second inputs, length ``N_inputs``.
        :param outputs: 1D int32 array of CoolProp output keys (e.g.
            ``[iDmolar, iHmolar]``), length ``N_outputs``.
        :param out: 2D contiguous (C-order) output array, shape
            ``(N_inputs, N_outputs)``. Populated in place.
        :param status: 1D int32 status array, length ``N_inputs``. ``0`` on
            success, nonzero ``fast_evaluate_status`` on per-point failure.
        :param imposed_phase: Optional phase hint (default: detect).

        .. note::
            Unlike :meth:`update`, this fast path does not perform saturation-
            curve cell bumping for tabular backends. PT_INPUTS points that
            land near the saturation curve may evaluate using a cell that
            straddles the two-phase notch, giving values that diverge from
            the cached path. To get strict single-phase results near
            saturation, either:

            * use ``HmolarP_INPUTS`` (which natively distinguishes liquid
              vs vapor by enthalpy), or
            * pass ``imposed_phase=iphase_liquid`` / ``iphase_gas``.
        """

    def set_mole_fractions(self, z: list[float]):
        """ Set the mole fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_mole_fractions` """

    def set_mass_fractions(self, z: list[float]):
        """ Set the mass fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_mass_fractions` """

    def set_volu_fractions(self, z: list[float]):
        """ Set the volume fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_volu_fractions` """

    def get_mole_fractions(self):
        """ Get the mole fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_mole_fractions` """

    def get_mass_fractions(self):
        """ Get the mass fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_mass_fractions` """

    def Tmin(self) -> float:
        """ Set the minimum temperature in K- wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmin` """

    def Tmax(self) -> float:
        """ Set the maximum temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmax` """

    def pmax(self) -> float:
        """ Set the maximum pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::pmax` """

    def Ttriple(self) -> float:
        """ Set the triple point temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Ttriple` """

    def T_critical(self) -> float:
        """ Gets the critical temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T_critical` """

    def rhomass_critical(self) -> float:
        """ Gets the critical density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass_critical` """

    def rhomolar_critical(self) -> float:
        """ Gets the critical density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar_critical` """

    def p_critical(self) -> float:
        """ Gets the critical pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::p_critical` """

    def all_critical_points(self) -> list:
        """ Calculate all the critical points - wrapper of c++ function :cpapi:`CoolProp::AbstractState::all_critical_points` """

    def criticality_contour_values(self) -> tuple:
        """
        Gets the criticality matrix values L1* and M1* - wrapper of c++ function :cpapi:`CoolProp::AbstractState::criticality_contour_values`
        Returns a tuple of (L1*, M1*)
        """

    def build_spinodal(self):
        """ Calculate the spinodal - wrapper of c++ function :cpapi:`CoolProp::AbstractState::build_spinodal` """

    def get_spinodal_data(self) -> PySpinodalData:
        """ Get the data from the spinodal - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_spinodal_data` """

    def T_reducing(self) -> float:
        """ Gets the reducing temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T_reducing` """

    def rhomolar_reducing(self) -> float:
        """ Gets the reducing density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar_reducing` """

    def rhomass_reducing(self) -> float:
        """ Gets the reducing density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass_reducing` """

    def tangent_plane_distance(self, T: float, p: float, w: list[float], rhomolar_guess: float=-1) -> float:
        """ Gets the tangent_plane_distance - wrapper of c++ function :cpapi:`CoolProp::AbstractState::tangent_plane_distance` """

    def get_fluid_constant(self, i: int, param: int) -> float:
        """ Get a constant for a fluid in the mixture :cpapi:`CoolProp::AbstractState::get_fluid_constant` """

    def keyed_output(self, iOutput: int) -> float:
        """ Get a keyed output :cpapi:`CoolProp::AbstractState::keyed_output(parameters key)` """

    def trivial_keyed_output(self, iOutput: int) -> float:
        """ Get a trivial keyed output not requiring any iteration :cpapi:`CoolProp::AbstractState::trivial_keyed_output(parameters key)` """

    def saturated_liquid_keyed_output(self, iOutput: int) -> float:
        """ Get a trivial output for the saturated liquid :cpapi:`CoolProp::AbstractState::saturated_liquid_keyed_output(parameters key)` """

    def saturated_vapor_keyed_output(self, iOutput: int) -> float:
        """ Get a trivial output for the saturated vapor :cpapi:`CoolProp::AbstractState::saturated_vapor_keyed_output(parameters key)` """

    def T(self) -> float:
        """ Get the temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T(void)` """

    def p(self) -> float:
        """ Get the pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::p(void)` """

    def compressibility_factor(self) -> float:
        """ Get the compressibility factor Z=p/(rho*R*T) - wrapper of c++ function :cpapi:`CoolProp::AbstractState::compressibility_factor(void)` """

    def Q(self) -> float:
        """ Get the vapor quality in mol/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Q(void)` """

    def Qmass(self) -> float:
        """ Get the vapor quality on a mass basis in kg/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Qmass(void)` """

    def rhomolar(self) -> float:
        """ Get the density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar(void)` """

    def rhomass(self) -> float:
        """ Get the density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass(void)` """

    def hmolar(self) -> float:
        """ Get the enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar(void)` """

    def hmass(self) -> float:
        """ Get the enthalpy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmass(void)` """

    def umolar(self) -> float:
        """ Get the internal energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umolar(void)` """

    def umass(self) -> float:
        """ Get the internal energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umass(void)` """

    def smolar(self) -> float:
        """ Get the entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar(void)` """

    def smass(self) -> float:
        """ Get the entropy in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smass(void)` """

    def cpmolar(self) -> float:
        """ Get the constant pressure specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cpmolar(void)` """

    def cpmass(self) -> float:
        """ Get the constant pressure specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cpmass(void)` """

    def cp0molar(self) -> float:
        """ Get the ideal gas constant pressure specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cp0molar(void)` """

    def cp0mass(self) -> float:
        """ Get the ideal gas constant pressure specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cp0mass(void)` """

    def cvmolar(self) -> float:
        """ Get the constant volume specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cvmolar(void)` """

    def cvmass(self) -> float:
        """ Get the constant volume specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cvmass(void)` """

    def gibbsmass(self) -> float:
        """ Get the mass-specific Gibbs energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmass(void)` """

    def gibbsmolar(self) -> float:
        """ Get the mole-specific Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar(void)` """

    def helmholtzmass(self) -> float:
        """ Get the mass-specific Helmholtz energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmass(void)` """

    def helmholtzmolar(self) -> float:
        """ Get the mole-specific Helmholtz energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmolar(void)` """

    def tau(self) -> float:
        """ Get the reciprocal reduced temperature - wrapper of c++ function :cpapi:`CoolProp::AbstractState::tau(void)` """

    def delta(self) -> float:
        """ Get the reduced density - wrapper of c++ function :cpapi:`CoolProp::AbstractState::delta(void)` """

    def speed_sound(self) -> float:
        """ Get the speed of sound in m/s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::speed_sound(void)` """

    def molar_mass(self) -> float:
        """ Get the molar mass in kg/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::molar_mass(void)` """

    def acentric_factor(self) -> float:
        """ Get the acentric factor - wrapper of c++ function :cpapi:`CoolProp::AbstractState::acentric_factor(void)` """

    def gas_constant(self) -> float:
        """ Get the gas constant in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gas_constant(void)` """

    def viscosity(self) -> float:
        """ Get the viscosity in Pa-s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::viscosity(void)` """

    def conductivity(self) -> float:
        """ Get the thermal conductivity in W/m/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conductivity(void)` """

    def surface_tension(self) -> float:
        """ Get the surface tension N/m - wrapper of c++ function :cpapi:`CoolProp::AbstractState::surface_tension(void)` """

    def Prandtl(self) -> float:
        """ Get the Prandtl number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Prandtl(void)` """

    def Bvirial(self) -> float:
        """ Get the B virial coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Bvirial(void)` """

    def Cvirial(self) -> float:
        """ Get the C virial coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Cvirial(void)` """

    def fundamental_derivative_of_gas_dynamics(self) -> float:
        """ Get the fundamental derivative of gas dynamics - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fundamental_derivative_of_gas_dynamics(void)` """

    def PIP(self) -> float:
        """ Get the phase identification parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::PIP` """

    def isobaric_expansion_coefficient(self) -> float:
        """ Get the isobaric expansion coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::isobaric_expansion_coefficient(void)` """

    def isothermal_compressibility(self) -> float:
        """ Get the isothermal_compressibility - wrapper of c++ function :cpapi:`CoolProp::AbstractState::isothermal_compressibility(void)` """

    def fugacity(self, i: int) -> float:
        """ Get the fugacity of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fugacity(std::size_t)` """

    def fugacity_coefficient(self, i: int) -> float:
        """ Get the fugacity coefficient of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fugacity_coefficient(std::size_t)` """

    def chemical_potential(self, i: int) -> float:
        """ Get the chemical potential of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::chemical_potential(std::size_t)` """

    def mole_fractions_liquid(self):
        """ Get the mole fractions of the liquid phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_liquid(void)` """

    def mole_fractions_vapor(self):
        """ Get the mole fractions of the vapor phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_vapor(void)` """

    def true_critical_point(self) -> tuple:
        """ Get the "true" critical point where dp/drho|T = 0 & d2p/drho^2|T = 0 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::true_critical_point` """

    def conformal_state(self, reference_fluid: str, T: float, rho: float) -> dict:
        """ Solve for conformal state used in extended corresponding states - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conformal_state` """

    def conductivity_contributions(self) -> dict:
        """ Retrieve each of the contributions to the conductivity, each in W/m/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conductivity_contributions` """

    def viscosity_contributions(self) -> dict:
        """ Retrieve each of the contributions to the viscosity, each in Pa-s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::viscosity_contributions` """

    def helmholtzmolar_excess(self) -> float:
        """ Get the mole-specific excess Helmholtz energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmolar_excess(void)` """

    def helmholtzmass_excess(self) -> float:
        """ Get the mass-specific excess Helmholtz energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmass_excess(void)` """

    def gibbsmolar_excess(self) -> float:
        """ Get the mole-specific excess Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar_excess(void)` """

    def gibbsmass_excess(self) -> float:
        """ Get the mass-specific excess Gibbs energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmass_excess(void)` """

    def umolar_excess(self) -> float:
        """ Get the mole-specific excess internal energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umolar_excess(void)` """

    def umass_excess(self) -> float:
        """ Get the mass-specific excess internal energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umass_excess(void)` """

    def hmolar_excess(self) -> float:
        """ Get the mole-specific excess enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar_excess(void)` """

    def hmass_excess(self) -> float:
        """ Get the mass-specific excess enthalpy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmass_excess(void)` """

    def smolar_excess(self) -> float:
        """ Get the mole-specific excess entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar_excess(void)` """

    def smass_excess(self) -> float:
        """ Get the mass-specific excess entropy in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smass_excess(void)` """

    def volumemolar_excess(self) -> float:
        """ Get the mole-specific excess volume in m^3/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::volumemolar_excess(void)` """

    def volumemass_excess(self) -> float:
        """ Get the mass-specific excess volume in m^3/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::volumemass_excess(void)` """

    def gibbsmolar_residual(self) -> float:
        """ Get the mole-specific residual Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar_residual(void)` """

    def hmolar_residual(self) -> float:
        """ Get the mole-specific residual enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar_residual(void)` """

    def smolar_residual(self) -> float:
        """ Get the mole-specific residual entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar_residual(void)` """

    def neff(self) -> float:
        """ Get the effective hardness of interaction - wrapper of c++ function :cpapi:`CoolProp::AbstractState::neff(void)` """

    def umolar_idealgas(self) -> float:
        """ Get the mole-specific ideal gas internal energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umolar_idealgas(void)` """

    def umass_idealgas(self) -> float:
        """ Get the mass-specific ideal gas internal energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umass_idealgas(void)` """

    def hmolar_idealgas(self) -> float:
        """ Get the mole-specific ideal gas enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar_idealgas(void)` """

    def hmass_idealgas(self) -> float:
        """ Get the mass-specific ideal gas enthalpy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmass_idealgas(void)` """

    def smolar_idealgas(self) -> float:
        """ Get the mole-specific ideal gas entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar_idealgas(void)` """

    def smass_idealgas(self) -> float:
        """ Get the mass-specific ideal gas entropy in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smass_idealgas(void)` """

    def first_partial_deriv(self, OF: int, WRT: int, CONSTANT: int) -> float:
        """ Get the first partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_partial_deriv` """

    def second_partial_deriv(self, OF: int, WRT1: int, CONSTANT1: int, WRT2: int, CONSTANT2: int) -> float:
        """ Get the second partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_partial_deriv` """

    def first_saturation_deriv(self, OF: int, WRT: int) -> float:
        """ Get the first derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_saturation_deriv` """

    def second_saturation_deriv(self, OF1: int, WRT1: int, WRT2: int) -> float:
        """ Get the second derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_saturation_deriv` """

    def first_two_phase_deriv(self, Of: int, Wrt: int, Constant: int) -> float:
        """ Get the first two-phase derivative - wrapper of C++ function :cpapi:`CoolProp::AbstractState::first_two_phase_deriv` """

    def second_two_phase_deriv(self, Of1: int, Wrt1: int, Constant1: int, Wrt2: int, Constant2: int) -> float:
        """ Get the second two-phase derivative - wrapper of C++ function :cpapi:`CoolProp::AbstractState::second_two_phase_deriv` """

    def first_two_phase_deriv_splined(self, Of: int, Wrt: int, Constant: int, x_end: float) -> float:
        """ Get the first two-phase derivative using splines - wrapper of C++ function :cpapi:`CoolProp::AbstractState::first_two_phase_deriv_splined` """

    def has_melting_line(self) -> bool:
        """ Check if the fluid has a melting line - True if is does, False otherwise - wrapper of c++ function :cpapi:`CoolProp::AbstractState::has_melting_line` """

    def melting_line(self, param: int, given: int, value: float) -> float:
        """ Get values from the melting line - wrapper of c++ function :cpapi:`CoolProp::AbstractState::melting_line` """

    def saturation_ancillary(self, param: int, Q: int, given: int, value: float) -> float:
        """ Get values from the saturation_ancillary - wrapper of c++ function :cpapi:`CoolProp::AbstractState::saturation_ancillary` """

    def build_phase_envelope(self, type: str):
        """ Build the phase envelope - wrapper of c++ function :cpapi:`CoolProp::AbstractState::build_phase_envelope` """

    def get_phase_envelope_data(self) -> PyPhaseEnvelopeData:
        """ Get the phase envelope data - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_phase_envelope_data` """

    def ideal_curve(self, type: str) -> tuple:
        """ Get an ideal curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::ideal_curve` """

    def alpha0(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alpha0` """

    def dalpha0_dDelta(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dDelta` """

    def dalpha0_dTau(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dTau` """

    def d2alpha0_dDelta2(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta2` """

    def d2alpha0_dDelta_dTau(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta_dTau` """

    def d2alpha0_dTau2(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dTau2` """

    def d3alpha0_dTau3(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dTau3` """

    def d3alpha0_dDelta_dTau2(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta_dTau2` """

    def d3alpha0_dDelta2_dTau(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta2_dTau` """

    def d3alpha0_dDelta3(self) -> float:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta3` """

    def alphar(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alphar` """

    def dalphar_dDelta(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dDelta` """

    def dalphar_dTau(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dTau` """

    def d2alphar_dDelta2(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta2` """

    def d2alphar_dDelta_dTau(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta_dTau` """

    def d2alphar_dTau2(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dTau2` """

    def d3alphar_dTau3(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dTau3` """

    def d3alphar_dDelta_dTau2(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta_dTau2` """

    def d3alphar_dDelta2_dTau(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta2_dTau` """

    def d3alphar_dDelta3(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta3` """

    def d4alphar_dTau4(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dTau4` """

    def d4alphar_dDelta_dTau3(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta_dTau3` """

    def d4alphar_dDelta2_dTau2(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta2_dTau2` """

    def d4alphar_dDelta3_dTau(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta3_dTau` """

    def d4alphar_dDelta4(self) -> float:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta4` """

class State:
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
    pAS: AbstractState

    def __init__(self, _Fluid: object, StateDict: dict, phase: object=None, backend=None):
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

    def set_Fluid(self, Fluid, backend):
        ...

    def update_ph(self, p: float, h: float):
        """
        Use the pressure and enthalpy directly

        Parameters
        ----------
        p: float
            Pressure (absolute) [kPa]
        h: float
            Enthalpy [kJ/kg]
        """

    def update_Trho(self, T: float, rho: float):
        """
        Just use the temperature and density directly for speed

        Parameters
        ----------
        T: float
            Temperature [K]
        rho: float
            Density [kg/m^3]
        """

    def update(self, params: dict):
        """
        Parameters
        params, dictionary
            A dictionary of terms to be updated, with keys equal to single-char inputs to the Props function,
            for instance ``dict(T=298, P = 101.325)`` would be one standard atmosphere
        """

    def Phase(self) -> int:
        """
        Returns an integer flag for the phase of the fluid, where the flag value
        is one of iLiquid, iSupercritical, iGas, iTwoPhase

        These constants are defined in the phase_constants module, and are imported
        into this module
        """

    def Props(self, iOutput: int) -> float:
        ...

    def get_Q(self) -> float:
        """ Get the quality [-] """

    def get_MM(self) -> float:
        """ Get the mole mass [kg/kmol] or [g/mol] """

    def get_rho(self) -> float:
        """ Get the density [kg/m^3] """

    def get_p(self) -> float:
        """ Get the pressure [kPa] """

    def get_T(self) -> float:
        """ Get the temperature [K] """

    def get_h(self) -> float:
        """ Get the specific enthalpy [kJ/kg] """

    def get_u(self) -> float:
        """ Get the specific internal energy [kJ/kg] """

    def get_s(self) -> float:
        """ Get the specific enthalpy [kJ/kg/K] """

    def get_cp0(self) -> float:
        """ Get the specific heat at constant pressure for the ideal gas [kJ/kg/K] """

    def get_cp(self) -> float:
        """ Get the specific heat at constant pressure  [kJ/kg/K] """

    def get_cv(self) -> float:
        """ Get the specific heat at constant volume  [kJ/kg/K] """

    def get_speed_sound(self) -> float:
        """ Get the speed of sound  [m/s] """

    def get_visc(self) -> float:
        """ Get the viscosity, in [Pa-s]"""

    def get_cond(self) -> float:
        """ Get the thermal conductivity, in [kW/m/K]"""

    def get_Tsat(self, Q: float=1):
        """
        Get the saturation temperature, in [K]

        Returns ``None`` if pressure is not within the two-phase pressure range
        """

    def get_superheat(self):
        """
        Get the amount of superheat above the saturation temperature corresponding to the pressure, in [K]

        Returns ``None`` if pressure is not within the two-phase pressure range
        """

    def get_subcooling(self):
        """
        Get the amount of subcooling below the saturation temperature corresponding to the pressure, in [K]

        Returns ``None`` if pressure is not within the two-phase pressure range
        """

    def get_dpdT(self) -> float:
        ...

    def speed_test(self, N: int):
        ...

    def __str__(self):
        """
        Return a string representation of the state
        """

    def copy(self) -> State:
        """
        Make a copy of this State class
        """

def HAProps_Aux(OutputName: str, T: float, p: float, w: float) -> tuple:
    """
    Allows low-level access to some of the routines employed in HumidAirProps

    Returns tuples of the form ``(Value, Units)`` where ``Value`` is the actual value and ``Units`` is a string that describes the units

    The list of possible inputs is

    * Baa [First virial air-air coefficient]
    * Caaa [Second virial air coefficient]
    * Bww [First virial water-water coefficient]
    * Cwww [Second virial water coefficient]
    * Baw [First cross virial coefficient]
    * Caww [Second air-water-water virial coefficient]
    * Caaw [Second air-air-water virial coefficient]
    * beta_H 
    * kT
    * vbar_ws [Molar saturated volume of water vapor]
    * p_ws [Saturated vapor pressure of pure water (>=0.01C) or ice (<0.01 C)]
    * f [Enhancement factor]
    """

def cair_sat(T: float) -> float:
    """
    The derivative of the saturation enthalpy cair_sat = d(hsat)/dT
    """

def set_reference_state(FluidName, *args):
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

def generate_update_pair(key1: int, value1: float, key2: int, value2: float) -> tuple:
    """
    This function will generate an input pair to the update() function given the key, value pairs for both inputs
    """

def get_config_as_json_string() -> str:
    """
    Obtain a json formulation of the internal configuration in CoolProp

    Values can be set by passing a modified json library (converted to string) to set_config_as_json_string
    """

def config_key_description(key: str) -> str:
    """
    Obtain the string description for a configuration key.  Python wrapper of C++ function :cpapi:`CoolProp::config_key_description`
    """

def set_config_as_json_string(s: str):
    """
    Set the internal configuration in CoolProp from a json data string

    Current state can be obtained by calling get_config_as_json_string
    """

def set_config_double(key: int, value: float):
    """ Set configuration key that is a double-precision float;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_double` """

def set_config_string(key: int, value: str):
    """ Set a configuration key that is a string;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_string` """

def set_config_bool(key: int, value: bool):
    """ Set a configuration key that is a boolean;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_bool` """

def set_config_int(key: int, value: int):
    """ Set a configuration key that is an integer;  wrapper of wrapper of C++ function :cpapi:`CoolProp::set_config_int` """

def get_config_double(key: int) -> float:
    """ Get a configuration key that is a double-precision float;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_double` """

def get_config_string(key: int) -> str:
    """ Get a configuration key that is a string;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_string` """

def get_config_bool(key: int) -> bool:
    """ Get a configuration key that is a boolean;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_bool` """

def get_config_int(key: int) -> int:
    """ Get a configuration key that is an integer;  wrapper of wrapper of C++ function :cpapi:`CoolProp::get_config_int` """

def get_parameter_index(key: str) -> int:
    ...

def get_phase_index(key: str) -> int:
    ...

def get_parameter_information(key: int, info: str) -> str:
    ...

def get_mixture_binary_pair_data(CAS1, CAS2, key) -> str:
    """
    Obtain mixture interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::get_mixture_binary_pair_data`
    """

def set_mixture_binary_pair_data(CAS1, CAS2, key, val):
    """
    Set mixture interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::set_mixture_binary_pair_data`
    """

def get_mixture_binary_pair_pcsaft(CAS1, CAS2, key) -> str:
    """
    Obtain mixture PC-SAFT interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::get_mixture_binary_pair_pcsaft`
    """

def set_mixture_binary_pair_pcsaft(CAS1, CAS2, key, val):
    """
    Set mixture PC-SAFT interaction parameter.  Python wrapper of C++ function :cpapi:`CoolProp::set_mixture_binary_pair_pcsaft`
    """

def add_fluids_as_JSON(backend, JSONstring):
    """
    Add fluids in a JSON-formatted string format. Python wrapper of C++ function :cpapi:`CoolProp::add_fluids_as_JSON`
    """

def get_global_param_string(param):
    ...

def is_trivial_parameter(key: int):
    ...

def get_fluid_param_string(fluid, param):
    ...

def apply_simple_mixing_rule(CAS1, CAS2, rule):
    """
    Apply simple mixing rule.  Currently linear or Lorentz-Berthelot.  Python wrapper of C++ function :cpapi:`CoolProp::apply_simple_mixing_rule`
    """

def set_departure_functions(functions):
    """
    Specify the departure terms as JSON. Python wrapper of C++ function :cpapi:`CoolProp::set_departure_functions`
    """

def set_interaction_parameters(data):
    """
    Specify the binary interaction terms as JSON. Python wrapper of C++ function :cpapi:`CoolProp::set_interaction_parameters`
    """

def set_predefined_mixtures(data):
    """
    Specify predefined mixtures as JSON. Python wrapper of C++ function :cpapi:`CoolProp::set_predefined_mixtures`
    """

def saturation_ancillary(name: str, output: str, Q: int, input: str, value: float) -> float:
    """
    Return a value from the saturation ancillary equations; python wrapper of :cpapi:`CoolProp::saturation_ancillary`
    """

def FluidsList() -> list:
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

def get_aliases(Fluid):
    """
    Return a list of aliases for the given fluid.
    Uses regex to properly handle cases like '1,2-dichloroethane' where the comma
    is part of the chemical name rather than a separator.
    """

def get_REFPROPname(Fluid) -> str:
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

def get_BibTeXKey(Fluid, key) -> str:
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

def get_errstr() -> str:
    """
    Return the current error string
    """

def set_debug_level(level: int):
    """
    Set the current debug level as integer in the range [0,10]

    Parameters
    ----------
    level : int
        If level is 0, no output will be written to screen, if >0,
        some output will be written to screen.  The larger level is,
        the more verbose the output will be
    """

def get_debug_level() -> int:
    """
    Return the current debug level as integer

    Returns
    -------
    level : int
        If level is 0, no output will be written to screen, if >0,
        some output will be written to screen.  The larger level is,
        the more verbose the output will be
    """

def extract_backend(in_str):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::extract_backend` .
    """

def extract_fractions(flds):
    """
    A Python wrapper of C++ function :cpapi:`CoolProp::extract_fractions` .
    """

def rebuildState(d):
    ...


# --- hand-written overloads (do not regenerate; spliced by postprocess.py) ---
#
# The high-level property functions are typed by hand because stubgen cannot
# express their scalar/array polymorphism (overloads are a stub-only construct,
# and Cython annotation_typing would coerce the args at runtime).
#
# NOTE on overload overlap: PropsSI/HAPropsSI/HAProps return ndarray iff ANY
# numeric input is array-like — which cannot be expressed without the scalar
# overload (-> float) and the array overload (whose `_Num` includes float, to
# admit mixed scalar+array calls -> ndarray) overlapping on an all-float call.
# Ordering (scalar first) makes the all-scalar case resolve to float under both
# pyright and mypy; the array overload carries `# type: ignore[overload-overlap]`
# so mypy does not flag the unavoidable overlap.  Per a .pyi stub, NO
# implementation line follows the overloads (that is itself a mypy error).

from typing import Sequence

_Scalar = float
_Array = Sequence[float] | ndarray
_Num = _Scalar | _Array


# PropsSI: two-argument trivial-property lookup -> float; six-argument form is
# scalar->float or array->ndarray.  (in1/Output and in3/in5 are required for the
# six-arg form; the two-arg form is _Props1SI(Output, FluidName).)
@overload
def PropsSI(Output: str, FluidName: str) -> float: ...
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> float: ...  # type: ignore[overload-overlap]
@overload
def PropsSI(Output: str, Name1: str, Prop1: _Num, Name2: str, Prop2: _Num, FluidName: str) -> ndarray: ...


# PhaseSI does NOT vectorize (per its docstring) — scalar inputs, returns a
# single phase string.
def PhaseSI(Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> str: ...


# Props is the deprecated scalar-only predecessor of PropsSI (no array path).
def Props(Output: str, Name1: str, Prop1: _Scalar, Name2: str, Prop2: _Scalar, FluidName: str) -> float: ...


# HAPropsSI / HAProps DO vectorize (scalar->float, array->ndarray).
@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Scalar, Name2: str, Value2: _Scalar, Name3: str, Value3: _Scalar) -> float: ...  # type: ignore[overload-overlap]
@overload
def HAPropsSI(Output: str, Name1: str, Value1: _Num, Name2: str, Value2: _Num, Name3: str, Value3: _Num) -> ndarray: ...


@overload
def HAProps(Output: str, Name1: str, Value1: _Scalar, Name2: str, Value2: _Scalar, Name3: str, Value3: _Scalar) -> float: ...  # type: ignore[overload-overlap]
@overload
def HAProps(Output: str, Name1: str, Value1: _Num, Name2: str, Value2: _Num, Name3: str, Value3: _Num) -> ndarray: ...
