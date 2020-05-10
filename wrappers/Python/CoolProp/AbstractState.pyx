# This file is embedded directly in CoolProp.pyx

from . cimport constants_header

cdef class PyPhaseEnvelopeData:
    pass

cdef class PyCriticalState:
    pass

cdef class PyGuessesStructure:
    def __init__(self):
        self.T = get_HUGE()
        self.p = get_HUGE()
        self.rhomolar = get_HUGE()
        self.hmolar = get_HUGE()
        self.smolar = get_HUGE()
        self.rhomolar_liq = get_HUGE()
        self.rhomolar_vap = get_HUGE()
        self.x = []
        self.y = []

cdef class PySpinodalData:
    pass

cdef class AbstractState:
    """
    This class is a one-to-one python wrapper of the :cpapi:`AbstractState` class
    """

    def __cinit__(self, string backend, string fluid):
        self.thisptr = cAbstractState.factory(backend, fluid)

    def __dealloc__(self):
        del self.thisptr

    cpdef fluid_param_string(self, string key):
        """ Get a fluid parameter string - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fluid_param_string` """
        return self.thisptr.fluid_param_string(key)

    cpdef name(self):
        """ Get the fluid name - wrapper of c++ function :cpapi:`CoolProp::AbstractState::name` """
        return self.thisptr.name()
    cpdef backend_name(self):
        """ Get the backend name - wrapper of c++ function :cpapi:`CoolProp::AbstractState::backend_name` """
        return self.thisptr.backend_name()
    cpdef fluid_names(self):
        """ Get the list of fluid names - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fluid_names` """
        return self.thisptr.fluid_names()

    cpdef constants_header.phases phase(self) except *:
        """ Get the phase as key value- wrapper of c++ function :cpapi:`CoolProp::AbstractState::phase` """
        return self.thisptr.phase()

    cpdef specify_phase(self, constants_header.phases phase):
        """ Specify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::specify_phase` """
        self.thisptr.specify_phase(phase)
    cpdef unspecify_phase(self):
        """ Unspecify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::unspecify_phase` """
        self.thisptr.unspecify_phase()

    cpdef change_EOS(self, size_t i, string EOS_name):
        """ Change the EOS for one component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::change_EOS` """
        self.thisptr.change_EOS(i, EOS_name)

    cpdef apply_simple_mixing_rule(self, size_t i, size_t j, string model):
        """ Apply a simple mixing rule - wrapper of c++ function :cpapi:`CoolProp::AbstractState::apply_simple_mixing_rule` """
        self.thisptr.apply_simple_mixing_rule(i, j, model)

    cpdef set_binary_interaction_double(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter, double val):
        """ Set a double precision interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_binary_interaction_double` """
        if string_or_size_t in cython.integral:
            self.thisptr.set_binary_interaction_double(<size_t>CAS1, <size_t>CAS2, parameter, val)
        else:
            self.thisptr.set_binary_interaction_double(<string>CAS1, <string>CAS2, parameter, val)
    cpdef double get_binary_interaction_double(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter) except *:
        """ Get a double precision interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_binary_interaction_double` """
        if string_or_size_t in cython.integral:
            return self.thisptr.get_binary_interaction_double(<size_t>CAS1, <size_t>CAS2, parameter)
        else:
            return self.thisptr.get_binary_interaction_double(<string>CAS1, <string>CAS2, parameter)

    cpdef set_binary_interaction_string(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter, string val):
        """ Set a string interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_binary_interaction_string` """
        if string_or_size_t in cython.integral:
            self.thisptr.set_binary_interaction_string(<size_t>CAS1, <size_t>CAS2, parameter, val)
        else:
            self.thisptr.set_binary_interaction_string(<string>CAS1, <string>CAS2, parameter, val)

    cpdef string get_binary_interaction_string(self, string CAS1, string CAS2, string parameter) except *:
        """ Get a string interaction parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_binary_interaction_string` """
        return self.thisptr.get_binary_interaction_string(CAS1, CAS2, parameter)

    cpdef set_fluid_parameter_double(self, size_t i, string parameter, double val):
        """ Set a fluid parameter that is a double-precision number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_fluid_parameter_double` """
        self.thisptr.set_fluid_parameter_double(i, parameter, val)

    cpdef double get_fluid_parameter_double(self, size_t i, string parameter) except *:
        """ Get a fluid parameter that is a double-precision number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_fluid_parameter_double` """
        return self.thisptr.get_fluid_parameter_double(i, parameter)

    cpdef update(self, constants_header.input_pairs ipair, double Value1, double Value2):
        """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """
        self.thisptr.update(ipair, Value1, Value2)
    cpdef update_with_guesses(self, constants_header.input_pairs ipair, double Value1, double Value2, PyGuessesStructure guesses):
        """ Update function - wrapper of c++ function :cpapi:`CoolProp::AbstractState::update` """
        cdef cAbstractState.GuessesStructure _guesses
        _guesses.T = guesses.T
        _guesses.p = guesses.p
        _guesses.rhomolar = guesses.rhomolar
        _guesses.rhomolar_liq = guesses.rhomolar_liq
        _guesses.rhomolar_vap = guesses.rhomolar_vap
        _guesses.x = guesses.x
        _guesses.y = guesses.y
        self.thisptr.update_with_guesses(ipair, Value1, Value2, _guesses)

    cpdef set_mole_fractions(self, vector[double] z):
        """ Set the mole fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_mole_fractions` """
        self.thisptr.set_mole_fractions(z)
    cpdef set_mass_fractions(self, vector[double] z):
        """ Set the mass fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_mass_fractions` """
        self.thisptr.set_mass_fractions(z)
    cpdef set_volu_fractions(self, vector[double] z):
        """ Set the volume fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::set_volu_fractions` """
        self.thisptr.set_volu_fractions(z)
    cpdef get_mole_fractions(self):
        """ Get the mole fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_mole_fractions` """
        return self.thisptr.get_mole_fractions()
    cpdef get_mass_fractions(self):
        """ Get the mass fractions - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_mass_fractions` """
        return self.thisptr.get_mass_fractions()

    ## ----------------------------------------
    ##        Limits
    ## ----------------------------------------
    cpdef double Tmin(self) except *:
        """ Set the minimum temperature in K- wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmin` """
        return self.thisptr.Tmin()
    cpdef double Tmax(self) except *:
        """ Set the maximum temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmax` """
        return self.thisptr.Tmax()
    cpdef double pmax(self) except *:
        """ Set the maximum pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::pmax` """
        return self.thisptr.pmax()
    cpdef double Ttriple(self) except *:
        """ Set the triple point temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Ttriple` """
        return self.thisptr.Ttriple()


    ## Critical point
    cpdef double T_critical(self) except *:
        """ Gets the critical temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T_critical` """
        return self.thisptr.T_critical()
    cpdef double rhomass_critical(self) except *:
        """ Gets the critical density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass_critical` """
        return self.thisptr.rhomass_critical()
    cpdef double rhomolar_critical(self) except *:
        """ Gets the critical density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar_critical` """
        return self.thisptr.rhomolar_critical()
    cpdef double p_critical(self) except *:
        """ Gets the critical pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::p_critical` """
        return self.thisptr.p_critical()

    cpdef list all_critical_points(self):
        """ Calculate all the critical points - wrapper of c++ function :cpapi:`CoolProp::AbstractState::all_critical_points` """
        # Get all the critical points
        cdef vector[cAbstractState.CriticalState] critpts = self.thisptr.all_critical_points()
        cdef PyCriticalState pypt
        cdef list collection = []
        # Convert to python
        for pt in critpts:
            pypt = PyCriticalState()
            pypt.stable = pt.stable
            pypt.T = pt.T
            pypt.p = pt.p
            pypt.rhomolar = pt.rhomolar
            collection.append(pypt)
        return collection
    cpdef tuple criticality_contour_values(self):
        """
        Gets the criticality matrix values L1* and M1* - wrapper of c++ function :cpapi:`CoolProp::AbstractState::criticality_contour_values`
        Returns a tuple of (L1*, M1*)
        """
        cdef CoolPropDbl L1star = 0, M1star = 0
        self.thisptr.criticality_contour_values(L1star, M1star)
        return L1star, M1star

    cpdef build_spinodal(self):
        """ Calculate the spinodal - wrapper of c++ function :cpapi:`CoolProp::AbstractState::build_spinodal` """
        self.thisptr.build_spinodal()
    cpdef PySpinodalData get_spinodal_data(self):
        """ Get the data from the spinodal - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_spinodal_data` """
        cdef cAbstractState.SpinodalData data = self.thisptr.get_spinodal_data()
        cdef PySpinodalData out = PySpinodalData()
        out.tau = data.tau
        out.delta = data.delta
        out.M1 = data.M1
        return out

    ## Reducing point
    cpdef double T_reducing(self) except *:
        """ Gets the reducing temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T_reducing` """
        return self.thisptr.T_reducing()
    cpdef double rhomolar_reducing(self) except *:
        """ Gets the reducing density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar_reducing` """
        return self.thisptr.rhomolar_reducing()
    cpdef double rhomass_reducing(self) except *:
        """ Gets the reducing density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass_reducing` """
        return self.thisptr.rhomass_reducing()


    cpdef double tangent_plane_distance(self, double T, double p, vector[double] w, double rhomolar_guess = -1) except *:
        """ Gets the tangent_plane_distance - wrapper of c++ function :cpapi:`CoolProp::AbstractState::tangent_plane_distance` """
        return self.thisptr.tangent_plane_distance(T, p, w, rhomolar_guess)

    ## ----------------------------------------
    ##        Fluid property accessors
    ## ----------------------------------------

    cpdef double get_fluid_constant(self, size_t i,constants_header.parameters param) except *:
        """ Get a constant for a fluid in the mixture :cpapi:`CoolProp::AbstractState::get_fluid_constant` """
        return self.thisptr.get_fluid_constant(i, param)

    cpdef double keyed_output(self, parameters iOutput) except *:
        """ Get a keyed output :cpapi:`CoolProp::AbstractState::keyed_output(parameters key)` """
        return self.thisptr.keyed_output(iOutput)
    cpdef double trivial_keyed_output(self, parameters iOutput) except *:
        """ Get a trivial keyed output not requiring any iteration :cpapi:`CoolProp::AbstractState::trivial_keyed_output(parameters key)` """
        return self.thisptr.trivial_keyed_output(iOutput)
    cpdef double saturated_liquid_keyed_output(self, parameters iOutput) except *:
        """ Get a trivial output for the saturated liquid :cpapi:`CoolProp::AbstractState::saturated_liquid_keyed_output(parameters key)` """
        return self.thisptr.saturated_liquid_keyed_output(iOutput)
    cpdef double saturated_vapor_keyed_output(self, parameters iOutput) except *:
        """ Get a trivial output for the saturated vapor :cpapi:`CoolProp::AbstractState::saturated_vapor_keyed_output(parameters key)` """
        return self.thisptr.saturated_vapor_keyed_output(iOutput)

    cpdef double T(self) except *:
        """ Get the temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T(void)` """
        return self.thisptr.T()
    cpdef double p(self) except *:
        """ Get the pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::p(void)` """
        return self.thisptr.p()
    cpdef double compressibility_factor(self) except *:
        """ Get the compressibility factor Z=p/(rho*R*T) - wrapper of c++ function :cpapi:`CoolProp::AbstractState::compressibility_factor(void)` """
        return self.thisptr.compressibility_factor()
    cpdef double Q(self) except *:
        """ Get the vapor quality in mol/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Q(void)` """
        return self.thisptr.Q()
    cpdef double rhomolar(self) except *:
        """ Get the density in mol/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomolar(void)` """
        return self.thisptr.rhomolar()
    cpdef double rhomass(self) except *:
        """ Get the density in kg/m^3 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::rhomass(void)` """
        return self.thisptr.rhomass()
    cpdef double hmolar(self) except *:
        """ Get the enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar(void)` """
        return self.thisptr.hmolar()
    cpdef double hmass(self) except *:
        """ Get the enthalpy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmass(void)` """
        return self.thisptr.hmass()
    cpdef double umolar(self) except *:
        """ Get the internal energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umolar(void)` """
        return self.thisptr.umolar()
    cpdef double umass(self) except *:
        """ Get the internal energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umass(void)` """
        return self.thisptr.umass()
    cpdef double smolar(self) except *:
        """ Get the entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar(void)` """
        return self.thisptr.smolar()
    cpdef double smass(self) except *:
        """ Get the entropy in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smass(void)` """
        return self.thisptr.smass()
    cpdef double cpmolar(self) except *:
        """ Get the constant pressure specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cpmolar(void)` """
        return self.thisptr.cpmolar()
    cpdef double cpmass(self) except *:
        """ Get the constant pressure specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cpmass(void)` """
        return self.thisptr.cpmass()
    cpdef double cp0molar(self) except *:
        """ Get the ideal gas constant pressure specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cp0molar(void)` """
        return self.thisptr.cp0molar()
    cpdef double cp0mass(self) except *:
        """ Get the ideal gas constant pressure specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cp0mass(void)` """
        return self.thisptr.cp0mass()
    cpdef double cvmolar(self) except *:
        """ Get the constant volume specific heat in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cvmolar(void)` """
        return self.thisptr.cvmolar()
    cpdef double cvmass(self) except *:
        """ Get the constant volume specific heat in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::cvmass(void)` """
        return self.thisptr.cvmass()
    cpdef double gibbsmass(self) except *:
        """ Get the mass-specific Gibbs energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmass(void)` """
        return self.thisptr.gibbsmass()
    cpdef double gibbsmolar(self) except *:
        """ Get the mole-specific Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar(void)` """
        return self.thisptr.gibbsmolar()
    cpdef double helmholtzmass(self) except *:
        """ Get the mass-specific Helmholtz energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmass(void)` """
        return self.thisptr.helmholtzmass()
    cpdef double helmholtzmolar(self) except *:
        """ Get the mole-specific Helmholtz energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmolar(void)` """
        return self.thisptr.helmholtzmolar()
    cpdef double tau(self) except *:
        """ Get the reciprocal reduced temperature - wrapper of c++ function :cpapi:`CoolProp::AbstractState::tau(void)` """
        return self.thisptr.tau()
    cpdef double delta(self) except *:
        """ Get the reduced density - wrapper of c++ function :cpapi:`CoolProp::AbstractState::delta(void)` """
        return self.thisptr.delta()
    cpdef double speed_sound(self) except *:
        """ Get the speed of sound in m/s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::speed_sound(void)` """
        return self.thisptr.speed_sound()
    cpdef double molar_mass(self) except *:
        """ Get the molar mass in kg/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::molar_mass(void)` """
        return self.thisptr.molar_mass()
    cpdef double acentric_factor(self) except *:
        """ Get the acentric factor - wrapper of c++ function :cpapi:`CoolProp::AbstractState::acentric_factor(void)` """
        return self.thisptr.acentric_factor()
    cpdef double gas_constant(self) except *:
        """ Get the gas constant in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gas_constant(void)` """
        return self.thisptr.gas_constant()
    cpdef double viscosity(self) except *:
        """ Get the viscosity in Pa-s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::viscosity(void)` """
        return self.thisptr.viscosity()
    cpdef double conductivity(self) except *:
        """ Get the thermal conductivity in W/m/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conductivity(void)` """
        return self.thisptr.conductivity()
    cpdef double surface_tension(self) except *:
        """ Get the surface tension N/m - wrapper of c++ function :cpapi:`CoolProp::AbstractState::surface_tension(void)` """
        return self.thisptr.surface_tension()
    cpdef double Prandtl(self) except *:
        """ Get the Prandtl number - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Prandtl(void)` """
        return self.thisptr.Prandtl()
    cpdef double Bvirial(self) except *:
        """ Get the B virial coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Bvirial(void)` """
        return self.thisptr.Bvirial()
    cpdef double Cvirial(self) except *:
        """ Get the C virial coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Cvirial(void)` """
        return self.thisptr.Cvirial()
    cpdef double fundamental_derivative_of_gas_dynamics(self) except *:
        """ Get the fundamental derivative of gas dynamics - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fundamental_derivative_of_gas_dynamics(void)` """
        return self.thisptr.fundamental_derivative_of_gas_dynamics()
    cpdef double PIP(self) except *:
        """ Get the phase identification parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::PIP` """
        return self.thisptr.PIP()
    cpdef double isobaric_expansion_coefficient(self) except *:
        """ Get the isobaric expansion coefficient - wrapper of c++ function :cpapi:`CoolProp::AbstractState::isobaric_expansion_coefficient(void)` """
        return self.thisptr.isobaric_expansion_coefficient()
    cpdef double isothermal_compressibility(self) except *:
        """ Get the isothermal_compressibility - wrapper of c++ function :cpapi:`CoolProp::AbstractState::isothermal_compressibility(void)` """
        return self.thisptr.isothermal_compressibility()
    cpdef double fugacity(self, size_t i) except *:
        """ Get the fugacity of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fugacity(std::size_t)` """
        return self.thisptr.fugacity(i)
    cpdef double fugacity_coefficient(self, size_t i) except *:
        """ Get the fugacity coefficient of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::fugacity_coefficient(std::size_t)` """
        return self.thisptr.fugacity_coefficient(i)
    cpdef double chemical_potential(self, size_t i) except *:
        """ Get the chemical potential of the i-th component - wrapper of c++ function :cpapi:`CoolProp::AbstractState::chemical_potential(std::size_t)` """
        return self.thisptr.chemical_potential(i)

    cpdef mole_fractions_liquid(self):
        """ Get the mole fractions of the liquid phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_liquid(void)` """
        return self.thisptr.mole_fractions_liquid()
    cpdef mole_fractions_vapor(self):
        """ Get the mole fractions of the vapor phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_vapor(void)` """
        return self.thisptr.mole_fractions_vapor()

    cpdef tuple true_critical_point(self):
        """ Get the "true" critical point where dp/drho|T = 0 & d2p/drho^2|T = 0 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::true_critical_point` """
        cdef double T = 1e99, rho = 1e99
        self.thisptr.true_critical_point(T, rho)
        return T, rho
    cpdef dict conformal_state(self, string reference_fluid, CoolPropDbl T, CoolPropDbl rho):
        """ Solve for conformal state used in extended corresponding states - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conformal_state` """
        cdef CoolPropDbl T0 = T, rho0 = rho
        self.thisptr.conformal_state(reference_fluid, T0, rho0)
        return dict(T = T0, rhomolar = rho0)
    cpdef dict conductivity_contributions(self):
        """ Retrieve each of the contributions to the conductivity, each in W/m/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::conductivity_contributions` """
        cdef CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0
        self.thisptr.conductivity_contributions(dilute, initial_density, residual, critical)
        return dict(dilute = dilute, initial_density = initial_density, residual = residual, critical = critical)
    cpdef dict viscosity_contributions(self):
        """ Retrieve each of the contributions to the viscosity, each in Pa-s - wrapper of c++ function :cpapi:`CoolProp::AbstractState::viscosity_contributions` """
        cdef CoolPropDbl dilute = 0, initial_density = 0, residual = 0, critical = 0
        self.thisptr.viscosity_contributions(dilute, initial_density, residual, critical)
        return dict(dilute = dilute, initial_density = initial_density, residual = residual, critical = critical)


    cpdef double helmholtzmolar_excess(self) except *:
        """ Get the mole-specific excess Helmholtz energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmolar_excess(void)` """
        return self.thisptr.helmholtzmolar_excess()
    cpdef double helmholtzmass_excess(self) except *:
        """ Get the mass-specific excess Helmholtz energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::helmholtzmass_excess(void)` """
        return self.thisptr.helmholtzmass_excess()
    cpdef double gibbsmolar_excess(self) except *:
        """ Get the mole-specific excess Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar_excess(void)` """
        return self.thisptr.gibbsmolar_excess()
    cpdef double gibbsmass_excess(self) except *:
        """ Get the mass-specific excess Gibbs energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmass_excess(void)` """
        return self.thisptr.gibbsmass_excess()
    cpdef double umolar_excess(self) except *:
        """ Get the mole-specific excess internal energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umolar_excess(void)` """
        return self.thisptr.umolar_excess()
    cpdef double umass_excess(self) except *:
        """ Get the mass-specific excess internal energy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::umass_excess(void)` """
        return self.thisptr.umass_excess()
    cpdef double hmolar_excess(self) except *:
        """ Get the mole-specific excess enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar_excess(void)` """
        return self.thisptr.hmolar_excess()
    cpdef double hmass_excess(self) except *:
        """ Get the mass-specific excess enthalpy in J/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmass_excess(void)` """
        return self.thisptr.hmass_excess()
    cpdef double smolar_excess(self) except *:
        """ Get the mole-specific excess entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar_excess(void)` """
        return self.thisptr.smolar_excess()
    cpdef double smass_excess(self) except *:
        """ Get the mass-specific excess entropy in J/kg/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smass_excess(void)` """
        return self.thisptr.smass_excess()
    cpdef double volumemolar_excess(self) except *:
        """ Get the mole-specific excess volume in m^3/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::volumemolar_excess(void)` """
        return self.thisptr.volumemolar_excess()
    cpdef double volumemass_excess(self) except *:
        """ Get the mass-specific excess volume in m^3/kg - wrapper of c++ function :cpapi:`CoolProp::AbstractState::volumemass_excess(void)` """
        return self.thisptr.volumemass_excess()

    cpdef double gibbsmolar_residual(self) except *:
        """ Get the mole-specific residual Gibbs energy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::gibbsmolar_residual(void)` """
        return self.thisptr.gibbsmolar_residual()
    cpdef double hmolar_residual(self) except *:
        """ Get the mole-specific residual enthalpy in J/mol - wrapper of c++ function :cpapi:`CoolProp::AbstractState::hmolar_residual(void)` """
        return self.thisptr.hmolar_residual()
    cpdef double smolar_residual(self) except *:
        """ Get the mole-specific residual entropy in J/mol/K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::smolar_residual(void)` """
        return self.thisptr.smolar_residual()


    ## ----------------------------------------
    ##        Derivatives
    ## ----------------------------------------

    cpdef CoolPropDbl first_partial_deriv(self, constants_header.parameters OF , constants_header.parameters WRT, constants_header.parameters CONSTANT) except *:
        """ Get the first partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_partial_deriv` """
        return self.thisptr.first_partial_deriv(OF, WRT, CONSTANT)
    cpdef CoolPropDbl second_partial_deriv(self, constants_header.parameters OF , constants_header.parameters WRT1, constants_header.parameters CONSTANT1, constants_header.parameters WRT2, constants_header.parameters CONSTANT2) except *:
        """ Get the second partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_partial_deriv` """
        return self.thisptr.second_partial_deriv(OF, WRT1, CONSTANT1, WRT2, CONSTANT2)
    cpdef CoolPropDbl first_saturation_deriv(self, constants_header.parameters OF , constants_header.parameters WRT) except *:
        """ Get the first derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_saturation_deriv` """
        return self.thisptr.first_saturation_deriv(OF, WRT)
    cpdef CoolPropDbl second_saturation_deriv(self, constants_header.parameters OF1 , constants_header.parameters WRT1, constants_header.parameters WRT2) except *:
        """ Get the second derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_saturation_deriv` """
        return self.thisptr.second_saturation_deriv(OF1, WRT1, WRT2)
    cpdef double first_two_phase_deriv(self, constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant) except *:
        """ Get the first two-phase derivative - wrapper of C++ function :cpapi:`CoolProp::AbstractState::first_two_phase_deriv` """
        return self.thisptr.first_two_phase_deriv(Of, Wrt, Constant)
    cpdef double second_two_phase_deriv(self, constants_header.parameters Of1, constants_header.parameters Wrt1, constants_header.parameters Constant1, constants_header.parameters Wrt2, constants_header.parameters Constant2) except *:
        """ Get the second two-phase derivative - wrapper of C++ function :cpapi:`CoolProp::AbstractState::second_two_phase_deriv` """
        return self.thisptr.second_two_phase_deriv(Of1, Wrt1, Constant1, Wrt2, Constant2)
    cpdef double first_two_phase_deriv_splined(self, constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant, double x_end) except *:
        """ Get the first two-phase derivative using splines - wrapper of C++ function :cpapi:`CoolProp::AbstractState::first_two_phase_deriv_splined` """
        return self.thisptr.first_two_phase_deriv_splined(Of, Wrt, Constant, x_end)

    ## ----------------------------------------
    ##        Ancillary curves
    ## ----------------------------------------

    cpdef bint has_melting_line(self) except *:
        """ Check if the fluid has a melting line - True if is does, False otherwise - wrapper of c++ function :cpapi:`CoolProp::AbstractState::has_melting_line` """
        return self.thisptr.has_melting_line()
    cpdef double melting_line(self, int param, int given, double value) except *:
        """ Get values from the melting line - wrapper of c++ function :cpapi:`CoolProp::AbstractState::melting_line` """
        return self.thisptr.melting_line(param, given, value)
    cpdef double saturation_ancillary(self, constants_header.parameters param, int Q, constants_header.parameters given, double value) except *:
        """ Get values from the saturation_ancillary - wrapper of c++ function :cpapi:`CoolProp::AbstractState::saturation_ancillary` """
        return self.thisptr.saturation_ancillary(param, Q, given, value)

    ## ----------------------------------------
    ##        Phase envelope
    ## ----------------------------------------

    cpdef build_phase_envelope(self, string type):
        """ Build the phase envelope - wrapper of c++ function :cpapi:`CoolProp::AbstractState::build_phase_envelope` """
        self.thisptr.build_phase_envelope(type)
    cpdef PyPhaseEnvelopeData get_phase_envelope_data(self):
        """ Get the phase envelope data - wrapper of c++ function :cpapi:`CoolProp::AbstractState::get_phase_envelope_data` """
        cdef cAbstractState.PhaseEnvelopeData pe_data = self.thisptr.get_phase_envelope_data()
        cdef PyPhaseEnvelopeData pe_out = PyPhaseEnvelopeData()
        pe_out.T = pe_data.T
        pe_out.p = pe_data.p
        pe_out.Q = pe_data.Q
        pe_out.rhomolar_liq = pe_data.rhomolar_liq
        pe_out.rhomolar_vap = pe_data.rhomolar_vap
        pe_out.hmolar_liq = pe_data.hmolar_liq
        pe_out.hmolar_vap = pe_data.hmolar_vap
        pe_out.smolar_liq = pe_data.smolar_liq
        pe_out.smolar_vap = pe_data.smolar_vap
        pe_out.iTsat_max = pe_data.iTsat_max
        pe_out.ipsat_max = pe_data.ipsat_max
        pe_out.TypeI = pe_data.TypeI
        pe_out.x = pe_data.x
        pe_out.y = pe_data.y
        pe_out.K = pe_data.K
        return pe_out

    ## -----------------------------------------
    ##   Ideal curves
    ## -----------------------------------------

    cpdef tuple ideal_curve(self, string type):
        """ Get an ideal curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::ideal_curve` """
        cdef vector[double] T, p
        self.thisptr.ideal_curve(type, T, p)
        return T, p

    ## -----------------------------------------
    ##   Helmholtz energy derivatives
    ## -----------------------------------------

    cpdef CoolPropDbl alpha0(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alpha0` """
        return self.thisptr.alpha0()
    cpdef CoolPropDbl dalpha0_dDelta(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dDelta` """
        return self.thisptr.dalpha0_dDelta()
    cpdef CoolPropDbl dalpha0_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dTau` """
        return self.thisptr.dalpha0_dTau()
    cpdef CoolPropDbl d2alpha0_dDelta2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta2` """
        return self.thisptr.d2alpha0_dDelta2()
    cpdef CoolPropDbl d2alpha0_dDelta_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta_dTau` """
        return self.thisptr.d2alpha0_dDelta_dTau()
    cpdef CoolPropDbl d2alpha0_dTau2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dTau2` """
        return self.thisptr.d2alpha0_dTau2()
    cpdef CoolPropDbl d3alpha0_dTau3(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dTau3` """
        return self.thisptr.d3alpha0_dTau3()
    cpdef CoolPropDbl d3alpha0_dDelta_dTau2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta_dTau2` """
        return self.thisptr.d3alpha0_dDelta_dTau2()
    cpdef CoolPropDbl d3alpha0_dDelta2_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta2_dTau` """
        return self.thisptr.d3alpha0_dDelta2_dTau()
    cpdef CoolPropDbl d3alpha0_dDelta3(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta3` """
        return self.thisptr.d3alpha0_dDelta3()

    cpdef CoolPropDbl alphar(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alphar` """
        return self.thisptr.alphar()
    cpdef CoolPropDbl dalphar_dDelta(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dDelta` """
        return self.thisptr.dalphar_dDelta()
    cpdef CoolPropDbl dalphar_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dTau` """
        return self.thisptr.dalphar_dTau()
    cpdef CoolPropDbl d2alphar_dDelta2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta2` """
        return self.thisptr.d2alphar_dDelta2()
    cpdef CoolPropDbl d2alphar_dDelta_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta_dTau` """
        return self.thisptr.d2alphar_dDelta_dTau()
    cpdef CoolPropDbl d2alphar_dTau2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dTau2` """
        return self.thisptr.d2alphar_dTau2()
    cpdef CoolPropDbl d3alphar_dTau3(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dTau3` """
        return self.thisptr.d3alphar_dTau3()
    cpdef CoolPropDbl d3alphar_dDelta_dTau2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta_dTau2` """
        return self.thisptr.d3alphar_dDelta_dTau2()
    cpdef CoolPropDbl d3alphar_dDelta2_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta2_dTau` """
        return self.thisptr.d3alphar_dDelta2_dTau()
    cpdef CoolPropDbl d3alphar_dDelta3(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta3` """
        return self.thisptr.d3alphar_dDelta3()
    cpdef CoolPropDbl d4alphar_dTau4(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dTau4` """
        return self.thisptr.d4alphar_dTau4()
    cpdef CoolPropDbl d4alphar_dDelta_dTau3(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta_dTau3` """
        return self.thisptr.d4alphar_dDelta_dTau3()
    cpdef CoolPropDbl d4alphar_dDelta2_dTau2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta2_dTau2` """
        return self.thisptr.d4alphar_dDelta2_dTau2()
    cpdef CoolPropDbl d4alphar_dDelta3_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta3_dTau` """
        return self.thisptr.d4alphar_dDelta3_dTau()
    cpdef CoolPropDbl d4alphar_dDelta4(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d4alphar_dDelta4` """
        return self.thisptr.d4alphar_dDelta4()
