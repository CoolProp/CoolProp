# This file is embedded directly in CoolProp.pyx

cimport constants_header
        
cdef class PyPhaseEnvelopeData:
    pass

cdef class PyGuessesStructure:
    pass
    
cdef class AbstractState:
    """
    This class is a one-to-one python wrapper of the :cpapi:`AbstractState` class
    """
    
    def __cinit__(self, string backend, string fluid):
        self.thisptr = cAbstractState.factory(backend, fluid)
        
    def __dealloc__(self):
        del self.thisptr
        
    cpdef name(self):
        """ Get the backend name - wrapper of c++ function :cpapi:`CoolProp::AbstractState::name` """
        return self.thisptr.name()
        
    cpdef constants_header.phases phase(self) except *:
        """ Get the phase as key value- wrapper of c++ function :cpapi:`CoolProp::AbstractState::phase` """
        return self.thisptr.phase()
        
    cpdef specify_phase(self, constants_header.phases phase):
        """ Specify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::specify_phase` """
        self.thisptr.specify_phase(phase)
    cpdef unspecify_phase(self):
        """ Unspecify the phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::unspecify_phase` """
        self.thisptr.unspecify_phase()
    
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
    
    ## ----------------------------------------	
    ##        Limits
    ## ----------------------------------------
    cpdef double Tmin(self):
        """ Set the minimum temperature in K- wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmin` """
        return self.thisptr.Tmin()
    cpdef double Tmax(self):
        """ Set the maximum temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Tmax` """
        return self.thisptr.Tmax()
    cpdef double pmax(self):
        """ Set the maximum pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::pmax` """
        return self.thisptr.pmax()
    cpdef double Ttriple(self):
        """ Set the triple point temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::Ttriple` """
        return self.thisptr.Ttriple()

    ## ----------------------------------------	
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double keyed_output(self, parameters iOutput) except *: 
        """ Update :cpapi:`CoolProp::AbstractState::keyed_output(parameters key)` """
        return self.thisptr.keyed_output(iOutput)
    cpdef double trivial_keyed_output(self, parameters iOutput) except *: 
        """ Update :cpapi:`CoolProp::AbstractState::trivial_keyed_output(parameters key)` """
        return self.thisptr.trivial_keyed_output(iOutput)
    cpdef double saturated_liquid_keyed_output(self, parameters iOutput) except *: 
        """ Update :cpapi:`CoolProp::AbstractState::saturated_liquid_keyed_output(parameters key)` """
        return self.thisptr.saturated_liquid_keyed_output(iOutput)
    cpdef double saturated_vapor_keyed_output(self, parameters iOutput) except *: 
        """ Update :cpapi:`CoolProp::AbstractState::saturated_vapor_keyed_output(parameters key)` """
        return self.thisptr.saturated_vapor_keyed_output(iOutput)
    
    cpdef double T(self) except *: 
        """ Get the temperature in K - wrapper of c++ function :cpapi:`CoolProp::AbstractState::T(void)` """
        return self.thisptr.T()
    cpdef double p(self) except *: 
        """ Get the pressure in Pa - wrapper of c++ function :cpapi:`CoolProp::AbstractState::p(void)` """
        return self.thisptr.p()
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
        return self.thisptr.molar_mass()
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
    cpdef double PIP(self) except *: 
        """ Get the phase identification parameter - wrapper of c++ function :cpapi:`CoolProp::AbstractState::PIP(void)` """
        return self.thisptr.PIP()
    
    cpdef mole_fractions_liquid(self):
        """ Get the mole fractions of the liquid phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_liquid(void)` """
        return self.thisptr.mole_fractions_liquid()
    cpdef mole_fractions_vapor(self):
        """ Get the mole fractions of the vapor phase - wrapper of c++ function :cpapi:`CoolProp::AbstractState::mole_fractions_vapor(void)` """
        return self.thisptr.mole_fractions_vapor()
    
    cpdef tuple true_critical_point(self):
        """ Get the "true" critical point where dp/drho|T = 0 & d2p/drho^2|T = 0 - wrapper of c++ function :cpapi:`CoolProp::AbstractState::true_critical_point(void)` """
        cdef double T = 1e99, rho = 1e99
        self.thisptr.true_critical_point(T, rho)
        return T, rho
        
    ## ----------------------------------------	
    ##        Derivatives
    ## ----------------------------------------
    
    cpdef long double first_partial_deriv(self, constants_header.parameters OF , constants_header.parameters WRT, constants_header.parameters CONSTANT) except *: 
        """ Get the first partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_partial_deriv` """
        return self.thisptr.first_partial_deriv(OF, WRT, CONSTANT)
    cpdef long double second_partial_deriv(self, constants_header.parameters OF , constants_header.parameters WRT1, constants_header.parameters CONSTANT1, constants_header.parameters WRT2, constants_header.parameters CONSTANT2) except *: 
        """ Get the second partial derivative - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_partial_deriv` """
        return self.thisptr.second_partial_deriv(OF, WRT1, CONSTANT1, WRT2, CONSTANT2)
    cpdef long double first_saturation_deriv(self, constants_header.parameters OF , constants_header.parameters WRT) except *: 
        """ Get the first derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::first_saturation_deriv` """
        return self.thisptr.first_saturation_deriv(OF, WRT)
    cpdef long double second_saturation_deriv(self, constants_header.parameters OF1 , constants_header.parameters WRT1, constants_header.parameters OF2, constants_header.parameters WRT2) except *: 
        """ Get the second derivative along the saturation curve - wrapper of c++ function :cpapi:`CoolProp::AbstractState::second_saturation_deriv` """
        return self.thisptr.second_saturation_deriv(OF1, WRT1, OF2, WRT2)
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
        return pe_out
        
    ## -----------------------------------------
    ##   Helmholtz energy derivatives
    ## -----------------------------------------
    
    cpdef long double alpha0(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alpha0` """
        return self.thisptr.alpha0()
    cpdef long double dalpha0_dDelta(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dDelta` """
        return self.thisptr.dalpha0_dDelta()
    cpdef long double dalpha0_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalpha0_dTau` """
        return self.thisptr.dalpha0_dTau()
    cpdef long double d2alpha0_dDelta2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta2` """
        return self.thisptr.d2alpha0_dDelta2()
    cpdef long double d2alpha0_dDelta_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dDelta_dTau` """
        return self.thisptr.d2alpha0_dDelta_dTau()
    cpdef long double d2alpha0_dTau2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alpha0_dTau2` """
        return self.thisptr.d2alpha0_dTau2()
    cpdef long double d3alpha0_dTau3(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dTau3` """
        return self.thisptr.d3alpha0_dTau3()
    cpdef long double d3alpha0_dDelta_dTau2(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta_dTau2` """
        return self.thisptr.d3alpha0_dDelta_dTau2()
    cpdef long double d3alpha0_dDelta2_dTau(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta2_dTau` """
        return self.thisptr.d3alpha0_dDelta2_dTau()
    cpdef long double d3alpha0_dDelta3(self) except *:
        """ Get the ideal-gas reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alpha0_dDelta3` """
        return self.thisptr.d3alpha0_dDelta3()
        
    cpdef long double alphar(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::alphar` """
        return self.thisptr.alphar()
    cpdef long double dalphar_dDelta(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dDelta` """
        return self.thisptr.dalphar_dDelta()
    cpdef long double dalphar_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::dalphar_dTau` """
        return self.thisptr.dalphar_dTau()
    cpdef long double d2alphar_dDelta2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta2` """
        return self.thisptr.d2alphar_dDelta2()
    cpdef long double d2alphar_dDelta_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dDelta_dTau` """
        return self.thisptr.d2alphar_dDelta_dTau()
    cpdef long double d2alphar_dTau2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d2alphar_dTau2` """
        return self.thisptr.d2alphar_dTau2()
    cpdef long double d3alphar_dTau3(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dTau3` """
        return self.thisptr.d3alphar_dTau3()
    cpdef long double d3alphar_dDelta_dTau2(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta_dTau2` """
        return self.thisptr.d3alphar_dDelta_dTau2()
    cpdef long double d3alphar_dDelta2_dTau(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta2_dTau` """
        return self.thisptr.d3alphar_dDelta2_dTau()
    cpdef long double d3alphar_dDelta3(self) except *:
        """ Get the residual reduced Helmholtz energy - wrapper of c++ function :cpapi:`CoolProp::AbstractState::d3alphar_dDelta3` """
        return self.thisptr.d3alphar_dDelta3()