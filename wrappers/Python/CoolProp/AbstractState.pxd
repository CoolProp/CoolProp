from libcpp cimport bool 
from libcpp.string cimport string

# A header defining the AbstractState class
cimport cAbstractState

cimport constants_header

cdef class PyPhaseEnvelopeData:
    cpdef public bool TypeI
    cpdef public size_t iTsat_max, ipsat_max, icrit
    cpdef public list T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q
    
cdef class PyGuessesStructure:
    cpdef public long double T, p, rhomolar, hmolar, smolar
    cpdef public long double rhomolar_liq, rhomolar_vap
    cpdef public vector[long double] x, y

cdef class AbstractState:
    cdef cAbstractState.AbstractState *thisptr     # hold a C++ instance which we're wrapping
    cpdef update(self, constants_header.input_pairs iInput1, double Value1, double Value2)
    cpdef update_with_guesses(self, constants_header.input_pairs iInput1, double Value1, double Value2, PyGuessesStructure guesses)
    cpdef set_mole_fractions(self, vector[double] z)
    cpdef set_mass_fractions(self, vector[double] z)
    cpdef set_volu_fractions(self, vector[double] z)
    
    cpdef name(self)
    
    cpdef constants_header.phases phase(self) except *  
    cpdef specify_phase(self, constants_header.phases phase)
    cpdef unspecify_phase(self)
    
    ## Limits
    cpdef double Tmin(self)
    cpdef double Tmax(self)
    cpdef double pmax(self)
    cpdef double Ttriple(self)
        
    ## ---------------------------------------- 
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self) except *
    cpdef double p(self) except *
    cpdef double Q(self) except *
    cpdef double rhomolar(self) except *
    cpdef double hmolar(self) except *
    cpdef double smolar(self) except *
    cpdef double umolar(self) except *
    cpdef double cpmolar(self) except *
    cpdef double cp0molar(self) except *
    cpdef double cvmolar(self) except *
    cpdef double rhomass(self) except *
    cpdef double hmass(self) except *
    cpdef double smass(self) except *
    cpdef double umass(self) except *
    cpdef double cpmass(self) except *
    cpdef double cp0mass(self) except *
    cpdef double cvmass(self) except *
    cpdef double speed_sound(self) except *
    cpdef double gas_constant(self) except *
    cpdef double tau(self) except *
    cpdef double delta(self) except *
    cpdef double viscosity(self) except *
    cpdef double conductivity(self) except *
    cpdef double surface_tension(self) except *
    cpdef double Prandtl(self) except *
    
    cpdef double molar_mass(self) except *
    cpdef double acentric_factor(self) except*
    cpdef double keyed_output(self, constants_header.parameters) except *
    cpdef double trivial_keyed_output(self, constants_header.parameters) except *
    cpdef double saturated_liquid_keyed_output(self, constants_header.parameters) except *
    cpdef double saturated_vapor_keyed_output(self, constants_header.parameters) except *
    
    ## ----------------------------------------	
    ##        Derivatives
    ## ----------------------------------------
    
    cpdef long double first_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef long double second_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef long double first_saturation_deriv(self, constants_header.parameters, constants_header.parameters) except *
    cpdef long double second_saturation_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    
    cpdef double first_two_phase_deriv(self, constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant) except *
    cpdef double second_two_phase_deriv(self, constants_header.parameters Of, constants_header.parameters Wrt1, constants_header.parameters Constant1, constants_header.parameters Wrt2, constants_header.parameters Constant2) except *
    cpdef double first_two_phase_deriv_splined(self ,constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant, double x_end) except *
    
    cpdef double melting_line(self, int, int, double) except *
    cpdef bool has_melting_line(self) except *
    cpdef double saturation_ancillary(self, constants_header.parameters, int, constants_header.parameters, double) except *
        
    cpdef build_phase_envelope(self, string)
    cpdef PyPhaseEnvelopeData get_phase_envelope_data(self)
    
    cpdef mole_fractions_liquid(self)
    cpdef mole_fractions_vapor(self)
    
    cpdef long double alpha0(self) except *
    cpdef long double dalpha0_dDelta(self) except *
    cpdef long double dalpha0_dTau(v) except *
    cpdef long double d2alpha0_dDelta2(self) except *
    cpdef long double d2alpha0_dDelta_dTau(self) except *
    cpdef long double d2alpha0_dTau2(self) except *
    cpdef long double d3alpha0_dTau3(self) except *
    cpdef long double d3alpha0_dDelta_dTau2(self) except *
    cpdef long double d3alpha0_dDelta2_dTau(self) except *
    cpdef long double d3alpha0_dDelta3(self) except *

    cpdef long double alphar(self) except *
    cpdef long double dalphar_dDelta(self) except *
    cpdef long double dalphar_dTau(self) except *
    cpdef long double d2alphar_dDelta2(self) except *
    cpdef long double d2alphar_dDelta_dTau(self) except *
    cpdef long double d2alphar_dTau2(self) except *
    cpdef long double d3alphar_dDelta3(self) except *
    cpdef long double d3alphar_dDelta2_dTau(self) except *
    cpdef long double d3alphar_dDelta_dTau2(self) except *
    cpdef long double d3alphar_dTau3(self) except *