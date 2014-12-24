from libcpp cimport bool 
from libcpp.string cimport string

# A header defining the AbstractState class
cimport cAbstractState

cimport constants_header

cdef class PyPhaseEnvelopeData:
    cpdef public bool TypeI
    cpdef public size_t iTsat_max, ipsat_max, icrit
    cpdef public list T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q

cdef class AbstractState:
    cdef cAbstractState.AbstractState *thisptr     # hold a C++ instance which we're wrapping
    cpdef update(self, constants_header.input_pairs iInput1, double Value1, double Value2)
    cpdef set_mole_fractions(self, vector[double] z)
    cpdef set_mass_fractions(self, vector[double] z)
    cpdef set_volu_fractions(self, vector[double] z)
    
    cpdef specify_phase(self, constants_header.phases phase)
    cpdef unspecify_phase(self)
    
    ## ---------------------------------------- 
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self)
    cpdef double p(self)
    cpdef double rhomolar(self)
    cpdef double hmolar(self)
    cpdef double smolar(self)
    cpdef double umolar(self)
    cpdef double cpmolar(self)
    cpdef double cp0molar(self)
    cpdef double cvmolar(self)
    cpdef double rhomass(self)
    cpdef double hmass(self) except *
    cpdef double smass(self) except *
    cpdef double umass(self) except *
    cpdef double cpmass(self) except *
    cpdef double cp0mass(self) except *
    cpdef double cvmass(self) except *
    cpdef double speed_sound(self) except *
    cpdef double gas_constant(self) except *
    
    cpdef double molar_mass(self) except *
    cpdef double keyed_output(self, constants_header.parameters) except *
    
    cpdef long double first_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef long double second_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef long double first_saturation_deriv(self, constants_header.parameters, constants_header.parameters) except *
    cpdef long double second_saturation_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    
    cpdef double melting_line(self, int, int, double) except *
    cpdef bool has_melting_line(self) except *
    
    cpdef build_phase_envelope(self, string)
    cpdef PyPhaseEnvelopeData get_phase_envelope_data(self)
    
    cpdef mole_fractions_liquid(self)
    cpdef mole_fractions_vapor(self)