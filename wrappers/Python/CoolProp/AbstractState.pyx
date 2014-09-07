# This file is embedded directly in CoolProp.pyx

cimport constants_header
        
cdef class PyPhaseEnvelopeData:
    pass
    
cdef class AbstractState:
    """
    This class is a one-to-one python wrapper of the :cpapi:`AbstractState` class
    """
    
    def __cinit__(self, string backend, string fluid):
        self.thisptr = cAbstractState.factory(backend, fluid)
        
    def __dealloc__(self):
        del self.thisptr
    
    cpdef update(self, constants_header.input_pairs ipair, double Value1, double Value2):
        """ Update function - mirrors c++ function :cpapi:`AbstractState::update` """
        self.thisptr.update(ipair, Value1, Value2)
    
    cpdef set_mole_fractions(self, vector[double] z): 
        """ Set the mole fractions - wrapper of c++ function :cpapi:`AbstractState::set_mole_fractions` """
        self.thisptr.set_mole_fractions(z)
        
    ## ----------------------------------------	
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double keyed_output(self, parameters iOutput) except *: 
        """ Update :cpapi:`AbstractState::update` """
        return self.thisptr.keyed_output(iOutput)
    
    cpdef double T(self) except *: 
        """ Get the temperature in K - wrapper of c++ function :cpapi:`AbstractState::T` """
        return self.thisptr.T()
    cpdef double p(self) except *: 
        """ Get the pressure in Pa - wrapper of c++ function :cpapi:`AbstractState::p` """
        return self.thisptr.p()
    cpdef double rhomolar(self) except *: 
        """ Get the density in mol/m^3 - wrapper of c++ function :cpapi:`AbstractState::rhomolar` """
        return self.thisptr.rhomolar()
    cpdef double rhomass(self) except *: 
        """ Get the density in kg/m^3 - wrapper of c++ function :cpapi:`AbstractState::rhomass` """
        return self.thisptr.rhomass()
    cpdef double hmolar(self) except *: 
        """ Get the enthalpy in J/mol - wrapper of c++ function :cpapi:`AbstractState::hmolar` """
        return self.thisptr.hmolar()
    cpdef double smolar(self) except *: 
        """ Get the entropy in J/mol/K - wrapper of c++ function :cpapi:`AbstractState::smolar` """
        return self.thisptr.smolar()
    cpdef double cpmolar(self) except *: 
        """ Get the constant pressure specific heat in J/mol/K - wrapper of c++ function :cpapi:`AbstractState::cpmolar` """
        return self.thisptr.cpmolar()
    cpdef double cvmolar(self) except *: 
        """ Get the constant volume specific heat in J/mol/K - wrapper of c++ function :cpapi:`AbstractState::cvmolar` """
        return self.thisptr.cvmolar()
    cpdef double hmass(self) except *: 
        """ Get the enthalpy in J/kg - wrapper of c++ function :cpapi:`AbstractState::hmass` """
        return self.thisptr.hmass()
    cpdef double smass(self) except *: 
        """ Get the entropy in J/kg/K - wrapper of c++ function :cpapi:`AbstractState::smass` """
        return self.thisptr.smass()
    cpdef double cpmass(self) except *: 
        """ Get the constant pressure specific heat in J/kg/K - wrapper of c++ function :cpapi:`AbstractState::cpmass` """
        return self.thisptr.cpmass()
    cpdef double cvmass(self) except *: 
        """ Get the constant volume specific heat in J/kg/K - wrapper of c++ function :cpapi:`AbstractState::cvmass` """
        return self.thisptr.cvmass()
    cpdef double speed_sound(self) except *: 
        """ Get the speed of sound in m/s - wrapper of c++ function :cpapi:`AbstractState::speed_sound` """
        return self.thisptr.speed_sound()
    cpdef double molar_mass(self) except *: 
        """ Get the molar mass in kg/mol - wrapper of c++ function :cpapi:`AbstractState::molar_mass` """
        return self.thisptr.molar_mass()
        
    ## ----------------------------------------	
    ##        Melting Line
    ## ----------------------------------------
    
    cpdef double melting_line(self, int param, int given, double value) except *: 
        """ Get values from the melting line - wrapper of c++ function :cpapi:`AbstractState::melting_line` """
        return self.thisptr.melting_line(param, given, value)
    cpdef bint has_melting_line(self) except *: 
        """ Check if the fluid has a melting line - True if is does, False otherwise - wrapper of c++ function :cpapi:`AbstractState::has_melting_line` """
        return self.thisptr.has_melting_line()
    
    ## ----------------------------------------	
    ##        Phase envelope
    ## ----------------------------------------
    
    cpdef build_phase_envelope(self, string type):
        """ Build the phase envelope - wrapper of c++ function :cpapi:`AbstractState::build_phase_envelope` """
        self.thisptr.build_phase_envelope(type)
    cpdef PyPhaseEnvelopeData get_phase_envelope_data(self):
        """ Get the phase envelope data - wrapper of c++ function :cpapi:`AbstractState::get_phase_envelope_data` """
        cdef cAbstractState.PhaseEnvelopeData pe_data = self.thisptr.get_phase_envelope_data()
        cdef PyPhaseEnvelopeData pe_out = PyPhaseEnvelopeData()
        pe_out.T = pe_data.T
        pe_out.p = pe_data.p
        pe_out.rhomolar_liq = pe_data.rhomolar_liq
        pe_out.rhomolar_vap = pe_data.rhomolar_vap
        pe_out.hmolar_liq = pe_data.hmolar_liq
        pe_out.hmolar_vap = pe_data.hmolar_vap
        
        return pe_out
        