#cython: embedsignature = True, c_string_type=str, c_string_encoding=ascii
        
cdef class AbstractState:
    """
    This class is a one-to-one python wrapper of the :cpapi:`AbstractState` class
    """
    
    def __cinit__(self, string backend, string fluid):
        self.thisptr = cAbstractState.factory(backend, fluid)
        
    def __dealloc__(self):
        del self.thisptr
    
    cpdef update(self, long ipair, double Value1, double Value2):
        """ Update function - mirrors c++ function :cpapi:`AbstractState::update` """
        self.thisptr.update(ipair, Value1, Value2)
    
    ## ----------------------------------------	
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double keyed_output(self, long iOutput) except *: 
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
        """ Get the molar mass kg/mol - wrapper of c++ function :cpapi:`AbstractState::molar_mass` """
        return self.thisptr.molar_mass()   
    
    cpdef double melting_line(self, int param, int given, double value) except *: 
        """ Get values from the melting line - wrapper of c++ function :cpapi:`AbstractState::melting_line` """
        return self.thisptr.melting_line(param, given, value)