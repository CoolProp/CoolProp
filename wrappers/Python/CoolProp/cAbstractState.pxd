from libcpp cimport bool 
from libcpp.string cimport string

cdef extern from "AbstractState.h" namespace "CoolProp":
    cdef cppclass AbstractState:
        
        ## Nullary Constructor
        AbstractState() except +ValueError
        
        ## Constructor with fluid name
        AbstractState(string FluidName) except +ValueError

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(long iInput1, double Value1, double Value2) except +ValueError

        ## Bulk properties accessors - temperature and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T() except +ValueError
        double rhomolar() except +ValueError
        double p() except +ValueError
        double hmolar() except +ValueError
        double smolar() except +ValueError
        double cpmolar() except +ValueError
        double cvmolar() except +ValueError
        double speed_sound() except +ValueError
        double rhomass() except +ValueError
        double hmass() except +ValueError
        double smass() except +ValueError
        double cpmass() except +ValueError
        double cvmass() except +ValueError
        
        double keyed_output(long) except+ValueError
        double molar_mass() except+ValueError
        double gas_constant() except+ValueError
        double build_phase_envelope() except+ValueError
        double viscosity() except+ValueError
        double conductivity() except+ValueError
        double surface_tension() except+ValueError
        
        double melting_line(int,int,double) except+ValueError
        bool has_melting_line() except+ValueError

# The static factory method for the AbstractState
cdef extern from "AbstractState.h" namespace "CoolProp::AbstractState":
    AbstractState* factory(const string &backend, const string &fluid_string) except+ValueError
