from libcpp cimport bool 
from libcpp.string cimport string

cdef extern from "AbstractState.h" namespace "CoolProp":
    cdef cppclass AbstractState:
        
        ## Nullary Constructor
        AbstractState() except +
        
        ## Constructor with fluid name
        AbstractState(string FluidName) except +

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(long iInput1, double Value1, double Value2) except +

        ## Bulk properties accessors - temperature and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T() except +
        double rhomolar() except +
        double p() except +
        double hmolar() except +
        double smolar() except +
        double cpmolar() except +
        double cvmolar() except +
        double speed_sound() except +
        double rhomass() except +
        double hmass() except +
        double smass() except +
        double cpmass() except +
        double cvmass() except +
        
        double keyed_output(long) except+
        double molar_mass() except+
        double gas_constant() except+
        double build_phase_envelope() except+
        double viscosity() except+
        double conductivity() except+
        double surface_tension() except+

# The static factory method for the AbstractState
cdef extern from "AbstractState.h" namespace "CoolProp::AbstractState":
    AbstractState* factory(const string &backend, const string &fluid_string) except+
