from libcpp cimport bool 
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport constants_header

cdef extern from "PhaseEnvelope.h" namespace "CoolProp":
    cdef cppclass PhaseEnvelopeData:
        bool TypeI
        size_t iTsat_max, ipsat_max, icrit
        vector[long double] T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q
    
cdef extern from "AbstractState.h" namespace "CoolProp":
    cdef cppclass AbstractState:
        
        ## Nullary Constructor
        AbstractState() except +ValueError
        
        ## Constructor with fluid name
        AbstractState(string FluidName) except +ValueError
        
        void specify_phase(constants_header.phases phase) except +ValueError
        
        void unspecify_phase() except +ValueError
        
        bool clear()

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(constants_header.input_pairs iInput1, double Value1, double Value2) except +ValueError

        ## Bulk properties accessors - temperature, pressure and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T() except +ValueError
        double rhomolar() except +ValueError
        double rhomass() except +ValueError
        double p() except +ValueError
        double hmolar() except +ValueError
        double hmass() except +ValueError
        double smolar() except +ValueError
        double smass() except +ValueError
        double umolar() except +ValueError
        double umass() except +ValueError
        double cpmolar() except +ValueError
        double cpmass() except +ValueError
        double cp0molar() except +ValueError
        double cp0mass() except +ValueError
        double cvmolar() except +ValueError
        double cvmass() except +ValueError
        double speed_sound() except +ValueError
        
        double keyed_output(constants_header.parameters) except+ValueError
        double molar_mass() except+ValueError
        double gas_constant() except+ValueError
        double build_phase_envelope() except+ValueError
        
        double viscosity() except+ValueError
        double conductivity() except+ValueError
        double surface_tension() except+ValueError
        
        long double first_partial_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        long double second_partial_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        long double first_saturation_deriv(constants_header.parameters, constants_header.parameters) except+ValueError
        long double second_saturation_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        
        void set_mole_fractions(vector[double]) except+ValueError
        void set_mass_fractions(vector[double]) except+ValueError
        void set_volu_fractions(vector[double]) except+ValueError
        
        double melting_line(int,int,double) except+ValueError
        bool has_melting_line() except+ValueError
        
        void build_phase_envelope(string) except+ValueError
        PhaseEnvelopeData get_phase_envelope_data() except+ValueError
        
        vector[long double] mole_fractions_liquid() except +ValueError
        vector[long double] mole_fractions_vapor() except +ValueError

# The static factory method for the AbstractState
cdef extern from "AbstractState.h" namespace "CoolProp::AbstractState":
    AbstractState* factory(const string &backend, const string &fluid_string) except+ValueError
