from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from . cimport constants_header

from .typedefs cimport CoolPropDbl

cdef extern from "PhaseEnvelope.h" namespace "CoolProp":
    cdef cppclass PhaseEnvelopeData:
        bool TypeI
        size_t iTsat_max, ipsat_max, icrit
        vector[double] T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q
        vector[vector[double]] x, y, K

cdef extern from "DataStructures.h" namespace "CoolProp":
    cdef cppclass CriticalState:
        double T, p, rhomolar, hmolar, smolar
        bool stable

cdef extern from "AbstractState.h" namespace "CoolProp":

    cdef cppclass GuessesStructure:
        double T, p, rhomolar, hmolar, smolar
        double rhomolar_liq, rhomolar_vap
        double _rhomolar_liq, _rhomolar_vap
        vector[double] x, y

    cdef cppclass SpinodalData:
        vector[double] tau, delta, M1

    cdef cppclass AbstractState:

        ## Nullary Constructor
        AbstractState() except +ValueError

        ## Constructor with fluid name
        AbstractState(string FluidName) except +ValueError

        void set_mole_fractions(vector[double]) except+ValueError
        void set_mass_fractions(vector[double]) except+ValueError
        void set_volu_fractions(vector[double]) except+ValueError

        vector[CoolPropDbl] mole_fractions_liquid() except +ValueError
        vector[CoolPropDbl] mole_fractions_vapor() except +ValueError
        vector[CoolPropDbl] get_mole_fractions() except +ValueError
        vector[CoolPropDbl] get_mass_fractions() except +ValueError

        constants_header.phases phase() except +ValueError
        void specify_phase(constants_header.phases phase) except +ValueError
        void unspecify_phase() except +ValueError

        void change_EOS(const size_t, const string &) except +ValueError

        void set_binary_interaction_double(const string, const string &, const string &, const double s) except +ValueError
        void set_binary_interaction_string(const string &, const string &, const string &, const string &) except +ValueError
        double get_binary_interaction_double(const string &, const string &, const string &) except +ValueError
        string get_binary_interaction_string(const string &, const string &, const string &) except +ValueError
        void set_binary_interaction_string(const size_t, const size_t, const string &, const string &) except +ValueError
        void apply_simple_mixing_rule(size_t, size_t, const string &) except +ValueError

        double get_binary_interaction_double(const size_t, const size_t, const string &) except +ValueError
        void set_binary_interaction_double(const size_t, const size_t, const string &, const double s) except +ValueError

        string name() except +ValueError
        string backend_name() except +ValueError
        vector[string] fluid_names() except +ValueError
        string fluid_param_string(const string &) except +ValueError
        void set_cubic_alpha_C(const size_t, const string&, const double, const double, const double) except +ValueError
        void set_fluid_parameter_double(const size_t, const string&, const double) except +ValueError
        double get_fluid_parameter_double(const size_t, const string&) except +ValueError

        bool clear()

        ## Limits
        double Tmin() except +ValueError
        double Tmax() except +ValueError
        double pmax() except +ValueError
        double Ttriple() except +ValueError

        ## Critical point
        double T_critical() except +ValueError
        double rhomass_critical() except +ValueError
        double rhomolar_critical() except +ValueError
        double p_critical() except +ValueError
        vector[CriticalState] all_critical_points() except +ValueError
        void criticality_contour_values(double &L1star, double &M1star) except +ValueError

        ## Spinodal curve
        void build_spinodal() except +ValueError
        SpinodalData get_spinodal_data() except +ValueError

        ## Tangent plane analysis
        double tangent_plane_distance(const double, const double, const vector[double], const double) except +ValueError

        ## Reducing point
        double T_reducing() except +ValueError
        double rhomolar_reducing() except +ValueError
        double rhomass_reducing() except +ValueError

        void ideal_curve(const string &, vector[double] &T, vector[double] &p) except +ValueError

        ## Property updater
        ## Uses the indices in CoolProp for the input parameters
        void update(constants_header.input_pairs iInput1, double Value1, double Value2) except +ValueError
        ## Uses the indices in CoolProp for the input parameters
        void update_with_guesses(constants_header.input_pairs iInput1, double Value1, double Value2, GuessesStructure) except +ValueError

        ## Bulk properties accessors - temperature, pressure and density are directly calculated every time
        ## All other parameters are calculated on an as-needed basis
        ## If single-phase, just plug into the EOS, otherwise need to do two-phase analysis
        double T() except +ValueError
        double rhomolar() except +ValueError
        double rhomass() except +ValueError
        double p() except +ValueError
        double compressibility_factor() except +ValueError
        double Q() except +ValueError
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
        double gibbsmolar() except +ValueError
        double gibbsmass() except +ValueError
        double helmholtzmolar() except +ValueError
        double helmholtzmass() except +ValueError
        double speed_sound() except +ValueError
        double tau() except +ValueError
        double delta() except +ValueError
        double viscosity() except+ValueError
        double conductivity() except+ValueError
        void conformal_state(const string &, CoolPropDbl &, CoolPropDbl &) except +ValueError
        void conductivity_contributions(CoolPropDbl &dilute, CoolPropDbl &initial_density, CoolPropDbl &residual, CoolPropDbl &critical) except +ValueError
        void viscosity_contributions(CoolPropDbl &dilute, CoolPropDbl &initial_density, CoolPropDbl &residual, CoolPropDbl &critical) except +ValueError

        double gibbsmolar_excess() except +ValueError
        double gibbsmass_excess() except +ValueError
        double hmolar_excess() except +ValueError
        double hmass_excess() except +ValueError
        double smolar_excess() except +ValueError
        double smass_excess() except +ValueError
        double umolar_excess() except +ValueError
        double umass_excess() except +ValueError
        double volumemolar_excess() except +ValueError
        double volumemass_excess() except +ValueError
        double helmholtzmolar_excess() except +ValueError
        double helmholtzmass_excess() except +ValueError

        double gibbsmolar_residual() except +ValueError
        double hmolar_residual() except +ValueError
        double smolar_residual() except +ValueError
        double neff() except +ValueError

        double surface_tension() except +ValueError
        double Prandtl() except +ValueError
        double Bvirial() except +ValueError
        double Cvirial() except +ValueError
        double PIP() except +ValueError
        double fundamental_derivative_of_gas_dynamics() except +ValueError
        double isothermal_compressibility() except +ValueError
        double isobaric_expansion_coefficient() except +ValueError
        double fugacity(size_t) except +ValueError
        double fugacity_coefficient(size_t) except +ValueError
        double chemical_potential(size_t) except +ValueError

        double get_fluid_constant(size_t,constants_header.parameters) except+ValueError
        double keyed_output(constants_header.parameters) except+ValueError
        double trivial_keyed_output(constants_header.parameters) except+ValueError
        double saturated_liquid_keyed_output(constants_header.parameters) except+ValueError
        double saturated_vapor_keyed_output(constants_header.parameters) except+ValueError

        double molar_mass() except+ValueError
        double acentric_factor() except+ValueError
        double gas_constant() except+ValueError

        CoolPropDbl first_partial_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        CoolPropDbl second_partial_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        CoolPropDbl first_saturation_deriv(constants_header.parameters, constants_header.parameters) except+ValueError
        CoolPropDbl second_saturation_deriv(constants_header.parameters, constants_header.parameters, constants_header.parameters) except+ValueError
        double first_two_phase_deriv(constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant) except+ValueError
        double second_two_phase_deriv(constants_header.parameters Of, constants_header.parameters Wrt1, constants_header.parameters Constant1, constants_header.parameters Wrt2, constants_header.parameters Constant2) except+ValueError
        double first_two_phase_deriv_splined(constants_header.parameters Of, constants_header.parameters Wrt, constants_header.parameters Constant, double x_end) except+ValueError
        void true_critical_point(double &T, double &rho) except +ValueError

        double melting_line(int,int,double) except+ValueError
        bool has_melting_line() except+ValueError
        double saturation_ancillary(constants_header.parameters, int, constants_header.parameters, double) except +ValueError

        double build_phase_envelope() except+ValueError
        void build_phase_envelope(string) except+ValueError
        PhaseEnvelopeData get_phase_envelope_data() except+ValueError

        CoolPropDbl alpha0() except+ValueError
        CoolPropDbl dalpha0_dDelta() except+ValueError
        CoolPropDbl dalpha0_dTau() except+ValueError
        CoolPropDbl d2alpha0_dDelta2() except+ValueError
        CoolPropDbl d2alpha0_dDelta_dTau() except+ValueError
        CoolPropDbl d2alpha0_dTau2() except+ValueError
        CoolPropDbl d3alpha0_dTau3() except+ValueError
        CoolPropDbl d3alpha0_dDelta_dTau2() except+ValueError
        CoolPropDbl d3alpha0_dDelta2_dTau() except+ValueError
        CoolPropDbl d3alpha0_dDelta3() except+ValueError

        CoolPropDbl alphar() except+ValueError
        CoolPropDbl dalphar_dDelta() except+ValueError
        CoolPropDbl dalphar_dTau() except+ValueError
        CoolPropDbl d2alphar_dDelta2() except+ValueError
        CoolPropDbl d2alphar_dDelta_dTau() except+ValueError
        CoolPropDbl d2alphar_dTau2() except+ValueError
        CoolPropDbl d3alphar_dDelta3() except+ValueError
        CoolPropDbl d3alphar_dDelta2_dTau() except+ValueError
        CoolPropDbl d3alphar_dDelta_dTau2() except+ValueError
        CoolPropDbl d3alphar_dTau3() except+ValueError
        CoolPropDbl d4alphar_dDelta4() except+ValueError
        CoolPropDbl d4alphar_dDelta3_dTau() except+ValueError
        CoolPropDbl d4alphar_dDelta2_dTau2() except+ValueError
        CoolPropDbl d4alphar_dDelta_dTau3() except+ValueError
        CoolPropDbl d4alphar_dTau4() except+ValueError


# The static factory method for the AbstractState
cdef extern from "AbstractState.h" namespace "CoolProp::AbstractState":
    AbstractState* factory(const string &backend, const string &fluid_string) except+ValueError
