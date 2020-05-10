from libcpp cimport bool
from libcpp.string cimport string

# A header defining the AbstractState class
from . cimport cAbstractState

from .typedefs cimport *
from . cimport constants_header

#ctypedef fused string_t:
#    cython.p_char
#    bytes
#    unicode
#    string
#
#ctypedef fused string_or_size_t:
#    string_t
#    cython.integral

ctypedef fused string_or_size_t:
    cython.p_char
    bytes
    unicode
    string
    short
    int
    long

cdef class PyPhaseEnvelopeData:
    cpdef public bool TypeI
    cpdef public size_t iTsat_max, ipsat_max, icrit
    cpdef public list T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap, Q
    cpdef public list x, y, K

cdef class PyGuessesStructure:
    cpdef public double T, p, rhomolar, hmolar, smolar
    cpdef public double rhomolar_liq, rhomolar_vap
    cpdef public list x, y

cdef class PyCriticalState:
    cpdef public double T, p, rhomolar, hmolar, smolar
    cpdef public bool stable

cdef class PySpinodalData:
    cpdef public vector[double] tau, delta, M1

cdef class AbstractState:
    cdef cAbstractState.AbstractState *thisptr     # hold a C++ instance which we're wrapping
    cpdef update(self, constants_header.input_pairs iInput1, double Value1, double Value2)
    cpdef update_with_guesses(self, constants_header.input_pairs iInput1, double Value1, double Value2, PyGuessesStructure guesses)
    cpdef set_mole_fractions(self, vector[double] z)
    cpdef set_mass_fractions(self, vector[double] z)
    cpdef set_volu_fractions(self, vector[double] z)

    cpdef set_binary_interaction_double(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter, double val)
    cpdef double get_binary_interaction_double(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter) except *
    cpdef set_binary_interaction_string(self, string_or_size_t CAS1, string_or_size_t CAS2, string parameter, string val)
    cpdef string get_binary_interaction_string(self, string CAS1, string CAS2, string parameter) except *
    cpdef apply_simple_mixing_rule(self, size_t, size_t, string)

    cpdef name(self)
    cpdef backend_name(self)
    cpdef fluid_names(self)
    cpdef fluid_param_string(self, string key)
    cpdef set_fluid_parameter_double(self, size_t i, string parameter, double value)
    cpdef double get_fluid_parameter_double(self, size_t i, string parameter) except *
    cpdef change_EOS(self, size_t, string)

    cpdef constants_header.phases phase(self) except *
    cpdef specify_phase(self, constants_header.phases phase)
    cpdef unspecify_phase(self)

    ## Limits
    cpdef double Tmin(self) except *
    cpdef double Tmax(self) except *
    cpdef double pmax(self) except *
    cpdef double Ttriple(self) except *

    ## Critical point
    cpdef double T_critical(self) except *
    cpdef double rhomass_critical(self) except *
    cpdef double rhomolar_critical(self) except *
    cpdef double p_critical(self) except *
    cpdef list all_critical_points(self)
    cpdef tuple criticality_contour_values(self)

    ## Spinodal curve(s)
    cpdef build_spinodal(self)
    cpdef PySpinodalData get_spinodal_data(self)

    ## Reducing point
    cpdef double T_reducing(self) except *
    cpdef double rhomolar_reducing(self) except *
    cpdef double rhomass_reducing(self) except *

    ## Tangent plane distance
    cpdef double tangent_plane_distance(self, double T, double p, vector[double] w, double rhomolar_guess=*) except *

    ## ----------------------------------------
    ##        Fluid property accessors
    ## ----------------------------------------

    cpdef double T(self) except *
    cpdef double p(self) except *
    cpdef double Q(self) except *
    cpdef double compressibility_factor(self) except *
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
    cpdef double gibbsmass(self) except *
    cpdef double gibbsmolar(self) except *
    cpdef double helmholtzmass(self) except *
    cpdef double helmholtzmolar(self) except *
    cpdef double speed_sound(self) except *
    cpdef double gas_constant(self) except *
    cpdef double tau(self) except *
    cpdef double delta(self) except *
    cpdef double viscosity(self) except *
    cpdef double conductivity(self) except *
    cpdef dict conformal_state(self, string, CoolPropDbl, CoolPropDbl)
    cpdef dict conductivity_contributions(self)
    cpdef dict viscosity_contributions(self)
    cpdef double surface_tension(self) except *
    cpdef double Prandtl(self) except *
    cpdef double Bvirial(self) except *
    cpdef double Cvirial(self) except *
    cpdef double PIP(self) except *
    cpdef double fundamental_derivative_of_gas_dynamics(self) except *
    cpdef double isothermal_compressibility(self) except *
    cpdef double isobaric_expansion_coefficient(self) except *
    cpdef double fugacity(self, size_t) except *
    cpdef double fugacity_coefficient(self, size_t) except *
    cpdef double chemical_potential(self, size_t) except *

    cpdef double gibbsmolar_excess(self) except *
    cpdef double gibbsmass_excess(self) except *
    cpdef double hmolar_excess(self) except *
    cpdef double hmass_excess(self) except *
    cpdef double smolar_excess(self) except *
    cpdef double smass_excess(self) except *
    cpdef double umolar_excess(self) except *
    cpdef double umass_excess(self) except *
    cpdef double volumemolar_excess(self) except *
    cpdef double volumemass_excess(self) except *
    cpdef double helmholtzmolar_excess(self) except *
    cpdef double helmholtzmass_excess(self) except *

    cpdef double gibbsmolar_residual(self) except *
    cpdef double hmolar_residual(self) except *
    cpdef double smolar_residual(self) except *




    cpdef double molar_mass(self) except *
    cpdef double acentric_factor(self) except*
    cpdef tuple true_critical_point(self)
    cpdef double get_fluid_constant(self,size_t,constants_header.parameters) except*
    cpdef double keyed_output(self, constants_header.parameters) except *
    cpdef double trivial_keyed_output(self, constants_header.parameters) except *
    cpdef double saturated_liquid_keyed_output(self, constants_header.parameters) except *
    cpdef double saturated_vapor_keyed_output(self, constants_header.parameters) except *

    cpdef tuple ideal_curve(self, string)

    ## ----------------------------------------
    ##        Derivatives
    ## ----------------------------------------

    cpdef CoolPropDbl first_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef CoolPropDbl second_partial_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *
    cpdef CoolPropDbl first_saturation_deriv(self, constants_header.parameters, constants_header.parameters) except *
    cpdef CoolPropDbl second_saturation_deriv(self, constants_header.parameters, constants_header.parameters, constants_header.parameters) except *

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
    cpdef get_mass_fractions(self)
    cpdef get_mole_fractions(self)

    cpdef CoolPropDbl alpha0(self) except *
    cpdef CoolPropDbl dalpha0_dDelta(self) except *
    cpdef CoolPropDbl dalpha0_dTau(v) except *
    cpdef CoolPropDbl d2alpha0_dDelta2(self) except *
    cpdef CoolPropDbl d2alpha0_dDelta_dTau(self) except *
    cpdef CoolPropDbl d2alpha0_dTau2(self) except *
    cpdef CoolPropDbl d3alpha0_dTau3(self) except *
    cpdef CoolPropDbl d3alpha0_dDelta_dTau2(self) except *
    cpdef CoolPropDbl d3alpha0_dDelta2_dTau(self) except *
    cpdef CoolPropDbl d3alpha0_dDelta3(self) except *

    cpdef CoolPropDbl alphar(self) except *
    cpdef CoolPropDbl dalphar_dDelta(self) except *
    cpdef CoolPropDbl dalphar_dTau(self) except *
    cpdef CoolPropDbl d2alphar_dDelta2(self) except *
    cpdef CoolPropDbl d2alphar_dDelta_dTau(self) except *
    cpdef CoolPropDbl d2alphar_dTau2(self) except *
    cpdef CoolPropDbl d3alphar_dDelta3(self) except *
    cpdef CoolPropDbl d3alphar_dDelta2_dTau(self) except *
    cpdef CoolPropDbl d3alphar_dDelta_dTau2(self) except *
    cpdef CoolPropDbl d3alphar_dTau3(self) except *
    cpdef CoolPropDbl d4alphar_dDelta4(self) except *
    cpdef CoolPropDbl d4alphar_dDelta3_dTau(self) except *
    cpdef CoolPropDbl d4alphar_dDelta2_dTau2(self) except *
    cpdef CoolPropDbl d4alphar_dDelta_dTau3(self) except *
    cpdef CoolPropDbl d4alphar_dTau4(self) except *
