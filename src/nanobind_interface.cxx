#ifdef NANOBIND

#    include "CoolProp/CoolProp.h"
#    include "CoolProp/AbstractState.h"
#    include "CoolProp/Configuration.h"
#    include "CoolProp/HumidAirProp.h"
#    include "CoolProp/DataStructures.h"
#    include "CoolProp/detail/state_capi.h"
#    include "Backends/Helmholtz/MixtureParameters.h"

#    include <nanobind/nanobind.h>
#    include <nanobind/stl/string.h>
#    include <nanobind/stl/string_view.h>  // get_parameter_information takes std::string_view
#    include <nanobind/stl/vector.h>
#    include <nanobind/stl/tuple.h>
#    include <nanobind/stl/pair.h>
#    include <nanobind/stl/map.h>
#    include <array>
namespace nb = nanobind;

CoolProp::AbstractState* factory(const std::string& backend, const std::string& fluid_names) {
    return CoolProp::AbstractState::factory(backend, fluid_names);
}

// ---- State C-ABI capsule (bridge for the frozen Cython `State` compat shim) --
// The opaque handle carries a shared_ptr<AbstractState> directly.
namespace {
using _CAPI_SP = std::shared_ptr<CoolProp::AbstractState>;
thread_local std::string g_capi_error;  // last C++ error message; "" == none
void* capi_make(const char* backend, const char* fluids) {
    try {
        auto* p = new _CAPI_SP(CoolProp::AbstractState::factory(backend, fluids));
        g_capi_error.clear();
        return p;
    } catch (const std::exception& e) {
        g_capi_error = e.what();
        return nullptr;
    }
}
void capi_destroy(void* h) {
    delete static_cast<_CAPI_SP*>(h);
}
void capi_update(void* h, long input_pair, double v1, double v2) {
    try {
        (*static_cast<_CAPI_SP*>(h))->update(static_cast<CoolProp::input_pairs>(input_pair), v1, v2);
        g_capi_error.clear();
    } catch (const std::exception& e) {
        g_capi_error = e.what();
    }
}
double capi_keyed_output(void* h, long key) {
    try {
        double v = (*static_cast<_CAPI_SP*>(h))->keyed_output(static_cast<CoolProp::parameters>(key));
        g_capi_error.clear();
        return v;
    } catch (const std::exception& e) {
        g_capi_error = e.what();
        return NAN;
    }
}
double capi_first_partial_deriv(void* h, long Of, long Wrt, long Constant) {
    try {
        double v = (*static_cast<_CAPI_SP*>(h))
                     ->first_partial_deriv(static_cast<CoolProp::parameters>(Of), static_cast<CoolProp::parameters>(Wrt),
                                           static_cast<CoolProp::parameters>(Constant));
        g_capi_error.clear();
        return v;
    } catch (const std::exception& e) {
        g_capi_error = e.what();
        return NAN;
    }
}
const char* capi_last_error() {
    return g_capi_error.empty() ? nullptr : g_capi_error.c_str();
}
const CoolProp_StateCAPI g_state_capi = {capi_make, capi_destroy, capi_update, capi_keyed_output, capi_first_partial_deriv, capi_last_error};
}  // namespace

void init_CoolProp(nb::module_& m) {
    using namespace CoolProp;

    nb::class_<SimpleState>(m, "SimpleState")
      .def(nb::init<>())
      .def_rw("T", &SimpleState::T)
      .def_rw("p", &SimpleState::p)
      .def_rw("rhomolar", &SimpleState::rhomolar);

    nb::class_<GuessesStructure>(m, "GuessesStructure")
      .def(nb::init<>())
      .def_rw("T", &GuessesStructure::T)
      .def_rw("p", &GuessesStructure::p)
      .def_rw("rhomolar", &GuessesStructure::rhomolar)
      .def_rw("hmolar", &GuessesStructure::hmolar)
      .def_rw("smolar", &GuessesStructure::smolar)
      .def_rw("rhomolar_liq", &GuessesStructure::rhomolar_liq)
      .def_rw("rhomolar_vap", &GuessesStructure::rhomolar_vap)
      .def_rw("x", &GuessesStructure::x)
      .def_rw("y", &GuessesStructure::y)
      .def("clear", &GuessesStructure::clear);

    nb::class_<CriticalState, SimpleState>(m, "CriticalState").def_rw("stable", &CriticalState::stable);

    nb::class_<PhaseEnvelopeData>(m, "PhaseEnvelopeData")
      .def_rw("K", &PhaseEnvelopeData::K)
      .def_rw("lnK", &PhaseEnvelopeData::lnK)
      .def_rw("x", &PhaseEnvelopeData::x)
      .def_rw("y", &PhaseEnvelopeData::y)
      .def_rw("T", &PhaseEnvelopeData::T)
      .def_rw("p", &PhaseEnvelopeData::p)
      .def_rw("lnT", &PhaseEnvelopeData::lnT)
      .def_rw("lnp", &PhaseEnvelopeData::lnp)
      .def_rw("rhomolar_liq", &PhaseEnvelopeData::rhomolar_liq)
      .def_rw("rhomolar_vap", &PhaseEnvelopeData::rhomolar_vap)
      .def_rw("lnrhomolar_liq", &PhaseEnvelopeData::lnrhomolar_liq)
      .def_rw("lnrhomolar_vap", &PhaseEnvelopeData::lnrhomolar_vap)
      .def_rw("hmolar_liq", &PhaseEnvelopeData::hmolar_liq)
      .def_rw("hmolar_vap", &PhaseEnvelopeData::hmolar_vap)
      .def_rw("smolar_liq", &PhaseEnvelopeData::smolar_liq)
      .def_rw("smolar_vap", &PhaseEnvelopeData::smolar_vap)
      .def_rw("Q", &PhaseEnvelopeData::Q)
      .def_rw("cpmolar_liq", &PhaseEnvelopeData::cpmolar_liq)
      .def_rw("cpmolar_vap", &PhaseEnvelopeData::cpmolar_vap)
      .def_rw("cvmolar_liq", &PhaseEnvelopeData::cvmolar_liq)
      .def_rw("cvmolar_vap", &PhaseEnvelopeData::cvmolar_vap)
      .def_rw("viscosity_liq", &PhaseEnvelopeData::viscosity_liq)
      .def_rw("viscosity_vap", &PhaseEnvelopeData::viscosity_vap)
      .def_rw("conductivity_liq", &PhaseEnvelopeData::conductivity_liq)
      .def_rw("conductivity_vap", &PhaseEnvelopeData::conductivity_vap)
      .def_rw("speed_sound_vap", &PhaseEnvelopeData::speed_sound_vap);

    // See http://stackoverflow.com/a/148610 and http://stackoverflow.com/questions/147267/easy-way-to-use-variables-of-enum-types-as-string-in-c#202511
    nb::enum_<configuration_keys>(m, "configuration_keys", nb::is_arithmetic())
#    define X(Enum, String, Default, Desc) .value(String, configuration_keys::Enum)
      CONFIGURATION_KEYS_ENUM
#    undef X
        .export_values();

    nb::enum_<parameters>(m, "parameters", nb::is_arithmetic())
      .value("igas_constant", parameters::igas_constant)
      .value("imolar_mass", parameters::imolar_mass)
      .value("iacentric_factor", parameters::iacentric_factor)
      .value("irhomolar_reducing", parameters::irhomolar_reducing)
      .value("irhomolar_critical", parameters::irhomolar_critical)
      .value("iT_reducing", parameters::iT_reducing)
      .value("iT_critical", parameters::iT_critical)
      .value("irhomass_reducing", parameters::irhomass_reducing)
      .value("irhomass_critical", parameters::irhomass_critical)
      .value("iP_critical", parameters::iP_critical)
      .value("iP_reducing", parameters::iP_reducing)
      .value("iT_triple", parameters::iT_triple)
      .value("iP_triple", parameters::iP_triple)
      .value("iT_min", parameters::iT_min)
      .value("iT_max", parameters::iT_max)
      .value("iP_max", parameters::iP_max)
      .value("iP_min", parameters::iP_min)
      .value("idipole_moment", parameters::idipole_moment)
      .value("iT", parameters::iT)
      .value("iP", parameters::iP)
      .value("iQ", parameters::iQ)
      .value("iTau", parameters::iTau)
      .value("iDelta", parameters::iDelta)
      .value("iDmolar", parameters::iDmolar)
      .value("iHmolar", parameters::iHmolar)
      .value("iSmolar", parameters::iSmolar)
      .value("iCpmolar", parameters::iCpmolar)
      .value("iCp0molar", parameters::iCp0molar)
      .value("iCvmolar", parameters::iCvmolar)
      .value("iUmolar", parameters::iUmolar)
      .value("iGmolar", parameters::iGmolar)
      .value("iHelmholtzmolar", parameters::iHelmholtzmolar)
      .value("iSmolar_residual", parameters::iSmolar_residual)
      .value("iHmolar_residual", parameters::iHmolar_residual)
      .value("iGmolar_residual", parameters::iGmolar_residual)
      .value("iDmass", parameters::iDmass)
      .value("iHmass", parameters::iHmass)
      .value("iSmass", parameters::iSmass)
      .value("iCpmass", parameters::iCpmass)
      .value("iCp0mass", parameters::iCp0mass)
      .value("iCvmass", parameters::iCvmass)
      .value("iUmass", parameters::iUmass)
      .value("iGmass", parameters::iGmass)
      .value("iHelmholtzmass", parameters::iHelmholtzmass)
      .value("iviscosity", parameters::iviscosity)
      .value("iconductivity", parameters::iconductivity)
      .value("isurface_tension", parameters::isurface_tension)
      .value("iPrandtl", parameters::iPrandtl)
      .value("ispeed_sound", parameters::ispeed_sound)
      .value("iisothermal_compressibility", parameters::iisothermal_compressibility)
      .value("iisobaric_expansion_coefficient", parameters::iisobaric_expansion_coefficient)
      .value("ifundamental_derivative_of_gas_dynamics", parameters::ifundamental_derivative_of_gas_dynamics)
      .value("ialphar", parameters::ialphar)
      .value("idalphar_ddelta_consttau", parameters::idalphar_ddelta_consttau)
      .value("idalpha0_dtau_constdelta", parameters::idalpha0_dtau_constdelta)
      .value("iBvirial", parameters::iBvirial)
      .value("iCvirial", parameters::iCvirial)
      .value("idBvirial_dT", parameters::idBvirial_dT)
      .value("idCvirial_dT", parameters::idCvirial_dT)
      .value("iZ", parameters::iZ)
      .value("iPIP", parameters::iPIP)
      .value("ifraction_min", parameters::ifraction_min)
      .value("ifraction_max", parameters::ifraction_max)
      .value("iT_freeze", parameters::iT_freeze)
      .value("iGWP20", parameters::iGWP20)
      .value("iGWP100", parameters::iGWP100)
      .value("iGWP500", parameters::iGWP500)
      .value("iFH", parameters::iFH)
      .value("iHH", parameters::iHH)
      .value("iPH", parameters::iPH)
      .value("iODP", parameters::iODP)
      .value("iPhase", parameters::iPhase)
      .value("iundefined_parameter", parameters::iundefined_parameter)
      .export_values();

    nb::enum_<input_pairs>(m, "input_pairs", nb::is_arithmetic())
      .value("QT_INPUTS", input_pairs::QT_INPUTS)
      .value("PQ_INPUTS", input_pairs::PQ_INPUTS)
      .value("QSmolar_INPUTS", input_pairs::QSmolar_INPUTS)
      .value("QSmass_INPUTS", input_pairs::QSmass_INPUTS)
      .value("HmolarQ_INPUTS", input_pairs::HmolarQ_INPUTS)
      .value("HmassQ_INPUTS", input_pairs::HmassQ_INPUTS)
      .value("DmolarQ_INPUTS", input_pairs::DmolarQ_INPUTS)
      .value("DmassQ_INPUTS", input_pairs::DmassQ_INPUTS)
      .value("PT_INPUTS", input_pairs::PT_INPUTS)
      .value("DmassT_INPUTS", input_pairs::DmassT_INPUTS)
      .value("DmolarT_INPUTS", input_pairs::DmolarT_INPUTS)
      .value("HmolarT_INPUTS", input_pairs::HmolarT_INPUTS)
      .value("HmassT_INPUTS", input_pairs::HmassT_INPUTS)
      .value("SmolarT_INPUTS", input_pairs::SmolarT_INPUTS)
      .value("SmassT_INPUTS", input_pairs::SmassT_INPUTS)
      .value("TUmolar_INPUTS", input_pairs::TUmolar_INPUTS)
      .value("TUmass_INPUTS", input_pairs::TUmass_INPUTS)
      .value("DmassP_INPUTS", input_pairs::DmassP_INPUTS)
      .value("DmolarP_INPUTS", input_pairs::DmolarP_INPUTS)
      .value("HmassP_INPUTS", input_pairs::HmassP_INPUTS)
      .value("HmolarP_INPUTS", input_pairs::HmolarP_INPUTS)
      .value("PSmass_INPUTS", input_pairs::PSmass_INPUTS)
      .value("PSmolar_INPUTS", input_pairs::PSmolar_INPUTS)
      .value("PUmass_INPUTS", input_pairs::PUmass_INPUTS)
      .value("PUmolar_INPUTS", input_pairs::PUmolar_INPUTS)
      .value("HmassSmass_INPUTS", input_pairs::HmassSmass_INPUTS)
      .value("HmolarSmolar_INPUTS", input_pairs::HmolarSmolar_INPUTS)
      .value("SmassUmass_INPUTS", input_pairs::SmassUmass_INPUTS)
      .value("SmolarUmolar_INPUTS", input_pairs::SmolarUmolar_INPUTS)
      .value("DmassHmass_INPUTS", input_pairs::DmassHmass_INPUTS)
      .value("DmolarHmolar_INPUTS", input_pairs::DmolarHmolar_INPUTS)
      .value("DmassSmass_INPUTS", input_pairs::DmassSmass_INPUTS)
      .value("DmolarSmolar_INPUTS", input_pairs::DmolarSmolar_INPUTS)
      .value("DmassUmass_INPUTS", input_pairs::DmassUmass_INPUTS)
      .value("DmolarUmolar_INPUTS", input_pairs::DmolarUmolar_INPUTS)
      .export_values();

    nb::enum_<phases>(m, "phases", nb::is_arithmetic())
      .value("iphase_liquid", phases::iphase_liquid)
      .value("iphase_supercritical", phases::iphase_supercritical)
      .value("iphase_supercritical_gas", phases::iphase_supercritical_gas)
      .value("iphase_supercritical_liquid", phases::iphase_supercritical_liquid)
      .value("iphase_critical_point", phases::iphase_critical_point)
      .value("iphase_gas", phases::iphase_gas)
      .value("iphase_twophase", phases::iphase_twophase)
      .value("iphase_unknown", phases::iphase_unknown)
      .export_values();

    nb::class_<AbstractState>(m, "_AbstractState")
      .def("set_T", &AbstractState::set_T)
      .def("backend_name", &AbstractState::backend_name)
      .def("using_mole_fractions", &AbstractState::using_mole_fractions)
      .def("using_mass_fractions", &AbstractState::using_mass_fractions)
      .def("using_volu_fractions", &AbstractState::using_volu_fractions)
      .def("set_mole_fractions", &AbstractState::set_mole_fractions)
      .def("set_mass_fractions", &AbstractState::set_mass_fractions)
      .def("set_volu_fractions", &AbstractState::set_volu_fractions)
      .def("mole_fractions_liquid", &AbstractState::mole_fractions_liquid)
      .def("mole_fractions_liquid_double", &AbstractState::mole_fractions_liquid_double)
      .def("mole_fractions_vapor", &AbstractState::mole_fractions_vapor)
      .def("mole_fractions_vapor_double", &AbstractState::mole_fractions_vapor_double)
      .def("get_mole_fractions", &AbstractState::get_mole_fractions)
      .def("get_mass_fractions", &AbstractState::get_mass_fractions)
      .def("update", &AbstractState::update)
      .def("update_with_guesses", &AbstractState::update_with_guesses)
      .def("available_in_high_level", &AbstractState::available_in_high_level)
      .def("fluid_param_string", &AbstractState::fluid_param_string)
      .def("fluid_names", &AbstractState::fluid_names)
      .def("set_binary_interaction_double", (void(AbstractState::*)(const std::string&, const std::string&, const std::string&, const double))
                                              & AbstractState::set_binary_interaction_double)
      .def("set_binary_interaction_double", (void(AbstractState::*)(const std::size_t, const std::size_t, const std::string&, const double))
                                              & AbstractState::set_binary_interaction_double)
      .def("set_binary_interaction_string", (void(AbstractState::*)(const std::string&, const std::string&, const std::string&, const std::string&))
                                              & AbstractState::set_binary_interaction_string)
      .def("set_binary_interaction_string", (void(AbstractState::*)(const std::size_t, const std::size_t, const std::string&, const std::string&))
                                              & AbstractState::set_binary_interaction_string)
      .def("get_binary_interaction_double",
           (double(AbstractState::*)(const std::string&, const std::string&, const std::string&)) & AbstractState::get_binary_interaction_double)
      .def("get_binary_interaction_double",
           (double(AbstractState::*)(const std::size_t, const std::size_t, const std::string&)) & AbstractState::get_binary_interaction_double)
      .def("get_binary_interaction_string", &AbstractState::get_binary_interaction_string)
      .def("apply_simple_mixing_rule", &AbstractState::apply_simple_mixing_rule)
      .def("set_fluid_parameter_double", &AbstractState::set_fluid_parameter_double)
      .def("clear", &AbstractState::clear)
      .def("get_reducing_state", &AbstractState::get_reducing_state)
      .def("get_state", &AbstractState::get_state)
      .def("Tmin", &AbstractState::Tmin)
      .def("Tmax", &AbstractState::Tmax)
      .def("pmax", &AbstractState::pmax)
      .def("Ttriple", &AbstractState::Ttriple)
      .def("phase", &AbstractState::phase)
      .def("specify_phase", &AbstractState::specify_phase)
      .def("unspecify_phase", &AbstractState::unspecify_phase)
      .def("T_critical", &AbstractState::T_critical)
      .def("p_critical", &AbstractState::p_critical)
      .def("rhomolar_critical", &AbstractState::rhomolar_critical)
      .def("rhomass_critical", &AbstractState::rhomass_critical)
      .def("all_critical_points", &AbstractState::all_critical_points)
      .def("build_spinodal", &AbstractState::build_spinodal)
      .def("get_spinodal_data", &AbstractState::get_spinodal_data)
      .def("criticality_contour_values",
           [](AbstractState& AS) {
               double L, M;
               AS.criticality_contour_values(L, M);
               return nb::make_tuple(L, M);
           })
      .def("tangent_plane_distance", &AbstractState::tangent_plane_distance)
      .def("T_reducing", &AbstractState::T_reducing)
      .def("rhomolar_reducing", &AbstractState::rhomolar_reducing)
      .def("rhomass_reducing", &AbstractState::rhomass_reducing)
      .def("p_triple", &AbstractState::p_triple)
      .def("name", &AbstractState::name)
      .def("dipole_moment", &AbstractState::dipole_moment)
      .def("keyed_output", &AbstractState::keyed_output)
      .def("trivial_keyed_output", &AbstractState::trivial_keyed_output)
      .def("saturated_liquid_keyed_output", &AbstractState::saturated_liquid_keyed_output)
      .def("saturated_vapor_keyed_output", &AbstractState::saturated_vapor_keyed_output)
      .def("T", &AbstractState::T)
      .def("rhomolar", &AbstractState::rhomolar)
      .def("rhomass", &AbstractState::rhomass)
      .def("p", &AbstractState::p)
      .def("Q", &AbstractState::Q)
      .def("Qmass", &AbstractState::Qmass)
      .def("tau", &AbstractState::tau)
      .def("delta", &AbstractState::delta)
      .def("molar_mass", &AbstractState::molar_mass)
      .def("acentric_factor", &AbstractState::acentric_factor)
      .def("gas_constant", &AbstractState::gas_constant)
      .def("Bvirial", &AbstractState::Bvirial)
      .def("dBvirial_dT", &AbstractState::dBvirial_dT)
      .def("Cvirial", &AbstractState::Cvirial)
      .def("dCvirial_dT", &AbstractState::dCvirial_dT)
      .def("compressibility_factor", &AbstractState::compressibility_factor)
      .def("hmolar", &AbstractState::hmolar)
      .def("hmass", &AbstractState::hmass)
      .def("hmolar_excess", &AbstractState::hmolar_excess)
      .def("hmass_excess", &AbstractState::hmass_excess)
      .def("smolar", &AbstractState::smolar)
      .def("smass", &AbstractState::smass)
      .def("smolar_excess", &AbstractState::smolar_excess)
      .def("smass_excess", &AbstractState::smass_excess)
      .def("umolar", &AbstractState::umolar)
      .def("umass", &AbstractState::umass)
      .def("umolar_excess", &AbstractState::umolar_excess)
      .def("umass_excess", &AbstractState::umass_excess)
      .def("cpmolar", &AbstractState::cpmolar)
      .def("cpmass", &AbstractState::cpmass)
      .def("cp0molar", &AbstractState::cp0molar)
      .def("cp0mass", &AbstractState::cp0mass)
      .def("cvmolar", &AbstractState::cvmolar)
      .def("cvmass", &AbstractState::cvmass)
      .def("gibbsmolar", &AbstractState::gibbsmolar)
      .def("gibbsmass", &AbstractState::gibbsmass)
      .def("gibbsmolar_excess", &AbstractState::gibbsmolar_excess)
      .def("gibbsmass_excess", &AbstractState::gibbsmass_excess)
      .def("helmholtzmolar", &AbstractState::helmholtzmolar)
      .def("helmholtzmass", &AbstractState::helmholtzmass)
      .def("helmholtzmolar_excess", &AbstractState::helmholtzmolar_excess)
      .def("helmholtzmass_excess", &AbstractState::helmholtzmass_excess)
      .def("volumemolar_excess", &AbstractState::volumemolar_excess)
      .def("volumemass_excess", &AbstractState::volumemass_excess)
      .def("speed_sound", &AbstractState::speed_sound)
      .def("isothermal_compressibility", &AbstractState::isothermal_compressibility)
      .def("isobaric_expansion_coefficient", &AbstractState::isobaric_expansion_coefficient)
      .def("fugacity_coefficient", &AbstractState::fugacity_coefficient)
      .def("fugacity", &AbstractState::fugacity)
      .def("chemical_potential", &AbstractState::chemical_potential)
      .def("fundamental_derivative_of_gas_dynamics", &AbstractState::fundamental_derivative_of_gas_dynamics)
      .def("PIP", &AbstractState::PIP)
      .def("true_critical_point", &AbstractState::true_critical_point)
      .def("ideal_curve",
           [](AbstractState& AS, const std::string& name) {
               std::vector<double> T, p;
               AS.ideal_curve(name, T, p);
               return nb::make_tuple(T, p);
           })
      .def("first_partial_deriv", &AbstractState::first_partial_deriv)
      .def("second_partial_deriv", &AbstractState::second_partial_deriv)
      .def("first_saturation_deriv", &AbstractState::first_saturation_deriv)
      .def("second_saturation_deriv", &AbstractState::second_saturation_deriv)
      .def("first_two_phase_deriv", &AbstractState::first_two_phase_deriv)
      .def("second_two_phase_deriv", &AbstractState::second_two_phase_deriv)
      .def("first_two_phase_deriv_splined", &AbstractState::first_two_phase_deriv_splined)
      .def("build_phase_envelope", &AbstractState::build_phase_envelope)
      .def("get_phase_envelope_data", &AbstractState::get_phase_envelope_data)
      .def("has_melting_line", &AbstractState::has_melting_line)
      .def("melting_line", &AbstractState::melting_line)
      .def("saturation_ancillary", &AbstractState::saturation_ancillary)
      .def("viscosity", &AbstractState::viscosity)
      .def("viscosity_contributions", &AbstractState::viscosity_contributions)
      .def("conductivity", &AbstractState::conductivity)
      .def("conductivity_contributions", &AbstractState::conductivity_contributions)
      .def("surface_tension", &AbstractState::surface_tension)
      .def("Prandtl", &AbstractState::Prandtl)
      .def("conformal_state", &AbstractState::conformal_state)
      .def("change_EOS", &AbstractState::change_EOS)
      .def("alpha0", &AbstractState::alpha0)
      .def("dalpha0_dDelta", &AbstractState::dalpha0_dDelta)
      .def("dalpha0_dTau", &AbstractState::dalpha0_dTau)
      .def("d2alpha0_dDelta2", &AbstractState::d2alpha0_dDelta2)
      .def("d2alpha0_dDelta_dTau", &AbstractState::d2alpha0_dDelta_dTau)
      .def("d2alpha0_dTau2", &AbstractState::d2alpha0_dTau2)
      .def("d3alpha0_dTau3", &AbstractState::d3alpha0_dTau3)
      .def("d3alpha0_dDelta_dTau2", &AbstractState::d3alpha0_dDelta_dTau2)
      .def("d3alpha0_dDelta2_dTau", &AbstractState::d3alpha0_dDelta2_dTau)
      .def("d3alpha0_dDelta3", &AbstractState::d3alpha0_dDelta3)
      .def("alphar", &AbstractState::alphar)
      .def("dalphar_dDelta", &AbstractState::dalphar_dDelta)
      .def("dalphar_dTau", &AbstractState::dalphar_dTau)
      .def("d2alphar_dDelta2", &AbstractState::d2alphar_dDelta2)
      .def("d2alphar_dDelta_dTau", &AbstractState::d2alphar_dDelta_dTau)
      .def("d2alphar_dTau2", &AbstractState::d2alphar_dTau2)
      .def("d3alphar_dDelta3", &AbstractState::d3alphar_dDelta3)
      .def("d3alphar_dDelta2_dTau", &AbstractState::d3alphar_dDelta2_dTau)
      .def("d3alphar_dDelta_dTau2", &AbstractState::d3alphar_dDelta_dTau2)
      .def("d3alphar_dTau3", &AbstractState::d3alphar_dTau3)
      .def("d4alphar_dDelta4", &AbstractState::d4alphar_dDelta4)
      .def("d4alphar_dDelta3_dTau", &AbstractState::d4alphar_dDelta3_dTau)
      .def("d4alphar_dDelta2_dTau2", &AbstractState::d4alphar_dDelta2_dTau2)
      .def("d4alphar_dDelta_dTau3", &AbstractState::d4alphar_dDelta_dTau3)
      .def("d4alphar_dTau4", &AbstractState::d4alphar_dTau4);

    m.def("AbstractState", &factory);

    m.def("get_config_as_json_string", &get_config_as_json_string);
    m.def("set_config_as_json_string", &set_config_as_json_string);
    m.def("config_key_description", (std::string(*)(configuration_keys)) & config_key_description);
    m.def("config_key_description", (std::string(*)(const std::string&)) & config_key_description);
    m.def("set_config_string", &set_config_string);
    m.def("set_config_double", &set_config_double);
    m.def("set_departure_functions", &set_departure_functions);
    m.def("set_config_bool", &set_config_bool);
    m.def("get_config_string", &get_config_string);
    m.def("get_config_double", &get_config_double);
    m.def("get_config_bool", &get_config_bool);
    m.def("get_parameter_information", &get_parameter_information);
    m.def("get_parameter_index", &get_parameter_index);
    m.def("get_phase_index", &get_phase_index);
    m.def("is_trivial_parameter", &is_trivial_parameter);
    m.def("generate_update_pair", &generate_update_pair<double>);
    m.def("Props1SI", &Props1SI);
    m.def("PropsSI", &PropsSI);
    m.def("PhaseSI", &PhaseSI);
    m.def("PropsSImulti", &PropsSImulti);
    m.def("get_global_param_string", &get_global_param_string);
    m.def("get_debug_level", &get_debug_level);
    m.def("set_debug_level", &set_debug_level);
    m.def("get_fluid_param_string", &get_fluid_param_string);
    m.def("extract_backend", &extract_backend);
    m.def("extract_fractions", &extract_fractions);
    m.def("set_reference_stateS", &set_reference_stateS);
    m.def("set_reference_stateD", &set_reference_stateD);
    m.def("saturation_ancillary", &saturation_ancillary);
    m.def("add_fluids_as_JSON", &add_fluids_as_JSON);
    m.def("HAPropsSI", &HumidAir::HAPropsSI);
    m.def("HAProps", &HumidAir::HAProps);
    m.def("HAProps_Aux", [](std::string out_string, double T, double p, double psi_w) {
        std::array<char, 1000> units{};
        double out = HumidAir::HAProps_Aux(out_string.c_str(), T, p, psi_w, units.data());
        return nb::make_tuple(out, std::string(units.data()));
    });
    m.def("cair_sat", &HumidAir::cair_sat);
    m.def("get_mixture_binary_pair_data", &get_mixture_binary_pair_data);
    m.def("set_mixture_binary_pair_data", &set_mixture_binary_pair_data);
    m.def("apply_simple_mixing_rule", &apply_simple_mixing_rule);
}

#    if defined(COOLPROP_NANOBIND_MODULE)
NB_MODULE(CoolProp, m) {
    init_CoolProp(m);
    // Export the State C-ABI table so the frozen Cython `State` shim can forward
    // through it (PDSim cimport compatibility without a Cython AbstractState).
    m.attr("_capi") = nb::capsule(&g_state_capi, "CoolProp._capi");
}
#    endif

#endif
