
/// *********************************************************************************
/// *********************************************************************************
///                     EMSCRIPTEN (for javascript)
/// *********************************************************************************
/// *********************************************************************************

#ifdef EMSCRIPTEN

#    include "CoolProp.h"
#    include "AbstractState.h"
#    include "Configuration.h"
#    include "HumidAirProp.h"
#    include "DataStructures.h"
#    include "Backends/Helmholtz/MixtureParameters.h"
#    include "CoolPropLib.h"

/// *********************************************************************************
/// *********************************************************************************
///                     EMSCRIPTEN (for javascript)
/// *********************************************************************************
/// *********************************************************************************

#    include <emscripten/bind.h>
using namespace emscripten;

// Binding code
EMSCRIPTEN_BINDINGS(coolprop_bindings) {
    function("F2K", &F2K);
    function("Props1SI", &CoolProp::Props1SI);
    function("PropsSI", &CoolProp::PropsSI);
    function("get_global_param_string", &CoolProp::get_global_param_string);
    function("get_fluid_param_string", &CoolProp::get_fluid_param_string);
    function("apply_simple_mixing_rule", &CoolProp::apply_simple_mixing_rule);
    function("get_mixture_binary_pair_data", &CoolProp::get_mixture_binary_pair_data);
    function("add_fluids_as_JSON", &CoolProp::add_fluids_as_JSON);

    enum_<CoolProp::parameters>("parameters")
        .value("igas_constant", CoolProp::igas_constant)
        .value("imolar_mass", CoolProp::imolar_mass)
        .value("iacentric_factor", CoolProp::iacentric_factor)
        .value("irhomolar_reducing", CoolProp::irhomolar_reducing)
        .value("irhomolar_critical", CoolProp::irhomolar_critical)
        .value("iT_reducing", CoolProp::iT_reducing)
        .value("iT_critical", CoolProp::iT_critical)
        .value("irhomass_reducing", CoolProp::irhomass_reducing)
        .value("irhomass_critical", CoolProp::irhomass_critical)
        .value("iP_critical", CoolProp::iP_critical)
        .value("iP_reducing", CoolProp::iP_reducing)
        .value("iT_triple", CoolProp::iT_triple)
        .value("iP_triple", CoolProp::iP_triple)
        .value("iT_min", CoolProp::iT_min)
        .value("iT_max", CoolProp::iT_max)
        .value("iP_max", CoolProp::iP_max)
        .value("iP_min", CoolProp::iP_min)
        .value("idipole_moment", CoolProp::idipole_moment)
        .value("iT", CoolProp::iT)
        .value("iP", CoolProp::iP)
        .value("iQ", CoolProp::iQ)
        .value("iTau", CoolProp::iTau)
        .value("iDelta", CoolProp::iDelta)
        .value("iDmolar", CoolProp::iDmolar)
        .value("iHmolar", CoolProp::iHmolar)
        .value("iSmolar", CoolProp::iSmolar)
        .value("iCpmolar", CoolProp::iCpmolar)
        .value("iCp0molar", CoolProp::iCp0molar)
        .value("iCvmolar", CoolProp::iCvmolar)
        .value("iUmolar", CoolProp::iUmolar)
        .value("iGmolar", CoolProp::iGmolar)
        .value("iHelmholtzmolar", CoolProp::iHelmholtzmolar)
        .value("iHmolar_residual", CoolProp::iHmolar_residual)
        .value("iSmolar_residual", CoolProp::iSmolar_residual)
        .value("iGmolar_residual", CoolProp::iGmolar_residual)
        .value("iHmolar_idealgas", CoolProp::iHmolar_idealgas)
        .value("iSmolar_idealgas", CoolProp::iSmolar_idealgas)
        .value("iUmolar_idealgas", CoolProp::iUmolar_idealgas)
        .value("iDmass", CoolProp::iDmass)
        .value("iHmass", CoolProp::iHmass)
        .value("iSmass", CoolProp::iSmass)
        .value("iCpmass", CoolProp::iCpmass)
        .value("iCp0mass", CoolProp::iCp0mass)
        .value("iCvmass", CoolProp::iCvmass)
        .value("iUmass", CoolProp::iUmass)
        .value("iGmass", CoolProp::iGmass)
        .value("iHelmholtzmass", CoolProp::iHelmholtzmass)
        .value("iHmass_idealgas", CoolProp::iHmass_idealgas)
        .value("iSmass_idealgas", CoolProp::iSmass_idealgas)
        .value("iUmass_idealgas", CoolProp::iUmass_idealgas)
        .value("iviscosity", CoolProp::iviscosity)
        .value("iconductivity", CoolProp::iconductivity)
        .value("isurface_tension", CoolProp::isurface_tension)
        .value("iPrandtl", CoolProp::iPrandtl)
        .value("ispeed_sound", CoolProp::ispeed_sound)
        .value("iisothermal_compressibility", CoolProp::iisothermal_compressibility)
        .value("iisobaric_expansion_coefficient", CoolProp::iisobaric_expansion_coefficient)
        .value("iisentropic_expansion_coefficient", CoolProp::iisentropic_expansion_coefficient)
        .value("ifundamental_derivative_of_gas_dynamics", CoolProp::ifundamental_derivative_of_gas_dynamics)
        .value("ialphar", CoolProp::ialphar)
        .value("idalphar_dtau_constdelta", CoolProp::idalphar_dtau_constdelta)
        .value("idalphar_ddelta_consttau", CoolProp::idalphar_ddelta_consttau)
        .value("ialpha0", CoolProp::ialpha0)
        .value("idalpha0_dtau_constdelta", CoolProp::idalpha0_dtau_constdelta)
        .value("idalpha0_ddelta_consttau", CoolProp::idalpha0_ddelta_consttau)
        .value("id2alpha0_ddelta2_consttau", CoolProp::id2alpha0_ddelta2_consttau)
        .value("id3alpha0_ddelta3_consttau", CoolProp::id3alpha0_ddelta3_consttau)
        .value("iBvirial", CoolProp::iBvirial)
        .value("iCvirial", CoolProp::iCvirial)
        .value("idBvirial_dT", CoolProp::idBvirial_dT)
        .value("idCvirial_dT", CoolProp::idCvirial_dT)
        .value("iZ", CoolProp::iZ)
        .value("iPIP", CoolProp::iPIP)
        .value("ifraction_min", CoolProp::ifraction_min)
        .value("ifraction_max", CoolProp::ifraction_max)
        .value("iT_freeze", CoolProp::iT_freeze)
        .value("iGWP20", CoolProp::iGWP20)
        .value("iGWP100", CoolProp::iGWP100)
        .value("iGWP500", CoolProp::iGWP500)
        .value("iFH", CoolProp::iFH)
        .value("iHH", CoolProp::iHH)
        .value("iPH", CoolProp::iPH)
        .value("iODP", CoolProp::iODP)
        .value("iPhase", CoolProp::iPhase)
        .value("iundefined_parameter", CoolProp::iundefined_parameter)
        ;

    enum_<CoolProp::input_pairs>("input_pairs")
        .value("QT_INPUTS", CoolProp::QT_INPUTS)
        .value("PQ_INPUTS", CoolProp::PQ_INPUTS)
        .value("QSmolar_INPUTS", CoolProp::QSmolar_INPUTS)
        .value("QSmass_INPUTS", CoolProp::QSmass_INPUTS)
        .value("HmolarQ_INPUTS", CoolProp::HmolarQ_INPUTS)
        .value("HmassQ_INPUTS", CoolProp::HmassQ_INPUTS)
        .value("DmolarQ_INPUTS", CoolProp::DmolarQ_INPUTS)
        .value("DmassQ_INPUTS", CoolProp::DmassQ_INPUTS)
        .value("PT_INPUTS", CoolProp::PT_INPUTS)
        .value("DmassT_INPUTS", CoolProp::DmassT_INPUTS)
        .value("DmolarT_INPUTS", CoolProp::DmolarT_INPUTS)
        .value("HmolarT_INPUTS", CoolProp::HmolarT_INPUTS)
        .value("HmassT_INPUTS", CoolProp::HmassT_INPUTS)
        .value("SmolarT_INPUTS", CoolProp::SmolarT_INPUTS)
        .value("SmassT_INPUTS", CoolProp::SmassT_INPUTS)
        .value("TUmolar_INPUTS", CoolProp::TUmolar_INPUTS)
        .value("TUmass_INPUTS", CoolProp::TUmass_INPUTS)
        .value("DmassP_INPUTS", CoolProp::DmassP_INPUTS)
        .value("DmolarP_INPUTS", CoolProp::DmolarP_INPUTS)
        .value("HmassP_INPUTS", CoolProp::HmassP_INPUTS)
        .value("HmolarP_INPUTS", CoolProp::HmolarP_INPUTS)
        .value("PSmass_INPUTS", CoolProp::PSmass_INPUTS)
        .value("PSmolar_INPUTS", CoolProp::PSmolar_INPUTS)
        .value("PUmass_INPUTS", CoolProp::PUmass_INPUTS)
        .value("PUmolar_INPUTS", CoolProp::PUmolar_INPUTS)
        .value("HmassSmass_INPUTS", CoolProp::HmassSmass_INPUTS)
        .value("HmolarSmolar_INPUTS", CoolProp::HmolarSmolar_INPUTS)
        .value("SmassUmass_INPUTS", CoolProp::SmassUmass_INPUTS)
        .value("SmolarUmolar_INPUTS", CoolProp::SmolarUmolar_INPUTS)
        .value("DmassHmass_INPUTS", CoolProp::DmassHmass_INPUTS)
        .value("DmolarHmolar_INPUTS", CoolProp::DmolarHmolar_INPUTS)
        .value("DmassSmass_INPUTS", CoolProp::DmassSmass_INPUTS)
        .value("DmolarSmolar_INPUTS", CoolProp::DmolarSmolar_INPUTS)
        .value("DmassUmass_INPUTS", CoolProp::DmassUmass_INPUTS)
        .value("DmolarUmolar_INPUTS", CoolProp::DmolarUmolar_INPUTS)
        ;

    enum_<CoolProp::phases>("phases")
        .value("iphase_liquid", CoolProp::iphase_liquid)
        .value("iphase_supercritical", CoolProp::iphase_supercritical)
        .value("iphase_supercritical_gas", CoolProp::iphase_supercritical_gas)
        .value("iphase_supercritical_liquid", CoolProp::iphase_supercritical_liquid)
        .value("iphase_critical_point", CoolProp::iphase_critical_point)
        .value("iphase_gas", CoolProp::iphase_gas)
        .value("iphase_twophase", CoolProp::iphase_twophase)
        .value("iphase_unknown", CoolProp::iphase_unknown)
        .value("iphase_not_imposed", CoolProp::iphase_not_imposed)
        ;

    enum_<CoolProp::backend_families>("backend_families")
        .value("HEOS_BACKEND_FAMILY", CoolProp::HEOS_BACKEND_FAMILY)
        .value("REFPROP_BACKEND_FAMILY", CoolProp::REFPROP_BACKEND_FAMILY)
        .value("INCOMP_BACKEND_FAMILY", CoolProp::INCOMP_BACKEND_FAMILY)
        .value("IF97_BACKEND_FAMILY", CoolProp::IF97_BACKEND_FAMILY)
        .value("TREND_BACKEND_FAMILY", CoolProp::TREND_BACKEND_FAMILY)
        .value("TTSE_BACKEND_FAMILY", CoolProp::TTSE_BACKEND_FAMILY)
        .value("BICUBIC_BACKEND_FAMILY", CoolProp::BICUBIC_BACKEND_FAMILY)
        .value("SRK_BACKEND_FAMILY", CoolProp::SRK_BACKEND_FAMILY)
        .value("PR_BACKEND_FAMILY", CoolProp::PR_BACKEND_FAMILY)
        .value("VTPR_BACKEND_FAMILY", CoolProp::VTPR_BACKEND_FAMILY)
        .value("PCSAFT_BACKEND_FAMILY", CoolProp::PCSAFT_BACKEND_FAMILY)
        ;
}
// Binding code
EMSCRIPTEN_BINDINGS(humid_air_bindings) {
    function("HAPropsSI", &HumidAir::HAPropsSI);
}

CoolProp::AbstractState* factory(const std::string& backend, const std::string& fluid_names) {
    return CoolProp::AbstractState::factory(backend, strsplit(fluid_names, '&'));
}

// Binding code
EMSCRIPTEN_BINDINGS(abstract_state_bindings) {

    register_vector<double>("VectorDouble");
    register_vector<std::string>("VectorString");

    value_object<CoolProp::PhaseEnvelopeData>("CoolProp::PhaseEnvelopeData")
// Use X macros to auto-generate the variables;
// each will look something like: .field("T", &CoolProp::PhaseEnvelopeData::T);
#    define X(name) .field(#    name, &CoolProp::PhaseEnvelopeData::name)
      PHASE_ENVELOPE_VECTORS
#    undef X
      ;

    function("factory", &factory, allow_raw_pointers());

    class_<CoolProp::AbstractState>("AbstractState")
      .function("backend_name", &CoolProp::AbstractState::backend_name)
      .function("using_mole_fractions", &CoolProp::AbstractState::using_mole_fractions)
      .function("using_mass_fractions", &CoolProp::AbstractState::using_mass_fractions)
      .function("using_volu_fractions", &CoolProp::AbstractState::using_volu_fractions)
      .function("set_mass_fractions", &CoolProp::AbstractState::set_mass_fractions)
      .function("set_volu_fractions", &CoolProp::AbstractState::set_volu_fractions)
      .function("set_mole_fractions", &CoolProp::AbstractState::set_mole_fractions_double)
      .function("mole_fractions_liquid", &CoolProp::AbstractState::mole_fractions_liquid_double)
      .function("mole_fractions_vapor", &CoolProp::AbstractState::mole_fractions_vapor_double)
      .function("get_mole_fractions", &CoolProp::AbstractState::get_mole_fractions)
      .function("get_mass_fractions", &CoolProp::AbstractState::get_mass_fractions)

      .function("update", &CoolProp::AbstractState::update)

      .function("T", &CoolProp::AbstractState::T)
      .function("rhomolar", &CoolProp::AbstractState::rhomolar)
      .function("rhomass", &CoolProp::AbstractState::rhomass)
      .function("p", &CoolProp::AbstractState::p)
      .function("Q", &CoolProp::AbstractState::Q)
      .function("tau", &CoolProp::AbstractState::tau)
      .function("delta", &CoolProp::AbstractState::delta)
      .function("molar_mass", &CoolProp::AbstractState::molar_mass)
      .function("acentric_factor", &CoolProp::AbstractState::acentric_factor)
      .function("gas_constant", &CoolProp::AbstractState::gas_constant)
      .function("Bvirial", &CoolProp::AbstractState::Bvirial)
      .function("Cvirial", &CoolProp::AbstractState::Cvirial)
      .function("compressibility_factor", &CoolProp::AbstractState::compressibility_factor)
      .function("hmolar", &CoolProp::AbstractState::hmolar)
      .function("hmass", &CoolProp::AbstractState::hmass)
      .function("smolar", &CoolProp::AbstractState::smolar)
      .function("smass", &CoolProp::AbstractState::smass)
      .function("umolar", &CoolProp::AbstractState::umolar)
      .function("umass", &CoolProp::AbstractState::umass)
      .function("cpmolar", &CoolProp::AbstractState::cpmolar)
      .function("cpmass", &CoolProp::AbstractState::cpmass)
      .function("cvmolar", &CoolProp::AbstractState::cvmolar)
      .function("cvmass", &CoolProp::AbstractState::cvmass)
      .function("gibbsmolar", &CoolProp::AbstractState::gibbsmolar)
      .function("gibbsmass", &CoolProp::AbstractState::gibbsmass)
      .function("helmholtzmolar", &CoolProp::AbstractState::helmholtzmolar)
      .function("helmholtzmass", &CoolProp::AbstractState::helmholtzmass)

      .function("speed_sound", &CoolProp::AbstractState::speed_sound)
      .function("isothermal_compressibility", &CoolProp::AbstractState::isothermal_compressibility)
      .function("isobaric_expansion_coefficient", &CoolProp::AbstractState::isobaric_expansion_coefficient)
      .function("isentropic_expansion_coefficient", &CoolProp::AbstractState::isentropic_expansion_coefficient)
      .function("viscosity", &CoolProp::AbstractState::viscosity)
      .function("conductivity", &CoolProp::AbstractState::conductivity)
      .function("surface_tension", &CoolProp::AbstractState::surface_tension)
      .function("Prandtl", &CoolProp::AbstractState::Prandtl)

      .function("keyed_output", &CoolProp::AbstractState::keyed_output)
      .function("trivial_keyed_output", &CoolProp::AbstractState::trivial_keyed_output)
      .function("saturated_liquid_keyed_output", &CoolProp::AbstractState::saturated_liquid_keyed_output)
      .function("saturated_vapor_keyed_output", &CoolProp::AbstractState::saturated_vapor_keyed_output)

      .function("first_partial_deriv", &CoolProp::AbstractState::first_partial_deriv)
      .function("second_partial_deriv", &CoolProp::AbstractState::second_partial_deriv)
      .function("first_saturation_deriv", &CoolProp::AbstractState::first_saturation_deriv)
      .function("second_saturation_deriv", &CoolProp::AbstractState::second_saturation_deriv)
      .function("first_two_phase_deriv", &CoolProp::AbstractState::first_two_phase_deriv)
      .function("second_two_phase_deriv", &CoolProp::AbstractState::second_two_phase_deriv)

      .function("build_phase_envelope", &CoolProp::AbstractState::build_phase_envelope)
      .function("get_phase_envelope_data", &CoolProp::AbstractState::get_phase_envelope_data)

      .function("melting_line", &CoolProp::AbstractState::melting_line)
      .function("saturation_ancillary", &CoolProp::AbstractState::saturation_ancillary)

      .function("T_critical", &CoolProp::AbstractState::T_critical)
      .function("p_critical", &CoolProp::AbstractState::p_critical)
      .function("rhomolar_critical", &CoolProp::AbstractState::rhomolar_critical)
      .function("rhomass_critical", &CoolProp::AbstractState::rhomass_critical)

      .function("T_reducing", &CoolProp::AbstractState::T_reducing)
      .function("rhomolar_reducing", &CoolProp::AbstractState::rhomolar_reducing)
      .function("rhomass_reducing", &CoolProp::AbstractState::rhomass_reducing)

      .function("p_triple", &CoolProp::AbstractState::p_triple)
      .function("Ttriple", &CoolProp::AbstractState::Ttriple)
      .function("Tmin", &CoolProp::AbstractState::Tmin)
      .function("Tmax", &CoolProp::AbstractState::Tmax)
      .function("pmax", &CoolProp::AbstractState::pmax)

      .function("dipole_moment", &CoolProp::AbstractState::dipole_moment)
      ;
}

#endif
