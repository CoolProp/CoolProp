
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

    enum_<CoolProp::input_pairs>("input_pairs").value("PT_INPUTS", CoolProp::PT_INPUTS);
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
      .function("gas_constant", &CoolProp::AbstractState::gas_constant)
      .function("update", &CoolProp::AbstractState::update)
      .function("p", &CoolProp::AbstractState::p)
      .function("rhomass", &CoolProp::AbstractState::rhomass)
      .function("viscosity", &CoolProp::AbstractState::viscosity)
      .function("set_mole_fractions", &CoolProp::AbstractState::set_mole_fractions_double)
      .function("build_phase_envelope", &CoolProp::AbstractState::build_phase_envelope)
      .function("get_phase_envelope_data", &CoolProp::AbstractState::get_phase_envelope_data);
}

#endif