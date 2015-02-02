%module CoolProp

%ignore CoolProp::AbstractState::set_mole_fractions(const std::vector< long double> &);
%ignore CoolProp::AbstractState::set_mass_fractions(const std::vector< long double> &);
%ignore CoolProp::set_config_json(rapidjson::Document &);
%ignore CoolProp::get_config_as_json(rapidjson::Document &);

%include "std_string.i" // This %include allows the use of std::string natively
%include "std_vector.i" // This allows for the use of STL vectors natively(ish)
%include "exception.i" //

// Instantiate templates used by example
namespace std {
   %template(LongDoubleVector) vector<long double>;
   %template(DoubleVector) vector<double>;
}

%exception { 
    try {
        $action
    }  catch (std::exception &e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
    }
}

// 
%rename(SS) CoolProp::SimpleState;
%rename(PED) CoolProp::PhaseEnvelopeData;

// This stuff will get included verbatim in CoolProp_wrap
%{
#include "DataStructures.h"
#include "AbstractState.h"
#include "CoolProp.h"
#include "PhaseEnvelope.h"
#define SWIG
#include "Configuration.h"
#undef SWIG
#include "HumidAirProp.h"
%}

%include "DataStructures.h"
%include "AbstractState.h"
%include "CoolProp.h"
%include "PhaseEnvelope.h"
%include "Configuration.h"
%include "HumidAirProp.h"