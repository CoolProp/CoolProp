%module CoolProp

// Hack the AbstractState to return the mole fractions as a vector<double> which SWIG can wrap 
%ignore CoolProp::AbstractState::mole_fractions_liquid();
%ignore CoolProp::AbstractState::mole_fractions_vapor();
%rename (mole_fractions_liquid) CoolProp::AbstractState::mole_fractions_liquid_double();
%rename (mole_fractions_vapor) CoolProp::AbstractState::mole_fractions_vapor_double();

%ignore CoolProp::AbstractState::set_mole_fractions(const std::vector<CoolPropDbl> &);
%ignore CoolProp::AbstractState::set_mass_fractions(const std::vector<CoolPropDbl> &);
%ignore CoolProp::AbstractState::set_volu_fractions(const std::vector<CoolPropDbl> &);

%ignore CoolProp::set_config_json(rapidjson::Document &);
%ignore CoolProp::get_config_as_json(rapidjson::Document &);

%include "std_string.i" // This %include allows the use of std::string natively
%include "std_vector.i" // This allows for the use of STL vectors natively(ish)
%include "exception.i" //

// Instantiate templates used
%template(DoubleVector) std::vector<double>;
%template(VectorOfDoubleVector) std::vector<std::vector<double> >;
%template(StringVector) std::vector<std::string>;
%template(VectorOfStringVector) std::vector< std::vector<std::string> >;

%apply double { CoolPropDbl }; 

%exception { 
    try {
        $action
    }  catch (std::exception &e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
    }
}

#ifdef SWIGSCILAB
// Shorten some names to make Scilab a bit happier
%rename(SS) CoolProp::SimpleState;
%rename(PED) CoolProp::PhaseEnvelopeData;
%rename(GE) CoolProp::GuessesStructure;
#endif

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