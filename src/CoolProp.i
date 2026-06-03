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

// fast_evaluate is a C++-only raw-pointer batch API; SWIG's default
// typemaps for `const double*` / `int*` produce wrappers that don't
// compile on at least the R binding (int -> SEXP).  Scripted-language
// callers should keep using update()/keyed_output() or the Cython
// AbstractState.fast_evaluate wrapper in the Python binding.
%ignore CoolProp::AbstractState::fast_evaluate;

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
#include "CoolProp/DataStructures.h"
#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolProp.h"
#include "CoolProp/fluids/PhaseEnvelope.h"
#define SWIG
#include "CoolProp/Configuration.h"
#undef SWIG
#include "CoolProp/HumidAirProp.h"
#include "Backends/Helmholtz/MixtureParameters.h"
%}

%include "CoolProp/DataStructures.h"
%include "CoolProp/AbstractState.h"
%include "CoolProp/CoolProp.h"
%include "CoolProp/fluids/PhaseEnvelope.h"
%include "CoolProp/Configuration.h"
%include "CoolProp/HumidAirProp.h"
namespace CoolProp {
void apply_simple_mixing_rule(const std::string& identifier1, const std::string& identifier2, const std::string& rule);
}
