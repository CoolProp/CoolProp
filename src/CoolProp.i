%module CoolProp

%ignore CoolProp::AbstractState::set_mole_fractions(const std::vector< long double> &);
%ignore CoolProp::AbstractState::set_mass_fractions(const std::vector< long double> &);

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
%rename(CoolProp_SimpleState) CoolProp::SimpleState;

#ifdef SWIGOCTAVE
    %ignore DerivTerms(char *Term, double T, double rho, Fluid * pFluid, bool SinglePhase, bool TwoPhase);
#endif

// This stuff will get included verbatim in CoolProp_wrap
%{
#include "AbstractState.h"
#include "DataStructures.h"
#include "CoolProp.h"
%}

%include "AbstractState.h"
%include "DataStructures.h"
%include "CoolProp.h"