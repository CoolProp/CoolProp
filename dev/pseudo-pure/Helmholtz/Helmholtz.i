%module helmholtz

//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
//// *************************** EXCEPTION HANDLING ****************************
%include exception.i

// A generic exception handler.  Any exceptions thrown in C++ will be caught here
%exception {
	try {
		$action
	}
    catch(std::exception &e) {
		SWIG_exception(SWIG_RuntimeError,e.what());
	}
    catch(...) {
		SWIG_exception(SWIG_RuntimeError,"Unknown exception");
	}
}

// This allows for the use of STL vectors
%include "std_vector.i"
// This allows for the use of STL strings
%include "std_string.i"

namespace std {
   %template(vectord) vector<double>;
};

%ignore phir_power::phir_power();
%ignore phir_power::phir_power(std::vector<double>,std::vector< double>,std::vector< double>,int,int);
%ignore phir_power::phir_power(const double [],const double [], const double [],int,int,int);
%ignore phir_power::phir_power(double [],double [],double [],int,int,int);
%ignore phir_power::phir_power(double const [],double const [],double const [],double const [],int,int,int);
%ignore phir_power::phir_power(double [],double [],double [],double [],int,int,int);
%ignore phi0_Planck_Einstein::phi0_Planck_Einstein(double [],double [],int,int,int);
%ignore phi0_Planck_Einstein::phi0_Planck_Einstein(double const [],double const [],int,int,int);
%ignore phi0_Planck_Einstein::phi0_Planck_Einstein(double,double);

// This stuff will get included verbatim in CoolProp_wrap.cpp
%{
#include "../../../CoolProp/Helmholtz.h"
#include "../../../CoolProp/CoolPropTools.h"
%}

// This is where the parsing actually happens
%include "../../../CoolProp/Helmholtz.h"