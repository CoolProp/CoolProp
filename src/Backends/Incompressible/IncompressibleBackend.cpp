#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#include <string>
//#include "CoolProp.h"

#include "IncompressibleBackend.h"
#include "../Helmholtz/Fluids/FluidLibrary.h"
#include "Solvers.h"
#include "MatrixMath.h"

namespace CoolProp {

IncompressibleBackend::IncompressibleBackend(const std::string &fluid_name) {
	throw CoolProp::NotImplementedError("Mixture-style constructor is not implemented yet for incompressible fluids");
}

IncompressibleBackend::IncompressibleBackend(const std::vector<std::string> &component_names) {
	throw CoolProp::NotImplementedError("Mixture-style constructor is not implemented yet for incompressible fluids");
}

bool IncompressibleBackend::using_mole_fractions(){return false;}


void IncompressibleBackend::update(long input_pair, double value1, double value2) {
	switch (input_pair) {
	case PT_INPUTS: {

		break;
	}
//	case DmassP_INPUTS: {
//
//	}
//		break;
//	}
//	case HmassP_INPUTS: {
//		// Call again, but this time with molar units
//		// H: [J/kg] * [kg/mol] -> [J/mol]
//		update(HmolarP_INPUTS, value1 * (double) _molar_mass, value2);
//		return;
//	}
//	case PUmass_INPUTS: {
//		// Call again, but this time with molar units
//		// U: [J/kg] * [kg/mol] -> [J/mol]
//		update(PUmolar_INPUTS, value1, value2 * (double) _molar_mass);
//		return;
//	}
//	case PSmass_INPUTS: {
//		// Call again, but this time with molar units
//		// U: [J/kg] * [kg/mol] -> [J/mol]
//		update(PUmolar_INPUTS, value1, value2 * (double) _molar_mass);
//		return;
//	}
	default: {
		throw ValueError(
				format("This pair of inputs [%s] is not yet supported",
						get_input_pair_short_desc(input_pair).c_str()));
	}
	}
}

/// Set the mole fractions
/**
@param mole_fractions The vector of mole fractions of the components
*/
void IncompressibleBackend::set_mole_fractions(const std::vector<long double> &mole_fractions) {
	throw CoolProp::NotImplementedError("Cannot set mole fractions for incompressible fluid");
}

/// Set the mass fractions
/**
@param mass_fractions The vector of mass fractions of the components
*/
void IncompressibleBackend::set_mass_fractions(const std::vector<long double> &mass_fractions) {
	this->mass_fractions = mass_fractions;
}

/// Check if the mole fractions have been set, etc.
void IncompressibleBackend::check_status() {
	throw CoolProp::NotImplementedError("Cannot check status for incompressible fluid");
}

/// Get the viscosity [Pa-s]
long double IncompressibleBackend::calc_viscosity(void){
	throw NotImplementedError();
}
/// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
long double IncompressibleBackend::calc_conductivity(void){
	throw NotImplementedError();
}

}
