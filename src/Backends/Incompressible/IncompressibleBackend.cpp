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
#include "IncompressibleFluid.h"
#include "IncompressibleLibrary.h"
#include "Solvers.h"
#include "MatrixMath.h"

namespace CoolProp {

IncompressibleBackend::IncompressibleBackend(const std::string &fluid_name) {
    fluid = &get_incompressible_fluid(fluid_name);
}

IncompressibleBackend::IncompressibleBackend(const std::vector<std::string> &component_names) {
	throw NotImplementedError("Mixture-style constructor is not implemented yet for incompressible fluids");
}

void IncompressibleBackend::update(long input_pair, double value1, double value2) {
    //if (mass_fractions.empty()){
    //    throw ValueError("mass fractions have not been set");
    //}
	switch (input_pair) 
	{
        case PT_INPUTS: {
            _p = value1; _T = value2; 
            break;
        }
        case DmassP_INPUTS: {
            break;
        }
        case PUmass_INPUTS: {
            break;
        }
        case PSmass_INPUTS: {
            break;
        }
        case HmassP_INPUTS: {
            break;
        }
        default: {
            throw ValueError(
                    format("This pair of inputs [%s] is not yet supported",
                            get_input_pair_short_desc(input_pair).c_str()));
        }
	}
	if (_p < 0){ throw ValueError("p is less than zero");}
    if (!ValidNumber(_p)){ throw ValueError("p is not a valid number");}
    if (_T < 0){ throw ValueError("T is less than zero");}
    if (!ValidNumber(_T)){ throw ValueError("T is not a valid number");}
}

/// Set the mole fractions
/**
@param mole_fractions The vector of mole fractions of the components
*/
void IncompressibleBackend::set_mole_fractions(const std::vector<long double> &mole_fractions) {
	throw NotImplementedError("Cannot set mole fractions for incompressible fluid");
}

/// Set the mass fractions
/**
@param mass_fractions The vector of mass fractions of the components
*/
void IncompressibleBackend::set_mass_fractions(const std::vector<long double> &mass_fractions) {
	if (mass_fractions.size()!=1) throw ValueError(format("The incompressible backend only supports one entry in the mass fraction vector and not %d.",mass_fractions.size()));
	this->mass_fractions = mass_fractions;
}

/// Check if the mole fractions have been set, etc.
void IncompressibleBackend::check_status() {
	throw NotImplementedError("Cannot check status for incompressible fluid");
}

///// Calculate T given pressure and density
///**
//@param rhomass The mass density in kg/m^3
//@param p The pressure in Pa
//@returns T The temperature in K
//*/
//long double IncompressibleBackend::DmassP_flash(long double rhomass, long double p){
//
//}
///// Calculate T given pressure and enthalpy
///**
//@param hmass The mass enthalpy in J/kg
//@param p The pressure in Pa
//@returns T The temperature in K
//*/
//long double IncompressibleBackend::HmassP_flash(long double hmass, long double p);
///// Calculate T given pressure and entropy
///**
//@param smass The mass entropy in J/kg/K
//@param p The pressure in Pa
//@returns T The temperature in K
//*/
//long double IncompressibleBackend::PSmass_flash(long double p, long double smass);
//
///// Calculate T given pressure and internal energy
///**
//@param umass The mass internal energy in J/kg
//@param p The pressure in Pa
//@returns T The temperature in K
//*/
//long double IncompressibleBackend::PUmass_flash(long double p, long double umass);
//









}
