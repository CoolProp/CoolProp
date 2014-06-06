/*
 * Incompressible.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#include "Incompressible.h"


namespace CoolProp {

std::vector<std::vector<double> > Incompressible::changeAxis(const std::vector<double> &input){
	std::vector<std::vector<double> > output;
	std::vector<double> tmp;
	size_t sizeX = input.size();
	for (size_t i = 0; i < sizeX; ++i){
		tmp.clear();
		tmp.push_back(input[i]);
		output.push_back(tmp);
	}
	return output;
}

// /* All functions need T and p as input. Might not
//  * be necessary, but gives a clearer structure.
//  */
// /// Density as a function of temperature, pressure and composition.
// double Incompressible::rho(double T_K, double p) {
// 	return poly.polyval(cRho, getxInput(x), getTInput(T_K));
// }
// /// Heat capacities as a function of temperature, pressure and composition.
// double Incompressible::c(double T_K, double p) {
// 	return poly.polyval(cHeat, getxInput(x), getTInput(T_K));
// }
// /// Enthalpy as a function of temperature, pressure and composition.
// double Incompressible::h(double T_K, double p) {
// 	return h_u(T_K, p);
// }
// /// Entropy as a function of temperature, pressure and composition.
// double Incompressible::s(double T_K, double p) {
// 	return poly.polyfracintcentral(cHeat, getxInput(x), T_K, Tbase)
// 			- poly.polyfracintcentral(cHeat, getxInput(x), Tref, Tbase);
// }
// /// Viscosity as a function of temperature, pressure and composition.
// double Incompressible::visc(double T_K, double p) {
// 	return expo.expval(cVisc, getxInput(x), getTInput(T_K), 2) / 1e3;
// }
// /// Thermal conductivity as a function of temperature, pressure and composition.
// double Incompressible::cond(double T_K, double p) {
// 	return poly.polyval(cCond, getxInput(x), getTInput(T_K));
// }
// /// Internal energy as a function of temperature, pressure and composition.
// double Incompressible::u(double T_K, double p) {
// 	return poly.polyint(cHeat, getxInput(x), getTInput(T_K))
// 			- poly.polyint(cHeat, getxInput(x), getTInput(Tref));
// }

/// Saturation pressure as a function of temperature and composition.
double Incompressible::psat(double T_K          ){throw NotImplementedError("Psat is not available");};
/// Freezing temperature as a function of pressure and composition.
double Incompressible::Tfreeze(         double p){throw NotImplementedError("Tfreeze is not available");};


/*
 * Some more functions to provide a single implementation
 * of important routines.
 * We start with the check functions that can validate input
 * in terms of pressure p, temperature T and composition x.
 */

/// Check validity of temperature input.
/** Compares the given temperature T to the result of a
 *  freezing point calculation. This is not necessarily
 *  defined for all fluids, default values do not
 *  cause errors. */
bool Incompressible::checkT(double T_K, double p){
	if( Tmin < 0. ) {
		throw ValueError("Please specify the minimum temperature.");
	} else if( Tmax < 0.) {
		throw ValueError("Please specify the maximum temperature.");
	} else if ( (Tmin>T_K) || (T_K>Tmax) ) {
		throw ValueError(format("Your temperature %f is not between %f and %f.",T_K,Tmin,Tmax));
	} else if (T_K < Tfreeze(p)) {
		throw ValueError(format("Your temperature %f is below the freezing point of %f.",T_K,Tfreeze(p)));
	} else {
		return true;
	}
	return false;
}

/// Check validity of pressure input.
/** Compares the given pressure p to the saturation pressure at
 *  temperature T and throws and exception if p is lower than
 *  the saturation conditions.
 *  The default value for psat is -1 yielding true if psat
 *  is not redefined in the subclass.
 *  */
bool Incompressible::checkP(double T_K, double p) {
	double ps = psat(T_K);
	if (p<ps) {
		throw ValueError(format("Equations are valid for solution phase only: %f < %f. ",p,ps));
	} else {
		return true;
	}
}

/// Check validity of composition input.
/** Compares the given composition x to a stored minimum and
 *  maximum value. Enforces the redefinition of xmin and
 *  xmax since the default values cause an error. */
bool Incompressible::checkX(double x){
	if( xmin < 0. ) {
		throw ValueError("Please specify the minimum concentration.");
	} else if( xmax < 0.) {
		throw ValueError("Please specify the maximum concentration.");
	} else if ( (xmin>x) || (x>xmax) ) {
		throw ValueError(format("Your composition %f is not between %f and %f.",x,xmin,xmax));
	} else {
		return true;
	}
	return false;
}

/// Check validity of temperature, pressure and composition input.
bool Incompressible::checkTPX(double T, double p, double x) {
	return (checkT(T,p) && checkP(T,p) && checkX(x));
}

} /* namespace CoolProp */
