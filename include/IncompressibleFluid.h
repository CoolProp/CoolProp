/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INCOMPRESSIBLEFLUID_H_
#define INCOMPRESSIBLEFLUID_H_

#include "DataStructures.h"
#include "Helmholtz.h"
#include "Solvers.h"

#include <numeric>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>

#include <Eigen/Core>
#include "PolyMath.h"
#include "MatrixMath.h"

namespace CoolProp {

struct IncompressibleData {
	int type;
	enum IncompressibleTypeEnum {
		INCOMPRESSIBLE_POLYNOMIAL,
		INCOMPRESSIBLE_EXPONENTIAL,
		INCOMPRESSIBLE_EXPPOLYNOMIAL,
		INCOMPRESSIBLE_NOT_SET
	};
	Eigen::MatrixXd coeffs; //TODO: Can we store the Eigen::Matrix objects more efficiently?
	//std::vector<std::vector<double> > coeffs;
	IncompressibleData() {
		type = INCOMPRESSIBLE_NOT_SET;
	};
};

/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
class IncompressibleFluid{

protected:
	std::string name;
	std::string description;
	std::string reference;

	double Tmin, Tmax;
	double xmin, xmax;

	double TminPsat;
	double xref, Tref, pref;
	double href, sref;
	double uref, rhoref;
	double xbase, Tbase;

	IncompressibleData density;
	IncompressibleData specific_heat;
	IncompressibleData viscosity;
	IncompressibleData conductivity;
	IncompressibleData p_sat;
	IncompressibleData T_freeze;
	IncompressibleData volToMass;
	IncompressibleData massToMole;

	Polynomial2DFrac poly;

	// Forward declaration of the some internal functions
	//double h_u(double T, double p, double x);
	//double u_h(double T, double p, double x);

public:
	IncompressibleFluid(){};
	virtual ~IncompressibleFluid(){};

	std::string getName() const {return name;}
	std::string get_name() const {return getName();}// For backwards-compatibility.
	std::string getDescription() const {return description;}
	std::string getReference() const {return reference;}

	double getTmax() const {return Tmax;}
	double getTmin() const {return Tmin;}
	double getxmax() const {return xmax;}
	double getxmin() const {return xmin;}
	double getTminPsat() const {return TminPsat;}
	double getTref() const {return Tref;}
	double getpref() const {return pref;}
	double getxref() const {return xref;}
	double gethref() const {return href;}
	double getsref() const {return sref;}
	double getTbase() const {return Tbase;}
	double getxbase() const {return xbase;}

	void setName(std::string name) {this->name = name;}
	void setDescription(std::string description) {this->description = description;}
	void setReference(std::string reference) {this->reference = reference;}
	void setTmax(double Tmax) {this->Tmax = Tmax;}
	void setTmin(double Tmin) {this->Tmin = Tmin;}
	void setxmax(double xmax) {this->xmax = xmax;}
	void setxmin(double xmin) {this->xmin = xmin;}
	void setTminPsat(double TminPsat) {this->TminPsat = TminPsat;}
	//void setTref(double Tref) {this->Tref = Tref;}
	//void setpref(double pref) {this->pref = pref;}
	//void setxref(double xref) {this->xref = xref;}
	void set_reference_state(double T0, double p0, double x0, double h0, double s0);
	void setTbase(double Tbase) {this->Tbase = Tbase;}
	void setxbase(double xbase) {this->xbase = xbase;}

	/// Setters for the coefficients
	void setDensity(IncompressibleData density){this->density = density;}
	void setSpecificHeat(IncompressibleData specific_heat){this->specific_heat = specific_heat;}
	void setViscosity(IncompressibleData viscosity){this->viscosity = viscosity;}
	void setConductivity(IncompressibleData conductivity){this->conductivity = conductivity;}
	void setPsat(IncompressibleData p_sat){this->p_sat = p_sat;}
	void setTfreeze(IncompressibleData T_freeze){this->T_freeze = T_freeze;}
	void setVolToMass(IncompressibleData volToMass){this->volToMass = volToMass;}
	void setMassToMole(IncompressibleData massToMole){this->massToMole = massToMole;}

	/// A function to check coefficients and equation types.
	void validate();


protected:
	/// Base function that handles the custom data type, just a place holder to show the structure.
	double baseFunction(IncompressibleData data, double x_in, double y_in);

public:
	/* All functions need T and p as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
	/// Density as a function of temperature, pressure and composition.
	double rho (double T, double p, double x=0.0){return baseFunction(density, T, x);};
	/// Heat capacities as a function of temperature, pressure and composition.
	double c   (double T, double p, double x=0.0){return baseFunction(specific_heat, T, x);};
	double cp  (double T, double p, double x=0.0){return c(T,p,x);};
	double cv  (double T, double p, double x=0.0){return c(T,p,x);};
	/// Entropy as a function of temperature, pressure and composition.
	double s   (double T, double p, double x);
	/// Internal energy as a function of temperature, pressure and composition.
	double u   (double T, double p, double x);
	/// Enthalpy as a function of temperature, pressure and composition.
	double h   (double T, double p, double x=0.0){return h_u(T,p,x);};
	/// Viscosity as a function of temperature, pressure and composition.
	double visc(double T, double p, double x=0.0){return baseFunction(viscosity, T, x);};
	/// Thermal conductivity as a function of temperature, pressure and composition.
	double cond(double T, double p, double x=0.0){return baseFunction(conductivity, T, x);};
	/// Saturation pressure as a function of temperature and composition.
	double psat(double T,           double x=0.0){return baseFunction(p_sat, T, x);};
	/// Freezing temperature as a function of pressure and composition.
	double Tfreeze(       double p, double x);
	/// Conversion from volume-based to mass-based composition.
	double V2M (double T,           double y);
	/// Conversion from mass-based to mole-based composition.
	double M2M (double T,           double x);


protected:

	/* Define internal energy and enthalpy as functions of the
	 * other properties to provide data in case there are no
	 * coefficients.
	 */
	/// Enthalpy from u, p and rho.
	/** Calculate enthalpy as a function of temperature and
	 *  pressure employing functions for internal energy and
	 *  density. Provides consistent formulations. */
	double h_u(double T, double p, double x) {
		return u(T,p,x)+p/rho(T,p,x)-href;
	};

	/// Internal energy from h, p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Provides consistent formulations. */
	double u_h(double T, double p, double x) {
		return h(T,p,x)-p/rho(T,p,x)+href;
	};


	/*
	 * Some more functions to provide a single implementation
	 * of important routines.
	 * We start with the check functions that can validate input
	 * in terms of pressure p, temperature T and composition x.
	 */
	/// Check validity of temperature input.
	/** Compares the given temperature T to the result of a
	 *  freezing point calculation. This is not necessarily
	 *  defined for all fluids, default values do not cause errors. */
	bool checkT(double T, double p, double x);

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T, double p, double x);

	/// Check validity of composition input.
	/** Compares the given composition x to a stored minimum and
	 *  maximum value. Enforces the redefinition of xmin and
	 *  xmax since the default values cause an error. */
	bool checkX(double x);

	/// Check validity of temperature, pressure and composition input.
	bool checkTPX(double T, double p, double x);
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEFLUID_H_ */
