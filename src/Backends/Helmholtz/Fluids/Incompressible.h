/*
 * Incompressible.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INCOMPRESSIBLE_H_
#define INCOMPRESSIBLE_H_

#include <Eigen/Core>
#include "PolyMath.h"
#include "MatrixMath.h"

namespace CoolProp {

class Incompressible{

protected:
	std::string name;
	std::string description;
	std::string reference;

	double Tmin, Tmax;
	double xmin, xmax;

	double TminPsat;
	double xref, Tref;
	double xbase, Tbase;

	Eigen::MatrixXd cRho;
	Eigen::MatrixXd cHeat;
	Eigen::MatrixXd cVisc;
	Eigen::MatrixXd cCond;
	Eigen::MatrixXd cPsat;
	Eigen::MatrixXd cTfreeze;
	Eigen::MatrixXd cV2M;

	Polynomial2DFrac poly;

public:
	Incompressible();
	virtual ~Incompressible();

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
	double getxref() const {return xref;}
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
	void setTref(double Tref) {this->Tref = Tref;}
	void setxref(double xref) {this->xref = xref;}
	void setTbase(double Tbase) {this->Tbase = Tbase;}
	void setxbase(double xbase) {this->xbase = xbase;}

	// Setters for the coefficients
	void setcRho(Eigen::MatrixXd cRho){this->cRho = cRho;}
	void setcHeat(Eigen::MatrixXd cHeat){this->cHeat = cHeat;}
	void setcVisc(Eigen::MatrixXd cVisc){this->cVisc = cVisc;}
	void setcCond(Eigen::MatrixXd cCond){this->cCond = cCond;}
	void setcPsat(Eigen::MatrixXd cPsat){this->cPsat = cPsat;}
	void setcTfreeze(Eigen::MatrixXd cTfreeze){this->cTfreeze = cTfreeze;}
	void setcV2M(Eigen::MatrixXd cV2M){this->cV2M = cV2M;}

	double getTInput(double curTValue){return curTValue-Tbase;}
	double getxInput(double curxValue){return (curxValue-xbase)*100.0;}

	/* All functions need T and p as input. Might not
	 * be necessary, but gives a clearer structure.
	 */
	/// Density as a function of temperature, pressure and composition.
	virtual double rho (double T_K, double p, double x);
	/// Heat capacities as a function of temperature, pressure and composition.
	virtual double c   (double T_K, double p, double x);
	virtual double cp  (double T_K, double p, double x){return c(T_K,p,x);};
	virtual double cv  (double T_K, double p, double x){return c(T_K,p,x);};
	/// Entropy as a function of temperature, pressure and composition.
	virtual double s   (double T_K, double p, double x);
	/// Internal energy as a function of temperature, pressure and composition.
	virtual double u   (double T_K, double p, double x);
	/// Enthalpy as a function of temperature, pressure and composition.
	virtual double h   (double T_K, double p, double x);
	/// Viscosity as a function of temperature, pressure and composition.
	virtual double visc(double T_K, double p, double x);
	/// Thermal conductivity as a function of temperature, pressure and composition.
	virtual double cond(double T_K, double p, double x);
	/// Saturation pressure as a function of temperature and composition.
	virtual double psat(double T_K, double x          );
	/// Freezing temperature as a function of pressure and composition.
	virtual double Tfreeze(         double p, double x);
	/// Conversion from volume-based to mass-based composition.
	virtual double V2M(                       double x);


protected:
	/* Define internal energy and enthalpy as functions of the
	 * other properties to provide data in case there are no
	 * coefficients.
	 */

	/// Enthalpy from u, p and rho.
	/** Calculate enthalpy as a function of temperature and
	 *  pressure employing functions for internal energy and
	 *  density. Provides consistent formulations. */
	double h_u(double T_K, double p, double x) {
		return u(T_K,p,x)+p/rho(T_K,p,x);
	};

	/// Internal energy from h, p and rho.
	/** Calculate internal energy as a function of temperature
	 *  and pressure employing functions for enthalpy and
	 *  density. Provides consistent formulations. */
	double u_h(double T_K, double p, double x) {
		return h(T_K,p,x)-p/rho(T_K,p,x);
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
	bool checkT(double T_K, double p, double x);

	/// Check validity of pressure input.
	/** Compares the given pressure p to the saturation pressure at
	 *  temperature T and throws and exception if p is lower than
	 *  the saturation conditions.
	 *  The default value for psat is -1 yielding true if psat
	 *  is not redefined in the subclass.
	 *  */
	bool checkP(double T_K, double p, double x);

	/// Check validity of composition input.
	/** Compares the given composition x to a stored minimum and
	 *  maximum value. Enforces the redefinition of xmin and
	 *  xmax since the default values cause an error. */
	bool checkX(double x);

	/// Check validity of temperature, pressure and composition input.
	bool checkTPX(double T, double p, double x);


};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLE_H_ */
