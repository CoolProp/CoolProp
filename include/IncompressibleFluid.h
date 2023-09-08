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
#include <cassert>
#include <iterator>

#include <Eigen/Core>
#include "PolyMath.h"
#include "MatrixMath.h"

namespace CoolProp {

struct IncompressibleData
{
    enum IncompressibleTypeEnum
    {
        INCOMPRESSIBLE_NOT_SET,
        INCOMPRESSIBLE_POLYNOMIAL,
        INCOMPRESSIBLE_EXPPOLYNOMIAL,
        INCOMPRESSIBLE_EXPONENTIAL,
        INCOMPRESSIBLE_LOGEXPONENTIAL,
        INCOMPRESSIBLE_POLYOFFSET
    };
    IncompressibleTypeEnum type;
    Eigen::MatrixXd coeffs;  //TODO: Can we store the Eigen::Matrix objects more efficiently?
    //std::vector<std::vector<double> > coeffs;
    IncompressibleData() {
        type = INCOMPRESSIBLE_NOT_SET;
    };
};

#if !defined(NO_FMTLIB) && FMT_VERSION >= 90000
inline int format_as(IncompressibleData::IncompressibleTypeEnum type) {
    return fmt::underlying(type);
}
#endif

/// A property provider for incompressible solutions and pure fluids
/**
This fluid instance is populated using an entry from a JSON file
*/
class IncompressibleFluid
{

   protected:
    bool strict;

    std::string name;
    std::string description;
    std::string reference;

    double Tmin, Tmax;
    double xmin, xmax;
    composition_types xid;

    double TminPsat;
    double xbase, Tbase;

    /// These are the objects that hold the coefficients
    /** Note that all polynomials require a 2-dimensional array
     *  of coefficients. This array may have only one row or
     *  column, but the structure should be 2D. This behaviour is
     *  hard-coded in the JSON file reader that resides inside
     *  the IncompressibleLibrary.cpp
     *  All other functions, also polyoffset, can only handle 1D
     *  input and throw an error if you feed them other coefficients.
     */

    /// Density coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData density;
    /// Specific heat coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     *  Fails for all other forms than polynomial due to the automatic
     *  integration for internal energy and entropy.
     */
    IncompressibleData specific_heat;
    /// Viscosity coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData viscosity;
    /// Conductivity coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData conductivity;
    /// Saturation pressure coefficients
    /** If 2D, the rows are temperature and the columns are concentration.
     *  If 1D, should be a column vector of temperature coefficients
     */
    IncompressibleData p_sat;
    /// Freezing temperature coefficients
    /** If 2D, the rows are concentration and the columns are pressure.
     *  If 1D, should be a column vector of concentration coefficients
     */
    IncompressibleData T_freeze;

    /// Mass fraction conversion coefficients
    /** If the fluid type is mass-based, it does not do anything. Otherwise,
     *  it converts the mass fraction to the required input.
     */
    IncompressibleData mass2input;
    /// Volume fraction conversion coefficients
    /** If the fluid type is volume-based, it does not do anything. Otherwise,
     *  it converts the volume fraction to the required input.
     */
    IncompressibleData volume2input;
    /// Mole fraction conversion coefficients
    /** If the fluid type is mole-based, it does not do anything. Otherwise,
     *  it converts the mole fraction to the required input.
     */
    IncompressibleData mole2input;

    Polynomial2DFrac poly;

    // Forward declaration of the some internal functions
    //double h_u(double T, double p, double x);
    //double u_h(double T, double p, double x);

   public:
    IncompressibleFluid() : Tmin(_HUGE), Tmax(_HUGE), xmin(_HUGE), xmax(_HUGE), TminPsat(_HUGE), xbase(_HUGE), Tbase(_HUGE) {
        strict = true;
        xid = IFRAC_UNDEFINED;
    };
    virtual ~IncompressibleFluid(){};

    std::string getName() const {
        return name;
    }
    std::string get_name() const {
        return getName();
    }  // For backwards-compatibility.
    std::string getDescription() const {
        return description;
    }
    std::string getReference() const {
        return reference;
    }

    double getTmax() const {
        return Tmax;
    }
    double getTmin() const {
        return Tmin;
    }
    double getxmax() const {
        return xmax;
    }
    double getxmin() const {
        return xmin;
    }
    composition_types getxid() const {
        return xid;
    }
    double getTminPsat() const {
        return TminPsat;
    }
    double getTbase() const {
        return Tbase;
    }
    double getxbase() const {
        return xbase;
    }

    void setName(const std::string& name) {
        this->name = name;
    }
    void setDescription(const std::string& description) {
        this->description = description;
    }
    void setReference(const std::string& reference) {
        this->reference = reference;
    }
    void setTmax(double Tmax) {
        this->Tmax = Tmax;
    }
    void setTmin(double Tmin) {
        this->Tmin = Tmin;
    }
    void setxmax(double xmax) {
        this->xmax = xmax;
    }
    void setxmin(double xmin) {
        this->xmin = xmin;
    }
    void setxid(composition_types xid) {
        this->xid = xid;
    }
    void setTminPsat(double TminPsat) {
        this->TminPsat = TminPsat;
    }
    void setTbase(double Tbase) {
        this->Tbase = Tbase;
    }
    void setxbase(double xbase) {
        this->xbase = xbase;
    }

    /// Setters for the coefficients
    void setDensity(IncompressibleData density) {
        this->density = density;
    }
    void setSpecificHeat(IncompressibleData specific_heat) {
        this->specific_heat = specific_heat;
    }
    void setViscosity(IncompressibleData viscosity) {
        this->viscosity = viscosity;
    }
    void setConductivity(IncompressibleData conductivity) {
        this->conductivity = conductivity;
    }
    void setPsat(IncompressibleData p_sat) {
        this->p_sat = p_sat;
    }
    void setTfreeze(IncompressibleData T_freeze) {
        this->T_freeze = T_freeze;
    }

    /// Setters for the concentration conversion coefficients
    void setMass2input(IncompressibleData mass2input) {
        this->mass2input = mass2input;
    }
    void setVolume2input(IncompressibleData volume2input) {
        this->volume2input = volume2input;
    }
    void setMole2input(IncompressibleData mole2input) {
        this->mole2input = mole2input;
    }

    /// A function to check coefficients and equation types.
    void validate();
    /// A function to test the density coefficients for 1D or 2D
    bool is_pure();

   protected:
    /// Base functions that handle the custom function types
    double baseExponential(IncompressibleData data, double y, double ybase);
    double baseLogexponential(IncompressibleData data, double y, double ybase);
    double baseExponentialOffset(IncompressibleData data, double y);
    double basePolyOffset(IncompressibleData data, double y, double z = 0.0);

   public:
    /* All functions need T and p as input. Might not
     * be necessary, but gives a clearer structure.
     */
    /// Density as a function of temperature, pressure and composition.
    double rho(double T, double p, double x);
    /// Heat capacities as a function of temperature, pressure and composition.
    double c(double T, double p, double x);
    double cp(double T, double p, double x) {
        throw ValueError(format("%s (%d): Please use the c-function instead.", __FILE__, __LINE__));
    }
    double cv(double T, double p, double x) {
        throw ValueError(format("%s (%d): Please use the c-function instead.", __FILE__, __LINE__));
    }
    /// Entropy as a function of temperature, pressure and composition.
    double s(double T, double p, double x) {
        throw ValueError(format("%s (%d): The internal calculations have changed, use the backend to calculate entropy from the partial derivatives.",
                                __FILE__, __LINE__));
    }
    /// Internal energy as a function of temperature, pressure and composition.
    double u(double T, double p, double x) {
        throw ValueError(
          format("%s (%d): The internal calculations have changed, use the backend to calculate internal energy from enthalpy.", __FILE__, __LINE__));
    }
    /// Enthalpy as a function of temperature, pressure and composition.
    double h(double T, double p, double x) {
        throw ValueError(
          format("%s (%d): The internal calculations have changed, use the backend to calculate enthalpy from the partial derivatives.", __FILE__,
                 __LINE__));
    }
    /// Viscosity as a function of temperature, pressure and composition.
    double visc(double T, double p, double x);
    /// Thermal conductivity as a function of temperature, pressure and composition.
    double cond(double T, double p, double x);
    /// Saturation pressure as a function of temperature and composition.
    double psat(double T, double x);
    /// Freezing temperature as a function of pressure and composition.
    double Tfreeze(double p, double x);

    /* Below are direct calculations of the derivatives. Nothing
     * special is going on, we simply use the polynomial class to
     * derive the different functions with respect to temperature.
     */
    /// Partial derivative of density
    //  with respect to temperature at constant pressure and composition
    double drhodTatPx(double T, double p, double x);
    ///// Partial derivative of entropy
    ////  with respect to temperature at constant pressure and composition
    //double dsdTatPx  (double T, double p, double x){return c(T,p,x)/T;};
    ///// Partial derivative of enthalpy
    ////  with respect to temperature at constant pressure and composition
    //double dhdTatPx  (double T, double p, double x){return c(T,p,x);};
    /// Partial derivative of entropy
    //  with respect to temperature at constant pressure and composition
    //  integrated in temperature
    double dsdTatPxdT(double T, double p, double x);
    /// Partial derivative of enthalpy
    //  with respect to temperature at constant pressure and composition
    //  integrated in temperature
    double dhdTatPxdT(double T, double p, double x);

    /// Mass fraction conversion function
    /** If the fluid type is mass-based, it does not do anything. Otherwise,
     *  it converts the mass fraction to the required input. */
    double inputFromMass(double T, double x);
    /// Volume fraction conversion function
    /** If the fluid type is volume-based, it does not do anything. Otherwise,
     *  it converts the volume fraction to the required input. */
    double inputFromVolume(double T, double x);
    /// Mole fraction conversion function
    /** If the fluid type is mole-based, it does not do anything. Otherwise,
     *  it converts the mole fraction to the required input. */
    double inputFromMole(double T, double x);

    /* Some functions can be inverted directly, those are listed
     * here. It is also possible to solve for other quantities, but
     * that involves some more sophisticated processing and is not
     * done here, but in the backend, T(h,p) for example.
     */
    /// Temperature as a function of density, pressure and composition.
    double T_rho(double Dmass, double p, double x);
    /// Temperature as a function of heat capacities as a function of temperature, pressure and composition.
    double T_c(double Cmass, double p, double x);
    /// Temperature as a function of entropy as a function of temperature, pressure and composition.
    double T_s(double Smass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from entropy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Temperature as a function of internal energy as a function of temperature, pressure and composition.
    double T_u(double Umass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from internal energy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Temperature as a function of enthalpy, pressure and composition.
    double T_h(double Hmass, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from enthalpy is not implemented in the fluid, use the backend.", __FILE__, __LINE__));
    }
    /// Viscosity as a function of temperature, pressure and composition.
    double T_visc(double visc, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from viscosity is not implemented.", __FILE__, __LINE__));
    }
    /// Thermal conductivity as a function of temperature, pressure and composition.
    double T_cond(double cond, double p, double x) {
        throw NotImplementedError(format("%s (%d): T from conductivity is not implemented.", __FILE__, __LINE__));
    }
    /// Saturation pressure as a function of temperature and composition.
    double T_psat(double psat, double x) {
        throw NotImplementedError(format("%s (%d): T from psat is not implemented.", __FILE__, __LINE__));
    }
    /// Composition as a function of freezing temperature and pressure.
    double x_Tfreeze(double Tfreeze, double p) {
        throw NotImplementedError(format("%s (%d): x from T_freeze is not implemented.", __FILE__, __LINE__));
    }

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
        return u(T, p, x) + p / rho(T, p, x);
    };

    /// Internal energy from h, p and rho.
    /** Calculate internal energy as a function of temperature
     *  and pressure employing functions for enthalpy and
     *  density. Provides consistent formulations. */
    double u_h(double T, double p, double x) {
        return h(T, p, x) - p / rho(T, p, x);
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

   public:
    /// Check validity of composition input.
    /** Compares the given composition x to a stored minimum and
     *  maximum value. Enforces the redefinition of xmin and
     *  xmax since the default values cause an error. */
    bool checkX(double x);

    /// Check validity of temperature, pressure and composition input.
    bool checkTPX(double T, double p, double x) {
        return (checkT(T, p, x) && checkP(T, p, x) && checkX(x));
    };
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEFLUID_H_ */
