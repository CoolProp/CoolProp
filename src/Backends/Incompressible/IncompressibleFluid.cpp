
#include "DataStructures.h"
#include "IncompressibleFluid.h"
#include "IncompressibleLibrary.h"
#include "math.h"
#include "MatrixMath.h"
#include "PolyMath.h"
#include <Eigen/Core>

namespace CoolProp {

/// A thermophysical property provider for all properties
/**
This fluid instance is populated using an entry from a JSON file and uses
simplified polynomial and exponential functions to calculate thermophysical
and transport properties.
*/
//IncompressibleFluid::IncompressibleFluid();

void IncompressibleFluid::validate() {
    return;
    // TODO: Implement validation function

    // u and s have to be of the polynomial type!
    //throw NotImplementedError("TODO");
}

bool IncompressibleFluid::is_pure() {
    if (density.coeffs.cols() == 1) return true;
    return false;
}

/// Base exponential function
double IncompressibleFluid::baseExponential(IncompressibleData data, double y, double ybase) {
    Eigen::VectorXd coeffs = makeVector(data.coeffs);
    size_t r = coeffs.rows(), c = coeffs.cols();
    if (strict && (r != 3 || c != 1))
        throw ValueError(format("%s (%d): You have to provide a 3,1 matrix of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
    return exp((double)(coeffs[0] / ((y - ybase) + coeffs[1]) - coeffs[2]));
}
/// Base exponential function with logarithmic term
double IncompressibleFluid::baseLogexponential(IncompressibleData data, double y, double ybase) {
    Eigen::VectorXd coeffs = makeVector(data.coeffs);
    size_t r = coeffs.rows(), c = coeffs.cols();
    if (strict && (r != 3 || c != 1))
        throw ValueError(format("%s (%d): You have to provide a 3,1 matrix of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
    return exp(
      (double)(log((double)(1.0 / ((y - ybase) + coeffs[0]) + 1.0 / ((y - ybase) + coeffs[0]) / ((y - ybase) + coeffs[0]))) * coeffs[1] + coeffs[2]));
}

double IncompressibleFluid::basePolyOffset(IncompressibleData data, double y, double z) {
    size_t r = data.coeffs.rows(), c = data.coeffs.cols();
    double offset = 0.0;
    double in = 0.0;
    Eigen::MatrixXd coeffs;
    if (r > 0 && c > 0) {
        offset = data.coeffs(0, 0);
        if (r == 1 && c > 1) {  // row vector -> function of z
            coeffs = Eigen::MatrixXd(data.coeffs.block(0, 1, r, c - 1));
            in = z;
        } else if (r > 1 && c == 1) {  // column vector -> function of y
            coeffs = Eigen::MatrixXd(data.coeffs.block(1, 0, r - 1, c));
            in = y;
        } else {
            throw ValueError(format("%s (%d): You have to provide a vector (1D matrix) of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
        }
        return poly.evaluate(coeffs, in, 0, offset);
    }
    throw ValueError(format("%s (%d): You have to provide a vector (1D matrix) of coefficients, not  (%d,%d).", __FILE__, __LINE__, r, c));
}

/// Density as a function of temperature, pressure and composition.
double IncompressibleFluid::rho(double T, double p, double x) {
    switch (density.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(density.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(density, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(density, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(density.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(density, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, density.type));
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, density.type));
    }
}

/// Heat capacities as a function of temperature, pressure and composition.
double IncompressibleFluid::c(double T, double p, double x) {
    switch (specific_heat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            //throw NotImplementedError("Here you should implement the polynomial.");
            return poly.evaluate(specific_heat.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, specific_heat.type));
        default:
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for specific heat.", __FILE__, __LINE__,
                                    specific_heat.type));
    }
}

/// Viscosity as a function of temperature, pressure and composition.
double IncompressibleFluid::visc(double T, double p, double x) {
    switch (viscosity.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(viscosity.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(viscosity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(viscosity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(viscosity.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(viscosity, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, viscosity.type));
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, viscosity.type));
    }
}
/// Thermal conductivity as a function of temperature, pressure and composition.
double IncompressibleFluid::cond(double T, double p, double x) {
    switch (conductivity.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(conductivity.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(conductivity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(conductivity, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(conductivity.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(conductivity, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, conductivity.type));
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, conductivity.type));
    }
}
/// Saturation pressure as a function of temperature and composition.
double IncompressibleFluid::psat(double T, double x) {
    if (T <= this->TminPsat) return 0.0;
    switch (p_sat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(p_sat.coeffs, T, x, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(p_sat, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(p_sat, T, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(p_sat.coeffs, T, x, 0, 0, Tbase, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(p_sat, T, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, p_sat.type));
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, p_sat.type));
    }
}
/// Freezing temperature as a function of pressure and composition.
double IncompressibleFluid::Tfreeze(double p, double x) {
    switch (T_freeze.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.evaluate(T_freeze.coeffs, p, x, 0, 0, 0.0, xbase);
        case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
            return baseExponential(T_freeze, x, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
            return baseLogexponential(T_freeze, x, 0.0);
        case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return exp(poly.evaluate(T_freeze.coeffs, p, x, 0, 0, 0.0, xbase));
        case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
            return basePolyOffset(T_freeze, p, x);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, T_freeze.type));
        default:
            throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.", __FILE__, __LINE__, T_freeze.type));
    }
}

/* Below are direct calculations of the derivatives. Nothing
 * special is going on, we simply use the polynomial class to
 * derive the different functions with respect to temperature.
 */
/// Partial derivative of density with respect to temperature at constant pressure and composition
double IncompressibleFluid::drhodTatPx(double T, double p, double x) {
    switch (density.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.derivative(density.coeffs, T, x, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, density.type));
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for density.", __FILE__, __LINE__, density.type));
    }
}
/// Partial derivative of entropy
//  with respect to temperature at constant pressure and composition
//  integrated in temperature
double IncompressibleFluid::dsdTatPxdT(double T, double p, double x) {
    switch (specific_heat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.integral(specific_heat.coeffs, T, x, 0, -1, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, specific_heat.type));
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for entropy.", __FILE__, __LINE__, specific_heat.type));
    }
}
/// Partial derivative of enthalpy
//  with respect to temperature at constant pressure and composition
//  integrated in temperature
double IncompressibleFluid::dhdTatPxdT(double T, double p, double x) {
    switch (specific_heat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.integral(specific_heat.coeffs, T, x, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, specific_heat.type));
        default:
            throw ValueError(
              format("%s (%d): There is no predefined way to use this function type \"[%d]\" for entropy.", __FILE__, __LINE__, specific_heat.type));
    }
}

/// Mass fraction conversion function
/** If the fluid type is mass-based, it does not do anything. Otherwise,
 *  it converts the mass fraction to the required input. */
double IncompressibleFluid::inputFromMass(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_MASS) {
        return x;
    } else {
        throw NotImplementedError("Mass composition conversion has not been implemented.");
        //switch (mass2input.type) {
        //    case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
        //        return poly.evaluate(mass2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
        //        return baseExponential(mass2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
        //        return baseLogexponential(mass2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
        //        return exp(poly.evaluate(mass2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
        //        return basePolyOffset(mass2input, T, x);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
        //        throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,mass2input.type));
        //        break;
        //    default:
        //        throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,mass2input.type));
        //        break;
        //}
        //return _HUGE;
    }
}

/// Volume fraction conversion function
/** If the fluid type is volume-based, it does not do anything. Otherwise,
 *  it converts the volume fraction to the required input. */
double IncompressibleFluid::inputFromVolume(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_VOLUME) {
        return x;
    } else {
        throw NotImplementedError("Volume composition conversion has not been implemented.");
        //switch (volume2input.type) {
        //    case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
        //        return poly.evaluate(volume2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
        //        return baseExponential(volume2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
        //        return baseLogexponential(volume2input, x, 0.0);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
        //        return exp(poly.evaluate(volume2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
        //        return basePolyOffset(volume2input, T, x);
        //        break;
        //    case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
        //        throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,volume2input.type));
        //        break;
        //    default:
        //        throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,volume2input.type));
        //        break;
        //}
        //return _HUGE;
    }
}

/// Mole fraction conversion function
/** If the fluid type is mole-based, it does not do anything. Otherwise,
 *  it converts the mole fraction to the required input. */
double IncompressibleFluid::inputFromMole(double T, double x) {
    if (this->xid == IFRAC_PURE) {
        return _HUGE;
    } else if (this->xid == IFRAC_MOLE) {
        return x;
    } else {
        throw NotImplementedError("Mole composition conversion has not been implemented.");
        /*
        switch (mole2input.type) {
            case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
                return poly.evaluate(mole2input.coeffs, T, x, 0, 0, 0.0, 0.0); // TODO: make sure Tbase and xbase are defined in the correct way
                break;
            case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
                return baseExponential(mole2input, x, 0.0);
                break;
            case IncompressibleData::INCOMPRESSIBLE_LOGEXPONENTIAL:
                return baseLogexponential(mole2input, x, 0.0);
                break;
            case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
                return exp(poly.evaluate(mole2input.coeffs, T, x, 0, 0, 0.0, 0.0)); // TODO: make sure Tbase and xbase are defined in the correct way
                break;
            case IncompressibleData::INCOMPRESSIBLE_POLYOFFSET:
                return basePolyOffset(mole2input, T, x);
                break;
            case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
                throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,mole2input.type));
                break;
            default:
                throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,mole2input.type));
                break;
        }
        return _HUGE;
        */
    }
}

/* Some functions can be inverted directly, those are listed
 * here. It is also possible to solve for other quantities, but
 * that involves some more sophisticated processing and is not
 * done here, but in the backend, T(h,p) for example.
 */
/// Temperature as a function of density, pressure and composition.
double IncompressibleFluid::T_rho(double Dmass, double p, double x) {
    double d_raw = Dmass;  // No changes needed, no reference values...
    switch (density.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.solve_limits(density.coeffs, x, d_raw, Tmin, Tmax, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, specific_heat.type));
        default:
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for inverse density.", __FILE__, __LINE__,
                                    specific_heat.type));
    }
}
/// Temperature as a function of heat capacities as a function of temperature, pressure and composition.
double IncompressibleFluid::T_c(double Cmass, double p, double x) {
    double c_raw = Cmass;  // No changes needed, no reference values...
    switch (specific_heat.type) {
        case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
            return poly.solve_limits(specific_heat.coeffs, x, c_raw, Tmin, Tmax, 0, 0, 0, Tbase, xbase);
        case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
            throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?", __FILE__,
                                    __LINE__, specific_heat.type));
        default:
            throw ValueError(format("%s (%d): There is no predefined way to use this function type \"[%d]\" for inverse specific heat.", __FILE__,
                                    __LINE__, specific_heat.type));
    }
}

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
bool IncompressibleFluid::checkT(double T, double p, double x) {
    if (Tmin <= 0.) throw ValueError("Please specify the minimum temperature.");
    if (Tmax <= 0.) throw ValueError("Please specify the maximum temperature.");
    if ((Tmin > T) || (T > Tmax)) throw ValueError(format("Your temperature %f is not between %f and %f.", T, Tmin, Tmax));
    double TF = 0.0;
    if (T_freeze.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET) TF = Tfreeze(p, x);
    if (T < TF) throw ValueError(format("Your temperature %f is below the freezing point of %f.", T, TF));
    return true;
}

/// Check validity of pressure input.
/** Compares the given pressure p to the saturation pressure at
 *  temperature T and throws and exception if p is lower than
 *  the saturation conditions.
 *  The default value for psat is -1 yielding true if psat
 *  is not redefined in the subclass.
 *  */
bool IncompressibleFluid::checkP(double T, double p, double x) {
    double ps = 0.0;
    if (p_sat.type != IncompressibleData::INCOMPRESSIBLE_NOT_SET) ps = psat(T, x);
    if (p < 0.0) throw ValueError(format("You cannot use negative pressures: %f < %f. ", p, 0.0));
    if (ps > 0.0 && p < ps) throw ValueError(format("Equations are valid for liquid phase only: %f < %f (psat). ", p, ps));
    return true;
}

/// Check validity of composition input.
/** Compares the given composition x to a stored minimum and
 *  maximum value. Enforces the redefinition of xmin and
 *  xmax since the default values cause an error. */
bool IncompressibleFluid::checkX(double x) {
    if (xmin < 0.0 || xmin > 1.0) throw ValueError("Please specify the minimum concentration between 0 and 1.");
    if (xmax < 0.0 || xmax > 1.0) throw ValueError("Please specify the maximum concentration between 0 and 1.");
    if ((xmin > x) || (x > xmax)) throw ValueError(format("Your composition %f is not between %f and %f.", x, xmin, xmax));
    return true;
}

} /* namespace CoolProp */

// Testing still needs to be enhanced.
/* Below, I try to carry out some basic tests for both 2D and 1D
 * polynomials as well as the exponential functions for vapour
 * pressure etc.
 */
#ifdef ENABLE_CATCH
#    include <math.h>
#    include <iostream>
#    include <catch2/catch_all.hpp>
#    include "TestObjects.h"

Eigen::MatrixXd makeMatrix(const std::vector<double>& coefficients) {
    //IncompressibleClass::checkCoefficients(coefficients,18);
    std::vector<std::vector<double>> matrix;
    std::vector<double> tmpVector;

    tmpVector.clear();
    tmpVector.push_back(coefficients[0]);
    tmpVector.push_back(coefficients[6]);
    tmpVector.push_back(coefficients[11]);
    tmpVector.push_back(coefficients[15]);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[1] * 100.0);
    tmpVector.push_back(coefficients[7] * 100.0);
    tmpVector.push_back(coefficients[12] * 100.0);
    tmpVector.push_back(coefficients[16] * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[2] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[8] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[13] * 100.0 * 100.0);
    tmpVector.push_back(coefficients[17] * 100.0 * 100.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[3] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[9] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[14] * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[4] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(coefficients[10] * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    tmpVector.push_back(coefficients[5] * 100.0 * 100.0 * 100.0 * 100.0 * 100.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    tmpVector.push_back(0.0);
    matrix.push_back(tmpVector);

    tmpVector.clear();
    return CoolProp::vec_to_eigen(matrix).transpose();
}

TEST_CASE("Internal consistency checks and example use cases for the incompressible fluids", "[IncompressibleFluids]") {
    bool PRINT = false;
    std::string tmpStr;
    std::vector<double> tmpVector;
    std::vector<std::vector<double>> tmpMatrix;

    SECTION("Test case for \"SylthermXLT\" by Dow Chemicals") {

        std::vector<double> cRho;
        cRho.push_back(+1.1563685145E+03);
        cRho.push_back(-1.0269048032E+00);
        cRho.push_back(-9.3506079577E-07);
        cRho.push_back(+1.0368116627E-09);
        CoolProp::IncompressibleData density;
        density.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        density.coeffs = CoolProp::vec_to_eigen(cRho);

        std::vector<double> cHeat;
        cHeat.push_back(+1.1562261074E+03);
        cHeat.push_back(+2.0994549103E+00);
        cHeat.push_back(+7.7175381057E-07);
        cHeat.push_back(-3.7008444051E-20);
        CoolProp::IncompressibleData specific_heat;
        specific_heat.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        specific_heat.coeffs = CoolProp::vec_to_eigen(cHeat);

        std::vector<double> cCond;
        cCond.push_back(+1.6121957379E-01);
        cCond.push_back(-1.3023781944E-04);
        cCond.push_back(-1.4395238766E-07);
        CoolProp::IncompressibleData conductivity;
        conductivity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        conductivity.coeffs = CoolProp::vec_to_eigen(cCond);

        std::vector<double> cVisc;
        cVisc.push_back(+1.0337654989E+03);
        cVisc.push_back(-4.3322764383E+01);
        cVisc.push_back(+1.0715062356E+01);
        CoolProp::IncompressibleData viscosity;
        viscosity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL;
        viscosity.coeffs = CoolProp::vec_to_eigen(cVisc);

        CoolProp::IncompressibleFluid XLT;
        XLT.setName("XLT");
        XLT.setDescription("SylthermXLT");
        XLT.setReference("Dow Chemicals data sheet");
        XLT.setTmax(533.15);
        XLT.setTmin(173.15);
        XLT.setxmax(0.0);
        XLT.setxmin(0.0);
        XLT.setTminPsat(533.15);

        XLT.setTbase(0.0);
        XLT.setxbase(0.0);

        /// Setters for the coefficients
        XLT.setDensity(density);
        XLT.setSpecificHeat(specific_heat);
        XLT.setViscosity(viscosity);
        XLT.setConductivity(conductivity);
        //XLT.setPsat(parse_coefficients(fluid_json, "saturation_pressure", false));
        //XLT.setTfreeze(parse_coefficients(fluid_json, "T_freeze", false));
        //XLT.setVolToMass(parse_coefficients(fluid_json, "volume2mass", false));
        //XLT.setMassToMole(parse_coefficients(fluid_json, "mass2mole", false));

        /// A function to check coefficients and equation types.
        //XLT.validate();
        double acc = 0.0001;
        double val = 0;
        double res = 0;

        // Prepare the results and compare them to the calculated values
        double T = 273.15 + 50;
        double p = 10e5;
        double x = 0.0;

        // Compare density
        val = 824.4615702148608;
        res = XLT.rho(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Compare cp
        val = 1834.7455527670554;
        res = XLT.c(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Check property functions
        CHECK_THROWS(XLT.s(T, p, x));
        CHECK_THROWS(XLT.h(T, p, x));
        CHECK_THROWS(XLT.u(T, p, x));

        // Compare v
        val = 0.0008931435169681835;
        res = XLT.visc(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Compare l
        val = 0.10410086156049088;
        res = XLT.cond(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }
    }

    SECTION("Test case for Methanol from SecCool") {

        CoolProp::IncompressibleFluid CH3OH = CoolProp::get_incompressible_fluid("CH3OH");

        // Prepare the results and compare them to the calculated values
        double acc = 0.0001;
        double T = 273.15 + 10;
        double p = 10e5;
        double x = 0.25;
        double expected = 0;
        double actual = 0;

        // Compare density
        expected = 963.2886528091547;
        actual = CH3OH.rho(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare cp
        expected = 3993.9748117022423;
        actual = CH3OH.c(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            CHECK(check_abs(expected, actual, acc));
        }

        // Check property functions
        CHECK_THROWS(CH3OH.s(T, p, x));
        CHECK_THROWS(CH3OH.h(T, p, x));
        CHECK_THROWS(CH3OH.u(T, p, x));

        // Compare v
        expected = 0.0023970245009602097;
        actual = CH3OH.visc(T, p, x) / 1e3;
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare conductivity
        expected = 0.44791148414693727;
        actual = CH3OH.cond(T, p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }

        // Compare Tfreeze
        expected = -20.02 + 273.15;  // 253.1293105454671;
        actual = CH3OH.Tfreeze(p, x);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }
    }
}

#endif /* ENABLE_CATCH */
