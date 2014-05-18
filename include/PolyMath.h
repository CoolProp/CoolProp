#ifndef POLYMATH_H
#define POLYMATH_H

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
#include "Solvers.h"
//#include <numeric> // inner_product
//#include <sstream>
//#include "float.h"

namespace CoolProp{

/// The base class for Polynomials
class BasePolynomial{

protected:
	bool DEBUG;

public:
	// Constructor
	BasePolynomial();
	// Destructor.  No implementation
	virtual ~BasePolynomial(){};

public:
	/// Basic checks for coefficient vectors.
	/** Starts with only the first coefficient dimension
	 *  and checks the vector length against parameter n. */
	bool checkCoefficients(const std::vector<double> &coefficients, const unsigned int n);
	bool checkCoefficients(const std::vector< std::vector<double> > &coefficients, const unsigned int rows, const unsigned int columns);

	/** Integrating coefficients for polynomials is done by dividing the
	 *  original coefficients by (i+1) and elevating the order by 1
	 *  through adding a zero as first coefficient.
	 *  Some reslicing needs to be applied to integrate along the x-axis.
	 *  In the brine/solution equations, reordering of the parameters
	 *  avoids this expensive operation. However, it is included for the
	 *  sake of completeness.
	 */
	std::vector<double> integrateCoeffs(const std::vector<double> &coefficients);
	std::vector< std::vector<double> > integrateCoeffs(const std::vector< std::vector<double> > &coefficients, bool axis);

	/** Deriving coefficients for polynomials is done by multiplying the
	 *  original coefficients with i and lowering the order by 1.
	 *
	 *  It is not really deprecated, but untested and therefore a warning
	 *  is issued. Please check this method before you use it.
	 */
	std::vector<double> deriveCoeffs(const std::vector<double> &coefficients);
	std::vector< std::vector<double> > deriveCoeffs(const std::vector< std::vector<double> > &coefficients, unsigned int axis);

private:
	/** The core of the polynomial wrappers are the different
	 *  implementations that follow below. In case there are
	 *  new calculation schemes available, please do not delete
	 *  the implementations, but mark them as deprecated.
	 *  The old functions are good for debugging since the
	 *  structure is easier to read than the backward Horner-scheme
	 *  or the recursive Horner-scheme.
	 */

	/// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
	/** Base function to produce n-th order polynomials
	 *  based on the length of the coefficient vector.
	 *  Starts with only the first coefficient at x^0. */
	DEPRECATED(double simplePolynomial(const std::vector<double> &coefficients, double x));
	DEPRECATED(double simplePolynomial(const std::vector<std::vector<double> > &coefficients, double x, double y));

	/// Simple integrated polynomial function generator.
	/** Base function to produce integrals of n-th order polynomials based on
	 *  the length of the coefficient vector.
	 *  Starts with only the first coefficient at x^0 */
	///Indefinite integral in x-direction
	double simplePolynomialInt(const std::vector<double> &coefficients, double x);
	///Indefinite integral in y-direction only
	double simplePolynomialInt(const std::vector<std::vector<double> > &coefficients, double x, double y);

	/// Simple integrated polynomial function generator divided by independent variable.
	/** Base function to produce integrals of n-th order
	 *  polynomials based on the length of the coefficient
	 *  vector. Starts with only the first coefficient at x^0 */
	///Indefinite integral of a polynomial divided by its independent variable
	double simpleFracInt(const std::vector<double> &coefficients, double x);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	double simpleFracInt(const std::vector<std::vector<double> > &coefficients, double x, double y);

	/** Simple integrated centred(!) polynomial function generator divided by independent variable.
	 *  We need to rewrite some of the functions in order to
	 *  use central fit. Having a central temperature xbase
	 *  allows for a better fit, but requires a different
	 *  formulation of the fracInt function group. Other
	 *  functions are not affected.
	 *  Starts with only the first coefficient at x^0 */
	///Helper function to calculate the D vector:
	double factorial(double nValue);
	double binom(double nValue, double nValue2);
	std::vector<double> fracIntCentralDvector(int m, double x, double xbase);
	///Indefinite integral of a centred polynomial divided by its independent variable
	double fracIntCentral(const std::vector<double> &coefficients, double x, double xbase);

	/// Horner function generator implementations
	/** Represent polynomials according to Horner's scheme.
	 *  This avoids unnecessary multiplication and thus
	 *  speeds up calculation.
	 */
	double baseHorner(const std::vector<double> &coefficients, double x);
	double baseHorner(const std::vector< std::vector<double> > &coefficients, double x, double y);
	///Indefinite integral in x-direction
	double baseHornerInt(const std::vector<double> &coefficients, double x);
	///Indefinite integral in y-direction only
	double baseHornerInt(const std::vector<std::vector<double> > &coefficients, double x, double y);
	///Indefinite integral of a polynomial divided by its independent variable
	double baseHornerFracInt(const std::vector<double> &coefficients, double x);
	///Indefinite integral of a polynomial divided by its 2nd independent variable
	double baseHornerFracInt(const std::vector<std::vector<double> > &coefficients, double x, double y);

	/** Alternatives
	 *  Simple functions that heavily rely on other parts of this file.
	 *  We still need to check which combinations yield the best
	 *  performance.
	 */
	///Derivative in x-direction
	double deriveIn2Steps(const std::vector<double> &coefficients, double x); // TODO: Check results!
	///Derivative in terms of x(axis=true) or y(axis=false).
	double deriveIn2Steps(const std::vector< std::vector<double> > &coefficients, double x, double y, bool axis); // TODO: Check results!
	///Indefinite integral in x-direction
	double integrateIn2Steps(const std::vector<double> &coefficients, double x);
	///Indefinite integral in terms of x(axis=true) or y(axis=false).
	double integrateIn2Steps(const std::vector< std::vector<double> > &coefficients, double x, double y, bool axis);
	///Indefinite integral in x-direction of a polynomial divided by its independent variable
	double fracIntIn2Steps(const std::vector<double> &coefficients, double x);
	///Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
	double fracIntIn2Steps(const std::vector<std::vector<double> > &coefficients, double x, double y);
	///Indefinite integral of a centred polynomial divided by its 2nd independent variable
	double fracIntCentral2Steps(const std::vector<std::vector<double> > &coefficients, double x, double y, double ybase);

public:
	/** Here we define the functions that should be used by the
	 *  respective implementations. Please do no use any other
	 *  method since this would break the purpose of this interface.
	 *  Note that the functions below are supposed to be aliases
	 *  to implementations declared elsewhere in this file.
	 */

	/** Everything related to the normal polynomials goes in this
	 *  section, holds all the functions for evaluating polynomials.
	 */
	/// Evaluates a one-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	virtual inline double polyval(const std::vector<double> &coefficients, double x){
		return baseHorner(coefficients,x);
	}

	/// Evaluates a two-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyval(const std::vector< std::vector<double> > &coefficients, double x, double y){
		return baseHorner(coefficients,x,y);
	}


	/** Everything related to the integrated polynomials goes in this
	 *  section, holds all the functions for evaluating polynomials.
	 */
	/// Evaluates the indefinite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	virtual inline double polyint(const std::vector<double> &coefficients, double x){
		return baseHornerInt(coefficients,x);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyint(const std::vector< std::vector<double> > &coefficients, double x, double y){
		return baseHornerInt(coefficients,x,y);
	}


	/** Everything related to the derived polynomials goes in this
	 *  section, holds all the functions for evaluating polynomials.
	 */
	/// Evaluates the derivative of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	virtual inline double polyder(const std::vector<double> &coefficients, double x){
		return deriveIn2Steps(coefficients,x);
	}

	/// Evaluates the derivative of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyder(const std::vector< std::vector<double> > &coefficients, double x, double y){
		return deriveIn2Steps(coefficients,x,y,false);
	}


	/** Everything related to the polynomials divided by one variable goes in this
	 *  section, holds all the functions for evaluating polynomials.
	 */
	/// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current position
	virtual inline double polyfracval(const std::vector<double> &coefficients, double x){
		return baseHorner(coefficients,x)/x;
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double y){
		return baseHorner(coefficients,x,y)/y;
	}


	/** Everything related to the integrated polynomials divided by one variable goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Evaluates the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current position
	virtual inline double polyfracint(const std::vector<double> &coefficients, double x){
		return baseHornerFracInt(coefficients,x);
	}

	/// Evaluates the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double y){
		return baseHornerFracInt(coefficients,x,y);
	}

	/// Evaluates the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current position
	/// @param xbase central temperature for fitted function
	virtual inline double polyfracintcentral(const std::vector<double> &coefficients, double x, double xbase){
		return fracIntCentral(coefficients,x,xbase);
	}

	/// Evaluates the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	/// @param ybase central temperature for fitted function
	virtual inline double polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double y, double ybase){
		return fracIntCentral2Steps(coefficients,x,y,ybase);
	}


	/** Everything related to the derived polynomials divided by one variable goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Evaluates the derivative of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current position
	virtual inline double polyfracder(const std::vector<double> &coefficients, double x){
		throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracder1D
	}

	/// Evaluates the derivative of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	virtual inline double polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double y){
		throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracder2D
	}

	/// Evaluates the derivative of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current position
	/// @param xbase central temperature for fitted function
	virtual inline double polyfracdercentral(const std::vector<double> &coefficients, double x, double xbase){
		throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracdercentral1D
	}

	/// Evaluates the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	/// @param ybase central temperature for fitted function
	virtual inline double polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double y, double ybase){
		throw CoolProp::NotImplementedError("Derivatives of polynomials divided by their independent variable have not been implemented."); // TODO: Implement polyfracdercentral2D
	}
};




/** Implements the function wrapper interface and can be
 *  used by the solvers.
 *  TODO: Make multidimensional
 */
class PolyResidual : public FuncWrapper1D {
protected:
	enum dims {i1D, i2D};
	/// Object that evaluates the equation
	BasePolynomial poly;
	/// Current output value
	double output, firstDim;
	int dim;
	std::vector< std::vector<double> > coefficients;
private:
	PolyResidual();
public:
	PolyResidual(const std::vector<double> &coefficients, double y);
	PolyResidual(const std::vector< std::vector<double> > &coefficients, double x, double z);
	virtual ~PolyResidual(){};
	bool is2D(){return (this->dim==i2D);};
	virtual double call(double x);
	virtual double deriv(double x);
};
class PolyIntResidual : public PolyResidual {
public:
	PolyIntResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
	PolyIntResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
	virtual double call(double x);
	virtual double deriv(double x);
};
class PolyFracIntResidual : public PolyResidual {
public:
	PolyFracIntResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
	PolyFracIntResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
	virtual double call(double x);
	virtual double deriv(double x);
};
class PolyFracIntCentralResidual : public PolyResidual {
protected:
	double baseVal;
public:
	PolyFracIntCentralResidual(const std::vector<double> &coefficients, double y, double xBase):PolyResidual(coefficients, y){this->baseVal = xBase;};
	PolyFracIntCentralResidual(const std::vector< std::vector<double> > &coefficients, double x, double z, double yBase): PolyResidual(coefficients, x, z){this->baseVal = yBase;};
	virtual double call(double x);
	virtual double deriv(double x);
};
class PolyDerResidual : public PolyResidual {
public:
	PolyDerResidual(const std::vector<double> &coefficients, double y):PolyResidual(coefficients, y){};
	PolyDerResidual(const std::vector< std::vector<double> > &coefficients, double x, double z):PolyResidual(coefficients, x, z){};
	virtual double call(double x);
	virtual double deriv(double x);
};




/** Implements the same public functions as the
 *  but solves the polynomial for the given value
 *  instead of evaluating it.
 *  TODO: This class does not check for bijective
 *        polynomials and is therefore a little
 *        fragile.
 */
class PolynomialSolver : public BasePolynomial{
private:
	enum solvers {iNewton, iBrent};
	int uses;
	double guess, min, max;
	double macheps, tol;
	int maxiter;

public:
	// Constructor
	PolynomialSolver();
	// Destructor.  No implementation
	virtual ~PolynomialSolver(){};

public:
	/** Here we redefine the functions that solve the polynomials.
	 *  These implementations all use the base class to evaluate
	 *  the polynomial during the solution process.
	 */

	/** Everything related to the normal polynomials goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves a one-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current input
	virtual double polyval(const std::vector<double> &coefficients, double y);

	/// Solves a two-dimensional polynomial for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyval(const std::vector< std::vector<double> > &coefficients, double x, double z);


	/** Everything related to the integrated polynomials goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves the indefinite integral of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	virtual double polyint(const std::vector<double> &coefficients, double y);

	/// Solves the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyint(const std::vector< std::vector<double> > &coefficients, double x, double z);


	/** Everything related to the derived polynomials goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves the derivative of a one-dimensional polynomial
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	virtual double polyder(const std::vector<double> &coefficients, double y);

	/// Solves the derivative of a two-dimensional polynomial along the 2nd axis (y)
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyder(const std::vector< std::vector<double> > &coefficients, double x, double z);


	/** Everything related to the polynomials divided by one variable goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	virtual double polyfracval(const std::vector<double> &coefficients, double y);

	/// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double z);


	/** Everything related to the integrated polynomials divided by one variable goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	virtual double polyfracint(const std::vector<double> &coefficients, double y);

	/// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double z);

	/// Solves the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	/// @param xbase central x-value for fitted function
	virtual double polyfracintcentral(const std::vector<double> &coefficients, double y, double xbase);

	/// Solves the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	/// @param ybase central y-value for fitted function
	virtual double polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase);


	/** Everything related to the derived polynomials divided by one variable goes in this
	 *  section, holds all the functions for solving polynomials.
	 */
	/// Solves the derivative of a one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	virtual double polyfracder(const std::vector<double> &coefficients, double y);

	/// Solves the derivative of a two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	virtual double polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double z);

	/// Solves the derivative of a centred one-dimensional polynomial divided by its independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param y double value that represents the current output
	/// @param xbase central x-value for fitted function
	virtual double polyfracdercentral(const std::vector<double> &coefficients, double y, double xbase);

	/// Solves the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param z double value that represents the current output
	/// @param ybase central y-value for fitted function
	virtual double polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase);


	/** Set the solvers and updates either the guess values or the
	 *  boundaries for the variable to solve for.
	 */
	/// Sets the guess value for the Newton solver and enables it.
	/// @param guess double value that represents the guess value
	virtual void setGuess(double guess);
	/// Sets the limits for the Brent solver and enables it.
	/// @param min double value that represents the lower boundary
	/// @param max double value that represents the upper boundary
	virtual void setLimits(double min, double max);
	/// Solves the equations based on previously defined parameters.
	/// @param min double value that represents the lower boundary
	/// @param max double value that represents the upper boundary
	virtual double solve(PolyResidual &res);
};


/// The base class for exponential functions
class BaseExponential{

protected:
	BasePolynomial poly;
	bool DEBUG;

public:
	BaseExponential();
	virtual ~BaseExponential(){};

public:
	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input
	/// @param n int value that determines the kind of exponential function
	double expval(const std::vector<double> &coefficients, double x, int n);

	/// Evaluates an exponential function for the given coefficients
	/// @param coefficients vector containing the ordered coefficients
	/// @param x double value that represents the current input in the 1st dimension
	/// @param y double value that represents the current input in the 2nd dimension
	/// @param n int value that determines the kind of exponential function
	double expval(const std::vector< std::vector<double> > &coefficients, double x, double y, int n);
};


}; /* namespace CoolProp */
#endif
