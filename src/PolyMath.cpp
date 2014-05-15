
#include "PolyMath.h"

#include "CoolPropTools.h"
#include "Exceptions.h"

#include <vector>
#include <string>
//#include <sstream>
//#include <numeric>
#include <math.h>

#include "Solvers.h"

namespace CoolProp{

BasePolynomial::BasePolynomial(){
  this->DEBUG = false;
}


/// Basic checks for coefficient vectors.
/** Starts with only the first coefficient dimension
 *  and checks the vector length against parameter n. */
bool BasePolynomial::checkCoefficients(const std::vector<double> &coefficients, unsigned int n){
	if (coefficients.size() == n){
		return true;
	} else {
		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
	}
	return false;
}
bool BasePolynomial::checkCoefficients(std::vector< std::vector<double> > const& coefficients, unsigned int rows, unsigned int columns){
	if (coefficients.size() == rows){
		bool result = true;
		for(unsigned int i=0; i<rows; i++) {
			result = result && checkCoefficients(coefficients[i],columns);
		}
		return result;
	} else {
		throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),rows));
	}
	return false;
}


/** Integrating coefficients for polynomials is done by dividing the
 *  original coefficients by (i+1) and elevating the order by 1.
 *  Some reslicing needs to be applied to integrate along the x-axis.
 */
std::vector<double> BasePolynomial::integrateCoeffs(std::vector<double> const& coefficients){
	std::vector<double> newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
	// pushing a zero elevates the order by 1
	newCoefficients.push_back(0.0);
	for(unsigned int i=0; i<coefficients.size(); i++) {
		newCoefficients.push_back(coefficients[i]/(i+1.));
	}
	return newCoefficients;
}
std::vector< std::vector<double> > BasePolynomial::integrateCoeffs(std::vector< std::vector<double> > const& coefficients, bool axis){
	std::vector< std::vector<double> > newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));

	if (axis==true){
		std::vector< std::vector<double> > tmpCoefficients;
		tmpCoefficients = transpose(coefficients);
		unsigned int sizeY = tmpCoefficients.size();
		for(unsigned int i=0; i<sizeY; i++) {
			newCoefficients.push_back(integrateCoeffs(tmpCoefficients[i]));
		}
		return transpose(newCoefficients);
	} else if (axis==false){
		for(unsigned int i=0; i<sizeX; i++) {
			newCoefficients.push_back(integrateCoeffs(coefficients[i]));
		}
		return newCoefficients;
	} else {
		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for integration. %d is not a valid input. ",axis));
	}
	return newCoefficients;
}


/** Deriving coefficients for polynomials is done by multiplying the
 *  original coefficients with i and lowering the order by 1.
 */
std::vector<double> BasePolynomial::deriveCoeffs(std::vector<double> const& coefficients){
	std::vector<double> newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
	// skipping the first element lowers the order
	for(unsigned int i=1; i<coefficients.size(); i++) {
		newCoefficients.push_back(coefficients[i]*i);
	}
	return newCoefficients;
}
std::vector< std::vector<double> > BasePolynomial::deriveCoeffs(const std::vector< std::vector<double> > &coefficients, unsigned int axis){
	std::vector< std::vector<double> > newCoefficients;
	unsigned int sizeX = coefficients.size();
	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));

	if (axis==0){
		std::vector< std::vector<double> > tmpCoefficients;
		tmpCoefficients = transpose(coefficients);
		unsigned int sizeY = tmpCoefficients.size();
		for(unsigned int i=0; i<sizeY; i++) {
			newCoefficients.push_back(deriveCoeffs(tmpCoefficients[i]));
		}
		return transpose(newCoefficients);
	} else if (axis==1){
		for(unsigned int i=0; i<sizeX; i++) {
			newCoefficients.push_back(deriveCoeffs(coefficients[i]));
		}
		return newCoefficients;
	} else {
		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for derivation. %d is not a valid input. ",axis));
	}
	return newCoefficients;
}


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
double BasePolynomial::simplePolynomial(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomial(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += coefficients[i] * pow(x,(int)i);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}
double BasePolynomial::simplePolynomial(std::vector<std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomial(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += pow(x,(int)i) * simplePolynomial(coefficients[i], y);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/// Simple integrated polynomial function generator.
/** Base function to produce integrals of n-th order
 *  polynomials based on the length of the coefficient
 *  vector.
 *  Starts with only the first coefficient at x^0 */
///Indefinite integral in x-direction
double BasePolynomial::simplePolynomialInt(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += 1./(i+1.) * coefficients[i] * pow(x,(int)(i+1.));
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in y-direction only
double BasePolynomial::simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running simplePolynomialInt(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0.;
	for(unsigned int i=0; i<coefficients.size();i++) {
		result += pow(x,(int)i) * simplePolynomialInt(coefficients[i], y);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/// Simple integrated polynomial function generator divided by independent variable.
/** Base function to produce integrals of n-th order
 *  polynomials based on the length of the coefficient
 *  vector.
 *  Starts with only the first coefficient at x^0 */
double BasePolynomial::simpleFracInt(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << x << "): ";
	}
	double result = coefficients[0] * log(x);
	if (coefficients.size() > 1) {
		for (unsigned int i=1; i<coefficients.size(); i++){
			result += 1./(i) * coefficients[i] * pow(x,(int)(i));
		}
	}
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

double BasePolynomial::simpleFracInt(std::vector< std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running      simpleFracInt(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for (unsigned int i=0; i<coefficients.size(); i++){
		result += pow(x,(int)i) * polyfracint(coefficients[i],y);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/** Simple integrated centred(!) polynomial function generator divided by independent variable.
 *  We need to rewrite some of the functions in order to
 *  use central fit. Having a central temperature xbase
 *  allows for a better fit, but requires a different
 *  formulation of the fracInt function group. Other
 *  functions are not affected.
 *  Starts with only the first coefficient at x^0 */

///Helper functions to calculate binomial coefficients: http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C.2B.2B
//double BasePolynomial::factorial(double nValue){
//   double result = nValue;
//   double result_next;
//   double pc = nValue;
//   do {
//	   result_next = result*(pc-1);
//	   result = result_next;
//	   pc--;
//   } while(pc>2);
//   nValue = result;
//   return nValue;
//}
//double BasePolynomial::factorial(double nValue){
//	if (nValue == 0) return (1);
//	else return (nValue * factorial(nValue - 1));
//}
double BasePolynomial::factorial(double nValue){
    double value = 1;
    for(int i = 2; i <= nValue; i++){
        value = value * i;
    }
    return value;
}

double BasePolynomial::binom(double nValue, double nValue2){
   double result;
   if(nValue2 == 1) return nValue;
   result = (factorial(nValue)) / (factorial(nValue2)*factorial((nValue - nValue2)));
   nValue2 = result;
   return nValue2;
}

///Helper functions to calculate the D vector:
std::vector<double> BasePolynomial::fracIntCentralDvector(int m, double x, double xbase){
	std::vector<double> D;
	double tmp;
	if (m<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",m));
	for (int j=0; j<m; j++){ // loop through row
		tmp = pow(-1.0,j) * log(x) * pow(xbase,(int)j);
		for(int k=0; k<j; k++) { // internal loop for every entry
			tmp += binom(j,k) * pow(-1.0,k) / (j-k) * pow(x,j-k) * pow(xbase,k);
		}
		D.push_back(tmp);
	}
	return D;
}

///Indefinite integral of a centred polynomial divided by its independent variable
double BasePolynomial::fracIntCentral(std::vector<double> const& coefficients, double x, double xbase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral(std::vector, " << x << ", " << xbase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	int m = coefficients.size();
	std::vector<double> D = fracIntCentralDvector(m, x, xbase);
	double result = 0;
	for(int j=0; j<m; j++) {
		result += coefficients[j] * D[j];
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/// Horner function generator implementations
/** Represent polynomials according to Horner's scheme.
 *  This avoids unnecessary multiplication and thus
 *  speeds up calculation.
 */
double BasePolynomial::baseHorner(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running       baseHorner(std::vector, " << x << "): ";
	}
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += coefficients[i];
	}
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

double BasePolynomial::baseHorner(std::vector< std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running       baseHorner(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHorner(coefficients[i], y);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in x-direction
double BasePolynomial::baseHornerInt(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << x << "): ";
	}
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += coefficients[i]/(i+1.);
	}
	result = result * x;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in y-direction only
double BasePolynomial::baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running       baseHornerInt(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerInt(coefficients[i], y);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in x-direction of a polynomial divided by its independent variable
double BasePolynomial::baseHornerFracInt(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = 0;
	if (coefficients.size() > 1) {
		for(int i=coefficients.size()-1; i>=1; i--) {
			result *= x;
			result += coefficients[i]/(i);
		}
		result *= x;
	}
	result += coefficients[0] * log(x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
double BasePolynomial::baseHornerFracInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running      baseHornerFra(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;

	double result = 0;
	for(int i=coefficients.size()-1; i>=0; i--) {
		result *= x;
		result += baseHornerFracInt(coefficients[i], y);
	}

	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}


/** Alternatives
 *  Simple functions that heavily rely on other parts of this file.
 *  We still need to check which combinations yield the best
 *  performance.
 */
///Derivative in x-direction
double BasePolynomial::deriveIn2Steps(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running   deriveIn2Steps(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result =  polyval(deriveCoeffs(coefficients),x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Derivative in terms of x(axis=true) or y(axis=false).
double BasePolynomial::deriveIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double y, bool axis){
	if (this->DEBUG) {
		std::cout << "Running   deriveIn2Steps(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = polyval(deriveCoeffs(coefficients,axis),x,y);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in x-direction
double BasePolynomial::integrateIn2Steps(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result =  polyval(integrateCoeffs(coefficients),x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in terms of x(axis=true) or y(axis=false).
double BasePolynomial::integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double y, bool axis){
	if (this->DEBUG) {
		std::cout << "Running   integrateIn2Steps(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = polyval(integrateCoeffs(coefficients,axis),x,y);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in x-direction of a polynomial divided by its independent variable
double BasePolynomial::fracIntIn2Steps(std::vector<double> const& coefficients, double x){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	double result = coefficients[0] * log(x);
	if (coefficients.size() > 1) {
		std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
		result += polyint(newCoeffs,x);
	}
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
double BasePolynomial::fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double y){
	if (this->DEBUG) {
		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << ", " << y << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(polyfracint(coefficients[i],y));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}

///Indefinite integral in y-direction of a centred polynomial divided by its 2nd independent variable
double BasePolynomial::fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double y, double ybase){
	if (this->DEBUG) {
		std::cout << "Running    fracIntCentral2Steps(std::vector, " << x << ", " << y << ", " << ybase << "): ";
	}
	bool db = this->DEBUG;
	this->DEBUG = false;
	std::vector<double> newCoeffs;
	for (unsigned int i=0; i<coefficients.size(); i++){
		newCoeffs.push_back(fracIntCentral(coefficients[i], y, ybase));
	}
	double result = polyval(newCoeffs,x);
	this->DEBUG = db;
	if (this->DEBUG) {
		std::cout << result << std::endl;
	}
	return result;
}




/** Implements the same public functions as the
 *  but solves the polynomial for the given value
 *  instead of evaluating it.
 *  TODO: This class does not check for bijective
 *        polynomials and is therefore a little
 *        fragile.
 */
PolynomialSolver::PolynomialSolver(){
  this->DEBUG   = false;
  this->macheps = DBL_EPSILON;
  this->tol     = DBL_EPSILON*1e3;
  this->maxiter = 100;
}

/** Everything related to the normal polynomials goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves a one-dimensional polynomial for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current input
virtual double PolynomialSolver::polyval(const std::vector<double> &coefficients, double y) {

	BasePolynomial polynomial = BasePolynomial();
	PolyResidual residual = PolyResidual(coefficients, y);

	std::string errstring;
	double result = -1.0;

	switch (this->uses) {
	case iNewton: ///< Newton solver with derivative and guess value
		result = Newton(residual, this->guess, this->tol, this->maxiter, errstring);
		break;

	case iBrent: ///< Brent solver with bounds
		result = Brent(residual, this->min, this->max, this->macheps, this->tol, this->maxiter, errstring);
		break;

	default:
		throw CoolProp::NotImplementedError("This solver has not been implemented.");
	}
	return result;
}

/// Solves a two-dimensional polynomial for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyval(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Everything related to the integrated polynomials goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves the indefinite integral of a one-dimensional polynomial
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
virtual double PolynomialSolver::polyint(const std::vector<double> &coefficients, double y){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyint(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Everything related to the derived polynomials goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves the derivative of a one-dimensional polynomial
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
virtual double PolynomialSolver::polyder(const std::vector<double> &coefficients, double y){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the derivative of a two-dimensional polynomial along the 2nd axis (y)
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyder(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Everything related to the polynomials divided by one variable goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
virtual double PolynomialSolver::polyfracval(const std::vector<double> &coefficients, double y){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Everything related to the integrated polynomials divided by one variable goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
virtual double PolynomialSolver::polyfracint(const std::vector<double> &coefficients, double y){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
/// @param xbase central x-value for fitted function
virtual double PolynomialSolver::polyfracintcentral(const std::vector<double> &coefficients, double y, double xbase){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
/// @param ybase central y-value for fitted function
virtual double PolynomialSolver::polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Everything related to the derived polynomials divided by one variable goes in this
 *  section, holds all the functions for solving polynomials.
 */
/// Solves the derivative of a one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
virtual double PolynomialSolver::polyfracder(const std::vector<double> &coefficients, double y){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the derivative of a two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
virtual double PolynomialSolver::polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double z){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the derivative of a centred one-dimensional polynomial divided by its independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param y double value that represents the current output
/// @param xbase central x-value for fitted function
virtual double PolynomialSolver::polyfracdercentral(const std::vector<double> &coefficients, double y, double xbase){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}

/// Solves the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param z double value that represents the current output
/// @param ybase central y-value for fitted function
virtual double PolynomialSolver::polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase){
	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
}


/** Set the solvers and updates either the guess values or the
 *  boundaries for the variable to solve for.
 */
/// Sets the guess value for the Newton solver and enables it.
/// @param guess double value that represents the guess value
virtual void PolynomialSolver::setGuess(double guess){
	this->uses  = iNewton;
	this->guess = guess;
	this->min   = -1;
	this->max   = -1;
}
/// Sets the limits for the Brent solver and enables it.
/// @param min double value that represents the lower boundary
/// @param max double value that represents the upper boundary
virtual void PolynomialSolver::setLimits(double min, double max){
	this->uses  = iBrent;
	this->guess = -1;
	this->min   = min;
	this->max   = max;
}



/** Implements the function wrapper interface and can be
 *  used by the solvers. This is only an example and you should
 *  use local redefinitions of the class.
 *  TODO: Make multidimensional
 */
PolyResidual::PolyResidual(){
	this->dim = -1;
}

PolyResidual::PolyResidual(const std::vector<double> &coefficients, double y){
	this->output = y;
	this->firstDim = 0;
	this->coefficients.clear();
	this->coefficients.push_back(coefficients);
	this->dim = i1D;
}

PolyResidual::PolyResidual(const std::vector< std::vector<double> > &coefficients, double x, double z){
	this->output = z;
	this->firstDim = x;
	this->coefficients = coefficients;
	this->dim = i2D;
}

virtual double PolyResidual::call(double x){
	//throw CoolProp::NotImplementedError("Please redefine your classes locally.");
	double polyRes = -1;
	if (this->dim==i1D) {
		polyRes = this->poly.polyval(this->coefficients[0], x);
	} else if (this->dim==i2D) {
		polyRes = this->poly.polyval(this->coefficients, this->firstDim, x);
	} else {
		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
	}
	return polyRes - this->output;
}

virtual double PolyResidual::deriv(double x){
//	throw CoolProp::NotImplementedError("Please redefine your classes locally.");
	double polyRes = -1;
	if (this->dim==i1D) {
		polyRes = this->poly.polyder(this->coefficients[0], x);
	} else if (this->dim==i2D) {
		polyRes = this->poly.polyder(this->coefficients, this->firstDim, x);
	} else {
		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
	}
	return polyRes;
}

virtual double PolyIntResidual::call(double x){
	//throw CoolProp::NotImplementedError("Please redefine your classes locally.");
	double polyRes = -1;
	if (this->dim==i1D) {
		polyRes = this->poly.polyint(this->coefficients[0], x);
	} else if (this->dim==i2D) {
		polyRes = this->poly.polyint(this->coefficients, this->firstDim, x);
	} else {
		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
	}
	return polyRes - this->output;
}

virtual double PolyIntResidual::deriv(double x){
//	throw CoolProp::NotImplementedError("Please redefine your classes locally.");
	double polyRes = -1;
	if (this->dim==i1D) {
		polyRes = this->poly.polyval(this->coefficients[0], x);
	} else if (this->dim==i2D) {
		polyRes = this->poly.polyval(this->coefficients, this->firstDim, x);
	} else {
		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
	}
	return polyRes;
}

virtual double PolyDerResidual::call(double x){
	//throw CoolProp::NotImplementedError("Please redefine your classes locally.");
	double polyRes = -1;
	if (this->dim==i1D) {
		polyRes = this->poly.polyder(this->coefficients[0], x);
	} else if (this->dim==i2D) {
		polyRes = this->poly.polyder(this->coefficients, this->firstDim, x);
	} else {
		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
	}
	return polyRes - this->output;
}

virtual double PolyDerResidual::deriv(double x){
	throw CoolProp::NotImplementedError("2nd derivative of a polynomial is not defined.");
}


/** Here we define the functions that should be to evaluate exponential
 *  functions. Not really polynomials, I know...
 */

BaseExponential::BaseExponential(){
  this->DEBUG = false;
//  this->poly = new BaseExponential();
}
//
//BaseExponential::~BaseExponential(){
//  delete this->poly;
//}

/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input
/// @param n int value that determines the kind of exponential function
double BaseExponential::expval(const std::vector<double> &coefficients, double x, int n){
	double result = 0.;
	if (n==1) {
		this->poly.checkCoefficients(coefficients,3);
		result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
	} else if (n==2) {
		result = exp(this->poly.polyval(coefficients, x));
	} else {
		throw ValueError(format("There is no function defined for this input (%d). ",n));
	}
	return result;
}

/// Evaluates an exponential function for the given coefficients
/// @param coefficients vector containing the ordered coefficients
/// @param x double value that represents the current input in the 1st dimension
/// @param y double value that represents the current input in the 2nd dimension
/// @param n int value that determines the kind of exponential function
double BaseExponential::expval(const std::vector< std::vector<double> > &coefficients, double x, double y, int n){
	double result = 0.;
	if (n==2) {
		result = exp(this->poly.polyval(coefficients, x, y));
	} else {
		throw ValueError(format("There is no function defined for this input (%d). ",n));
	}
	return result;
}



//// Define the residual to be driven to zero
//    class solver_resid : public FuncWrapper1D
//    {
//    public:
//        int other;
//        double T, value, r, eos, rhomolar;
//        HelmholtzEOSMixtureBackend *HEOS;
//
//        solver_resid(HelmholtzEOSMixtureBackend *HEOS, double T, double value, int other){
//            this->HEOS = HEOS; this->T = T; this->value = value; this->other = other;
//        };
//        double call(double rhomolar){
//            this->rhomolar = rhomolar;
//            switch(other)
//            {
//            case iSmolar:
//                eos = HEOS->calc_smolar_nocache(T,rhomolar); break;
//            case iHmolar:
//                eos = HEOS->calc_hmolar_nocache(T,rhomolar); break;
//            case iUmolar:
//                eos = HEOS->calc_umolar_nocache(T,rhomolar); break;
//            default:
//                throw ValueError(format("Input not supported"));
//            }
//
//            r = eos-value;
//            return r;
//        };
//    };
//    solver_resid resid(this, T, value, other);
//    std::string errstring;



}

//int main() {
//
//	SimpleIncompressible* liquid = new DowthermQClass();
//	double AT      =  150.0 + 273.15;
//	double Ap      =  3e5;
//    liquid->testInputs(AT,Ap);
//
//
//	SecCoolSolution* obj = new MethanolSolution();
//    double x      =   0.25;
//    double T      =   5.0 + 273.15;
//    double p      =   3e5;
//
//	obj->testInputs(T+00,p,x);
//	obj->testInputs(T+05,p,x);
//	obj->testInputs(T+10,p,x);
//	obj->testInputs(T+15,p,x);
//
//
//}

int main() {

	std::vector<double> cHeat;
    cHeat.clear();
    cHeat.push_back(999.729);
    cHeat.push_back(2.87576);

	CoolProp::BasePolynomial base = CoolProp::BasePolynomial();

	double result = base.polyval(cHeat,273.15+50);

	printf("From object:      h = %3.3f \t kg/m3    \n",result);

}
