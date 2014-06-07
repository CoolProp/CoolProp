
#include "PolyMath.h"

#include "CoolPropTools.h"
#include "Exceptions.h"
#include "MatrixMath.h"

#include <vector>
#include <string>
//#include <sstream>
//#include <numeric>
#include <math.h>
#include <iostream>

#include "Solvers.h"

#include <unsupported/Eigen/Polynomials>

namespace CoolProp{

/// Set the coefficient matrix.
/// @param coefficients matrix containing the ordered coefficients
void Polynomial2D::setCoefficients(const Eigen::MatrixXd &coefficients){
	this->coefficients = coefficients;
}
void Polynomial2D::setCoefficients(const std::vector<std::vector<double> > &coefficients){
	std::size_t c = num_cols(coefficients), r = num_rows(coefficients);
	Eigen::MatrixXd tmp = Eigen::MatrixXd::Constant(r,c,0.0);
	convert(coefficients,tmp);
	this->setCoefficients(tmp);
}

/// Basic checks for coefficient vectors.
/** Starts with only the first coefficient dimension
 *  and checks the matrix size against the parameters rows and columns. */
/// @param rows unsigned integer value that represents the desired degree of the polynomial
bool Polynomial2D::checkCoefficients(const unsigned int rows);
/// @param rows unsigned integer value that represents the desired degree of the polynomial in the 1st dimension
/// @param columns unsigned integer value that represents the desired degree of the polynomial in the 2nd dimension
bool Polynomial2D::checkCoefficients(const unsigned int rows, const unsigned int columns);
/// @param coefficients vector containing the ordered coefficients
/// @param rows unsigned integer value that represents the desired degree of the polynomial
bool Polynomial2D::checkCoefficients(const Eigen::VectorXd &coefficients, const unsigned int rows){
	if (coefficients.rows() == rows){
		return true;
	} else {
		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.rows(),rows));
	}
	return false;
}
/// @param coefficients vector containing the ordered coefficients
/// @param rows unsigned integer value that represents the desired degree of the polynomial
bool Polynomial2D::checkCoefficients(const std::vector<double> &coefficients, const unsigned int rows){
	if (coefficients.size() == rows){
		return true;
	} else {
		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),rows));
	}
	return false;
}
/// @param coefficients matrix containing the ordered coefficients
/// @param rows unsigned integer value that represents the desired degree of the polynomial in the 1st dimension
/// @param columns unsigned integer value that represents the desired degree of the polynomial in the 2nd dimension
bool Polynomial2D::checkCoefficients(const Eigen::MatrixXd &coefficients, const unsigned int rows, const unsigned int columns);
/// @param coefficients matrix containing the ordered coefficients
/// @param rows unsigned integer value that represents the desired degree of the polynomial in the 1st dimension
/// @param columns unsigned integer value that represents the desired degree of the polynomial in the 2nd dimension
bool Polynomial2D::checkCoefficients(const std::vector<std::vector<double> > &coefficients, const unsigned int rows, const unsigned int columns);



///// Basic checks for coefficient vectors.
///** Starts with only the first coefficient dimension
// *  and checks the vector length against parameter n. */
//bool BasePolynomial::checkCoefficients(const std::vector<double> &coefficients, unsigned int n){

//bool BasePolynomial::checkCoefficients(std::vector< std::vector<double> > const& coefficients, unsigned int rows, unsigned int columns){
//	if (coefficients.size() == rows){
//		bool result = true;
//		for(unsigned int i=0; i<rows; i++) {
//			result = result && checkCoefficients(coefficients[i],columns);
//		}
//		return result;
//	} else {
//		throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),rows));
//	}
//	return false;
//}






///// Integration functions
///** Integrating coefficients for polynomials is done by dividing the
// *  original coefficients by (i+1) and elevating the order by 1
// *  through adding a zero as first coefficient.
// *  Some reslicing needs to be applied to integrate along the x-axis.
// *  In the brine/solution equations, reordering of the parameters
// *  avoids this expensive operation. However, it is included for the
// *  sake of completeness.
// */
///// @param coefficients matrix containing the ordered coefficients
///// @param axis unsigned integer value that represents the desired direction of integration
//Eigen::MatrixXd Polynomial2D::integrateCoeffs(const Eigen::MatrixXd &coefficients, unsigned int axis = 1);
///// @param coefficients matrix containing the ordered coefficients
///// @param axis unsigned integer value that represents the desired direction of integration
//std::vector<std::vector<double> > Polynomial2D::integrateCoeffs(const std::vector<std::vector<double> > &coefficients, unsigned int axis = 1);
//
///// Derivative coefficients calculation
///** Deriving coefficients for polynomials is done by multiplying the
// *  original coefficients with i and lowering the order by 1.
// *
// *  It is not really deprecated, but untested and therefore a warning
// *  is issued. Please check this method before you use it.
// */
///// @param coefficients matrix containing the ordered coefficients
///// @param axis unsigned integer value that represents the desired direction of derivation
//Eigen::MatrixXd Polynomial2D::deriveCoeffs(const Eigen::MatrixXd &coefficients, unsigned int axis = 1);
///// @param coefficients matrix containing the ordered coefficients
///// @param axis unsigned integer value that represents the desired direction of derivation
//std::vector<std::vector<double> > Polynomial2D::deriveCoeffs(const std::vector<std::vector<double> > &coefficients, unsigned int axis = 1);
//
///// The core functions to evaluate the polynomial
///** It is here we implement the different special
// *  functions that allow us to specify certain
// *  types of polynomials.
// *  The derivative might bee needed during the
// *  solution process of the solver. It could also
// *  be a protected function...
// */
///// @param x_in double value that represents the current input in the 1st dimension
///// @param y_in double value that represents the current input in the 2nd dimension
//double Polynomial2D::evaluate(const double &x_in, const double &y_in);
///// @param x_in double value that represents the current input in the 1st dimension
///// @param y_in double value that represents the current input in the 2nd dimension
///// @param axis unsigned integer value that represents the axis to solve for (1=x, 2=y)
//double Polynomial2D::derivative(const double &x_in, const double &y_in, unsigned int axis = 1);
///// @param in double value that represents the current input in x (1st dimension) or y (2nd dimension)
///// @param z_in double value that represents the current output in the 3rd dimension
///// @param axis unsigned integer value that represents the axis to solve for (1=x, 2=y)
//double Polynomial2D::solve(const double &in, const double &z_in, unsigned int axis = 1);


};


#ifdef ENABLE_CATCH
#include <math.h>
#include <iostream>
#include "catch.hpp"

TEST_CASE("Internal consistency checks and example use cases for PolyMath.cpp","[PolyMath]")
{
	/// Test case for "SylthermXLT" by "Dow Chemicals"
	std::vector<double> cHeat;
	cHeat.clear();
	cHeat.push_back(+1.1562261074E+03);
	cHeat.push_back(+2.0994549103E+00);
	cHeat.push_back(+7.7175381057E-07);
	cHeat.push_back(-3.7008444051E-20);

	std::vector<std::vector<double> > cHeat2D;
	cHeat2D.push_back(cHeat);
	cHeat2D.push_back(cHeat);

	CoolProp::Polynomial2D poly2D;

	SECTION("Coefficient parsing and setting") {
		poly2D.setCoefficients(cHeat2D);
	}
}


#endif /* ENABLE_CATCH */















/// Constructors for the base class for all Polynomials
//Polynomial1D::Polynomial1D();
//bool Polynomial2D::setCoefficients(const Eigen::MatrixXd &coefficients){
//	this.coefficients = coefficients;
//	return this.coefficients == coefficients;
//}
//bool Polynomial2D::setCoefficients(const std::vector< std::vector<double> > &coefficients){
//	return this->setCoefficients(convert(coefficients));
//}







//namespace CoolProp{
//
//BasePolynomial::BasePolynomial(){
//  this->POLYMATH_DEBUG = false;
//}
//
//
///// Basic checks for coefficient vectors.
///** Starts with only the first coefficient dimension
// *  and checks the vector length against parameter n. */
//bool BasePolynomial::checkCoefficients(const std::vector<double> &coefficients, unsigned int n){
//	if (coefficients.size() == n){
//		return true;
//	} else {
//		throw ValueError(format("The number of coefficients %d does not match with %d. ",coefficients.size(),n));
//	}
//	return false;
//}
//bool BasePolynomial::checkCoefficients(std::vector< std::vector<double> > const& coefficients, unsigned int rows, unsigned int columns){
//	if (coefficients.size() == rows){
//		bool result = true;
//		for(unsigned int i=0; i<rows; i++) {
//			result = result && checkCoefficients(coefficients[i],columns);
//		}
//		return result;
//	} else {
//		throw ValueError(format("The number of rows %d does not match with %d. ",coefficients.size(),rows));
//	}
//	return false;
//}
//
//
///** Integrating coefficients for polynomials is done by dividing the
// *  original coefficients by (i+1) and elevating the order by 1.
// *  Some reslicing needs to be applied to integrate along the x-axis.
// */
//std::vector<double> BasePolynomial::integrateCoeffs(std::vector<double> const& coefficients){
//	std::vector<double> newCoefficients;
//	unsigned int sizeX = coefficients.size();
//	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
//	// pushing a zero elevates the order by 1
//	newCoefficients.push_back(0.0);
//	for(unsigned int i=0; i<coefficients.size(); i++) {
//		newCoefficients.push_back(coefficients[i]/(i+1.));
//	}
//	return newCoefficients;
//}
//std::vector< std::vector<double> > BasePolynomial::integrateCoeffs(std::vector< std::vector<double> > const& coefficients, bool axis){
//	std::vector< std::vector<double> > newCoefficients;
//	unsigned int sizeX = coefficients.size();
//	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
//
//	if (axis==true){
//		std::vector< std::vector<double> > tmpCoefficients;
//		tmpCoefficients = transpose(coefficients);
//		unsigned int sizeY = tmpCoefficients.size();
//		for(unsigned int i=0; i<sizeY; i++) {
//			newCoefficients.push_back(integrateCoeffs(tmpCoefficients[i]));
//		}
//		return transpose(newCoefficients);
//	} else if (axis==false){
//		for(unsigned int i=0; i<sizeX; i++) {
//			newCoefficients.push_back(integrateCoeffs(coefficients[i]));
//		}
//		return newCoefficients;
//	} else {
//		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for integration. %d is not a valid input. ",axis));
//	}
//	return newCoefficients;
//}
//
//
///** Deriving coefficients for polynomials is done by multiplying the
// *  original coefficients with i and lowering the order by 1.
// */
//std::vector<double> BasePolynomial::deriveCoeffs(std::vector<double> const& coefficients){
//	std::vector<double> newCoefficients;
//	unsigned int sizeX = coefficients.size();
//	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
//	// skipping the first element lowers the order
//	for(unsigned int i=1; i<coefficients.size(); i++) {
//		newCoefficients.push_back(coefficients[i]*i);
//	}
//	return newCoefficients;
//}
//std::vector< std::vector<double> > BasePolynomial::deriveCoeffs(const std::vector< std::vector<double> > &coefficients, unsigned int axis){
//	std::vector< std::vector<double> > newCoefficients;
//	unsigned int sizeX = coefficients.size();
//	if (sizeX<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",sizeX));
//
//	if (axis==0){
//		std::vector< std::vector<double> > tmpCoefficients;
//		tmpCoefficients = transpose(coefficients);
//		unsigned int sizeY = tmpCoefficients.size();
//		for(unsigned int i=0; i<sizeY; i++) {
//			newCoefficients.push_back(deriveCoeffs(tmpCoefficients[i]));
//		}
//		return transpose(newCoefficients);
//	} else if (axis==1){
//		for(unsigned int i=0; i<sizeX; i++) {
//			newCoefficients.push_back(deriveCoeffs(coefficients[i]));
//		}
//		return newCoefficients;
//	} else {
//		throw ValueError(format("You can only use x-axis (0) and y-axis (1) for derivation. %d is not a valid input. ",axis));
//	}
//	return newCoefficients;
//}
//
//
///** The core of the polynomial wrappers are the different
// *  implementations that follow below. In case there are
// *  new calculation schemes available, please do not delete
// *  the implementations, but mark them as deprecated.
// *  The old functions are good for debugging since the
// *  structure is easier to read than the backward Horner-scheme
// *  or the recursive Horner-scheme.
// */
//
///// Simple polynomial function generator. <- Deprecated due to poor performance, use Horner-scheme instead
///** Base function to produce n-th order polynomials
// *  based on the length of the coefficient vector.
// *  Starts with only the first coefficient at x^0. */
//double BasePolynomial::simplePolynomial(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running simplePolynomial(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0.;
//	for(unsigned int i=0; i<coefficients.size();i++) {
//		result += coefficients[i] * pow(x,(int)i);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//double BasePolynomial::simplePolynomial(std::vector<std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running simplePolynomial(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0;
//	for(unsigned int i=0; i<coefficients.size();i++) {
//		result += pow(x,(int)i) * simplePolynomial(coefficients[i], y);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
///// Simple integrated polynomial function generator.
///** Base function to produce integrals of n-th order
// *  polynomials based on the length of the coefficient
// *  vector.
// *  Starts with only the first coefficient at x^0 */
/////Indefinite integral in x-direction
//double BasePolynomial::simplePolynomialInt(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running simplePolynomialInt(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0.;
//	for(unsigned int i=0; i<coefficients.size();i++) {
//		result += 1./(i+1.) * coefficients[i] * pow(x,(int)(i+1.));
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in y-direction only
//double BasePolynomial::simplePolynomialInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running simplePolynomialInt(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0.;
//	for(unsigned int i=0; i<coefficients.size();i++) {
//		result += pow(x,(int)i) * simplePolynomialInt(coefficients[i], y);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
///// Simple integrated polynomial function generator divided by independent variable.
///** Base function to produce integrals of n-th order
// *  polynomials based on the length of the coefficient
// *  vector.
// *  Starts with only the first coefficient at x^0 */
//double BasePolynomial::simpleFracInt(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running      simpleFracInt(std::vector, " << x << "): ";
//	}
//	double result = coefficients[0] * log(x);
//	if (coefficients.size() > 1) {
//		for (unsigned int i=1; i<coefficients.size(); i++){
//			result += 1./(i) * coefficients[i] * pow(x,(int)(i));
//		}
//	}
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//double BasePolynomial::simpleFracInt(std::vector< std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running      simpleFracInt(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0;
//	for (unsigned int i=0; i<coefficients.size(); i++){
//		result += pow(x,(int)i) * polyfracint(coefficients[i],y);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
///** Simple integrated centred(!) polynomial function generator divided by independent variable.
// *  We need to rewrite some of the functions in order to
// *  use central fit. Having a central temperature xbase
// *  allows for a better fit, but requires a different
// *  formulation of the fracInt function group. Other
// *  functions are not affected.
// *  Starts with only the first coefficient at x^0 */
//
/////Helper functions to calculate binomial coefficients: http://rosettacode.org/wiki/Evaluate_binomial_coefficients#C.2B.2B
////double BasePolynomial::factorial(double nValue){
////   double result = nValue;
////   double result_next;
////   double pc = nValue;
////   do {
////	   result_next = result*(pc-1);
////	   result = result_next;
////	   pc--;
////   } while(pc>2);
////   nValue = result;
////   return nValue;
////}
////double BasePolynomial::factorial(double nValue){
////	if (nValue == 0) return (1);
////	else return (nValue * factorial(nValue - 1));
////}
//double BasePolynomial::factorial(double nValue){
//    double value = 1;
//    for(int i = 2; i <= nValue; i++){
//        value = value * i;
//    }
//    return value;
//}
//
//double BasePolynomial::binom(double nValue, double nValue2){
//   double result;
//   if(nValue2 == 1) return nValue;
//   result = (factorial(nValue)) / (factorial(nValue2)*factorial((nValue - nValue2)));
//   nValue2 = result;
//   return nValue2;
//}
//
/////Helper functions to calculate the D vector:
//std::vector<double> BasePolynomial::fracIntCentralDvector(int m, double x, double xbase){
//	std::vector<double> D;
//	double tmp;
//	if (m<1) throw ValueError(format("You have to provide coefficients, a vector length of %d is not a valid. ",m));
//	for (int j=0; j<m; j++){ // loop through row
//		tmp = pow(-1.0,j) * log(x) * pow(xbase,(int)j);
//		for(int k=0; k<j; k++) { // internal loop for every entry
//			tmp += binom(j,k) * pow(-1.0,k) / (j-k) * pow(x,j-k) * pow(xbase,k);
//		}
//		D.push_back(tmp);
//	}
//	return D;
//}
//
/////Indefinite integral of a centred polynomial divided by its independent variable
//double BasePolynomial::fracIntCentral(std::vector<double> const& coefficients, double x, double xbase){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running    fracIntCentral(std::vector, " << x << ", " << xbase << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	int m = coefficients.size();
//	std::vector<double> D = fracIntCentralDvector(m, x, xbase);
//	double result = 0;
//	for(int j=0; j<m; j++) {
//		result += coefficients[j] * D[j];
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
///// Horner function generator implementations
///** Represent polynomials according to Horner's scheme.
// *  This avoids unnecessary multiplication and thus
// *  speeds up calculation.
// */
//double BasePolynomial::baseHorner(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running       baseHorner(std::vector, " << x << "): ";
//	}
//	double result = 0;
//	for(int i=coefficients.size()-1; i>=0; i--) {
//		result *= x;
//		result += coefficients[i];
//	}
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//double BasePolynomial::baseHorner(std::vector< std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running       baseHorner(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0;
//	for(int i=coefficients.size()-1; i>=0; i--) {
//		result *= x;
//		result += baseHorner(coefficients[i], y);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in x-direction
//double BasePolynomial::baseHornerInt(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running       baseHornerInt(std::vector, " << x << "): ";
//	}
//	double result = 0;
//	for(int i=coefficients.size()-1; i>=0; i--) {
//		result *= x;
//		result += coefficients[i]/(i+1.);
//	}
//	result = result * x;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in y-direction only
//double BasePolynomial::baseHornerInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running       baseHornerInt(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0;
//	for(int i=coefficients.size()-1; i>=0; i--) {
//		result *= x;
//		result += baseHornerInt(coefficients[i], y);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in x-direction of a polynomial divided by its independent variable
//double BasePolynomial::baseHornerFracInt(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running      baseHornerFra(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = 0;
//	if (coefficients.size() > 1) {
//		for(int i=coefficients.size()-1; i>=1; i--) {
//			result *= x;
//			result += coefficients[i]/(i);
//		}
//		result *= x;
//	}
//	result += coefficients[0] * log(x);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
//double BasePolynomial::baseHornerFracInt(std::vector<std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running      baseHornerFra(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//
//	double result = 0;
//	for(int i=coefficients.size()-1; i>=0; i--) {
//		result *= x;
//		result += baseHornerFracInt(coefficients[i], y);
//	}
//
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
///** Alternatives
// *  Simple functions that heavily rely on other parts of this file.
// *  We still need to check which combinations yield the best
// *  performance.
// */
/////Derivative in x-direction
//double BasePolynomial::deriveIn2Steps(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running   deriveIn2Steps(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result =  polyval(deriveCoeffs(coefficients),x);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Derivative in terms of x(axis=true) or y(axis=false).
//double BasePolynomial::deriveIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double y, bool axis){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running   deriveIn2Steps(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = polyval(deriveCoeffs(coefficients,axis),x,y);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in x-direction
//double BasePolynomial::integrateIn2Steps(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running   integrateIn2Steps(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result =  polyval(integrateCoeffs(coefficients),x);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in terms of x(axis=true) or y(axis=false).
//double BasePolynomial::integrateIn2Steps(std::vector< std::vector<double> > const& coefficients, double x, double y, bool axis){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running   integrateIn2Steps(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = polyval(integrateCoeffs(coefficients,axis),x,y);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in x-direction of a polynomial divided by its independent variable
//double BasePolynomial::fracIntIn2Steps(std::vector<double> const& coefficients, double x){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	double result = coefficients[0] * log(x);
//	if (coefficients.size() > 1) {
//		std::vector<double> newCoeffs(coefficients.begin() + 1, coefficients.end());
//		result += polyint(newCoeffs,x);
//	}
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in y-direction of a polynomial divided by its 2nd independent variable
//double BasePolynomial::fracIntIn2Steps(std::vector<std::vector<double> > const& coefficients, double x, double y){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running    fracIntIn2Steps(std::vector, " << x << ", " << y << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	std::vector<double> newCoeffs;
//	for (unsigned int i=0; i<coefficients.size(); i++){
//		newCoeffs.push_back(polyfracint(coefficients[i],y));
//	}
//	double result = polyval(newCoeffs,x);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
/////Indefinite integral in y-direction of a centred polynomial divided by its 2nd independent variable
//double BasePolynomial::fracIntCentral2Steps(std::vector<std::vector<double> > const& coefficients, double x, double y, double ybase){
//	if (this->POLYMATH_DEBUG) {
//		std::cout << "Running    fracIntCentral2Steps(std::vector, " << x << ", " << y << ", " << ybase << "): ";
//	}
//	bool db = this->POLYMATH_DEBUG;
//	this->POLYMATH_DEBUG = false;
//	std::vector<double> newCoeffs;
//	for (unsigned int i=0; i<coefficients.size(); i++){
//		newCoeffs.push_back(fracIntCentral(coefficients[i], y, ybase));
//	}
//	double result = polyval(newCoeffs,x);
//	this->POLYMATH_DEBUG = db;
//	if (this->POLYMATH_DEBUG) {
//		std::cout << result << std::endl;
//	}
//	return result;
//}
//
//
//
//
///** Implements the function wrapper interface and can be
// *  used by the solvers. This is only an example and you should
// *  use local redefinitions of the class.
// *  TODO: Make multidimensional
// */
//PolyResidual::PolyResidual(){
//	this->dim = -1;
//}
//
//PolyResidual::PolyResidual(const std::vector<double> &coefficients, double y){
//	this->output = y;
//	this->firstDim = 0;
//	this->coefficients.clear();
//	this->coefficients.push_back(coefficients);
//	this->dim = i1D;
//}
//
//PolyResidual::PolyResidual(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	this->output = z;
//	this->firstDim = x;
//	this->coefficients = coefficients;
//	this->dim = i2D;
//}
//
//double PolyResidual::call(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyval(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyval(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes - this->output;
//}
//
//double PolyResidual::deriv(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyder(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyder(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes;
//}
//
//double PolyIntResidual::call(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyint(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyint(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes - this->output;
//}
//
//double PolyIntResidual::deriv(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyval(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyval(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes;
//}
//
//double PolyFracIntResidual::call(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyfracint(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyfracint(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes - this->output;
//}
//
//double PolyFracIntResidual::deriv(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyfracval(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyfracval(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes;
//}
//
//double PolyFracIntCentralResidual::call(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyfracintcentral(this->coefficients[0], x, this->baseVal);
//		break;
//	case i2D:
//		polyRes = this->poly.polyfracintcentral(this->coefficients, this->firstDim, x, this->baseVal);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes - this->output;
//}
//
//double PolyFracIntCentralResidual::deriv(double x){
//	throw CoolProp::NotImplementedError("Derivative of a polynomial frac int is not defined.");
//}
//
//double PolyDerResidual::call(double x){
//	double polyRes = -1;
//	switch (this->dim) {
//	case i1D:
//		polyRes = this->poly.polyder(this->coefficients[0], x);
//		break;
//	case i2D:
//		polyRes = this->poly.polyder(this->coefficients, this->firstDim, x);
//		break;
//	default:
//		throw CoolProp::NotImplementedError("There are only 1D and 2D, a polynomial's live is not 3D.");
//	}
//	return polyRes - this->output;
//}
//
//double PolyDerResidual::deriv(double x){
//	throw CoolProp::NotImplementedError("2nd derivative of a polynomial is not defined.");
//}
//
//
//
//
///** Implements the same public functions as the BasePolynomial
// *  but solves the polynomial for the given value
// *  instead of evaluating it.
// *  TODO: This class does not check for bijective
// *        polynomials and is therefore a little
// *        fragile.
// */
//PolynomialSolver::PolynomialSolver(){
//  this->POLYMATH_DEBUG   = false;
//  this->macheps = DBL_EPSILON;
//  this->tol     = DBL_EPSILON*1e3;
//  this->maxiter = 50;
//}
//
///** Everything related to the normal polynomials goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves a one-dimensional polynomial for the given coefficients
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current input
//double PolynomialSolver::polyval(const std::vector<double> &coefficients, double y) {
//	PolyResidual residual = PolyResidual(coefficients, y);
//	return this->solve(residual);
//}
//
///// Solves a two-dimensional polynomial for the given coefficients
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyval(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	PolyResidual residual = PolyResidual(coefficients, x, z);
//	return this->solve(residual);
//}
//
//
///** Everything related to the integrated polynomials goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves the indefinite integral of a one-dimensional polynomial
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
//double PolynomialSolver::polyint(const std::vector<double> &coefficients, double y){
//	PolyIntResidual residual = PolyIntResidual(coefficients, y);
//	return this->solve(residual);
//}
//
///// Solves the indefinite integral of a two-dimensional polynomial along the 2nd axis (y)
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyint(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	PolyIntResidual residual = PolyIntResidual(coefficients, x, z);
//	return this->solve(residual);
//}
//
//
///** Everything related to the derived polynomials goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves the derivative of a one-dimensional polynomial
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
//double PolynomialSolver::polyder(const std::vector<double> &coefficients, double y){
//	PolyDerResidual residual = PolyDerResidual(coefficients, y);
//	return this->solve(residual);
//}
//
///// Solves the derivative of a two-dimensional polynomial along the 2nd axis (y)
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyder(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	PolyDerResidual residual = PolyDerResidual(coefficients, x, z);
//	return this->solve(residual);
//}
//
//
///** Everything related to the polynomials divided by one variable goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
//double PolynomialSolver::polyfracval(const std::vector<double> &coefficients, double y){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
///// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyfracval(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
//
///** Everything related to the integrated polynomials divided by one variable goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves the indefinite integral of a one-dimensional polynomial divided by its independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
//double PolynomialSolver::polyfracint(const std::vector<double> &coefficients, double y){
//	PolyFracIntResidual residual = PolyFracIntResidual(coefficients, y);
//	return this->solve(residual);
//}
//
///// Solves the indefinite integral of a two-dimensional polynomial divided by its 2nd independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyfracint(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	PolyFracIntResidual residual = PolyFracIntResidual(coefficients, x, z);
//	return this->solve(residual);
//}
//
///// Solves the indefinite integral of a centred one-dimensional polynomial divided by its independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
///// @param xbase central x-value for fitted function
//double PolynomialSolver::polyfracintcentral(const std::vector<double> &coefficients, double y, double xbase){
//	PolyFracIntCentralResidual residual = PolyFracIntCentralResidual(coefficients, y, xbase);
//	return this->solve(residual);
//}
//
///// Solves the indefinite integral of a centred two-dimensional polynomial divided by its 2nd independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
///// @param ybase central y-value for fitted function
//double PolynomialSolver::polyfracintcentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase){
//	PolyFracIntCentralResidual residual = PolyFracIntCentralResidual(coefficients, x, z, ybase);
//	return this->solve(residual);
//}
//
//
///** Everything related to the derived polynomials divided by one variable goes in this
// *  section, holds all the functions for solving polynomials.
// */
///// Solves the derivative of a one-dimensional polynomial divided by its independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
//double PolynomialSolver::polyfracder(const std::vector<double> &coefficients, double y){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
///// Solves the derivative of a two-dimensional polynomial divided by its 2nd independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
//double PolynomialSolver::polyfracder(const std::vector< std::vector<double> > &coefficients, double x, double z){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
///// Solves the derivative of a centred one-dimensional polynomial divided by its independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param y double value that represents the current output
///// @param xbase central x-value for fitted function
//double PolynomialSolver::polyfracdercentral(const std::vector<double> &coefficients, double y, double xbase){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
///// Solves the derivative of a centred two-dimensional polynomial divided by its 2nd independent variable
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param z double value that represents the current output
///// @param ybase central y-value for fitted function
//double PolynomialSolver::polyfracdercentral(const std::vector< std::vector<double> > &coefficients, double x, double z, double ybase){
//	throw CoolProp::NotImplementedError("This solver has not been implemented, yet."); // TODO: Implement function
//}
//
//
///** Set the solvers and updates either the guess values or the
// *  boundaries for the variable to solve for.
// */
///// Sets the guess value for the Newton solver and enables it.
///// @param guess double value that represents the guess value
//void PolynomialSolver::setGuess(double guess){
//	this->uses  = iNewton;
//	this->guess = guess;
//	this->min   = -1;
//	this->max   = -1;
//}
///// Sets the limits for the Brent solver and enables it.
///// @param min double value that represents the lower boundary
///// @param max double value that represents the upper boundary
//void PolynomialSolver::setLimits(double min, double max){
//	this->uses  = iBrent;
//	this->guess = -1;
//	this->min   = min;
//	this->max   = max;
//}
//
///// Solves the equations based on previously defined parameters.
///// @param min double value that represents the lower boundary
///// @param max double value that represents the upper boundary
//double PolynomialSolver::solve(PolyResidual &res){
//	std::string errstring;
//	double result = -1.0;
//	switch (this->uses) {
//	case iNewton: ///< Newton solver with derivative and guess value
//		if (res.is2D()) {
//			throw CoolProp::NotImplementedError("The Newton solver is not suitable for 2D polynomials, yet.");
//		}
//		result = Newton(res, this->guess, this->tol, this->maxiter, errstring);
//		break;
//
//	case iBrent: ///< Brent solver with bounds
//		result = Brent(res, this->min, this->max, this->macheps, this->tol, this->maxiter, errstring);
//		break;
//
//	default:
//		throw CoolProp::NotImplementedError("This solver has not been implemented or you forgot to select a solver...");
//	}
//	return result;
//}
//
//
///** Here we define the functions that should be to evaluate exponential
// *  functions. Not really polynomials, I know...
// */
//
//BaseExponential::BaseExponential(){
//  this->POLYMATH_DEBUG = false;
////  this->poly = new BaseExponential();
//}
////
////BaseExponential::~BaseExponential(){
////  delete this->poly;
////}
//
///// Evaluates an exponential function for the given coefficients
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input
///// @param n int value that determines the kind of exponential function
//double BaseExponential::expval(const std::vector<double> &coefficients, double x, int n){
//	double result = 0.;
//	if (n==1) {
//		this->poly.checkCoefficients(coefficients,3);
//		result = exp(coefficients[0]/(x+coefficients[1]) - coefficients[2]);
//	} else if (n==2) {
//		result = exp(this->poly.polyval(coefficients, x));
//	} else {
//		throw ValueError(format("There is no function defined for this input (%d). ",n));
//	}
//	return result;
//}
//
///// Evaluates an exponential function for the given coefficients
///// @param coefficients vector containing the ordered coefficients
///// @param x double value that represents the current input in the 1st dimension
///// @param y double value that represents the current input in the 2nd dimension
///// @param n int value that determines the kind of exponential function
//double BaseExponential::expval(const std::vector< std::vector<double> > &coefficients, double x, double y, int n){
//	double result = 0.;
//	if (n==2) {
//		result = exp(this->poly.polyval(coefficients, x, y));
//	} else {
//		throw ValueError(format("There is no function defined for this input (%d). ",n));
//	}
//	return result;
//}
//
//
//}
//
//
//#ifdef ENABLE_CATCH
//#include <math.h>
//#include "catch.hpp"
//
//class PolynomialConsistencyFixture {
//public:
//	CoolProp::BasePolynomial poly;
//	CoolProp::PolynomialSolver solver;
////	enum dims {i1D, i2D};
////	double firstDim;
////	int dim;
////	std::vector< std::vector<double> > coefficients;
////
////    void setInputs(const std::vector<double> &coefficients){
////    	this->firstDim = 0;
////    	this->coefficients.clear();
////    	this->coefficients.push_back(coefficients);
////    	this->dim = i1D;
////    }
////
////    void setInputs(const std::vector< std::vector<double> > &coefficients, double x){
////    	this->firstDim = x;
////    	this->coefficients = coefficients;
////    	this->dim = i2D;
////    }
//};
//
//
//TEST_CASE("Internal consistency checks with PolyMath objects","[PolyMath]")
//{
//	CoolProp::BasePolynomial poly;
//	CoolProp::PolynomialSolver solver;
//
//	/// Test case for "SylthermXLT" by "Dow Chemicals"
//	std::vector<double> cHeat;
//	cHeat.clear();
//	cHeat.push_back(+1.1562261074E+03);
//	cHeat.push_back(+2.0994549103E+00);
//	cHeat.push_back(+7.7175381057E-07);
//	cHeat.push_back(-3.7008444051E-20);
//
//	double deltaT = 0.1;
//	double Tmin   = 273.15- 50;
//	double Tmax   = 273.15+250;
//	double Tinc   = 15;
//
//	double val1,val2,val3,val4;
//
//	SECTION("DerFromVal1D") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyval(cHeat, T-deltaT);
//			val2 = poly.polyval(cHeat, T+deltaT);
//			val3 = (val2-val1)/2/deltaT;
//			val4 = poly.polyder(cHeat, T);
//			CAPTURE(T);
//			CAPTURE(val3);
//			CAPTURE(val4);
//			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
//		}
//	}
//	SECTION("ValFromInt1D") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyint(cHeat, T-deltaT);
//			val2 = poly.polyint(cHeat, T+deltaT);
//			val3 = (val2-val1)/2/deltaT;
//			val4 = poly.polyval(cHeat, T);
//			CAPTURE(T);
//			CAPTURE(val3);
//			CAPTURE(val4);
//			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
//		}
//	}
//
//	SECTION("Solve1DNewton") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyval(cHeat, T);
//			solver.setGuess(T+100);
//			val2 = solver.polyval(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyint(cHeat, T);
//			solver.setGuess(T+100);
//			val2 = solver.polyint(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
////			val1 = poly.polyder(cHeat, T);
////			solver.setGuess(T+100);
////			val2 = solver.polyder(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////
////			val1 = poly.polyfracint(cHeat, T);
////			solver.setGuess(T+100);
////			val2 = solver.polyfracint(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
//		}
//	}
//	SECTION("Solve1DBrent") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyval(cHeat, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyval(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyint(cHeat, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyint(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyder(cHeat, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyder(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyfracint(cHeat, T);
//			solver.setLimits(T-100,T+100);
//			val2 = solver.polyfracint(cHeat, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyfracintcentral(cHeat, T, 250.0);
//			solver.setLimits(T-100,T+100);
//			val2 = solver.polyfracintcentral(cHeat, val1, 250.0);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//		}
//	}
//
//	/// Test case for 2D
//	double xDim1 = 0.3;
//	std::vector< std::vector<double> > cHeat2D;
//	cHeat2D.clear();
//	cHeat2D.push_back(cHeat);
//	cHeat2D.push_back(cHeat);
//	cHeat2D.push_back(cHeat);
//
//	//setInputs(cHeat2D, 0.3);
//
//	SECTION("DerFromVal2D") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyval(cHeat2D, xDim1, T-deltaT);
//			val2 = poly.polyval(cHeat2D, xDim1, T+deltaT);
//			val3 = (val2-val1)/2/deltaT;
//			val4 = poly.polyder(cHeat2D, xDim1, T);
//			CAPTURE(T);
//			CAPTURE(val3);
//			CAPTURE(val4);
//			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
//		}
//	}
//
//	SECTION("ValFromInt2D") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyint(cHeat2D, xDim1, T-deltaT);
//			val2 = poly.polyint(cHeat2D, xDim1, T+deltaT);
//			val3 = (val2-val1)/2/deltaT;
//			val4 = poly.polyval(cHeat2D, xDim1, T);
//			CAPTURE(T);
//			CAPTURE(val3);
//			CAPTURE(val4);
//			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
//		}
//	}
//
////	SECTION("Solve2DNewton") {
////		for (double T = Tmin; T<Tmax; T+=Tinc) {
////			val1 = poly.polyval(cHeat2D, xDim1, T);
////			solver.setGuess(T+100);
////			val2 = solver.polyval(cHeat2D, xDim1, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////		}
////	}
//	SECTION("Solve2DBrent") {
//		for (double T = Tmin; T<Tmax; T+=Tinc) {
//			val1 = poly.polyval(cHeat2D, xDim1, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyval(cHeat2D, xDim1, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyint(cHeat2D, xDim1, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyint(cHeat2D, xDim1, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyder(cHeat2D, xDim1, T);
//			solver.setLimits(T-300,T+300);
//			val2 = solver.polyder(cHeat2D, xDim1, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyfracint(cHeat2D, xDim1, T);
//			solver.setLimits(T-100,T+100);
//			val2 = solver.polyfracint(cHeat2D, xDim1, val1);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//
//			val1 = poly.polyfracintcentral(cHeat2D, xDim1, T, 250);
//			solver.setLimits(T-100,T+100);
//			val2 = solver.polyfracintcentral(cHeat2D, xDim1, val1, 250);
//			CAPTURE(T);
//			CAPTURE(val1);
//			CAPTURE(val2);
//			CHECK(fabs(T-val2) < 1e-1);
//		}
//	}
//
//}
//
////TEST_CASE_METHOD(PolynomialConsistencyFixture,"Internal consistency checks","[PolyMath]")
////{
////	/// Test case for "SylthermXLT" by "Dow Chemicals"
////	std::vector<double> cHeat;
////	cHeat.clear();
////	cHeat.push_back(+1.1562261074E+03);
////	cHeat.push_back(+2.0994549103E+00);
////	cHeat.push_back(+7.7175381057E-07);
////	cHeat.push_back(-3.7008444051E-20);
////
////	//setInputs(cHeat);
////	double deltaT = 0.1;
////	double val1,val2,val3,val4;
////
////    SECTION("DerFromVal1D") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T-deltaT);
////			val2 = this->poly.polyval(cHeat, T+deltaT);
////			val3 = (val2-val1)/2/deltaT;
////			val4 = this->poly.polyder(cHeat, T);
////			CAPTURE(T);
////			CAPTURE(val3);
////			CAPTURE(val4);
////			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
////		}
////    }
////
////    SECTION("ValFromInt1D") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyint(cHeat, T-deltaT);
////			val2 = this->poly.polyint(cHeat, T+deltaT);
////			val3 = (val2-val1)/2/deltaT;
////			val4 = this->poly.polyval(cHeat, T);
////			CAPTURE(T);
////			CAPTURE(val3);
////			CAPTURE(val4);
////			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
////		}
////	}
////
////    SECTION("Solve1DNewton") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T);
////			this->solver.setGuess(T+100);
////			val2 = this->solver.polyval(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////		}
////	}
////    SECTION("Solve1DBrent") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T);
////			this->solver.setLimits(T-300,T+300);
////			val2 = this->solver.polyval(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////		}
////	}
////
////    /// Test case for 2D
////    std::vector< std::vector<double> > cHeat2D;
////	cHeat2D.clear();
////	cHeat2D.push_back(cHeat);
////	cHeat2D.push_back(cHeat);
////	cHeat2D.push_back(cHeat);
////
////	//setInputs(cHeat2D, 0.3);
////
////	SECTION("DerFromVal2D") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T-deltaT);
////			val2 = this->poly.polyval(cHeat, T+deltaT);
////			val3 = (val2-val1)/2/deltaT;
////			val4 = this->poly.polyder(cHeat, T);
////			CAPTURE(T);
////			CAPTURE(val3);
////			CAPTURE(val4);
////			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
////		}
////	}
////
////	SECTION("ValFromInt2D") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyint(cHeat, T-deltaT);
////			val2 = this->poly.polyint(cHeat, T+deltaT);
////			val3 = (val2-val1)/2/deltaT;
////			val4 = this->poly.polyval(cHeat, T);
////			CAPTURE(T);
////			CAPTURE(val3);
////			CAPTURE(val4);
////			CHECK( (1.0-fabs(val4/val3)) < 1e-1);
////		}
////	}
////
////    SECTION("Solve2DNewton") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T);
////			this->solver.setGuess(T+100);
////			val2 = this->solver.polyval(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////		}
////	}
////    SECTION("Solve2DBrent") {
////		for (double T = 273.15-50; T<273.15+250; T+=15) {
////			val1 = this->poly.polyval(cHeat, T);
////			this->solver.setLimits(T-300,T+300);
////			val2 = this->solver.polyval(cHeat, val1);
////			CAPTURE(T);
////			CAPTURE(val1);
////			CAPTURE(val2);
////			CHECK(fabs(T-val2) < 1e-1);
////		}
////	}
////
////}
////
////TEST_CASE("Check against hard coded data","[PolyMath]")
////{
////    CHECK(fabs(HumidAir::f_factor(-60+273.15,101325)/(1.00708)-1) < 1e-3);
////    CHECK(fabs(HumidAir::f_factor( 80+273.15,101325)/(1.00573)-1) < 1e-3);
////    CHECK(fabs(HumidAir::f_factor(-60+273.15,10000e3)/(2.23918)-1) < 1e-3);
////    CHECK(fabs(HumidAir::f_factor(300+273.15,10000e3)/(1.04804)-1) < 1e-3);
////}
//
//
//
////int main() {
////
////	Catch::ConfigData &config = session.configData();
////	config.testsOrTags.clear();
////	config.testsOrTags.push_back("[fast]");
////	session.useConfigData(config);
////	return session.run();
////
////}
//
//#endif /* CATCH_ENABLED */
//
//
////int main() {
////
////	std::vector<double> cHeat;
////	cHeat.clear();
////	cHeat.push_back(+1.1562261074E+03);
////	cHeat.push_back(+2.0994549103E+00);
////	cHeat.push_back(+7.7175381057E-07);
////	cHeat.push_back(-3.7008444051E-20);
////
////	CoolProp::BasePolynomial base = CoolProp::BasePolynomial();
////	CoolProp::PolynomialSolver solve = CoolProp::PolynomialSolver();
////
////	double T = 273.15+50;
////
////	double c = base.polyval(cHeat,T);
////	printf("Should be  :      c = %3.3f \t J/kg/K    \n",1834.746);
////	printf("From object:      c = %3.3f \t J/kg/K    \n",c);
////
////	T = 0.0;
////	solve.setGuess(75+273.15);
////	T = solve.polyval(cHeat,c);
////	printf("Should be  :      T = %3.3f \t K    \n",273.15+50.0);
////	printf("From object:      T = %3.3f \t K    \n",T);
////
////	T = 0.0;
////	solve.setLimits(273.15+10,273.15+100);
////	T = solve.polyval(cHeat,c);
////	printf("Should be  :      T = %3.3f \t K    \n",273.15+50.0);
////	printf("From object:      T = %3.3f \t K    \n",T);
////
////}
