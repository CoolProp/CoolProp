
#include "IncompressibleFluid.h"
#include "math.h"
#include "MatrixMath.h"
#include "PolyMath.h"

namespace CoolProp {



/// A thermophysical property provider for critical and reducing values as well as derivatives of Helmholtz energy
/**
This fluid instance is populated using an entry from a JSON file
*/
//IncompressibleFluid::IncompressibleFluid();

void IncompressibleFluid::set_reference_state(double T0, double p0, double x0=0.0, double h0=0.0, double s0=0.0){

	this->rhoref = rho(T0,p0,x0);
	this->pref = p0;

	this->uref = 0.0;
	this->uref = u(T0,p0,x0);//(href - pref/rhoref);

	this->href = h0; // set new reference value
	this->sref = s0; // set new reference value
	this->href = h(T0,p0,x0); // adjust offset to fit to equations
	this->sref = s(T0,p0,x0); // adjust offset to fit to equations
}

void IncompressibleFluid::validate(){throw NotImplementedError("TODO");}

/// Base function that handles the custom data type
double IncompressibleFluid::baseFunction(IncompressibleData data, double x_in, double y_in=0.0){
	switch (data.type) {
		case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
			//throw NotImplementedError("Here you should implement the polynomial.");
			return poly.evaluate(data.coeffs, x_in, y_in, 0, 0, Tbase, xbase);
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
			//throw NotImplementedError("Here you should implement the exponential.");
			poly.checkCoefficients(data.coeffs, 3,1);
			return exp( data.coeffs(0,0) / ( (x_in-Tbase)+data.coeffs(1,0) ) - data.coeffs(2,0) );
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
			//throw NotImplementedError("Here you should implement the exponential polynomial.");
			return exp(poly.evaluate(data.coeffs, x_in, y_in, 0, 0, Tbase, xbase));
			break;
		case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
			throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,data.type));
			break;
		default:
			throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,data.type));
			break;
	}
	return -_HUGE;
}

/// Entropy as a function of temperature, pressure and composition.
double IncompressibleFluid::s   (double T, double p, double x=0.0){
	IncompressibleData data = specific_heat;
	switch (data.type) {
		case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
			//throw NotImplementedError("Here you should implement the polynomial.");
			return poly.integral(data.coeffs, T, x, 0, -1, 0, Tbase, xbase) - sref;
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
			throw NotImplementedError(format("%s (%d): There is no automatic integration of the exponential function.",__FILE__,__LINE__));
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
			throw NotImplementedError(format("%s (%d): There is no automatic integration of the exponential polynomial function.",__FILE__,__LINE__));
			break;
		case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
			throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,data.type));
			break;
		default:
			throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,data.type));
			break;
	}
	return -_HUGE;
}

/// Internal energy as a function of temperature, pressure and composition.
double IncompressibleFluid::u   (double T, double p, double x=0.0){
	IncompressibleData data = specific_heat;
	switch (data.type) {
		case IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL:
			//throw NotImplementedError("Here you should implement the polynomial.");
			return poly.integral(data.coeffs, T, x, 0, 0, 0, Tbase, xbase) - uref;
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPONENTIAL:
			throw NotImplementedError(format("%s (%d): There is no automatic integration of the exponential function.",__FILE__,__LINE__));
			break;
		case IncompressibleData::INCOMPRESSIBLE_EXPPOLYNOMIAL:
			throw NotImplementedError(format("%s (%d): There is no automatic integration of the exponential polynomial function.",__FILE__,__LINE__));
			break;
		case IncompressibleData::INCOMPRESSIBLE_NOT_SET:
			throw ValueError(format("%s (%d): The function type is not specified (\"[%d]\"), are you sure the coefficients have been set?",__FILE__,__LINE__,data.type));
			break;
		default:
			throw ValueError(format("%s (%d): Your function type \"[%d]\" is unknown.",__FILE__,__LINE__,data.type));
			break;
	}
	return -_HUGE;
}

/// Freezing temperature as a function of pressure and composition.
double Tfreeze(       double p, double x){throw NotImplementedError("TODO");}
/// Define freezing point calculations
//double Tfreeze(double p, double x){
//	IncompressibleClass::checkCoefficients(cTfreeze,5);
//	std::vector<double> tmpVector(cTfreeze.begin()+1,cTfreeze.end());
//	return polyval(tmpVector, x*100.0-cTfreeze[0])+273.15;
//}
/// Conversion from volume-based to mass-based composition.
double V2M (double T,           double y){throw NotImplementedError("TODO");}
/// Conversion from mass-based to mole-based composition.
double M2M (double T,           double x){throw NotImplementedError("TODO");}

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
bool IncompressibleFluid::checkT(double T, double p, double x=0.0){throw NotImplementedError("TODO");}

/// Check validity of pressure input.
/** Compares the given pressure p to the saturation pressure at
 *  temperature T and throws and exception if p is lower than
 *  the saturation conditions.
 *  The default value for psat is -1 yielding true if psat
 *  is not redefined in the subclass.
 *  */
bool IncompressibleFluid::checkP(double T, double p, double x=0.0){throw NotImplementedError("TODO");}

/// Check validity of composition input.
/** Compares the given composition x to a stored minimum and
 *  maximum value. Enforces the redefinition of xmin and
 *  xmax since the default values cause an error. */
bool IncompressibleFluid::checkX(double x){throw NotImplementedError("TODO");}

/// Check validity of temperature, pressure and composition input.
bool IncompressibleFluid::checkTPX(double T, double p, double x=0.0){throw NotImplementedError("TODO");}

} /* namespace CoolProp */



// Testing still needs to be enhanced.
/* Below, I try to carry out some basic tests for both 2D and 1D
 * polynomials as well as the exponential functions for vapour
 * pressure etc.
 */
#ifdef ENABLE_CATCH
#include <math.h>
#include <iostream>
#include "catch.hpp"

TEST_CASE("Internal consistency checks and example use cases for the incompressible fluids","[IncompressibleFluids]")
{
	bool PRINT = false;
	std::string tmpStr;
	std::vector<double> tmpVector;
	std::vector< std::vector<double> > tmpMatrix;


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

		//XLT.set_reference_state(25+273.15, 1.01325e5, 0.0, 0.0, 0.0);
		double Tref = 25+273.15;
		double pref = 0.0;
		XLT.set_reference_state(Tref, pref, 0.0, 0.0, 0.0);

		/// A function to check coefficients and equation types.
		//XLT.validate();


		// Prepare the results and compare them to the calculated values
		double acc = 0.0001;
		double T = 273.15+50;
		double p = 10e5;
		double val = 0;
		double res = 0;

		// Compare density
		val = 824.4615702148608;
		res = XLT.rho(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		// Compare cp
		val = 1834.7455527670554;
		res = XLT.c(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		// Compare s
		val = 145.59157247249246;
		res = XLT.s(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		val = 0.0;
		res = XLT.s(Tref,pref);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( val==res );
		}

		// Compare u
		val = 45212.407309106304;
		res = XLT.u(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		val = 0.0;
		res = XLT.u(Tref,pref);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( val==res );
		}

		// Compare h
		val = 46425.32011926845;
		res = XLT.h(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		val = 0.0;
		res = XLT.h(Tref,pref);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( val==res );
		}

		// Compare v
		val = 0.0008931435169681835;
		res = XLT.visc(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}

		// Compare l
		val = 0.10410086156049088;
		res = XLT.cond(T,p);
		{
		CAPTURE(T);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(val,res,acc) );
		}




	}


	/// Test case for Methanol from SecCool


//	name = std::string("SecCoolSolution");
//	description = std::string("Test Methanol SecCool");
//	reference = std::string("Test");
//
//	Tmin     = -50 + 273.15;
//	Tmax     =  20 + 273.15;
//	TminPsat = Tmax;
//
//	xmin     = 0.0;
//	xmax     = 0.5;
//
//	Tbase    = -4.48 + 273.15;
//	xbase    = 31.57 / 100.0;
//
//	tmpVector.clear();
//	tmpVector.push_back( 960.24665800);
//	tmpVector.push_back(-1.2903839100);
//	tmpVector.push_back(-0.0161042520);
//	tmpVector.push_back(-0.0001969888);
//	tmpVector.push_back( 1.131559E-05);
//	tmpVector.push_back( 9.181999E-08);
//	tmpVector.push_back(-0.4020348270);
//	tmpVector.push_back(-0.0162463989);
//	tmpVector.push_back( 0.0001623301);
//	tmpVector.push_back( 4.367343E-06);
//	tmpVector.push_back( 1.199000E-08);
//	tmpVector.push_back(-0.0025204776);
//	tmpVector.push_back( 0.0001101514);
//	tmpVector.push_back(-2.320217E-07);
//	tmpVector.push_back( 7.794999E-08);
//	tmpVector.push_back( 9.937483E-06);
//	tmpVector.push_back(-1.346886E-06);
//	tmpVector.push_back( 4.141999E-08);
//	cRho.clear();
//	cRho = makeMatrix(tmpVector);
//
//
//
//
//	std::vector<double> tmpVector;
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[0]);
//	tmpVector.push_back(coefficients[6]);
//	tmpVector.push_back(coefficients[11]);
//	tmpVector.push_back(coefficients[15]);
//	matrix.push_back(tmpVector);
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[1]);
//	tmpVector.push_back(coefficients[7]);
//	tmpVector.push_back(coefficients[12]);
//	tmpVector.push_back(coefficients[16]);
//	matrix.push_back(tmpVector);
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[2]);
//	tmpVector.push_back(coefficients[8]);
//	tmpVector.push_back(coefficients[13]);
//	tmpVector.push_back(coefficients[17]);
//	matrix.push_back(tmpVector);
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[3]);
//	tmpVector.push_back(coefficients[9]);
//	tmpVector.push_back(coefficients[14]);
//	tmpVector.push_back(0.0);
//	matrix.push_back(tmpVector);
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[4]);
//	tmpVector.push_back(coefficients[10]);
//	tmpVector.push_back(0.0);
//	tmpVector.push_back(0.0);
//	matrix.push_back(tmpVector);
//
//	tmpVector.clear();
//	tmpVector.push_back(coefficients[5]);
//	tmpVector.push_back(0.0);
//	tmpVector.push_back(0.0);
//	tmpVector.push_back(0.0);
//	matrix.push_back(tmpVector);
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	cTfreeze.clear();
//	cTfreeze.push_back( 27.755555600); // reference concentration in per cent
//	cTfreeze.push_back(-22.973221700);
//	cTfreeze.push_back(-1.1040507200);
//	cTfreeze.push_back(-0.0120762281);
//	cTfreeze.push_back(-9.343458E-05);
//
//
//
//
//
//
//
//
//	double deltaT = 0.1;
//	double Tmin   = 273.15- 50;
//	double Tmax   = 273.15+250;
//	double Tinc   = 200;
//
//	std::vector<std::vector<double> > cHeat2D;
//	cHeat2D.push_back(cHeat);
//	cHeat2D.push_back(cHeat);
//	cHeat2D.push_back(cHeat);
//
//	Eigen::MatrixXd matrix2D = CoolProp::vec_to_eigen(cHeat2D);
//
//	Eigen::MatrixXd matrix2Dtmp;
//	std::vector<std::vector<double> > vec2Dtmp;
//
//	SECTION("Coefficient parsing") {
//		CoolProp::Polynomial2D poly;
//		CHECK_THROWS(poly.checkCoefficients(matrix2D,4,5));
//		CHECK( poly.checkCoefficients(matrix2D,3,4) );
//	}


}


#endif /* ENABLE_CATCH */
