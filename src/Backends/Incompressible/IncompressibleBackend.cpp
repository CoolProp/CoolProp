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

IncompressibleBackend::IncompressibleBackend(IncompressibleFluid* fluid) {
    this->fluid = fluid;
}

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
/// Set the mass fractions
/**
@param mass_fraction The mass fraction of the component other than water
*/
void IncompressibleBackend::set_mass_fractions(const long double &mass_fraction) {
	this->mass_fractions.clear();
	this->mass_fractions.push_back(mass_fraction);
}

/// Check if the mole fractions have been set, etc.
void IncompressibleBackend::check_status() {
	throw NotImplementedError("Cannot check status for incompressible fluid");
}

/// Calculate T given pressure and density
/**
@param rhomass The mass density in kg/m^3
@param p The pressure in Pa
@returns T The temperature in K
*/
long double IncompressibleBackend::DmassP_flash(long double rhomass, long double p){
	return fluid->T_rho(rhomass, p, mass_fractions[0]);
}
/// Calculate T given pressure and enthalpy
/**
@param hmass The mass enthalpy in J/kg
@param p The pressure in Pa
@returns T The temperature in K
*/
long double IncompressibleBackend::HmassP_flash(long double hmass, long double p){

	class HmassP_residual : public FuncWrapper1D {
	protected:
		double p,x,h_in;
		IncompressibleFluid* fluid;
	protected:
		HmassP_residual();
	public:
		HmassP_residual(IncompressibleFluid* fluid, const double &p,  const double &x, const double &h_in){
			this->p = p;
			this->x = x;
			this->h_in = h_in;
			this->fluid = fluid;
		}
		virtual ~HmassP_residual(){};
		double call(double target){
			return fluid->h(target,p,x) - h_in; //fluid.u(target,p,x)+ p / fluid.rho(target,p,x) - h_in;
		}
		//double deriv(double target);
	};

	//double T_tmp = this->PUmass_flash(p, hmass); // guess value from u=h

	HmassP_residual res = HmassP_residual(fluid, p, mass_fractions[0], hmass);

	std::string errstring;
	double macheps = DBL_EPSILON;
	double tol     = DBL_EPSILON*1e3;
	int    maxiter = 10;
	double result = Brent(res, fluid->getTmin(), fluid->getTmax(), macheps, tol, maxiter, errstring);
	//if (this->do_debug()) std::cout << "Brent solver message: " << errstring << std::endl;
	return result;
}
/// Calculate T given pressure and entropy
/**
@param smass The mass entropy in J/kg/K
@param p The pressure in Pa
@returns T The temperature in K
*/
long double IncompressibleBackend::PSmass_flash(long double p, long double smass){
	return fluid->T_s(smass, p, mass_fractions[0]);
}

/// Calculate T given pressure and internal energy
/**
@param umass The mass internal energy in J/kg
@param p The pressure in Pa
@returns T The temperature in K
*/
long double IncompressibleBackend::PUmass_flash(long double p, long double umass){
	return fluid->T_u(umass, p, mass_fractions[0]);
}


}


// Testing routines with fixed parameters and known results
/* These functions try to cover as much as possible, but
 * they still need some serious additions.
 */

#ifdef ENABLE_CATCH
#include <math.h>
#include <iostream>
#include "catch.hpp"

#include "TestObjects.h"

TEST_CASE("Internal consistency checks and example use cases for the incompressible backend","[IncompressibleBackend]")
{
	CoolProp::IncompressibleFluid fluid = CoolPropTesting::incompressibleFluidObject();
	CoolProp::IncompressibleBackend backend = CoolProp::IncompressibleBackend(&fluid);

	SECTION("Test case for Methanol from SecCool") {

		// Some basic functions
		// has to return false
		CHECK( backend.using_mole_fractions()==false );

		//void update(long input_pair, double value1, double value2);

		std::vector<long double> fractions;
		fractions.push_back(0.4);
		CHECK_THROWS( backend.set_mole_fractions(fractions) );
		CHECK_NOTHROW( backend.set_mass_fractions(fractions) );
		fractions.push_back(0.4);
		CHECK_THROWS( backend.set_mass_fractions(fractions) );
		CHECK_NOTHROW( backend.set_mass_fractions(0.4) );
		CHECK_THROWS( backend.check_status() );



		// Prepare the results and compare them to the calculated values
		double acc = 0.0001;
		double T   = 273.15+10;
		double p   = 10e5;
		double x   = 0.25;
		backend.set_mass_fractions(x);
		double val = 0;
		double res = 0;

		//CoolProp::set_debug_level(100);

		// Compare density flash
		val = fluid.rho(T,p,x);
		//res = backend.DmassP_flash(val, p);
		{
		CAPTURE(T);
		CAPTURE(p);
		CAPTURE(x);
		CAPTURE(val);
		CAPTURE(res);
		CHECK( check_abs(T,backend.DmassP_flash(val, p),acc) );
		}

//
//
//
//
//		/// Calculate T given pressure and density
//		/**
//		@param rhomass The mass density in kg/m^3
//		@param p The pressure in Pa
//		@returns T The temperature in K
//		*/
//		long double DmassP_flash(long double rhomass, long double p);
//		/// Calculate T given pressure and enthalpy
//		/**
//		@param hmass The mass enthalpy in J/kg
//		@param p The pressure in Pa
//		@returns T The temperature in K
//		*/
//		long double HmassP_flash(long double hmass, long double p);
//		/// Calculate T given pressure and entropy
//		/**
//		@param smass The mass entropy in J/kg/K
//		@param p The pressure in Pa
//		@returns T The temperature in K
//		*/
//		long double PSmass_flash(long double p, long double smass);
//
//		/// Calculate T given pressure and internal energy
//		/**
//		@param umass The mass internal energy in J/kg
//		@param p The pressure in Pa
//		@returns T The temperature in K
//		*/
//		long double PUmass_flash(long double p, long double umass);
//
//
//
////		/// Get the viscosity [Pa-s]
////		long double calc_viscosity(void){return fluid->visc(_T, _p, mass_fractions[0]);};
////		/// Get the thermal conductivity [W/m/K] (based on the temperature and pressure in the state class)
////		long double calc_conductivity(void){return fluid->cond(_T, _p, mass_fractions[0]);};
////
////		long double calc_rhomass(void){return fluid->rho(_T, _p, mass_fractions[0]);};
////		long double calc_hmass(void){return fluid->h(_T, _p, mass_fractions[0]);};
////		long double calc_smass(void){return fluid->s(_T, _p, mass_fractions[0]);};
////		long double calc_umass(void){return fluid->u(_T, _p, mass_fractions[0]);};
////		long double calc_cpmass(void){return fluid->cp(_T, _p, mass_fractions[0]);};
////		long double calc_cvmass(void){return fluid->cv(_T, _p, mass_fractions[0]);};
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
//
//
//
//
//		// Compare cp
//		val = 3993.9748117022423;
//		res = CH3OH.c(T,p,x);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		// Compare s
//		val = -206.62646783739274;
//		res = CH3OH.s(T,p,x);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		val = 0.0;
//		res = CH3OH.s(Tref,pref,xref);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( val==res );
//		}
//
//		// Compare u
//		val = -60043.78429641827;
//		res = CH3OH.u(T,p,x);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		val = href - pref/CH3OH.rho(Tref,pref,xref);
//		res = CH3OH.u(Tref,pref,xref);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( val==res );
//		}
//
//		// Compare h
//		val = -59005.67386390795;
//		res = CH3OH.h(T,p,x);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		val = 0.0;
//		res = CH3OH.h(Tref,pref,xref);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( val==res );
//		}
//
//		// Compare v
//		val = 0.0023970245009602097;
//		res = CH3OH.visc(T,p,x)/1e3;
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		// Compare l
//		val = 0.44791148414693727;
//		res = CH3OH.cond(T,p,x);
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}
//
//		// Compare Tfreeze
//		val = -20.02+273.15;// 253.1293105454671;
//		res = CH3OH.Tfreeze(p,x)+273.15;
//		{
//		CAPTURE(T);
//		CAPTURE(p);
//		CAPTURE(x);
//		CAPTURE(val);
//		CAPTURE(res);
//		CHECK( check_abs(val,res,acc) );
//		}


	}


}

#endif /* ENABLE_CATCH */
