
#if defined(_MSC_VER)
#    define _CRTDBG_MAP_ALLOC
#    ifndef _CRT_SECURE_NO_WARNINGS
#        define _CRT_SECURE_NO_WARNINGS
#    endif
#    include <crtdbg.h>
#    include <sys/stat.h>
#else
#    include <sys/stat.h>
#endif

#include <string>
#include "CoolProp/CoolProp.h"

#include "IncompressibleBackend.h"
#include "CoolProp/fluids/IncompressibleFluid.h"
#include "IncompressibleLibrary.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/numerics/Solvers.h"
#include "CoolProp/numerics/MatrixMath.h"

namespace CoolProp {

IncompressibleBackend::IncompressibleBackend() {
    throw NotImplementedError("Empty constructor is not implemented for incompressible fluids");
    //this->fluid = NULL;
}

IncompressibleBackend::IncompressibleBackend(IncompressibleFluid* fluid)
  : fluid(fluid), chebyshev_caloric(IncompressibleFluid::chebyshevCaloricPreferred()) {

    this->set_reference_state();
    if (this->fluid->is_pure()) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
    } else {
        this->set_fractions(std::vector<CoolPropDbl>(1, 0.0));
    }
    //} else if (ValidNumber(this->fluid->xmin)) {
    //	this->set_fractions(std::vector<CoolPropDbl>(1,this->fluid->getxmin()));
    //} else if (ValidNumber(this->fluid->xmax)) {
    //	this->set_fractions(std::vector<CoolPropDbl>(1,this->fluid->getxmax()));
    //}
}

IncompressibleBackend::IncompressibleBackend(const std::string& fluid_name)
  : fluid(&get_incompressible_fluid(fluid_name)), chebyshev_caloric(IncompressibleFluid::chebyshevCaloricPreferred()) {

    this->set_reference_state();
    if (this->fluid->is_pure()) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
    } else {
        this->set_fractions(std::vector<CoolPropDbl>(1, 0.0));
    }
}

IncompressibleBackend::IncompressibleBackend(const std::vector<std::string>& component_names) {
    throw NotImplementedError("Mixture-style constructor is not implemented yet for incompressible fluids");
}

void IncompressibleBackend::update(CoolProp::input_pairs input_pair, double value1, double value2) {
    //if (mass_fractions.empty()){
    //    throw ValueError("mass fractions have not been set");
    //}

    if (get_debug_level() >= 10) {
        //throw ValueError(format("%s (%d): You have to provide a dimension, 0 or 1, for the solver, %d is not valid. ",__FILE__,__LINE__,axis));
        std::cout << format("Incompressible backend: Called update with %d and %f, %f ", input_pair, value1, value2) << '\n';
    }

    clear();

    if (get_debug_level() >= 50) {
        std::cout << format("Incompressible backend: _fractions are %s ", vec_to_string(_fractions).c_str()) << '\n';
    }
    if (_fractions.size() != 1) {
        throw ValueError(format("%s is an incompressible fluid, mass fractions must be set to a vector with ONE entry, not %d.", this->name().c_str(),
                                _fractions.size()));
    }
    if (fluid->is_pure()) {
        this->_fluid_type = FLUID_TYPE_INCOMPRESSIBLE_LIQUID;
        if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Fluid type is  %d ", this->_fluid_type) << '\n';
        if (_fractions[0] != 1.0) {
            throw ValueError(format("%s is a pure fluid. The composition has to be set to a vector with one entry equal to 1.0. %s is not valid.",
                                    this->name().c_str(), vec_to_string(_fractions).c_str()));
        }
    } else {
        this->_fluid_type = FLUID_TYPE_INCOMPRESSIBLE_SOLUTION;
        if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Fluid type is  %d ", this->_fluid_type) << '\n';
        if (!is_in_closed_range(0.0, 1.0, static_cast<double>(_fractions[0]))) {  // rejects NaN too
            throw ValueError(
              format("%s is a solution or brine. Mass fractions must be set to a vector with one entry between 0 and 1. %s is not valid.",
                     this->name().c_str(), vec_to_string(_fractions).c_str()));
        }
    }

    this->_phase = iphase_liquid;
    if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Phase type is  %d ", this->_phase) << '\n';

    switch (input_pair) {
        case PT_INPUTS: {
            _p = value1;
            _T = value2;
            break;
        }
        case DmassP_INPUTS: {
            _p = value2;
            _T = this->DmassP_flash(value1, value2);
            break;
        }
            //        case PUmass_INPUTS: {
            //            _p = value1;
            //            _T = this->PUmass_flash(value1, value2);
            //            break;
            //        }
        case PSmass_INPUTS: {
            _p = value1;
            _T = this->PSmass_flash(value1, value2);
            break;
        }
        case HmassP_INPUTS: {
            _p = value2;
            _T = this->HmassP_flash(value1, value2);
            break;
        }
        case QT_INPUTS: {
            if (value1 != 0) {
                throw ValueError("Incompressible fluids can only handle saturated liquid, Q=0.");
            }
            _T = value2;
            _p = fluid->psat(value2, _fractions[0]);
            break;
        }
        default: {
            throw ValueError(format("This pair of inputs [%s] is not yet supported", get_input_pair_short_desc(input_pair).c_str()));
        }
    }
    if (_p < 0) {
        throw ValueError("p is less than zero");
    }
    if (!ValidNumber(_p)) {
        throw ValueError("p is not a valid number");
    }
    if (_T < 0) {
        throw ValueError("T is less than zero");
    }
    if (!ValidNumber(_T)) {
        throw ValueError("T is not a valid number");
    }
    if (get_debug_level() >= 50)
        std::cout << format("Incompressible backend: Update finished T=%f, p=%f, x=%s ", this->_T, this->_p, vec_to_string(_fractions).c_str())
                  << '\n';
    fluid->checkTPX(_T, _p, _fractions[0]);
}

/// Clear all the cached values
bool IncompressibleBackend::clear() {
    AbstractState::clear();  // Call the base class
                             /// Additional cached elements used for the partial derivatives
    this->_cmass.clear();
    this->_hmass.clear();
    this->_rhomass.clear();
    this->_smass.clear();
    this->_umass.clear();

    this->_drhodTatPx.clear();
    this->_dsdTatPx.clear();
    this->_dhdTatPx.clear();
    this->_dsdTatPxdT.clear();
    this->_dhdTatPxdT.clear();
    this->_dsdpatTx.clear();
    this->_dhdpatTx.clear();
    // Done
    return true;
}

/// Update the reference values and clear the state
void IncompressibleBackend::set_reference_state(double T0, double p0, double x0, double h0, double s0) {
    // Explicit scope: this method is reachable from the ctor
    // (IncompressibleBackend(IncompressibleFluid*) calls
    // set_reference_state()), so the virtual clear() would bypass
    // derived overrides anyway.  Make that explicit.
    IncompressibleBackend::clear();
    /// Reference values, no need to calculate them each time
    this->_hmass_ref.clear();
    this->_smass_ref.clear();
    //
    this->_T_ref = T0;
    this->_p_ref = p0;
    this->_x_ref = x0;
    this->_h_ref = h0;
    this->_s_ref = s0;
}

/// Set the fractions
/**
@param fractions The vector of fractions of the components converted to the correct input
*/
void IncompressibleBackend::set_fractions(const std::vector<CoolPropDbl>& fractions) {
    if (get_debug_level() >= 10)
        std::cout << format("Incompressible backend: Called set_fractions with %s ", vec_to_string(fractions).c_str()) << '\n';
    if (fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the fraction vector and not %d.", fractions.size()));
    if ((this->_fractions.size() != 1) || (this->_fractions[0] != fractions[0])) {  // Change it!
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Updating the fractions triggered a change in reference state %s -> %s",
                                vec_to_string(this->_fractions).c_str(), vec_to_string(fractions).c_str())
                      << '\n';
        this->_fractions = fractions;
        set_reference_state(T_ref(), p_ref(), this->_fractions[0], h_ref(), s_ref());
    }
}

/// Set the mole fractions
/**
@param mole_fractions The vector of mole fractions of the components
*/
void IncompressibleBackend::set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
    if (get_debug_level() >= 10)
        std::cout << format("Incompressible backend: Called set_mole_fractions with %s ", vec_to_string(mole_fractions).c_str()) << '\n';
    if (mole_fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the mole fraction vector and not %d.", mole_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  //( this->_fractions[0]!=1.0 )){
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(mole_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << '\n';
    } else if (fluid->getxid() == IFRAC_MOLE) {
        this->set_fractions(mole_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        tmp_fractions.reserve(mole_fractions.size());
        for (auto mole_fraction : mole_fractions) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromMole(0.0, mole_fraction));
        }
        this->set_fractions(tmp_fractions);
    }
}

/// Set the mass fractions
/**
@param mass_fractions The vector of mass fractions of the components
*/
void IncompressibleBackend::set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
    if (get_debug_level() >= 10)
        std::cout << format("Incompressible backend: Called set_mass_fractions with %s ", vec_to_string(mass_fractions).c_str()) << '\n';
    if (mass_fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the mass fraction vector and not %d.", mass_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  // ( this->_fractions[0]!=1.0 )) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(mass_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << '\n';
    } else if (fluid->getxid() == IFRAC_MASS) {
        this->set_fractions(mass_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        tmp_fractions.reserve(mass_fractions.size());
        for (auto mass_fraction : mass_fractions) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromMass(0.0, mass_fraction));
        }
        this->set_fractions(tmp_fractions);
    }
}

/// Set the volume fractions
/**
@param volu_fractions The vector of volume fractions of the components
*/
void IncompressibleBackend::set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions) {
    if (get_debug_level() >= 10)
        std::cout << format("Incompressible backend: Called set_volu_fractions with %s ", vec_to_string(volu_fractions).c_str()) << '\n';
    if (volu_fractions.size() != 1)
        throw ValueError(
          format("The incompressible backend only supports one entry in the volume fraction vector and not %d.", volu_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  // ( this->_fractions[0]!=1.0 )) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(volu_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << '\n';
    } else if (fluid->getxid() == IFRAC_VOLUME) {
        this->set_fractions(volu_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        tmp_fractions.reserve(volu_fractions.size());
        for (auto volu_fraction : volu_fractions) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromVolume(0.0, volu_fraction));
        }
        this->set_fractions(tmp_fractions);
    }
}

/// Check if the mole fractions have been set, etc.
void IncompressibleBackend::check_status() {
    throw NotImplementedError("Cannot check status for incompressible fluid");
}

/** We have to override some of the functions from the AbstractState.
 *  The incompressibles are only mass-based and do not support conversion
 *  from molar to specific quantities.
 *  We also have a few new chaced variables that we need.
 */
/// Return the mass density in kg/m^3
double IncompressibleBackend::rhomass() {
    if (!_rhomass) _rhomass = calc_rhomass();
    return _rhomass;
}
/// Return the mass enthalpy in J/kg
double IncompressibleBackend::hmass() {
    if (!_hmass) _hmass = calc_hmass();
    return _hmass;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::smass() {
    if (!_smass) _smass = calc_smass();
    return _smass;
}
/// Return the molar internal energy in J/mol
double IncompressibleBackend::umass() {
    if (!_umass) _umass = calc_umass();
    return _umass;
}
/// Return the mass constant pressure specific heat in J/kg/K
double IncompressibleBackend::cmass() {
    if (!_cmass) _cmass = calc_cmass();
    return _cmass;
}

double IncompressibleBackend::drhodTatPx() {
    if (!_drhodTatPx) _drhodTatPx = calc_drhodTatPx(_T, _p, _fractions[0]);
    return _drhodTatPx;
}
double IncompressibleBackend::dsdTatPx() {
    if (!_dsdTatPx) _dsdTatPx = calc_dsdTatPx(_T, _p, _fractions[0]);
    return _dsdTatPx;
}
double IncompressibleBackend::dhdTatPx() {
    if (!_dhdTatPx) _dhdTatPx = calc_dhdTatPx(_T, _p, _fractions[0]);
    return _dhdTatPx;
}
double IncompressibleBackend::dsdTatPxdT() {
    if (!_dsdTatPxdT) _dsdTatPxdT = calc_dsdTatPxdT(_T, _p, _fractions[0]);
    return _dsdTatPxdT;
}
double IncompressibleBackend::dhdTatPxdT() {
    if (!_dhdTatPxdT) _dhdTatPxdT = calc_dhdTatPxdT(_T, _p, _fractions[0]);
    return _dhdTatPxdT;
}
double IncompressibleBackend::dsdpatTx() {
    if (!_dsdpatTx) _dsdpatTx = calc_dsdpatTx(rhomass(), drhodTatPx());
    return _dsdpatTx;
}
double IncompressibleBackend::dhdpatTx() {
    if (!_dhdpatTx) _dhdpatTx = calc_dhdpatTx(_T, rhomass(), drhodTatPx());
    return _dhdpatTx;
}

/// Return the temperature in K
double IncompressibleBackend::T_ref() {
    if (!_T_ref) throw ValueError("Reference temperature is not set");
    return _T_ref;
}
/// Return the pressure in Pa
double IncompressibleBackend::p_ref() {
    if (!_p_ref) throw ValueError("Reference pressure is not set");
    return _p_ref;
}
/// Return the composition
double IncompressibleBackend::x_ref() {
    if (!_x_ref) throw ValueError("Reference composition is not set");
    return _x_ref;
}
/// Return the mass enthalpy in J/kg
double IncompressibleBackend::h_ref() {
    if (!_h_ref) throw ValueError("Reference enthalpy is not set");
    return _h_ref;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::s_ref() {
    if (!_s_ref) throw ValueError("Reference entropy is not set");
    return _s_ref;
}

/// Return the mass enthalpy in J/kg
double IncompressibleBackend::hmass_ref() {
    if (!_hmass_ref) _hmass_ref = raw_calc_hmass(T_ref(), p_ref(), x_ref());
    return _hmass_ref;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::smass_ref() {
    if (!_smass_ref) _smass_ref = raw_calc_smass(T_ref(), p_ref(), x_ref());
    return _smass_ref;
}

CoolPropDbl IncompressibleBackend::DmassP_flash(CoolPropDbl rhomass, CoolPropDbl p) {
    return fluid->T_rho(rhomass, p, _fractions[0], chebyshev_caloric);
}
/// Residual of the HmassP and PSmass flashes, with its derivative in T.
/**
h = Int c dT + p (v - T dv/dT) and s = Int c/T dT - p dv/dT, with v = 1/rho,
so at constant p the exact slopes are
    dh/dT = c - p T d2v/dT2,    ds/dT = c/T - p d2v/dT2,
    d2v/dT2 = 2 rho'^2 / rho^3 - rho'' / rho^2.
The term with d2v/dT2 grows with p, so at high pressure cp alone is a poor
slope. call() keeps rho and rho' at its T for deriv() at the same T.
*/
class IncompressibleBackend::CaloricResidual : public FuncWrapper1DWithDeriv
{
   public:
    CaloricResidual(IncompressibleBackend* backend, bool enthalpy, double p, double x, double target)
      : backend(backend), enthalpy(enthalpy), p(p), x(x), target(target) {}
    double call(double T) override {
        set_density_terms(T);
        return (enthalpy ? backend->raw_hmass_from(T, p, x, rho, drho) : backend->raw_smass_from(T, p, x, rho, drho)) - target;
    }
    double deriv(double T) override {
        const double c = backend->c_at(T, p, x);
        set_density_terms(T);
        const double d2rho = backend->d2rhodT2_at(T, p, x);
        const double d2v = 2 * drho * drho / (rho * rho * rho) - d2rho / (rho * rho);
        return enthalpy ? c - p * T * d2v : c / T - p * d2v;
    }

   private:
    IncompressibleBackend* backend;
    bool enthalpy;
    double p, x, target;
    double T_cached = NAN, rho = NAN, drho = NAN;
    void set_density_terms(double T) {
        if (T == T_cached) return;
        rho = backend->rho_at(T, p, x);
        drho = backend->calc_drhodTatPx(T, p, x);
        T_cached = T;
    }
};

/// First guess for the HmassP and PSmass flashes from the reference state and cp.
/**
h and s are exactly linear in p at fixed T in this model, so their values at
(T_ref, p) follow from the cached reference values without approximation:
    h(T_ref, p) = h(T_ref, p_ref) + (p - p_ref) dh/dp,   and the same for s.
With cp taken constant from T_ref, h - h(T_ref, p) = cp (T - T_ref) and
s - s(T_ref, p) = cp ln(T / T_ref). A guess outside the range is fine:
BoundedNewton then starts from the false-position point of its bracket.
*/
double IncompressibleBackend::flash_first_guess(bool enthalpy, double p, double raw_target) {
    const double x = _fractions[0], T_r = T_ref();
    // rho, rho' and cp at T_ref do not depend on p: keep them until the
    // composition, the reference temperature or the fluid's data change
    GuessCache& g = guess_cache;
    if (g.x != x || g.T_ref != T_r || g.data_version != fluid->getDataVersion()) {
        g.rho = rho_at(T_r, p, x);
        g.drho = calc_drhodTatPx(T_r, p, x);
        g.c = c_at(T_r, p, x);
        g.x = x;
        g.T_ref = T_r;
        g.data_version = fluid->getDataVersion();
    }
    const double rho_r = g.rho, drho_r = g.drho, c_r = g.c;
    if (enthalpy) {
        const double h_r = hmass_ref() + (p - p_ref()) * calc_dhdpatTx(T_r, rho_r, drho_r);
        return T_r + (raw_target - h_r) / c_r;
    }
    const double s_r = smass_ref() + (p - p_ref()) * calc_dsdpatTx(rho_r, drho_r);
    return T_r * std::exp((raw_target - s_r) / c_r);
}

/// Calculate T given pressure and enthalpy
/**
Newton on T inside a sign-change bracket over the liquid range at p (BoundedNewton,
see IncompressibleFluid::solveInLiquidRange),
with the exact slope dh/dT at constant p, starting from the cp-based guess of
flash_first_guess. The bracket keeps every step inside the range and raises
an error when h is outside it, as Brent did.
@param hmass The mass enthalpy in J/kg
@param p The pressure in Pa
@returns T The temperature in K
*/
CoolPropDbl IncompressibleBackend::HmassP_flash(CoolPropDbl hmass, CoolPropDbl p) {
    const double raw_target = hmass - h_ref() + hmass_ref();
    CaloricResidual res(this, true, p, _fractions[0], raw_target);
    const double tol = DBL_EPSILON * 1e3;
    const double T_guess = flash_first_guess(true, p, raw_target);
    return fluid->solveInLiquidRange(p, _fractions[0], [&](double T_lo, double T_hi) { return BoundedNewton(&res, T_guess, T_lo, T_hi, tol, 100); });
}
/// Calculate T given pressure and entropy
/**
Solved like HmassP_flash, with the exact slope ds/dT at constant p.
@param smass The mass entropy in J/kg/K
@param p The pressure in Pa
@returns T The temperature in K
*/
CoolPropDbl IncompressibleBackend::PSmass_flash(CoolPropDbl p, CoolPropDbl smass) {
    const double raw_target = smass - s_ref() + smass_ref();
    CaloricResidual res(this, false, p, _fractions[0], raw_target);
    // s on the classic fits carries a large ln(T) constant, so one ulp of s exceeds what an absolute 2e-13 K tolerance can resolve.
    const double tol = 1e-12 * fluid->getTmax();
    const double T_guess = flash_first_guess(false, p, raw_target);
    return fluid->solveInLiquidRange(p, _fractions[0], [&](double T_lo, double T_hi) { return BoundedNewton(&res, T_guess, T_lo, T_hi, tol, 100); });
}

///// Calculate T given pressure and internal energy
///**
//@param umass The mass internal energy in J/kg
//@param p The pressure in Pa
//@returns T The temperature in K
//*/
//CoolPropDbl IncompressibleBackend::PUmass_flash(CoolPropDbl p, CoolPropDbl umass){
//
//	class PUmass_residual : public FuncWrapper1D {
//	protected:
//		double p,x,u_in;
//		IncompressibleBackend* fluid;
//	protected:
//		PUmass_residual(){};
//	public:
//		PUmass_residual(IncompressibleBackend* fluid, const double &p,  const double &x, const double &u_in){
//			this->p = p;
//			this->x = x;
//			this->u_in = u_in;
//			this->fluid = fluid;
//		}
//		virtual ~PUmass_residual(){};
//		double call(double target){
//			return fluid->u(target,p,x) - u_in;
//		}
//	};
//
//	PUmass_residual res = PUmass_residual(*this, p, _fractions[0], umass);
//
//	std::string errstring;
//	double macheps = DBL_EPSILON;
//	double tol     = DBL_EPSILON*1e3;
//	int    maxiter = 15;
//	double result = Brent(&res, fluid->getTmin(), fluid->getTmax(), macheps, tol, maxiter, errstring);
//	//if (this->do_debug()) std::cout << "Brent solver message: " << errstring << std::endl;
//	return result;
//}

/// We start with the functions that do not need a reference state
CoolPropDbl IncompressibleBackend::calc_melting_line(int param, int given, CoolPropDbl value) {
    if (param == iT && given == iP) {
        return fluid->Tfreeze(value, _fractions[0]);
    } else {
        throw ValueError("For incompressibles, the only valid inputs to calc_melting_line are T(p)");
    }
};
CoolPropDbl IncompressibleBackend::calc_umass() {
    return hmass() - _p / rhomass();
};

/// ... and continue with the ones that depend on reference conditions.
CoolPropDbl IncompressibleBackend::calc_hmass() {
    return h_ref() + raw_calc_hmass(_T, _p, _fractions[0]) - hmass_ref();
};
CoolPropDbl IncompressibleBackend::calc_smass() {
    return s_ref() + raw_calc_smass(_T, _p, _fractions[0]) - smass_ref();
};

/// Functions that can be used with the solver, they miss the reference values!
CoolPropDbl IncompressibleBackend::raw_calc_hmass(double T, double p, double x) {
    return raw_hmass_from(T, p, x, rho_at(T, p, x), calc_drhodTatPx(T, p, x));
};
CoolPropDbl IncompressibleBackend::raw_calc_smass(double T, double p, double x) {
    return raw_smass_from(T, p, x, rho_at(T, p, x), calc_drhodTatPx(T, p, x));
};

/// Calculate the first partial derivative for the desired derivative
CoolPropDbl IncompressibleBackend::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
    // TODO: Can this be accelerated?
    if ((Of == iDmass) && (Wrt == iP)) return 0.0;  // incompressible!
    if ((Of == iDmass) && (Wrt == iHmass) && (Constant == iP)) return drhodTatPx() / dhdTatPx();
    if ((Of == iHmass) && (Wrt == iDmass) && (Constant == iP)) return dhdTatPx() / drhodTatPx();
    if ((Of == iDmass) && (Wrt == iSmass) && (Constant == iP)) return drhodTatPx() / dsdTatPx();
    if ((Of == iSmass) && (Wrt == iDmass) && (Constant == iP)) return dsdTatPx() / drhodTatPx();
    if ((Of == iDmass) && (Wrt == iT) && (Constant == iP)) return drhodTatPx();
    if ((Of == iT) && (Wrt == iDmass) && (Constant == iP)) return 1.0 / drhodTatPx();
    //
    if ((Of == iHmass) && (Wrt == iP) && (Constant == iT)) return dhdpatTx();
    if ((Of == iP) && (Wrt == iHmass) && (Constant == iT)) return 1.0 / dhdpatTx();
    if ((Of == iHmass) && (Wrt == iSmass) && (Constant == iT)) return dhdpatTx() / dsdpatTx();
    if ((Of == iSmass) && (Wrt == iHmass) && (Constant == iT)) return dsdpatTx() / dhdpatTx();
    if ((Of == iHmass) && (Wrt == iT) && (Constant == iP)) return dhdTatPx();
    if ((Of == iT) && (Wrt == iHmass) && (Constant == iP)) return 1.0 / dhdTatPx();
    //
    if ((Of == iSmass) && (Wrt == iP) && (Constant == iT)) return dsdpatTx();
    if ((Of == iP) && (Wrt == iSmass) && (Constant == iT)) return 1.0 / dsdpatTx();
    if ((Of == iSmass) && (Wrt == iT) && (Constant == iP)) return dsdTatPx();
    if ((Of == iT) && (Wrt == iSmass) && (Constant == iP)) return 1.0 / dsdTatPx();
    //if ( (Of==iHmass) && (Wrt==iP) && (Constant==iT) ) return dsdTatPxdT();
    //if ( (Of==iHmass) && (Wrt==iP) && (Constant==iT) ) return dhdTatPxdT();
    throw ValueError("Incompressible fluids only support a limited subset of partial derivatives.");
}

/* Other useful derivatives
 */
/// Partial derivative of entropy with respect to pressure at constant temperature and composition
//  \f[ \left( \frac{\partial s}{\partial p} \right)_T = - \left( \frac{\partial v}{\partial T} \right)_p = \rho^{-2} \left( \frac{\partial \rho}{\partial T} \right)_p \right) \f]
double IncompressibleBackend::calc_dsdpatTx(double rho, double drhodTatPx) {
    return 1 / rho / rho * drhodTatPx;
}
/// Partial derivative of enthalpy with respect to pressure at constant temperature and composition
//  \f[ \left( \frac{\partial h}{\partial p} \right)_T = v - T \left( \frac{\partial v}{\partial T} \right)_p = \rho^{-1} \left( 1 + T \rho^{-1} \left( \frac{\partial \rho}{\partial T} \right)_p \right) \f]
double IncompressibleBackend::calc_dhdpatTx(double T, double rho, double drhodTatPx) {
    return 1 / rho * (1 + T / rho * drhodTatPx);
}

}  // namespace CoolProp

// Testing routines with fixed parameters and known results
/* These functions try to cover as much as possible, but
 * they still need some serious additions.
 */

#ifdef ENABLE_CATCH
#    include <math.h>
#    include <iostream>
#    include <catch2/catch_all.hpp>

#    include "Tests/TestObjects.h"
#    include "CoolProp/Configuration.h"
#    include <memory>
#    include "CoolProp/AbstractState.h"

TEST_CASE("Internal consistency checks and example use cases for the incompressible backend", "[IncompressibleBackend]") {
    CoolProp::IncompressibleFluid fluid = CoolProp::get_incompressible_fluid("MMA2");
    CoolProp::IncompressibleBackend backend = CoolProp::IncompressibleBackend(&fluid);

    SECTION("Test case for Methanol from SecCool") {

        // Some basic functions
        // has to return false
        CHECK(backend.using_mole_fractions() == false);

        //void update(long input_pair, double value1, double value2);

        std::vector<CoolPropDbl> fractions;
        fractions.push_back(0.4);
        CHECK_THROWS(backend.set_mole_fractions(fractions));
        CHECK_NOTHROW(backend.set_mass_fractions(fractions));
        fractions.push_back(0.4);
        CHECK_THROWS(backend.set_mass_fractions(fractions));
        CHECK_THROWS(backend.check_status());

        // Prepare the results and compare them to the calculated values
        double acc = 0.0001;
        double T = 273.15 + 10;
        double p = 10e5;
        double x = 0.25;
        backend.set_mass_fractions(std::vector<CoolPropDbl>(1, x));
        double val = 0;
        double res = 0;

        //CoolProp::set_debug_level(100);

        // Compare density flash
        val = fluid.rho(T, p, x);
        res = backend.DmassP_flash(val, p);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(T, res, acc));
        }

        // call the update function to set internal variables,
        // concentration has been set before.
        CHECK_THROWS(backend.update(CoolProp::DmassT_INPUTS, val, T));  // First with wrong parameters
        CHECK_NOTHROW(backend.update(CoolProp::PT_INPUTS, p, T));       // ... and then with the correct ones.

        // Compare enthalpy flash
        val = backend.hmass();
        res = backend.HmassP_flash(val, p);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(T, res, acc));
        }

        // Compare entropy flash
        val = backend.smass();
        res = backend.PSmass_flash(p, val);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(T, res, acc));
        }

        /// Get the viscosity [Pa-s]
        val = fluid.visc(T, p, x);
        res = backend.calc_viscosity();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        /// Get the thermal conductivity [W/m/K] (based on the temperature and pressure in the state class)
        val = fluid.cond(T, p, x);
        res = backend.calc_conductivity();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        val = fluid.rho(T, p, x);
        res = backend.rhomass();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        //        val = fluid.h(T, p, x);
        //        res = backend.hmass();
        //        {
        //        CAPTURE(T);
        //        CAPTURE(p);
        //        CAPTURE(x);
        //        CAPTURE(val);
        //        CAPTURE(res);
        //        CHECK( check_abs(val,res,acc) );
        //        }
        //
        //        val = fluid.s(T, p, x);
        //        res = backend.smass();
        //        {
        //        CAPTURE(T);
        //        CAPTURE(p);
        //        CAPTURE(x);
        //        CAPTURE(val);
        //        CAPTURE(res);
        //        CHECK( check_abs(val,res,acc) );
        //        }
        // Make sure the result does not change -> reference state...
        val = backend.smass();
        backend.update(CoolProp::PT_INPUTS, p, T);
        res = backend.smass();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        val = fluid.c(T, p, x);
        res = backend.cpmass();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        res = backend.cvmass();
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Compare Tfreeze
        val = fluid.Tfreeze(p, x);      //-20.02+273.15;// 253.1293105454671;
        res = backend.calc_T_freeze();  // -20.02+273.15;
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        // Testing the reference state function
        double Tref = 20 + 273.15 - 20;
        double pref = 101325 * 10;
        double xref = x;
        backend.set_reference_state(Tref, pref, xref, 0.5, 0.5);
        backend.set_mass_fractions(std::vector<CoolPropDbl>(1, x));
        backend.update(CoolProp::PT_INPUTS, pref, Tref);
        val = 0.5;
        res = backend.hmass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        val = 0.5;
        res = backend.smass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        backend.set_reference_state(Tref, pref, xref, 123, 456);
        backend.update(CoolProp::PT_INPUTS, pref, Tref);
        val = 123;
        res = backend.hmass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        val = 456;
        res = backend.smass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        backend.set_reference_state(Tref, pref, xref, 789, 123);
        backend.update(CoolProp::PT_INPUTS, pref, Tref);
        val = 789;
        res = backend.hmass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }

        val = 123;
        res = backend.smass();
        {
            CAPTURE(Tref);
            CAPTURE(pref);
            CAPTURE(val);
            CAPTURE(res);
            CHECK(check_abs(val, res, acc));
        }
    }

    SECTION("Tests for the full implementation using PropsSI") {

        // Prepare the results and compare them to the calculated values
        std::string fluid("ExampleMelinder");
        double acc = 0.0001;
        double T = -5 + 273.15;
        double p = 10e5;
        double x = 0.3;
        double expected = 0;
        double actual = 0;

        // Compare different inputs
        // ... as vector
        expected = 9.6212e+02;
        std::vector<std::string> fluid_Melinder(1, fluid);
        std::vector<std::vector<double>> IO = CoolProp::PropsSImulti(std::vector<std::string>(1, "D"), "T", std::vector<double>(1, T), "P",
                                                                     std::vector<double>(1, p), "INCOMP", fluid_Melinder, std::vector<double>(1, x));
        REQUIRE(!IO.empty());
        actual = IO[0][0];
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            CHECK(check_abs(expected, actual, acc));
        }
        // ... as %
        actual = CoolProp::PropsSI("D", "T", T, "P", p, "INCOMP::" + fluid + format("-%f%%", x * 100.0));
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
        // ... as mass fraction
        actual = CoolProp::PropsSI("D", "T", T, "P", p, "INCOMP::" + fluid + format("[%f]", x));
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string name = "INCOMP::" + fluid + format("[%f]", x);
            CAPTURE(name);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }
        // entropy reference state problems
        //CoolProp::set_debug_level(51);
        expected = CoolProp::PropsSI("S", "T", T, "P", p, "INCOMP::" + fluid + format("[%f]", x));
        actual = CoolProp::PropsSI("S", "T", T, "P", p, "INCOMP::" + fluid + format("[%f]", x));
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(expected);
            CAPTURE(actual);
            std::string name = "INCOMP::" + fluid + format("[%f]", x);
            CAPTURE(name);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(expected, actual, acc));
        }
    }
    SECTION("SecCool example") {
        double acc = 0.0001;
        std::string backend = "INCOMP";
        std::vector<std::string> fluids(1, "ExampleSecCool");
        double T = -5 + 273.15;
        double p = 10e5;
        double x = 0.4;
        std::vector<double> x_vec = std::vector<double>(1, x);
        std::vector<double> T_vec = std::vector<double>(1, T);
        std::vector<double> p_vec = std::vector<double>(1, p);

        // Compare d
        double dexpected = 9.4844e+02;
        std::vector<std::vector<double>> IO =
          CoolProp::PropsSImulti(std::vector<std::string>(1, "D"), "T", T_vec, "P", p_vec, backend, fluids, x_vec);
        REQUIRE(!IO.empty());
        double dactual = IO[0][0];
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(dexpected);
            CAPTURE(dactual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(dexpected, dactual, acc));
        }

        // Compare cp
        double cpexpected = 3.6304e+03;
        double cpactual = CoolProp::PropsSImulti(std::vector<std::string>(1, "C"), "T", T_vec, "P", p_vec, backend, fluids, x_vec)[0][0];
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(x);
            CAPTURE(cpexpected);
            CAPTURE(cpactual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(cpexpected, cpactual, acc));
        }
    }
    SECTION("INCOMP::ExamplePure") {
        // The golden values below were computed from the centered-polynomial
        // fit of this example fluid's tabulated data. The Chebyshev caloric
        // fits (preferred since the *_cheb entries shipped) reproduce the
        // same data with an independently chosen order, so the two fits may
        // legitimately differ at fit-quality level between the data points
        // (~2e-4 relative for cp). The tolerance is fit-level, not
        // bit-level; a wrong-fit-class regression is orders of magnitude
        // larger.
        double acc = 0.001;
        std::string fluid = std::string("INCOMP::ExamplePure");
        double T = +55 + 273.15;
        double p = 10e5;

        // Compare d
        double dexpected = 7.3646e+02;
        double dactual = CoolProp::PropsSI("D", "T", T, "P", p, fluid);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(dexpected);
            CAPTURE(dactual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(dexpected, dactual, acc));
        }

        // Compare cp
        double cpexpected = 2.2580e+03;
        double cpactual = CoolProp::PropsSI("C", "T", T, "P", p, fluid);
        {
            CAPTURE(T);
            CAPTURE(p);
            CAPTURE(cpexpected);
            CAPTURE(cpactual);
            std::string errmsg = CoolProp::get_global_param_string("errstring");
            CAPTURE(errmsg);
            CHECK(check_abs(cpexpected, cpactual, acc));
        }
    }

    //    SECTION("Tests for the hardcoded fluids") {
    //
    //        // Prepare the results and compare them to the calculated values
    //        std::string fluid("INCOMP::LiBr");
    //        double acc = 0.0001;
    //        double T   =  50 + 273.15;
    //        double p   = 10e5;
    //        double x   = 0.3;
    //        double val = 0;
    //        double res = 0;
    //
    //        // Compare different inputs
    //        // ... as vector
    //        val = 9.6212e+02;
    //        res = CoolProp::PropsSI("D","T",T,"P",p,fluid,std::vector<double>(1,x));
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //        // ... as %
    //        res = CoolProp::PropsSI("D","T",T,"P",p,fluid+format("-%f%s",x*100.0,"%"));
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //        // ... as mass fraction
    //        res = CoolProp::PropsSI("D","T",T,"P",p,fluid+format("[%f]",x));
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //
    //
    //        fluid = std::string("INCOMP::ExampleSecCool");
    //        T   = -5  + 273.15;
    //        p   = 10e5;
    //        x   = 0.4;
    //        std::vector<double> x_vec = std::vector<double>(1,x);
    //
    //        // Compare d
    //        val = 9.4844e+02;
    //        res = CoolProp::PropsSI("D","T",T,"P",p,fluid,x_vec);
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //
    //        // Compare cp
    //        val = 3.6304e+03;
    //        res = CoolProp::PropsSI("C","T",T,"P",p,fluid,x_vec);
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //
    //        fluid = std::string("INCOMP::ExamplePure");
    //        T   = +55  + 273.15;
    //        p   = 10e5;
    //
    //        // Compare d
    //        val = 7.3646e+02;
    //        res = CoolProp::PropsSI("D","T",T,"P",p,fluid);
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //
    //        // Compare cp
    //        val = 2.2580e+03;
    //        res = CoolProp::PropsSI("C","T",T,"P",p,fluid);
    //        {
    //        CAPTURE(T);
    //        CAPTURE(p);
    //        CAPTURE(x);
    //        CAPTURE(val);
    //        CAPTURE(res);
    //        CHECK( check_abs(val,res,acc) );
    //        }
    //    }
}

TEST_CASE("Fits with the composition put in give exactly the full evaluation", "[IncompressibleBackend]") {
    // The backend reduces each fit to a series in T once per composition
    // (CompositionFits). Both evaluations run the same sums in the same order,
    // so they agree to the last bit on x86-64. The check allows a few units in
    // the last place, because a compiler may fuse multiply-adds differently in
    // the two inlined copies (aarch64 does by default).
    // A few ulp relative, plus a margin far below any physical resolution for
    // values that are themselves the small difference of larger terms
    const auto same_value = [](double a, double b) { return a == Catch::Approx(b).epsilon(8 * DBL_EPSILON).margin(1e-9); };
    // The fast path must actually be taken: a fit left unreduced would also
    // pass the comparisons below, through the full evaluation.
    for (const bool chebyshev : {false, true}) {
        CAPTURE(chebyshev);
        CoolProp::CompositionFits reduced;
        CoolProp::get_incompressible_fluid("MEG").fixComposition(0.35, chebyshev, reduced);
        CHECK(reduced.rho.kind != CoolProp::FitAtComposition::NONE);
        CHECK(reduced.c.kind != CoolProp::FitAtComposition::NONE);
        CHECK(reduced.int_c.kind != CoolProp::FitAtComposition::NONE);
        CHECK(reduced.int_c_divT.kind != CoolProp::FitAtComposition::NONE);
        CHECK(reduced.visc.kind != CoolProp::FitAtComposition::NONE);
    }
    struct Case
    {
        const char* name;
        double x;
    };
    for (const Case& c : {Case{"MEG", 0.35}, Case{"LiBr", 0.5}, Case{"Water", 1.0}, Case{"DowQ", 1.0}}) {
        for (const bool chebyshev : {false, true}) {
            CAPTURE(c.name, chebyshev);
            CoolProp::IncompressibleFluid& fluid = CoolProp::get_incompressible_fluid(c.name);
            std::unique_ptr<CoolProp::IncompressibleBackend> made;
            {
                ChebyshevCaloricSetting setting(chebyshev);  // read once, by the constructor
                made = std::make_unique<CoolProp::IncompressibleBackend>(&fluid);
            }
            CoolProp::IncompressibleBackend& backend = *made;
            if (!fluid.is_pure()) backend.set_mass_fractions(std::vector<CoolPropDbl>(1, c.x));
            const double x = fluid.is_pure() ? 1.0 : c.x;
            for (double T : {300.0, 320.0, 345.0}) {
                CAPTURE(T);
                const double p = 1e5;
                backend.update(CoolProp::PT_INPUTS, p, T);
                CHECK(same_value(backend.rhomass(), fluid.rho(T, p, x, chebyshev)));
                CHECK(same_value(backend.cmass(), fluid.c(T, p, x, chebyshev)));
                CHECK(same_value(backend.viscosity(), fluid.visc(T, p, x)));
                CHECK(same_value(backend.conductivity(), fluid.cond(T, p, x)));
                const double rho = fluid.rho(T, p, x, chebyshev), drho = fluid.drhodTatPx(T, p, x, chebyshev);
                CHECK(same_value(backend.raw_calc_hmass(T, p, x), fluid.dhdTatPxdT(T, p, x, chebyshev) + p * (1 / rho * (1 + T / rho * drho))));
                CHECK(same_value(backend.raw_calc_smass(T, p, x), fluid.dsdTatPxdT(T, p, x, chebyshev) + p * (1 / rho / rho * drho)));
            }
        }
    }
}

TEST_CASE("A state sees changes made to its fluid through the setters", "[IncompressibleBackend]") {
    // CompositionFits are rebuilt when the fluid's data_version changes, which
    // every fit setter and setTbase / setxbase bump, not only validate().
    CoolProp::IncompressibleFluid fluid = CoolProp::get_incompressible_fluid("MEG");  // a private copy
    std::unique_ptr<CoolProp::IncompressibleBackend> made;
    {
        ChebyshevCaloricSetting setting(false);  // classic fits: their entropy depends on Tbase
        made = std::make_unique<CoolProp::IncompressibleBackend>(&fluid);
    }
    CoolProp::IncompressibleBackend& backend = *made;
    backend.set_mass_fractions(std::vector<CoolPropDbl>(1, 0.35));
    backend.update(CoolProp::PT_INPUTS, 1e5, 300.0);
    CHECK(std::isfinite(backend.smass()));

    SECTION("setTbase drops the entropy integral until validate() runs again") {
        fluid.setTbase(fluid.getTbase() + 10.0);
        backend.update(CoolProp::PT_INPUTS, 1e5, 310.0);
        CHECK_THROWS(backend.raw_calc_smass(310.0, 1e5, 0.35));
    }
    SECTION("New density and conductivity fits are used at once") {
        CoolProp::IncompressibleData density;
        density.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        density.coeffs = Eigen::MatrixXd::Constant(1, 1, 1234.0);
        fluid.setDensity(density);
        CoolProp::IncompressibleData conductivity = density;
        conductivity.coeffs = Eigen::MatrixXd::Constant(1, 1, 0.25);
        fluid.setConductivity(conductivity);
        backend.update(CoolProp::PT_INPUTS, 1e5, 310.0);
        CHECK(backend.rhomass() == 1234.0);
        CHECK(backend.conductivity() == 0.25);
    }
    SECTION("A new viscosity fit is used at once") {
        CoolProp::IncompressibleData viscosity;
        viscosity.type = CoolProp::IncompressibleData::INCOMPRESSIBLE_POLYNOMIAL;
        viscosity.coeffs = Eigen::MatrixXd::Constant(1, 1, 0.5);
        fluid.setViscosity(viscosity);
        backend.update(CoolProp::PT_INPUTS, 1e5, 310.0);
        CHECK(backend.viscosity() == 0.5);
    }
}

TEST_CASE("The flashes search only the liquid range at the given pressure", "[IncompressibleBackend]") {
    CoolProp::IncompressibleFluid& water = CoolProp::get_incompressible_fluid("Water");
    // Water boils at about 373 K at 1 atm, well below Tmax = 473 K
    const double T_boil = water.TmaxAt(101325.0, 0.0);
    CHECK(T_boil == Catch::Approx(373.1).margin(0.5));
    CHECK(water.psat(T_boil, 0.0) == Catch::Approx(101325.0).epsilon(1e-6));
    // At 20 bar it stays liquid up to Tmax, and without a saturation curve
    // (MEG) the limit is always Tmax
    CHECK(water.TmaxAt(2e6, 0.0) == water.getTmax());
    CoolProp::IncompressibleFluid& meg = CoolProp::get_incompressible_fluid("MEG");
    CHECK(meg.TmaxAt(101325.0, 0.35) == meg.getTmax());

    // An enthalpy that belongs to 400 K has no liquid state at 1 atm; the
    // error says that the liquid boils, not only that there is no bracket
    {
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("INCOMP", "Water"));
        AS->update(CoolProp::PT_INPUTS, 2e6, 400.0);
        const double h_400 = AS->hmass(), s_400 = AS->smass();  // read before the failing updates clear the state
        CHECK_THROWS_WITH(AS->update(CoolProp::HmassP_INPUTS, h_400, 101325.0), Catch::Matchers::ContainsSubstring("boils at"));
        CHECK_THROWS_WITH(AS->update(CoolProp::PSmass_INPUTS, 101325.0, s_400), Catch::Matchers::ContainsSubstring("boils at"));
    }
    // Too cold is not boiling: an enthalpy below the liquid range keeps the
    // bracket error, also at 1 atm where the liquid range ends at boiling
    {
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("INCOMP", "Water"));
        AS->update(CoolProp::PT_INPUTS, 101325.0, 280.0);
        const double h_cold = AS->hmass() - 1e6;
        CHECK_THROWS_WITH(AS->update(CoolProp::HmassP_INPUTS, h_cold, 101325.0),
                          !Catch::Matchers::ContainsSubstring("boils at") && Catch::Matchers::ContainsSubstring("do not bracket"));
    }

    // A liquid exactly at its boiling point is still inside the bracket
    for (const bool chebyshev : {false, true}) {
        CAPTURE(chebyshev);
        ChebyshevCaloricSetting setting(chebyshev);
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("INCOMP", "Water"));
        for (double T : {300.0, 350.0, 400.0, 450.0}) {
            CAPTURE(T);
            AS->update(CoolProp::QT_INPUTS, 0.0, T);
            const double p = AS->p(), h = AS->hmass(), s = AS->smass();
            AS->update(CoolProp::HmassP_INPUTS, h, p);
            CHECK(std::abs(AS->T() - T) < 1e-9);
            AS->update(CoolProp::PSmass_INPUTS, p, s);
            CHECK(std::abs(AS->T() - T) < 1e-9);
        }
    }
}

TEST_CASE("HmassP and PSmass flashes converge above 2048 K", "[IncompressibleBackend]") {
    // Above 1024 K neighbouring doubles are at least as far apart as the
    // HmassP tolerance of 2.2e-13 K, so a purely absolute tolerance could not
    // be met there (liquid sodium goes up to 2500 K, lead-bismuth to 1100 K).
    // The solver-level test in src/Solvers.cpp forces that case; this one
    // checks real fluid data at those temperatures.
    for (const bool chebyshev : {false, true}) {
        CAPTURE(chebyshev);
        ChebyshevCaloricSetting setting(chebyshev);
        std::shared_ptr<CoolProp::AbstractState> AS(CoolProp::AbstractState::factory("INCOMP", "LiqNa"));
        // Above 2048 K neighbouring doubles are 4.5e-13 apart. At 200 bar sodium
        // is liquid up to about 2360 K; above that the fits are extrapolated
        // and h can fall with T, so the flash must stop at the boiling point
        // (TmaxAt) to keep its bracket.
        const double p = 2e7;
        for (int i = 0; i <= 200; ++i) {
            const double T = 2050.0 + 1.5 * i + 0.123456789;
            CAPTURE(T);
            AS->update(CoolProp::PT_INPUTS, p, T);
            const double h = AS->hmass(), s = AS->smass();
            AS->update(CoolProp::HmassP_INPUTS, h, p);
            CHECK(std::abs(AS->T() - T) < 1e-9);
            AS->update(CoolProp::PSmass_INPUTS, p, s);
            CHECK(std::abs(AS->T() - T) < 1e-9);
        }
    }
}

#endif /* ENABLE_CATCH */
