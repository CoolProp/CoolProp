
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
//#include "CoolProp.h"

#include "IncompressibleBackend.h"
#include "IncompressibleFluid.h"
#include "IncompressibleLibrary.h"
#include "DataStructures.h"
#include "Solvers.h"
#include "MatrixMath.h"

namespace CoolProp {

IncompressibleBackend::IncompressibleBackend() {
    throw NotImplementedError("Empty constructor is not implemented for incompressible fluids");
    //this->fluid = NULL;
}

IncompressibleBackend::IncompressibleBackend(IncompressibleFluid* fluid) {
    this->fluid = fluid;
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

IncompressibleBackend::IncompressibleBackend(const std::string& fluid_name) {
    this->fluid = &get_incompressible_fluid(fluid_name);
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
        std::cout << format("Incompressible backend: Called update with %d and %f, %f ", input_pair, value1, value2) << std::endl;
    }

    clear();

    if (get_debug_level() >= 50) {
        std::cout << format("Incompressible backend: _fractions are %s ", vec_to_string(_fractions).c_str()) << std::endl;
    }
    if (_fractions.size() != 1) {
        throw ValueError(format("%s is an incompressible fluid, mass fractions must be set to a vector with ONE entry, not %d.", this->name().c_str(),
                                _fractions.size()));
    }
    if (fluid->is_pure()) {
        this->_fluid_type = FLUID_TYPE_INCOMPRESSIBLE_LIQUID;
        if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Fluid type is  %d ", this->_fluid_type) << std::endl;
        if (_fractions[0] != 1.0) {
            throw ValueError(format("%s is a pure fluid. The composition has to be set to a vector with one entry equal to 1.0. %s is not valid.",
                                    this->name().c_str(), vec_to_string(_fractions).c_str()));
        }
    } else {
        this->_fluid_type = FLUID_TYPE_INCOMPRESSIBLE_SOLUTION;
        if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Fluid type is  %d ", this->_fluid_type) << std::endl;
        if ((_fractions[0] < 0.0) || (_fractions[0] > 1.0)) {
            throw ValueError(
              format("%s is a solution or brine. Mass fractions must be set to a vector with one entry between 0 and 1. %s is not valid.",
                     this->name().c_str(), vec_to_string(_fractions).c_str()));
        }
    }

    this->_phase = iphase_liquid;
    if (get_debug_level() >= 50) std::cout << format("Incompressible backend: Phase type is  %d ", this->_phase) << std::endl;

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
                  << std::endl;
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
    this->clear();
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
        std::cout << format("Incompressible backend: Called set_fractions with %s ", vec_to_string(fractions).c_str()) << std::endl;
    if (fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the fraction vector and not %d.", fractions.size()));
    if ((this->_fractions.size() != 1) || (this->_fractions[0] != fractions[0])) {  // Change it!
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Updating the fractions triggered a change in reference state %s -> %s",
                                vec_to_string(this->_fractions).c_str(), vec_to_string(fractions).c_str())
                      << std::endl;
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
        std::cout << format("Incompressible backend: Called set_mole_fractions with %s ", vec_to_string(mole_fractions).c_str()) << std::endl;
    if (mole_fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the mole fraction vector and not %d.", mole_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  //( this->_fractions[0]!=1.0 )){
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(mole_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << std::endl;
    } else if (fluid->getxid() == IFRAC_MOLE) {
        this->set_fractions(mole_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        for (std::size_t i = 0; i < mole_fractions.size(); i++) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromMole(0.0, mole_fractions[i]));
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
        std::cout << format("Incompressible backend: Called set_mass_fractions with %s ", vec_to_string(mass_fractions).c_str()) << std::endl;
    if (mass_fractions.size() != 1)
        throw ValueError(format("The incompressible backend only supports one entry in the mass fraction vector and not %d.", mass_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  // ( this->_fractions[0]!=1.0 )) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(mass_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << std::endl;
    } else if (fluid->getxid() == IFRAC_MASS) {
        this->set_fractions(mass_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        for (std::size_t i = 0; i < mass_fractions.size(); i++) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromMass(0.0, mass_fractions[i]));
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
        std::cout << format("Incompressible backend: Called set_volu_fractions with %s ", vec_to_string(volu_fractions).c_str()) << std::endl;
    if (volu_fractions.size() != 1)
        throw ValueError(
          format("The incompressible backend only supports one entry in the volume fraction vector and not %d.", volu_fractions.size()));
    if ((fluid->getxid() == IFRAC_PURE) && true) {  // ( this->_fractions[0]!=1.0 )) {
        this->set_fractions(std::vector<CoolPropDbl>(1, 1.0));
        if (get_debug_level() >= 20)
            std::cout << format("Incompressible backend: Overwriting fractions for pure fluid with %s -> %s", vec_to_string(volu_fractions).c_str(),
                                vec_to_string(this->_fractions).c_str())
                      << std::endl;
    } else if (fluid->getxid() == IFRAC_VOLUME) {
        this->set_fractions(volu_fractions);
    } else {
        std::vector<CoolPropDbl> tmp_fractions;
        for (std::size_t i = 0; i < volu_fractions.size(); i++) {
            tmp_fractions.push_back((CoolPropDbl)fluid->inputFromVolume(0.0, volu_fractions[i]));
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
double IncompressibleBackend::rhomass(void) {
    if (!_rhomass) _rhomass = calc_rhomass();
    return _rhomass;
}
/// Return the mass enthalpy in J/kg
double IncompressibleBackend::hmass(void) {
    if (!_hmass) _hmass = calc_hmass();
    return _hmass;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::smass(void) {
    if (!_smass) _smass = calc_smass();
    return _smass;
}
/// Return the molar internal energy in J/mol
double IncompressibleBackend::umass(void) {
    if (!_umass) _umass = calc_umass();
    return _umass;
}
/// Return the mass constant pressure specific heat in J/kg/K
double IncompressibleBackend::cmass(void) {
    if (!_cmass) _cmass = calc_cmass();
    return _cmass;
}

double IncompressibleBackend::drhodTatPx(void) {
    if (!_drhodTatPx) _drhodTatPx = calc_drhodTatPx(_T, _p, _fractions[0]);
    return _drhodTatPx;
}
double IncompressibleBackend::dsdTatPx(void) {
    if (!_dsdTatPx) _dsdTatPx = calc_dsdTatPx(_T, _p, _fractions[0]);
    return _dsdTatPx;
}
double IncompressibleBackend::dhdTatPx(void) {
    if (!_dhdTatPx) _dhdTatPx = calc_dhdTatPx(_T, _p, _fractions[0]);
    return _dhdTatPx;
}
double IncompressibleBackend::dsdTatPxdT(void) {
    if (!_dsdTatPxdT) _dsdTatPxdT = calc_dsdTatPxdT(_T, _p, _fractions[0]);
    return _dsdTatPxdT;
}
double IncompressibleBackend::dhdTatPxdT(void) {
    if (!_dhdTatPxdT) _dhdTatPxdT = calc_dhdTatPxdT(_T, _p, _fractions[0]);
    return _dhdTatPxdT;
}
double IncompressibleBackend::dsdpatTx(void) {
    if (!_dsdpatTx) _dsdpatTx = calc_dsdpatTx(rhomass(), drhodTatPx());
    return _dsdpatTx;
}
double IncompressibleBackend::dhdpatTx(void) {
    if (!_dhdpatTx) _dhdpatTx = calc_dhdpatTx(_T, rhomass(), drhodTatPx());
    return _dhdpatTx;
}

/// Return the temperature in K
double IncompressibleBackend::T_ref(void) {
    if (!_T_ref) throw ValueError("Reference temperature is not set");
    return _T_ref;
}
/// Return the pressure in Pa
double IncompressibleBackend::p_ref(void) {
    if (!_p_ref) throw ValueError("Reference pressure is not set");
    return _p_ref;
}
/// Return the composition
double IncompressibleBackend::x_ref(void) {
    if (!_x_ref) throw ValueError("Reference composition is not set");
    return _x_ref;
}
/// Return the mass enthalpy in J/kg
double IncompressibleBackend::h_ref(void) {
    if (!_h_ref) throw ValueError("Reference enthalpy is not set");
    return _h_ref;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::s_ref(void) {
    if (!_s_ref) throw ValueError("Reference entropy is not set");
    return _s_ref;
}

/// Return the mass enthalpy in J/kg
double IncompressibleBackend::hmass_ref(void) {
    if (!_hmass_ref) _hmass_ref = raw_calc_hmass(T_ref(), p_ref(), x_ref());
    return _hmass_ref;
}
/// Return the molar entropy in J/mol/K
double IncompressibleBackend::smass_ref(void) {
    if (!_smass_ref) _smass_ref = raw_calc_smass(T_ref(), p_ref(), x_ref());
    return _smass_ref;
}

/// Calculate T given pressure and density
/**
@param rhomass The mass density in kg/m^3
@param p The pressure in Pa
@returns T The temperature in K
*/
CoolPropDbl IncompressibleBackend::DmassP_flash(CoolPropDbl rhomass, CoolPropDbl p) {
    return fluid->T_rho(rhomass, p, _fractions[0]);
}
/// Calculate T given pressure and enthalpy
/**
@param hmass The mass enthalpy in J/kg
@param p The pressure in Pa
@returns T The temperature in K
*/
CoolPropDbl IncompressibleBackend::HmassP_flash(CoolPropDbl hmass, CoolPropDbl p) {

    class HmassP_residual : public FuncWrapper1D
    {
       protected:
        double p, x, h_in;
        IncompressibleBackend* backend;

       public:
        HmassP_residual(IncompressibleBackend* backend, const double& p, const double& x, const double& h_in)
          : p(p), x(x), h_in(h_in), backend(backend) {}
        double call(double target) {
            return backend->raw_calc_hmass(target, p, x) - h_in;  //fluid.u(target,p,x)+ p / fluid.rho(target,p,x) - h_in;
        }
        //double deriv(double target);
    };

    HmassP_residual res = HmassP_residual(this, p, _fractions[0], hmass - h_ref() + hmass_ref());

    double macheps = DBL_EPSILON;
    double tol = DBL_EPSILON * 1e3;
    int maxiter = 10;
    double result = Brent(&res, fluid->getTmin(), fluid->getTmax(), macheps, tol, maxiter);
    //if (this->do_debug()) std::cout << "Brent solver message: " << errstring << std::endl;
    return result;
}
/// Calculate T given pressure and entropy
/**
@param smass The mass entropy in J/kg/K
@param p The pressure in Pa
@returns T The temperature in K
*/
CoolPropDbl IncompressibleBackend::PSmass_flash(CoolPropDbl p, CoolPropDbl smass) {

    class PSmass_residual : public FuncWrapper1D
    {
       protected:
        double p, x, s_in;
        IncompressibleBackend* backend;

       public:
        PSmass_residual(IncompressibleBackend* backend, const double& p, const double& x, const double& s_in)
          : p(p), x(x), s_in(s_in), backend(backend) {}
        double call(double target) {
            return backend->raw_calc_smass(target, p, x) - s_in;
        }
    };

    PSmass_residual res = PSmass_residual(this, p, _fractions[0], smass - s_ref() + smass_ref());

    double macheps = DBL_EPSILON;
    double tol = DBL_EPSILON * 1e3;
    int maxiter = 10;
    double result = Brent(&res, fluid->getTmin(), fluid->getTmax(), macheps, tol, maxiter);
    //if (this->do_debug()) std::cout << "Brent solver message: " << errstring << std::endl;
    return result;
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
//	int    maxiter = 10;
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
CoolPropDbl IncompressibleBackend::calc_umass(void) {
    return hmass() - _p / rhomass();
};

/// ... and continue with the ones that depend on reference conditions.
CoolPropDbl IncompressibleBackend::calc_hmass(void) {
    return h_ref() + raw_calc_hmass(_T, _p, _fractions[0]) - hmass_ref();
};
CoolPropDbl IncompressibleBackend::calc_smass(void) {
    return s_ref() + raw_calc_smass(_T, _p, _fractions[0]) - smass_ref();
};

/// Functions that can be used with the solver, they miss the reference values!
CoolPropDbl IncompressibleBackend::raw_calc_hmass(double T, double p, double x) {
    return calc_dhdTatPxdT(T, p, x) + p * calc_dhdpatTx(T, fluid->rho(T, p, x), calc_drhodTatPx(T, p, x));
};
CoolPropDbl IncompressibleBackend::raw_calc_smass(double T, double p, double x) {
    return calc_dsdTatPxdT(T, p, x) + p * calc_dsdpatTx(fluid->rho(T, p, x), calc_drhodTatPx(T, p, x));
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

#    include "TestObjects.h"

TEST_CASE("Internal consistency checks and example use cases for the incompressible backend", "[IncompressibleBackend]") {
    CoolProp::IncompressibleFluid fluid = CoolProp::get_incompressible_fluid("Methanol");
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
        double acc = 0.0001;
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

#endif /* ENABLE_CATCH */
