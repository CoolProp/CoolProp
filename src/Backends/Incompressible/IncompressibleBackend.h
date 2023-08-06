
#ifndef INCOMPRESSIBLEBACKEND_H_
#define INCOMPRESSIBLEBACKEND_H_

#include "DataStructures.h"
#include "IncompressibleFluid.h"
#include "AbstractState.h"
#include "Exceptions.h"

#include <vector>

namespace CoolProp {

class IncompressibleBackend : public AbstractState
{

   protected:
    /// Bulk values, state variables
    //double _T, _p; // From AbstractState
    std::vector<CoolPropDbl> _fractions;

    /// Reference values, no need to calculate them each time
    CachedElement _T_ref, _p_ref, _x_ref, _h_ref, _s_ref;
    CachedElement _hmass_ref, _smass_ref;

    /// Additional cached elements used for the partial derivatives
    CachedElement _cmass, _hmass, _rhomass, _smass, _umass;
    CachedElement _drhodTatPx, _dsdTatPx, _dhdTatPx, _dsdTatPxdT, _dhdTatPxdT, _dsdpatTx, _dhdpatTx;

    IncompressibleFluid* fluid;

    /// Set the fractions
    /**
    @param fractions The vector of fractions of the components converted to the correct input
    */
    void set_fractions(const std::vector<CoolPropDbl>& fractions);

   public:
    IncompressibleBackend();
    virtual ~IncompressibleBackend(){};
    std::string backend_name(void) {
        return get_backend_string(INCOMP_BACKEND);
    }

    /// The instantiator
    /// @param fluid object, mostly for testing purposes
    IncompressibleBackend(IncompressibleFluid* fluid);
    /// The instantiator
    /// @param fluid_name the string with the fluid name
    IncompressibleBackend(const std::string& fluid_name);
    /// The instantiator
    /// @param component_names The vector of strings of the fluid components, without file ending
    IncompressibleBackend(const std::vector<std::string>& component_names);

    // Incompressible backend uses different compositions
    bool using_mole_fractions(void) {
        return this->fluid->getxid() == IFRAC_MOLE;
    };
    bool using_mass_fractions(void) {
        return (this->fluid->getxid() == IFRAC_MASS || this->fluid->getxid() == IFRAC_PURE);
    };
    bool using_volu_fractions(void) {
        return this->fluid->getxid() == IFRAC_VOLUME;
    };

    /// Updating function for incompressible fluid
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update all the internal variables that we can.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs input_pair, double value1, double value2);

    std::string fluid_param_string(const std::string& ParamName) {
        if (!ParamName.compare("long_name")) {
            return calc_name();
        } else {
            throw ValueError(format("Input value [%s] is invalid.", ParamName.c_str()));
        }
    }

    /// Clear all the cached values
    bool clear();

    /// Update the reference values and clear the state
    void set_reference_state(double T0 = 20 + 273.15, double p0 = 101325, double x0 = 0.0, double h0 = 0.0, double s0 = 0.0);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions);
    const std::vector<CoolPropDbl>& get_mole_fractions(void) {
        throw NotImplementedError("get_mole_fractions not implemented for this backend");
    };

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions);

    /// Set the volume fractions
    /**
    @param volu_fractions The vector of volume fractions of the components
    */
    void set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions);

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /** We have to override some of the functions from the AbstractState.
	 *  The incompressibles are only mass-based and do not support conversion
	 *  from molar to specific quantities.
	 *  We also have a few new chaced variables that we need.
	 */
    /// Return the mass density in kg/m^3
    double rhomass(void);
    /// Return the mass enthalpy in J/kg
    double hmass(void);
    /// Return the molar entropy in J/mol/K
    double smass(void);
    /// Return the molar internal energy in J/mol
    double umass(void);
    /// Return the mass constant pressure specific heat in J/kg/K
    double cmass(void);

    double drhodTatPx(void);
    double dsdTatPx(void);
    double dhdTatPx(void);
    double dsdTatPxdT(void);
    double dhdTatPxdT(void);
    double dsdpatTx(void);
    double dhdpatTx(void);

    /// Return the temperature in K
    double T_ref(void);
    /// Return the pressure in Pa
    double p_ref(void);
    /// Return the composition
    double x_ref(void);
    /// Return the mass enthalpy in J/kg
    double h_ref(void);
    /// Return the molar entropy in J/mol/K
    double s_ref(void);

    /// Return the mass enthalpy in J/kg
    double hmass_ref(void);
    /// Return the molar entropy in J/mol/K
    double smass_ref(void);

    /** These functions should be protected, but that requires new tests.
     *  I'll leave that as a TODO item for now.
     */
    /// Calculate T given pressure and density
    /**
    @param rhomass The mass density in kg/m^3
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    CoolPropDbl DmassP_flash(CoolPropDbl rhomass, CoolPropDbl p);
    /// Calculate T given pressure and enthalpy
    /**
    @param hmass The mass enthalpy in J/kg
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    CoolPropDbl HmassP_flash(CoolPropDbl hmass, CoolPropDbl p);
    /// Calculate T given pressure and entropy
    /**
    @param smass The mass entropy in J/kg/K
    @param p The pressure in Pa
    @returns T The temperature in K
    */
    CoolPropDbl PSmass_flash(CoolPropDbl p, CoolPropDbl smass);

    //    /// Calculate T given pressure and internal energy
    //    /**
    //    @param umass The mass internal energy in J/kg
    //    @param p The pressure in Pa
    //    @returns T The temperature in K
    //    */
    //    CoolPropDbl PUmass_flash(CoolPropDbl p, CoolPropDbl umass);

    /// We start with the functions that do not need a reference state
    CoolPropDbl calc_rhomass(void) {
        return fluid->rho(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_cmass(void) {
        return fluid->c(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_cpmass(void) {
        return cmass();
    };
    CoolPropDbl calc_cvmass(void) {
        return cmass();
    };
    CoolPropDbl calc_viscosity(void) {
        return fluid->visc(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_conductivity(void) {
        return fluid->cond(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_T_freeze(void) {
        // No update is called - T_freeze is a trivial output
        fluid->checkX(_fractions[0]);
        return fluid->Tfreeze(_p, _fractions[0]);
    };
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value);
    CoolPropDbl calc_umass(void);

    /// ... and continue with the ones that depend on reference conditions.
    CoolPropDbl calc_hmass(void);
    CoolPropDbl calc_smass(void);

   public:
    /// Functions that can be used with the solver, they miss the reference values!
    CoolPropDbl raw_calc_hmass(double T, double p, double x);
    CoolPropDbl raw_calc_smass(double T, double p, double x);

   protected:
    /// Calculate the first partial derivative for the desired derivative
    CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);

    /* Below are direct calculations of the derivatives. Nothing
	 * special is going on, we simply use the polynomial class to
	 * derive the different functions with respect to temperature.
	 */
    /// Partial derivative of density with respect to temperature at constant pressure and composition
    double calc_drhodTatPx(double T, double p, double x) {
        return fluid->drhodTatPx(T, p, x);
    };
    /// Partial derivative of entropy with respect to temperature at constant pressure and composition
    double calc_dsdTatPx(double T, double p, double x) {
        return fluid->c(T, p, x) / T;
    };
    /// Partial derivative of enthalpy with respect to temperature at constant pressure and composition
    double calc_dhdTatPx(double T, double p, double x) {
        return fluid->c(T, p, x);
    };
    /// Partial derivative of entropy
    ///  with respect to temperature at constant pressure and composition
    ///  integrated in temperature
    double calc_dsdTatPxdT(double T, double p, double x) {
        return fluid->dsdTatPxdT(T, p, x);
    };
    /// Partial derivative of enthalpy
    ///  with respect to temperature at constant pressure and composition
    ///  integrated in temperature
    double calc_dhdTatPxdT(double T, double p, double x) {
        return fluid->dhdTatPxdT(T, p, x);
    };

    /* Other useful derivatives
	 */
    /// Partial derivative of entropy with respect to pressure at constant temperature and composition
    ///  \f[ \left( \frac{\partial s}{\partial p} \right)_T = - \left( \frac{\partial v}{\partial T} \right)_p = \rho^{-2} \left( \frac{\partial \rho}{\partial T} \right)_p \right) \f]
    double calc_dsdpatTx(double rho, double drhodTatPx);
    /// Partial derivative of enthalpy with respect to pressure at constant temperature and composition
    ///  \f[ \left( \frac{\partial h}{\partial p} \right)_T = v - T \left( \frac{\partial v}{\partial T} \right)_p = \rho^{-1} \left( 1 + T \rho^{-1} \left( \frac{\partial \rho}{\partial T} \right)_p \right) \f]
    double calc_dhdpatTx(double T, double rho, double drhodTatPx);

   public:
    /// Constants from the fluid object
    CoolPropDbl calc_Tmax(void) {
        return fluid->getTmax();
    };
    CoolPropDbl calc_Tmin(void) {
        return fluid->getTmin();
    };
    CoolPropDbl calc_fraction_min(void) {
        return fluid->getxmin();
    };
    CoolPropDbl calc_fraction_max(void) {
        return fluid->getxmax();
    };
    std::string calc_name(void) {
        return fluid->getName();
    };
    std::string calc_description(void) {
        return fluid->getDescription();
    };
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEBACKEND_H_ */
