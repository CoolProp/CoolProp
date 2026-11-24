
#ifndef INCOMPRESSIBLEBACKEND_H_
#define INCOMPRESSIBLEBACKEND_H_

#include "CoolProp/DataStructures.h"
#include "CoolProp/fluids/IncompressibleFluid.h"
#include "CoolProp/AbstractState.h"
#include "CoolProp/Exceptions.h"

#include <algorithm>
#include <cmath>
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
    /// data_version of the fluid that _hmass_ref and _smass_ref were computed from
    unsigned long ref_data_version = 0;
    /// True while _T_ref is default_T_ref_K(), which depends on the fluid's range
    bool T_ref_is_default = false;
    void drop_stale_reference_values() {
        if (ref_data_version != fluid->getDataVersion()) {
            _hmass_ref.clear();
            _smass_ref.clear();
            // a replaced fluid can have another range, which moves the default
            if (T_ref_is_default) _T_ref = default_T_ref_K();
            ref_data_version = fluid->getDataVersion();
        }
    }

    /// Additional cached elements used for the partial derivatives
    CachedElement _cmass, _hmass, _rhomass, _smass, _umass;
    CachedElement _drhodTatPx, _dsdTatPx, _dhdTatPx, _dsdTatPxdT, _dhdTatPxdT, _dsdpatTx, _dhdpatTx;

    IncompressibleFluid* fluid;

    /// Evaluate density and specific heat from the fluid's Chebyshev fits
    /// (true) or its classic fits (false), see IncompressibleFluid::useChebyshevFit.
    /// Read once from INCOMP_PREFER_CHEBYSHEV_CALORIC when the state is
    /// created and never changed afterwards: the reference enthalpy and
    /// entropy cached in _hmass_ref and _smass_ref were computed from the
    /// fits selected here, and mixing them with the other fits would shift
    /// h and s by a constant.
    bool chebyshev_caloric = true;

    /// The fits of this state with its composition put in (see FitAtComposition),
    /// built when the composition changes or the fluid's data is replaced.
    CompositionFits fits;
    /// The reduced fits if x is this state's composition (the reference state
    /// uses the same one), else nullptr: other compositions use the full
    /// evaluation.
    const CompositionFits* fits_at(double x) {
        if (_fractions.size() != 1 || x != _fractions[0]) return nullptr;
        if (fits.x != x || fits.chebyshev != chebyshev_caloric || fits.data_version != fluid->getDataVersion()) {
            fluid->fixComposition(x, chebyshev_caloric, fits);
        }
        return &fits;
    }
    double rho_at(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        return (f && f->rho.kind != FitAtComposition::NONE) ? f->rho.eval(T) : fluid->rho(T, p, x, chebyshev_caloric);
    }
    double d2rhodT2_at(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        return (f && f->d2rhodT2.kind != FitAtComposition::NONE) ? f->d2rhodT2.eval(T) : fluid->d2rhodT2atPx(T, p, x, chebyshev_caloric);
    }
    double c_at(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        return (f && f->c.kind != FitAtComposition::NONE) ? f->c.eval(T) : fluid->c(T, p, x, chebyshev_caloric);
    }

    /// Set the fractions
    /**
    @param fractions The vector of fractions of the components converted to the correct input
    */
    void set_fractions(const std::vector<CoolPropDbl>& fractions);

   public:
    IncompressibleBackend();
    virtual ~IncompressibleBackend() = default;
    std::string backend_name() override {
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
    bool using_mole_fractions() override {
        return this->fluid->getxid() == IFRAC_MOLE;
    };
    bool using_mass_fractions() override {
        return (this->fluid->getxid() == IFRAC_MASS || this->fluid->getxid() == IFRAC_PURE);
    };
    bool using_volu_fractions() override {
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
    void update(CoolProp::input_pairs input_pair, double value1, double value2) override;

    std::string fluid_param_string(const std::string& ParamName) override {
        if (!ParamName.compare("long_name")) {
            return calc_name();
        } else {
            throw ValueError(format("Input value [%s] is invalid.", ParamName.c_str()));
        }
    }

    /// Clear all the cached values
    bool clear() override;

    /// Update the reference values and clear the state. An explicit T0 is kept
    /// as given; the constructors use default_T_ref_K() instead of the 293.15 K
    /// default argument.
    void set_reference_state(double T0 = 20 + 273.15, double p0 = 101325, double x0 = 0.0, double h0 = 0.0, double s0 = 0.0);
    /// Default reference temperature: 293.15 K, or the nearest end of the
    /// fluid's range when 293.15 K lies outside it. The fits are not valid
    /// outside the range and a Chebyshev fit grows without bound there, so a
    /// reference value taken there is meaningless and costs precision in every
    /// h and s (IceNA's reference entropy was 2.4e11 J/kg/K, COO-203).
    /// A range that is not set yet (infinite or NaN) keeps 293.15 K.
    [[nodiscard]] double default_T_ref_K() const {
        const double T_default_K = 20 + 273.15, T_min_K = fluid->getTmin(), T_max_K = fluid->getTmax();
        if (!std::isfinite(T_min_K) || !std::isfinite(T_max_K) || T_min_K > T_max_K) return T_default_K;
        return std::min(std::max(T_default_K, T_min_K), T_max_K);
    }
    /// Set the default reference state (see default_T_ref_K)
    void set_default_reference_state() {
        set_reference_state(default_T_ref_K());
        T_ref_is_default = true;
    }

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) override;
    const std::vector<CoolPropDbl>& get_mole_fractions() override {
        throw NotImplementedError("get_mole_fractions not implemented for this backend");
    };

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) override;

    /// Set the volume fractions
    /**
    @param volu_fractions The vector of volume fractions of the components
    */
    void set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions) override;

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /** We have to override some of the functions from the AbstractState.
	 *  The incompressibles are only mass-based and do not support conversion
	 *  from molar to specific quantities.
	 *  We also have a few new chaced variables that we need.
	 */
    /// Return the mass density in kg/m^3
    double rhomass();
    /// Return the mass enthalpy in J/kg
    double hmass();
    /// Return the molar entropy in J/mol/K
    double smass();
    /// Return the molar internal energy in J/mol
    double umass();
    /// Return the mass constant pressure specific heat in J/kg/K
    double cmass();

    double drhodTatPx();
    double dsdTatPx();
    double dhdTatPx();
    double dsdTatPxdT();
    double dhdTatPxdT();
    double dsdpatTx();
    double dhdpatTx();

    /// Return the temperature in K
    double T_ref();
    /// Return the pressure in Pa
    double p_ref();
    /// Return the composition
    double x_ref();
    /// Return the mass enthalpy in J/kg
    double h_ref();
    /// Return the molar entropy in J/mol/K
    double s_ref();

    /// Return the mass enthalpy in J/kg
    double hmass_ref();
    /// Return the molar entropy in J/mol/K
    double smass_ref();

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

    // CoolPropDbl PUmass_flash(CoolPropDbl p, CoolPropDbl umass);  // not implemented

    // Molar quantities are not meaningful for incompressible solutions,
    // which are formulated entirely on a mass basis with no fixed
    // composition / molecular weight (#1908). Throw a controlled
    // exception with a helpful pointer at the mass-basis equivalents
    // instead of silently returning the AbstractState default
    // (-_HUGE) which surfaced as "PropsSI failed ungracefully".
    CoolPropDbl calc_molar_mass() override {
        throw NotImplementedError("Molar mass is not defined for the INCOMP (incompressible) backend; INCOMP fluids are mass-based.");
    }
    CoolPropDbl calc_rhomolar() override {
        throw NotImplementedError("Dmolar / rhomolar is not defined for the INCOMP backend; use Dmass / rhomass instead.");
    }
    CoolPropDbl calc_hmolar() override {
        throw NotImplementedError("Hmolar / hmolar is not defined for the INCOMP backend; use Hmass / hmass instead.");
    }
    CoolPropDbl calc_smolar() override {
        throw NotImplementedError("Smolar / smolar is not defined for the INCOMP backend; use Smass / smass instead.");
    }
    CoolPropDbl calc_umolar() override {
        throw NotImplementedError("Umolar / umolar is not defined for the INCOMP backend; use Umass / umass instead.");
    }
    CoolPropDbl calc_cpmolar() override {
        throw NotImplementedError("Cpmolar / cpmolar is not defined for the INCOMP backend; use Cpmass / cpmass instead.");
    }
    CoolPropDbl calc_cvmolar() override {
        throw NotImplementedError("Cvmolar / cvmolar is not defined for the INCOMP backend; use Cvmass / cvmass instead.");
    }

    /// We start with the functions that do not need a reference state
    CoolPropDbl calc_rhomass() override {
        return rho_at(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_cmass() {
        return c_at(_T, _p, _fractions[0]);
    };
    /// Specific heat at constant pressure, the slope of h in T at constant p.
    /// It differs from the fitted heat capacity c by -(p - p_fit) T d2v/dT2, the
    /// pressure term of h (see calc_dhdTatPx), so it agrees with h and s at any
    /// pressure and equals c at the pressure of the data.
    CoolPropDbl calc_cpmass() override {
        return dhdTatPx();
    };
    /// Specific heat at constant volume, returned as the fitted heat capacity c.
    /// With v a function of T only, the isothermal compressibility is zero and
    /// cp - cv = T (dv/dT)^2 / (v kappa_T) cannot be evaluated: cv is not defined
    /// by this model. c is the fitted value, which cp equals at p_fit_Pa.
    CoolPropDbl calc_cvmass() override {
        return cmass();
    };
    CoolPropDbl calc_viscosity() override {
        const CompositionFits* f = fits_at(_fractions[0]);
        return (f && f->visc.kind != FitAtComposition::NONE) ? f->visc.eval(_T) : fluid->visc(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_conductivity() override {
        const CompositionFits* f = fits_at(_fractions[0]);
        return (f && f->cond.kind != FitAtComposition::NONE) ? f->cond.eval(_T) : fluid->cond(_T, _p, _fractions[0]);
    };
    CoolPropDbl calc_T_freeze() override {
        // No update is called - T_freeze is a trivial output
        fluid->checkX(_fractions[0]);
        return fluid->Tfreeze(_p, _fractions[0]);
    };
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value) override;
    CoolPropDbl calc_umass() override;

    /// ... and continue with the ones that depend on reference conditions.
    CoolPropDbl calc_hmass() override;
    CoolPropDbl calc_smass() override;

    /// Pressure of the fitted data, 1 atm. The pressure terms of h and s are
    /// taken from here, so cp equals the fitted c at this pressure (#3577).
    static constexpr double p_fit_Pa = 101325.0;
    /// h and s without the reference offsets, from rho and drho/dT at T:
    /// h = Int c dT + (p - p_fit) dh/dp|T and s = Int c/T dT + (p - p_fit) ds/dp|T.
    /// The one place these are written; raw_calc_hmass/smass and the flash
    /// residual use them.
    double raw_hmass_from(double T, double p, double x, double rho, double drhodT) {
        return calc_dhdTatPxdT(T, p, x) + (p - p_fit_Pa) * calc_dhdpatTx(T, rho, drhodT);
    }
    double raw_smass_from(double T, double p, double x, double rho, double drhodT) {
        return calc_dsdTatPxdT(T, p, x) + (p - p_fit_Pa) * calc_dsdpatTx(rho, drhodT);
    }
    /// Residual of the HmassP and PSmass flashes, see IncompressibleBackend.cpp
    class CaloricResidual;
    double flash_first_guess(bool enthalpy, double p, double raw_target);
    /// Values at T_ref that flash_first_guess needs, kept per composition
    struct GuessCache
    {
        double x = NAN, T_ref = NAN, rho = NAN, drho = NAN, c = NAN;
        unsigned long data_version = 0;
    } guess_cache;

   public:
    /// Functions that can be used with the solver, they miss the reference values!
    CoolPropDbl raw_calc_hmass(double T, double p, double x);
    CoolPropDbl raw_calc_smass(double T, double p, double x);

   protected:
    /// Calculate the first partial derivative for the desired derivative
    CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) override;

    /* Below are direct calculations of the derivatives. Nothing
	 * special is going on, we simply use the polynomial class to
	 * derive the different functions with respect to temperature.
	 */
    /// Partial derivative of density with respect to temperature at constant pressure and composition
    double calc_drhodTatPx(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        return (f && f->drhodT.kind != FitAtComposition::NONE) ? f->drhodT.eval(T) : fluid->drhodTatPx(T, p, x, chebyshev_caloric);
    };
    /// Second derivative of v = 1/rho in T at constant p, from rho, drho/dT and d2rho/dT2
    static double d2vdT2_from(double rho, double drhodT, double d2rhodT2) {
        return 2 * drhodT * drhodT / (rho * rho * rho) - d2rhodT2 / (rho * rho);
    }
    /// Slopes of h and s in T at constant p and composition, given d2v/dT2 at T.
    /// h = Int c dT + (p - p_fit) (v - T dv/dT) and s = Int c/T dT - (p - p_fit) dv/dT, so
    ///   dh/dT = c - (p - p_fit) T d2v/dT2  and  ds/dT = c/T - (p - p_fit) d2v/dT2.
    /// The one place these are written; cpmass, the partial derivatives and the
    /// flash residual all use them.
    static double dhdT_from(double T, double p, double c, double d2vdT2) {
        return c - (p - p_fit_Pa) * T * d2vdT2;
    }
    static double dsdT_from(double T, double p, double c, double d2vdT2) {
        return c / T - (p - p_fit_Pa) * d2vdT2;
    }
    double d2vdT2_at(double T, double p, double x) {
        return d2vdT2_from(rho_at(T, p, x), calc_drhodTatPx(T, p, x), d2rhodT2_at(T, p, x));
    }
    /// Partial derivative of entropy with respect to temperature at constant pressure and composition
    double calc_dsdTatPx(double T, double p, double x) {
        return dsdT_from(T, p, c_at(T, p, x), d2vdT2_at(T, p, x));
    };
    /// Partial derivative of enthalpy with respect to temperature at constant pressure and composition
    double calc_dhdTatPx(double T, double p, double x) {
        return dhdT_from(T, p, c_at(T, p, x), d2vdT2_at(T, p, x));
    };
    /// Partial derivative of entropy
    ///  with respect to temperature at constant pressure and composition
    ///  integrated in temperature
    double calc_dsdTatPxdT(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        if (f && f->int_c_divT.kind != FitAtComposition::NONE) {
            // Polynomial fits add a ln(T) term, as IncompressibleFluid::dsdTatPxdT does
            double value = f->int_c_divT.eval(T);
            if (f->int_c_divT_log.kind != FitAtComposition::NONE) value += std::log(T) * f->int_c_divT_log.eval(T);
            return value;
        }
        return fluid->dsdTatPxdT(T, p, x, chebyshev_caloric);
    };
    /// Partial derivative of enthalpy
    ///  with respect to temperature at constant pressure and composition
    ///  integrated in temperature
    double calc_dhdTatPxdT(double T, double p, double x) {
        const CompositionFits* f = fits_at(x);
        return (f && f->int_c.kind != FitAtComposition::NONE) ? f->int_c.eval(T) : fluid->dhdTatPxdT(T, p, x, chebyshev_caloric);
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
    CoolPropDbl calc_Tmax() override {
        return fluid->getTmax();
    };
    CoolPropDbl calc_Tmin() override {
        return fluid->getTmin();
    };
    CoolPropDbl calc_fraction_min() override {
        return fluid->getxmin();
    };
    CoolPropDbl calc_fraction_max() override {
        return fluid->getxmax();
    };
    std::string calc_name() override {
        return fluid->getName();
    };
    std::string calc_description() override {
        return fluid->getDescription();
    };
};

} /* namespace CoolProp */
#endif /* INCOMPRESSIBLEBACKEND_H_ */
