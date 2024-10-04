/*
 * AbstractState.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef ABSTRACTSTATE_H_
#define ABSTRACTSTATE_H_

#include "CachedElement.h"
#include "Exceptions.h"
#include "DataStructures.h"
#include "PhaseEnvelope.h"
#include "crossplatform_shared_ptr.h"

#include <numeric>

namespace CoolProp {

/// This structure holds values obtained while tracing the spinodal curve
/// (most often in the process of finding critical points, but not only)
class SpinodalData
{
   public:
    std::vector<double> tau,  ///< The reciprocal reduced temperature (\f$\tau=T_r/T\f$)
      delta,                  ///< The reduced density (\f$\delta=\rho/\rho_r\f$)
      M1;                     ///< The determinant of the scaled matrix for the second criticality condition
};

/// This simple class holds the values for guesses for use in some solvers
/// that have the ability to use guess values intelligently
class GuessesStructure
{
   public:
    double T,               ///< temperature in K
      p,                    ///< pressure in Pa
      rhomolar,             ///< molar density in mol/m^3
      hmolar,               ///< molar enthalpy in J/mol
      smolar,               ///< molar entropy in J/mol/K
      rhomolar_liq,         ///< molar density of the liquid phase in mol/m^3
      rhomolar_vap;         ///< molar density of the vapor phase in mol/m^3
    std::vector<double> x,  ///< molar composition of the liquid phase
      y;                    ///< molar composition of the vapor phase
    GuessesStructure() {
        clear();
    };
    void clear() {
        T = _HUGE;
        p = _HUGE;
        rhomolar = _HUGE;
        hmolar = _HUGE;
        smolar = _HUGE;
        rhomolar_liq = _HUGE;
        rhomolar_vap = _HUGE;
        x.clear(), y.clear();
    }
};

//! The mother of all state classes
/*!
This class provides the basic properties based on interrelations of the
properties, their derivatives and the Helmholtz energy terms. It does not
provide the mechanism to update the values. This has to be implemented in
a subclass. Most functions are defined as virtual functions allowing us
redefine them later, for example to implement the TTSE technique. The
functions defined here are always used as a fall-back.

This base class does not perform any checks on the two-phase conditions and
alike. Most of the functions defined here only apply to compressible single
state substances. Make sure you are aware of all the assumptions we made
when using this class.

Add build table function to Abstract State
Interpolator inherit AS implemented by TTSE BICUBIC

*/
class AbstractState
{
   protected:
    /// Some administrative variables
    long _fluid_type;
    phases _phase;               ///< The key for the phase from CoolProp::phases enum
    phases imposed_phase_index;  ///< If the phase is imposed, the imposed phase index

    bool isSupercriticalPhase(void) {
        return (this->_phase == iphase_supercritical || this->_phase == iphase_supercritical_liquid || this->_phase == iphase_supercritical_gas);
    }

    bool isHomogeneousPhase(void) {
        return (this->_phase == iphase_liquid || this->_phase == iphase_gas || isSupercriticalPhase() || this->_phase == iphase_critical_point);
    }

    bool isTwoPhase(void) {
        return (this->_phase == iphase_twophase);
    }

    /// Two important points
    SimpleState _critical, _reducing;

    /// Molar mass [mol/kg]
    CachedElement _molar_mass;

    /// Universal gas constant [J/mol/K]
    CachedElement _gas_constant;

    /// Bulk values
    double _rhomolar, _T, _p, _Q, _R;

    CachedElement _tau, _delta;

    /// Transport properties
    CachedElement _viscosity, _conductivity, _surface_tension;

    CachedElement _hmolar, _smolar, _umolar, _logp, _logrhomolar, _cpmolar, _cp0molar, _cvmolar, _speed_sound, _gibbsmolar, _helmholtzmolar;

    /// Residual properties
    CachedElement _hmolar_residual, _smolar_residual, _gibbsmolar_residual;

    /// Excess properties
    CachedElement _hmolar_excess, _smolar_excess, _gibbsmolar_excess, _umolar_excess, _volumemolar_excess, _helmholtzmolar_excess;

    /// Ancillary values
    CachedElement _rhoLanc, _rhoVanc, _pLanc, _pVanc, _TLanc, _TVanc;

    CachedElement _fugacity_coefficient;

    /// Smoothing values
    CachedElement _rho_spline, _drho_spline_dh__constp, _drho_spline_dp__consth;

    /// Cached low-level elements for in-place calculation of other properties
    CachedElement _alpha0, _dalpha0_dTau, _dalpha0_dDelta, _d2alpha0_dTau2, _d2alpha0_dDelta_dTau, _d2alpha0_dDelta2, _d3alpha0_dTau3,
      _d3alpha0_dDelta_dTau2, _d3alpha0_dDelta2_dTau, _d3alpha0_dDelta3, _alphar, _dalphar_dTau, _dalphar_dDelta, _d2alphar_dTau2,
      _d2alphar_dDelta_dTau, _d2alphar_dDelta2, _d3alphar_dTau3, _d3alphar_dDelta_dTau2, _d3alphar_dDelta2_dTau, _d3alphar_dDelta3, _d4alphar_dTau4,
      _d4alphar_dDelta_dTau3, _d4alphar_dDelta2_dTau2, _d4alphar_dDelta3_dTau, _d4alphar_dDelta4;

    CachedElement _dalphar_dDelta_lim, _d2alphar_dDelta2_lim, _d2alphar_dDelta_dTau_lim, _d3alphar_dDelta2_dTau_lim;

    /// Two-Phase variables
    CachedElement _rhoLmolar, _rhoVmolar;

    // ----------------------------------------
    // Property accessors to be optionally implemented by the backend
    // for properties that are not always calculated
    // ----------------------------------------
    /// Using this backend, calculate the molar enthalpy in J/mol
    virtual CoolPropDbl calc_hmolar(void) {
        throw NotImplementedError("calc_hmolar is not implemented for this backend");
    };
    /// Using this backend, calculate the residual molar enthalpy in J/mol
    virtual CoolPropDbl calc_hmolar_residual(void) {
        throw NotImplementedError("calc_hmolar_residual is not implemented for this backend");
    };
    /// Using this backend, calculate the molar entropy in J/mol/K
    virtual CoolPropDbl calc_smolar(void) {
        throw NotImplementedError("calc_smolar is not implemented for this backend");
    };
    /// Using this backend, calculate the residual molar entropy in J/mol/K
    virtual CoolPropDbl calc_smolar_residual(void) {
        throw NotImplementedError("calc_smolar_residual is not implemented for this backend");
    };
    /// Using this backend, calculate effective hardness of interaction
    virtual CoolPropDbl calc_neff(void) {
        throw NotImplementedError("calc_neff is not implemented for this backend");
    };
    /// Using this backend, calculate the molar internal energy in J/mol
    virtual CoolPropDbl calc_umolar(void) {
        throw NotImplementedError("calc_umolar is not implemented for this backend");
    };
    /// Using this backend, calculate the molar constant-pressure specific heat in J/mol/K
    virtual CoolPropDbl calc_cpmolar(void) {
        throw NotImplementedError("calc_cpmolar is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal gas molar constant-pressure specific heat in J/mol/K
    virtual CoolPropDbl calc_cpmolar_idealgas(void) {
        throw NotImplementedError("calc_cpmolar_idealgas is not implemented for this backend");
    };
    /// Using this backend, calculate the molar constant-volume specific heat in J/mol/K
    virtual CoolPropDbl calc_cvmolar(void) {
        throw NotImplementedError("calc_cvmolar is not implemented for this backend");
    };
    /// Using this backend, calculate the molar Gibbs function in J/mol
    virtual CoolPropDbl calc_gibbsmolar(void) {
        throw NotImplementedError("calc_gibbsmolar is not implemented for this backend");
    };
    /// Using this backend, calculate the residual molar Gibbs function in J/mol
    virtual CoolPropDbl calc_gibbsmolar_residual(void) {
        throw NotImplementedError("calc_gibbsmolar_residual is not implemented for this backend");
    };
    /// Using this backend, calculate the molar Helmholtz energy in J/mol
    virtual CoolPropDbl calc_helmholtzmolar(void) {
        throw NotImplementedError("calc_helmholtzmolar is not implemented for this backend");
    };
    /// Using this backend, calculate the speed of sound in m/s
    virtual CoolPropDbl calc_speed_sound(void) {
        throw NotImplementedError("calc_speed_sound is not implemented for this backend");
    };
    /// Using this backend, calculate the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    virtual CoolPropDbl calc_isothermal_compressibility(void) {
        throw NotImplementedError("calc_isothermal_compressibility is not implemented for this backend");
    };
    /// Using this backend, calculate the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    virtual CoolPropDbl calc_isobaric_expansion_coefficient(void) {
        throw NotImplementedError("calc_isobaric_expansion_coefficient is not implemented for this backend");
    };
    /// Using this backend, calculate the isentropic expansion coefficient \f$ \kappa_s = -\frac{c_p}{c_v}\frac{v}{p}\left.\frac{\partial p}{\partial v}\right|_T = \frac{\rho}{p}\left.\frac{\partial p}{\partial \rho}\right|_s\f$
    virtual CoolPropDbl calc_isentropic_expansion_coefficient(void) {
        throw NotImplementedError("calc_isentropic_expansion_coefficient is not implemented for this backend");
    };
    /// Using this backend, calculate the viscosity in Pa-s
    virtual CoolPropDbl calc_viscosity(void) {
        throw NotImplementedError("calc_viscosity is not implemented for this backend");
    };
    /// Using this backend, calculate the thermal conductivity in W/m/K
    virtual CoolPropDbl calc_conductivity(void) {
        throw NotImplementedError("calc_conductivity is not implemented for this backend");
    };
    /// Using this backend, calculate the surface tension in N/m
    virtual CoolPropDbl calc_surface_tension(void) {
        throw NotImplementedError("calc_surface_tension is not implemented for this backend");
    };
    /// Using this backend, calculate the molar mass in kg/mol
    virtual CoolPropDbl calc_molar_mass(void) {
        throw NotImplementedError("calc_molar_mass is not implemented for this backend");
    };
    /// Using this backend, calculate the acentric factor
    virtual CoolPropDbl calc_acentric_factor(void) {
        throw NotImplementedError("calc_acentric_factor is not implemented for this backend");
    };
    /// Using this backend, calculate the pressure in Pa
    virtual CoolPropDbl calc_pressure(void) {
        throw NotImplementedError("calc_pressure is not implemented for this backend");
    };
    /// Using this backend, calculate the universal gas constant \f$R_u\f$ in J/mol/K
    virtual CoolPropDbl calc_gas_constant(void) {
        throw NotImplementedError("calc_gas_constant is not implemented for this backend");
    };
    /// Using this backend, calculate the fugacity coefficient (dimensionless)
    virtual CoolPropDbl calc_fugacity_coefficient(std::size_t i) {
        throw NotImplementedError("calc_fugacity_coefficient is not implemented for this backend");
    };
    /// Using this backend, calculate the fugacity in Pa
    virtual std::vector<CoolPropDbl> calc_fugacity_coefficients() {
        throw NotImplementedError("calc_fugacity_coefficients is not implemented for this backend");
    };
    /// Using this backend, calculate the fugacity in Pa
    virtual CoolPropDbl calc_fugacity(std::size_t i) {
        throw NotImplementedError("calc_fugacity is not implemented for this backend");
    };
    /// Using this backend, calculate the chemical potential in J/mol
    virtual CoolPropDbl calc_chemical_potential(std::size_t i) {
        throw NotImplementedError("calc_chemical_potential is not implemented for this backend");
    };
    /// Using this backend, calculate the phase identification parameter (PIP)
    virtual CoolPropDbl calc_PIP(void) {
        throw NotImplementedError("calc_PIP is not implemented for this backend");
    };

    // Excess properties
    /// Using this backend, calculate and cache the excess properties
    virtual void calc_excess_properties(void) {
        throw NotImplementedError("calc_excess_properties is not implemented for this backend");
    };

    // Derivatives of residual helmholtz energy
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    virtual CoolPropDbl calc_alphar(void) {
        throw NotImplementedError("calc_alphar is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalphar_dDelta(void) {
        throw NotImplementedError("calc_dalphar_dDelta is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalphar_dTau(void) {
        throw NotImplementedError("calc_dalphar_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dDelta2(void) {
        throw NotImplementedError("calc_d2alphar_dDelta2 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dDelta_dTau(void) {
        throw NotImplementedError("calc_d2alphar_dDelta_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dTau2(void) {
        throw NotImplementedError("calc_d2alphar_dTau2 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta3(void) {
        throw NotImplementedError("calc_d3alphar_dDelta3 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta2_dTau(void) {
        throw NotImplementedError("calc_d3alphar_dDelta2_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta_dTau2(void) {
        throw NotImplementedError("calc_d3alphar_dDelta_dTau2 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dTau3(void) {
        throw NotImplementedError("calc_d3alphar_dTau3 is not implemented for this backend");
    };

    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d4alphar_dDelta4(void) {
        throw NotImplementedError("calc_d4alphar_dDelta4 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d4alphar_dDelta3_dTau(void) {
        throw NotImplementedError("calc_d4alphar_dDelta3_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d4alphar_dDelta2_dTau2(void) {
        throw NotImplementedError("calc_d4alphar_dDelta2_dTau2 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d4alphar_dDelta_dTau3(void) {
        throw NotImplementedError("calc_d4alphar_dDelta_dTau3 is not implemented for this backend");
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d4alphar_dTau4(void) {
        throw NotImplementedError("calc_d4alphar_dTau4 is not implemented for this backend");
    };

    // Derivatives of ideal-gas helmholtz energy
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0\f$ (dimensionless)
    virtual CoolPropDbl calc_alpha0(void) {
        throw NotImplementedError("calc_alpha0 is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalpha0_dDelta(void) {
        throw NotImplementedError("calc_dalpha0_dDelta is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalpha0_dTau(void) {
        throw NotImplementedError("calc_dalpha0_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dDelta_dTau(void) {
        throw NotImplementedError("calc_d2alpha0_dDelta_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dDelta2(void) {
        throw NotImplementedError("calc_d2alpha0_dDelta2 is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dTau2(void) {
        throw NotImplementedError("calc_d2alpha0_dTau2 is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta3(void) {
        throw NotImplementedError("calc_d3alpha0_dDelta3 is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta2_dTau(void) {
        throw NotImplementedError("calc_d3alpha0_dDelta2_dTau is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta_dTau2(void) {
        throw NotImplementedError("calc_d3alpha0_dDelta_dTau2 is not implemented for this backend");
    };
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dTau3(void) {
        throw NotImplementedError("calc_d3alpha0_dTau3 is not implemented for this backend");
    };

    virtual void calc_reducing_state(void) {
        throw NotImplementedError("calc_reducing_state is not implemented for this backend");
    };

    /// Using this backend, calculate the maximum temperature in K
    virtual CoolPropDbl calc_Tmax(void) {
        throw NotImplementedError("calc_Tmax is not implemented for this backend");
    };
    /// Using this backend, calculate the minimum temperature in K
    virtual CoolPropDbl calc_Tmin(void) {
        throw NotImplementedError("calc_Tmin is not implemented for this backend");
    };
    /// Using this backend, calculate the maximum pressure in Pa
    virtual CoolPropDbl calc_pmax(void) {
        throw NotImplementedError("calc_pmax is not implemented for this backend");
    };

    /// Using this backend, calculate the 20-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP20(void) {
        throw NotImplementedError("calc_GWP20 is not implemented for this backend");
    };
    /// Using this backend, calculate the 100-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP100(void) {
        throw NotImplementedError("calc_GWP100 is not implemented for this backend");
    };
    /// Using this backend, calculate the 500-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP500(void) {
        throw NotImplementedError("calc_GWP500 is not implemented for this backend");
    };
    /// Using this backend, calculate the ozone depletion potential (ODP)
    virtual CoolPropDbl calc_ODP(void) {
        throw NotImplementedError("calc_ODP is not implemented for this backend");
    };
    /// Using this backend, calculate the flame hazard
    virtual CoolPropDbl calc_flame_hazard(void) {
        throw NotImplementedError("calc_flame_hazard is not implemented for this backend");
    };
    /// Using this backend, calculate the health hazard
    virtual CoolPropDbl calc_health_hazard(void) {
        throw NotImplementedError("calc_health_hazard is not implemented for this backend");
    };
    /// Using this backend, calculate the physical hazard
    virtual CoolPropDbl calc_physical_hazard(void) {
        throw NotImplementedError("calc_physical_hazard is not implemented for this backend");
    };
    /// Using this backend, calculate the dipole moment in C-m (1 D = 3.33564e-30 C-m)
    virtual CoolPropDbl calc_dipole_moment(void) {
        throw NotImplementedError("calc_dipole_moment is not implemented for this backend");
    };

    /// Calculate the first partial derivative for the desired derivative
    virtual CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);
    /// Calculate the second partial derivative using the given backend
    virtual CoolPropDbl calc_second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2);

    /// Using this backend, calculate the reduced density (rho/rhoc)
    virtual CoolPropDbl calc_reduced_density(void) {
        throw NotImplementedError("calc_reduced_density is not implemented for this backend");
    };
    /// Using this backend, calculate the reciprocal reduced temperature (Tc/T)
    virtual CoolPropDbl calc_reciprocal_reduced_temperature(void) {
        throw NotImplementedError("calc_reciprocal_reduced_temperature is not implemented for this backend");
    };

    /// Using this backend, calculate the second virial coefficient
    virtual CoolPropDbl calc_Bvirial(void) {
        throw NotImplementedError("calc_Bvirial is not implemented for this backend");
    };
    /// Using this backend, calculate the third virial coefficient
    virtual CoolPropDbl calc_Cvirial(void) {
        throw NotImplementedError("calc_Cvirial is not implemented for this backend");
    };
    /// Using this backend, calculate the derivative dB/dT
    virtual CoolPropDbl calc_dBvirial_dT(void) {
        throw NotImplementedError("calc_dBvirial_dT is not implemented for this backend");
    };
    /// Using this backend, calculate the derivative dC/dT
    virtual CoolPropDbl calc_dCvirial_dT(void) {
        throw NotImplementedError("calc_dCvirial_dT is not implemented for this backend");
    };
    /// Using this backend, calculate the compressibility factor Z \f$ Z = p/(\rho R T) \f$
    virtual CoolPropDbl calc_compressibility_factor(void) {
        throw NotImplementedError("calc_compressibility_factor is not implemented for this backend");
    };

    /// Using this backend, get the name of the fluid
    virtual std::string calc_name(void) {
        throw NotImplementedError("calc_name is not implemented for this backend");
    };
    /// Using this backend, get the description of the fluid
    virtual std::string calc_description(void) {
        throw NotImplementedError("calc_description is not implemented for this backend");
    };

    /// Using this backend, get the triple point temperature in K
    virtual CoolPropDbl calc_Ttriple(void) {
        throw NotImplementedError("calc_Ttriple is not implemented for this backend");
    };
    /// Using this backend, get the triple point pressure in Pa
    virtual CoolPropDbl calc_p_triple(void) {
        throw NotImplementedError("calc_p_triple is not implemented for this backend");
    };

    /// Using this backend, get the critical point temperature in K
    virtual CoolPropDbl calc_T_critical(void) {
        throw NotImplementedError("calc_T_critical is not implemented for this backend");
    };
    /// Using this backend, get the reducing point temperature in K
    virtual CoolPropDbl calc_T_reducing(void) {
        throw NotImplementedError("calc_T_reducing is not implemented for this backend");
    };
    /// Using this backend, get the critical point pressure in Pa
    virtual CoolPropDbl calc_p_critical(void) {
        throw NotImplementedError("calc_p_critical is not implemented for this backend");
    };
    /// Using this backend, get the reducing point pressure in Pa
    virtual CoolPropDbl calc_p_reducing(void) {
        throw NotImplementedError("calc_p_reducing is not implemented for this backend");
    };
    /// Using this backend, get the critical point molar density in mol/m^3
    virtual CoolPropDbl calc_rhomolar_critical(void) {
        throw NotImplementedError("calc_rhomolar_critical is not implemented for this backend");
    };
    /// Using this backend, get the critical point mass density in kg/m^3 - Added for IF97Backend which is mass based
    virtual CoolPropDbl calc_rhomass_critical(void) {
        throw NotImplementedError("calc_rhomass_critical is not implemented for this backend");
    };
    /// Using this backend, get the reducing point molar density in mol/m^3
    virtual CoolPropDbl calc_rhomolar_reducing(void) {
        throw NotImplementedError("calc_rhomolar_reducing is not implemented for this backend");
    };
    /// Using this backend, construct the phase envelope, the variable type describes the type of phase envelope to be built.
    virtual void calc_phase_envelope(const std::string& type) {
        throw NotImplementedError("calc_phase_envelope is not implemented for this backend");
    };
    ///
    virtual CoolPropDbl calc_rhomass(void) {
        return rhomolar() * molar_mass();
    }
    virtual CoolPropDbl calc_hmass(void) {
        return hmolar() / molar_mass();
    }
    virtual CoolPropDbl calc_hmass_excess(void) {
        return hmolar_excess() / molar_mass();
    }
    virtual CoolPropDbl calc_smass(void) {
        return smolar() / molar_mass();
    }
    virtual CoolPropDbl calc_smass_excess(void) {
        return smolar_excess() / molar_mass();
    }
    virtual CoolPropDbl calc_cpmass(void) {
        return cpmolar() / molar_mass();
    }
    virtual CoolPropDbl calc_cp0mass(void) {
        return cp0molar() / molar_mass();
    }
    virtual CoolPropDbl calc_cvmass(void) {
        return cvmolar() / molar_mass();
    }
    virtual CoolPropDbl calc_umass(void) {
        return umolar() / molar_mass();
    }
    virtual CoolPropDbl calc_umass_excess(void) {
        return umolar_excess() / molar_mass();
    }
    virtual CoolPropDbl calc_gibbsmass(void) {
        return gibbsmolar() / molar_mass();
    }
    virtual CoolPropDbl calc_gibbsmass_excess(void) {
        return gibbsmolar_excess() / molar_mass();
    }
    virtual CoolPropDbl calc_helmholtzmass(void) {
        return helmholtzmolar() / molar_mass();
    }
    virtual CoolPropDbl calc_helmholtzmass_excess(void) {
        return helmholtzmolar_excess() / molar_mass();
    }
    virtual CoolPropDbl calc_volumemass_excess(void) {
        return volumemolar_excess() / molar_mass();
    }

    /// Update the states after having changed the reference state for enthalpy and entropy
    virtual void update_states(void) {
        throw NotImplementedError("This backend does not implement update_states function");
    };

    virtual CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value) {
        throw NotImplementedError("This backend does not implement calc_melting_line function");
    };

    /// @param param The key for the parameter to be returned
    /// @param Q The quality for the parameter that is given (0 = saturated liquid, 1 = saturated vapor)
    /// @param given The key for the parameter that is given
    /// @param value The value for the parameter that is given
    virtual CoolPropDbl calc_saturation_ancillary(parameters param, int Q, parameters given, double value) {
        throw NotImplementedError("This backend does not implement calc_saturation_ancillary");
    };

    /// Using this backend, calculate the phase
    virtual phases calc_phase(void) {
        throw NotImplementedError("This backend does not implement calc_phase function");
    };
    /// Using this backend, specify the phase to be used for all further calculations
    virtual void calc_specify_phase(phases phase) {
        throw NotImplementedError("This backend does not implement calc_specify_phase function");
    };
    /// Using this backend, unspecify the phase
    virtual void calc_unspecify_phase(void) {
        throw NotImplementedError("This backend does not implement calc_unspecify_phase function");
    };
    /// Using this backend, get a vector of fluid names
    virtual std::vector<std::string> calc_fluid_names(void) {
        throw NotImplementedError("This backend does not implement calc_fluid_names function");
    };
    /// Using this backend, calculate a phase given by the state string
    /// @param state A string that describes the state desired, one of "hs_anchor", "critical"/"crit", "reducing"
    virtual const CoolProp::SimpleState& calc_state(const std::string& state) {
        throw NotImplementedError("calc_state is not implemented for this backend");
    };

    virtual const CoolProp::PhaseEnvelopeData& calc_phase_envelope_data(void) {
        throw NotImplementedError("calc_phase_envelope_data is not implemented for this backend");
    };

    virtual std::vector<CoolPropDbl> calc_mole_fractions_liquid(void) {
        throw NotImplementedError("calc_mole_fractions_liquid is not implemented for this backend");
    };
    virtual std::vector<CoolPropDbl> calc_mole_fractions_vapor(void) {
        throw NotImplementedError("calc_mole_fractions_vapor is not implemented for this backend");
    };
    virtual const std::vector<CoolPropDbl> calc_mass_fractions(void) {
        throw NotImplementedError("calc_mass_fractions is not implemented for this backend");
    };

    /// Get the minimum fraction (mole, mass, volume) for incompressible fluid
    virtual CoolPropDbl calc_fraction_min(void) {
        throw NotImplementedError("calc_fraction_min is not implemented for this backend");
    };
    /// Get the maximum fraction (mole, mass, volume) for incompressible fluid
    virtual CoolPropDbl calc_fraction_max(void) {
        throw NotImplementedError("calc_fraction_max is not implemented for this backend");
    };
    virtual CoolPropDbl calc_T_freeze(void) {
        throw NotImplementedError("calc_T_freeze is not implemented for this backend");
    };

    virtual CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1) {
        throw NotImplementedError("calc_first_saturation_deriv is not implemented for this backend");
    };
    virtual CoolPropDbl calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Wrt2) {
        throw NotImplementedError("calc_second_saturation_deriv is not implemented for this backend");
    };
    virtual CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
        throw NotImplementedError("calc_first_two_phase_deriv is not implemented for this backend");
    };
    virtual CoolPropDbl calc_second_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant, parameters Wrt2, parameters Constant2) {
        throw NotImplementedError("calc_second_two_phase_deriv is not implemented for this backend");
    };
    virtual CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) {
        throw NotImplementedError("calc_first_two_phase_deriv_splined is not implemented for this backend");
    };

    virtual CoolPropDbl calc_saturated_liquid_keyed_output(parameters key) {
        throw NotImplementedError("calc_saturated_liquid_keyed_output is not implemented for this backend");
    };
    virtual CoolPropDbl calc_saturated_vapor_keyed_output(parameters key) {
        throw NotImplementedError("calc_saturated_vapor_keyed_output is not implemented for this backend");
    };
    virtual void calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) {
        throw NotImplementedError("calc_ideal_curve is not implemented for this backend");
    };

    /// Using this backend, get the temperature
    virtual CoolPropDbl calc_T(void) {
        return _T;
    }
    /// Using this backend, get the molar density in mol/m^3
    virtual CoolPropDbl calc_rhomolar(void) {
        return _rhomolar;
    }

    /// Using this backend, calculate the tangent plane distance for a given trial composition
    virtual double calc_tangent_plane_distance(const double T, const double p, const std::vector<double>& w, const double rhomolar_guess) {
        throw NotImplementedError("calc_tangent_plane_distance is not implemented for this backend");
    };

    /// Using this backend, return true critical point where dp/drho|T = 0 and d2p/drho^2|T = 0
    virtual void calc_true_critical_point(double& T, double& rho) {
        throw NotImplementedError("calc_true_critical_point is not implemented for this backend");
    };

    virtual void calc_conformal_state(const std::string& reference_fluid, CoolPropDbl& T, CoolPropDbl& rhomolar) {
        throw NotImplementedError("calc_conformal_state is not implemented for this backend");
    };

    virtual void calc_viscosity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) {
        throw NotImplementedError("calc_viscosity_contributions is not implemented for this backend");
    };
    virtual void calc_conductivity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) {
        throw NotImplementedError("calc_conductivity_contributions is not implemented for this backend");
    };
    virtual std::vector<CriticalState> calc_all_critical_points(void) {
        throw NotImplementedError("calc_all_critical_points is not implemented for this backend");
    };
    virtual void calc_build_spinodal() {
        throw NotImplementedError("calc_build_spinodal is not implemented for this backend");
    };
    virtual SpinodalData calc_get_spinodal_data() {
        throw NotImplementedError("calc_get_spinodal_data is not implemented for this backend");
    };
    virtual void calc_criticality_contour_values(double& L1star, double& M1star) {
        throw NotImplementedError("calc_criticality_contour_values is not implemented for this backend");
    };

    /// Convert mass-based input pair to molar-based input pair;  If molar-based, do nothing
    virtual void mass_to_molar_inputs(CoolProp::input_pairs& input_pair, CoolPropDbl& value1, CoolPropDbl& value2);

    /// Change the equation of state for a given component to a specified EOS
    virtual void calc_change_EOS(const std::size_t i, const std::string& EOS_name) {
        throw NotImplementedError("calc_change_EOS is not implemented for this backend");
    };

   public:
    AbstractState() : _fluid_type(FLUID_TYPE_UNDEFINED), _phase(iphase_unknown) {
        clear();
    }
    virtual ~AbstractState(){};

    /// A factory function to return a pointer to a new-allocated instance of one of the backends.
    /**
     * @brief This is a convenience function to allow for the use of '&' delimited fluid names.  Slightly less computationally efficient than the
     * @param backend The backend in use, one of "HEOS", "REFPROP", etc.
     * @param fluid_names Fluid names as a '&' delimited string
     * @return
     */
    static AbstractState* factory(const std::string& backend, const std::string& fluid_names) {
        return factory(backend, strsplit(fluid_names, '&'));
    };

    /**
     * @brief A factory function to return a pointer to a new-allocated instance of one of the backends.
     * @param backend The backend in use, "HEOS", "REFPROP", etc.
     * @param fluid_names A vector of strings of the fluid names
     * @return A pointer to the instance generated
     *
     * Several backends are possible:
     *
     * 1. "?" : The backend is unknown, we will parse the fluid string to determine the backend to be used.  Probably will use HEOS backend (see below)
     * 2. "HEOS" : The Helmholtz Equation of State backend for use with pure and pseudo-pure fluids, and mixtures, all of which are based on multi-parameter Helmholtz Energy equations of state.  The fluid part of the string should then either be
     *    1. A pure or pseudo-pure fluid name (eg. "PROPANE" or "R410A"), yielding a HelmholtzEOSBackend instance.
     *    2. A string that encodes the components of the mixture with a "&" between them (e.g. "R32&R125"), yielding a HelmholtzEOSMixtureBackend instance.
     *
     * 3. "REFPROP" : The REFPROP backend will be used.  The fluid part of the string should then either be
     *    1. A pure or pseudo-pure fluid name (eg. "PROPANE" or "R410A"), yielding a REFPROPBackend instance.
     *    2. A string that encodes the components of the mixture with a "&" between them (e.g. "R32&R125"), yielding a REFPROPMixtureBackend instance.
     *
     * 4. "INCOMP": The incompressible backend will be used
     * 5. "TTSE&XXXX": The TTSE backend will be used, and the tables will be generated using the XXXX backend where XXXX is one of the base backends("HEOS", "REFPROP", etc. )
     * 6. "BICUBIC&XXXX": The Bicubic backend will be used, and the tables will be generated using the XXXX backend where XXXX is one of the base backends("HEOS", "REFPROP", etc. )
     *
     * Very Important!! : Use a smart pointer to manage the pointer returned.  In older versions of C++, you can use std::tr1::smart_ptr. In C++2011 you can use std::shared_ptr
     */
    static AbstractState* factory(const std::string& backend, const std::vector<std::string>& fluid_names);

    /// Set the internal variable T without a flash call (expert use only!)
    void set_T(CoolPropDbl T) {
        _T = T;
    }

    /// Get a string representation of the backend - for instance "HelmholtzEOSMixtureBackend"
    /// for the core mixture model in CoolProp
    ///
    /// Must be overloaded by the backend to provide the backend's name
    virtual std::string backend_name(void) = 0;

    // The derived classes must implement this function to define whether they use mole fractions (true) or mass fractions (false)
    virtual bool using_mole_fractions(void) = 0;
    virtual bool using_mass_fractions(void) = 0;
    virtual bool using_volu_fractions(void) = 0;

    virtual void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) = 0;
    virtual void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) = 0;
    virtual void set_volu_fractions(const std::vector<CoolPropDbl>& mass_fractions) {
        throw NotImplementedError("Volume composition has not been implemented.");
    }

    /**
    \brief Set the reference state based on a string representation

    @param reference_state The reference state to use, one of

    Reference State | Description
    -------------   | -------------------
    "IIR"           | h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
    "ASHRAE"        | h = 0, s = 0 @ -40C saturated liquid
    "NBP"           | h = 0, s = 0 @ 1.0 bar saturated liquid
    "DEF"           | Reset to the default reference state for the fluid
    "RESET"         | Remove the offset

    The offset in the ideal gas Helmholtz energy can be obtained from
    \f[
    \displaystyle\frac{\Delta s}{R_u/M}+\frac{\Delta h}{(R_u/M)T}\tau
    \f]
    where \f$ \Delta s = s-s_{spec} \f$ and \f$ \Delta h = h-h_{spec} \f$
    */
    virtual void set_reference_stateS(const std::string& reference_state) {
        throw NotImplementedError(
          "Setting reference state has not been implemented for this backend. Try using CoolProp::set_reference_stateD instead.");
    }

    /// Set the reference state based on a thermodynamic state point specified by temperature and molar density
    /// @param T Temperature at reference state [K]
    /// @param rhomolar Molar density at reference state [mol/m^3]
    /// @param hmolar0 Molar enthalpy at reference state [J/mol]
    /// @param smolar0 Molar entropy at reference state [J/mol/K]
    virtual void set_reference_stateD(double T, double rhomolar, double hmolar0, double smolar0) {
        throw NotImplementedError(
          "Setting reference state has not been implemented for this backend. Try using CoolProp::set_reference_stateD instead.");
    }

#ifndef COOLPROPDBL_MAPS_TO_DOUBLE
    void set_mole_fractions(const std::vector<double>& mole_fractions) {
        set_mole_fractions(std::vector<CoolPropDbl>(mole_fractions.begin(), mole_fractions.end()));
    };
    void set_mass_fractions(const std::vector<double>& mass_fractions) {
        set_mass_fractions(std::vector<CoolPropDbl>(mass_fractions.begin(), mass_fractions.end()));
    };
    void set_volu_fractions(const std::vector<double>& volu_fractions) {
        set_volu_fractions(std::vector<CoolPropDbl>(volu_fractions.begin(), volu_fractions.end()));
    };
#endif

#ifdef EMSCRIPTEN
    void set_mole_fractions_double(const std::vector<double>& mole_fractions) {
        set_mole_fractions(std::vector<CoolPropDbl>(mole_fractions.begin(), mole_fractions.end()));
    };
#endif

    /// Get the mole fractions of the equilibrium liquid phase
    std::vector<CoolPropDbl> mole_fractions_liquid(void) {
        return calc_mole_fractions_liquid();
    };
    /// Get the mole fractions of the equilibrium liquid phase (but as a double for use in SWIG wrapper)
    std::vector<double> mole_fractions_liquid_double(void) {
        std::vector<CoolPropDbl> x = calc_mole_fractions_liquid();
        return std::vector<double>(x.begin(), x.end());
    };

    /// Get the mole fractions of the equilibrium vapor phase
    std::vector<CoolPropDbl> mole_fractions_vapor(void) {
        return calc_mole_fractions_vapor();
    };
    /// Get the mole fractions of the equilibrium vapor phase (but as a double for use in SWIG wrapper)
    std::vector<double> mole_fractions_vapor_double(void) {
        std::vector<CoolPropDbl> y = calc_mole_fractions_vapor();
        return std::vector<double>(y.begin(), y.end());
    };

    /// Get the mole fractions of the fluid
    virtual const std::vector<CoolPropDbl>& get_mole_fractions(void) = 0;
    /// Get the mass fractions of the fluid
    virtual const std::vector<CoolPropDbl> get_mass_fractions(void) {
        return this->calc_mass_fractions();
    };

    /// Update the state using two state variables
    virtual void update(CoolProp::input_pairs input_pair, double Value1, double Value2) = 0;

    /// Update the state using two state variables and providing guess values
    /// Some or all of the guesses will be used - this is backend dependent
    virtual void update_with_guesses(CoolProp::input_pairs input_pair, double Value1, double Value2, const GuessesStructure& guesses) {
        throw NotImplementedError("update_with_guesses is not implemented for this backend");
    };

    /// A function that says whether the backend instance can be instantiated in the high-level interface
    /// In general this should be true, except for some other backends (especially the tabular backends)
    /// To disable use in high-level interface, implement this function and return false
    virtual bool available_in_high_level(void) {
        return true;
    }

    /// Return a string from the backend for the mixture/fluid - backend dependent - could be CAS #, name, etc.
    virtual std::string fluid_param_string(const std::string&) {
        throw NotImplementedError("fluid_param_string has not been implemented for this backend");
    }

    /// Return a vector of strings of the fluid names that are in use
    std::vector<std::string> fluid_names(void);

    /** Get a constant for one of the fluids forming this mixture
     *  @param i Index (0-based) of the fluid
     *  @param param parameter you want to obtain (probably one that is a trivial parameter)
     */
    virtual const double get_fluid_constant(std::size_t i, parameters param) const {
        throw NotImplementedError("get_fluid_constant is not implemented for this backend");
    };
    ;

    /// Set binary mixture floating point parameter (EXPERT USE ONLY!!!)
    virtual void set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter, const double value) {
        throw NotImplementedError("set_binary_interaction_double is not implemented for this backend");
    };
    /// Set binary mixture floating point parameter (EXPERT USE ONLY!!!)
    virtual void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, const double value) {
        throw NotImplementedError("set_binary_interaction_double is not implemented for this backend");
    };
    /// Set binary mixture string parameter (EXPERT USE ONLY!!!)
    virtual void set_binary_interaction_string(const std::string& CAS1, const std::string& CAS2, const std::string& parameter,
                                               const std::string& value) {
        throw NotImplementedError("set_binary_interaction_string is not implemented for this backend");
    };
    /// Set binary mixture string parameter (EXPERT USE ONLY!!!)
    virtual void set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter, const std::string& value) {
        throw NotImplementedError("set_binary_interaction_string is not implemented for this backend");
    };
    /// Get binary mixture double value (EXPERT USE ONLY!!!)
    virtual double get_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) {
        throw NotImplementedError("get_binary_interaction_double is not implemented for this backend");
    };
    /// Get binary mixture double value (EXPERT USE ONLY!!!)
    virtual double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) {
        throw NotImplementedError("get_binary_interaction_double is not implemented for this backend");
    };
    /// Get binary mixture string value (EXPERT USE ONLY!!!)
    virtual std::string get_binary_interaction_string(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) {
        throw NotImplementedError("get_binary_interaction_string is not implemented for this backend");
    };
    /// Apply a simple mixing rule (EXPERT USE ONLY!!!)
    virtual void apply_simple_mixing_rule(std::size_t i, std::size_t j, const std::string& model) {
        throw NotImplementedError("apply_simple_mixing_rule is not implemented for this backend");
    };
    /// Set the cubic alpha function's constants:
    virtual void set_cubic_alpha_C(const size_t i, const std::string& parameter, const double c1, const double c2, const double c3) {
        throw ValueError("set_cubic_alpha_C only defined for cubic backends");
    };
    /// Set fluid parameter (currently the volume translation parameter for cubic)
    virtual void set_fluid_parameter_double(const size_t i, const std::string& parameter, const double value) {
        throw ValueError("set_fluid_parameter_double only defined for cubic backends");
    };
    /// Double fluid parameter (currently the volume translation parameter for cubic)
    virtual double get_fluid_parameter_double(const size_t i, const std::string& parameter) {
        throw ValueError("get_fluid_parameter_double only defined for cubic backends");
    };

    /// Clear all the cached values
    virtual bool clear();
    /// When the composition changes, clear all cached values that are only dependent on composition, but not the thermodynamic state
    virtual bool clear_comp_change();

    /// Get the state that is used in the equation of state or mixture model
    /// to reduce the state.  For pure fluids this is usually, but not always,
    /// the critical point.  For mixture models, it is usually composition dependent
    virtual const CoolProp::SimpleState& get_reducing_state() {
        return _reducing;
    };

    /// Get a desired state point - backend dependent
    const CoolProp::SimpleState& get_state(const std::string& state) {
        return calc_state(state);
    };

    /// Get the minimum temperature in K
    double Tmin(void);
    /// Get the maximum temperature in K
    double Tmax(void);
    /// Get the maximum pressure in Pa
    double pmax(void);
    /// Get the triple point temperature in K
    double Ttriple(void);

    /// Get the phase of the state
    phases phase(void) {
        return calc_phase();
    };
    /// Specify the phase for all further calculations with this state class
    void specify_phase(phases phase) {
        calc_specify_phase(phase);
    };
    /// Unspecify the phase and go back to calculating it based on the inputs
    void unspecify_phase(void) {
        calc_unspecify_phase();
    };

    /// Return the critical temperature in K
    double T_critical(void);
    /// Return the critical pressure in Pa
    double p_critical(void);
    /// Return the critical molar density in mol/m^3
    double rhomolar_critical(void);
    /// Return the critical mass density in kg/m^3
    double rhomass_critical(void);

    /// Return the vector of critical points, including points that are unstable or correspond to negative pressure
    std::vector<CriticalState> all_critical_points(void) {
        return calc_all_critical_points();
    };

    /// Construct the spinodal curve for the mixture (or pure fluid)
    void build_spinodal() {
        calc_build_spinodal();
    };

    /// Get the data from the spinodal curve constructed in the call to build_spinodal()
    SpinodalData get_spinodal_data() {
        return calc_get_spinodal_data();
    };

    /// Calculate the criticality contour values \f$\mathcal{L}_1^*\f$ and \f$\mathcal{M}_1^*\f$
    void criticality_contour_values(double& L1star, double& M1star) {
        return calc_criticality_contour_values(L1star, M1star);
    }

    /// Return the tangent plane distance for a given trial composition w
    /// @param T Temperature (K)
    /// @param p Pressure (Pa)
    /// @param w The trial composition
    /// @param rhomolar_guess (mol/m^3) The molar density guess value (if <0 (default), not used; if >0, guess value will be used in flash evaluation)
    ///
    /// \f[
    /// tpd(w) = \sum_i w_i(\ln w_i + \ln \phi_i(w) - d_i)
    /// \f]
    /// with
    /// \f[ d_i = \ln z_i + \ln \phi_i(z) \f]
    /// Or you can express the \f$ tpd \f$ in terms of fugacity (See Table 7.3 from GERG 2004 monograph)
    /// since \f$ \ln \phi_i = \ln f_i - \ln p -\ln z_i\f$
    /// thus
    /// \f[ d_i = \ln f_i(z) - \ln p\f]
    /// and
    /// \f[
    /// tpd(w) = \sum_i w_i(\ln f_i(w) - \ln p - d_i)
    /// \f]
    /// and the \f$ \ln p \f$ cancel, leaving
    /// \f[
    /// tpd(w) = \sum_i w_i(\ln f_i(w) - \ln f_i(z))
    /// \f]
    double tangent_plane_distance(const double T, const double p, const std::vector<double>& w, const double rhomolar_guess = -1) {
        return calc_tangent_plane_distance(T, p, w, rhomolar_guess);
    };

    /// Return the reducing point temperature in K
    double T_reducing(void);
    /// Return the molar density at the reducing point in mol/m^3
    double rhomolar_reducing(void);
    /// Return the mass density at the reducing point in kg/m^3
    double rhomass_reducing(void);

    /// Return the triple point pressure in Pa
    double p_triple(void);

    /// Return the name - backend dependent
    std::string name() {
        return calc_name();
    };
    /// Return the description - backend dependent
    std::string description() {
        return calc_description();
    };

    /// Return the dipole moment in C-m (1 D = 3.33564e-30 C-m)
    double dipole_moment() {
        return calc_dipole_moment();
    }

    // ----------------------------------------
    // Bulk properties - temperature and density are directly calculated every time
    // All other parameters are calculated on an as-needed basis
    // ----------------------------------------
    /// Retrieve a value by key
    double keyed_output(parameters key);
    /// A trivial keyed output like molar mass that does not depend on the state
    double trivial_keyed_output(parameters key);
    /// Get an output from the saturated liquid state by key
    double saturated_liquid_keyed_output(parameters key) {
        return calc_saturated_liquid_keyed_output(key);
    };
    /// Get an output from the saturated vapor state by key
    double saturated_vapor_keyed_output(parameters key) {
        return calc_saturated_vapor_keyed_output(key);
    };

    /// Return the temperature in K
    double T(void) {
        return calc_T();
    };
    /// Return the molar density in mol/m^3
    double rhomolar(void) {
        return calc_rhomolar();
    };
    /// Return the mass density in kg/m^3
    double rhomass(void) {
        return calc_rhomass();
    };
    /// Return the pressure in Pa
    double p(void) {
        return _p;
    };
    /// Return the vapor quality (mol/mol); Q = 0 for saturated liquid
    double Q(void) {
        return _Q;
    };
    /// Return the reciprocal of the reduced temperature (\f$\tau = T_c/T\f$)
    double tau(void);
    /// Return the reduced density (\f$\delta = \rho/\rho_c\f$)
    double delta(void);
    /// Return the molar mass in kg/mol
    double molar_mass(void);
    /// Return the acentric factor
    double acentric_factor(void);
    /// Return the mole-fraction weighted gas constant in J/mol/K
    double gas_constant(void);
    /// Return the B virial coefficient
    double Bvirial(void);
    /// Return the derivative of the B virial coefficient with respect to temperature
    double dBvirial_dT(void);
    /// Return the C virial coefficient
    double Cvirial(void);
    /// Return the derivative of the C virial coefficient with respect to temperature
    double dCvirial_dT(void);
    /// Return the compressibility factor \f$ Z = p/(rho R T) \f$
    double compressibility_factor(void);
    /// Return the molar enthalpy in J/mol
    double hmolar(void);
    /// Return the residual molar enthalpy in J/mol
    double hmolar_residual(void);
    /// Return the mass enthalpy in J/kg
    double hmass(void) {
        return calc_hmass();
    };
    /// Return the excess molar enthalpy in J/mol
    double hmolar_excess(void);
    /// Return the excess mass enthalpy in J/kg
    double hmass_excess(void) {
        return calc_hmass_excess();
    };
    /// Return the molar entropy in J/mol/K
    double smolar(void);
    /// Return the residual molar entropy (as a function of temperature and density) in J/mol/K
    double smolar_residual(void);
    /// Return the effective hardness of interaction
    double neff(void);
    /// Return the molar entropy in J/kg/K
    double smass(void) {
        return calc_smass();
    };
    /// Return the molar entropy in J/mol/K
    double smolar_excess(void);
    /// Return the molar entropy in J/kg/K
    double smass_excess(void) {
        return calc_smass_excess();
    };
    /// Return the molar internal energy in J/mol
    double umolar(void);
    /// Return the mass internal energy in J/kg
    double umass(void) {
        return calc_umass();
    };
    /// Return the excess internal energy in J/mol
    double umolar_excess(void);
    /// Return the excess internal energy in J/kg
    double umass_excess(void) {
        return calc_umass_excess();
    };
    /// Return the molar constant pressure specific heat in J/mol/K
    double cpmolar(void);
    /// Return the mass constant pressure specific heat in J/kg/K
    double cpmass(void) {
        return calc_cpmass();
    };
    /// Return the molar constant pressure specific heat for ideal gas part only in J/mol/K
    double cp0molar(void);
    /// Return the mass constant pressure specific heat for ideal gas part only in J/kg/K
    double cp0mass(void) {
        return calc_cp0mass();
    };
    /// Return the molar constant volume specific heat in J/mol/K
    double cvmolar(void);
    /// Return the mass constant volume specific heat in J/kg/K
    double cvmass(void) {
        return calc_cvmass();
    };
    /// Return the Gibbs energy in J/mol
    double gibbsmolar(void);
    /// Return the residual Gibbs energy in J/mol
    double gibbsmolar_residual(void);
    /// Return the Gibbs energy in J/kg
    double gibbsmass(void) {
        return calc_gibbsmass();
    };
    /// Return the excess Gibbs energy in J/mol
    double gibbsmolar_excess(void);
    /// Return the excess Gibbs energy in J/kg
    double gibbsmass_excess(void) {
        return calc_gibbsmass_excess();
    };
    /// Return the Helmholtz energy in J/mol
    double helmholtzmolar(void);
    /// Return the Helmholtz energy in J/kg
    double helmholtzmass(void) {
        return calc_helmholtzmass();
    };
    /// Return the excess Helmholtz energy in J/mol
    double helmholtzmolar_excess(void);
    /// Return the excess Helmholtz energy in J/kg
    double helmholtzmass_excess(void) {
        return calc_helmholtzmass_excess();
    };
    /// Return the excess volume in m^3/mol
    double volumemolar_excess(void);
    /// Return the excess volume in m^3/kg
    double volumemass_excess(void) {
        return calc_volumemass_excess();
    };
    /// Return the speed of sound in m/s
    double speed_sound(void);
    /// Return the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    double isothermal_compressibility(void);
    /// Return the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    double isobaric_expansion_coefficient(void);
    /// Return the isentropic expansion coefficient \f$ \kappa_s = -\frac{c_p}{c_v}\frac{v}{p}\left.\frac{\partial p}{\partial v}\right|_T = \frac{\rho}{p}\left.\frac{\partial p}{\partial \rho}\right|_s\f$
    double isentropic_expansion_coefficient(void);
    /// Return the fugacity coefficient of the i-th component of the mixture
    double fugacity_coefficient(std::size_t i);
    /// Return a vector of the fugacity coefficients for all components in the mixture
    std::vector<double> fugacity_coefficients();
    /// Return the fugacity of the i-th component of the mixture
    double fugacity(std::size_t i);
    /// Return the chemical potential of the i-th component of the mixture
    double chemical_potential(std::size_t i);
    /** \brief Return the fundamental derivative of gas dynamics \f$ \Gamma \f$
     *
     * see also Colonna et al, FPE, 2010
     *
     * \f[ \Gamma = 1+\frac{\rho}{c}\left(\frac{\partial c}{\partial \rho}\right)_{s} = 1+\frac{\rho}{2c^2}\left(\frac{\partial^2 p}{\partial \rho^2}\right)_{s} = \frac{v^3}{2c^2}\left(\frac{\partial^2 p}{\partial v^2}\right)_{s}\f]
     *
     * Note: densities are mass-based densities, not mole-based densities
     */
    double fundamental_derivative_of_gas_dynamics(void);
    /// Return the phase identification parameter (PIP) of G. Venkatarathnam and L.R. Oellrich, "Identification of the phase of a fluid using partial derivatives of pressure, volume, and temperature without reference to saturation properties: Applications in phase equilibria calculations"
    double PIP() {
        return calc_PIP();
    };

    /// Calculate the "true" critical point for pure fluids where dpdrho|T and d2p/drho2|T are equal to zero
    void true_critical_point(double& T, double& rho) {
        calc_true_critical_point(T, rho);
    }

    /**
     * \brief Calculate an ideal curve for a pure fluid
     *
     * @param type The type of ideal curve you would like to calculate - "Ideal", "Boyle", "Joule-Thomson", "Joule Inversion", etc.
     * @param T The temperatures along the curve in K
     * @param p The pressures along the curve in Pa
    */
    void ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) {
        calc_ideal_curve(type, T, p);
    };

    // ----------------------------------------
    //    Partial derivatives
    // ----------------------------------------

    /** \brief The first partial derivative in homogeneous phases
     *
     * \f[ \left(\frac{\partial A}{\partial B}\right)_C = \frac{\left(\frac{\partial A}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial A}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta}{\left(\frac{\partial B}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial B}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta} = \frac{N}{D}\f]
     */
    CoolPropDbl first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
        return calc_first_partial_deriv(Of, Wrt, Constant);
    };

    /** \brief The second partial derivative in homogeneous phases
     *
     * The first partial derivative (\ref CoolProp::AbstractState::first_partial_deriv) can be expressed as
     *
     * \f[ \left(\frac{\partial A}{\partial B}\right)_C = \frac{\left(\frac{\partial A}{\partial T}\right)_\rho\left(\frac{\partial C}{\partial \rho}\right)_T-\left(\frac{\partial A}{\partial \rho}\right)_T\left(\frac{\partial C}{\partial T}\right)_\rho}{\left(\frac{\partial B}{\partial T}\right)_\rho\left(\frac{\partial C}{\partial \rho}\right)_T-\left(\frac{\partial B}{\partial \rho}\right)_T\left(\frac{\partial C}{\partial T}\right)_\rho} = \frac{N}{D}\f]
     *
     * and the second derivative can be expressed as
     *
     * \f[
     * \frac{\partial}{\partial D}\left(\left(\frac{\partial A}{\partial B}\right)_C\right)_E = \frac{\frac{\partial}{\partial T}\left( \left(\frac{\partial A}{\partial B}\right)_C \right)_\rho\left(\frac{\partial E}{\partial \rho}\right)_T-\frac{\partial}{\partial \rho}\left(\left(\frac{\partial A}{\partial B}\right)_C\right)_T\left(\frac{\partial E}{\partial T}\right)_\rho}{\left(\frac{\partial D}{\partial T}\right)_\rho\left(\frac{\partial E}{\partial \rho}\right)_T-\left(\frac{\partial D}{\partial \rho}\right)_T\left(\frac{\partial E}{\partial T}\right)_\rho}
     * \f]
     *
     * which can be expressed in parts as
     *
     * \f[\left(\frac{\partial N}{\partial \rho}\right)_{T} = \left(\frac{\partial A}{\partial T}\right)_\rho\left(\frac{\partial^2 C}{\partial \rho^2}\right)_{T}+\left(\frac{\partial^2 A}{\partial T\partial\rho}\right)\left(\frac{\partial C}{\partial \rho}\right)_{T}-\left(\frac{\partial A}{\partial \rho}\right)_T\left(\frac{\partial^2 C}{\partial T\partial\rho}\right)-\left(\frac{\partial^2 A}{\partial \rho^2}\right)_{T}\left(\frac{\partial C}{\partial T}\right)_\rho\f]
     * \f[\left(\frac{\partial D}{\partial \rho}\right)_{T} = \left(\frac{\partial B}{\partial T}\right)_\rho\left(\frac{\partial^2 C}{\partial \rho^2}\right)_{T}+\left(\frac{\partial^2 B}{\partial T\partial\rho}\right)\left(\frac{\partial C}{\partial \rho}\right)_{T}-\left(\frac{\partial B}{\partial \rho}\right)_T\left(\frac{\partial^2 C}{\partial T\partial\rho}\right)-\left(\frac{\partial^2 B}{\partial \rho^2}\right)_{T}\left(\frac{\partial C}{\partial T}\right)_\rho\f]
     * \f[\left(\frac{\partial N}{\partial T}\right)_{\rho} = \left(\frac{\partial A}{\partial T}\right)_\rho\left(\frac{\partial^2 C}{\partial \rho\partial T}\right)+\left(\frac{\partial^2 A}{\partial T^2}\right)_\rho\left(\frac{\partial C}{\partial \rho}\right)_{T}-\left(\frac{\partial A}{\partial \rho}\right)_T\left(\frac{\partial^2 C}{\partial T^2}\right)_\rho-\left(\frac{\partial^2 A}{\partial \rho\partial T}\right)\left(\frac{\partial C}{\partial T}\right)_\rho\f]
     * \f[\left(\frac{\partial D}{\partial T}\right)_{\rho} = \left(\frac{\partial B}{\partial T}\right)_\rho\left(\frac{\partial^2 C}{\partial \rho\partial T}\right)+\left(\frac{\partial^2 B}{\partial T^2}\right)_\rho\left(\frac{\partial C}{\partial \rho}\right)_{T}-\left(\frac{\partial B}{\partial \rho}\right)_T\left(\frac{\partial^2 C}{\partial T^2}\right)_\rho-\left(\frac{\partial^2 B}{\partial \rho\partial T}\right)\left(\frac{\partial C}{\partial T}\right)_\rho\f]
     * \f[\frac{\partial}{\partial \rho}\left( \left(\frac{\partial A}{\partial B}\right)_C \right)_T = \frac{D\left(\frac{\partial N}{\partial \rho}\right)_{T}-N\left(\frac{\partial D}{\partial \rho}\right)_{\tau}}{D^2}\f]
     * \f[\frac{\partial}{\partial T}\left( \left(\frac{\partial A}{\partial B}\right)_C \right)_\rho = \frac{D\left(\frac{\partial N}{\partial T}\right)_{\rho}-N\left(\frac{\partial D}{\partial T}\right)_{\rho}}{D^2}\f]
     *
     * The terms \f$ N \f$ and \f$ D \f$ are the numerator and denominator from \ref CoolProp::AbstractState::first_partial_deriv respectively
     */
    CoolPropDbl second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) {
        return calc_second_partial_deriv(Of1, Wrt1, Constant1, Wrt2, Constant2);
    };

    /** \brief The first partial derivative along the saturation curve
     *
     * Implementing the algorithms and ideas of:
     * Matthis Thorade, Ali Saadat, "Partial derivatives of thermodynamic state properties for dynamic simulation",
     * Environmental Earth Sciences, December 2013, Volume 70, Issue 8, pp 3497-3503
     *
     * Basically the idea is that the p-T derivative is given by Clapeyron relations:
     *
     * \f[ \left(\frac{\partial T}{\partial p}\right)_{\sigma} = T\left(\frac{v'' - v'}{h'' - h'}\right)_{\sigma} \f]
     *
     * and then other derivatives can be obtained along the saturation curve from
     *
     * \f[ \left(\frac{\partial y}{\partial p}\right)_{\sigma} = \left(\frac{\partial y}{\partial p}\right)+\left(\frac{\partial y}{\partial T}\right)\left(\frac{\partial T}{\partial p}\right)_{\sigma} \f]
     *
     * \f[ \left(\frac{\partial y}{\partial T}\right)_{\sigma} = \left(\frac{\partial y}{\partial T}\right)+\left(\frac{\partial y}{\partial p}\right)\left(\frac{\partial p}{\partial T}\right)_{\sigma} \f]
     *
     * where derivatives without the \f$ \sigma \f$ are homogeneous (conventional) derivatives.
     *
     * @param Of1 The parameter that the derivative is taken of
     * @param Wrt1 The parameter that the derivative is taken with respect to
     */
    CoolPropDbl first_saturation_deriv(parameters Of1, parameters Wrt1) {
        return calc_first_saturation_deriv(Of1, Wrt1);
    };

    /** \brief The second partial derivative along the saturation curve
     *
     * Implementing the algorithms and ideas of:
     * Matthis Thorade, Ali Saadat, "Partial derivatives of thermodynamic state properties for dynamic simulation",
     * Environmental Earth Sciences, December 2013, Volume 70, Issue 8, pp 3497-3503
     *
     * Like with \ref first_saturation_deriv, we can express the derivative as
     * \f[ \left(\frac{\partial y}{\partial T}\right)_{\sigma} = \left(\frac{\partial y}{\partial T}\right)+\left(\frac{\partial y}{\partial p}\right)\left(\frac{\partial p}{\partial T}\right)_{\sigma} \f]
     *
     * where \f$ y \f$ is already a saturation derivative. So you might end up with something like
     *
     * \f[ \left(\frac{\partial \left(\frac{\partial T}{\partial p}\right)_{\sigma}}{\partial T}\right)_{\sigma} = \left(\frac{\partial \left(\frac{\partial T}{\partial p}\right)_{\sigma}}{\partial T}\right)+\left(\frac{\partial \left(\frac{\partial T}{\partial p}\right)_{\sigma}}{\partial p}\right)\left(\frac{\partial p}{\partial T}\right)_{\sigma} \f]
     *
     * @param Of1 The parameter that the first derivative is taken of
     * @param Wrt1 The parameter that the first derivative is taken with respect to
     * @param Wrt2 The parameter that the second derivative is taken with respect to
     * */
    CoolPropDbl second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Wrt2) {
        return calc_second_saturation_deriv(Of1, Wrt1, Wrt2);
    };

    /**
     * @brief Calculate the first "two-phase" derivative as described by Thorade and Sadaat, EAS, 2013
     *
     * Implementing the algorithms and ideas of:
     * Matthis Thorade, Ali Saadat, "Partial derivatives of thermodynamic state properties for dynamic simulation",
     * Environmental Earth Sciences, December 2013, Volume 70, Issue 8, pp 3497-3503
     *
     * Spline evaluation is as described in:
     * S Quoilin, I Bell, A Desideri, P Dewallef, V Lemort,
     * "Methods to increase the robustness of finite-volume flow models in thermodynamic systems",
     * Energies 7 (3), 1621-1640
     *
     * \note Not all derivatives are supported!
     *
     * @param Of The parameter to be derived
     * @param Wrt The parameter that the derivative is taken with respect to
     * @param Constant The parameter that is held constant
     * @return
     */
    double first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
        return calc_first_two_phase_deriv(Of, Wrt, Constant);
    };

    /**
    * @brief Calculate the second "two-phase" derivative as described by Thorade and Sadaat, EAS, 2013
    *
    * Implementing the algorithms and ideas of:
    * Matthis Thorade, Ali Saadat, "Partial derivatives of thermodynamic state properties for dynamic simulation",
    * Environmental Earth Sciences, December 2013, Volume 70, Issue 8, pp 3497-3503
    *
    * \note Not all derivatives are supported!
    *
    * @param Of The parameter to be derived
    * @param Wrt1 The parameter that the derivative is taken with respect to in the first derivative
    * @param Constant1 The parameter that is held constant in the first derivative
    * @param Wrt2 The parameter that the derivative is taken with respect to in the second derivative
    * @param Constant2 The parameter that is held constant in the second derivative
    * @return
    */
    double second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) {
        return calc_second_two_phase_deriv(Of, Wrt1, Constant1, Wrt2, Constant2);
    };

    /**
    * @brief Calculate the first "two-phase" derivative as described by Thorade and Sadaat, EAS, 2013
    *
    * Implementing the algorithms and ideas of:
    * Matthis Thorade, Ali Saadat, "Partial derivatives of thermodynamic state properties for dynamic simulation",
    * Environmental Earth Sciences, December 2013, Volume 70, Issue 8, pp 3497-3503
    *
    * Spline evaluation is as described in:
    * S Quoilin, I Bell, A Desideri, P Dewallef, V Lemort,
    * "Methods to increase the robustness of finite-volume flow models in thermodynamic systems",
    * Energies 7 (3), 1621-1640
    *
    * \note Not all derivatives are supported! If you need all three currently supported values (drho_dh__p, drho_dp__h and rho_spline), you should calculate drho_dp__h first to avoid duplicate calculations.
    *
    * @param Of The parameter to be derived
    * @param Wrt The parameter that the derivative is taken with respect to
    * @param Constant The parameter that is held constant
    * @param x_end The end vapor quality at which the spline is defined (spline is active in [0, x_end])
    * @return
    */
    double first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, double x_end) {
        return calc_first_two_phase_deriv_splined(Of, Wrt, Constant, x_end);
    };

    // ----------------------------------------
    //    Phase envelope for mixtures
    // ----------------------------------------

    /**
     * \brief Construct the phase envelope for a mixture
     *
     * @param type currently a dummy variable that is not used
     */
    void build_phase_envelope(const std::string& type = "");
    /**
     * \brief After having calculated the phase envelope, return the phase envelope data
     */
    const CoolProp::PhaseEnvelopeData& get_phase_envelope_data() {
        return calc_phase_envelope_data();
    };

    // ----------------------------------------
    //    Ancillary equations
    // ----------------------------------------

    /// Return true if the fluid has a melting line - default is false, but can be re-implemented by derived class
    virtual bool has_melting_line(void) {
        return false;
    };
    /// Return a value from the melting line
    /// @param param The key for the parameter to be returned
    /// @param given The key for the parameter that is given
    /// @param value The value for the parameter that is given
    double melting_line(int param, int given, double value);
    /// Return the value from a saturation ancillary curve (if the backend implements it)
    /// @param param The key for the parameter to be returned
    /// @param Q The quality for the parameter that is given (0 = saturated liquid, 1 = saturated vapor)
    /// @param given The key for the parameter that is given
    /// @param value The value for the parameter that is given
    double saturation_ancillary(parameters param, int Q, parameters given, double value);

    // ----------------------------------------
    // Transport properties
    // ----------------------------------------
    /// Return the viscosity in Pa-s
    double viscosity(void);
    /// Return the viscosity contributions, each in Pa-s
    void viscosity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) {
        calc_viscosity_contributions(dilute, initial_density, residual, critical);
    };
    /// Return the thermal conductivity in W/m/K
    double conductivity(void);
    /// Return the thermal conductivity contributions, each in W/m/K
    void conductivity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) {
        calc_conductivity_contributions(dilute, initial_density, residual, critical);
    };
    /// Return the surface tension in N/m
    double surface_tension(void);
    /// Return the Prandtl number (dimensionless)
    double Prandtl(void) {
        return cpmass() * viscosity() / conductivity();
    };
    /**
     * @brief Find the conformal state needed for ECS
     * @param reference_fluid The reference fluid for which the conformal state will be calculated
     * @param T Temperature (initial guess must be provided, or < 0 to start with unity shape factors)
     * @param rhomolar Molar density (initial guess must be provided, or < 0 to start with unity shape factors)
     */
    void conformal_state(const std::string& reference_fluid, CoolPropDbl& T, CoolPropDbl& rhomolar) {
        return calc_conformal_state(reference_fluid, T, rhomolar);
    };

    /// \brief Change the equation of state for a given component to a specified EOS
    /// @param i Index of the component to change (if a pure fluid, i=0)
    /// @param EOS_name Name of the EOS to use (something like "SRK", "PR", "XiangDeiters", but backend-specific)
    /// \note Calls the calc_change_EOS function of the implementation
    void change_EOS(const std::size_t i, const std::string& EOS_name) {
        calc_change_EOS(i, EOS_name);
    }

    // ----------------------------------------
    // Helmholtz energy and derivatives
    // ----------------------------------------
    /// Return the term \f$ \alpha^0 \f$
    CoolPropDbl alpha0(void) {
        if (!_alpha0) _alpha0 = calc_alpha0();
        return _alpha0;
    };
    /// Return the term \f$ \alpha^0_{\delta} \f$
    CoolPropDbl dalpha0_dDelta(void) {
        if (!_dalpha0_dDelta) _dalpha0_dDelta = calc_dalpha0_dDelta();
        return _dalpha0_dDelta;
    };
    /// Return the term \f$ \alpha^0_{\tau} \f$
    CoolPropDbl dalpha0_dTau(void) {
        if (!_dalpha0_dTau) _dalpha0_dTau = calc_dalpha0_dTau();
        return _dalpha0_dTau;
    };
    /// Return the term \f$ \alpha^0_{\delta\delta} \f$
    CoolPropDbl d2alpha0_dDelta2(void) {
        if (!_d2alpha0_dDelta2) _d2alpha0_dDelta2 = calc_d2alpha0_dDelta2();
        return _d2alpha0_dDelta2;
    };
    /// Return the term \f$ \alpha^0_{\delta\tau} \f$
    CoolPropDbl d2alpha0_dDelta_dTau(void) {
        if (!_d2alpha0_dDelta_dTau) _d2alpha0_dDelta_dTau = calc_d2alpha0_dDelta_dTau();
        return _d2alpha0_dDelta_dTau;
    };
    /// Return the term \f$ \alpha^0_{\tau\tau} \f$
    CoolPropDbl d2alpha0_dTau2(void) {
        if (!_d2alpha0_dTau2) _d2alpha0_dTau2 = calc_d2alpha0_dTau2();
        return _d2alpha0_dTau2;
    };
    /// Return the term \f$ \alpha^0_{\tau\tau\tau} \f$
    CoolPropDbl d3alpha0_dTau3(void) {
        if (!_d3alpha0_dTau3) _d3alpha0_dTau3 = calc_d3alpha0_dTau3();
        return _d3alpha0_dTau3;
    };
    /// Return the term \f$ \alpha^0_{\delta\tau\tau} \f$
    CoolPropDbl d3alpha0_dDelta_dTau2(void) {
        if (!_d3alpha0_dDelta_dTau2) _d3alpha0_dDelta_dTau2 = calc_d3alpha0_dDelta_dTau2();
        return _d3alpha0_dDelta_dTau2;
    };
    /// Return the term \f$ \alpha^0_{\delta\delta\tau} \f$
    CoolPropDbl d3alpha0_dDelta2_dTau(void) {
        if (!_d3alpha0_dDelta2_dTau) _d3alpha0_dDelta2_dTau = calc_d3alpha0_dDelta2_dTau();
        return _d3alpha0_dDelta2_dTau;
    };
    /// Return the term \f$ \alpha^0_{\delta\delta\delta} \f$
    CoolPropDbl d3alpha0_dDelta3(void) {
        if (!_d3alpha0_dDelta3) _d3alpha0_dDelta3 = calc_d3alpha0_dDelta3();
        return _d3alpha0_dDelta3;
    };

    /// Return the term \f$ \alpha^r \f$
    CoolPropDbl alphar(void) {
        if (!_alphar) _alphar = calc_alphar();
        return _alphar;
    };
    /// Return the term \f$ \alpha^r_{\delta} \f$
    CoolPropDbl dalphar_dDelta(void) {
        if (!_dalphar_dDelta) _dalphar_dDelta = calc_dalphar_dDelta();
        return _dalphar_dDelta;
    };
    /// Return the term \f$ \alpha^r_{\tau} \f$
    CoolPropDbl dalphar_dTau(void) {
        if (!_dalphar_dTau) _dalphar_dTau = calc_dalphar_dTau();
        return _dalphar_dTau;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta} \f$
    CoolPropDbl d2alphar_dDelta2(void) {
        if (!_d2alphar_dDelta2) _d2alphar_dDelta2 = calc_d2alphar_dDelta2();
        return _d2alphar_dDelta2;
    };
    /// Return the term \f$ \alpha^r_{\delta\tau} \f$
    CoolPropDbl d2alphar_dDelta_dTau(void) {
        if (!_d2alphar_dDelta_dTau) _d2alphar_dDelta_dTau = calc_d2alphar_dDelta_dTau();
        return _d2alphar_dDelta_dTau;
    };
    /// Return the term \f$ \alpha^r_{\tau\tau} \f$
    CoolPropDbl d2alphar_dTau2(void) {
        if (!_d2alphar_dTau2) _d2alphar_dTau2 = calc_d2alphar_dTau2();
        return _d2alphar_dTau2;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta\delta} \f$
    CoolPropDbl d3alphar_dDelta3(void) {
        if (!_d3alphar_dDelta3) _d3alphar_dDelta3 = calc_d3alphar_dDelta3();
        return _d3alphar_dDelta3;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta\tau} \f$
    CoolPropDbl d3alphar_dDelta2_dTau(void) {
        if (!_d3alphar_dDelta2_dTau) _d3alphar_dDelta2_dTau = calc_d3alphar_dDelta2_dTau();
        return _d3alphar_dDelta2_dTau;
    };
    /// Return the term \f$ \alpha^r_{\delta\tau\tau} \f$
    CoolPropDbl d3alphar_dDelta_dTau2(void) {
        if (!_d3alphar_dDelta_dTau2) _d3alphar_dDelta_dTau2 = calc_d3alphar_dDelta_dTau2();
        return _d3alphar_dDelta_dTau2;
    };
    /// Return the term \f$ \alpha^r_{\tau\tau\tau} \f$
    CoolPropDbl d3alphar_dTau3(void) {
        if (!_d3alphar_dTau3) _d3alphar_dTau3 = calc_d3alphar_dTau3();
        return _d3alphar_dTau3;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta\delta\delta} \f$
    CoolPropDbl d4alphar_dDelta4(void) {
        if (!_d4alphar_dDelta4) _d4alphar_dDelta4 = calc_d4alphar_dDelta4();
        return _d4alphar_dDelta4;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta\delta\tau} \f$
    CoolPropDbl d4alphar_dDelta3_dTau(void) {
        if (!_d4alphar_dDelta3_dTau) _d4alphar_dDelta3_dTau = calc_d4alphar_dDelta3_dTau();
        return _d4alphar_dDelta3_dTau;
    };
    /// Return the term \f$ \alpha^r_{\delta\delta\tau\tau} \f$
    CoolPropDbl d4alphar_dDelta2_dTau2(void) {
        if (!_d4alphar_dDelta2_dTau2) _d4alphar_dDelta2_dTau2 = calc_d4alphar_dDelta2_dTau2();
        return _d4alphar_dDelta2_dTau2;
    };
    /// Return the term \f$ \alpha^r_{\delta\tau\tau\tau} \f$
    CoolPropDbl d4alphar_dDelta_dTau3(void) {
        if (!_d4alphar_dDelta_dTau3) _d4alphar_dDelta_dTau3 = calc_d4alphar_dDelta_dTau3();
        return _d4alphar_dDelta_dTau3;
    };
    /// Return the term \f$ \alpha^r_{\tau\tau\tau\tau} \f$
    CoolPropDbl d4alphar_dTau4(void) {
        if (!_d4alphar_dTau4) _d4alphar_dTau4 = calc_d4alphar_dTau4();
        return _d4alphar_dTau4;
    };
};

/** An abstract AbstractState generator class
 *
 *  This class should be derived and statically initialized in a C++ file.  In the initializer,
 *  the register_backend function should be called.  This will register the backend family, and
 *  when this generator is looked up in the map, the get_AbstractState function will be used
 *  to return an initialized instance
 */
class AbstractStateGenerator
{
   public:
    virtual AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) = 0;
    virtual ~AbstractStateGenerator(){};
};

/** Register a backend in the backend library (statically defined in AbstractState.cpp and not
 *  publicly accessible)
 */
void register_backend(const backend_families& bf, shared_ptr<AbstractStateGenerator> gen);

template <class T>
class GeneratorInitializer
{
   public:
    GeneratorInitializer(backend_families bf) {
        register_backend(bf, shared_ptr<AbstractStateGenerator>(new T()));
    };
};

} /* namespace CoolProp */
#endif /* ABSTRACTSTATE_H_ */
