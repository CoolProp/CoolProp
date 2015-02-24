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

#include <numeric>

namespace CoolProp {

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
class AbstractState {
protected:

    /// Some administrative variables
    long _fluid_type;
    phases _phase; ///< The key for the phase from CoolProp::phases enum

    bool isSupercriticalPhase(void){
        return (this->_phase == iphase_supercritical || this->_phase == iphase_supercritical_liquid || this->_phase == iphase_supercritical_gas);
    }
    
    bool isHomogeneousPhase(void){
        return (this->_phase==iphase_liquid || this->_phase==iphase_gas || isSupercriticalPhase() || this->_phase == iphase_critical_point);
    }

    bool isTwoPhase(void){
        return (this->_phase==iphase_twophase);
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

    CachedElement _hmolar, _smolar, _umolar, _logp, _logrhomolar, _cpmolar, _cp0molar, _cvmolar, _speed_sound, _gibbsmolar;

    /// Ancillary values
    CachedElement _rhoLanc, _rhoVanc, _pLanc, _pVanc, _TLanc, _TVanc;

    CachedElement _fugacity_coefficient;

    /// Smoothing values
    double _rhospline, _dsplinedp, _dsplinedh;

    /// Cached low-level elements for in-place calculation of other properties
    CachedElement _alpha0, _dalpha0_dTau, _dalpha0_dDelta, _d2alpha0_dTau2, _d2alpha0_dDelta_dTau,
            _d2alpha0_dDelta2, _d3alpha0_dTau3, _d3alpha0_dDelta_dTau2, _d3alpha0_dDelta2_dTau,
            _d3alpha0_dDelta3, _alphar, _dalphar_dTau, _dalphar_dDelta, _d2alphar_dTau2, _d2alphar_dDelta_dTau,
            _d2alphar_dDelta2, _d3alphar_dTau3, _d3alphar_dDelta_dTau2, _d3alphar_dDelta2_dTau,
            _d3alphar_dDelta3;

    CachedElement _dalphar_dDelta_lim, _d2alphar_dDelta2_lim,
            _d2alphar_dDelta_dTau_lim, _d3alphar_dDelta2_dTau_lim;

    /// Two-Phase variables
    CachedElement _rhoLmolar, _rhoVmolar;

    // ----------------------------------------
    // Property accessors to be optionally implemented by the backend
    // for properties that are not always calculated
    // ----------------------------------------
    /// Using this backend, calculate the molar enthalpy in J/mol
    virtual CoolPropDbl calc_hmolar(void){throw NotImplementedError("calc_hmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar entropy in J/mol/K
    virtual CoolPropDbl calc_smolar(void){throw NotImplementedError("calc_smolar is not implemented for this backend");};
    /// Using this backend, calculate the molar internal energy in J/mol
    virtual CoolPropDbl calc_umolar(void){throw NotImplementedError("calc_umolar is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-pressure specific heat in J/mol/K
    virtual CoolPropDbl calc_cpmolar(void){throw NotImplementedError("calc_cpmolar is not implemented for this backend");};
    /// Using this backend, calculate the ideal gas molar constant-pressure specific heat in J/mol/K
    virtual CoolPropDbl calc_cpmolar_idealgas(void){throw NotImplementedError("calc_cpmolar_idealgas is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-volume specific heat in J/mol/K
    virtual CoolPropDbl calc_cvmolar(void){throw NotImplementedError("calc_cvmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar Gibbs function in J/mol
    virtual CoolPropDbl calc_gibbsmolar(void){throw NotImplementedError("calc_gibbsmolar is not implemented for this backend");};
    /// Using this backend, calculate the speed of sound in m/s
    virtual CoolPropDbl calc_speed_sound(void){throw NotImplementedError("calc_speed_sound is not implemented for this backend");};
    /// Using this backend, calculate the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    virtual CoolPropDbl calc_isothermal_compressibility(void){throw NotImplementedError("calc_isothermal_compressibility is not implemented for this backend");};
    /// Using this backend, calculate the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    virtual CoolPropDbl calc_isobaric_expansion_coefficient(void){throw NotImplementedError("calc_isobaric_expansion_coefficient is not implemented for this backend");};
    /// Using this backend, calculate the viscosity in Pa-s
    virtual CoolPropDbl calc_viscosity(void){throw NotImplementedError("calc_viscosity is not implemented for this backend");};
    /// Using this backend, calculate the thermal conductivity in W/m/K
    virtual CoolPropDbl calc_conductivity(void){throw NotImplementedError("calc_conductivity is not implemented for this backend");};
    /// Using this backend, calculate the surface tension in N/m
    virtual CoolPropDbl calc_surface_tension(void){throw NotImplementedError("calc_surface_tension is not implemented for this backend");};
    /// Using this backend, calculate the molar mass in kg/mol
    virtual CoolPropDbl calc_molar_mass(void){throw NotImplementedError("calc_molar_mass is not implemented for this backend");};
    /// Using this backend, calculate the acentric factor
    virtual CoolPropDbl calc_acentric_factor(void){throw NotImplementedError("calc_acentric_factor is not implemented for this backend");};   
    /// Using this backend, calculate the pressure in Pa
    virtual CoolPropDbl calc_pressure(void){throw NotImplementedError("calc_pressure is not implemented for this backend");};
    /// Using this backend, calculate the universal gas constant \f$R_u\f$ in J/mol/K
    virtual CoolPropDbl calc_gas_constant(void){throw NotImplementedError("calc_gas_constant is not implemented for this backend");};
    /// Using this backend, calculate the fugacity coefficient (dimensionless)
    virtual CoolPropDbl calc_fugacity_coefficient(int i){throw NotImplementedError("calc_fugacity_coefficient is not implemented for this backend");};


    // Derivatives of residual helmholtz energy
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    virtual CoolPropDbl calc_alphar(void){throw NotImplementedError("calc_alphar is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalphar_dDelta(void){throw NotImplementedError("calc_dalphar_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalphar_dTau(void){throw NotImplementedError("calc_dalphar_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dDelta2(void){throw NotImplementedError("calc_d2alphar_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dDelta_dTau(void){throw NotImplementedError("calc_d2alphar_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alphar_dTau2(void){throw NotImplementedError("calc_d2alphar_dTau2 is not implemented for this backend");};
	/// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta3(void){throw NotImplementedError("calc_d3alpha0_dDelta3 is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta2_dTau(void){throw NotImplementedError("calc_d3alpha0_dDelta2_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dDelta_dTau2(void){throw NotImplementedError("calc_d3alpha0_dDelta_dTau2 is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alphar_dTau3(void){throw NotImplementedError("calc_d3alpha0_dTau3 is not implemented for this backend");};
	
    // Derivatives of ideal-gas helmholtz energy
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0\f$ (dimensionless)
    virtual CoolPropDbl calc_alpha0(void){throw NotImplementedError("calc_alpha0 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalpha0_dDelta(void){throw NotImplementedError("calc_dalpha0_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_dalpha0_dTau(void){throw NotImplementedError("calc_dalpha0_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dDelta_dTau(void){throw NotImplementedError("calc_d2alpha0_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dDelta2(void){throw NotImplementedError("calc_d2alpha0_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d2alpha0_dTau2(void){throw NotImplementedError("calc_d2alpha0_dTau2 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta\delta}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta3(void){throw NotImplementedError("calc_d3alpha0_dDelta3 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta2_dTau(void){throw NotImplementedError("calc_d3alpha0_dDelta2_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dDelta_dTau2(void){throw NotImplementedError("calc_d3alpha0_dDelta_dTau2 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau\tau}\f$ (dimensionless)
    virtual CoolPropDbl calc_d3alpha0_dTau3(void){throw NotImplementedError("calc_d3alpha0_dTau3 is not implemented for this backend");};

    virtual void calc_reducing_state(void){throw NotImplementedError("calc_reducing_state is not implemented for this backend");};

    /// Using this backend, calculate the maximum temperature in K
    virtual CoolPropDbl calc_Tmax(void){throw NotImplementedError("calc_Tmax is not implemented for this backend");};
    /// Using this backend, calculate the minimum temperature in K
    virtual CoolPropDbl calc_Tmin(void){throw NotImplementedError("calc_Tmin is not implemented for this backend");};
    /// Using this backend, calculate the maximum pressure in Pa
    virtual CoolPropDbl calc_pmax(void){throw NotImplementedError("calc_pmax is not implemented for this backend");};

    /// Using this backend, calculate the 20-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP20(void){throw NotImplementedError("calc_GWP20 is not implemented for this backend");};
    /// Using this backend, calculate the 100-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP100(void){throw NotImplementedError("calc_GWP100 is not implemented for this backend");};
    /// Using this backend, calculate the 500-year global warming potential (GWP)
    virtual CoolPropDbl calc_GWP500(void){throw NotImplementedError("calc_GWP500 is not implemented for this backend");};
    /// Using this backend, calculate the ozone depletion potential (ODP)
    virtual CoolPropDbl calc_ODP(void){throw NotImplementedError("calc_ODP is not implemented for this backend");};
    /// Using this backend, calculate the flame hazard
    virtual CoolPropDbl calc_flame_hazard(void){throw NotImplementedError("calc_flame_hazard is not implemented for this backend");};
    /// Using this backend, calculate the health hazard
    virtual CoolPropDbl calc_health_hazard(void){throw NotImplementedError("calc_health_hazard is not implemented for this backend");};
    /// Using this backend, calculate the physical hazard
    virtual CoolPropDbl calc_physical_hazard(void){throw NotImplementedError("calc_physical_hazard is not implemented for this backend");};

    /// Calculate the first partial derivative for the desired derivative
    virtual CoolPropDbl calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant);
    /// Calculate the second partial derivative using the given backend
    virtual CoolPropDbl calc_second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Of2, parameters Constant2);

    /// Using this backend, calculate the reduced density (rho/rhoc)
    virtual CoolPropDbl calc_reduced_density(void){throw NotImplementedError("calc_reduced_density is not implemented for this backend");};
    /// Using this backend, calculate the reciprocal reduced temperature (Tc/T)
    virtual CoolPropDbl calc_reciprocal_reduced_temperature(void){throw NotImplementedError("calc_reciprocal_reduced_temperature is not implemented for this backend");};

    /// Using this backend, calculate the second virial coefficient
    virtual CoolPropDbl calc_Bvirial(void){throw NotImplementedError("calc_Bvirial is not implemented for this backend");};
    /// Using this backend, calculate the third virial coefficient
    virtual CoolPropDbl calc_Cvirial(void){throw NotImplementedError("calc_Cvirial is not implemented for this backend");};
    /// Using this backend, calculate the derivative dB/dT
    virtual CoolPropDbl calc_dBvirial_dT(void){throw NotImplementedError("calc_dBvirial_dT is not implemented for this backend");};
    /// Using this backend, calculate the derivative dC/dT
    virtual CoolPropDbl calc_dCvirial_dT(void){throw NotImplementedError("calc_dCvirial_dT is not implemented for this backend");};
    /// Using this backend, calculate the compressibility factor Z \f$ Z = p/(\rho R T) \f$
    virtual CoolPropDbl calc_compressibility_factor(void){throw NotImplementedError("calc_compressibility_factor is not implemented for this backend");};

    /// Using this backend, get the name of the fluid
    virtual std::string calc_name(void){throw NotImplementedError("calc_name is not implemented for this backend");};

    /// Using this backend, get the triple point temperature in K
    virtual CoolPropDbl calc_Ttriple(void){throw NotImplementedError("calc_Ttriple is not implemented for this backend");};
    /// Using this backend, get the triple point pressure in Pa
    virtual CoolPropDbl calc_p_triple(void){throw NotImplementedError("calc_p_triple is not implemented for this backend");};

    /// Using this backend, get the critical point temperature in K
    virtual CoolPropDbl calc_T_critical(void){throw NotImplementedError("calc_T_critical is not implemented for this backend");};
	/// Using this backend, get the reducing point temperature in K
    virtual CoolPropDbl calc_T_reducing(void){throw NotImplementedError("calc_T_reducing is not implemented for this backend");};
    /// Using this backend, get the critical point pressure in Pa
    virtual CoolPropDbl calc_p_critical(void){throw NotImplementedError("calc_p_critical is not implemented for this backend");};
    /// Using this backend, get the reducing point pressure in Pa
    virtual CoolPropDbl calc_p_reducing(void){throw NotImplementedError("calc_p_reducing is not implemented for this backend");};
    /// Using this backend, get the critical point molar density in mol/m^3
    virtual CoolPropDbl calc_rhomolar_critical(void){throw NotImplementedError("calc_rhomolar_critical is not implemented for this backend");};
	/// Using this backend, get the reducing point molar density in mol/m^3
    virtual CoolPropDbl calc_rhomolar_reducing(void){throw NotImplementedError("calc_rhomolar_reducing is not implemented for this backend");};
    /// Using this backend, construct the phase envelope, the variable type describes the type of phase envelope to be built.
    virtual void calc_phase_envelope(const std::string &type){throw NotImplementedError("calc_phase_envelope is not implemented for this backend");};
    /// 
    virtual CoolPropDbl calc_rhomass(void){return _rhomolar*molar_mass();}
    virtual CoolPropDbl calc_hmass(void){return hmolar()/molar_mass();}
    virtual CoolPropDbl calc_smass(void){return smolar()/molar_mass();}
    virtual CoolPropDbl calc_cpmass(void){return cpmolar()/molar_mass();}
    virtual CoolPropDbl calc_cp0mass(void){return cp0molar()/molar_mass();}
    virtual CoolPropDbl calc_cvmass(void){return cvmolar()/molar_mass();}
    virtual CoolPropDbl calc_umass(void){return umolar()/molar_mass();}
    
    /// Update the states after having changed the reference state for enthalpy and entropy
    virtual void update_states(void){throw NotImplementedError("This backend does not implement update_states function");};
    
    virtual CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value){throw NotImplementedError("This backend does not implement calc_melting_line function");};
    
    /// @param param The key for the parameter to be returned
    /// @param Q The quality for the parameter that is given (0 = saturated liquid, 1 = saturated vapor)
    /// @param given The key for the parameter that is given
    /// @param value The value for the parameter that is given
    virtual CoolPropDbl calc_saturation_ancillary(parameters param, int Q, parameters given, double value){throw NotImplementedError("This backend does not implement calc_saturation_ancillary");};
    
    /// Using this backend, calculate the phase
    virtual phases calc_phase(void){throw NotImplementedError("This backend does not implement calc_phase function");};
    /// Using this backend, specify the phase to be used for all further calculations
    virtual void calc_specify_phase(phases phase){throw NotImplementedError("This backend does not implement calc_specify_phase function");};
    /// Using this backend, unspecify the phase
    virtual void calc_unspecify_phase(void){throw NotImplementedError("This backend does not implement calc_unspecify_phase function");};
    /// Using this backend get a vector of fluid names
	virtual std::vector<std::string> calc_fluid_names(void){throw NotImplementedError("This backend does not implement calc_fluid_names function");};
    /// Using this backend, calculate a phase given by the state string
    /// @param state A string that describes the state desired, one of "hs_anchor", "critical"/"crit", "reducing"
    virtual const CoolProp::SimpleState & calc_state(const std::string &state){throw NotImplementedError("calc_state is not implemented for this backend");};
    
    virtual const CoolProp::PhaseEnvelopeData & calc_phase_envelope_data(void){throw NotImplementedError("calc_phase_envelope_data is not implemented for this backend");};
    
    virtual std::vector<CoolPropDbl> calc_mole_fractions_liquid(void){throw NotImplementedError("calc_mole_fractions_liquid is not implemented for this backend");};
    virtual std::vector<CoolPropDbl> calc_mole_fractions_vapor(void){throw NotImplementedError("calc_mole_fractions_vapor is not implemented for this backend");};

    /// Get the minimum fraction (mole, mass, volume) for incompressible fluid
    virtual CoolPropDbl calc_fraction_min(void){throw NotImplementedError("calc_fraction_min is not implemented for this backend");};
    /// Get the maximum fraction (mole, mass, volume) for incompressible fluid
    virtual CoolPropDbl calc_fraction_max(void){throw NotImplementedError("calc_fraction_max is not implemented for this backend");};
    virtual CoolPropDbl calc_T_freeze(void){throw NotImplementedError("calc_T_freeze is not implemented for this backend");};
	
	virtual CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1){throw NotImplementedError("calc_first_saturation_deriv is not implemented for this backend");};
	virtual CoolPropDbl calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Of2, parameters Wrt2){throw NotImplementedError("calc_second_saturation_deriv is not implemented for this backend");};
    virtual CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant){throw NotImplementedError("calc_first_two_phase_deriv is not implemented for this backend");};
    virtual CoolPropDbl calc_second_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant, parameters Wrt2, parameters Constant2){throw NotImplementedError("calc_second_two_phase_deriv is not implemented for this backend");};
    virtual CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end){throw NotImplementedError("calc_first_two_phase_deriv_splined is not implemented for this backend");};
    
    virtual CoolPropDbl calc_saturated_liquid_keyed_output(parameters key){throw NotImplementedError("calc_saturated_liquid_keyed_output is not implemented for this backend");};
    virtual CoolPropDbl calc_saturated_vapor_keyed_output(parameters key){throw NotImplementedError("calc_saturated_vapor_keyed_output is not implemented for this backend");};

public:

    AbstractState(){clear();};
    virtual ~AbstractState(){};

    /// A factory function to return a pointer to a new-allocated instance of one of the backends.
    /**
     * @brief This is a convenience function to allow for the use of '&' delimited fluid names.  Slightly less computationally efficient than the 
     * @param backend The backend in use, one of "HEOS", "REFPROP", etc.
     * @param fluid_names Fluid names as a '&' delimited string
     * @return 
     */
    static AbstractState * factory(const std::string &backend, const std::string &fluid_names)
    {
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
    static AbstractState * factory(const std::string &backend, const std::vector<std::string> &fluid_names);

    // The derived classes must implement this function to define whether they use mole fractions (true) or mass fractions (false)
    virtual bool using_mole_fractions(void) = 0;
    virtual bool using_mass_fractions(void) = 0;
    virtual bool using_volu_fractions(void) = 0;

    virtual void update(CoolProp::input_pairs input_pair, double Value1, double Value2) = 0;
    virtual void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions) = 0;
    virtual void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions) = 0;
    virtual void set_volu_fractions(const std::vector<CoolPropDbl> &mass_fractions){throw NotImplementedError("Volume composition has not been implemented.");}
    /// Get the mole fractions of the components
    virtual const std::vector<CoolPropDbl> & get_mole_fractions(void) = 0;
    
	std::vector<std::string> fluid_names(void);
    
    /// Clear all the cached values
    bool clear();

    void set_mole_fractions(const std::vector<double> &mole_fractions){set_mole_fractions(std::vector<CoolPropDbl>(mole_fractions.begin(), mole_fractions.end()));};
    void set_mass_fractions(const std::vector<double> &mass_fractions){set_mass_fractions(std::vector<CoolPropDbl>(mass_fractions.begin(), mass_fractions.end()));};
    void set_volu_fractions(const std::vector<double> &volu_fractions){set_volu_fractions(std::vector<CoolPropDbl>(volu_fractions.begin(), volu_fractions.end()));};

    virtual const CoolProp::SimpleState & get_reducing_state(){return _reducing;};
    const CoolProp::SimpleState & get_state(const std::string &state){return calc_state(state);};

    // Limits
    double Tmin(void);
    double Tmax(void);
    double pmax(void);
    double Ttriple(void);
    
    /// Get the phase of the state
    phases phase(void){return calc_phase();};
    /// Specify the phase for all further calculations with this state class
    void specify_phase(phases phase){calc_specify_phase(phase);};
    /// Unspecify the phase and go back to calculating it based on the inputs
    void unspecify_phase(void){calc_unspecify_phase();};

    /// Return the critical temperature in K
    /**
    For pure fluids, this is the critical point temperature
    For mixtures, it is the exact critical point temperature calculated by the methods of Michelsen( \todo fill in reference)
    */
    double T_critical(void);
	/// Return the reducing point temperature in K
    double T_reducing(void);
    /// Return the critical pressure in Pa
    /**
    For pure fluids, this is the critical point pressure as defined by the author of the EOS
    For mixtures, it is the exact critical point pressure calculated by the methods of Michelsen( \todo fill in reference)
    */
    double p_critical(void);
    /// 
    /** \brief Return the critical molar density in mol/m^3
	 * 
	 * For pure fluids, this is the critical point molar density
     * For mixtures, it is the exact critical point molar density calculated by the methods of Michelsen( \todo fill in reference)
     */
    double rhomolar_critical(void);
	/** \brief Return the critical molar density in kg/m^3
	 * 
	 * For pure fluids, this is the critical point molar density
     * For mixtures, it is the exact critical point molar density calculated by the methods of Michelsen( \todo fill in reference)
     */
    double rhomass_critical(void);
	/** \brief Return the molar density at the reducing point in mol/m^3
	 * 
	 * For pure fluids, this is the critical point molar density
     * For mixtures, it is the exact critical point molar density calculated by the methods of Michelsen( \todo fill in reference)
     */
    double rhomolar_reducing(void);
	/** \brief Return the mass density at the reducing point in kg/m^3
	 * 
	 * For pure fluids, this is the critical point molar density
     * For mixtures, it is the exact critical point molar density calculated by the methods of Michelsen( \todo fill in reference)
     */
    double rhomass_reducing(void);
    
    /// Return the triple point pressure 
    double p_triple(void);

    std::string name(){return calc_name();};

    // ----------------------------------------
    // Bulk properties - temperature and density are directly calculated every time
    // All other parameters are calculated on an as-needed basis
    // ----------------------------------------
    /// Retrieve a value by key
    double keyed_output(int key);
    /// A trivial keyed output like molar mass that does not depend on the state
    double trivial_keyed_output(int key);
    /// Get an output from the saturated liquid state by key
    double saturated_liquid_keyed_output(parameters key){return calc_saturated_liquid_keyed_output(key);};
    /// Get an output from the saturated vapor state by key
    double saturated_vapor_keyed_output(parameters key){return calc_saturated_vapor_keyed_output(key);};
    
    /// Return the temperature in K
    double T(void)  {return _T;};
    /// Return the molar density in mol/m^3
    double rhomolar(void){return _rhomolar;};
    /// Return the mass density in kg/m^3
    double rhomass(void){return calc_rhomass();};
    /// Return the pressure in Pa
    double p(void)  {return _p;};
    /// Return the vapor quality (mol/mol); Q = 0 for saturated liquid
    double Q(void)  {return _Q;};
    /// Return the reciprocal of the reduced temperature (\f$\tau = T_c/T\f$)
    double tau(void);
    /// Return the reduced density (\f$\delta = \rho/\rho_c\f$)
    double delta(void);
    /// Return the molar mass in kg/mol
    double molar_mass(void);
    double acentric_factor(void);
    /// Return the mole-fraction weighted gas constant in J/mol/K
    
    double gas_constant(void);
    double Bvirial(void);
    double dBvirial_dT(void);
    double Cvirial(void);
    double dCvirial_dT(void);
    /// Return the compressibility factor \f$ Z = p/(rho R T) \f$
    double compressibility_factor(void);
    /// Return the molar enthalpy in J/mol
    double hmolar(void);
    /// Return the mass enthalpy in J/kg
    double hmass(void){return calc_hmass();};
    /// Return the molar entropy in J/mol/K
    double smolar(void);
    /// Return the molar entropy in J/kg/K
    double smass(void){return calc_smass();};
    /// Return the molar internal energy in J/mol
    double umolar(void);
    /// Return the mass internal energy in J/kg
    double umass(void){return calc_umass();};
    /// Return the molar constant pressure specific heat in J/mol/K
    double cpmolar(void);
    /// Return the mass constant pressure specific heat in J/kg/K
    double cpmass(void){return calc_cpmass();};
    /// Return the molar constant pressure specific heat for ideal gas part only in J/mol/K
    double cp0molar(void);
    /// Return the mass constant pressure specific heat for ideal gas part only in J/kg/K
    double cp0mass(void){return calc_cp0mass();};
    /// Return the molar constant volume specific heat in J/mol/K
    double cvmolar(void);
    /// Return the mass constant volume specific heat in J/kg/K
    double cvmass(void){return calc_cvmass();};
    /// Return the Gibbs function in J/mol
    double gibbsmolar(void){return calc_gibbsmolar();};
    /// Return the speed of sound in m/s
    double speed_sound(void);
    /// Return the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    double isothermal_compressibility(void);
    /// Return the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    double isobaric_expansion_coefficient(void);
    double fugacity_coefficient(int i);
    //double fundamental_derivative_of_gas_dynamics(void);
    
    std::vector<CoolPropDbl> mole_fractions_liquid(void){return calc_mole_fractions_liquid();};
    std::vector<CoolPropDbl> mole_fractions_vapor(void){return calc_mole_fractions_vapor();};
    
    // ----------------------------------------
    //    Partial derivatives
    // ----------------------------------------
    
    /** \brief The first partial derivative in homogeneous phases
     * 
     * \f[ \left(\frac{\partial A}{\partial B}\right)_C = \frac{\left(\frac{\partial A}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial A}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta}{\left(\frac{\partial B}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial B}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta} = \frac{N}{D}\f]
     */
    CoolPropDbl first_partial_deriv(parameters Of, parameters Wrt, parameters Constant){return calc_first_partial_deriv(Of, Wrt, Constant);};
    
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
    CoolPropDbl second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Of2, parameters Constant2){return calc_second_partial_deriv(Of1,Wrt1,Constant1,Of2,Constant2);};
    
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
	CoolPropDbl first_saturation_deriv(parameters Of1, parameters Wrt1){return calc_first_saturation_deriv(Of1,Wrt1);};
	
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
	 * @param Of2 The parameter that the second derivative is taken of
	 * @param Wrt2 The parameter that the second derivative is taken with respect to
	 * */
	CoolPropDbl second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Of2, parameters Wrt2){return calc_second_saturation_deriv(Of1,Wrt1,Of2,Wrt2);};
    
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
     * @param type_flag A flag describing how the derivative should be calculated, either normal or using splines
     * @return 
     */
    double first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant){
        return calc_first_two_phase_deriv(Of, Wrt, Constant);
    };
    
    double second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2){
        return calc_second_two_phase_deriv(Of, Wrt1, Constant1, Wrt2, Constant2);
    };
    
    double first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, double x_end){
        return calc_first_two_phase_deriv_splined(Of, Wrt, Constant, x_end);
    };
    
    // ----------------------------------------
    //    Phase envelope for mixtures
    // ----------------------------------------
    
    void build_phase_envelope(const std::string &type);
    const CoolProp::PhaseEnvelopeData &get_phase_envelope_data(){return calc_phase_envelope_data();};
    
    // ----------------------------------------
    //    Ancillary equations
    // ----------------------------------------
    
    /// Return true if the fluid has a melting line - default is false, but can be re-implemented by derived class
    virtual bool has_melting_line(void){return false;};
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
    /// Return the thermal conductivity in W/m/K
    double conductivity(void);
    /// Return the surface tension in N/m
    double surface_tension(void);
    /// Return the Prandtl number (dimensionless)
    double Prandtl(void){return cpmass()*viscosity()/conductivity();};

    // ----------------------------------------
    // Helmholtz energy and derivatives
    // ----------------------------------------
    /// Return the term \f$ \alpha^0 \f$
    CoolPropDbl alpha0(void){
        if (!_alpha0) _alpha0 = calc_alpha0();
        return _alpha0;
    };
    CoolPropDbl dalpha0_dDelta(void){
        if (!_dalpha0_dDelta) _dalpha0_dDelta = calc_dalpha0_dDelta();
        return _dalpha0_dDelta;
    };
    CoolPropDbl dalpha0_dTau(void){
        if (!_dalpha0_dTau) _dalpha0_dTau = calc_dalpha0_dTau();
        return _dalpha0_dTau;
    };
    CoolPropDbl d2alpha0_dDelta2(void){
        if (!_d2alpha0_dDelta2) _d2alpha0_dDelta2 = calc_d2alpha0_dDelta2();
        return _d2alpha0_dDelta2;
    };
    CoolPropDbl d2alpha0_dDelta_dTau(void){
        if (!_d2alpha0_dDelta_dTau) _d2alpha0_dDelta_dTau = calc_d2alpha0_dDelta_dTau();
        return _d2alpha0_dDelta_dTau;
    };
    CoolPropDbl d2alpha0_dTau2(void){
        if (!_d2alpha0_dTau2) _d2alpha0_dTau2 = calc_d2alpha0_dTau2();
        return _d2alpha0_dTau2;
    };
    CoolPropDbl d3alpha0_dTau3(void){
        if (!_d3alpha0_dTau3) _d3alpha0_dTau3 = calc_d3alpha0_dTau3();
        return _d3alpha0_dTau3;
    };
    CoolPropDbl d3alpha0_dDelta_dTau2(void){
        if (!_d3alpha0_dDelta_dTau2) _d3alpha0_dDelta_dTau2 = calc_d3alpha0_dDelta_dTau2();
        return _d3alpha0_dDelta_dTau2;
    };
    CoolPropDbl d3alpha0_dDelta2_dTau(void){
        if (!_d3alpha0_dDelta2_dTau) _d3alpha0_dDelta2_dTau = calc_d3alpha0_dDelta2_dTau();
        return _d3alpha0_dDelta2_dTau;
    };
    CoolPropDbl d3alpha0_dDelta3(void){
        if (!_d3alpha0_dDelta3) _d3alpha0_dDelta3 = calc_d3alpha0_dDelta3();
        return _d3alpha0_dDelta3;
    };

    CoolPropDbl alphar(void){
        if (!_alphar) _alphar = calc_alphar();
        return _alphar;
    };
    CoolPropDbl dalphar_dDelta(void){
        if (!_dalphar_dDelta) _dalphar_dDelta = calc_dalphar_dDelta();
        return _dalphar_dDelta;
    };
    CoolPropDbl dalphar_dTau(void){
        if (!_dalphar_dTau) _dalphar_dTau = calc_dalphar_dTau();
        return _dalphar_dTau;
    };
    CoolPropDbl d2alphar_dDelta2(void){
        if (!_d2alphar_dDelta2) _d2alphar_dDelta2 = calc_d2alphar_dDelta2();
        return _d2alphar_dDelta2;
    };
    CoolPropDbl d2alphar_dDelta_dTau(void){
        if (!_d2alphar_dDelta_dTau) _d2alphar_dDelta_dTau = calc_d2alphar_dDelta_dTau();
        return _d2alphar_dDelta_dTau;
    };
    CoolPropDbl d2alphar_dTau2(void){
        if (!_d2alphar_dTau2) _d2alphar_dTau2 = calc_d2alphar_dTau2();
        return _d2alphar_dTau2;
    };
	CoolPropDbl d3alphar_dDelta3(void){
        if (!_d3alphar_dDelta3) _d3alphar_dDelta3 = calc_d3alphar_dDelta3();
        return _d3alphar_dDelta3;
    };
	CoolPropDbl d3alphar_dDelta2_dTau(void){
        if (!_d3alphar_dDelta2_dTau) _d3alphar_dDelta2_dTau = calc_d3alphar_dDelta2_dTau();
        return _d3alphar_dDelta2_dTau;
    };
	CoolPropDbl d3alphar_dDelta_dTau2(void){
        if (!_d3alphar_dDelta_dTau2) _d3alphar_dDelta_dTau2 = d3alphar_dDelta_dTau2();
        return _d3alphar_dDelta_dTau2;
    };
	CoolPropDbl d3alphar_dTau3(void){
        if (!_d3alphar_dTau3) _d3alphar_dTau3 = calc_d3alphar_dTau3();
        return _d3alphar_dTau3;
    };
	
    /*
    virtual double dalphar_dDelta_lim(void) = 0;
    virtual double d2alphar_dDelta2_lim(void) = 0;
    virtual double d2alphar_dDelta_dTau_lim(void) = 0;
    virtual double d3alphar_dDelta2_dTau_lim(void) = 0;
    */
};

} /* namespace CoolProp */
#endif /* ABSTRACTSTATE_H_ */
