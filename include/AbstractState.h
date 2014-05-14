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
    long _phase;
    bool _forceSinglePhase, _forceTwoPhase;

    //~ bool isCompressibleFluid(void){
        //~ return !(_fluid_type == FLUID_TYPE_INCOMPRESSIBLE_LIQUID
              //~ || _fluid_type == FLUID_TYPE_INCOMPRESSIBLE_SOLUTION);
    //~ }

    //~ bool checkCompressible(void){
        //~ if (!this->isCompressibleFluid()){throw ValueError(ERR_NOT_COMPRESSIBLE);}
        //~ return true;
    //~ }

    bool isHomogeneousPhase(void){
        return (this->_phase==iphase_liquid || this->_phase==iphase_gas || this->_phase == iphase_supercritical);
    }

    bool isTwoPhase(void){
        return (this->_phase==iphase_twophase);
    }

    //~ bool checkTwoPhase(void){
        //~ if (!this->isCompressibleFluid()){throw ValueError(ERR_NOT_A_TWO_PHASE_FLUID);}
        //~ if (!this->isTwoPhase()&&!_forceTwoPhase){throw ValueError(ERR_NOT_A_TWO_PHASE_STATE);}
        //~ return true;
    //~ }

    //~ bool checkSinglePhase(void){
        //~ if (!this->isHomogeneousPhase()||!_forceSinglePhase){throw ValueError(ERR_NOT_A_TWO_PHASE_FUNCTION);}
        //~ return true;
    //~ }

    
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

    CachedElement _hmolar, _smolar, _umolar, _logp, _logrhomolar, _cpmolar, _cvmolar, _speed_sound;

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
    virtual long double calc_hmolar(void){throw NotImplementedError("calc_hmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar entropy in J/mol/K
    virtual long double calc_smolar(void){throw NotImplementedError("calc_smolar is not implemented for this backend");};
    /// Using this backend, calculate the molar internal energy in J/mol
    virtual long double calc_umolar(void){throw NotImplementedError("calc_umolar is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-pressure specific heat in J/mol/K
    virtual long double calc_cpmolar(void){throw NotImplementedError("calc_cpmolar is not implemented for this backend");};
    /// Using this backend, calculate the molar constant-volume specific heat in J/mol/K
    virtual long double calc_cvmolar(void){throw NotImplementedError("calc_cvmolar is not implemented for this backend");};
    /// Using this backend, calculate the speed of sound in m/s
    virtual long double calc_speed_sound(void){throw NotImplementedError("calc_speed_sound is not implemented for this backend");};
    /// Using this backend, calculate the isothermal compressibility \f$ \kappa = -\frac{1}{v}\left.\frac{\partial v}{\partial p}\right|_T=\frac{1}{\rho}\left.\frac{\partial \rho}{\partial p}\right|_T\f$  in 1/Pa
    virtual long double calc_isothermal_compressibility(void){throw NotImplementedError("calc_isothermal_compressibility is not implemented for this backend");};
    /// Using this backend, calculate the isobaric expansion coefficient \f$ \beta = \frac{1}{v}\left.\frac{\partial v}{\partial T}\right|_p = -\frac{1}{\rho}\left.\frac{\partial \rho}{\partial T}\right|_p\f$  in 1/K
    virtual long double calc_isobaric_expansion_coefficient(void){throw NotImplementedError("calc_isobaric_expansion_coefficient is not implemented for this backend");};
    /// Using this backend, calculate the viscosity in Pa-s
    virtual long double calc_viscosity(void){throw NotImplementedError("calc_viscosity is not implemented for this backend");};
    /// Using this backend, calculate the thermal conductivity in W/m/K
    virtual long double calc_conductivity(void){throw NotImplementedError("calc_conductivity is not implemented for this backend");};
    /// Using this backend, calculate the surface tension in N/m
    virtual long double calc_surface_tension(void){throw NotImplementedError("calc_surface_tension is not implemented for this backend");};
    /// Using this backend, calculate the molar mass in kg/mol
    virtual long double calc_molar_mass(void){throw NotImplementedError("calc_molar_mass is not implemented for this backend");};
    /// Using this backend, calculate the pressure in Pa
    virtual long double calc_pressure(void){throw NotImplementedError("calc_pressure is not implemented for this backend");};
    /// Using this backend, calculate the universal gas constant \f$R_u\f$ in J/mol/K
    virtual long double calc_gas_constant(void){throw NotImplementedError("calc_gas_constant is not implemented for this backend");};
    /// Using this backend, calculate the fugacity coefficient (dimensionless)
    virtual long double calc_fugacity_coefficient(int i){throw NotImplementedError("calc_fugacity_coefficient is not implemented for this backend");};


    // Derivatives of residual helmholtz energy
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    virtual long double calc_alphar(void){throw NotImplementedError("calc_alphar is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    virtual long double calc_dalphar_dDelta(void){throw NotImplementedError("calc_dalphar_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    virtual long double calc_dalphar_dTau(void){throw NotImplementedError("calc_dalphar_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    virtual long double calc_d2alphar_dDelta2(void){throw NotImplementedError("calc_d2alphar_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    virtual long double calc_d2alphar_dDelta_dTau(void){throw NotImplementedError("calc_d2alphar_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    virtual long double calc_d2alphar_dTau2(void){throw NotImplementedError("calc_d2alphar_dTau2 is not implemented for this backend");};
    
    // Derivatives of ideal-gas helmholtz energy
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0\f$ (dimensionless)
    virtual long double calc_alpha0(void){throw NotImplementedError("calc_alpha0 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta}\f$ (dimensionless)
    virtual long double calc_dalpha0_dDelta(void){throw NotImplementedError("calc_dalpha0_dDelta is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau}\f$ (dimensionless)
    virtual long double calc_dalpha0_dTau(void){throw NotImplementedError("calc_dalpha0_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\delta}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dDelta_dTau(void){throw NotImplementedError("calc_d2alpha0_dDelta_dTau is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\delta\tau}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dDelta2(void){throw NotImplementedError("calc_d2alpha0_dDelta2 is not implemented for this backend");};
    /// Using this backend, calculate the ideal-gas Helmholtz energy term \f$\alpha^0_{\tau\tau}\f$ (dimensionless)
    virtual long double calc_d2alpha0_dTau2(void){throw NotImplementedError("calc_d2alpha0_dTau2 is not implemented for this backend");};

    virtual void calc_reducing_state(void){throw NotImplementedError("calc_reducing_state is not implemented for this backend");};

    /// Using this backend, calculate the maximum temperature in K
    virtual long double calc_Tmax(void){throw NotImplementedError("calc_Tmax is not implemented for this backend");};
    /// Using this backend, calculate the maximum pressure in Pa
    virtual long double calc_pmax(void){throw NotImplementedError("calc_pmax is not implemented for this backend");};

    /// Using this backend, calculate the 20-year global warming potential (GWP)
    virtual long double calc_GWP20(void){throw NotImplementedError("calc_GWP20 is not implemented for this backend");};
    /// Using this backend, calculate the 100-year global warming potential (GWP)
    virtual long double calc_GWP100(void){throw NotImplementedError("calc_GWP100 is not implemented for this backend");};
    /// Using this backend, calculate the 500-year global warming potential (GWP)
    virtual long double calc_GWP500(void){throw NotImplementedError("calc_GWP500 is not implemented for this backend");};
    /// Using this backend, calculate the ozone depletion potential (ODP)
    virtual long double calc_ODP(void){throw NotImplementedError("calc_ODP is not implemented for this backend");};
    /// Using this backend, calculate the flame hazard
    virtual long double calc_flame_hazard(void){throw NotImplementedError("calc_flame_hazard is not implemented for this backend");};
    /// Using this backend, calculate the health hazard
    virtual long double calc_health_hazard(void){throw NotImplementedError("calc_health_hazard is not implemented for this backend");};
    /// Using this backend, calculate the physical hazard
    virtual long double calc_physical_hazard(void){throw NotImplementedError("calc_physical_hazard is not implemented for this backend");};
    
    /// Calculate the first partial derivative for the desired derivative
    virtual long double calc_first_partial_deriv(int Of, int Wrt, int Constant){throw NotImplementedError("calc_first_partial_deriv is not implemented for this backend");};

    /// Using this backend, calculate the reduced density (rho/rhoc)
    virtual long double calc_reduced_density(void){throw NotImplementedError("calc_reduced_density is not implemented for this backend");};
    /// Using this backend, calculate the reciprocal reduced temperature (Tc/T)
    virtual long double calc_reciprocal_reduced_temperature(void){throw NotImplementedError("calc_reciprocal_reduced_temperature is not implemented for this backend");};

    /// Using this backend, calculate the second virial coefficient
    virtual long double calc_Bvirial(void){throw NotImplementedError("calc_Bvirial is not implemented for this backend");};
    /// Using this backend, calculate the third virial coefficient
    virtual long double calc_Cvirial(void){throw NotImplementedError("calc_Cvirial is not implemented for this backend");};
    /// Using this backend, calculate the derivative dB/dT
    virtual long double calc_dBvirial_dT(void){throw NotImplementedError("calc_dBvirial_dT is not implemented for this backend");};
    /// Using this backend, calculate the derivative dC/dT
    virtual long double calc_dCvirial_dT(void){throw NotImplementedError("calc_dCvirial_dT is not implemented for this backend");};

    /// Using this backend, get the name of the fluid
    virtual std::string calc_name(void){throw NotImplementedError("calc_name is not implemented for this backend");};
    
    /// Using this backend, get the triple point temperature in K
    virtual long double calc_Ttriple(void){throw NotImplementedError("calc_Ttriple is not implemented for this backend");};

public:

    virtual long double calc_melt_p_T(long double T){throw NotImplementedError("calc_melt_p_T is not implemented for this backend");};
    virtual long double calc_melt_T_p(long double p){throw NotImplementedError("calc_melt_T_p is not implemented for this backend");};
    virtual long double calc_melt_rho_T(long double T){throw NotImplementedError("calc_melt_rho_T is not implemented for this backend");};

    AbstractState(){};
    virtual ~AbstractState(){};

    
    
    /// A factory function to return a pointer to a new-allocated instance of one of the backends.
    /**
    Very Important!! : You must ensure to delete the backend instance that is created, otherwise there will be a memory leak

    The backend that is selected is based on the string passed in:
    
    1. If it starts with "REFPROP-", or no backend specification is provided, the function will assume that the backend is the CORE backend (for backwards-compatibility reasons)
    2. If it starts with "REFPROP-", the REFPROP backend will be used.  The remaining part of the string should then 
       either be
       1. A pure or pseudo-pure fluid name (eg. "PROPANE" or "R410A"), yielding a REFPROPBackend instance.
       2. A string that encodes the components of the mixture with a vertical bar between them (e.g. "R32|R125"), yielding a REFPROPMixtureBackend instance.
    3. If it starts with "TTSE", the TTSE backend will be used, yielding a TTSEBackend instance
    4. If it starts with "BICUBIC", the BICUBIC backend will be used, yielding a BICUBICBackend instance

    */
    static AbstractState * factory(const std::string &backend, const std::string &fluid_string);
    
    bool clear();
    virtual void update(long input_pair, double Value1, double Value2) = 0;
    virtual void set_mole_fractions(const std::vector<long double> &mole_fractions) = 0;
    virtual void set_mass_fractions(const std::vector<long double> &mass_fractions) = 0;

    void set_mole_fractions(const std::vector<double> &mole_fractions){set_mole_fractions(std::vector<long double>(mole_fractions.begin(), mole_fractions.end()));};
    void set_mass_fractions(const std::vector<double> &mass_fractions){set_mass_fractions(std::vector<long double>(mass_fractions.begin(), mass_fractions.end()));};
    
    // The derived classes must implement this function to define whether they use mole fractions (true) or mass fractions (false)
    virtual bool using_mole_fractions(void) = 0;

    const CoolProp::SimpleState & get_reducing(){return _reducing;};

    double keyed_output(int key);

    long double first_partial_deriv(int Of, int Wrt, int Constant){return calc_first_partial_deriv(Of,Wrt,Constant);};

    // Limits
    double Tmax(void);
    double pmax(void);
    double Ttriple(void);

    std::string name(){return calc_name();};

    // ----------------------------------------
    // Bulk properties - temperature and density are directly calculated every time
    // All other parameters are calculated on an as-needed basis
    // ----------------------------------------
    double T(void)  {return _T;};
    double rhomolar(void){return _rhomolar;};
    double p(void)  {return _p;};
    double Q(void)  {return _Q;};

    double tau(void);
    double delta(void);

    double molar_mass(void);
    double gas_constant(void);

    double Bvirial(void);
    double dBvirial_dT(void);
    double Cvirial(void);
    double dCvirial_dT(void);

    double hmolar(void);
    double smolar(void);
    double umolar(void);
    double cpmolar(void);
    double cvmolar(void);
    double speed_sound(void);
    double isothermal_compressibility(void);
    double isobaric_expansion_coefficient(void);
    double fugacity_coefficient(int i);
    //double fundamental_derivative_of_gas_dynamics(void);

    // ----------------------------------------
    // Transport properties
    // ----------------------------------------
    double viscosity(void);
    double conductivity(void);
    double surface_tension(void);

    // ----------------------------------------
    // Helmholtz energy and derivatives
    // ----------------------------------------
    /// Return the derivative \f$ \alpha^0 \f$
    long double alpha0(void){
        if (!_alpha0) _alpha0 = calc_alpha0();
        return _alpha0;
    };
    long double dalpha0_dDelta(void){
        if (!_dalpha0_dDelta) _dalpha0_dDelta = calc_dalpha0_dDelta();
        return _dalpha0_dDelta;
    };
    long double dalpha0_dTau(void){
        if (!_dalpha0_dTau) _dalpha0_dTau = calc_dalpha0_dTau();
        return _dalpha0_dTau;
    };
    long double d2alpha0_dDelta2(void){
        if (!_d2alpha0_dDelta2) _d2alpha0_dDelta2 = calc_d2alpha0_dDelta2();
        return _d2alpha0_dDelta2;
    };
    long double d2alpha0_dDelta_dTau(void){
        if (!_d2alpha0_dDelta_dTau) _d2alpha0_dDelta_dTau = calc_d2alpha0_dDelta_dTau();
        return _d2alpha0_dDelta_dTau;
    };
    long double d2alpha0_dTau2(void){
        if (!_d2alpha0_dTau2) _d2alpha0_dTau2 = calc_d2alpha0_dTau2();
        return _d2alpha0_dTau2;
    };
    /*
    virtual double d3alpha0_dDelta3(void) = 0;
    virtual double d3alpha0_dDelta2_dTau(void) = 0;
    virtual double d3alpha0_dDelta_dTau2(void) = 0;
    virtual double d3alpha0_dTau3(void) = 0;
    */

    long double alphar(void){
        if (!_alphar) _alphar = calc_alphar();
        return _alphar;
    };
    long double dalphar_dDelta(void){
        if (!_dalphar_dDelta) _dalphar_dDelta = calc_dalphar_dDelta();
        return _dalphar_dDelta;
    };
    long double dalphar_dTau(void){
        if (!_dalphar_dTau) _dalphar_dTau = calc_dalphar_dTau();
        return _dalphar_dTau;
    };
    long double d2alphar_dDelta2(void){
        if (!_d2alphar_dDelta2) _d2alphar_dDelta2 = calc_d2alphar_dDelta2();
        return _d2alphar_dDelta2;
    };
    long double d2alphar_dDelta_dTau(void){
        if (!_d2alphar_dDelta_dTau) _d2alphar_dDelta_dTau = calc_d2alphar_dDelta_dTau();
        return _d2alphar_dDelta_dTau;
    };
    long double d2alphar_dTau2(void){
        if (!_d2alphar_dTau2) _d2alphar_dTau2 = calc_d2alphar_dTau2();
        return _d2alphar_dTau2;
    };
    /*
    virtual double d3alphar_dDelta3(void) = 0;
    virtual double d3alphar_dDelta2_dTau(void) = 0;
    virtual double d3alphar_dDelta_dTau2(void) = 0;
    virtual double d3alphar_dTau3(void) = 0;

    virtual double dalphar_dDelta_lim(void) = 0;
    virtual double d2alphar_dDelta2_lim(void) = 0;
    virtual double d2alphar_dDelta_dTau_lim(void) = 0;
    virtual double d3alphar_dDelta2_dTau_lim(void) = 0;
    */
};

class AbstractStateWrapper
{
protected:
    AbstractState *p;
public:
    AbstractStateWrapper(){this->p = NULL;};
    AbstractStateWrapper(const std::string &backend, const std::string &fluid_string){
        this->p = AbstractState::factory(backend, fluid_string);
    };
    ~AbstractStateWrapper(){delete this->p;};
    void update(long input_pair, double Value1, double Value2){ this->p->update(input_pair,Value1,Value2); };
    double keyed_output(int key) { return this->p->keyed_output(key); }
};

} /* namespace CoolProp */
#endif /* ABSTRACTSTATE_H_ */
