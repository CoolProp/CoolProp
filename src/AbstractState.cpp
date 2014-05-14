/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#include "math.h"
#include "AbstractState.h"
#include "Backends\REFPROPBackend.h"
#include "Backends\HelmholtzEOSBackend.h"
#include "Backends\IncompressibleBackend.h"
#include "Fluids\FluidLibrary.h"

namespace CoolProp {

AbstractState * AbstractState::factory(const std::string &backend, const std::string &fluid_string)
{
    static std::string HEOS_string = "HEOS";
    if (!backend.compare("HEOS"))
    {
        if (fluid_string.find('&') == -1){
            return new HelmholtzEOSBackend(&get_fluid(fluid_string));
        }
        else{
            // Split at the '&'
            std::vector<std::string> components = strsplit(fluid_string,'&');

            return new HelmholtzEOSMixtureBackend(components);
        }
    }
    else if (!backend.compare("REFPROP"))
    {
        if (fluid_string.find('&') == -1){
            return new REFPROPBackend(fluid_string);
        }
        else{
            // Split at the '&'
            std::vector<std::string> components = strsplit(fluid_string,'&');

            return new REFPROPMixtureBackend(components);
        }
    }
    else if (!backend.compare("INCOMP"))
    {
        return new IncompressibleBackend(fluid_string);
    }
    else if (!backend.compare("BRINE"))
    {
        throw ValueError("BRINE backend not yet implemented");
    }
    else if (!backend.compare("TREND"))
    {
        throw ValueError("TREND backend not yet implemented");
    }
    else if (!backend.compare("?"))
    {
        std::size_t idel = fluid_string.find("::");
        // Backend has not been specified, and we have to figure out what the backend is by parsing the string
        if (idel == -1) // No '::' found, no backend specified, try HEOS, otherwise a failure
        {
            // Figure out what backend to use
            return factory(HEOS_string, fluid_string);
        }
        else
        {
            // Split string at the '::' into two std::string, call again
            return factory(std::string(fluid_string.begin(), fluid_string.begin() + idel), std::string(fluid_string.begin()+idel+2, fluid_string.end()));
        }

    }
    else
    {
        throw ValueError(format("Invalid backend name [%s] to factory function",backend.c_str()));
    }
}

bool AbstractState::clear() {
    // Reset all instances of CachedElement and overwrite
    // the internal double values with -_HUGE
    this->_fluid_type = FLUID_TYPE_UNDEFINED;
    this->_phase = iphase_unknown;
    this->_forceSinglePhase = false;
    this->_forceTwoPhase = false;
    this->_R = _HUGE;

    /// Ancillary curve values
    this->_rhoLanc.clear();
    this->_rhoVanc.clear();
    this->_pVanc.clear();
    this->_pLanc.clear();
    this->_TVanc.clear();
    this->_TLanc.clear();

    this->_critical.T = -_HUGE;
    this->_critical.hmolar = -_HUGE;
    this->_critical.p = -_HUGE;
    this->_critical.rhomolar = -_HUGE;
    this->_critical.smolar = -_HUGE;

    this->_reducing.T = -_HUGE;
    this->_reducing.hmolar = -_HUGE;
    this->_reducing.p = -_HUGE;
    this->_reducing.rhomolar = -_HUGE;
    this->_reducing.smolar = -_HUGE;

    /// Bulk values
    this->_rhomolar = -_HUGE;
    this->_T = -_HUGE;
    this->_p = -_HUGE;
    this->_Q = -_HUGE;
    this->_tau.clear();
    this->_delta.clear();

    this->_umolar.clear();
    this->_cpmolar.clear();
    this->_cvmolar.clear();
    this->_speed_sound.clear();
    this->_hmolar.clear();
    this->_smolar.clear();
    this->_logp.clear();
    this->_logrhomolar.clear();

    ///// Smoothing values
    //this->rhospline = -_HUGE;
    //this->dsplinedp = -_HUGE;
    //this->dsplinedh = -_HUGE;

    /// Cached low-level elements for in-place calculation of other properties
    this->_alpha0.clear();
    this->_dalpha0_dTau.clear();
    this->_dalpha0_dDelta.clear();
    this->_d2alpha0_dTau2.clear();
    this->_d2alpha0_dDelta_dTau.clear();
    this->_d2alpha0_dDelta2.clear();
    this->_d3alpha0_dTau3.clear();
    this->_d3alpha0_dDelta_dTau2.clear();
    this->_d3alpha0_dDelta2_dTau.clear();
    this->_d3alpha0_dDelta3.clear();
    this->_alphar.clear();
    this->_dalphar_dTau.clear();
    this->_dalphar_dDelta.clear();
    this->_d2alphar_dTau2.clear();
    this->_d2alphar_dDelta_dTau.clear();
    this->_d2alphar_dDelta2.clear();
    this->_d3alphar_dTau3.clear();
    this->_d3alphar_dDelta_dTau2.clear();
    this->_d3alphar_dDelta2_dTau.clear();
    this->_d3alphar_dDelta3.clear();

    this->_dalphar_dDelta_lim.clear();
    this->_d2alphar_dDelta2_lim.clear();
    this->_d2alphar_dDelta_dTau_lim.clear();
    this->_d3alphar_dDelta2_dTau_lim.clear();

    return true;
}

double AbstractState::keyed_output(int key)
{
    switch (key)
    {
    case iQ:
        return Q();
    case iT:
        return T();
    case iP:
        return p();
    case iDmolar:
        return rhomolar();
    case iDmass:
        return rhomolar()*molar_mass();
    case iHmolar:
        return hmolar();
    case iHmass:
        return hmolar()/molar_mass();
    case iSmolar:
        return smolar();
    case iSmass:
        return smolar()/molar_mass();
    case iUmolar:
        return umolar();
    case iUmass:
        return umolar()/molar_mass();
    case imolar_mass:
        return molar_mass();
    case iT_reducing:
        return get_reducing().T;
    case irhomolar_reducing:
        return get_reducing().rhomolar;
    //case iT_critical:
    //    return get_critical().T;
    //case irhomolar_critical:
    //    return get_critical().rhomolar; // TODO
    case ialpha0:
        return alpha0();
    case idalpha0_ddelta_consttau:
        return dalpha0_dDelta();
    case idalpha0_dtau_constdelta:
        return dalpha0_dTau();
    case iBvirial:
        return Bvirial();
    case idBvirial_dT:
        return dBvirial_dT();
    case iCvirial:
        return Cvirial();
    case idCvirial_dT:
        return dCvirial_dT();
    case iisothermal_compressibility:
        return isothermal_compressibility();
    default:
        throw ValueError(format("This input [%d: \"%s\"] is not valid for keyed_output",key,get_parameter_information(key,"short").c_str()));
    }
}

double AbstractState::tau(void){
    if (!_tau) _tau = calc_reciprocal_reduced_temperature();
    return _tau;
}
double AbstractState::delta(void){
    if (!_delta) _delta = calc_reduced_density();
    return _delta;
}
double AbstractState::Tmax(void){
    return calc_Tmax();
}
double AbstractState::Ttriple(void){
    return calc_Ttriple();
}
double AbstractState::pmax(void){
    return calc_pmax();
}

double AbstractState::hmolar(void){
    if (!_hmolar) _hmolar = calc_hmolar();
    return _hmolar;
}
double AbstractState::smolar(void){
    if (!_smolar) _smolar = calc_smolar();
    return _smolar;
}
double AbstractState::umolar(void){
    if (!_umolar) _umolar = calc_umolar();
    return _umolar;
}
double AbstractState::cpmolar(void){
    if (!_cpmolar) _cpmolar = calc_cpmolar();
    return _cpmolar;
}
double AbstractState::cvmolar(void){
    if (!_cvmolar) _cvmolar = calc_cvmolar();
    return _cvmolar;
}
double AbstractState::speed_sound(void){
    if (!_speed_sound) _speed_sound = calc_speed_sound();
    return _speed_sound;
}
double AbstractState::viscosity(void){
    if (!_viscosity) _viscosity = calc_viscosity();
    return _viscosity;
}
double AbstractState::conductivity(void){
    if (!_conductivity) _conductivity = calc_conductivity();
    return _conductivity;
}
double AbstractState::surface_tension(void){
    if (!_surface_tension) _surface_tension = calc_surface_tension();
    return _surface_tension;
}
double AbstractState::molar_mass(void){
    if (!_molar_mass) _molar_mass = calc_molar_mass();
    return _molar_mass;
}
double AbstractState::gas_constant(void){
    if (!_gas_constant) _gas_constant = calc_gas_constant();
    return _gas_constant;
}
double AbstractState::fugacity_coefficient(int i){
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity_coefficient(i);
}
double AbstractState::isothermal_compressibility(void){
	return 1.0/_rhomolar*first_partial_deriv(iDmolar, iP, iT);
}
double AbstractState::isobaric_expansion_coefficient(void){
	return -1.0/pow(_rhomolar,2)*first_partial_deriv(iDmolar, iT, iP);
}
double AbstractState::Bvirial(void){ return calc_Bvirial(); }
double AbstractState::Cvirial(void){ return calc_Cvirial(); }
double AbstractState::dBvirial_dT(void){ return calc_dBvirial_dT(); }
double AbstractState::dCvirial_dT(void){ return calc_dCvirial_dT(); }

//	// ----------------------------------------
//	// Smoothing functions for density
//	// ----------------------------------------
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodh_constp_smoothed(double xend);
//	/// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//	virtual double AbstractState::drhodp_consth_smoothed(double xend);
//	/// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
//	virtual void AbstractState::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);
//
//
//	// ----------------------------------------
//	// Transport properties // TODO: Fix it!
//	// ----------------------------------------
//
//	virtual double AbstractState::surface_tension(void);

} /* namespace CoolProp */
