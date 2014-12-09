/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#include <stdlib.h>
#include "math.h"
#include "AbstractState.h"
#include "Backends/REFPROP/REFPROPBackend.h"
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Incompressible/IncompressibleBackend.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "Backends/Tabular/TabularBackends.h"

namespace CoolProp {

AbstractState * AbstractState::factory(const std::string &backend, const std::string &fluid_string)
{
    static std::string HEOS_string = "HEOS";
    if (!backend.compare("HEOS"))
    {
        if (fluid_string.find('&') == std::string::npos){
            return new HelmholtzEOSBackend(fluid_string);
        }
        else{
            // Split at the '&'
            std::vector<std::string> components = strsplit(fluid_string,'&');

            return new HelmholtzEOSMixtureBackend(components);
        }
    }
    else if (!backend.compare("REFPROP"))
    {
        if (fluid_string.find('&') == std::string::npos){
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
    else if (backend.find("TTSE&") == 0)
    {
        // Will throw if there is a problem with this backend
        shared_ptr<AbstractState> AS(factory(backend.substr(5), fluid_string));
        return new TTSEBackend(*AS.get());
    }
    else if (!backend.compare("TREND"))
    {
        throw ValueError("TREND backend not yet implemented");
    }
    else if (!backend.compare("?"))
    {
        std::size_t idel = fluid_string.find("::");
        // Backend has not been specified, and we have to figure out what the backend is by parsing the string
        if (idel == std::string::npos) // No '::' found, no backend specified, try HEOS, otherwise a failure
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
    this->_cp0molar.clear();
    this->_cvmolar.clear();
    this->_speed_sound.clear();
    this->_hmolar.clear();
    this->_smolar.clear();
    this->_gibbsmolar.clear();
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
double AbstractState::trivial_keyed_output(int key)
{
    if (get_debug_level()>=50) std::cout << format("AbstractState: keyed_output called for %s ",get_parameter_information(key,"short").c_str()) << std::endl;
    switch (key)
    {
    case imolar_mass:
        return molar_mass();
    case iT_min:
        return Tmin();
    case iT_triple:
        return Ttriple();
    case iT_max:
        return Tmax();
    case iP_max:
        return pmax();
    case iP_min:
    case iP_triple:
        return this->p_triple();
    case iT_reducing:
        return get_reducing_state().T;
    case irhomolar_reducing:
        return get_reducing_state().rhomolar;
	case iP_reducing:
		return get_reducing_state().p;
    case iP_critical:
        return this->p_critical();
    case iT_critical:
        return this->T_critical();
    case irhomolar_critical:
        return this->rhomolar_critical();
    case irhomass_critical:
        return this->rhomass_critical();
    case iODP:
        return this->calc_ODP();
    case iGWP100:
        return this->calc_GWP100();
    case iGWP20:
        return this->calc_GWP20();
    case iGWP500:
        return this->calc_GWP500();
    case ifraction_min:
        return this->calc_fraction_min();
    case ifraction_max:
        return this->calc_fraction_max();
    case iT_freeze:
        return this->calc_T_freeze();
    default:
        throw ValueError(format("This input [%d: \"%s\"] is not valid for trivial_keyed_output",key,get_parameter_information(key,"short").c_str()));
    }
}
double AbstractState::keyed_output(int key)
{
    if (get_debug_level()>=50) std::cout << format("AbstractState: keyed_output called for %s ",get_parameter_information(key,"short").c_str()) << std::endl;
    // Handle trivial inputs
    if (is_trivial_parameter(key))
    {
        return trivial_keyed_output(key);
    }
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
        return rhomass();
    case iHmolar:
        return hmolar();
    case iHmass:
        return hmass();
    case iSmolar:
        return smolar();
    case iSmass:
        return smass();
    case iUmolar:
        return umolar();
    case iUmass:
        return umass();
    case iCvmolar:
        return cvmolar();
    case iCvmass:
        return cvmass();
    case iCpmolar:
        return cpmolar();
    case iCp0molar:
        return cp0molar();
    case iCpmass:
        return cpmass();
    case iCp0mass:
        return cp0mass();
    case imolar_mass:
        return molar_mass();
    case iT_reducing:
        return get_reducing_state().T;
    case irhomolar_reducing:
        return get_reducing_state().rhomolar;
    case ispeed_sound:
        return speed_sound();
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
    case iviscosity:
        return viscosity();
    case iconductivity:
        return conductivity();
	case isurface_tension:
		return surface_tension();
    case iPhase:
        return phase();
    case iZ:
        return compressibility_factor();
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
double AbstractState::Tmin(void){
    return calc_Tmin();
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
double AbstractState::T_critical(void){
    return calc_T_critical();
}
double AbstractState::T_reducing(void){
    return calc_T_reducing();
}
double AbstractState::p_critical(void){
    return calc_p_critical();
}
double AbstractState::p_triple(void){
    return calc_p_triple();
}
double AbstractState::rhomolar_critical(void){
    return calc_rhomolar_critical();
}
double AbstractState::rhomass_critical(void){
    return calc_rhomolar_critical()*molar_mass();
}
double AbstractState::rhomolar_reducing(void){
    return calc_rhomolar_reducing();
}
double AbstractState::rhomass_reducing(void){
    return calc_rhomolar_reducing()*molar_mass();
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
double AbstractState::cp0molar(void){
    return calc_cpmolar_idealgas();
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
double AbstractState::melting_line(int param, int given, double value){
    return calc_melting_line(param, given, value);
}
double AbstractState::saturation_ancillary(parameters param, int Q, parameters given, double value){
    return calc_saturation_ancillary(param, Q, given, value);
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
void AbstractState::build_phase_envelope(const std::string &type)
{
    calc_phase_envelope(type);
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
double AbstractState::compressibility_factor(void){ return calc_compressibility_factor(); }

// Get the derivatives of the parameters in the partial derivative with respect to T and rho
void get_dT_drho(AbstractState &AS, parameters index, long double &dT, long double &drho)
{
    long double T = AS.T(),
                rho = AS.rhomolar(),
                rhor = AS.rhomolar_reducing(),
                Tr = AS.T_reducing(),
                dT_dtau = -pow(T, 2)/Tr,
                R = AS.gas_constant(),
                delta = rho/rhor,
                tau = Tr/T;
    
    switch (index)
    {
    case iT:
        dT = 1; drho = 0; break;
    case iDmolar:
        dT = 0; drho = 1; break;
    case iDmass:
        dT = 0; drho = AS.molar_mass(); break;
    case iP:
    {
        // dp/drho|T
        drho = R*T*(1+2*delta*AS.dalphar_dDelta()+pow(delta, 2)*AS.d2alphar_dDelta2());
        // dp/dT|rho
        dT = rho*R*(1+delta*AS.dalphar_dDelta() - tau*delta*AS.d2alphar_dDelta_dTau());
        break;
    }
    case iHmass:
    case iHmolar:
    {
        // dh/dT|rho
        dT = R*(-pow(tau,2)*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()) + (1+delta*AS.dalphar_dDelta()-tau*delta*AS.d2alphar_dDelta_dTau()));
        // dh/drhomolar|T
        drho = T*R/rho*(tau*delta*AS.d2alphar_dDelta_dTau()+delta*AS.dalphar_dDelta()+pow(delta,2)*AS.d2alphar_dDelta2());
        if (index == iHmass){
            // dhmolar/drhomolar|T * dhmass/dhmolar where dhmass/dhmolar = 1/mole mass
            drho /= AS.molar_mass();
            dT /= AS.molar_mass();
        }
        break;
    }
    case iSmass:
    case iSmolar:
    {
        // ds/dT|rho
        dT = R/T*(-pow(tau,2)*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()));
        // ds/drho|T
        drho = R/rho*(-(1+delta*AS.dalphar_dDelta()-tau*delta*AS.d2alphar_dDelta_dTau()));
        if (index == iSmass){
            // ds/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
            drho /= AS.molar_mass();
            dT /= AS.molar_mass();
        }
        break;
    }
    case iUmass:
    case iUmolar:
    {
        // du/dT|rho
        dT = R*(-pow(tau,2)*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()));
        // du/drho|T
        drho = AS.T()*R/rho*(tau*delta*AS.d2alphar_dDelta_dTau());
        if (index == iUmass){
            // du/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
            drho /= AS.molar_mass();
            dT /= AS.molar_mass();
        }
        break;
    }
    case iTau:
        dT = 1/dT_dtau; drho = 0; break;
    case iDelta:
        dT = 0; drho = 1/rhor; break;
    default:
        throw ValueError(format("input to get_dT_drho[%s] is invalid",get_parameter_information(index,"short").c_str()));
    }
}
long double AbstractState::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant)
{
	long double dOf_dT, dOf_drho, dWrt_dT, dWrt_drho, dConstant_dT, dConstant_drho;

    get_dT_drho(*this, Of, dOf_dT, dOf_drho);
    get_dT_drho(*this, Wrt, dWrt_dT, dWrt_drho);
    get_dT_drho(*this, Constant, dConstant_dT, dConstant_drho);

    return (dOf_dT*dConstant_drho-dOf_drho*dConstant_dT)/(dWrt_dT*dConstant_drho-dWrt_drho*dConstant_dT);
}
//    // ----------------------------------------
//    // Smoothing functions for density
//    // ----------------------------------------
//    /// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//    virtual double AbstractState::drhodh_constp_smoothed(double xend);
//    /// A smoothed version of the derivative using a spline curve in the region of x=0 to x=xend
//    virtual double AbstractState::drhodp_consth_smoothed(double xend);
//    /// Density corresponding to the smoothed derivatives in the region of x=0 to x=xend
//    virtual void AbstractState::rho_smoothed(double xend, double *rho_spline, double *dsplinedh, double *dsplinedp);

} /* namespace CoolProp */


#ifdef ENABLE_CATCH

#include "catch.hpp"

TEST_CASE("Check AbstractState","[AbstractState]")
{
    SECTION("bad backend")
    {
        CHECK_THROWS(shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("DEFINITELY_A_BAD_BACKEND", "Water")));
    }
    SECTION("good backend - bad fluid")
    {
        CHECK_THROWS(shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("HEOS", "DEFINITELY_A_BAD_FLUID")));
    }
    SECTION("good backend - helmholtz")
    {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("HEOS", "Water")));
    }
    SECTION("good backend - incomp")
    {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("INCOMP", "DEB")));
    }
    SECTION("good backend - REFPROP")
    {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("REFPROP", "Water")));
    }

}

#endif
