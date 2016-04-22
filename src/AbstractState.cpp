/*
 * AbstractState.cpp
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>
#include "math.h"
#include "AbstractState.h"
#include "Backends/REFPROP/REFPROPBackend.h"
#include "Backends/Helmholtz/HelmholtzEOSBackend.h"
#include "Backends/Incompressible/IncompressibleBackend.h"
#include "Backends/Helmholtz/Fluids/FluidLibrary.h"
#include "Backends/IF97/IF97Backend.h"
#include "Backends/Cubics/CubicBackend.h"

#if !defined(NO_TABULAR_BACKENDS)
    #include "Backends/Tabular/TTSEBackend.h"
    #include "Backends/Tabular/BicubicBackend.h"
#endif

namespace CoolProp {

AbstractState * AbstractState::factory(const std::string &backend, const std::vector<std::string> &fluid_names)
{
    if (get_debug_level() > 0){
        std::cout << "AbstractState::factory(" << backend << "," << stringvec_to_string(fluid_names) << ")" << std::endl;
    }
    static const std::string HEOS_string = "HEOS";
    if (!backend.compare(HEOS_string))
    {
        if (fluid_names.size() == 1){
            return new HelmholtzEOSBackend(fluid_names[0]);
        }
        else{
            return new HelmholtzEOSMixtureBackend(fluid_names);
        }
    }
    else if (!backend.compare("REFPROP"))
    {
        if (fluid_names.size() == 1){
            return new REFPROPBackend(fluid_names[0]);
        }
        else{
            return new REFPROPMixtureBackend(fluid_names);
        }
    }
    else if (!backend.compare("INCOMP"))
    {
        if (fluid_names.size() != 1){throw ValueError(format("For INCOMP backend, name vector must be one element long"));}
        return new IncompressibleBackend(fluid_names[0]);
    }
    else if (!backend.compare("IF97"))
    {
        return new IF97Backend();
    }
    #if !defined(NO_TABULAR_BACKENDS)
    else if (backend.find("TTSE&") == 0)
    {
        // Will throw if there is a problem with this backend
        shared_ptr<AbstractState> AS(factory(backend.substr(5), fluid_names));
        return new TTSEBackend(AS);
    }
    else if (backend.find("BICUBIC&") == 0)
    {
        // Will throw if there is a problem with this backend
        shared_ptr<AbstractState> AS(factory(backend.substr(8), fluid_names));
        return new BicubicBackend(AS);
    }
    #endif
    else if (!backend.compare("TREND"))
    {
        throw ValueError("TREND backend not yet implemented");
    }
    else if (backend == "SRK"){
        return new SRKBackend(fluid_names, get_config_double(R_U_CODATA));
    }
    else if (backend == "PR" || backend == "Peng-Robinson"){
        return new PengRobinsonBackend(fluid_names, get_config_double(R_U_CODATA));
    }
    else if (!backend.compare("?") || backend.empty())
    {
        std::size_t idel = fluid_names[0].find("::");
        // Backend has not been specified, and we have to figure out what the backend is by parsing the string
        if (idel == std::string::npos) // No '::' found, no backend specified, try HEOS, otherwise a failure
        {
            // Figure out what backend to use
            return factory(HEOS_string, fluid_names);
        }
        else
        {
            // Split string at the '::' into two std::string, call again
            return factory(std::string(fluid_names[0].begin(), fluid_names[0].begin() + idel), std::string(fluid_names[0].begin()+(idel+2), fluid_names[0].end()));
        }
    }
    else
    {
        throw ValueError(format("Invalid backend name [%s] to factory function",backend.c_str()));
    }
}
std::vector<std::string> AbstractState::fluid_names(void)
{
    return calc_fluid_names();
}

bool AbstractState::clear() {
    // Reset all instances of CachedElement and overwrite
    // the internal double values with -_HUGE
    this->_R = _HUGE;
    this->_gas_constant.clear();
    this->_molar_mass.clear();

    /// Ancillary curve values
    this->_rhoLanc.clear();
    this->_rhoVanc.clear();
    this->_pVanc.clear();
    this->_pLanc.clear();
    this->_TVanc.clear();
    this->_TLanc.clear();

    this->_critical.fill(_HUGE);
    this->_reducing.fill(_HUGE);

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

    /// Two-Phase variables
    this->_rhoLmolar.clear();
    this->_rhoVmolar.clear();

    /// Transport properties
    this->_viscosity.clear();
    this->_conductivity.clear();
    this->_surface_tension.clear();

    return true;
}
void AbstractState::mass_to_molar_inputs(CoolProp::input_pairs &input_pair, CoolPropDbl &value1, CoolPropDbl &value2)
{
    // Check if a mass based input, convert it to molar units

    switch (input_pair)
    {
    case DmassT_INPUTS: ///< Mass density in kg/m^3, Temperature in K
        //case HmassT_INPUTS: ///< Enthalpy in J/kg, Temperature in K (NOT CURRENTLY IMPLEMENTED)
    case SmassT_INPUTS: ///< Entropy in J/kg/K, Temperature in K
        //case TUmass_INPUTS: ///< Temperature in K, Internal energy in J/kg (NOT CURRENTLY IMPLEMENTED)
    case DmassP_INPUTS: ///< Mass density in kg/m^3, Pressure in Pa
    case DmassQ_INPUTS: ///< Mass density in kg/m^3, molar quality
    case HmassP_INPUTS: ///< Enthalpy in J/kg, Pressure in Pa
    case PSmass_INPUTS: ///< Pressure in Pa, Entropy in J/kg/K
    case PUmass_INPUTS: ///< Pressure in Pa, Internal energy in J/kg
    case HmassSmass_INPUTS: ///< Enthalpy in J/kg, Entropy in J/kg/K
    case SmassUmass_INPUTS: ///< Entropy in J/kg/K, Internal energy in J/kg
    case DmassHmass_INPUTS: ///< Mass density in kg/m^3, Enthalpy in J/kg
    case DmassSmass_INPUTS: ///< Mass density in kg/m^3, Entropy in J/kg/K
    case DmassUmass_INPUTS: ///< Mass density in kg/m^3, Internal energy in J/kg
    {
        // Set the cache value for the molar mass if it hasn't been set yet
        molar_mass();

        // Molar mass (just for compactness of the following switch)
        CoolPropDbl mm = static_cast<CoolPropDbl>(_molar_mass);

        switch (input_pair)
        {
        case DmassT_INPUTS: input_pair = DmolarT_INPUTS; value1 /= mm;  break;
            //case HmassT_INPUTS: input_pair = HmolarT_INPUTS; value1 *= mm;  break; (NOT CURRENTLY IMPLEMENTED)
        case SmassT_INPUTS: input_pair = SmolarT_INPUTS; value1 *= mm;  break;
            //case TUmass_INPUTS: input_pair = TUmolar_INPUTS; value2 *= mm;  break; (NOT CURRENTLY IMPLEMENTED)
        case DmassP_INPUTS: input_pair = DmolarP_INPUTS; value1 /= mm;  break;
        case DmassQ_INPUTS: input_pair = DmolarQ_INPUTS; value1 /= mm;  break;
        case HmassP_INPUTS: input_pair = HmolarP_INPUTS; value1 *= mm;  break;
        case PSmass_INPUTS: input_pair = PSmolar_INPUTS; value2 *= mm;  break;
        case PUmass_INPUTS: input_pair = PUmolar_INPUTS; value2 *= mm;  break;
        case HmassSmass_INPUTS: input_pair = HmolarSmolar_INPUTS; value1 *= mm; value2 *= mm;  break;
        case SmassUmass_INPUTS: input_pair = SmolarUmolar_INPUTS; value1 *= mm; value2 *= mm;  break;
        case DmassHmass_INPUTS: input_pair = DmolarHmolar_INPUTS; value1 /= mm; value2 *= mm;  break;
        case DmassSmass_INPUTS: input_pair = DmolarSmolar_INPUTS; value1 /= mm; value2 *= mm;  break;
        case DmassUmass_INPUTS: input_pair = DmolarUmolar_INPUTS; value1 /= mm; value2 *= mm;  break;
        default: break;
        }
        break;
    }
    default:
        return;
    }
}
double AbstractState::trivial_keyed_output(parameters key)
{
    if (get_debug_level()>=50) std::cout << format("AbstractState: keyed_output called for %s ",get_parameter_information(key,"short").c_str()) << std::endl;
    switch (key)
    {
    case imolar_mass:
        return molar_mass();
    case iacentric_factor:
        return acentric_factor();
	case igas_constant:
		return gas_constant();
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
        return calc_T_reducing();
    case irhomolar_reducing:
        return calc_rhomolar_reducing();
	case iP_reducing:
		return calc_p_reducing();
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
    case iFH:
        return this->calc_flame_hazard();
    case iHH:
        return this->calc_health_hazard();
    case iPH:
        return this->calc_physical_hazard();
    case idipole_moment:
        return this->calc_dipole_moment();
    default:
        throw ValueError(format("This input [%d: \"%s\"] is not valid for trivial_keyed_output",key,get_parameter_information(key,"short").c_str()));
    }
}
double AbstractState::keyed_output(parameters key)
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
    case iSmolar_residual:
        return gas_constant()*(tau()*dalphar_dTau() - alphar());
    case iSmass:
        return smass();
    case iUmolar:
        return umolar();
    case iUmass:
        return umass();
    case iGmolar:
        return gibbsmolar();
    case iGmass:
        return gibbsmass();
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
    case ialphar:
        return alphar(); 
    case ialpha0:
        return alpha0();
    case idalpha0_ddelta_consttau:
        return dalpha0_dDelta();
    case idalpha0_dtau_constdelta:
        return dalpha0_dTau();
    case idalphar_ddelta_consttau:
        return dalphar_dDelta();
    case idalphar_dtau_constdelta:
        return dalphar_dTau();
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
    case iisobaric_expansion_coefficient:
        return isobaric_expansion_coefficient();
    case iviscosity:
        return viscosity();
    case iconductivity:
        return conductivity();
    case iPrandtl:
        return Prandtl();
	case isurface_tension:
		return surface_tension();
    case iPhase:
        return phase();
    case iZ:
        return compressibility_factor();
    case iPIP:
        return PIP();
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
double AbstractState::gibbsmolar(void){
    if (!_gibbsmolar) _gibbsmolar = calc_gibbsmolar();
    return _gibbsmolar;
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
double AbstractState::acentric_factor(){
    return calc_acentric_factor();
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
double AbstractState::fugacity_coefficient(std::size_t i){
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity_coefficient(i);
}
double AbstractState::fugacity(std::size_t i){
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity(i);
}
double AbstractState::chemical_potential(std::size_t i) {
    // TODO: Cache the chemical potential for each component
    return calc_chemical_potential(i);
}
void AbstractState::build_phase_envelope(const std::string &type)
{
    calc_phase_envelope(type);
}
double AbstractState::isothermal_compressibility(void){
    return 1.0/_rhomolar*first_partial_deriv(iDmolar, iP, iT);
}
double AbstractState::isobaric_expansion_coefficient(void){
    return -1.0/_rhomolar*first_partial_deriv(iDmolar, iT, iP);
}
double AbstractState::Bvirial(void){ return calc_Bvirial(); }
double AbstractState::Cvirial(void){ return calc_Cvirial(); }
double AbstractState::dBvirial_dT(void){ return calc_dBvirial_dT(); }
double AbstractState::dCvirial_dT(void){ return calc_dCvirial_dT(); }
double AbstractState::compressibility_factor(void){ return calc_compressibility_factor(); }

// Get the derivatives of the parameters in the partial derivative with respect to T and rho
void get_dT_drho(AbstractState &AS, parameters index, CoolPropDbl &dT, CoolPropDbl &drho)
{
    CoolPropDbl T = AS.T(),
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
    case iGmass:
    case iGmolar:
    {
        // dg/dT|rho
        double dTau_dT = 1/dT_dtau;
        dT = R*AS.T()*(AS.dalpha0_dTau()+AS.dalphar_dTau()+AS.delta()*AS.d2alphar_dDelta_dTau())*dTau_dT + R*(1+AS.alpha0() + AS.alphar() + AS.delta()*AS.dalphar_dDelta());
        // dg/drho|T
        double dDelta_drho = 1/rhor;
        drho = AS.T()*R*(AS.dalpha0_dDelta()+AS.dalphar_dDelta()+AS.delta()*AS.d2alphar_dDelta2() + AS.dalphar_dDelta())*dDelta_drho;
        if (index == iGmass){
            // dg/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
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
void get_dT_drho_second_derivatives(AbstractState &AS, int index, CoolPropDbl &dT2, CoolPropDbl &drho_dT, CoolPropDbl &drho2)
{
	CoolPropDbl T = AS.T(),
                rho = AS.rhomolar(),
                rhor = AS.rhomolar_reducing(),
                Tr = AS.T_reducing(),
                R = AS.gas_constant(),
                delta = rho/rhor,
                tau = Tr/T;

    // Here we use T and rho as independent variables since derivations are already done by Thorade, 2013,
    // Partial derivatives of thermodynamic state propertiesfor dynamic simulation, DOI 10.1007/s12665-013-2394-z

    switch (index)
    {
    case iT:
	case iDmass:
    case iDmolar:
        dT2 = 0; // d2rhomolar_dtau2
        drho2 = 0;
        drho_dT = 0;
        break;
    case iTau:
        dT2 = 2*Tr/pow(T,3); drho_dT = 0; drho2 = 0; break;
    case iDelta:
        dT2 = 0; drho_dT = 0; drho2 = 0; break;
    case iP:
    {
        drho2 = T*R/rho*(2*delta*AS.dalphar_dDelta()+4*pow(delta,2)*AS.d2alphar_dDelta2()+pow(delta,3)*AS.d3alphar_dDelta3());
        dT2 = rho*R/T*(pow(tau,2)*delta*AS.d3alphar_dDelta_dTau2());
        drho_dT = R*(1+2*delta*AS.dalphar_dDelta() +pow(delta,2)*AS.d2alphar_dDelta2() - 2*delta*tau*AS.d2alphar_dDelta_dTau() - tau*pow(delta, 2)*AS.d3alphar_dDelta2_dTau());
        break;
    }
    case iHmass:
    case iHmolar:
    {
        // d2h/drho2|T
        drho2 = R*T*pow(delta/rho,2)*(tau*AS.d3alphar_dDelta2_dTau() + 2*AS.d2alphar_dDelta2() + delta*AS.d3alphar_dDelta3());
        // d2h/dT2|rho
        dT2 = R/T*pow(tau, 2)*(tau*(AS.d3alpha0_dTau3()+AS.d3alphar_dTau3()) + 2*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()) + delta*AS.d3alphar_dDelta_dTau2());
        // d2h/drho/dT
        drho_dT = R/rho*delta*(delta*AS.d2alphar_dDelta2() - pow(tau,2)*AS.d3alphar_dDelta_dTau2() + AS.dalphar_dDelta() - tau*delta*AS.d3alphar_dDelta2_dTau() - tau*AS.d2alphar_dDelta_dTau());
        if (index == iHmass){
            drho2 /= AS.molar_mass();
            drho_dT /= AS.molar_mass();
            dT2 /= AS.molar_mass();
        }
        break;
    }
    case iSmass:
    case iSmolar:
    {
        // d2s/rho2|T
        drho2 = R/pow(rho,2)*(1-pow(delta,2)*AS.d2alphar_dDelta2() + tau*pow(delta,2)*AS.d3alphar_dDelta2_dTau());
        // d2s/dT2|rho
        dT2 = R*pow(tau/T, 2)*(tau*(AS.d3alpha0_dTau3()+AS.d3alphar_dTau3())+3*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()));
        // d2s/drho/dT
        drho_dT = R/(T*rho)*(-pow(tau,2)*delta*AS.d3alphar_dDelta_dTau2());
        if (index == iSmass){
            drho2 /= AS.molar_mass();
            drho_dT /= AS.molar_mass();
            dT2 /= AS.molar_mass();
        }
        break;
    }
    case iUmass:
    case iUmolar:
    {
        // d2u/rho2|T
        drho2 = R*T*tau*pow(delta/rho,2)*AS.d3alphar_dDelta2_dTau();
        // d2u/dT2|rho
        dT2 = R/T*pow(tau, 2)*(tau*(AS.d3alpha0_dTau3()+AS.d3alphar_dTau3())+2*(AS.d2alpha0_dTau2()+AS.d2alphar_dTau2()));
        // d2u/drho/dT
        drho_dT = R/rho*(-pow(tau,2)*delta*AS.d3alphar_dDelta_dTau2());
        if (index == iUmass){
            drho2 /= AS.molar_mass();
            drho_dT /= AS.molar_mass();
            dT2 /= AS.molar_mass();
        }
        break;
    }
    default:
        throw ValueError(format("input to get_dT_drho_second_derivatives[%s] is invalid", get_parameter_information(index,"short").c_str()));
    }
}
CoolPropDbl AbstractState::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant)
{
	CoolPropDbl dOf_dT, dOf_drho, dWrt_dT, dWrt_drho, dConstant_dT, dConstant_drho;

    get_dT_drho(*this, Of, dOf_dT, dOf_drho);
    get_dT_drho(*this, Wrt, dWrt_dT, dWrt_drho);
    get_dT_drho(*this, Constant, dConstant_dT, dConstant_drho);

    return (dOf_dT*dConstant_drho-dOf_drho*dConstant_dT)/(dWrt_dT*dConstant_drho-dWrt_drho*dConstant_dT);
}
CoolPropDbl AbstractState::calc_second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2)
{
    CoolPropDbl dOf1_dT, dOf1_drho, dWrt1_dT, dWrt1_drho, dConstant1_dT, dConstant1_drho, d2Of1_dT2, d2Of1_drhodT,
                d2Of1_drho2, d2Wrt1_dT2, d2Wrt1_drhodT, d2Wrt1_drho2, d2Constant1_dT2, d2Constant1_drhodT, d2Constant1_drho2,
                dWrt2_dT, dWrt2_drho, dConstant2_dT, dConstant2_drho, N, D, dNdrho__T, dDdrho__T, dNdT__rho, dDdT__rho,
                dderiv1_drho, dderiv1_dT, second;

    // First and second partials needed for terms involved in first derivative
    get_dT_drho(*this, Of1, dOf1_dT, dOf1_drho);
    get_dT_drho(*this, Wrt1, dWrt1_dT, dWrt1_drho);
    get_dT_drho(*this, Constant1, dConstant1_dT, dConstant1_drho);
    get_dT_drho_second_derivatives(*this, Of1, d2Of1_dT2, d2Of1_drhodT, d2Of1_drho2);
    get_dT_drho_second_derivatives(*this, Wrt1, d2Wrt1_dT2, d2Wrt1_drhodT, d2Wrt1_drho2);
    get_dT_drho_second_derivatives(*this, Constant1, d2Constant1_dT2, d2Constant1_drhodT, d2Constant1_drho2);

    // First derivatives of terms involved in the second derivative
    get_dT_drho(*this, Wrt2, dWrt2_dT, dWrt2_drho);
    get_dT_drho(*this, Constant2, dConstant2_dT, dConstant2_drho);

    // Numerator and denominator of first partial derivative term
    N = dOf1_dT*dConstant1_drho - dOf1_drho*dConstant1_dT;
    D = dWrt1_dT*dConstant1_drho - dWrt1_drho*dConstant1_dT;

    // Derivatives of the numerator and denominator of the first partial derivative term with respect to rho, T held constant
    // They are of similar form, with Of1 and Wrt1 swapped
    dNdrho__T = dOf1_dT*d2Constant1_drho2 + d2Of1_drhodT*dConstant1_drho - dOf1_drho*d2Constant1_drhodT - d2Of1_drho2*dConstant1_dT;
    dDdrho__T = dWrt1_dT*d2Constant1_drho2 + d2Wrt1_drhodT*dConstant1_drho - dWrt1_drho*d2Constant1_drhodT - d2Wrt1_drho2*dConstant1_dT;

    // Derivatives of the numerator and denominator of the first partial derivative term with respect to T, rho held constant
    // They are of similar form, with Of1 and Wrt1 swapped
    dNdT__rho = dOf1_dT*d2Constant1_drhodT + d2Of1_dT2*dConstant1_drho - dOf1_drho*d2Constant1_dT2 - d2Of1_drhodT*dConstant1_dT;
    dDdT__rho = dWrt1_dT*d2Constant1_drhodT + d2Wrt1_dT2*dConstant1_drho - dWrt1_drho*d2Constant1_dT2 - d2Wrt1_drhodT*dConstant1_dT;

    // First partial of first derivative term with respect to T
    dderiv1_drho = (D*dNdrho__T - N*dDdrho__T)/pow(D, 2);

    // First partial of first derivative term with respect to rho
    dderiv1_dT = (D*dNdT__rho - N*dDdT__rho)/pow(D, 2);

    // Complete second derivative
    second = (dderiv1_dT*dConstant2_drho - dderiv1_drho*dConstant2_dT)/(dWrt2_dT*dConstant2_drho - dWrt2_drho*dConstant2_dT);

    return second;
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

TEST_CASE("Check derivatives in first_partial_deriv","[derivs_in_first_partial_deriv]")
{
    shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> WaterplusT(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> WaterminusT(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> Waterplusrho(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> Waterminusrho(CoolProp::AbstractState::factory("HEOS", "Water"));

    double dT = 1e-3, drho = 1;
    Water->update(CoolProp::PT_INPUTS, 101325, 300);
    WaterplusT->update(CoolProp::DmolarT_INPUTS, Water->rhomolar(), 300+dT);
    WaterminusT->update(CoolProp::DmolarT_INPUTS, Water->rhomolar(), 300-dT);
    Waterplusrho->update(CoolProp::DmolarT_INPUTS, Water->rhomolar()+drho, 300);
    Waterminusrho->update(CoolProp::DmolarT_INPUTS, Water->rhomolar()-drho, 300);

    // Numerical derivatives
    CoolPropDbl dP_dT_num = (WaterplusT->p() - WaterminusT->p())/(2*dT);
    CoolPropDbl dP_drho_num = (Waterplusrho->p() - Waterminusrho->p())/(2*drho);
    
    CoolPropDbl dHmolar_dT_num = (WaterplusT->hmolar() - WaterminusT->hmolar())/(2*dT);
    CoolPropDbl dHmolar_drho_num = (Waterplusrho->hmolar() - Waterminusrho->hmolar())/(2*drho);
    CoolPropDbl dHmass_dT_num = (WaterplusT->hmass() - WaterminusT->hmass())/(2*dT);
    CoolPropDbl dHmass_drho_num = (Waterplusrho->hmass() - Waterminusrho->hmass())/(2*drho);

    CoolPropDbl dSmolar_dT_num = (WaterplusT->smolar() - WaterminusT->smolar())/(2*dT);
    CoolPropDbl dSmolar_drho_num = (Waterplusrho->smolar() - Waterminusrho->smolar())/(2*drho);
    CoolPropDbl dSmass_dT_num = (WaterplusT->smass() - WaterminusT->smass())/(2*dT);
    CoolPropDbl dSmass_drho_num = (Waterplusrho->smass() - Waterminusrho->smass())/(2*drho);

    CoolPropDbl dUmolar_dT_num = (WaterplusT->umolar() - WaterminusT->umolar())/(2*dT);
    CoolPropDbl dUmolar_drho_num = (Waterplusrho->umolar() - Waterminusrho->umolar())/(2*drho);
    CoolPropDbl dUmass_dT_num = (WaterplusT->umass() - WaterminusT->umass())/(2*dT);
    CoolPropDbl dUmass_drho_num = (Waterplusrho->umass() - Waterminusrho->umass())/(2*drho);

    CoolPropDbl dGmolar_dT_num = (WaterplusT->gibbsmolar() - WaterminusT->gibbsmolar())/(2*dT);
    CoolPropDbl dGmolar_drho_num = (Waterplusrho->gibbsmolar() - Waterminusrho->gibbsmolar())/(2*drho);
    CoolPropDbl dGmass_dT_num = (WaterplusT->gibbsmass() - WaterminusT->gibbsmass())/(2*dT);
    CoolPropDbl dGmass_drho_num = (Waterplusrho->gibbsmass() - Waterminusrho->gibbsmass())/(2*drho);

    // Pressure
    CoolPropDbl dP_dT_analyt, dP_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iP, dP_dT_analyt, dP_drho_analyt);
    // Enthalpy
    CoolPropDbl dHmolar_dT_analyt, dHmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iHmolar, dHmolar_dT_analyt, dHmolar_drho_analyt);
    CoolPropDbl dHmass_dT_analyt, dHmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iHmass, dHmass_dT_analyt, dHmass_drho_analyt);
    // Entropy
    CoolPropDbl dSmolar_dT_analyt, dSmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iSmolar, dSmolar_dT_analyt, dSmolar_drho_analyt);
    CoolPropDbl dSmass_dT_analyt, dSmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iSmass, dSmass_dT_analyt, dSmass_drho_analyt);
    // Internal energy
    CoolPropDbl dUmolar_dT_analyt, dUmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iUmolar, dUmolar_dT_analyt, dUmolar_drho_analyt);
    CoolPropDbl dUmass_dT_analyt, dUmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iUmass, dUmass_dT_analyt, dUmass_drho_analyt);
    // Gibbs
    CoolPropDbl dGmolar_dT_analyt, dGmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iGmolar, dGmolar_dT_analyt, dGmolar_drho_analyt);
    CoolPropDbl dGmass_dT_analyt, dGmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iGmass, dGmass_dT_analyt, dGmass_drho_analyt);

    double eps = 1e-3;

    CHECK( std::abs(dP_dT_analyt/dP_dT_num-1) < eps);
    CHECK( std::abs(dP_drho_analyt/dP_drho_num-1) < eps);
    
    CHECK( std::abs(dHmolar_dT_analyt/dHmolar_dT_num-1) < eps);
    CHECK( std::abs(dHmolar_drho_analyt/dHmolar_drho_num-1) < eps);
    CHECK( std::abs(dHmass_dT_analyt/dHmass_dT_num-1) < eps);
    CHECK( std::abs(dHmass_drho_analyt/dHmass_drho_num-1) < eps);

    CHECK( std::abs(dSmolar_dT_analyt/dSmolar_dT_num-1) < eps);
    CHECK( std::abs(dSmolar_drho_analyt/dSmolar_drho_num-1) < eps);
    CHECK( std::abs(dSmass_dT_analyt/dSmass_dT_num-1) < eps);
    CHECK( std::abs(dSmass_drho_analyt/dSmass_drho_num-1) < eps);

    CHECK( std::abs(dUmolar_dT_analyt/dUmolar_dT_num-1) < eps);
    CHECK( std::abs(dUmolar_drho_analyt/dUmolar_drho_num-1) < eps);
    CHECK( std::abs(dUmass_dT_analyt/dUmass_dT_num-1) < eps);
    CHECK( std::abs(dUmass_drho_analyt/dUmass_drho_num-1) < eps);

    CHECK( std::abs(dGmolar_dT_analyt/dGmolar_dT_num-1) < eps);
    CHECK( std::abs(dGmolar_drho_analyt/dGmolar_drho_num-1) < eps);
    CHECK( std::abs(dGmass_dT_analyt/dGmass_dT_num-1) < eps);
    CHECK( std::abs(dGmass_drho_analyt/dGmass_drho_num-1) < eps);

}

#endif


/// *********************************************************************************
/// *********************************************************************************
///                     EMSCRIPTEN (for javascript)
/// *********************************************************************************
/// *********************************************************************************

#ifdef EMSCRIPTEN

#include <emscripten/bind.h>
using namespace emscripten;

CoolProp::AbstractState * factory(const std::string &backend, const std::string &fluid_names)
{
    return CoolProp::AbstractState::factory(backend, strsplit(fluid_names, '&'));
}

// Binding code
EMSCRIPTEN_BINDINGS(abstract_state_bindings) {

    register_vector<double>("VectorDouble");
    register_vector<std::string>("VectorString");

    //value_object<CoolProp::PhaseEnvelopeData>("CoolProp::PhaseEnvelopeData")
    //    // Use X macros to auto-generate the variables; 
    //    // each will look something like: .field("T", &CoolProp::PhaseEnvelopeData::T);
    //    #define X(name) .field("name", &CoolProp::PhaseEnvelopeData::name)
    //        PHASE_ENVELOPE_VECTORS
    //    #undef X
    //    ;

    function("factory", &factory, allow_raw_pointers());

    class_<CoolProp::AbstractState>("AbstractState")
        .function("gas_constant", &CoolProp::AbstractState::gas_constant)
        .function("update", &CoolProp::AbstractState::update)
        .function("set_mole_fractions", &CoolProp::AbstractState::set_mole_fractions_double)
        .function("build_phase_envelope", &CoolProp::AbstractState::build_phase_envelope)
        .function("get_phase_envelope_data", &CoolProp::AbstractState::get_phase_envelope_data);

}

#endif