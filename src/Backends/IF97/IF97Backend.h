
#ifndef IF97BACKEND_H_
#define IF97BACKEND_H_

#include "DataStructures.h"
#include "externals/IF97/IF97.h"
#include "AbstractState.h"
#include "Exceptions.h"
#include <vector>

namespace CoolProp {

class IF97Backend : public AbstractState  {

protected:

	/// Additional cached elements used only in this backend since the "normal"
    /// backends use only molar units while IF97 uses mass-based units
	CachedElement  _cpmass, _cvmass, _hmass, _rhomass, _smass, _umass;

public:
    /// The name of the backend being used
    std::string backend_name(void) { return get_backend_string(IF97_BACKEND); }

    // REQUIRED BUT NOT USED IN IF97 FUNCTIONS
    bool using_mole_fractions(void){return false;};
    bool using_mass_fractions(void){return true;}; // But actually it doesn't matter since it is only a pure fluid
    bool using_volu_fractions(void){return false;};
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions){throw NotImplementedError("Mole composition has not been implemented.");};
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){}; // Not implemented, but don't throw any errors
    void set_volu_fractions(const std::vector<CoolPropDbl> &volu_fractions){throw NotImplementedError("Volume composition has not been implemented.");};
    const std::vector<CoolPropDbl> & get_mole_fractions(void){throw NotImplementedError("get_mole_fractions composition has not been implemented.");};

    /// Override clear() function of IF97 Water
    bool clear() {
    // Reset all instances of CachedElement and overwrite
    // the internal double values with -_HUGE
        this->_rhomolar = -_HUGE;
        this->_T = -_HUGE;
        this->_p = -_HUGE;
        this->_Q = -_HUGE;
        this->_phase = iphase_not_imposed; // Default condition is no phse imposed

        return true;
    };

    /// Updating function for IF97 Water
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update the bare minimum values needed to calculate the other values.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs input_pair, double value1, double value2){
        
        clear();  //clear the few cached values we are using

        switch(input_pair){
            case PT_INPUTS: _p = value1; _T = value2; break;
            case PQ_INPUTS: 
                _p = value1; 
                _Q = value2; 
                _T = IF97::Tsat97(_p); 
                if (std::abs(_Q) < 1e-10){
                    _phase = iphase_liquid; // bubble point
                }
                else if (std::abs(value2-1) < 1e-10){
                    _phase = iphase_gas; // dew point
                }
                else
                    _phase = iphase_twophase;
                break;
            case QT_INPUTS: 
                _Q = value1; 
                _T = value2; 
                _p = IF97::psat97(_T); 
                if (std::abs(_Q) < 1e-10){
                    _phase = iphase_liquid; // bubble point
                }
                else if (std::abs(value2-1) < 1e-10){
                    _phase = iphase_gas; // dew point
                }
                else
                    _phase = iphase_twophase;
                break;
            default:
                throw ValueError("Bad input_pair");
        }
    };

    /*  We have to override some of the functions from the AbstractState.
	 *  IF97 is only mass-based and does not support conversion
	 *  from mass- to molar-specific quantities.
	 */
    // ************************************************************************* //
    //                   Basic Thermodynamic Functions                           //
    // ************************************************************************* //
    //
    double calc_Liquid(parameters iCalc) {
        switch(iCalc){
        case iDmass:  return IF97::rholiq_p(_p); break; ///< Mass-based density
        case iHmass:  return IF97::hliq_p(_p); break;   ///< Mass-based enthalpy
        case iSmass:  return IF97::sliq_p(_p); break;   ///< Mass-based entropy
        case iCpmass: return IF97::cpliq_p(_p); break;  ///< Mass-based constant-pressure specific heat
        case iCvmass: return IF97::cvliq_p(_p); break;  ///< Mass-based constant-volume specific heat
        case iUmass:  return IF97::uliq_p(_p); break;   ///< Mass-based internal energy
        case ispeed_sound: return IF97::speed_soundliq_p(_p); break;  ///< Speed of Sound
        default: return -_HUGE;
        };
    }
    double calc_Vapor(parameters iCalc) {
        switch(iCalc){
        case iDmass:  return IF97::rhovap_p(_p); break; ///< Mass-based density
        case iHmass:  return IF97::hvap_p(_p); break;   ///< Mass-based enthalpy
        case iSmass:  return IF97::svap_p(_p); break;   ///< Mass-based entropy
        case iCpmass: return IF97::cpvap_p(_p); break;  ///< Mass-based constant-pressure specific heat
        case iCvmass: return IF97::cvvap_p(_p); break;  ///< Mass-based constant-volume specific heat
        case iUmass:  return IF97::uvap_p(_p); break;   ///< Mass-based internal energy
        case ispeed_sound: return IF97::speed_soundvap_p(_p); break;  ///< Speed of Sound
        default: return -_HUGE;
        };
    }
    double calc_Flash(parameters iCalc) {
        switch(_phase){
            case iphase_liquid: 
                return calc_Liquid(iCalc); 
                break;
            case iphase_gas:    
                return calc_Vapor(iCalc); 
                break;
            case iphase_twophase:
                switch(iCalc){ 
                case iDmass:
                    return 1.0/(_Q/calc_Vapor(iCalc) + (1.0-_Q)/calc_Liquid(iCalc)); break;
                case iCpmass:
                    throw NotImplementedError(format("Isobaric Specific Heat not implemented in two phase region")); break;
                case iCvmass:
                    throw NotImplementedError(format("Isochoric Specific Heat not implemented in two phase region")); break;
                case ispeed_sound:
                    throw NotImplementedError(format("Speed of Sound not implemented in two phase region")); break;
                default:
                    return _Q*calc_Vapor(iCalc) + (1.0-_Q)*calc_Liquid(iCalc);
                }; break;
            default: 
                switch(iCalc){
                case iDmass:  return IF97::rhomass_Tp(_T, _p); break; ///< Mass-based density
                case iHmass:  return IF97::hmass_Tp(_T, _p); break;   ///< Mass-based enthalpy
                case iSmass:  return IF97::smass_Tp(_T, _p); break;   ///< Mass-based entropy
                case iCpmass: return IF97::cpmass_Tp(_T, _p); break;  ///< Mass-based constant-pressure specific heat
                case iCvmass: return IF97::cvmass_Tp(_T, _p); break;  ///< Mass-based constant-volume specific heat
                case iUmass:  return IF97::umass_Tp(_T, _p); break;   ///< Mass-based internal energy
                case ispeed_sound: return IF97::speed_sound_Tp(_T, _p); break;  ///< speed of sound
                default: throw NotImplementedError(format("Output variable not implemented in IF97 Backend"));
                };
        }
    }
	/// Return the mass density in kg/m³
    double rhomass(void){ return calc_rhomass(); };
    double calc_rhomass(void){ return calc_Flash(iDmass); };
    /// Return the molar density in mol/m³
    double rhomolar(void){ return calc_rhomolar(); };
    double calc_rhomolar(void){ return rhomass()/molar_mass(); };    /// kg/m³ * mol/kg = mol/m³

	/// Return the mass enthalpy in J/kg
    double hmass(void){ return calc_hmass(); };
    double calc_hmass(void){ return calc_Flash(iHmass); };
	/// Return the molar enthalpy in J/mol
    double hmolar(void){ return calc_hmolar(); };
    double calc_hmolar(void){ return hmass()*molar_mass(); };        /// J/kg * kg/mol = J/mol

	/// Return the mass entropy in J/kg/K
    double smass(void){ return calc_smass(); };
    double calc_smass(void){ return calc_Flash(iSmass); };
	/// Return the molar entropy in J/mol/K
    double smolar(void){ return calc_smolar(); };
    double calc_smolar(void){ return smass()*molar_mass(); };        /// J/kg-K * kg/mol = J/mol-K

	/// Return the mass internal energy in J/kg
    double umass(void){ return calc_umass(); };
    double calc_umass(void){ return calc_Flash(iUmass); };
	/// Return the molar internal energy in J/mol
    double umolar(void){ return calc_umolar(); };
    double calc_umolar(void){ return umass()*molar_mass(); };        /// J/kg * kg/mol = J/mol

	/// Return the mass-based constant pressure specific heat in J/kg/K
    double cpmass(void){ return calc_cpmass(); };
    double calc_cpmass(void){ return calc_Flash(iCpmass); };
	/// Return the molar-based constant pressure specific heat in J/mol/K
    double cpmolar(void){ return calc_cpmolar(); };
    double calc_cpmolar(void){ return cpmass()*molar_mass(); };      /// J/kg-K * kg/mol = J/mol-K

    /// Return the mass-based constant volume specific heat in J/kg/K
    double cvmass(void){ return calc_cvmass(); };
    double calc_cvmass(void){ return calc_Flash(iCvmass); };
    /// Return the molar-based constant volume specific heat in J/mol/K
    double cvmolar(void){ return calc_cvmolar(); };
    double calc_cvmolar(void){ return cvmass()*molar_mass(); };      /// J/kg-K * kg/mol = J/mol-K

    /// Return the speed of sound in m/s
    double speed_sound(void){ return calc_speed_sound(); };
    double calc_speed_sound(void) { return calc_Flash(ispeed_sound); };
    //
    // ************************************************************************* //
    //                         Trivial Functions                                 //
    // ************************************************************************* //
    //
    /// Using this backend, get the triple point temperature in K
    double calc_Ttriple(void){ return IF97::get_Ttrip(); };
    /// Using this backend, get the triple point pressure in Pa
    double calc_p_triple(void){ return IF97::get_ptrip(); };
    /// Using this backend, get the critical point temperature in K
    double calc_T_critical(void){ return IF97::get_Tcrit(); };
    /// Using this backend, get the critical point pressure in Pa
    double calc_p_critical(void){ return IF97::get_pcrit(); };
    /// Using this backend, get the ideal gas constant in J/mol*K
    /// ==> multiplies IF97 Rgas by molar_mass() to put on molar basis per CoolProp convention
    double calc_gas_constant(void){ return IF97::get_Rgas()*molar_mass(); };
    /// Using this backend, get the molar mass in kg/mol
    double calc_molar_mass(void){ return IF97::get_MW(); };
    /// Using this backend, get the acentric factor (unitless)
    double calc_acentric_factor(void){ return IF97::get_Acentric(); };
    /// Using this backend, get the high pressure limit in Pa
    // TODO: May want to adjust this based on _T, since Region 5
    //       is limited to 50 MPa, instead of 100 MPa elsewhere.
    double calc_pmax(void){ return IF97::get_Pmax(); };
    /// Note: Pmin not implemented in Abstract State or CoolProp
    /// Using this backend, get the high temperature limit in K
    double calc_Tmax(void){ return IF97::get_Tmax(); };
    /// Using this backend, get the high pressure limit in K
    double calc_Tmin(void){ return IF97::get_Tmin(); };
    /// Using this backend, get the critical point density in kg/m³
    /// Replace molar-based AbstractState functions since IF97 is mass based only
    double rhomolar_critical(void){
        return calc_rhomass_critical()/molar_mass();
    }
    double rhomass_critical(void){
        return calc_rhomass_critical();
    }
    // Overwrite the virtual calc_ functions for density
    double calc_rhomolar_critical(void){ return rhomass_critical()/molar_mass(); };
    double calc_rhomass_critical(void){ return IF97::get_rhocrit(); };
    //
    // ************************************************************************* //
    //                      Saturation Functions                                 //
    // ************************************************************************* //
    //
    double calc_pressure(void){ return _p; };
    //
    // ************************************************************************* //
    //                      Saturation Functions                                 //
    // ************************************************************************* //
    //
};

} /* namespace CoolProp */
#endif /* IF97BACKEND_H_ */
