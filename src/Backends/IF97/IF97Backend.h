
#ifndef IF97BACKEND_H_
#define IF97BACKEND_H_

#include "DataStructures.h"
#include "externals/IF97/IF97.h"
#include "AbstractState.h"
#include "Exceptions.h"
#include <vector>
#include <cmath>

namespace CoolProp {

class IF97Backend : public AbstractState
{

   protected:
    /// Additional cached elements used only in this backend since the "normal"
    /// backends use only molar units while IF97 uses mass-based units
    CachedElement _hmass, _rhomass, _smass;
    /// CachedElement  _hVmass, _hLmass, _sVmass, sLmass;

   public:
    /// The name of the backend being used
    std::string backend_name(void) {
        return get_backend_string(IF97_BACKEND);
    }

    // REQUIRED BUT NOT USED IN IF97 FUNCTIONS
    bool using_mole_fractions(void) {
        return false;
    };
    bool using_mass_fractions(void) {
        return true;
    };  // But actually it doesn't matter since it is only a pure fluid
    bool using_volu_fractions(void) {
        return false;
    };
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        throw NotImplementedError("Mole composition has not been implemented.");
    };
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions){};  // Not implemented, but don't throw any errors
    void set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions) {
        throw NotImplementedError("Volume composition has not been implemented.");
    };
    const std::vector<CoolPropDbl>& get_mole_fractions(void) {
        throw NotImplementedError("get_mole_fractions composition has not been implemented.");
    };

    /// Override clear() function of IF97 Water
    bool clear() {
        // Reset all instances of CachedElement and overwrite
        // the internal double values with -_HUGE
        // Default phase condition is no phase imposed
        // IF97 will make phase/region determination
        this->_T = -_HUGE;
        this->_p = -_HUGE;
        this->_Q = -_HUGE;
        this->_rhomass.clear();
        this->_hmass.clear();
        this->_smass.clear();
        this->_phase = iphase_not_imposed;
        return true;
    };

    void set_phase() {
        double epsilon = 3.3e-5;                              // IAPWS-IF97 RMS saturated pressure inconsistency
        if ((std::abs(_T - IF97::Tcrit) < epsilon / 10.0) &&  //            RMS temperature inconsistency ~ epsilon/10
            (std::abs(_p - IF97::Pcrit) < epsilon)) {         // within epsilon of [Tcrit,Pcrit]
            _phase = iphase_critical_point;                   //     at critical point
        } else if (_T >= IF97::Tcrit) {                       // to the right of the critical point
            if (_p >= IF97::Pcrit) {                          //     above the critical point
                _phase = iphase_supercritical;
            } else {  //     below the critical point
                _phase = iphase_supercritical_gas;
            }
        } else {                      // to the left of the critical point
            if (_p >= IF97::Pcrit) {  //     above the critical point
                _phase = iphase_supercritical_liquid;
            } else {  //     below critical point
                double psat = IF97::psat97(_T);
                if (_p > psat * (1.0 + epsilon)) {  //         above the saturation curve
                    _phase = iphase_liquid;
                } else if (_p < psat * (1.0 - epsilon)) {  //         below the saturation curve
                    _phase = iphase_gas;
                } else  //         exactly on saturation curve (within 1e-4 %)
                    _phase = iphase_twophase;
            }
        }
    };

    /// Updating function for IF97 Water
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update the bare minimum values needed to calculate the other values.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs input_pair, double value1, double value2) {

        double H, S, hLmass, hVmass, sLmass, sVmass;

        clear();  //clear the few cached values we are using

        switch (input_pair) {
            case PT_INPUTS:
                _p = value1;
                _T = value2;
                _Q = -1;
                set_phase();
                //Two-Phase Check, with PT Inputs:
                if (_phase == iphase_twophase)
                    throw ValueError(format("Saturation pressure [%g Pa] corresponding to T [%g K] is within 3.3e-3 %% of given p [%Lg Pa]",
                                            IF97::psat97(_T), _T, _p));
                break;
            case PQ_INPUTS:
                _p = value1;
                _Q = value2;
                if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
                _T = IF97::Tsat97(_p);  // ...will throw exception if _P not on saturation curve
                _phase = iphase_twophase;
                break;
            case QT_INPUTS:
                _Q = value1;
                _T = value2;
                if ((_Q < 0) || (_Q > 1)) throw CoolProp::OutOfRangeError("Input vapor quality [Q] must be between 0 and 1");
                _p = IF97::psat97(_T);  // ...will throw exception if _P not on saturation curve
                _phase = iphase_twophase;
                break;
            case HmolarP_INPUTS:
                // IF97 is mass based so convert hmolar input to hmass
                _hmass = value1 / molar_mass();  // Convert to mass basis : (J/mol) / (kg/mol) = J/kg
                _p = value2;
                // Fall thru to mass basis inputs
            case HmassP_INPUTS:
                if (!(_hmass)) _hmass = value1;  // don't set if already set above
                _p = value2;
                _T = IF97::T_phmass(_p, _hmass);
                // ...if in the vapor dome (Region 4), calculate Quality...
                if (IF97::BackwardRegion(_p, _hmass, IF97_HMASS) == 4) {
                    H = _hmass;
                    hVmass = IF97::hvap_p(_p);
                    hLmass = IF97::hliq_p(_p);
                    _Q = (H - hLmass) / (hVmass - hLmass);
                    _phase = iphase_twophase;
                } else {
                    _Q = -1;
                    set_phase();
                };
                break;
            case PSmolar_INPUTS:
                // IF97 is mass based so convert smolar input to smass
                _p = value1;
                _smass = value2 / molar_mass();  // Convert to mass basis : (J/mol-K) / (kg/mol) = J/kg-K
                // Fall thru to mass basis inputs
            case PSmass_INPUTS:
                _p = value1;
                if (!(_smass)) _smass = value2;
                _T = IF97::T_psmass(_p, _smass);
                if (IF97::BackwardRegion(_p, _smass, IF97_SMASS) == 4) {
                    S = _smass;
                    sVmass = IF97::svap_p(_p);
                    sLmass = IF97::sliq_p(_p);
                    _Q = (S - sLmass) / (sVmass - sLmass);
                    _phase = iphase_twophase;
                } else {
                    _Q = -1;
                    set_phase();
                };
                break;
            case HmolarSmolar_INPUTS:
                // IF97 is mass based so convert smolar input to smass
                _hmass = value1 / molar_mass();  // Convert to mass basis : (J/mol) / (kg/mol) = J/kg
                _smass = value2 / molar_mass();  // Convert to mass basis : (J/mol-K) / (kg/mol) = J/kg-K
                // Fall thru to mass basis inputs
            case HmassSmass_INPUTS:
                _hmass = value1;
                _smass = value2;
                _p = IF97::p_hsmass(_hmass, _smass);
                _T = IF97::T_phmass(_p, _hmass);
                // ...if in the vapor dome (Region 4), calculate Quality...
                if (IF97::BackwardRegion(_p, _hmass, IF97_HMASS) == 4) {
                    H = _hmass;
                    hVmass = IF97::hvap_p(_p);
                    hLmass = IF97::hliq_p(_p);
                    _Q = (H - hLmass) / (hVmass - hLmass);
                    _phase = iphase_twophase;
                } else {
                    _Q = -1;
                    set_phase();
                };
                break;
            default:
                throw ValueError("This pair of inputs is not yet supported");
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
    double calc_SatLiquid(parameters iCalc) {
        switch (iCalc) {
            case iDmass:
                return IF97::rholiq_p(_p);
                break;  ///< Mass-based density
            case iHmass:
                return IF97::hliq_p(_p);
                break;  ///< Mass-based enthalpy
            case iSmass:
                return IF97::sliq_p(_p);
                break;  ///< Mass-based entropy
            case iCpmass:
                return IF97::cpliq_p(_p);
                break;  ///< Mass-based constant-pressure specific heat
            case iCvmass:
                return IF97::cvliq_p(_p);
                break;  ///< Mass-based constant-volume specific heat
            case iUmass:
                return IF97::uliq_p(_p);
                break;  ///< Mass-based internal energy
            case ispeed_sound:
                return IF97::speed_soundliq_p(_p);
                break;  ///< Speed of Sound
            case iviscosity:
                return IF97::viscliq_p(_p);
                break;  ///< Viscosity function
            case iconductivity:
                return IF97::tcondliq_p(_p);
                break;  ///< Thermal conductivity
            case isurface_tension:
                return IF97::sigma97(_T);
                break;  ///< Surface Tension
            case iPrandtl:
                return IF97::prandtlliq_p(_p);
                break;  ///< Prandtl number
            default:
                return -_HUGE;
        };
    }
    double calc_SatVapor(parameters iCalc) {
        switch (iCalc) {
            case iDmass:
                return IF97::rhovap_p(_p);
                break;  ///< Mass-based density
            case iHmass:
                return IF97::hvap_p(_p);
                break;  ///< Mass-based enthalpy
            case iSmass:
                return IF97::svap_p(_p);
                break;  ///< Mass-based entropy
            case iCpmass:
                return IF97::cpvap_p(_p);
                break;  ///< Mass-based constant-pressure specific heat
            case iCvmass:
                return IF97::cvvap_p(_p);
                break;  ///< Mass-based constant-volume specific heat
            case iUmass:
                return IF97::uvap_p(_p);
                break;  ///< Mass-based internal energy
            case ispeed_sound:
                return IF97::speed_soundvap_p(_p);
                break;  ///< Speed of Sound
            case iviscosity:
                return IF97::viscvap_p(_p);
                break;  ///< Viscosity function
            case iconductivity:
                return IF97::tcondvap_p(_p);
                break;  ///< Thermal conductivity
            case isurface_tension:
                return IF97::sigma97(_T);
                break;  ///< Surface Tension
            case iPrandtl:
                return IF97::prandtlvap_p(_p);
                break;  ///< Prandtl number
            default:
                return -_HUGE;
        };
    }
    double calc_Flash(parameters iCalc) {
        switch (_phase) {
            case iphase_twophase:  // In saturation envelope
                if (std::abs(_Q) < 1e-10) {
                    return calc_SatLiquid(iCalc);  // bubble point (Q == 0) on Sat. Liquid curve
                } else if (std::abs(_Q - 1) < 1e-10) {
                    return calc_SatVapor(iCalc);  // dew point (Q == 1) on Sat. Vapor curve
                } else {                          // else "inside" bubble ( 0 < Q < 1 )
                    switch (iCalc) {
                        case iDmass:
                            // Density is an inverse phase weighted property, since it's the inverse of specific volume
                            return 1.0 / (_Q / calc_SatVapor(iDmass) + (1.0 - _Q) / calc_SatLiquid(iDmass));
                            break;
                        case iCpmass:
                            throw NotImplementedError(format("Isobaric Specific Heat not valid in two phase region"));
                            break;
                        case iCvmass:
                            throw NotImplementedError(format("Isochoric Specific Heat not valid in two phase region"));
                            break;
                        case ispeed_sound:
                            throw NotImplementedError(format("Speed of Sound not valid in two phase region"));
                            break;
                        case iviscosity:
                            throw NotImplementedError(format("Viscosity not valid in two phase region"));
                            break;
                        case iconductivity:
                            throw NotImplementedError(format("Viscosity not valid in two phase region"));
                            break;
                        case isurface_tension:
                            return IF97::sigma97(_T);
                            break;  // Surface Tension is not a phase-weighted property
                        case iPrandtl:
                            throw NotImplementedError(format("Prandtl number is not valid in two phase region"));
                            break;
                        default:
                            return _Q * calc_SatVapor(iCalc) + (1.0 - _Q) * calc_SatLiquid(iCalc);  // Phase weighted combination
                    };
                }
                break;
            default:  // Outside saturation envelope (iphase_not_imposed), let IF97 determine phase/region
                switch (iCalc) {
                    case iDmass:
                        return IF97::rhomass_Tp(_T, _p);
                        break;  ///< Mass-based density
                    case iHmass:
                        return IF97::hmass_Tp(_T, _p);
                        break;  ///< Mass-based enthalpy
                    case iSmass:
                        return IF97::smass_Tp(_T, _p);
                        break;  ///< Mass-based entropy
                    case iCpmass:
                        return IF97::cpmass_Tp(_T, _p);
                        break;  ///< Mass-based constant-pressure specific heat
                    case iCvmass:
                        return IF97::cvmass_Tp(_T, _p);
                        break;  ///< Mass-based constant-volume specific heat
                    case iUmass:
                        return IF97::umass_Tp(_T, _p);
                        break;  ///< Mass-based internal energy
                    case ispeed_sound:
                        return IF97::speed_sound_Tp(_T, _p);
                        break;  ///< Speed of sound
                    case iviscosity:
                        return IF97::visc_Tp(_T, _p);
                        break;  ///< Viscosity function
                    case iconductivity:
                        return IF97::tcond_Tp(_T, _p);
                        break;  ///< Thermal conductivity
                    case isurface_tension:
                        throw NotImplementedError(format("Surface Tension is only valid within the two phase region; Try PQ or QT inputs"));
                        break;
                    case iPrandtl:
                        return IF97::prandtl_Tp(_T, _p);
                        break;  ///< Prandtl number
                    default:
                        throw NotImplementedError(format("Output variable not implemented in IF97 Backend"));
                };
        }
    }
    /// Return the mass density in kg/m³
    double rhomass(void) {
        return calc_rhomass();
    };
    double calc_rhomass(void) {
        return calc_Flash(iDmass);
    };
    /// Return the molar density in mol/m³
    double rhomolar(void) {
        return calc_rhomolar();
    };
    double calc_rhomolar(void) {
        return rhomass() / molar_mass();
    };  /// kg/m³ * mol/kg = mol/m³

    /// Return the mass enthalpy in J/kg
    double hmass(void) {
        return calc_hmass();
    };
    double calc_hmass(void) {
        return calc_Flash(iHmass);
    };
    /// Return the molar enthalpy in J/mol
    double hmolar(void) {
        return calc_hmolar();
    };
    double calc_hmolar(void) {
        return hmass() * molar_mass();
    };  /// J/kg * kg/mol = J/mol

    /// Return the mass entropy in J/kg/K
    double smass(void) {
        return calc_smass();
    };
    double calc_smass(void) {
        return calc_Flash(iSmass);
    };
    /// Return the molar entropy in J/mol/K
    double smolar(void) {
        return calc_smolar();
    };
    double calc_smolar(void) {
        return smass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the mass internal energy in J/kg
    double umass(void) {
        return calc_umass();
    };
    double calc_umass(void) {
        return calc_Flash(iUmass);
    };
    /// Return the molar internal energy in J/mol
    double umolar(void) {
        return calc_umolar();
    };
    double calc_umolar(void) {
        return umass() * molar_mass();
    };  /// J/kg * kg/mol = J/mol

    /// Return the mass-based constant pressure specific heat in J/kg/K
    double cpmass(void) {
        return calc_cpmass();
    };
    double calc_cpmass(void) {
        return calc_Flash(iCpmass);
    };
    /// Return the molar-based constant pressure specific heat in J/mol/K
    double cpmolar(void) {
        return calc_cpmolar();
    };
    double calc_cpmolar(void) {
        return cpmass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the mass-based constant volume specific heat in J/kg/K
    double cvmass(void) {
        return calc_cvmass();
    };
    double calc_cvmass(void) {
        return calc_Flash(iCvmass);
    };
    /// Return the molar-based constant volume specific heat in J/mol/K
    double cvmolar(void) {
        return calc_cvmolar();
    };
    double calc_cvmolar(void) {
        return cvmass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the speed of sound in m/s
    double speed_sound(void) {
        return calc_speed_sound();
    };
    double calc_speed_sound(void) {
        return calc_Flash(ispeed_sound);
    };

    // Return the phase
    phases calc_phase(void) {
        return _phase;
    };

    //
    // ************************************************************************* //
    //                         Trivial Functions                                 //
    // ************************************************************************* //
    //
    /// Using this backend, get the triple point temperature in K
    double calc_Ttriple(void) {
        return IF97::get_Ttrip();
    };
    /// Using this backend, get the triple point pressure in Pa
    double calc_p_triple(void) {
        return IF97::get_ptrip();
    };
    /// Using this backend, get the critical point temperature in K
    double calc_T_critical(void) {
        return IF97::get_Tcrit();
    };
    /// Using this backend, get the critical point pressure in Pa
    double calc_p_critical(void) {
        return IF97::get_pcrit();
    };
    /// Using this backend, get the ideal gas constant in J/mol*K
    /// ==> multiplies IF97 Rgas by molar_mass() to put on molar basis per CoolProp convention
    double calc_gas_constant(void) {
        return IF97::get_Rgas() * molar_mass();
    };
    /// Using this backend, get the molar mass in kg/mol
    double calc_molar_mass(void) {
        return IF97::get_MW();
    };
    /// Using this backend, get the acentric factor (unitless)
    double calc_acentric_factor(void) {
        return IF97::get_Acentric();
    };
    /// Using this backend, get the high pressure limit in Pa
    // TODO: May want to adjust this based on _T, since Region 5
    //       is limited to 50 MPa, instead of 100 MPa elsewhere.
    double calc_pmax(void) {
        return IF97::get_Pmax();
    };
    /// Note: Pmin not implemented in Abstract State or CoolProp
    /// Using this backend, get the high temperature limit in K
    double calc_Tmax(void) {
        return IF97::get_Tmax();
    };
    /// Using this backend, get the high pressure limit in K
    double calc_Tmin(void) {
        return IF97::get_Tmin();
    };
    /// Using this backend, get the critical point density in kg/m³
    /// Replace molar-based AbstractState functions since IF97 is mass based only
    double rhomolar_critical(void) {
        return calc_rhomass_critical() / molar_mass();
    }
    double rhomass_critical(void) {
        return calc_rhomass_critical();
    }
    // Overwrite the virtual calc_ functions for density
    double calc_rhomolar_critical(void) {
        return rhomass_critical() / molar_mass();
    };
    double calc_rhomass_critical(void) {
        return IF97::get_rhocrit();
    };
    //
    // ************************************************************************* //
    //                      Saturation Functions                                 //
    // ************************************************************************* //
    //
    double calc_pressure(void) {
        return _p;
    };
    //
    // ************************************************************************* //
    //                 Transport Property Functions                              //
    // ************************************************************************* //
    //
    // Return viscosity in [Pa-s]
    double viscosity(void) {
        return calc_viscosity();
    };
    double calc_viscosity(void) {
        return calc_Flash(iviscosity);
    };
    // Return thermal conductivity in [W/m-K]
    double conductivity(void) {
        return calc_conductivity();
    };
    double calc_conductivity(void) {
        return calc_Flash(iconductivity);
    };
    // Return surface tension in [N/m]
    double surface_tension(void) {
        return calc_surface_tension();
    };
    double calc_surface_tension(void) {
        return calc_Flash(isurface_tension);
    };
    // Return Prandtl number (mu*Cp/k) [dimensionless]
    double Prandtl(void) {
        return calc_Flash(iPrandtl);
    };
};

} /* namespace CoolProp */
#endif /* IF97BACKEND_H_ */
