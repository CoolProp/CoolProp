
#ifndef IF97BACKEND_H_
#define IF97BACKEND_H_

#include "DataStructures.h"
#include <IF97.h>
#include "AbstractState.h"
#include "Exceptions.h"
#include <vector>

namespace CoolProp {

class IF97Backend : public AbstractState
{

   protected:
    /// Additional cached elements used only in this backend since the "normal"
    /// backends use only molar units while IF97 uses mass-based units
    CachedElement _hmass, _rhomass, _smass,
      _reverse;  // Also need a way to flag using the IF97 reverse calcs for h(p,s) and s(p,h)
                 /// CachedElement  _hVmass, _hLmass, _sVmass, sLmass;

   public:
    /// The name of the backend being used
    std::string backend_name() {
        return get_backend_string(IF97_BACKEND);
    }

    /// The fluid this backend represents (always Water by definition).
    /// AbstractState's default calc_fluid_names() throws — override
    /// here so callers like the SVDSBTL adapter (which expects
    /// fluid_names() to work on any source-of-truth backend) don't
    /// trip on it.
    std::vector<std::string> calc_fluid_names() override {
        return {"Water"};
    }

    // REQUIRED BUT NOT USED IN IF97 FUNCTIONS
    bool using_mole_fractions() {
        return false;
    };
    bool using_mass_fractions() {
        return true;
    };  // But actually it doesn't matter since it is only a pure fluid
    bool using_volu_fractions() {
        return false;
    };
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) {
        throw NotImplementedError("Mole composition has not been implemented.");
    };
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) {};  // Not implemented, but don't throw any errors
    void set_volu_fractions(const std::vector<CoolPropDbl>& volu_fractions) {
        throw NotImplementedError("Volume composition has not been implemented.");
    };
    const std::vector<CoolPropDbl>& get_mole_fractions() {
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
        this->_reverse.clear();
        this->_phase = iphase_not_imposed;
        return true;
    };

    // Set phase based on cached values in _p & _T or, if reverse, from (_p,_smass) or (_p,_hmass)
    void set_phase() {
        double epsilon = 3.3e-5;  // IAPWS-IF97 RMS saturated pressure inconsistency
        if (!_reverse) {
            if ((std::abs(_T - IF97::Tcrit) < epsilon / 10.0) &&  //            RMS temperature inconsistency ~ epsilon/10
                (std::abs(_p - IF97::Pcrit) < epsilon)) {         // within epsilon of [Tcrit,Pcrit]
                _phase = iphase_critical_point;                   //     at critical point
            } else if (_T > IF97::Tcrit) {                        // to the right of the critical point
                if (_p > IF97::Pcrit) {                           //     above the critical point pressure
                    _phase = iphase_supercritical;
                } else {  //     below the critical point pressure
                    _phase = iphase_supercritical_gas;
                }
            } else {                     // to the left of the critical point
                if (_p > IF97::Pcrit) {  //     above the critical point pressure
                    _phase = iphase_supercritical_liquid;
                } else {  //     at or below critical point pressure
                    double psat = IF97::psat97(_T);
                    if (_p > psat * (1.0 + epsilon)) {  //         above the saturation curve
                        _phase = iphase_liquid;
                    } else if (_p < psat * (1.0 - epsilon)) {  //         below the saturation curve
                        _phase = iphase_gas;
                    } else  //         exactly on saturation curve (within 1e-4 %)
                        _phase = iphase_twophase;
                }
            }
        } else {  // Backwards: Determine phase from _p & _smass or _hmass
            int IF97Region;
            if (_smass) {                                                   // Get IF97 Region
                IF97Region = IF97::BackwardRegion(_p, _smass, IF97_SMASS);  //     using p, s
            } else {
                IF97Region = IF97::BackwardRegion(_p, _hmass, IF97_HMASS);  //     using p, h
            }
            switch (IF97Region) {  // Convert IF97 Region to CP iPhase
                case 1:            // IF97::REGION1
                    if (_p <= IF97::Pcrit) {
                        _phase = iphase_liquid;
                    } else {
                        _phase = iphase_supercritical_liquid;
                    };
                    break;
                case 2:  // IF97::REGION2
                    if (_T <= IF97::Tcrit) {
                        _phase = iphase_gas;
                    } else {
                        _phase = iphase_supercritical_gas;
                    };
                    break;
                case 3:  // IF97::REGION3
                    if (_T < IF97::Tsat97(_p)) {
                        if (_p <= IF97::Pcrit) {
                            _phase = iphase_liquid;
                        } else {
                            _phase = iphase_supercritical_liquid;
                        };
                    } else {
                        if (_T <= IF97::Tcrit) {
                            _phase = iphase_gas;
                        } else {
                            _phase = iphase_supercritical_gas;
                        }
                    };
                    break;
                case 4:                        // IF97::REGION4 (Saturation
                    _phase = iphase_twophase;  //   already handled but here for completeness
                    break;
                case 5:  // IF97::REGION5 or O.B.
                default:
                    throw CoolProp::OutOfRangeError("Outside of IF97 Reverse Function Bounds");
                    break;
            };
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
                _reverse = 1.0;  // only ever set for HP and SP Inputs
                _T = IF97::T_phmass(_p, _hmass);
                _Q = -1;  // Default if not set in 2-Phase Region
                // ...if in the vapor dome (Region 4), calculate Quality...
                if (IF97::BackwardRegion(_p, _hmass, IF97_HMASS) == 4) {
                    H = _hmass;
                    hVmass = IF97::hvap_p(_p);
                    hLmass = IF97::hliq_p(_p);
                    _Q = std::min(1.0, std::max(0.0, (H - hLmass) / (hVmass - hLmass)));  //bound between 0 and 1
                    _phase = iphase_twophase;
                } else {
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
                _reverse = 1.0;
                _T = IF97::T_psmass(_p, _smass);
                _Q = -1;
                // ...if in the vapor dome (Region 4), calculate Quality...
                if (IF97::BackwardRegion(_p, _smass, IF97_SMASS) == 4) {
                    S = _smass;
                    sVmass = IF97::svap_p(_p);
                    sLmass = IF97::sliq_p(_p);
                    _Q = std::min(1.0, std::max(0.0, (S - sLmass) / (sVmass - sLmass)));  //bound between 0 and 1
                    _phase = iphase_twophase;
                } else {
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
                    _Q = std::min(1.0, std::max(0.0, (H - hLmass) / (hVmass - hLmass)));  //bount between 0 and 1
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
    double rhomass() {
        return calc_rhomass();
    };
    double calc_rhomass() {
        return calc_Flash(iDmass);
    };
    /// Return the molar density in mol/m³
    double rhomolar() {
        return calc_rhomolar();
    };
    double calc_rhomolar() {
        return rhomass() / molar_mass();
    };  /// kg/m³ * mol/kg = mol/m³

    /// Return the mass enthalpy in J/kg
    double hmass() {
        return calc_hmass();
    };
    double calc_hmass() {
        if (_reverse && _smass)
            return IF97::hmass_psmass(_p, _smass);  // Special IF97 function for handling reverse h(p,s) evaluation
        else
            return calc_Flash(iHmass);
    };
    /// Return the molar enthalpy in J/mol
    double hmolar() {
        return calc_hmolar();
    };
    double calc_hmolar() {
        return hmass() * molar_mass();
    };  /// J/kg * kg/mol = J/mol

    /// Return the mass entropy in J/kg/K
    double smass() {
        return calc_smass();
    };
    double calc_smass() {
        if (_reverse && _hmass)
            return IF97::smass_phmass(_p, _hmass);  // Special IF97 function for handling reverse s(p,h) evaluation
        else
            return calc_Flash(iSmass);
    };
    /// Return the molar entropy in J/mol/K
    double smolar() {
        return calc_smolar();
    };
    double calc_smolar() {
        return smass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the mass internal energy in J/kg
    double umass() {
        return calc_umass();
    };
    double calc_umass() {
        return calc_Flash(iUmass);
    };
    /// Return the molar internal energy in J/mol
    double umolar() {
        return calc_umolar();
    };
    double calc_umolar() {
        return umass() * molar_mass();
    };  /// J/kg * kg/mol = J/mol

    /// Return the mass-based constant pressure specific heat in J/kg/K
    double cpmass() {
        return calc_cpmass();
    };
    double calc_cpmass() {
        return calc_Flash(iCpmass);
    };
    /// Return the molar-based constant pressure specific heat in J/mol/K
    double cpmolar() {
        return calc_cpmolar();
    };
    double calc_cpmolar() {
        return cpmass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the mass-based constant volume specific heat in J/kg/K
    double cvmass() {
        return calc_cvmass();
    };
    double calc_cvmass() {
        return calc_Flash(iCvmass);
    };
    /// Return the molar-based constant volume specific heat in J/mol/K
    double cvmolar() {
        return calc_cvmolar();
    };
    double calc_cvmolar() {
        return cvmass() * molar_mass();
    };  /// J/kg-K * kg/mol = J/mol-K

    /// Return the speed of sound in m/s
    double speed_sound() {
        return calc_speed_sound();
    };
    double calc_speed_sound() {
        return calc_Flash(ispeed_sound);
    };

    // Return the phase
    phases calc_phase() {
        return _phase;
    };

    //
    // ************************************************************************* //
    //                         Trivial Functions                                 //
    // ************************************************************************* //
    //
    /// Using this backend, get the triple point temperature in K
    double calc_Ttriple() {
        return IF97::get_Ttrip();
    };
    /// Using this backend, get the triple point pressure in Pa
    double calc_p_triple() {
        return IF97::get_ptrip();
    };
    /// Using this backend, get the critical point temperature in K
    double calc_T_critical() {
        return IF97::get_Tcrit();
    };
    /// Using this backend, get the critical point pressure in Pa
    double calc_p_critical() {
        return IF97::get_pcrit();
    };
    /// Using this backend, get the ideal gas constant in J/mol*K
    /// ==> multiplies IF97 Rgas by molar_mass() to put on molar basis per CoolProp convention
    double calc_gas_constant() {
        return IF97::get_Rgas() * molar_mass();
    };
    /// Using this backend, get the molar mass in kg/mol
    double calc_molar_mass() {
        return IF97::get_MW();
    };
    /// Using this backend, get the acentric factor (unitless)
    double calc_acentric_factor() {
        return IF97::get_Acentric();
    };
    /// Using this backend, get the high pressure limit in Pa
    // TODO: May want to adjust this based on _T, since Region 5
    //       is limited to 50 MPa, instead of 100 MPa elsewhere.
    double calc_pmax() {
        return IF97::get_Pmax();
    };
    /// Note: Pmin not implemented in Abstract State or CoolProp
    /// Using this backend, get the high temperature limit in K
    double calc_Tmax() {
        return IF97::get_Tmax();
    };
    /// Using this backend, get the high pressure limit in K
    double calc_Tmin() {
        return IF97::get_Tmin();
    };
    /// Using this backend, get the critical point density in kg/m³
    /// Replace molar-based AbstractState functions since IF97 is mass based only
    double rhomolar_critical() {
        return calc_rhomass_critical() / molar_mass();
    }
    double rhomass_critical() {
        return calc_rhomass_critical();
    }
    // Overwrite the virtual calc_ functions for density
    double calc_rhomolar_critical() {
        return rhomass_critical() / molar_mass();
    };
    double calc_rhomass_critical() {
        return IF97::get_rhocrit();
    };
    //
    // ************************************************************************* //
    //                      Saturation Functions                                 //
    // ************************************************************************* //
    //
    double calc_pressure() {
        return _p;
    };
    //
    // ************************************************************************* //
    //                 Transport Property Functions                              //
    // ************************************************************************* //
    //
    // Return viscosity in [Pa-s]
    double viscosity() {
        return calc_viscosity();
    };
    double calc_viscosity() {
        return calc_Flash(iviscosity);
    };
    // Return thermal conductivity in [W/m-K]
    double conductivity() {
        return calc_conductivity();
    };
    double calc_conductivity() {
        return calc_Flash(iconductivity);
    };
    // Return surface tension in [N/m]
    double surface_tension() {
        return calc_surface_tension();
    };
    double calc_surface_tension() {
        return calc_Flash(isurface_tension);
    };
    // Return Prandtl number (mu*Cp/k) [dimensionless]
    double Prandtl() {
        return calc_Flash(iPrandtl);
    };

    /// Vectorized direct evaluation. See AbstractState::fast_evaluate for contract.
    /// IF97 is stateless under the hood; this entry point just dispatches each
    /// (T,p) or (h,p) point straight to the IF97 namespace evaluators with no
    /// allocations and no writes to the AbstractState cache.
    void fast_evaluate(CoolProp::input_pairs input_pair, const double* val1, const double* val2, std::size_t N_inputs,
                       const CoolProp::parameters* outputs, std::size_t N_outputs, double* out_buffer, std::size_t out_buffer_size, int* status_flags,
                       std::size_t status_flags_size, CoolProp::phases imposed_phase = CoolProp::iphase_not_imposed) override {
        if (N_outputs != 0 && N_inputs > std::numeric_limits<std::size_t>::max() / N_outputs) {
            throw ValueError(format("fast_evaluate: N_inputs * N_outputs would overflow size_t (N_inputs=%zu, N_outputs=%zu)", N_inputs, N_outputs));
        }
        const std::size_t required_out = N_inputs * N_outputs;
        if (out_buffer_size < required_out) {
            throw ValueError(format("fast_evaluate: out_buffer_size=%zu < required %zu (N_inputs * N_outputs)", out_buffer_size, required_out));
        }
        if (status_flags_size < N_inputs) {
            throw ValueError(format("fast_evaluate: status_flags_size=%zu < required %zu (N_inputs)", status_flags_size, N_inputs));
        }
        if (N_inputs == 0) return;
        // Null-pointer check BEFORE the N_outputs==0 early-write path: that path
        // dereferences status_flags, so the check has to dominate it.
        if (val1 == nullptr || val2 == nullptr || outputs == nullptr || out_buffer == nullptr || status_flags == nullptr) {
            throw ValueError("fast_evaluate: null pointer argument");
        }
        if (N_outputs == 0) {
            for (std::size_t k = 0; k < N_inputs; ++k)
                status_flags[k] = CoolProp::fast_evaluate_ok;
            return;
        }

        const bool is_PT = (input_pair == PT_INPUTS);
        const bool is_HmassP = (input_pair == HmassP_INPUTS);
        const bool is_HmolarP = (input_pair == HmolarP_INPUTS);
        if (!is_PT && !is_HmassP && !is_HmolarP) {
            throw ValueError(format("fast_evaluate (IF97): input_pair %s not supported (use PT_INPUTS, HmassP_INPUTS, HmolarP_INPUTS)",
                                    get_input_pair_short_desc(input_pair).c_str()));
        }

        // Validate output keys up-front. IF97 covers everything its calc_Flash supports.
        constexpr std::size_t MAX_FE_OUTPUTS = 64;
        if (N_outputs > MAX_FE_OUTPUTS) {
            throw ValueError(format("fast_evaluate: N_outputs=%zu exceeds compile-time limit %zu", N_outputs, MAX_FE_OUTPUTS));
        }
        for (std::size_t o = 0; o < N_outputs; ++o) {
            switch (outputs[o]) {
                case iT:
                case iP:
                case iDmass:
                case iDmolar:
                case iHmass:
                case iHmolar:
                case iSmass:
                case iSmolar:
                case iUmass:
                case iUmolar:
                case iCpmass:
                case iCpmolar:
                case iCvmass:
                case iCvmolar:
                case ispeed_sound:
                case iviscosity:
                case iconductivity:
                case iPrandtl:
                    break;
                default:
                    throw ValueError(
                      format("fast_evaluate (IF97): output[%zu]=%s not supported", o, get_parameter_information(outputs[o], "short").c_str()));
            }
        }

        const double M = calc_molar_mass();  // kg/mol, constant for water
        const double NaN = std::numeric_limits<double>::quiet_NaN();

        // imposed_phase: IF97 derives phase from (T,p) — we honour the hint by
        // skipping the two-phase check when caller asserts single-phase.
        const bool skip_twophase_check = (imposed_phase != iphase_not_imposed && imposed_phase != iphase_twophase && imposed_phase != iphase_unknown);

        auto fill_nan_row = [&](std::size_t k) {
            for (std::size_t o = 0; o < N_outputs; ++o) {
                out_buffer[k * N_outputs + o] = NaN;
            }
        };

        for (std::size_t k = 0; k < N_inputs; ++k) {
            const double v1 = val1[k];
            const double v2 = val2[k];

            double T_k, p_k;
            if (is_PT) {
                p_k = v1;
                T_k = v2;
            } else {
                // (h, p) input — convert h to mass basis and use IF97's
                // backward T_phmass to recover T. Both branches must produce
                // a clean (T, p) for the forward evaluators below.
                p_k = v2;
                double hmass_k = is_HmolarP ? (v1 / M) : v1;
                try {
                    T_k = IF97::T_phmass(p_k, hmass_k);
                } catch (const std::exception&) {
                    status_flags[k] = CoolProp::fast_evaluate_out_of_range;
                    fill_nan_row(k);
                    continue;
                }
            }

            // Range check before evaluation. IF97 has explicit T/p limits; we
            // reject points outside and let in-range points go to the kernel,
            // which itself will throw on Region-5 etc. (caught below).
            if (!(p_k > 0 && T_k > 0 && p_k <= IF97::get_Pmax() && T_k >= IF97::get_Tmin() && T_k <= IF97::get_Tmax())) {
                status_flags[k] = CoolProp::fast_evaluate_out_of_range;
                fill_nan_row(k);
                continue;
            }

            // Two-phase rejection: PT_INPUTS landing exactly on the saturation
            // curve is ambiguous; caller can use imposed_phase or QT/PQ inputs
            // (which we don't support here yet). For (h,p) inputs the backward
            // solver places us in a specific region so this is less common.
            if (!skip_twophase_check && is_PT && T_k <= IF97::Tcrit) {
                const double psat = IF97::psat97(T_k);
                const double eps = 3.3e-5;
                if (std::abs(p_k - psat) <= psat * eps) {
                    status_flags[k] = CoolProp::fast_evaluate_two_phase_disallowed;
                    fill_nan_row(k);
                    continue;
                }
            }

            bool eval_failed = false;
            for (std::size_t o = 0; o < N_outputs; ++o) {
                const parameters out_key = outputs[o];
                try {
                    double val;
                    switch (out_key) {
                        case iT:
                            val = T_k;
                            break;
                        case iP:
                            val = p_k;
                            break;
                        case iDmass:
                            val = IF97::rhomass_Tp(T_k, p_k);
                            break;
                        case iDmolar:
                            val = IF97::rhomass_Tp(T_k, p_k) / M;
                            break;
                        case iHmass:
                            val = IF97::hmass_Tp(T_k, p_k);
                            break;
                        case iHmolar:
                            val = IF97::hmass_Tp(T_k, p_k) * M;
                            break;
                        case iSmass:
                            val = IF97::smass_Tp(T_k, p_k);
                            break;
                        case iSmolar:
                            val = IF97::smass_Tp(T_k, p_k) * M;
                            break;
                        case iUmass:
                            val = IF97::umass_Tp(T_k, p_k);
                            break;
                        case iUmolar:
                            val = IF97::umass_Tp(T_k, p_k) * M;
                            break;
                        case iCpmass:
                            val = IF97::cpmass_Tp(T_k, p_k);
                            break;
                        case iCpmolar:
                            val = IF97::cpmass_Tp(T_k, p_k) * M;
                            break;
                        case iCvmass:
                            val = IF97::cvmass_Tp(T_k, p_k);
                            break;
                        case iCvmolar:
                            val = IF97::cvmass_Tp(T_k, p_k) * M;
                            break;
                        case ispeed_sound:
                            val = IF97::speed_sound_Tp(T_k, p_k);
                            break;
                        case iviscosity:
                            val = IF97::visc_Tp(T_k, p_k);
                            break;
                        case iconductivity:
                            val = IF97::tcond_Tp(T_k, p_k);
                            break;
                        case iPrandtl:
                            val = IF97::prandtl_Tp(T_k, p_k);
                            break;
                        default:
                            val = NaN;
                            eval_failed = true;
                            break;
                    }
                    out_buffer[k * N_outputs + o] = val;
                } catch (const std::exception&) {
                    eval_failed = true;
                    out_buffer[k * N_outputs + o] = NaN;
                }
            }
            // API contract: when status_flags[k] is non-zero, the entire row is
            // NaN. The per-output catch above only NaNs the specific failing
            // output; finish the job here so callers don't see partial rows.
            if (eval_failed) {
                fill_nan_row(k);
                status_flags[k] = CoolProp::fast_evaluate_internal_error;
            } else {
                status_flags[k] = CoolProp::fast_evaluate_ok;
            }
        }
        // IF97 is stateless externally; clear() resets internal cached values to
        // be consistent with the TabularBackend contract (state untouched).
        clear();
    };
};

} /* namespace CoolProp */
#endif /* IF97BACKEND_H_ */
