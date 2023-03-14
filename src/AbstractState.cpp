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
#include "DataStructures.h"
#include "Backends/IF97/IF97Backend.h"
#include "Backends/Cubics/CubicBackend.h"
#include "Backends/Cubics/VTPRBackend.h"
#include "Backends/Incompressible/IncompressibleBackend.h"
#include "Backends/PCSAFT/PCSAFTBackend.h"

#if !defined(NO_TABULAR_BACKENDS)
#include "Backends/Tabular/TTSEBackend.h"
#include "Backends/Tabular/BicubicBackend.h"
#endif

namespace CoolProp {

/// This tiny class holds pointers to generators for the backends and can be used to look up
/// generators at runtime.  This class should be populated through the use of static initialized

class BackendLibrary
{
   private:
    std::map<backend_families, shared_ptr<AbstractStateGenerator>> backends;

   public:
    void add_backend(const backend_families& bg, const shared_ptr<AbstractStateGenerator>& asg) {
        backends[bg] = asg;
    };
    void get_generator_iterators(const backend_families& bg,
                                 std::map<backend_families, shared_ptr<AbstractStateGenerator>>::const_iterator& generator,
                                 std::map<backend_families, shared_ptr<AbstractStateGenerator>>::const_iterator& end) {
        generator = backends.find(bg);
        end = backends.end();
    };
    std::size_t size() {
        return backends.size();
    };
};
inline BackendLibrary& get_backend_library() {
    static BackendLibrary the_library;
    return the_library;
}

void register_backend(const backend_families& bf, shared_ptr<AbstractStateGenerator> gen) {
    get_backend_library().add_backend(bf, gen);
};

class IF97BackendGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        if (fluid_names.size() == 1) {         // Check that fluid_names[0] has only one component
            std::string str = fluid_names[0];  // Check that the fluid name is an alias for "Water"
            if ((upper(str) == "WATER") || (upper(str) == "H2O")) {
                return new IF97Backend();
            } else {
                throw ValueError(format("The IF97 backend returns Water props only; fluid name [%s] not allowed", fluid_names[0].c_str()));
            }
        } else {
            throw ValueError(format("The IF97 backend does not support mixtures, only Water"));
        };
    };
};
// This static initialization will cause the generator to register
static GeneratorInitializer<IF97BackendGenerator> if97_gen(IF97_BACKEND_FAMILY);
class SRKGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        return new SRKBackend(fluid_names, get_config_double(R_U_CODATA));
    };
};
static GeneratorInitializer<SRKGenerator> srk_gen(CoolProp::SRK_BACKEND_FAMILY);
class PRGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        return new PengRobinsonBackend(fluid_names, get_config_double(R_U_CODATA));
    };
};
static GeneratorInitializer<PRGenerator> pr_gen(CoolProp::PR_BACKEND_FAMILY);
class IncompressibleBackendGenerator : public AbstractStateGenerator
{
   public:
    AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        if (fluid_names.size() != 1) {
            throw ValueError(format("For INCOMP backend, name vector must be one element long"));
        }
        return new IncompressibleBackend(fluid_names[0]);
    };
};
// This static initialization will cause the generator to register
static GeneratorInitializer<IncompressibleBackendGenerator> incomp_gen(INCOMP_BACKEND_FAMILY);
class VTPRGenerator : public CoolProp::AbstractStateGenerator
{
   public:
    CoolProp::AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        return new CoolProp::VTPRBackend(fluid_names, CoolProp::get_config_double(R_U_CODATA));
    };
};
// This static initialization will cause the generator to register
static CoolProp::GeneratorInitializer<VTPRGenerator> vtpr_gen(CoolProp::VTPR_BACKEND_FAMILY);

class PCSAFTGenerator : public CoolProp::AbstractStateGenerator
{
   public:
    CoolProp::AbstractState* get_AbstractState(const std::vector<std::string>& fluid_names) {
        return new CoolProp::PCSAFTBackend(fluid_names);
    };
};
// This static initialization will cause the generator to register
static CoolProp::GeneratorInitializer<PCSAFTGenerator> pcsaft_gen(CoolProp::PCSAFT_BACKEND_FAMILY);

AbstractState* AbstractState::factory(const std::string& backend, const std::vector<std::string>& fluid_names) {
    if (get_debug_level() > 0) {
        std::cout << "AbstractState::factory(" << backend << "," << stringvec_to_string(fluid_names) << ")" << std::endl;
    }

    backend_families f1;
    std::string f2;
    extract_backend_families_string(backend, f1, f2);

    std::map<backend_families, shared_ptr<AbstractStateGenerator>>::const_iterator gen, end;
    get_backend_library().get_generator_iterators(f1, gen, end);

    if (get_debug_level() > 0) {
        std::cout << "AbstractState::factory backend_library size: " << get_backend_library().size() << std::endl;
    }

    if (gen != end) {
        // One of the registered backends was able to match the given backend family
        return gen->second->get_AbstractState(fluid_names);
    }
#if !defined(NO_TABULAR_BACKENDS)
    else if (f1 == TTSE_BACKEND_FAMILY) {
        // Will throw if there is a problem with this backend
        shared_ptr<AbstractState> AS(factory(f2, fluid_names));
        return new TTSEBackend(AS);
    } else if (f1 == BICUBIC_BACKEND_FAMILY) {
        // Will throw if there is a problem with this backend
        shared_ptr<AbstractState> AS(factory(f2, fluid_names));
        return new BicubicBackend(AS);
    }
#endif
    else if (!backend.compare("?") || backend.empty()) {
        std::size_t idel = fluid_names[0].find("::");
        // Backend has not been specified, and we have to figure out what the backend is by parsing the string
        if (idel == std::string::npos)  // No '::' found, no backend specified, try HEOS, otherwise a failure
        {
            // Figure out what backend to use
            return factory("HEOS", fluid_names);
        } else {
            // Split string at the '::' into two std::string, call again
            return factory(std::string(fluid_names[0].begin(), fluid_names[0].begin() + idel),
                           std::string(fluid_names[0].begin() + (idel + 2), fluid_names[0].end()));
        }
    } else {
        throw ValueError(format("Invalid backend name [%s] to factory function", backend.c_str()));
    }
}
std::vector<std::string> AbstractState::fluid_names(void) {
    return calc_fluid_names();
}
bool AbstractState::clear_comp_change() {
    // Reset all instances of CachedElement and overwrite
    // the internal double values with -_HUGE
    this->_R = _HUGE;
    this->_gas_constant.clear();
    this->_molar_mass.clear();
    this->_critical.fill(_HUGE);
    this->_reducing.fill(_HUGE);
    return true;
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
    this->_helmholtzmolar.clear();
    this->_logp.clear();
    this->_logrhomolar.clear();

    this->_hmolar_excess.clear();
    this->_smolar_excess.clear();
    this->_gibbsmolar_excess.clear();
    this->_volumemolar_excess.clear();
    this->_umolar_excess.clear();
    this->_helmholtzmolar_excess.clear();

    this->_hmolar_residual.clear();
    this->_smolar_residual.clear();
    this->_gibbsmolar_residual.clear();

    /// Smoothing values
    this->_rho_spline.clear();
    this->_drho_spline_dh__constp.clear();
    this->_drho_spline_dp__consth.clear();

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
void AbstractState::mass_to_molar_inputs(CoolProp::input_pairs& input_pair, CoolPropDbl& value1, CoolPropDbl& value2) {
    // Check if a mass based input, convert it to molar units

    switch (input_pair) {
        case DmassT_INPUTS:  ///< Mass density in kg/m^3, Temperature in K
          //case HmassT_INPUTS: ///< Enthalpy in J/kg, Temperature in K (NOT CURRENTLY IMPLEMENTED)
        case SmassT_INPUTS:  ///< Entropy in J/kg/K, Temperature in K
          //case TUmass_INPUTS: ///< Temperature in K, Internal energy in J/kg (NOT CURRENTLY IMPLEMENTED)
        case DmassP_INPUTS:      ///< Mass density in kg/m^3, Pressure in Pa
        case DmassQ_INPUTS:      ///< Mass density in kg/m^3, molar quality
        case HmassP_INPUTS:      ///< Enthalpy in J/kg, Pressure in Pa
        case PSmass_INPUTS:      ///< Pressure in Pa, Entropy in J/kg/K
        case PUmass_INPUTS:      ///< Pressure in Pa, Internal energy in J/kg
        case HmassSmass_INPUTS:  ///< Enthalpy in J/kg, Entropy in J/kg/K
        case SmassUmass_INPUTS:  ///< Entropy in J/kg/K, Internal energy in J/kg
        case DmassHmass_INPUTS:  ///< Mass density in kg/m^3, Enthalpy in J/kg
        case DmassSmass_INPUTS:  ///< Mass density in kg/m^3, Entropy in J/kg/K
        case DmassUmass_INPUTS:  ///< Mass density in kg/m^3, Internal energy in J/kg
        {
            // Set the cache value for the molar mass if it hasn't been set yet
            molar_mass();

            // Molar mass (just for compactness of the following switch)
            CoolPropDbl mm = static_cast<CoolPropDbl>(_molar_mass);

            switch (input_pair) {
                case DmassT_INPUTS:
                    input_pair = DmolarT_INPUTS;
                    value1 /= mm;
                    break;
                    //case HmassT_INPUTS: input_pair = HmolarT_INPUTS; value1 *= mm;  break; (NOT CURRENTLY IMPLEMENTED)
                case SmassT_INPUTS:
                    input_pair = SmolarT_INPUTS;
                    value1 *= mm;
                    break;
                    //case TUmass_INPUTS: input_pair = TUmolar_INPUTS; value2 *= mm;  break; (NOT CURRENTLY IMPLEMENTED)
                case DmassP_INPUTS:
                    input_pair = DmolarP_INPUTS;
                    value1 /= mm;
                    break;
                case DmassQ_INPUTS:
                    input_pair = DmolarQ_INPUTS;
                    value1 /= mm;
                    break;
                case HmassP_INPUTS:
                    input_pair = HmolarP_INPUTS;
                    value1 *= mm;
                    break;
                case PSmass_INPUTS:
                    input_pair = PSmolar_INPUTS;
                    value2 *= mm;
                    break;
                case PUmass_INPUTS:
                    input_pair = PUmolar_INPUTS;
                    value2 *= mm;
                    break;
                case HmassSmass_INPUTS:
                    input_pair = HmolarSmolar_INPUTS;
                    value1 *= mm;
                    value2 *= mm;
                    break;
                case SmassUmass_INPUTS:
                    input_pair = SmolarUmolar_INPUTS;
                    value1 *= mm;
                    value2 *= mm;
                    break;
                case DmassHmass_INPUTS:
                    input_pair = DmolarHmolar_INPUTS;
                    value1 /= mm;
                    value2 *= mm;
                    break;
                case DmassSmass_INPUTS:
                    input_pair = DmolarSmolar_INPUTS;
                    value1 /= mm;
                    value2 *= mm;
                    break;
                case DmassUmass_INPUTS:
                    input_pair = DmolarUmolar_INPUTS;
                    value1 /= mm;
                    value2 *= mm;
                    break;
                default:
                    break;
            }
            break;
        }
        default:
            return;
    }
}
double AbstractState::trivial_keyed_output(parameters key) {
    if (get_debug_level() >= 50)
        std::cout << format("AbstractState: trivial_keyed_output called for %s ", get_parameter_information(key, "short").c_str()) << std::endl;
    switch (key) {
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
            throw ValueError(
              format("This input [%d: \"%s\"] is not valid for trivial_keyed_output", key, get_parameter_information(key, "short").c_str()));
    }
}
double AbstractState::keyed_output(parameters key) {
    if (get_debug_level() >= 50)
        std::cout << format("AbstractState: keyed_output called for %s ", get_parameter_information(key, "short").c_str()) << std::endl;
    // Handle trivial inputs
    if (is_trivial_parameter(key)) {
        return trivial_keyed_output(key);
    }
    switch (key) {
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
        case iHmolar_residual:
            return hmolar_residual();
        case iHmass:
            return hmass();
        case iSmolar:
            return smolar();
        case iSmolar_residual:
            return smolar_residual();
        case iSmass:
            return smass();
        case iUmolar:
            return umolar();
        case iUmass:
            return umass();
        case iGmolar:
            return gibbsmolar();
        case iGmolar_residual:
            return gibbsmolar_residual();
        case iGmass:
            return gibbsmass();
        case iHelmholtzmolar:
            return helmholtzmolar();
        case iHelmholtzmass:
            return helmholtzmass();
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
        case id2alpha0_ddelta2_consttau:
            return d2alpha0_dDelta2();
        case id3alpha0_ddelta3_consttau:
            return d3alpha0_dDelta3();
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
        case iisentropic_expansion_coefficient:
            return isentropic_expansion_coefficient();
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
        case ifundamental_derivative_of_gas_dynamics:
            return fundamental_derivative_of_gas_dynamics();
        default:
            throw ValueError(format("This input [%d: \"%s\"] is not valid for keyed_output", key, get_parameter_information(key, "short").c_str()));
    }
}

double AbstractState::tau(void) {
    if (!_tau) _tau = calc_reciprocal_reduced_temperature();
    return _tau;
}
double AbstractState::delta(void) {
    if (!_delta) _delta = calc_reduced_density();
    return _delta;
}
double AbstractState::Tmin(void) {
    return calc_Tmin();
}
double AbstractState::Tmax(void) {
    return calc_Tmax();
}
double AbstractState::Ttriple(void) {
    return calc_Ttriple();
}
double AbstractState::pmax(void) {
    return calc_pmax();
}
double AbstractState::T_critical(void) {
    return calc_T_critical();
}
double AbstractState::T_reducing(void) {
    if (!ValidNumber(_reducing.T)) {
        calc_reducing_state();
    }
    return _reducing.T;
}
double AbstractState::p_critical(void) {
    return calc_p_critical();
}
double AbstractState::p_triple(void) {
    return calc_p_triple();
}
double AbstractState::rhomolar_critical(void) {
    return calc_rhomolar_critical();
}
double AbstractState::rhomass_critical(void) {
    return calc_rhomolar_critical() * molar_mass();
}
double AbstractState::rhomolar_reducing(void) {
    if (!ValidNumber(_reducing.rhomolar)) {
        calc_reducing_state();
    }
    return _reducing.rhomolar;
}
double AbstractState::rhomass_reducing(void) {
    return rhomolar_reducing() * molar_mass();
}
double AbstractState::hmolar(void) {
    if (!_hmolar) _hmolar = calc_hmolar();
    return _hmolar;
}
double AbstractState::hmolar_residual(void) {
    if (!_hmolar_residual) _hmolar_residual = calc_hmolar_residual();
    return _hmolar_residual;
}
double AbstractState::hmolar_excess(void) {
    if (!_hmolar_excess) calc_excess_properties();
    return _hmolar_excess;
}
double AbstractState::smolar(void) {
    if (!_smolar) _smolar = calc_smolar();
    return _smolar;
}
double AbstractState::smolar_residual(void) {
    if (!_smolar_residual) _smolar_residual = calc_smolar_residual();
    return _smolar_residual;
}
double AbstractState::smolar_excess(void) {
    if (!_smolar_excess) calc_excess_properties();
    return _smolar_excess;
}
double AbstractState::umolar(void) {
    if (!_umolar) _umolar = calc_umolar();
    return _umolar;
}
double AbstractState::umolar_excess(void) {
    if (!_umolar_excess) calc_excess_properties();
    return _umolar_excess;
}
double AbstractState::gibbsmolar(void) {
    if (!_gibbsmolar) _gibbsmolar = calc_gibbsmolar();
    return _gibbsmolar;
}
double AbstractState::gibbsmolar_residual(void) {
    if (!_gibbsmolar_residual) _gibbsmolar_residual = calc_gibbsmolar_residual();
    return _gibbsmolar_residual;
}
double AbstractState::gibbsmolar_excess(void) {
    if (!_gibbsmolar_excess) calc_excess_properties();
    return _gibbsmolar_excess;
}
double AbstractState::helmholtzmolar(void) {
    if (!_helmholtzmolar) _helmholtzmolar = calc_helmholtzmolar();
    return _helmholtzmolar;
}
double AbstractState::helmholtzmolar_excess(void) {
    if (!_helmholtzmolar_excess) calc_excess_properties();
    return _helmholtzmolar_excess;
}
double AbstractState::volumemolar_excess(void) {
    if (!_volumemolar_excess) calc_excess_properties();
    return _volumemolar_excess;
}
double AbstractState::cpmolar(void) {
    if (!_cpmolar) _cpmolar = calc_cpmolar();
    return _cpmolar;
}
double AbstractState::cp0molar(void) {
    return calc_cpmolar_idealgas();
}
double AbstractState::cvmolar(void) {
    if (!_cvmolar) _cvmolar = calc_cvmolar();
    return _cvmolar;
}
double AbstractState::speed_sound(void) {
    if (!_speed_sound) _speed_sound = calc_speed_sound();
    return _speed_sound;
}
double AbstractState::viscosity(void) {
    if (!_viscosity) _viscosity = calc_viscosity();
    return _viscosity;
}
double AbstractState::conductivity(void) {
    if (!_conductivity) _conductivity = calc_conductivity();
    return _conductivity;
}
double AbstractState::melting_line(int param, int given, double value) {
    return calc_melting_line(param, given, value);
}
double AbstractState::acentric_factor() {
    return calc_acentric_factor();
}
double AbstractState::saturation_ancillary(parameters param, int Q, parameters given, double value) {
    return calc_saturation_ancillary(param, Q, given, value);
}
double AbstractState::surface_tension(void) {
    if (!_surface_tension) _surface_tension = calc_surface_tension();
    return _surface_tension;
}
double AbstractState::molar_mass(void) {
    if (!_molar_mass) _molar_mass = calc_molar_mass();
    return _molar_mass;
}
double AbstractState::gas_constant(void) {
    if (!_gas_constant) _gas_constant = calc_gas_constant();
    return _gas_constant;
}
double AbstractState::fugacity_coefficient(std::size_t i) {
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity_coefficient(i);
}
std::vector<double> AbstractState::fugacity_coefficients() {
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity_coefficients();
}
double AbstractState::fugacity(std::size_t i) {
    // TODO: Cache the fug. coeff for each component
    return calc_fugacity(i);
}
double AbstractState::chemical_potential(std::size_t i) {
    // TODO: Cache the chemical potential for each component
    return calc_chemical_potential(i);
}
void AbstractState::build_phase_envelope(const std::string& type) {
    calc_phase_envelope(type);
}
double AbstractState::isothermal_compressibility(void) {
    return 1.0 / _rhomolar * first_partial_deriv(iDmolar, iP, iT);
}
double AbstractState::isobaric_expansion_coefficient(void) {
    return -1.0 / _rhomolar * first_partial_deriv(iDmolar, iT, iP);
}
double AbstractState::isentropic_expansion_coefficient(void) {
    return _rhomolar / _p * first_partial_deriv(iP, iDmolar, iSmolar);
}
double AbstractState::Bvirial(void) {
    return calc_Bvirial();
}
double AbstractState::Cvirial(void) {
    return calc_Cvirial();
}
double AbstractState::dBvirial_dT(void) {
    return calc_dBvirial_dT();
}
double AbstractState::dCvirial_dT(void) {
    return calc_dCvirial_dT();
}
double AbstractState::compressibility_factor(void) {
    return calc_compressibility_factor();
}

double AbstractState::fundamental_derivative_of_gas_dynamics() {
    // See Colonna, FPE, 2010, Eq. 1
    return 1 + this->second_partial_deriv(iP, iDmass, iSmolar, iDmass, iSmolar) * this->rhomass() / (2 * powInt(speed_sound(), 2));
};

// Get the derivatives of the parameters in the partial derivative with respect to T and rho
void get_dT_drho(AbstractState& AS, parameters index, CoolPropDbl& dT, CoolPropDbl& drho) {
    CoolPropDbl T = AS.T(), rho = AS.rhomolar(), rhor = AS.rhomolar_reducing(), Tr = AS.T_reducing(), dT_dtau = -pow(T, 2) / Tr,
                R = AS.gas_constant(), delta = rho / rhor, tau = Tr / T;

    switch (index) {
        case iT:
            dT = 1;
            drho = 0;
            break;
        case iDmolar:
            dT = 0;
            drho = 1;
            break;
        case iDmass:
            dT = 0;
            drho = AS.molar_mass();
            break;
        case iP: {
            // dp/drho|T
            drho = R * T * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2());
            // dp/dT|rho
            dT = rho * R * (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau());
            break;
        }
        case iHmass:
        case iHmolar: {
            // dh/dT|rho
            dT = R
                 * (-pow(tau, 2) * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2())
                    + (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()));
            // dh/drhomolar|T
            drho = T * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau() + delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2());
            if (index == iHmass) {
                // dhmolar/drhomolar|T * dhmass/dhmolar where dhmass/dhmolar = 1/mole mass
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case iSmass:
        case iSmolar: {
            // ds/dT|rho
            dT = R / T * (-pow(tau, 2) * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()));
            // ds/drho|T
            drho = R / rho * (-(1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()));
            if (index == iSmass) {
                // ds/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case iUmass:
        case iUmolar: {
            // du/dT|rho
            dT = R * (-pow(tau, 2) * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()));
            // du/drho|T
            drho = AS.T() * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau());
            if (index == iUmass) {
                // du/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case iGmass:
        case iGmolar: {
            // dg/dT|rho
            double dTau_dT = 1 / dT_dtau;
            dT = R * AS.T() * (AS.dalpha0_dTau() + AS.dalphar_dTau() + AS.delta() * AS.d2alphar_dDelta_dTau()) * dTau_dT
                 + R * (1 + AS.alpha0() + AS.alphar() + AS.delta() * AS.dalphar_dDelta());
            // dg/drho|T
            double dDelta_drho = 1 / rhor;
            drho = AS.T() * R * (AS.dalpha0_dDelta() + AS.dalphar_dDelta() + AS.delta() * AS.d2alphar_dDelta2() + AS.dalphar_dDelta()) * dDelta_drho;
            if (index == iGmass) {
                // dg/drho|T / drhomass/drhomolar where drhomass/drhomolar = mole mass
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case iTau:
            dT = 1 / dT_dtau;
            drho = 0;
            break;
        case iDelta:
            dT = 0;
            drho = 1 / rhor;
            break;
        case iCvmolar:
        case iCvmass: {
            // use the second order derivative of internal energy
            // make it cleaner by using the function get_dT_drho_second_derivatives directly?
            // dcvdT|rho = d2u/dT2|rho
            dT = R / T * pow(tau, 2) * (tau * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 2 * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()));
            // dcvdrho|T = d2u/dT/drho
            drho = R / rho * (-pow(tau, 2) * delta * AS.d3alphar_dDelta_dTau2());
            if (index == iCvmass) {
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case iCpmolar:
        case iCpmass: {
            // dcp/dT|rho = d2h/dT2 + dh/drho * dP/dT * d2P/drhodT / ( dp/drho )^2 - ( d2h/dTdrho * dP/dT + dh/drho * d2P/dT2 ) / ( dP/drho )
            dT = R / T * pow(tau, 2)
                 * (tau * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 2 * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2())
                    + delta * AS.d3alphar_dDelta_dTau2());
            dT += (T * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau() + delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()))
                  * (rho * R * (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()))
                  * (R
                     * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2() - 2 * delta * tau * AS.d2alphar_dDelta_dTau()
                        - tau * pow(delta, 2) * AS.d3alphar_dDelta2_dTau()))
                  / pow(R * T * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()), 2);
            dT -= ((R / rho * delta
                    * (delta * AS.d2alphar_dDelta2() - pow(tau, 2) * AS.d3alphar_dDelta_dTau2() + AS.dalphar_dDelta()
                       - tau * delta * AS.d3alphar_dDelta2_dTau() - tau * AS.d2alphar_dDelta_dTau()))
                     * (rho * R * (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()))
                   + (T * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau() + delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()))
                       * (rho * R / T * (pow(tau, 2) * delta * AS.d3alphar_dDelta_dTau2())))
                  / (R * T * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()));
            // dcpdrho|T = d2h/dTdrho + dh/drho * dP/dT * d2P/drho2 / ( dp/drho )^2 - ( d2h/drho2 * dP/dT + dh/drho * d2P/dTdrho ) / ( dP/drho )
            drho = R / rho * delta
                   * (delta * AS.d2alphar_dDelta2() - pow(tau, 2) * AS.d3alphar_dDelta_dTau2() + AS.dalphar_dDelta()
                      - tau * delta * AS.d3alphar_dDelta2_dTau() - tau * AS.d2alphar_dDelta_dTau());  //d2h/dTdrho
            drho +=
              (T * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau() + delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()))
              * (rho * R * (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()))
              * (T * R / rho * (2 * delta * AS.dalphar_dDelta() + 4 * pow(delta, 2) * AS.d2alphar_dDelta2() + pow(delta, 3) * AS.d3alphar_dDelta3()))
              / pow(R * T * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()), 2);
            drho -= ((R * T * pow(delta / rho, 2) * (tau * AS.d3alphar_dDelta2_dTau() + 2 * AS.d2alphar_dDelta2() + delta * AS.d3alphar_dDelta3()))
                       * (rho * R * (1 + delta * AS.dalphar_dDelta() - tau * delta * AS.d2alphar_dDelta_dTau()))
                     + (T * R / rho * (tau * delta * AS.d2alphar_dDelta_dTau() + delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()))
                         * (R
                            * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()
                               - 2 * delta * tau * AS.d2alphar_dDelta_dTau() - tau * pow(delta, 2) * AS.d3alphar_dDelta2_dTau())))
                    / (R * T * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2()));
            if (index == iCpmass) {
                drho /= AS.molar_mass();
                dT /= AS.molar_mass();
            }
            break;
        }
        case ispeed_sound: {
            //dwdT
            double aa = 1.0 + delta * AS.dalphar_dDelta() - delta * tau * AS.d2alphar_dDelta_dTau();
            double bb = pow(tau, 2) * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2());
            double daa_dTau = -delta * tau * AS.d3alphar_dDelta_dTau2();
            double dbb_dTau = pow(tau, 2) * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 2.0 * tau * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2());
            double w = AS.speed_sound();
            dT = 1.0 / 2.0 / w / T
                 * (pow(w, 2)
                    - R * Tr / AS.molar_mass()
                        * (2.0 * delta * AS.d2alphar_dDelta_dTau() + pow(delta, 2) * AS.d3alphar_dDelta2_dTau()
                           - (2 * aa / bb * daa_dTau - pow(aa / bb, 2) * dbb_dTau)));
            //dwdrho
            double daa_dDelta =
              AS.dalphar_dDelta() + delta * AS.d2alphar_dDelta2() - tau * (AS.d2alphar_dDelta_dTau() + delta * AS.d3alphar_dDelta2_dTau());
            double dbb_dDelta = pow(tau, 2) * (AS.d3alpha0_dDelta_dTau2() + AS.d3alphar_dDelta_dTau2());
            drho = R * T / 2.0 / AS.molar_mass() / w / rhor
                   * (2.0 * (AS.dalphar_dDelta() + delta * AS.d2alphar_dDelta2())
                      + (2.0 * delta * AS.d2alphar_dDelta2() + pow(delta, 2) * AS.d3alphar_dDelta3())
                      - (2 * aa / bb * daa_dDelta - pow(aa / bb, 2) * dbb_dDelta));
            break;
        }
        default:
            throw ValueError(format("input to get_dT_drho[%s] is invalid", get_parameter_information(index, "short").c_str()));
    }
}
void get_dT_drho_second_derivatives(AbstractState& AS, int index, CoolPropDbl& dT2, CoolPropDbl& drho_dT, CoolPropDbl& drho2) {
    CoolPropDbl T = AS.T(), rho = AS.rhomolar(), rhor = AS.rhomolar_reducing(), Tr = AS.T_reducing(), R = AS.gas_constant(), delta = rho / rhor,
                tau = Tr / T;

    // Here we use T and rho as independent variables since derivations are already done by Thorade, 2013,
    // Partial derivatives of thermodynamic state propertiesfor dynamic simulation, DOI 10.1007/s12665-013-2394-z

    switch (index) {
        case iT:
        case iDmass:
        case iDmolar:
            dT2 = 0;  // d2rhomolar_dtau2
            drho2 = 0;
            drho_dT = 0;
            break;
        case iTau:
            dT2 = 2 * Tr / pow(T, 3);
            drho_dT = 0;
            drho2 = 0;
            break;
        case iDelta:
            dT2 = 0;
            drho_dT = 0;
            drho2 = 0;
            break;
        case iP: {
            drho2 =
              T * R / rho * (2 * delta * AS.dalphar_dDelta() + 4 * pow(delta, 2) * AS.d2alphar_dDelta2() + pow(delta, 3) * AS.d3alphar_dDelta3());
            dT2 = rho * R / T * (pow(tau, 2) * delta * AS.d3alphar_dDelta_dTau2());
            drho_dT = R
                      * (1 + 2 * delta * AS.dalphar_dDelta() + pow(delta, 2) * AS.d2alphar_dDelta2() - 2 * delta * tau * AS.d2alphar_dDelta_dTau()
                         - tau * pow(delta, 2) * AS.d3alphar_dDelta2_dTau());
            break;
        }
        case iHmass:
        case iHmolar: {
            // d2h/drho2|T
            drho2 = R * T * pow(delta / rho, 2) * (tau * AS.d3alphar_dDelta2_dTau() + 2 * AS.d2alphar_dDelta2() + delta * AS.d3alphar_dDelta3());
            // d2h/dT2|rho
            dT2 = R / T * pow(tau, 2)
                  * (tau * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 2 * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2())
                     + delta * AS.d3alphar_dDelta_dTau2());
            // d2h/drho/dT
            drho_dT = R / rho * delta
                      * (delta * AS.d2alphar_dDelta2() - pow(tau, 2) * AS.d3alphar_dDelta_dTau2() + AS.dalphar_dDelta()
                         - tau * delta * AS.d3alphar_dDelta2_dTau() - tau * AS.d2alphar_dDelta_dTau());
            if (index == iHmass) {
                drho2 /= AS.molar_mass();
                drho_dT /= AS.molar_mass();
                dT2 /= AS.molar_mass();
            }
            break;
        }
        case iSmass:
        case iSmolar: {
            // d2s/rho2|T
            drho2 = R / pow(rho, 2) * (1 - pow(delta, 2) * AS.d2alphar_dDelta2() + tau * pow(delta, 2) * AS.d3alphar_dDelta2_dTau());
            // d2s/dT2|rho
            dT2 = R * pow(tau / T, 2) * (tau * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 3 * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()));
            // d2s/drho/dT
            drho_dT = R / (T * rho) * (-pow(tau, 2) * delta * AS.d3alphar_dDelta_dTau2());
            if (index == iSmass) {
                drho2 /= AS.molar_mass();
                drho_dT /= AS.molar_mass();
                dT2 /= AS.molar_mass();
            }
            break;
        }
        case iUmass:
        case iUmolar: {
            // d2u/rho2|T
            drho2 = R * T * tau * pow(delta / rho, 2) * AS.d3alphar_dDelta2_dTau();
            // d2u/dT2|rho
            dT2 = R / T * pow(tau, 2) * (tau * (AS.d3alpha0_dTau3() + AS.d3alphar_dTau3()) + 2 * (AS.d2alpha0_dTau2() + AS.d2alphar_dTau2()));
            // d2u/drho/dT
            drho_dT = R / rho * (-pow(tau, 2) * delta * AS.d3alphar_dDelta_dTau2());
            if (index == iUmass) {
                drho2 /= AS.molar_mass();
                drho_dT /= AS.molar_mass();
                dT2 /= AS.molar_mass();
            }
            break;
        }
        default:
            throw ValueError(format("input to get_dT_drho_second_derivatives[%s] is invalid", get_parameter_information(index, "short").c_str()));
    }
}
CoolPropDbl AbstractState::calc_first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
    CoolPropDbl dOf_dT, dOf_drho, dWrt_dT, dWrt_drho, dConstant_dT, dConstant_drho;

    get_dT_drho(*this, Of, dOf_dT, dOf_drho);
    get_dT_drho(*this, Wrt, dWrt_dT, dWrt_drho);
    get_dT_drho(*this, Constant, dConstant_dT, dConstant_drho);

    return (dOf_dT * dConstant_drho - dOf_drho * dConstant_dT) / (dWrt_dT * dConstant_drho - dWrt_drho * dConstant_dT);
}
CoolPropDbl AbstractState::calc_second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) {
    CoolPropDbl dOf1_dT, dOf1_drho, dWrt1_dT, dWrt1_drho, dConstant1_dT, dConstant1_drho, d2Of1_dT2, d2Of1_drhodT, d2Of1_drho2, d2Wrt1_dT2,
      d2Wrt1_drhodT, d2Wrt1_drho2, d2Constant1_dT2, d2Constant1_drhodT, d2Constant1_drho2, dWrt2_dT, dWrt2_drho, dConstant2_dT, dConstant2_drho, N, D,
      dNdrho__T, dDdrho__T, dNdT__rho, dDdT__rho, dderiv1_drho, dderiv1_dT, second;

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
    N = dOf1_dT * dConstant1_drho - dOf1_drho * dConstant1_dT;
    D = dWrt1_dT * dConstant1_drho - dWrt1_drho * dConstant1_dT;

    // Derivatives of the numerator and denominator of the first partial derivative term with respect to rho, T held constant
    // They are of similar form, with Of1 and Wrt1 swapped
    dNdrho__T = dOf1_dT * d2Constant1_drho2 + d2Of1_drhodT * dConstant1_drho - dOf1_drho * d2Constant1_drhodT - d2Of1_drho2 * dConstant1_dT;
    dDdrho__T = dWrt1_dT * d2Constant1_drho2 + d2Wrt1_drhodT * dConstant1_drho - dWrt1_drho * d2Constant1_drhodT - d2Wrt1_drho2 * dConstant1_dT;

    // Derivatives of the numerator and denominator of the first partial derivative term with respect to T, rho held constant
    // They are of similar form, with Of1 and Wrt1 swapped
    dNdT__rho = dOf1_dT * d2Constant1_drhodT + d2Of1_dT2 * dConstant1_drho - dOf1_drho * d2Constant1_dT2 - d2Of1_drhodT * dConstant1_dT;
    dDdT__rho = dWrt1_dT * d2Constant1_drhodT + d2Wrt1_dT2 * dConstant1_drho - dWrt1_drho * d2Constant1_dT2 - d2Wrt1_drhodT * dConstant1_dT;

    // First partial of first derivative term with respect to T
    dderiv1_drho = (D * dNdrho__T - N * dDdrho__T) / pow(D, 2);

    // First partial of first derivative term with respect to rho
    dderiv1_dT = (D * dNdT__rho - N * dDdT__rho) / pow(D, 2);

    // Complete second derivative
    second = (dderiv1_dT * dConstant2_drho - dderiv1_drho * dConstant2_dT) / (dWrt2_dT * dConstant2_drho - dWrt2_drho * dConstant2_dT);

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

#include <catch2/catch_all.hpp>

TEST_CASE("Check AbstractState", "[AbstractState]") {
    SECTION("bad backend") {
        CHECK_THROWS(shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("DEFINITELY_A_BAD_BACKEND", "Water")));
    }
    SECTION("good backend - bad fluid") {
        CHECK_THROWS(shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "DEFINITELY_A_BAD_FLUID")));
    }
    SECTION("good backend - helmholtz") {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", "Water")));
    }
    SECTION("good backend - incomp") {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("INCOMP", "DEB")));
    }
    SECTION("good backend - REFPROP") {
        CHECK_NOTHROW(shared_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("REFPROP", "Water")));
    }
}

TEST_CASE("Check derivatives in first_partial_deriv", "[derivs_in_first_partial_deriv]") {
    shared_ptr<CoolProp::AbstractState> Water(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> WaterplusT(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> WaterminusT(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> Waterplusrho(CoolProp::AbstractState::factory("HEOS", "Water"));
    shared_ptr<CoolProp::AbstractState> Waterminusrho(CoolProp::AbstractState::factory("HEOS", "Water"));

    double dT = 1e-3, drho = 1;
    Water->update(CoolProp::PT_INPUTS, 101325, 300);
    WaterplusT->update(CoolProp::DmolarT_INPUTS, Water->rhomolar(), 300 + dT);
    WaterminusT->update(CoolProp::DmolarT_INPUTS, Water->rhomolar(), 300 - dT);
    Waterplusrho->update(CoolProp::DmolarT_INPUTS, Water->rhomolar() + drho, 300);
    Waterminusrho->update(CoolProp::DmolarT_INPUTS, Water->rhomolar() - drho, 300);

    // Numerical derivatives
    CoolPropDbl dP_dT_num = (WaterplusT->p() - WaterminusT->p()) / (2 * dT);
    CoolPropDbl dP_drho_num = (Waterplusrho->p() - Waterminusrho->p()) / (2 * drho);

    CoolPropDbl dHmolar_dT_num = (WaterplusT->hmolar() - WaterminusT->hmolar()) / (2 * dT);
    CoolPropDbl dHmolar_drho_num = (Waterplusrho->hmolar() - Waterminusrho->hmolar()) / (2 * drho);
    CoolPropDbl dHmass_dT_num = (WaterplusT->hmass() - WaterminusT->hmass()) / (2 * dT);
    CoolPropDbl dHmass_drho_num = (Waterplusrho->hmass() - Waterminusrho->hmass()) / (2 * drho);

    CoolPropDbl dSmolar_dT_num = (WaterplusT->smolar() - WaterminusT->smolar()) / (2 * dT);
    CoolPropDbl dSmolar_drho_num = (Waterplusrho->smolar() - Waterminusrho->smolar()) / (2 * drho);
    CoolPropDbl dSmass_dT_num = (WaterplusT->smass() - WaterminusT->smass()) / (2 * dT);
    CoolPropDbl dSmass_drho_num = (Waterplusrho->smass() - Waterminusrho->smass()) / (2 * drho);

    CoolPropDbl dUmolar_dT_num = (WaterplusT->umolar() - WaterminusT->umolar()) / (2 * dT);
    CoolPropDbl dUmolar_drho_num = (Waterplusrho->umolar() - Waterminusrho->umolar()) / (2 * drho);
    CoolPropDbl dUmass_dT_num = (WaterplusT->umass() - WaterminusT->umass()) / (2 * dT);
    CoolPropDbl dUmass_drho_num = (Waterplusrho->umass() - Waterminusrho->umass()) / (2 * drho);

    CoolPropDbl dGmolar_dT_num = (WaterplusT->gibbsmolar() - WaterminusT->gibbsmolar()) / (2 * dT);
    CoolPropDbl dGmolar_drho_num = (Waterplusrho->gibbsmolar() - Waterminusrho->gibbsmolar()) / (2 * drho);
    CoolPropDbl dGmass_dT_num = (WaterplusT->gibbsmass() - WaterminusT->gibbsmass()) / (2 * dT);
    CoolPropDbl dGmass_drho_num = (Waterplusrho->gibbsmass() - Waterminusrho->gibbsmass()) / (2 * drho);

    CoolPropDbl dCvmolar_dT_num = (WaterplusT->cvmolar() - WaterminusT->cvmolar()) / (2 * dT);
    CoolPropDbl dCvmolar_drho_num = (Waterplusrho->cvmolar() - Waterminusrho->cvmolar()) / (2 * drho);
    CoolPropDbl dCvmass_dT_num = (WaterplusT->cvmass() - WaterminusT->cvmass()) / (2 * dT);
    CoolPropDbl dCvmass_drho_num = (Waterplusrho->cvmass() - Waterminusrho->cvmass()) / (2 * drho);

    CoolPropDbl dCpmolar_dT_num = (WaterplusT->cpmolar() - WaterminusT->cpmolar()) / (2 * dT);
    CoolPropDbl dCpmolar_drho_num = (Waterplusrho->cpmolar() - Waterminusrho->cpmolar()) / (2 * drho);
    CoolPropDbl dCpmass_dT_num = (WaterplusT->cpmass() - WaterminusT->cpmass()) / (2 * dT);
    CoolPropDbl dCpmass_drho_num = (Waterplusrho->cpmass() - Waterminusrho->cpmass()) / (2 * drho);

    CoolPropDbl dspeed_sound_dT_num = (WaterplusT->speed_sound() - WaterminusT->speed_sound()) / (2 * dT);
    CoolPropDbl dspeed_sound_drho_num = (Waterplusrho->speed_sound() - Waterminusrho->speed_sound()) / (2 * drho);

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
    // Isochoric heat capacity
    CoolPropDbl dCvmolar_dT_analyt, dCvmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iCvmolar, dCvmolar_dT_analyt, dCvmolar_drho_analyt);
    CoolPropDbl dCvmass_dT_analyt, dCvmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iCvmass, dCvmass_dT_analyt, dCvmass_drho_analyt);
    // Isobaric heat capacity
    CoolPropDbl dCpmolar_dT_analyt, dCpmolar_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iCpmolar, dCpmolar_dT_analyt, dCpmolar_drho_analyt);
    CoolPropDbl dCpmass_dT_analyt, dCpmass_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::iCpmass, dCpmass_dT_analyt, dCpmass_drho_analyt);
    // Speed of sound
    CoolPropDbl dspeed_sound_dT_analyt, dspeed_sound_drho_analyt;
    CoolProp::get_dT_drho(*Water, CoolProp::ispeed_sound, dspeed_sound_dT_analyt, dspeed_sound_drho_analyt);

    double eps = 1e-3;

    CHECK(std::abs(dP_dT_analyt / dP_dT_num - 1) < eps);
    CHECK(std::abs(dP_drho_analyt / dP_drho_num - 1) < eps);

    CHECK(std::abs(dHmolar_dT_analyt / dHmolar_dT_num - 1) < eps);
    CHECK(std::abs(dHmolar_drho_analyt / dHmolar_drho_num - 1) < eps);
    CHECK(std::abs(dHmass_dT_analyt / dHmass_dT_num - 1) < eps);
    CHECK(std::abs(dHmass_drho_analyt / dHmass_drho_num - 1) < eps);

    CHECK(std::abs(dSmolar_dT_analyt / dSmolar_dT_num - 1) < eps);
    CHECK(std::abs(dSmolar_drho_analyt / dSmolar_drho_num - 1) < eps);
    CHECK(std::abs(dSmass_dT_analyt / dSmass_dT_num - 1) < eps);
    CHECK(std::abs(dSmass_drho_analyt / dSmass_drho_num - 1) < eps);

    CHECK(std::abs(dUmolar_dT_analyt / dUmolar_dT_num - 1) < eps);
    CHECK(std::abs(dUmolar_drho_analyt / dUmolar_drho_num - 1) < eps);
    CHECK(std::abs(dUmass_dT_analyt / dUmass_dT_num - 1) < eps);
    CHECK(std::abs(dUmass_drho_analyt / dUmass_drho_num - 1) < eps);

    CHECK(std::abs(dGmolar_dT_analyt / dGmolar_dT_num - 1) < eps);
    CHECK(std::abs(dGmolar_drho_analyt / dGmolar_drho_num - 1) < eps);
    CHECK(std::abs(dGmass_dT_analyt / dGmass_dT_num - 1) < eps);
    CHECK(std::abs(dGmass_drho_analyt / dGmass_drho_num - 1) < eps);

    CHECK(std::abs(dCvmolar_dT_analyt / dCvmolar_dT_num - 1) < eps);
    CHECK(std::abs(dCvmolar_drho_analyt / dCvmolar_drho_num - 1) < eps);
    CHECK(std::abs(dCvmass_dT_analyt / dCvmass_dT_num - 1) < eps);
    CHECK(std::abs(dCvmass_drho_analyt / dCvmass_drho_num - 1) < eps);

    CHECK(std::abs(dCpmolar_dT_analyt / dCpmolar_dT_num - 1) < eps);
    CHECK(std::abs(dCpmolar_drho_analyt / dCpmolar_drho_num - 1) < eps);
    CHECK(std::abs(dCpmass_dT_analyt / dCpmass_dT_num - 1) < eps);
    CHECK(std::abs(dCpmass_drho_analyt / dCpmass_drho_num - 1) < eps);

    CHECK(std::abs(dspeed_sound_dT_analyt / dspeed_sound_dT_num - 1) < eps);
    CHECK(std::abs(dspeed_sound_drho_analyt / dspeed_sound_drho_num - 1) < eps);
}

#endif
