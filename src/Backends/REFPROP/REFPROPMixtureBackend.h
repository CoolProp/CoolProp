/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTUREBACKEND_H_
#define REFPROPMIXTUREBACKEND_H_

#include "AbstractState.h"
#include "DataStructures.h"

#include <vector>

namespace CoolProp {

class REFPROPMixtureBackend : public AbstractState
{
   private:
    std::string cached_component_string;

   protected:
    std::size_t Ncomp;
    bool _mole_fractions_set;

    static std::size_t instance_counter;
    static bool _REFPROP_supported;
    std::vector<CoolPropDbl> mole_fractions_long_double;  // read-only
    std::vector<double> mole_fractions, mass_fractions;
    std::vector<double> mole_fractions_liq, mole_fractions_vap;
    std::vector<std::string> fluid_names;

    /// Call the PHIXdll function in the dll
    CoolPropDbl call_phixdll(int itau, int idelta);
    /// Call the PHI0dll function in the dll
    CoolPropDbl call_phi0dll(int itau, int idelta);

   public:
    REFPROPMixtureBackend() : Ncomp(0), _mole_fractions_set(false) {
        instance_counter++;
    }

    /// The instantiator
    /// @param fluid_names The vector of strings of the fluid components, without file ending
    REFPROPMixtureBackend(const std::vector<std::string>& fluid_names) {
        construct(fluid_names);
    };

    /// A function to actually do the initialization to allow it to be called in derived classes
    void construct(const std::vector<std::string>& fluid_names);

    std::string backend_name(void) {
        return get_backend_string(REFPROP_BACKEND_MIX);
    }
    virtual ~REFPROPMixtureBackend();

    static std::string version();

    std::vector<std::string> calc_fluid_names() {
        return fluid_names;
    };
    PhaseEnvelopeData PhaseEnvelope;

    /// Set binary mixture floating point parameter
    void set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter, const double value);
    /// Get binary mixture double value
    double get_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter);

    /// Get binary mixture string value
    std::string get_binary_interaction_string(const std::string& CAS1, const std::string& CAS2, const std::string& parameter);
    /// Set binary mixture string value
    void set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter, const std::string& value);

    /// Set binary mixture string parameter (EXPERT USE ONLY!!!)
    void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, const double value);
    /// Get binary mixture double value (EXPERT USE ONLY!!!)
    double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter);

    /// Find the index (1-based for FORTRAN) of the fluid with the given CAS number
    int match_CAS(const std::string& CAS);

    // REFPROP backend uses mole fractions
    bool using_mole_fractions() {
        return true;
    }
    bool using_mass_fractions() {
        return false;
    }
    bool using_volu_fractions() {
        return false;
    }

    /// Calculate the name of the fluid
    std::string calc_name() {
        return fluid_param_string("name");
    }

    // Get _phase for pure fluids only
    phases calc_phase(void) {
        if (this->Ncomp > 1) {
            throw NotImplementedError("The REFPROP backend does not implement calc_phase function for mixtures.");
        } else {
            return _phase;
        }
    };

    // Utility function to determine the phase from quality value return from REFPROP
    phases GetRPphase();

    /** \brief Specify the phase - this phase will always be used in calculations
     *
     * @param phase_index The index from CoolProp::phases
     */
    void calc_specify_phase(phases phase_index) {
        imposed_phase_index = phase_index;
        _phase = phase_index;
    }
    /**\brief Unspecify the phase - the phase is no longer imposed, different solvers can do as they like
     */
    void calc_unspecify_phase() {
        imposed_phase_index = iphase_not_imposed;
    }

    /// Updating function for REFPROP
    /**
    In this function we take a pair of thermodynamic states, those defined in the input_pairs
    enumeration and update all the internal variables that we can.  REFPROP calculates
    a lot of other state variables each time you use a flash routine so we cache all the
    outputs that we can, which saves on computational time.

    @param input_pair Integer key from CoolProp::input_pairs to the two inputs that will be passed to the function
    @param value1 First input value
    @param value2 Second input value
    */
    void update(CoolProp::input_pairs, double value1, double value2);

    /**
     * @brief Update the state, while providing guess values
     */
    void update_with_guesses(CoolProp::input_pairs, double value1, double value2, const GuessesStructure& guesses);

    CoolPropDbl calc_molar_mass(void);

    void check_loaded_fluid(void);

    void calc_excess_properties();

    /// Returns true if REFPROP is supported on this platform
    static bool REFPROP_supported(void);

    std::string fluid_param_string(const std::string& ParamName);

    CoolPropDbl calc_PIP(void);

    CoolPropDbl calc_cpmolar_idealgas(void);

    /// Set the fluids in REFPROP DLL by calling the SETUPdll function
    /**
    @param fluid_names The vector of strings of the fluid components, without file ending
    */
    void set_REFPROP_fluids(const std::vector<std::string>& fluid_names);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions);

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions);

    const std::vector<CoolPropDbl>& get_mole_fractions() {
        return mole_fractions_long_double;
    };

    const std::vector<CoolPropDbl> calc_mass_fractions();

    void calc_phase_envelope(const std::string& type);

    CoolPropDbl calc_compressibility_factor(void) {
        return _p / (_rhomolar * gas_constant() * _T);
    };

    const CoolProp::PhaseEnvelopeData& calc_phase_envelope_data() {
        return PhaseEnvelope;
    };

    std::vector<CoolPropDbl> calc_mole_fractions_liquid(void) {
        return std::vector<CoolPropDbl>(mole_fractions_liq.begin(), mole_fractions_liq.begin() + this->Ncomp);
    }
    std::vector<CoolPropDbl> calc_mole_fractions_vapor(void) {
        return std::vector<CoolPropDbl>(mole_fractions_vap.begin(), mole_fractions_vap.begin() + this->Ncomp);
    }

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /// Get the viscosity [Pa-s] (based on the temperature and density in the state class)
    CoolPropDbl calc_viscosity(void);
    /// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
    CoolPropDbl calc_conductivity(void);
    /// Get the surface tension [N/m] (based on the temperature in the state class).  Invalid for temperatures above critical point or below triple point temperature
    CoolPropDbl calc_surface_tension(void);
    /// Calc the B virial coefficient
    CoolPropDbl calc_Bvirial(void);
    /// Calc the temperature derivative of the second virial coefficient
    CoolPropDbl calc_dBvirial_dT(void);
    /// Calc the C virial coefficient
    CoolPropDbl calc_Cvirial(void);

    CoolPropDbl calc_fugacity_coefficient(std::size_t i);
    CoolPropDbl calc_fugacity(std::size_t i);
    CoolPropDbl calc_chemical_potential(std::size_t i);
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value);
    bool has_melting_line();
    double calc_melt_Tmax();
    CoolPropDbl calc_T_critical(void);
    CoolPropDbl calc_T_reducing(void);
    void calc_reducing_state(void);
    CoolPropDbl calc_p_critical(void);
    CoolPropDbl calc_p_triple(void);
    CoolPropDbl calc_p_min(void) {
        return calc_p_triple();
    };
    CoolPropDbl calc_rhomolar_critical(void);
    CoolPropDbl calc_rhomolar_reducing(void);
    CoolPropDbl calc_Ttriple(void);
    CoolPropDbl calc_acentric_factor(void);
    CoolPropDbl calc_gas_constant(void);
    CoolPropDbl calc_dipole_moment(void);

    /// Calculate the "true" critical point where dp/drho|T and d2p/drho2|T are zero
    void calc_true_critical_point(double& T, double& rho);

    /// Calculate the saturation properties
    CoolPropDbl calc_saturated_liquid_keyed_output(parameters key);
    CoolPropDbl calc_saturated_vapor_keyed_output(parameters key);

    /// Calculate an ideal curve
    void calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p);

    /// A wrapper function to calculate the limits for the EOS
    void limits(double& Tmin, double& Tmax, double& rhomolarmax, double& pmax);
    /// Calculate the maximum pressure
    CoolPropDbl calc_pmax(void);
    /// Calculate the maximum temperature
    CoolPropDbl calc_Tmax(void);
    /// Calculate the minimum temperature
    CoolPropDbl calc_Tmin(void);

    /// Calculate the residual entropy in J/mol/K (should be a uniquely negative quantity)
    CoolPropDbl calc_smolar_residual(void) {
        return (tau() * calc_dalphar_dTau() - calc_alphar()) * gas_constant();
    }

    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    CoolPropDbl calc_alphar(void) {
        return call_phixdll(0, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dDelta(void) {
        return call_phixdll(0, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dTau(void) {
        return call_phixdll(1, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta2(void) {
        return call_phixdll(0, 2);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta_dTau(void) {
        return call_phixdll(1, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dTau2(void) {
        return call_phixdll(2, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta3(void) {
        return call_phixdll(0, 3);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta2_dTau(void) {
        return call_phixdll(1, 2);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta_dTau2(void) {
        return call_phixdll(2, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dTau3(void) {
        return call_phixdll(3, 0);
    };

    CoolPropDbl calc_alpha0(void) {
        return call_phi0dll(0, 0);
    };
    CoolPropDbl calc_dalpha0_dDelta(void) {
        return call_phi0dll(0, 1);
    };
    CoolPropDbl calc_dalpha0_dTau(void) {
        return call_phi0dll(1, 0);
    };
    CoolPropDbl calc_d2alpha0_dDelta2(void) {
        return call_phi0dll(0, 2);
    };
    CoolPropDbl calc_d2alpha0_dDelta_dTau(void) {
        return call_phi0dll(1, 1);
    };
    CoolPropDbl calc_d2alpha0_dTau2(void) {
        return call_phi0dll(2, 0);
    };
    CoolPropDbl calc_d3alpha0_dDelta3(void) {
        return call_phi0dll(0, 3);
    };
    CoolPropDbl calc_d3alpha0_dDelta2_dTau(void) {
        return call_phi0dll(1, 2);
    };
    CoolPropDbl calc_d3alpha0_dDelta_dTau2(void) {
        return call_phi0dll(2, 1);
    };
    CoolPropDbl calc_d3alpha0_dTau3(void) {
        return call_phi0dll(3, 0);
    };
};

bool force_load_REFPROP();
bool force_unload_REFPROP();
void REFPROP_SETREF(char hrf[3], int ixflag, double x0[1], double& h0, double& s0, double& T0, double& p0, int& ierr, char herr[255], int l1, int l2);

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
