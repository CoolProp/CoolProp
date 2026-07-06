/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTUREBACKEND_H_
#define REFPROPMIXTUREBACKEND_H_

#include "CoolProp/AbstractState.h"
#include "CoolProp/CoolPropFluid.h"
#include "CoolProp/DataStructures.h"

#include <string>
#include <vector>

namespace CoolProp {

/// Return the REFPROP .FLD stem for a CoolPropFluid.
/// Falls back to `fallback` when REFPROPname is absent or the sentinel "N/A",
/// so the original user-supplied name is preserved rather than fluid.name
/// (which may differ, e.g. "R1336mzz(E)" vs the file "R1336MZZE.FLD").
inline std::string refprop_stem(const CoolPropFluid& fluid, const std::string& fallback) {
    if (!fluid.REFPROPname.empty() && fluid.REFPROPname != "N/A") {
        return fluid.REFPROPname;
    }
    return fallback;
}

struct THERM0dllOutputs
{
    double p_kPa;        /// Pressure [kPa]
    double umol_Jmol;    /// Internal energy [J/mol]
    double hmol_Jmol;    /// Enthalpy [J/mol]
    double smol_JmolK;   /// Entropy [J/mol-K]
    double cvmol_JmolK;  /// Isochoric heat capacity [J/mol-K]
    double cpmol_JmolK;  /// Isobaric heat capacity [J/mol-K]
    double w_ms;         /// Speed of sound [m/s]
    double amol_Jmol;    ///Helmholtz energy [J/mol]
    double gmol_Jmol;    /// Gibbs free energy [J/mol]
};

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

    /// Configure *this* as a shim linked to `host`'s already-loaded REFPROP fluids (no SETUPdll).
    void link_to_loaded_fluids(const REFPROPMixtureBackend& host);

    /// dP/dT [Pa/K] along the pure-component saturation line via DPTSATKdll. kph: 1=liquid, 2=vapor.
    double dpdT_along_saturation_pure(int kph);
    /// Saturation pressure [Pa] at temperature T for the bubble (Q=0) / dew (Q=1) branch via SATTdll.
    double saturation_pressure_at_T(double T, int Q);

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

    std::string backend_name() override {
        return get_backend_string(REFPROP_BACKEND_MIX);
    }
    virtual ~REFPROPMixtureBackend();

    static std::string version();

    std::vector<std::string> calc_fluid_names() override {
        return fluid_names;
    };
    PhaseEnvelopeData PhaseEnvelope;

    /// Set binary mixture floating point parameter
    void set_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter, const double value) override;
    /// Get binary mixture double value
    double get_binary_interaction_double(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) override;

    /// Get binary mixture string value
    std::string get_binary_interaction_string(const std::string& CAS1, const std::string& CAS2, const std::string& parameter) override;
    /// Set binary mixture string value
    void set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter, const std::string& value) override;

    /// Set binary mixture string parameter (EXPERT USE ONLY!!!)
    void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, const double value) override;
    /// Get binary mixture double value (EXPERT USE ONLY!!!)
    double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) override;

    /// Find the index (1-based for FORTRAN) of the fluid with the given CAS number
    int match_CAS(const std::string& CAS);

    // REFPROP backend uses mole fractions
    bool using_mole_fractions() override {
        return true;
    }
    bool using_mass_fractions() override {
        return false;
    }
    bool using_volu_fractions() override {
        return false;
    }

    /// Calculate the name of the fluid
    std::string calc_name() override {
        return fluid_param_string("name");
    }

    // Get _phase for pure fluids only
    phases calc_phase() override {
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
    void calc_specify_phase(phases phase_index) override {
        imposed_phase_index = phase_index;
        _phase = phase_index;
    }
    /**\brief Unspecify the phase - the phase is no longer imposed, different solvers can do as they like
     */
    void calc_unspecify_phase() override {
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
    void update(CoolProp::input_pairs, double value1, double value2) override;

    void update_Qmass_pair(CoolProp::input_pairs pair, double v1, double v2) override;

    /**
     * @brief Update the state, while providing guess values
     */
    void update_with_guesses(CoolProp::input_pairs, double value1, double value2, const GuessesStructure& guesses) override;

    /// Set the state directly from (rho_molar, T) using THERMdll, no flash.
    /// Used internally to drive saturation/two-phase derivative shims.
    void update_DmolarT_direct(CoolPropDbl rhomolar, CoolPropDbl T);

    /// Build a fully-populated saturated-phase shim. `Q` must be 0 (liquid) or 1 (vapor).
    /// The shim reuses the fluids already loaded by this instance (no SETUPdll) and
    /// carries the correct per-phase composition (mole_fractions_liq / mole_fractions_vap).
    shared_ptr<REFPROPMixtureBackend> build_saturation_shim(int Q);

    CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1) override;

    CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) override;
    CoolPropDbl calc_second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) override;
    CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) override;

    CoolPropDbl calc_molar_mass() override;

    PhaseMolarMasses calc_phase_molar_masses() override;

    void check_loaded_fluid();

    void calc_excess_properties() override;

    /// Returns true if REFPROP is supported on this platform
    static bool REFPROP_supported();

    std::string fluid_param_string(const std::string& ParamName) override;

    CoolPropDbl calc_PIP() override;

    CoolPropDbl calc_cpmolar_idealgas() override;

    /// Set the fluids in REFPROP DLL by calling the SETUPdll function
    /**
    @param fluid_names The vector of strings of the fluid components, without file ending
    */
    void set_REFPROP_fluids(const std::vector<std::string>& fluid_names);

    /// Set the mole fractions
    /**
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<CoolPropDbl>& mole_fractions) override;

    /// Set the mass fractions
    /**
    @param mass_fractions The vector of mass fractions of the components
    */
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) override;

    const std::vector<CoolPropDbl>& get_mole_fractions() override {
        return mole_fractions_long_double;
    };

    const std::vector<CoolPropDbl> calc_mass_fractions() override;

    void calc_phase_envelope(const std::string& type) override;

    CoolPropDbl calc_compressibility_factor() override {
        return _p / (_rhomolar * gas_constant() * _T);
    };

    const CoolProp::PhaseEnvelopeData& calc_phase_envelope_data() override {
        return PhaseEnvelope;
    };

    std::vector<CoolPropDbl> calc_mole_fractions_liquid() override {
        return {mole_fractions_liq.begin(), mole_fractions_liq.begin() + this->Ncomp};
    }
    std::vector<CoolPropDbl> calc_mole_fractions_vapor() override {
        return {mole_fractions_vap.begin(), mole_fractions_vap.begin() + this->Ncomp};
    }

    /// Check if the mole fractions have been set, etc.
    void check_status();

    /// Get the viscosity [Pa-s] (based on the temperature and density in the state class)
    CoolPropDbl calc_viscosity() override;
    /// Get the thermal conductivity [W/m/K] (based on the temperature and density in the state class)
    CoolPropDbl calc_conductivity() override;
    /// Get the surface tension [N/m] (based on the temperature in the state class).  Invalid for temperatures above critical point or below triple point temperature
    CoolPropDbl calc_surface_tension() override;
    /// Calc the B virial coefficient
    CoolPropDbl calc_Bvirial() override;
    /// Calc the temperature derivative of the second virial coefficient
    CoolPropDbl calc_dBvirial_dT() override;
    /// Calc the C virial coefficient
    CoolPropDbl calc_Cvirial() override;

    CoolPropDbl calc_fugacity_coefficient(std::size_t i) override;
    CoolPropDbl calc_fugacity(std::size_t i) override;
    CoolPropDbl calc_chemical_potential(std::size_t i) override;
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value) override;
    bool has_melting_line() override;
    double calc_melt_Tmax();
    CoolPropDbl calc_T_critical() override;
    CoolPropDbl calc_T_reducing() override;
    void calc_reducing_state() override;
    CoolPropDbl calc_p_critical() override;
    CoolPropDbl calc_p_triple() override;
    CoolPropDbl calc_p_min() {
        return calc_p_triple();
    };
    CoolPropDbl calc_rhomolar_critical() override;
    CoolPropDbl calc_rhomolar_reducing() override;
    CoolPropDbl calc_Ttriple() override;
    CoolPropDbl calc_acentric_factor() override;
    CoolPropDbl calc_gas_constant() override;
    CoolPropDbl calc_dipole_moment() override;

    /// Calculate the "true" critical point where dp/drho|T and d2p/drho2|T are zero
    void calc_true_critical_point(double& T, double& rho) override;

    /// Calculate the saturation properties
    CoolPropDbl calc_saturated_liquid_keyed_output(parameters key) override;
    CoolPropDbl calc_saturated_vapor_keyed_output(parameters key) override;

    /// Calculate an ideal curve
    void calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) override;

    /// A wrapper function to calculate the limits for the EOS
    void limits(double& Tmin, double& Tmax, double& rhomolarmax, double& pmax);
    /// Calculate the maximum pressure
    CoolPropDbl calc_pmax() override;
    /// Calculate the maximum temperature
    CoolPropDbl calc_Tmax() override;
    /// Calculate the minimum temperature
    CoolPropDbl calc_Tmin() override;

    /// Call into the THERM0dll method and return outputs as a struct. REFPROP must already be setup
    THERM0dllOutputs call_THERM0dll(double T, double rho_mol_dm3, const std::vector<double>& mole_fractions);

    /// Calculate the residual entropy in J/mol/K (should be a uniquely negative quantity)
    CoolPropDbl calc_smolar_residual() override {
        return (tau() * calc_dalphar_dTau() - calc_alphar()) * gas_constant();
    }

    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    CoolPropDbl calc_alphar() override {
        return call_phixdll(0, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dDelta() override {
        return call_phixdll(0, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dTau() override {
        return call_phixdll(1, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta2() override {
        return call_phixdll(0, 2);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta_dTau() override {
        return call_phixdll(1, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dTau2() override {
        return call_phixdll(2, 0);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta3() override {
        return call_phixdll(0, 3);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta2_dTau() override {
        return call_phixdll(1, 2);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta_dTau2() override {
        return call_phixdll(2, 1);
    };
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dTau3() override {
        return call_phixdll(3, 0);
    };

    CoolPropDbl calc_alpha0() override {
        return call_phi0dll(0, 0);
    };
    CoolPropDbl calc_dalpha0_dDelta() override {
        return call_phi0dll(0, 1);
    };
    CoolPropDbl calc_dalpha0_dTau() override {
        return call_phi0dll(1, 0);
    };
    CoolPropDbl calc_d2alpha0_dDelta2() override {
        return call_phi0dll(0, 2);
    };
    CoolPropDbl calc_d2alpha0_dDelta_dTau() override {
        return call_phi0dll(1, 1);
    };
    CoolPropDbl calc_d2alpha0_dTau2() override {
        return call_phi0dll(2, 0);
    };
    CoolPropDbl calc_d3alpha0_dDelta3() override {
        return call_phi0dll(0, 3);
    };
    CoolPropDbl calc_d3alpha0_dDelta2_dTau() override {
        return call_phi0dll(1, 2);
    };
    CoolPropDbl calc_d3alpha0_dDelta_dTau2() override {
        return call_phi0dll(2, 1);
    };
    CoolPropDbl calc_d3alpha0_dTau3() override {
        return call_phi0dll(3, 0);
    };
};

bool force_load_REFPROP();
bool force_unload_REFPROP();
// Free-function wrapper around REFPROP's SETREFdll. The char/double pointer
// parameters were originally written as `char hrf[3]` / `double x0[1]` /
// `char herr[255]` — those decay to plain pointers (the [N] is doc-only,
// not enforced). Use the explicit `*` form per modernize-avoid-c-arrays
// (#2869) to make the decay visible at the call site.
void REFPROP_SETREF(char* hrf, int ixflag, double* x0, double& h0, double& s0, double& T0, double& p0, int& ierr, char* herr, int l1, int l2);

} /* namespace CoolProp */
#endif /* REFPROPMIXTUREBACKEND_H_ */
