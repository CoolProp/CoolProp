
#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "AbstractState.h"
#include "CoolPropFluid.h"
#include "ReducingFunctions.h"
#include "ExcessHEFunction.h"
#include "Solvers.h"
#include "PhaseEnvelope.h"
#include "DataStructures.h"
#include "Configuration.h"

#include <vector>

namespace CoolProp {

class FlashRoutines;

class ResidualHelmholtz;

// This class contains the mole fractions for a given mixture.
class MoleFractions
{
   private:
    std::vector<CoolPropDbl> mole_fractions;  ///< The bulk mole fractions of the mixture
    template <typename T>
    bool verify_mole_fractions_set(T i) const {
        if (i >= mole_fractions.size()) {
            throw CoolProp::ValueError("mole fractions are not set for all components");
        }
        return true;
    }

   public:
    template <typename T>
    void resize(T N) {
        return mole_fractions.resize(N);
    }
    std::size_t size() const {
        return mole_fractions.size();
    }
    void clear() {
        mole_fractions.clear();
    }
    // operator overloads
    template <typename T>
    MoleFractions& operator=(const std::vector<T>& values) {
        mole_fractions = values;
        return *this;
    }
    template <typename T>
    CoolPropDbl operator[](T i) const {
        verify_mole_fractions_set(i);
        return mole_fractions[i];
    }
    operator std::vector<CoolPropDbl>&() {
        return mole_fractions;
    }
};

class HelmholtzEOSMixtureBackend : public AbstractState
{

   protected:
    void pre_update(CoolProp::input_pairs& input_pair, CoolPropDbl& value1, CoolPropDbl& value2);
    void post_update(bool optional_checks = true);
    std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>
      linked_states;  ///< States that are linked to this one, and should be updated (BIP, reference state, etc.)
    shared_ptr<HelmholtzEOSMixtureBackend> transient_pure_state;  ///< A temporary state used for calculations of pure fluid properties
    shared_ptr<HelmholtzEOSMixtureBackend> TPD_state;             ///< A temporary state used for calculations of the tangent-plane-distance
    shared_ptr<HelmholtzEOSMixtureBackend> critical_state;        ///< A temporary state used for calculations of the critical point(s)
    /// Update the state class used to calculate the tangent-plane-distance
    virtual void add_TPD_state() {
        if (TPD_state.get() == nullptr) {
            bool sat_states = false;
            TPD_state.reset(get_copy(sat_states));
            linked_states.push_back(TPD_state);
        }
    };
    /// Update the state class used to calculate the critical point(s)
    virtual void add_critical_state() {
        if (critical_state.get() == nullptr) {
            bool sat_states = true;
            critical_state.reset(get_copy(sat_states));
            linked_states.push_back(critical_state);
        }
    };
    /// Update the state class used to calculate the critical point(s)
    virtual void add_transient_pure_state() {
        if (transient_pure_state.get() == nullptr) {
            bool sat_states = true;
            transient_pure_state.reset(get_copy(sat_states));
            linked_states.push_back(transient_pure_state);
        }
    };

    std::vector<CoolPropFluid> components;  ///< The components that are in use
    bool is_pure_or_pseudopure;             ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    MoleFractions mole_fractions;           ///< The bulk mole fractions of the mixture
    std::vector<CoolPropDbl> K,             ///< The K factors for the components
      lnK;                                  ///< The natural logarithms of the K factors of the components

    SimpleState _crit;
    std::size_t N;  ///< Number of components

    /// This overload is protected because it doesn't follow the base class definition, since this function is needed for constructing spinodals
    std::vector<CoolProp::CriticalState> _calc_all_critical_points(bool find_critical_points = true);

    static void set_fluid_enthalpy_entropy_offset(CoolPropFluid& component, double delta_a1, double delta_a2, const std::string& ref);

   public:
    /// Return the SuperAncillary expansion for the (pure) fluid, or
    /// nullptr if no SuperAncillary ships with this fluid.  Throws
    /// ValueError on mixtures / pseudo-pure fluids.
    std::shared_ptr<EquationOfState::SuperAncillary_t> get_superanc();

    /// Lazily build saturation superancillaries for caloric properties (h_sat, s_sat,
    /// both liquid and vapor branches) using the existing rho_sat superancillary's
    /// piece structure. Each entry in the new ChebyshevApproximation1D is constructed
    /// by sampling (T, rho_sat(T)) -> EOS at the rho-superancillary's Chebyshev nodes.
    /// At construction time the ChebyshevApproximation1D's companion-matrix extrema
    /// finder identifies non-monotonic regions, so subsequent inversion via
    /// get_x_for_y enumerates all roots. Reference state is captured at first build —
    /// see GitHub #2773. This is a no-op for fluids without a superancillary.
    void ensure_caloric_superancillaries();

    HelmholtzEOSMixtureBackend();
    HelmholtzEOSMixtureBackend(const std::vector<CoolPropFluid>& components, bool generate_SatL_and_SatV = true);
    HelmholtzEOSMixtureBackend(const std::vector<std::string>& component_names, bool generate_SatL_and_SatV = true);
    virtual HelmholtzEOSMixtureBackend* get_copy(bool generate_SatL_and_SatV = true);

    // Copy over the reducing and departure terms to all linked states (recursively)
    void sync_linked_states(const HelmholtzEOSMixtureBackend* const);

    virtual ~HelmholtzEOSMixtureBackend() {};
    std::string backend_name() override {
        return get_backend_string(HEOS_BACKEND_MIX);
    }
    shared_ptr<ReducingFunction> Reducing;
    shared_ptr<ResidualHelmholtz> residual_helmholtz;
    PhaseEnvelopeData PhaseEnvelope;
    SimpleState hsat_max;
    SsatSimpleState ssat_max;
    SpinodalData spinodal_values;

    bool clear() override {
        // Clear the locally cached values for the derivatives of the Helmholtz energy
        // in each component
        for (std::vector<CoolPropFluid>::iterator it = components.begin(); it != components.end(); ++it) {
            (*it).EOS().alphar.clear();
            (*it).EOS().alpha0.clear();
        }
        return AbstractState::clear();
    };

    friend class
      FlashRoutines;  // Allows the static methods in the FlashRoutines class to have access to all the protected members and methods of this class
    friend class
      TransportRoutines;  // Allows the static methods in the TransportRoutines class to have access to all the protected members and methods of this class
    friend class
      MixtureDerivatives;  // Allows the static methods in the MixtureDerivatives class to have access to all the protected members and methods of this class
    friend class
      PhaseEnvelopeRoutines;  // Allows the static methods in the PhaseEnvelopeRoutines class to have access to all the protected members and methods of this class
    friend class
      MixtureParameters;  // Allows the static methods in the MixtureParameters class to have access to all the protected members and methods of this class
    friend class
      CorrespondingStatesTerm;  // // Allows the methods in the CorrespondingStatesTerm class to have access to all the protected members and methods of this class

    // Helmholtz EOS backend uses mole fractions
    bool using_mole_fractions() override {
        return true;
    }
    bool using_mass_fractions() override {
        return false;
    }
    bool using_volu_fractions() override {
        return false;
    }
    bool is_pure() {
        return components.size() == 1 && !components[0].EOS().pseudo_pure;
    }
    bool has_melting_line() override {
        return is_pure_or_pseudopure && components[0].ancillaries.melting_line.enabled();
    };
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value) override;
    /// Return a string from the backend for the mixture/fluid
    std::string fluid_param_string(const std::string&) override;

    /// brief Set the reference state based on a string representation
    void set_reference_stateS(const std::string& reference_state) override;

    /// Set the reference state based on a thermodynamic state point specified by temperature and molar density
    void set_reference_stateD(double T, double rhomolar, double hmolar0, double smolar0) override;

    /// Set binary mixture floating point parameter
    void set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter, const double value) override;
    /// Get binary mixture double value
    double get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) override;
    ///// Get binary mixture string value
    //virtual std::string get_binary_interaction_string(const std::size_t &i, const std::size_t &j, const std::string &parameter);
    /// Set a binary interaction string
    void set_binary_interaction_string(const std::size_t i, const std::size_t j, const std::string& parameter, const std::string& value) override;
    /// Apply a simple mixing rule
    void apply_simple_mixing_rule(std::size_t i, std::size_t j, const std::string& model) override;

    // Set the cubic alpha function's constants:
    void set_cubic_alpha_C(const size_t i, const std::string& parameter, const double c1, const double c2, const double c3) override {
        throw ValueError("set_cubic_alpha_C only defined for cubic backends");
    };

    // Set fluid parameter (currently the volume translation parameter for cubic)
    void set_fluid_parameter_double(const size_t i, const std::string& parameter, const double value) override {
        throw ValueError("set_fluid_parameter_double only defined for cubic backends");
    };
    double get_fluid_parameter_double(const size_t i, const std::string& parameter) override;

    phases calc_phase() override {
        return _phase;
    };

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
    CoolPropDbl calc_saturation_ancillary(parameters param, int Q, parameters given, double value) override;
    void calc_ssat_max();
    void calc_hsat_max();
    CoolPropDbl calc_GWP20() override;
    CoolPropDbl calc_GWP500() override;
    CoolPropDbl calc_GWP100() override;
    CoolPropDbl calc_ODP() override;

    CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1) override;
    CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1, HelmholtzEOSMixtureBackend& SatL, HelmholtzEOSMixtureBackend& SatV);
    CoolPropDbl calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Wrt2) override;
    CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) override;
    CoolPropDbl calc_second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) override;
    CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end) override;

    CriticalState calc_critical_point(double rho0, double T0);
    /// An overload to make the compiler (clang in this case) happy
    std::vector<CoolProp::CriticalState> calc_all_critical_points() override {
        bool find_critical_points = true;
        return _calc_all_critical_points(find_critical_points);
    };

    virtual void get_critical_point_starting_values(double& delta0, double& tau0) {
        delta0 = get_config_double(SPINODAL_MINIMUM_DELTA);  // The value of delta where we start searching for crossing with Lstar=0 contour
        tau0 = 0.66;                                         // The value of tau where we start searching at delta=delta0
    }
    /// Get the search radius in delta and tau for the tracer
    virtual void get_critical_point_search_radii(double& R_delta, double& R_tau);
    /// Checking function to see if we should stop the tracing of the critical contour
    virtual bool get_critical_is_terminated(double& delta, double& tau) {
        return delta > 5 || tau > 5;
    }

    /// Build the spinodal curve
    void calc_build_spinodal() override;

    /// Get the data from the spinodal curve
    SpinodalData calc_get_spinodal_data() override {
        return spinodal_values;
    };

    /// Calculate the values \f$\mathcal{L}_1^*\f$ and \f$\mathcal{M}_1^*\f$
    void calc_criticality_contour_values(double& L1star, double& M1star) override;

    /// Calculate tangent plane distance for given trial composition w
    double calc_tangent_plane_distance(const double T, const double p, const std::vector<double>& w, const double rhomolar_guess) override;

    /// Calculate the phase once the state is fully calculated but you aren't sure if it is liquid or gas or ...
    void recalculate_singlephase_phase();

    /// Change the equation of state for one component
    void calc_change_EOS(const std::size_t i, const std::string& EOS_name) override;

    const CoolProp::SimpleState& calc_state(const std::string& state) override;

    const double get_fluid_constant(std::size_t i, parameters param) const override {
        const CoolPropFluid& fld = components[i];
        switch (param) {
            case iP_critical:
                return fld.crit.p;
            case iT_critical:
                return fld.crit.T;
            case iT_reducing:
                return fld.EOS().reduce.T;
            case irhomolar_reducing:
                return fld.EOS().reduce.rhomolar;
            case irhomolar_critical:
                return fld.crit.rhomolar;
            case iacentric_factor:
                return fld.EOS().acentric;
            case imolar_mass:
                return fld.EOS().molar_mass;
            case iT_triple:
                return fld.EOS().sat_min_liquid.T;
            case iP_triple:
                return fld.EOS().sat_min_liquid.p;
            case igas_constant:
                return fld.EOS().R_u;
            default:
                throw ValueError(format("I don't know what to do with this fluid constant: %s", get_parameter_information(param, "short").c_str()));
        }
    }

    const std::vector<CoolPropFluid>& get_components() const {
        return components;
    }
    std::vector<CoolPropFluid>& get_components() {
        return components;
    }
    std::vector<CoolPropDbl>& get_K() {
        return K;
    };
    std::vector<CoolPropDbl>& get_lnK() {
        return lnK;
    };
    HelmholtzEOSMixtureBackend& get_SatL() {
        return *SatL;
    };
    HelmholtzEOSMixtureBackend& get_SatV() {
        return *SatV;
    };

    std::vector<CoolPropDbl> calc_mole_fractions_liquid() override {
        // SatL/SatV retain composition vectors from the most recent VLE flash
        // (or phase-envelope build) — those are not meaningful when the
        // current state is single-phase, and returning them silently surfaced
        // as bug #2308. Require the current state to be two-phase before
        // returning anything.
        if (_phase != iphase_twophase) {
            throw ValueError("mole_fractions_liquid is only defined in the two-phase region (current state is single-phase)");
        }
        return SatL->get_mole_fractions();
    };
    std::vector<CoolPropDbl> calc_mole_fractions_vapor() override {
        if (_phase != iphase_twophase) {
            throw ValueError("mole_fractions_vapor is only defined in the two-phase region (current state is single-phase)");
        }
        return SatV->get_mole_fractions();
    };

    const std::vector<CoolPropDbl> calc_mass_fractions() override;

    const CoolProp::PhaseEnvelopeData& calc_phase_envelope_data() override {
        return PhaseEnvelope;
    };

    /// Calculate the conformal state (unity shape factors starting point if T < 0 and rhomolar < 0)
    void calc_conformal_state(const std::string& reference_fluid, CoolPropDbl& T, CoolPropDbl& rhomolar) override;

    void resize(std::size_t N);
    shared_ptr<HelmholtzEOSMixtureBackend> SatL, SatV;  ///<

    /** \brief The standard update function
     * @param input_pair The pair of inputs that will be provided
     * @param value1 The first input value
     * @param value2 The second input value
     */
    void update(CoolProp::input_pairs input_pair, double value1, double value2) override;

    /** \brief Update the state using guess values
	 *
	 */
    void update_with_guesses(CoolProp::input_pairs input_pair, double Value1, double Value2, const GuessesStructure& guesses) override;

    /** \brief Update all the internal variables for a state by copying from another state
     */
    void update_internal(HelmholtzEOSMixtureBackend& HEOS);

    /** \brief Update with TP and a guess for rho
     * @param T Temperature in K
     * @param p Pressure in Pa
     * @param rho_guess Density in mol/m^3 guessed
     */
    void update_TP_guessrho(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess);
    void update_DmolarT_direct(CoolPropDbl rhomolar, CoolPropDbl T);
    void update_TDmolarP_unchecked(CoolPropDbl T, CoolPropDbl rhomolarL, CoolPropDbl p);
    void update_QT_pure_superanc(CoolPropDbl Q, CoolPropDbl T) override;
    void update_HmolarQ_with_guessT(CoolPropDbl hmolar, CoolPropDbl Q, CoolPropDbl Tguess);

    /** \brief Set the components of the mixture
     *
     * @param components The components that are to be used in this mixture
     * @param generate_SatL_and_SatV true if SatL and SatV classes should be added, false otherwise.  Added so that saturation classes can be added without infinite recursion of adding saturation classes
     */
    virtual void set_components(const std::vector<CoolPropFluid>& components, bool generate_SatL_and_SatV = true);

    /** \brief Set the mixture parameters - binary pair reducing functions, departure functions, F_ij, etc.
     */
    void set_mixture_parameters();

    /** \brief Set the mole fractions
     *
     * @param mf The vector of mole fractions of the components
     */
    void set_mole_fractions(const std::vector<CoolPropDbl>& mf) override;

    const std::vector<CoolPropDbl>& get_mole_fractions() override {
        return mole_fractions;
    };
    std::vector<CoolPropDbl>& get_mole_fractions_ref() {
        return mole_fractions;
    };
    std::vector<double>& get_mole_fractions_doubleref() {
        return mole_fractions;
    }

    /** \brief Set the mass fractions
     *
     * @param mass_fractions The vector of mass fractions of the components
     */
    void set_mass_fractions(const std::vector<CoolPropDbl>& mass_fractions) override;

    void calc_ideal_curve(const std::string& type, std::vector<double>& T, std::vector<double>& p) override;

    CoolPropDbl calc_molar_mass() override;
    PhaseMolarMasses calc_phase_molar_masses() override;
    CoolPropDbl calc_gas_constant() override;
    CoolPropDbl calc_acentric_factor() override;

    CoolPropDbl calc_Bvirial() override;
    CoolPropDbl calc_Cvirial() override;
    CoolPropDbl calc_dBvirial_dT() override;
    CoolPropDbl calc_dCvirial_dT() override;

    CoolPropDbl calc_pressure() override;
    CoolPropDbl calc_cvmolar() override;
    CoolPropDbl calc_cpmolar() override;
    CoolPropDbl calc_gibbsmolar() override;
    CoolPropDbl calc_gibbsmolar_residual() override {
        return gas_constant() * _T * (alphar() + delta() * dalphar_dDelta());
    }
    CoolPropDbl calc_gibbsmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);

    CoolPropDbl calc_helmholtzmolar() override;
    CoolPropDbl calc_cpmolar_idealgas() override;
    CoolPropDbl calc_pressure_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_smolar() override;
    CoolPropDbl calc_smolar_residual() override {
        return gas_constant() * (tau() * dalphar_dTau() - alphar());
    }
    CoolPropDbl calc_smolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);

    CoolPropDbl calc_hmolar() override;
    CoolPropDbl calc_hmolar_residual() override {
        return gas_constant() * _T * (tau() * dalphar_dTau() + delta() * dalphar_dDelta());
    }
    CoolPropDbl calc_hmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);

    CoolPropDbl calc_umolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_umolar() override;
    CoolPropDbl calc_speed_sound() override;

    void calc_excess_properties() override;

    /** \brief The phase identification parameter of Venkatarathnam et al., FPE, 2011
     *
     *
     */

    CoolPropDbl calc_phase_identification_parameter();
    CoolPropDbl calc_fugacity(std::size_t i) override;
    CoolPropDbl calc_fugacity_coefficient(std::size_t i) override;
    CoolPropDbl calc_chemical_potential(std::size_t i) override;

    /// Using this backend, calculate the flame hazard
    CoolPropDbl calc_flame_hazard() override {
        return components[0].environment.FH;
    };
    /// Using this backend, calculate the health hazard
    CoolPropDbl calc_health_hazard() override {
        return components[0].environment.HH;
    };
    /// Using this backend, calculate the physical hazard
    CoolPropDbl calc_physical_hazard() override {
        return components[0].environment.PH;
    };

    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    CoolPropDbl calc_alphar() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dDelta() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
    CoolPropDbl calc_dalphar_dTau() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta2() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dDelta_dTau() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d2alphar_dTau2() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta3() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta2_dTau() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dDelta_dTau2() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d3alphar_dTau3() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta\delta}\f$ (dimensionless)
    CoolPropDbl calc_d4alphar_dDelta4() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta\tau}\f$ (dimensionless)
    CoolPropDbl calc_d4alphar_dDelta3_dTau() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d4alphar_dDelta2_dTau2() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d4alphar_dDelta_dTau3() override;
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau\tau}\f$ (dimensionless)
    CoolPropDbl calc_d4alphar_dTau4() override;

    CoolPropDbl calc_alpha0() override;
    CoolPropDbl calc_dalpha0_dDelta() override;
    CoolPropDbl calc_dalpha0_dTau() override;
    CoolPropDbl calc_d2alpha0_dDelta2() override;
    CoolPropDbl calc_d2alpha0_dDelta_dTau() override;
    CoolPropDbl calc_d2alpha0_dTau2() override;
    CoolPropDbl calc_d3alpha0_dDelta3() override;
    CoolPropDbl calc_d3alpha0_dDelta2_dTau() override;
    CoolPropDbl calc_d3alpha0_dDelta_dTau2() override;
    CoolPropDbl calc_d3alpha0_dTau3() override;

    CoolPropDbl calc_surface_tension() override;
    CoolPropDbl calc_viscosity() override;
    CoolPropDbl calc_viscosity_dilute();
    CoolPropDbl calc_viscosity_background();
    CoolPropDbl calc_viscosity_background(CoolPropDbl eta_dilute, CoolPropDbl& initial_density, CoolPropDbl& residual);
    CoolPropDbl calc_conductivity() override;
    CoolPropDbl calc_conductivity_background();

    /**
     * \brief Calculate each of the contributions to the viscosity
     *
     * If the viscosity model is hardcoded or ECS is being used, there will only be one entry in critical and all others will be zero
     */
    void calc_viscosity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) override;
    /**
     * \brief Calculate each of the contributions to the conductivity
     *
     * If the conductivity model is hardcoded or ECS is being used, there will only be one entry in initial_density and all others will be zero
     */
    void calc_conductivity_contributions(CoolPropDbl& dilute, CoolPropDbl& initial_density, CoolPropDbl& residual, CoolPropDbl& critical) override;

    CoolPropDbl calc_saturated_liquid_keyed_output(parameters key) override;
    CoolPropDbl calc_saturated_vapor_keyed_output(parameters key) override;

    CoolPropDbl calc_Tmin() override;
    CoolPropDbl calc_Tmax() override;
    CoolPropDbl calc_pmax() override;
    CoolPropDbl calc_Ttriple() override;
    CoolPropDbl calc_p_triple() override;
    CoolPropDbl calc_pmax_sat();
    CoolPropDbl calc_Tmax_sat();
    void calc_Tmin_sat(CoolPropDbl& Tmin_satL, CoolPropDbl& Tmin_satV);
    void calc_pmin_sat(CoolPropDbl& pmin_satL, CoolPropDbl& pmin_satV);

    CoolPropDbl calc_T_critical() override;
    CoolPropDbl calc_p_critical() override;
    CoolPropDbl calc_rhomolar_critical() override;

    CoolPropDbl calc_T_reducing() override {
        return get_reducing_state().T;
    };
    CoolPropDbl calc_rhomolar_reducing() override {
        return get_reducing_state().rhomolar;
    };
    CoolPropDbl calc_p_reducing() override {
        return get_reducing_state().p;
    };

    // Calculate the phase identification parameter of Venkatarathnam et al, Fluid Phase Equilibria
    CoolPropDbl calc_PIP() override {
        return 2
               - rhomolar()
                   * (second_partial_deriv(iP, iDmolar, iT, iT, iDmolar) / first_partial_deriv(iP, iT, iDmolar)
                      - second_partial_deriv(iP, iDmolar, iT, iDmolar, iT) / first_partial_deriv(iP, iDmolar, iT));
    };

    std::string calc_name() override;
    std::vector<std::string> calc_fluid_names() override;

    void calc_all_alphar_deriv_cache(const std::vector<CoolPropDbl>& mole_fractions, const CoolPropDbl& tau, const CoolPropDbl& delta);
    virtual CoolPropDbl calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl>& mole_fractions,
                                                  const CoolPropDbl& tau, const CoolPropDbl& delta);

    /**
    \brief Take derivatives of the ideal-gas part of the Helmholtz energy, don't use any cached values, or store any cached values

    @param nTau How many derivatives with respect to \f$\tau\f$ to take
    @param nDelta How many derivatives with respect to \f$\delta\f$ to take
    @param mole_fractions Mole fractions
    @param tau Reciprocal reduced temperature where \f$\tau=T_r / T\f$
    @param delta Reduced density where \f$\delta = \rho / \rho_r \f$
    @param Tr Reducing temperature of the mixture [K]
    @param rhor Reducing molar density of the mixture [mol/m^3]

    \f[
    \alpha^0 = \displaystyle\sum_{i=1}^{N}x_i[\alpha^0_{oi}(\rho,T) + \ln x_i]
    \f]
    where in this case, we use the \f$\alpha^0\f$ for the given fluid, which uses the inputs \f$\tau_i\f$ and \f$\delta_i\f$, so we do the conversion between mixture and component reduced states with
    \f[
    \tau_i = \frac{T_{c,i}}{T} = \frac{\tau T_{c,i}}{T_r}
    \f]
    \f[
    \delta_i = \frac{\rho}{\rho_{c,i}} = \frac{\delta\rho_r}{\rho_{c,i}}
    \f]

    \sa Table B5, GERG 2008 from Kunz Wagner, JCED, 2012
    */
    CoolPropDbl calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl>& mole_fractions, const CoolPropDbl& tau,
                                          const CoolPropDbl& delta, const CoolPropDbl& Tr, const CoolPropDbl& rhor);

    HelmholtzDerivatives calc_all_alpha0_derivs_nocache(const std::vector<CoolPropDbl>& mole_fractions, const CoolPropDbl& tau,
                                                        const CoolPropDbl& delta, const CoolPropDbl& Tr, const CoolPropDbl& rhor);

    void calc_reducing_state() override;
    virtual SimpleState calc_reducing_state_nocache(const std::vector<CoolPropDbl>& mole_fractions);

    const CoolProp::SimpleState& get_reducing_state() override {
        calc_reducing_state();
        return _reducing;
    };

    void update_states() override;

    CoolPropDbl calc_compressibility_factor() override {
        return 1 + delta() * dalphar_dDelta();
    };

    void calc_phase_envelope(const std::string& type) override;

    CoolPropDbl SRK_covolume();

    // ***************************************************************
    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    // ***************************************************************
    void T_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value);
    void p_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value, bool& saturation_called);
    void DmolarP_phase_determination();

    // ***************************************************************
    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
    // ***************************************************************

    virtual CoolPropDbl solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess = -1);
    virtual CoolPropDbl solver_rho_Tp_SRK(CoolPropDbl T, CoolPropDbl p, phases phase);
    enum StationaryPointReturnFlag
    {
        ZERO_STATIONARY_POINTS,
        ONE_STATIONARY_POINT_FOUND,
        TWO_STATIONARY_POINTS_FOUND
    };
    virtual StationaryPointReturnFlag solver_dpdrho0_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomax, CoolPropDbl& light, CoolPropDbl& heavy);
    virtual CoolPropDbl solver_rho_Tp_global(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rhomax);
};

class CorrespondingStatesTerm
{
   public:
    /// Calculate all the derivatives that do not involve any composition derivatives
    virtual HelmholtzDerivatives all(HelmholtzEOSMixtureBackend& HEOS, double tau, double delta, const std::vector<CoolPropDbl>& x,
                                     bool cache_values = false) {
        HelmholtzDerivatives summer;
        std::size_t N = x.size();
        for (std::size_t i = 0; i < N; ++i) {
            HelmholtzDerivatives derivs = HEOS.components[i].EOS().alphar.all(tau, delta, cache_values);
            summer = summer + derivs * x[i];
        }
        return summer;
    }
    CoolPropDbl dalphar_dxi(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().baser(HEOS.tau(), HEOS.delta());
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().baser(HEOS.tau(), HEOS.delta()) - HEOS.components[N - 1].EOS().baser(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().dalphar_dTau(HEOS._tau, HEOS._delta);
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().dalphar_dTau(HEOS._tau, HEOS._delta) - HEOS.components[N - 1].EOS().dalphar_dTau(HEOS._tau, HEOS._delta);
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().dalphar_dDelta(HEOS.tau(), HEOS.delta());
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().dalphar_dDelta(HEOS.tau(), HEOS.delta())
                   - HEOS.components[N - 1].EOS().dalphar_dDelta(HEOS._tau, HEOS._delta);
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d3alphar_dxi_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d2alphar_dDelta2(HEOS.tau(), HEOS.delta());
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().d2alphar_dDelta2(HEOS.tau(), HEOS.delta())
                   - HEOS.components[N - 1].EOS().d2alphar_dDelta2(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d3alphar_dxi_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d2alphar_dTau2(HEOS.tau(), HEOS.delta());
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().d2alphar_dTau2(HEOS.tau(), HEOS.delta())
                   - HEOS.components[N - 1].EOS().d2alphar_dTau2(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d3alphar_dxi_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d2alphar_dDelta_dTau(HEOS.tau(), HEOS.delta());
        } else if (xN_flag == XN_DEPENDENT) {
            std::size_t N = x.size();
            if (i == N - 1) return 0;
            return HEOS.components[i].EOS().d2alphar_dDelta_dTau(HEOS.tau(), HEOS.delta())
                   - HEOS.components[N - 1].EOS().d2alphar_dDelta_dTau(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }

    CoolPropDbl d2alphardxidxj(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                               x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return 0;
        } else if (xN_flag == XN_DEPENDENT) {
            return 0;
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }

    CoolPropDbl d4alphar_dxi_dDelta3(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d3alphar_dDelta3(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d4alphar_dxi_dTau3(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d3alphar_dTau3(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d4alphar_dxi_dDelta_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d3alphar_dDelta_dTau2(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }
    CoolPropDbl d4alphar_dxi_dDelta2_dTau(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, x_N_dependency_flag xN_flag) {
        if (xN_flag == XN_INDEPENDENT) {
            return HEOS.components[i].EOS().d3alphar_dDelta2_dTau(HEOS.tau(), HEOS.delta());
        } else {
            throw ValueError(format("xN_flag is invalid"));
        }
    }

    CoolPropDbl d3alphardxidxjdxk(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j, std::size_t k,
                                  x_N_dependency_flag xN_flag) {
        return 0;
    }
    CoolPropDbl d3alphar_dxi_dxj_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                        x_N_dependency_flag xN_flag) {
        return 0;
    }
    CoolPropDbl d3alphar_dxi_dxj_dTau(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                      x_N_dependency_flag xN_flag) {
        return 0;
    }
    CoolPropDbl d4alphar_dxi_dxj_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                         x_N_dependency_flag xN_flag) {
        return 0;
    }
    CoolPropDbl d4alphar_dxi_dxj_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                             x_N_dependency_flag xN_flag) {
        return 0;
    }
    CoolPropDbl d4alphar_dxi_dxj_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::vector<CoolPropDbl>& x, std::size_t i, std::size_t j,
                                       x_N_dependency_flag xN_flag) {
        return 0;
    }
};

/// This class contains the two primary contributions to the residual Helmholtz energy - a corresponding states
/// contribution, sometimes (incorrectly) referred to as ideal mixing, and an excess term
///
/// It delegates the calls to the corresponding states and excess contributions
/// The entire class can be replaced with a derived class
class ResidualHelmholtz
{
   public:
    ExcessTerm Excess;
    CorrespondingStatesTerm CS;

    ResidualHelmholtz() {};
    ResidualHelmholtz(const ExcessTerm& E, const CorrespondingStatesTerm& C) : Excess(E), CS(C) {};
    virtual ~ResidualHelmholtz() = default;

    ResidualHelmholtz copy() {
        return ResidualHelmholtz(Excess.copy(), CS);
    }
    ResidualHelmholtz* copy_ptr() {
        return new ResidualHelmholtz(Excess.copy(), CS);
    }

    virtual HelmholtzDerivatives all(HelmholtzEOSMixtureBackend& HEOS, const std::vector<CoolPropDbl>& mole_fractions, double tau, double delta,
                                     bool cache_values = false) {
        HelmholtzDerivatives a = CS.all(HEOS, tau, delta, mole_fractions, cache_values) + Excess.all(tau, delta, mole_fractions, cache_values);
        a.delta_x_dalphar_ddelta = delta * a.dalphar_ddelta;
        a.tau_x_dalphar_dtau = tau * a.dalphar_dtau;

        a.delta2_x_d2alphar_ddelta2 = POW2(delta) * a.d2alphar_ddelta2;
        a.deltatau_x_d2alphar_ddelta_dtau = delta * tau * a.d2alphar_ddelta_dtau;
        a.tau2_x_d2alphar_dtau2 = POW2(tau) * a.d2alphar_dtau2;

        return a;
    }
    virtual CoolPropDbl dalphar_dxi(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.dalphar_dxi(HEOS, mole_fractions, i, xN_flag) + Excess.dalphar_dxi(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d2alphardxidxj(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d2alphardxidxj(HEOS, mole_fractions, i, j, xN_flag) + Excess.d2alphardxidxj(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d2alphar_dxi_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d2alphar_dxi_dTau(HEOS, mole_fractions, i, xN_flag) + Excess.d2alphar_dxi_dTau(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d2alphar_dxi_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d2alphar_dxi_dDelta(HEOS, mole_fractions, i, xN_flag) + Excess.d2alphar_dxi_dDelta(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d3alphar_dxi_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphar_dxi_dTau2(HEOS, mole_fractions, i, xN_flag) + Excess.d3alphar_dxi_dTau2(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphar_dxi_dDelta_dTau(HEOS, mole_fractions, i, xN_flag) + Excess.d3alphar_dxi_dDelta_dTau(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d3alphar_dxi_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphar_dxi_dDelta2(HEOS, mole_fractions, i, xN_flag) + Excess.d3alphar_dxi_dDelta2(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphar_dxi_dxj_dTau(HEOS, mole_fractions, i, j, xN_flag) + Excess.d3alphar_dxi_dxj_dTau(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d3alphar_dxi_dxj_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphar_dxi_dxj_dDelta(HEOS, mole_fractions, i, j, xN_flag) + Excess.d3alphar_dxi_dxj_dDelta(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d3alphardxidxjdxk(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                          x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d3alphardxidxjdxk(HEOS, mole_fractions, i, j, k, xN_flag) + Excess.d3alphardxidxjdxk(mole_fractions, i, j, k, xN_flag);
    }

    virtual CoolPropDbl d4alphar_dxi_dTau3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dTau3(HEOS, mole_fractions, i, xN_flag) + Excess.d4alphar_dxi_dTau3(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta2_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dDelta2_dTau(HEOS, mole_fractions, i, xN_flag) + Excess.d4alphar_dxi_dDelta2_dTau(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dDelta_dTau2(HEOS, mole_fractions, i, xN_flag) + Excess.d4alphar_dxi_dDelta_dTau2(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dDelta3(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dDelta3(HEOS, mole_fractions, i, xN_flag) + Excess.d4alphar_dxi_dDelta3(mole_fractions, i, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dTau2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dxj_dTau2(HEOS, mole_fractions, i, j, xN_flag) + Excess.d4alphar_dxi_dxj_dTau2(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dxj_dDelta_dTau(HEOS, mole_fractions, i, j, xN_flag)
               + Excess.d4alphar_dxi_dxj_dDelta_dTau(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dDelta2(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, x_N_dependency_flag xN_flag) {
        std::vector<CoolPropDbl>& mole_fractions = HEOS.get_mole_fractions_ref();
        return CS.d4alphar_dxi_dxj_dDelta2(HEOS, mole_fractions, i, j, xN_flag) + Excess.d4alphar_dxi_dxj_dDelta2(mole_fractions, i, j, xN_flag);
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dxk_dDelta(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                    x_N_dependency_flag xN_flag) {
        return 0;
    }
    virtual CoolPropDbl d4alphar_dxi_dxj_dxk_dTau(HelmholtzEOSMixtureBackend& HEOS, std::size_t i, std::size_t j, std::size_t k,
                                                  x_N_dependency_flag xN_flag) {
        return 0;
    }
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
