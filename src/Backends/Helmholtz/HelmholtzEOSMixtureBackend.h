
#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "AbstractState.h"
#include "CoolPropFluid.h"
#include "ReducingFunctions.h"
#include "ExcessHEFunction.h"
#include "Solvers.h"
#include "PhaseEnvelope.h"

#include <vector>

namespace CoolProp {

class FlashRoutines;

class HelmholtzEOSMixtureBackend : public AbstractState {
    
private:
    void pre_update(CoolProp::input_pairs &input_pair, CoolPropDbl &value1, CoolPropDbl &value2 );
    void post_update();
protected:
    std::vector<CoolPropFluid*> components; ///< The components that are in use

    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<CoolPropDbl> mole_fractions; ///< The bulk mole fractions of the mixture
    std::vector<CoolPropDbl> K, ///< The K factors for the components
                             lnK; ///< The natural logarithms of the K factors of the components

    SimpleState _crit;
    phases imposed_phase_index;
    std::size_t N; ///< Number of components
    
public:
    HelmholtzEOSMixtureBackend(){
        imposed_phase_index = iphase_not_imposed; _phase = iphase_unknown;};
    HelmholtzEOSMixtureBackend(const std::vector<CoolPropFluid*> &components, bool generate_SatL_and_SatV = true);
    HelmholtzEOSMixtureBackend(const std::vector<std::string> &component_names, bool generate_SatL_and_SatV = true);
    virtual ~HelmholtzEOSMixtureBackend(){};
    shared_ptr<ReducingFunction> Reducing;
    ExcessTerm Excess;
    PhaseEnvelopeData PhaseEnvelope;
    SimpleState hsat_max;
    SsatSimpleState ssat_max;

    friend class FlashRoutines; // Allows the static methods in the FlashRoutines class to have access to all the protected members and methods of this class
    friend class TransportRoutines; // Allows the static methods in the TransportRoutines class to have access to all the protected members and methods of this class
    friend class MixtureDerivatives; // Allows the static methods in the MixtureDerivatives class to have access to all the protected members and methods of this class
    friend class PhaseEnvelopeRoutines; // Allows the static methods in the PhaseEnvelopeRoutines class to have access to all the protected members and methods of this class
    friend class MixtureParameters; // Allows the static methods in the MixtureParameters class to have access to all the protected members and methods of this class

    // Helmholtz EOS backend uses mole fractions
    bool using_mole_fractions(){return true;}
    bool using_mass_fractions(){return false;}
    bool using_volu_fractions(){return false;}
    bool is_pure(){ return components.size() == 1 && !components[0]->EOSVector[0].pseudo_pure; }
    bool has_melting_line(){ return is_pure_or_pseudopure && components[0]->ancillaries.melting_line.enabled();};
    CoolPropDbl calc_melting_line(int param, int given, CoolPropDbl value);
    phases calc_phase(void){return _phase;};
    void calc_specify_phase(phases phase){ specify_phase(phase); }
    void calc_unspecify_phase(){ unspecify_phase(); }
    CoolPropDbl calc_saturation_ancillary(parameters param, int Q, parameters given, double value);
    void calc_ssat_max(void);
    void calc_hsat_max(void);
    CoolPropDbl calc_GWP20();
    CoolPropDbl calc_GWP500();
    CoolPropDbl calc_GWP100();
    CoolPropDbl calc_ODP();
	
	CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1);
    CoolPropDbl calc_first_saturation_deriv(parameters Of1, parameters Wrt1, HelmholtzEOSMixtureBackend &SatL, HelmholtzEOSMixtureBackend &SatV);
	CoolPropDbl calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Of2, parameters Wrt2);
    CoolPropDbl calc_first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant);
    CoolPropDbl calc_second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2);
    CoolPropDbl calc_first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, CoolPropDbl x_end);
    
	/// Calculate the phase once the state is fully calculated but you aren't sure if it is liquid or gas or ...
	void recalculate_singlephase_phase();

    const CoolProp::SimpleState &calc_state(const std::string &state);

    const std::vector<CoolPropFluid*> &get_components(){return components;};
    std::vector<CoolPropDbl> &get_K(){return K;};
    std::vector<CoolPropDbl> &get_lnK(){return lnK;};
    HelmholtzEOSMixtureBackend &get_SatL(){return *SatL;};
    HelmholtzEOSMixtureBackend &get_SatV(){return *SatV;};
    
    std::vector<CoolPropDbl> calc_mole_fractions_liquid(void){return SatL->get_mole_fractions();};
    std::vector<CoolPropDbl> calc_mole_fractions_vapor(void){return SatV->get_mole_fractions();};
    
    const CoolProp::PhaseEnvelopeData &calc_phase_envelope_data(){return PhaseEnvelope;};

    void resize(std::size_t N);
    shared_ptr<HelmholtzEOSMixtureBackend> SatL, SatV; ///<

    /** \brief The standard update function
     * @param input_pair The pair of inputs that will be provided
     * @param value1 The first input value
     * @param value2 The second input value
     */
    void update(CoolProp::input_pairs input_pair, double value1, double value2);
    
    /** \brief Update all the internal variables for a state by copying from another state
     */
    void update_internal(HelmholtzEOSMixtureBackend &HEOS);

    /** \brief Update with TP and a guess for rho
     * @param T Temperature in K
     * @param p Pressure in Pa
     * @param rho_guess Density in mol/m^3 guessed
     */
    void update_TP_guessrho(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess);
    void update_DmolarT_direct(CoolPropDbl rhomolar, CoolPropDbl T);
    void update_HmolarQ_with_guessT(CoolPropDbl hmolar, CoolPropDbl Q, CoolPropDbl Tguess);

    /** \brief Set the components of the mixture
     * 
     * @param components The components that are to be used in this mixture
     * @param generate_SatL_and_SatV true if SatL and SatV classes should be added, false otherwise.  Added so that saturation classes can be added without infinite recursion of adding saturation classes
     */
    void set_components(const std::vector<CoolPropFluid*> &components, bool generate_SatL_and_SatV = true);

    /** \brief Specify the phase - this phase will always be used in calculations
     * 
     * @param phase_index The index from CoolProp::phases
     */
    void specify_phase(phases phase_index){imposed_phase_index = phase_index; _phase = phase_index;};
    
    /**\brief Unspecify the phase - the phase is no longer imposed, different solvers can do as they like
     */
    void unspecify_phase(){imposed_phase_index = iphase_not_imposed;};

    /** \brief Set the mixture parameters - binary pair reducing functions, departure functions, F_ij, etc.
     */
    void set_mixture_parameters();

    /** \brief Set the mole fractions
     * 
     * @param mole_fractions The vector of mole fractions of the components
     */
    void set_mole_fractions(const std::vector<CoolPropDbl> &mole_fractions);

    const std::vector<CoolPropDbl> &get_mole_fractions(){return mole_fractions;};
    std::vector<CoolPropDbl> &get_mole_fractions_ref(){return mole_fractions;};

    /** \brief Set the mass fractions
     * 
     * @param mass_fractions The vector of mass fractions of the components
     */
    void set_mass_fractions(const std::vector<CoolPropDbl> &mass_fractions){throw std::exception();};

    CoolPropDbl calc_molar_mass(void);
    CoolPropDbl calc_gas_constant(void);
    CoolPropDbl calc_acentric_factor(void);

    CoolPropDbl calc_Bvirial(void);
    CoolPropDbl calc_Cvirial(void);
    CoolPropDbl calc_dBvirial_dT(void);
    CoolPropDbl calc_dCvirial_dT(void);

    CoolPropDbl calc_pressure(void);
    CoolPropDbl calc_cvmolar(void);
    CoolPropDbl calc_cpmolar(void);
    CoolPropDbl calc_gibbsmolar(void);
    CoolPropDbl calc_cpmolar_idealgas(void);
    CoolPropDbl calc_pressure_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_smolar(void);
    CoolPropDbl calc_smolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    
    CoolPropDbl calc_hmolar(void);
    CoolPropDbl calc_hmolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    
    CoolPropDbl calc_umolar_nocache(CoolPropDbl T, CoolPropDbl rhomolar);
    CoolPropDbl calc_umolar(void);
    CoolPropDbl calc_speed_sound(void);
    /** \brief The phase identification parameter of Venkatarathnam et al., FPE, 2011
     * 
     * 
     */
    CoolPropDbl calc_fugacity_coefficient(int i);
    CoolPropDbl calc_phase_identification_parameter(void);

	/// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    CoolPropDbl calc_alphar(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
	CoolPropDbl calc_dalphar_dDelta(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
	CoolPropDbl calc_dalphar_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dDelta2(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dDelta_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d2alphar_dTau2(void);    
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dDelta3(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dDelta2_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dDelta_dTau2(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
	CoolPropDbl calc_d3alphar_dTau3(void);

    CoolPropDbl calc_alpha0(void);
    CoolPropDbl calc_dalpha0_dDelta(void);
    CoolPropDbl calc_dalpha0_dTau(void);
    CoolPropDbl calc_d2alpha0_dDelta2(void);
    CoolPropDbl calc_d2alpha0_dDelta_dTau(void);
    CoolPropDbl calc_d2alpha0_dTau2(void);
    CoolPropDbl calc_d3alpha0_dDelta3(void);
    CoolPropDbl calc_d3alpha0_dDelta2_dTau(void);
    CoolPropDbl calc_d3alpha0_dDelta_dTau2(void);
    CoolPropDbl calc_d3alpha0_dTau3(void);

    CoolPropDbl calc_surface_tension(void);
    CoolPropDbl calc_viscosity(void);
    CoolPropDbl calc_viscosity_dilute(void);
    CoolPropDbl calc_viscosity_background(void);
    CoolPropDbl calc_viscosity_background(CoolPropDbl eta_dilute);
    CoolPropDbl calc_conductivity(void);
    CoolPropDbl calc_conductivity_background(void);

    CoolPropDbl calc_saturated_liquid_keyed_output(parameters key) { return SatL->keyed_output(key); }
    CoolPropDbl calc_saturated_vapor_keyed_output(parameters key) { return SatV->keyed_output(key); }

    CoolPropDbl calc_Tmin(void);
    CoolPropDbl calc_Tmax(void);
    CoolPropDbl calc_pmax(void);
    CoolPropDbl calc_Ttriple(void);
    CoolPropDbl calc_p_triple(void);
    CoolPropDbl calc_pmax_sat(void);
    CoolPropDbl calc_Tmax_sat(void);
    void calc_Tmin_sat(CoolPropDbl &Tmin_satL, CoolPropDbl &Tmin_satV);
    void calc_pmin_sat(CoolPropDbl &pmin_satL, CoolPropDbl &pmin_satV);
	
    CoolPropDbl calc_T_critical(void);
    CoolPropDbl calc_p_critical(void);
    CoolPropDbl calc_rhomolar_critical(void);
	
	CoolPropDbl calc_T_reducing(void){return get_reducing_state().T;};
    CoolPropDbl calc_rhomolar_reducing(void){return get_reducing_state().rhomolar;};
    CoolPropDbl calc_p_reducing(void){return get_reducing_state().p;};

    std::string calc_name(void);
	std::vector<std::string> calc_fluid_names(void);

    void calc_all_alphar_deriv_cache(const std::vector<CoolPropDbl> &mole_fractions, const CoolPropDbl &tau, const CoolPropDbl &delta);
    CoolPropDbl calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl> & mole_fractions, const CoolPropDbl &tau, const CoolPropDbl &delta);

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
    CoolPropDbl calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl> & mole_fractions, const CoolPropDbl &tau, const CoolPropDbl &delta, const CoolPropDbl &Tr, const CoolPropDbl &rhor);

    void calc_reducing_state(void);
    SimpleState calc_reducing_state_nocache(const std::vector<CoolPropDbl> & mole_fractions);

    const CoolProp::SimpleState & get_reducing_state(){calc_reducing_state(); return _reducing;};

    void update_states();
    
    CoolPropDbl calc_compressibility_factor(void){return 1+delta()*dalphar_dDelta();};
    
    void calc_phase_envelope(const std::string &type);
    
    void mass_to_molar_inputs(CoolProp::input_pairs &input_pair, CoolPropDbl &value1, CoolPropDbl &value2);

    // ***************************************************************
    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    // ***************************************************************
    void T_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value);
    void p_phase_determination_pure_or_pseudopure(int other, CoolPropDbl value, bool &saturation_called);
    void DmolarP_phase_determination();
    

    // ***************************************************************
    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
    // ***************************************************************

    CoolPropDbl solver_rho_Tp(CoolPropDbl T, CoolPropDbl p, CoolPropDbl rho_guess = -1);
    CoolPropDbl solver_rho_Tp_SRK(CoolPropDbl T, CoolPropDbl p, int phase);
    CoolPropDbl solver_for_rho_given_T_oneof_HSU(CoolPropDbl T, CoolPropDbl value, int other);

};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
