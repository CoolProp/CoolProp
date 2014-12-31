
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
    void pre_update(CoolProp::input_pairs &input_pair, long double &value1, long double &value2 );
    void post_update();
protected:
    std::vector<CoolPropFluid*> components; ///< The components that are in use

    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<long double> mole_fractions; ///< The bulk mole fractions of the mixture
    std::vector<long double> K, ///< The K factors for the components
                             lnK; ///< The natural logarithms of the K factors of the components

    SimpleState _crit;
    phases imposed_phase_index;
    std::size_t N; ///< Number of components
    
public:
    HelmholtzEOSMixtureBackend(){
        imposed_phase_index = iphase_not_imposed; _phase = iphase_unknown;};
    HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);
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
    bool is_pure(){ return is_pure_or_pseudopure; }
    bool has_melting_line(){ return is_pure_or_pseudopure && components[0]->ancillaries.melting_line.enabled();};
    long double calc_melting_line(int param, int given, long double value);
    phases calc_phase(void){return _phase;};
    void calc_specify_phase(phases phase){ specify_phase(phase); }
    void calc_unspecify_phase(){ unspecify_phase(); }
    long double calc_saturation_ancillary(parameters param, int Q, parameters given, double value);
    void calc_ssat_max(void);
    void calc_hsat_max(void);
    long double calc_GWP20();
    long double calc_GWP500();
    long double calc_GWP100();
    long double calc_ODP();
	
	long double calc_first_saturation_deriv(parameters Of1, parameters Wrt1);
	long double calc_second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Of2, parameters Wrt2);
    
	/// Calculate the phase once the state is fully calculated but you aren't sure if it is liquid or gas or ...
	void recalculate_singlephase_phase();

    const CoolProp::SimpleState &calc_state(const std::string &state);

    const std::vector<CoolPropFluid*> &get_components(){return components;};
    std::vector<long double> &get_K(){return K;};
    std::vector<long double> &get_lnK(){return lnK;};
    HelmholtzEOSMixtureBackend &get_SatL(){return *SatL;};
    HelmholtzEOSMixtureBackend &get_SatV(){return *SatV;};
    
    std::vector<long double> calc_mole_fractions_liquid(void){return SatL->get_mole_fractions();};
    std::vector<long double> calc_mole_fractions_vapor(void){return SatV->get_mole_fractions();};
    
    const CoolProp::PhaseEnvelopeData &calc_phase_envelope_data(){return PhaseEnvelope;};

    void resize(unsigned int N);
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
    void update_TP_guessrho(long double T, long double p, long double rho_guess);
    void update_DmolarT_direct(long double rhomolar, long double T);
    void update_HmolarQ_with_guessT(long double hmolar, long double Q, long double Tguess);

    /** \brief Set the components of the mixture
     * 
     * @param components The components that are to be used in this mixture
     * @param generate_SatL_and_SatV true if SatL and SatV classes should be added, false otherwise.  Added so that saturation classes can be added without infinite recursion of adding saturation classes
     */
    void set_components(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);

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
    void set_mole_fractions(const std::vector<long double> &mole_fractions);

    const std::vector<long double> &get_mole_fractions(){return mole_fractions;};
    std::vector<long double> &get_mole_fractions_ref(){return mole_fractions;};

    /** \brief Set the mass fractions
     * 
     * @param mass_fractions The vector of mass fractions of the components
     */
    void set_mass_fractions(const std::vector<long double> &mass_fractions){throw std::exception();};

    long double calc_molar_mass(void);
    long double calc_gas_constant(void);

    long double calc_Bvirial(void);
    long double calc_Cvirial(void);
    long double calc_dBvirial_dT(void);
    long double calc_dCvirial_dT(void);

    long double calc_pressure(void);
    long double calc_cvmolar(void);
    long double calc_cpmolar(void);
    long double calc_gibbsmolar(void);
    long double calc_cpmolar_idealgas(void);
    long double calc_pressure_nocache(long double T, long double rhomolar);
    long double calc_smolar(void);
    long double calc_smolar_nocache(long double T, long double rhomolar);
    
    long double calc_hmolar(void);
    long double calc_hmolar_nocache(long double T, long double rhomolar);
    
    long double calc_umolar_nocache(long double T, long double rhomolar);
    long double calc_umolar(void);
    long double calc_speed_sound(void);
    /** \brief The phase identification parameter of Venkatarathnam et al., FPE, 2011
     * 
     * 
     */
    long double calc_fugacity_coefficient(int i);
    long double calc_phase_identification_parameter(void);

	/// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r\f$ (dimensionless)
    long double calc_alphar(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta}\f$ (dimensionless)
	long double calc_dalphar_dDelta(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau}\f$ (dimensionless)
	long double calc_dalphar_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta}\f$ (dimensionless)
	long double calc_d2alphar_dDelta2(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau}\f$ (dimensionless)
	long double calc_d2alphar_dDelta_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau}\f$ (dimensionless)
	long double calc_d2alphar_dTau2(void);    
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\delta}\f$ (dimensionless)
	long double calc_d3alphar_dDelta3(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\delta\tau}\f$ (dimensionless)
	long double calc_d3alphar_dDelta2_dTau(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\delta\tau\tau}\f$ (dimensionless)
	long double calc_d3alphar_dDelta_dTau2(void);
    /// Using this backend, calculate the residual Helmholtz energy term \f$\alpha^r_{\tau\tau\tau}\f$ (dimensionless)
	long double calc_d3alphar_dTau3(void);

    long double calc_alpha0(void);
    long double calc_dalpha0_dDelta(void);
    long double calc_dalpha0_dTau(void);
    long double calc_d2alpha0_dDelta2(void);
    long double calc_d2alpha0_dDelta_dTau(void);
    long double calc_d2alpha0_dTau2(void);
    long double calc_d3alpha0_dDelta3(void);
    long double calc_d3alpha0_dDelta2_dTau(void);
    long double calc_d3alpha0_dDelta_dTau2(void);
    long double calc_d3alpha0_dTau3(void);

    long double calc_surface_tension(void);
    long double calc_viscosity(void);
    long double calc_viscosity_dilute(void);
    long double calc_viscosity_background(void);
    long double calc_viscosity_background(long double eta_dilute);
    long double calc_conductivity(void);
    long double calc_conductivity_background(void);

    long double calc_saturated_liquid_keyed_output(parameters key) { return SatL->keyed_output(key); }
    long double calc_saturated_vapor_keyed_output(parameters key) { return SatV->keyed_output(key); }

    long double calc_Tmin(void);
    long double calc_Tmax(void);
    long double calc_pmax(void);
    long double calc_Ttriple(void);
    long double calc_p_triple(void);
    long double calc_pmax_sat(void);
    long double calc_Tmax_sat(void);
    void calc_Tmin_sat(long double &Tmin_satL, long double &Tmin_satV);
    void calc_pmin_sat(long double &pmin_satL, long double &pmin_satV);
	
    long double calc_T_critical(void);
    long double calc_p_critical(void);
    long double calc_rhomolar_critical(void);
	
	long double calc_T_reducing(void){return get_reducing_state().T;};
    long double calc_rhomolar_reducing(void){return get_reducing_state().rhomolar;};

    std::string calc_name(void);
	std::vector<std::string> calc_fluid_names(void);

    void calc_all_alphar_deriv_cache(const std::vector<long double> &mole_fractions, const long double &tau, const long double &delta);
    long double calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<long double> & mole_fractions, const long double &tau, const long double &delta);

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
    long double calc_alpha0_deriv_nocache(const int nTau, const int nDelta, const std::vector<long double> & mole_fractions, const long double &tau, const long double &delta, const long double &Tr, const long double &rhor);

    void calc_reducing_state(void);
    SimpleState calc_reducing_state_nocache(const std::vector<long double> & mole_fractions);

    const CoolProp::SimpleState & get_reducing_state(){calc_reducing_state(); return _reducing;};

    void update_states();
    
    long double calc_compressibility_factor(void){return 1+delta()*dalphar_dDelta();};
    
    void calc_phase_envelope(const std::string &type);
    
    void mass_to_molar_inputs(CoolProp::input_pairs &input_pair, long double &value1, long double &value2);

    // ***************************************************************
    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    // ***************************************************************
    void T_phase_determination_pure_or_pseudopure(int other, long double value);
    void p_phase_determination_pure_or_pseudopure(int other, long double value, bool &saturation_called);
    void DmolarP_phase_determination();
    

    // ***************************************************************
    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
    // ***************************************************************

    long double solver_rho_Tp(long double T, long double p, long double rho_guess = -1);
    long double solver_rho_Tp_SRK(long double T, long double p, int phase);
    long double solver_for_rho_given_T_oneof_HSU(long double T, long double value, int other);

};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
