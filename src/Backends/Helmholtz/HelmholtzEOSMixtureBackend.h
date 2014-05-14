
#ifndef HELMHOLTZEOSMIXTUREBACKEND_H_
#define HELMHOLTZEOSMIXTUREBACKEND_H_

#include "AbstractState.h"
#include "../Fluids/CoolPropFluid.h"
#include "ReducingFunctions.h"
#include "ExcessHEFunction.h"
#include "Solvers.h"

#include <vector>

namespace CoolProp {

class FlashRoutines;

class HelmholtzEOSMixtureBackend : public AbstractState {
protected:
    std::vector<CoolPropFluid*> components; ///< The components that are in use
    
    bool is_pure_or_pseudopure; ///< A flag for whether the substance is a pure or pseudo-pure fluid (true) or a mixture (false)
    std::vector<long double> mole_fractions; ///< The mole fractions of the components
    std::vector<long double> mole_fractions_liq, ///< The mole fractions of the saturated liquid 
                             mole_fractions_vap, ///< The mole fractions of the saturated vapor
                             K, ///< The K factors for the components
                             lnK; ///< The natural logarithms of the K factors of the components

    SimpleState _crit;
    int imposed_phase_index;
public:
    HelmholtzEOSMixtureBackend(){SatL = NULL; SatV = NULL; imposed_phase_index = -1;};
    HelmholtzEOSMixtureBackend(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);
    HelmholtzEOSMixtureBackend(std::vector<std::string> &component_names, bool generate_SatL_and_SatV = true);
    virtual ~HelmholtzEOSMixtureBackend(){};
    ReducingFunctionContainer Reducing;
    ExcessTerm Excess;

    friend class FlashRoutines; // Allows the routines in the FlashRoutines class to have access to all the protected members and methods of this class

    // Helmholtz EOS backend uses mole fractions
    bool using_mole_fractions(){return true;}

    const std::vector<CoolPropFluid*> &get_components(){return components;};
    std::vector<long double> &get_K(){return K;};
    std::vector<long double> &get_lnK(){return lnK;};

    HelmholtzEOSMixtureBackend *SatL, *SatV; ///< 

    void update(long input_pair, double value1, double value2);

    void update_TP_guessrho(long double T, long double p, long double rho_guess);

    /// Set the components of the mixture
    /**
    @param components The components that are to be used in this mixture
    @param generate_SatL_and_SatV true if SatL and SatV classes should be added, false otherwise.  Added so that saturation classes can be added without infinite recursion of adding saturation classes
    */
    void set_components(std::vector<CoolPropFluid*> components, bool generate_SatL_and_SatV = true);

    /**
    \brief Specify the phase - this phase will always be used in calculations
    @param phase_index The index from CoolProp::phases
    */
    void specify_phase(int phase_index){imposed_phase_index = phase_index;};

    void set_reducing_function();
    void set_excess_term();

    /// Set the mole fractions
    /** 
    @param mole_fractions The vector of mole fractions of the components
    */
    void set_mole_fractions(const std::vector<long double> &mole_fractions);

    const std::vector<long double> &get_mole_fractions(){return mole_fractions;};
    
    /// Set the mass fractions
    /** 
    @param mass_fractions The vector of mass fractions of the components
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
    long double calc_hmolar(void);
    long double calc_smolar(void);
    long double calc_pressure_nocache(long double T, long double rhomolar);
    long double calc_smolar_nocache(long double T, long double rhomolar);
    long double calc_hmolar_nocache(long double T, long double rhomolar);
    long double calc_umolar_nocache(long double T, long double rhomolar);
    long double calc_umolar(void);
    long double calc_speed_sound(void);
    long double calc_fugacity_coefficient(int i);

    long double calc_alphar(void);
    long double calc_dalphar_dDelta(void);
    long double calc_dalphar_dTau(void);
    long double calc_d2alphar_dDelta2(void);
    long double calc_d2alphar_dDelta_dTau(void);
    long double calc_d2alphar_dTau2(void);

    long double calc_alpha0(void);
    long double calc_dalpha0_dDelta(void);
    long double calc_dalpha0_dTau(void);
    long double calc_d2alpha0_dDelta2(void);
    long double calc_d2alpha0_dDelta_dTau(void);
    long double calc_d2alpha0_dTau2(void);

    long double calc_surface_tension(void);

    long double calc_Tmax(void);
    long double calc_pmax(void);
    long double calc_Ttriple(void);

    std::string calc_name(void);

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

    /**
    In a general way we can calculate any first partial derivative based on calculating derivatives with respect to the fundamental variables of the EOS, \f$\tau\f$ and \f$\delta\f$
    \f[
    \left(\frac{\partial A}{\partial B}\right)_C = \frac{\left(\frac{\partial A}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial A}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta}{\left(\frac{\partial B}{\partial \tau}\right)_\delta\left(\frac{\partial C}{\partial \delta}\right)_\tau-\left(\frac{\partial B}{\partial \delta}\right)_\tau\left(\frac{\partial C}{\partial \tau}\right)_\delta}
    \f]
    */
    long double calc_first_partial_deriv(int Of, int Wrt, int Constant);

    void mass_to_molar_inputs(long &input_pair, double &value1, double &value2);

    // *************************************************************** 
    // ***************************************************************
    // *************  PHASE DETERMINATION ROUTINES  ******************
    // ***************************************************************
    // *************************************************************** 
    void T_phase_determination_pure_or_pseudopure(int other, long double value);
    void p_phase_determination_pure_or_pseudopure(int other, long double value);
    void DmolarP_phase_determination();


    // ***************************************************************
    // ***************************************************************
    // *******************  SOLVER ROUTINES  *************************
    // ***************************************************************
    // ***************************************************************        
    
    long double solver_rho_Tp(long double T, long double p, long double rho_guess = -1);
    long double solver_rho_Tp_SRK(long double T, long double p, int phase);
    long double solver_for_rho_given_T_oneof_HSU(long double T, long double value, int other);








    // ***************************************************************
    // ***************************************************************
    // *****************  MIXTURE DERIVATIVES  ***********************
    // ***************************************************************
    // ***************************************************************


    long double mixderiv_dalphar_dxi(int i);
    long double mixderiv_d2alphar_dxi_dTau(int i);
    long double mixderiv_d2alphar_dxi_dDelta(int i);
    long double mixderiv_d2alphardxidxj(int i, int j);

    /*! The derivative term
	\f[
	\left(\frac{\partial p}{\partial T} \right)_{V,\bar n} = \rho R(1+\delta \alpha_{\delta}^r-\delta \tau \alpha^r_{\delta\tau})
	\f]
	GERG 2004 Monograph equation 7.61
	*/
	long double mixderiv_dpdT__constV_n();

    long double mixderiv_dpdrho__constT_n();

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial V} \right)_{T,\bar n} = -\rho^2 RT(1+2\delta \alpha_{\delta}^r+\delta^2\alpha^r_{\delta\delta})
	\f]
	GERG 2004 Monograph equation 7.62
	*/
	long double mixderiv_ndpdV__constT_n();

	/*! The derivative term
	\f[
	n\left(\frac{\partial p}{\partial n_i} \right)_{T,V,n_j} = \rho RT\left[1+\delta\alpha_{\delta}^r\left[2- \frac{1}{\rho_r}\cdot n\left( \frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] +\delta\cdot n\left(\frac{\partial\alpha_{\delta}^r}{\partial n_i}\right)_{T,V,n_j}\right]
	\f]
	GERG 2004 Monograph equation 7.63
	*/
	long double mixderiv_ndpdni__constT_V_nj(int i);

    /// GERG 2004 monograph Eqn. 7.32
    /*! The partial molar volume
	\f[
	\hat v_i = \left( \frac{\partial V}{\partial n_i}\right)_{T,p,n_j} = \frac{-\left(\dfrac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\dfrac{\partial p}{\partial V}\right)_{T,\bar n}}
	\f]
	*/
	long double mixderiv_partial_molar_volume(int i);

    /*!
	Natural logarithm of the fugacity coefficient
	*/
    long double mixderiv_ln_fugacity_coefficient(int i);

	/*!
	Derivative of the natural logarithm of the fugacity coefficient with respect to T
	*/
	long double mixderiv_dln_fugacity_coefficient_dT__constrho_n(int i);

    /*!
	Derivative of the natural logarithm of the fugacity coefficient with respect to T
	*/
	long double mixderiv_dln_fugacity_coefficient_drho__constT_n(int i);

    /// GERG 2004 Monograph Eqn. 7.29
	/** The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial T} \right)_{p,\bar n} = \left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} + \frac{1}{T}-\frac{\hat v}{RT}\left(\frac{\partial p}{\partial T}\right)_{V,\bar n}
	\f]
	*/
	long double mixderiv_dln_fugacity_coefficient_dT__constp_n(int i);

    /// Table B4, Kunz, JCED, 2012 for the original term and the subsequent substitutions
    /*! The derivative term
	\f[
	n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j}
	\f]
	which is equal to
	\f{eqnarray*}{
	n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} &=& \delta \phi^r_{\delta}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
	&& +\tau \phi^r_{\tau}\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
	&& +\phi^r_{x_i}-\sum_{k=1}^{N}x_k\phi^r_{x_k}
	\f}
	*/
	long double mixderiv_ndalphar_dni__constT_V_nj(int i);

	/// GERG Equation 7.42
	long double mixderiv_dnalphar_dni__constT_V_nj(int i);

    /// GERG 2004 Monograph Eqn. 7.30
	/**
    The derivative term
	\f[
	\left(\frac{\partial \ln \phi_i}{\partial p} \right)_{T,\bar n} = \frac{\hat v_i}{RT}-\frac{1}{p}
	\f]
	*/
	long double mixderiv_dln_fugacity_coefficient_dp__constT_n(int i);

	///GERG 2004 Monograph Equation 7.31
    /** The derivative term
	\f[
	n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1+\frac{n}{RT}\frac{\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\left(\frac{\partial p}{\partial n_i}\right)_{T,V,n_j}}{\left(\frac{\partial p}{\partial V}\right)_{T,\bar n}}
	\f]
	which is also equal to 
	\f[
	n\left(\frac{\partial \ln \phi_i}{\partial n_j}\right)_{T,p} = n\left(\frac{\partial^2n\alpha^r}{\partial n_j \partial n_i} \right)_{T,V}+1-\frac{\hat v_i}{RT}\left[n\left(\frac{\partial p}{\partial n_j}\right)_{T,V,n_i}\right]
	\f]
	*/
	long double mixderiv_ndln_fugacity_coefficient_dnj__constT_p(int i, int j);

	/// Gernert Equation 3.115
	/// Catch test provided
	long double mixderiv_dln_fugacity_coefficient_dxj__constT_p_xi(int i, int j);
	
	/// Gernert Equation 3.130
	/// Catch test provided
	long double mixderiv_dpdxj__constT_V_xi(int j);

	/// Gernert Equation 3.117
    long double mixderiv_d2nalphar_dni_dxj__constT_V(int i, int j){ return mixderiv_d_ndalphardni_dxj__constT_V_xi(i, j) + mixderiv_dalphar_dxj__constT_V_xi(j);};

	/// Gernert Equation 3.119
	/// Catch test provided
	long double mixderiv_dalphar_dxj__constT_V_xi(int j);

	/// Gernert Equation 3.118
	/// Catch test provided
	long double mixderiv_d_ndalphardni_dxj__constT_V_xi(int i, int j);

	/// Gernert Equation 3.134
	/// Catch test provided
	long double mixderiv_d_dalpharddelta_dxj__constT_V_xi(int j);

	/// Gernert Equation 3.121
	/// Catch test provided
	long double mixderiv_ddelta_dxj__constT_V_xi(int j);

	/// Gernert Equation 3.122
	/// Catch test provided
	long double mixderiv_dtau_dxj__constT_V_xi(int j);

	///  GERG 2004 Monograph, equations 7.44 and 7.51
    /** The derivative term
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = \left( \frac{\partial}{\partial T}\left(\frac{\partial n \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right)_{V,\bar n}
	\f]
	\f[
	\left(\frac{\partial^2n\alpha^r}{\partial T\partial n_i} \right)_{V,n_j} = -\frac{\tau}{T}\left[\alpha_{\tau}^r +\left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\bar x}\right]
	\f]
	*/
	long double mixderiv_d2nalphar_dni_dT(int i);

	/// GERG 2004 Monograph Equation 7.51 and Table B4, Kunz, JCED, 2012
    /** The derivative term
	\f{eqnarray*}{
	\frac{\partial }{\partial \tau} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right) &=& \delta \phi^r_{\delta\tau}\left[ 1-\frac{1}{\rho_r}\left[\left(\frac{\partial \rho_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial \rho_r}{\partial x_k}\right)_{x_j}  \right]\right]\\
	&& +(\tau \phi^r_{\tau\tau}+\phi^r_{\tau})\frac{1}{T_r}\left[\left(\frac{\partial T_r}{\partial x_i}\right)_{x_j} - \sum_{k=1}^N x_k\left(\frac{\partial T_r}{\partial x_k}\right)_{x_j}  \right]\\
	&& +\phi^r_{x_i\tau}-\sum_{k=1}^{N}x_k\phi^r_{x_k\tau}
	\f}
	*/
	long double mixderiv_d_ndalphardni_dTau(int i);

	/// GERG 2004 Monograph Equation 7.50 and Table B4, Kunz, JCED, 2012
    /** The derivative term
	\f{eqnarray*}{
	\left(\frac{\partial }{\partial \delta} \left( n\left(\frac{\partial \phi^r}{\partial n_i} \right)_{T,V,n_j} \right)\right)_{\tau,\bar x} &=& (\alpha_{\delta}^r+\delta\alpha_{\delta\delta}^r)\left[1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j} \right] \\
	&+&\tau\alpha^r_{\delta\tau}\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
	&+&\phi^r_{\delta x_i}-\sum_{k=1}^{N}x_k\phi^r_{\delta x_k}
	\f}
	*/
	long double mixderiv_d_ndalphardni_dDelta(int i);

	/** GERG 2004 Monograph equation 7.41:
	\f[
	n\left(\frac{\partial^2n\alpha^r}{\partial n_i \partial n_j} \right)_{T,V} = n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i}
	\f]
	and
	GERG 2004 Monograph equation 7.46:
	\f[
	n\left( \frac{\partial}{\partial n_j}\left(\frac{\partial n\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)_{T,V,n_i} = n\left( \frac{\partial \alpha^r}{\partial n_j}\right)_{T,V,n_i} + n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i}
	\f]
	GERG 2004 Monograph equation 7.47:
	\f{eqnarray*}{
	n\left( \frac{\partial}{\partial n_j}\left(n\left(\frac{\partial \alpha^r}{\partial n_i}\right)_{T,V,n_j} \right) \right)_{T,V,n_i} &=& \left( \frac{\partial}{\partial \delta}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\delta}{\partial n_j}\right)_{T,V,n_i}\\ 
	&+& \left( \frac{\partial}{\partial \tau}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\tau,\bar x}\cdot n\left(\frac{\partial\tau}{\partial n_j}\right)_{T,V,n_i}\\
	&+& \left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}-\sum_{k=1}^{N}x_k \left( \frac{\partial}{\partial x_k}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i}\\
	\f}
	*/
	long double mixderiv_nd2nalphardnidnj__constT_V(int i, int j);

    /// GERG 2004 Monograph equation 7.48
	/** The derivative term
	\f[
	n\left(\frac{\partial \delta}{\partial n_i} \right)_{T,V,n_j} = \delta - \frac{\delta}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i} \right)_{n_j}
	\f]
	*/
	long double mixderiv_nddeltadni__constT_V_nj(int i);
	
    /// GERG 2004 Monograph equation 7.49
    /** The derivative term
	\f[
	n\left(\frac{\partial \tau}{\partial n_i} \right)_{T,V,n_j} = \frac{\tau}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i} \right)_{n_j}
	\f]
	*/
	long double mixderiv_ndtaudni__constT_V_nj(int i);

    /// \brief GERG 2004 Monograph equation 7.52
	/** The derivative term
	\f{eqnarray*}{
	\left( \frac{\partial}{\partial x_j}\left(n\left(\frac{\partial\alpha^r}{\partial n_i}\right)_{T,V,n_j}\right)\right)_{\delta,\tau,x_i} &=& \delta\alpha_{\delta x_j}^{r}\left[ 1-\frac{1}{\rho_r}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
	&-& \delta\alpha_{\delta}^{r}\frac{1}{\rho_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{\rho_r}\left(\frac{\partial \rho_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial \rho_r}{\partial n_i}\right)_{n_j}\right] \\
	&+& \tau\alpha_{\tau x_j}^r\frac{1}{T_r}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\\
	&+& \tau\alpha_{\tau}^{r}\frac{1}{T_r}\left[ \left(\frac{\partial}{\partial x_j}\left(n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right)\right)_{x_i}-\frac{1}{T_r}\left(\frac{\partial T_r}{\partial x_j}\right)_{x_i}\cdot n\left(\frac{\partial T_r}{\partial n_i}\right)_{n_j}\right] \\
	&+& \alpha_{x_ix_j}^r-\alpha_{x_j}^r-\sum_{m=1}^Nx_m\alpha_{x_jx_m}^r
	\f}
	*/
	long double mixderiv_d_ndalphardni_dxj__constdelta_tau_xi(int i, int j);
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSMIXTUREBACKEND_H_ */
