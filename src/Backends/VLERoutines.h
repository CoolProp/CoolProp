
#ifndef VLEROUTINES_H
#define VLEROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

namespace SaturationSolvers
{
    struct saturation_T_pure_Akasaka_options{
        bool use_guesses; ///< true to start off at the values specified by rhoL, rhoV
        long double omega, rhoL, rhoV, pL, pV;
    };
    struct saturation_T_pure_options{
        bool use_guesses; ///< true to start off at the values specified by rhoL, rhoV
        long double omega, rhoL, rhoV, pL, pV;
    };
    
    struct saturation_D_pure_options{
        enum imposed_rho_options{IMPOSED_RHOL, IMPOSED_RHOV};
        bool use_guesses, ///< True to start off at the values specified by rhoL, rhoV, T
             use_logdelta; ///< True to use partials with respect to log(delta) rather than delta
        long double omega, rhoL, rhoV, pL, pV;
        int imposed_rho;
        saturation_D_pure_options(){ use_logdelta = true; omega = 1.0;} // Defaults
    };

    enum sstype_enum {imposed_T, imposed_p};
    struct mixture_VLE_IO
    {
        int sstype, Nstep_max;
        long double rhomolar_liq, rhomolar_vap, p, T, beta;
        std::vector<long double> x,y,K;
    };

    /*! Returns the natural logarithm of K for component i using the method from Wilson as in
	\f[
	\ln K_i = \ln\left(\frac{p_{c,i}}{p}\right)+5.373(1+\omega_i)\left(1-\frac{T_{c,i}}{T}\right)
	\f]
    @param HEOS The Helmholtz EOS mixture backend
	@param T Temperature [K]
	@param p Pressure [Pa]
	@param i Index of component [-]
	*/
	static long double Wilson_lnK_factor(HelmholtzEOSMixtureBackend *HEOS, long double T, long double p, int i){ 
        EquationOfState *EOS = (HEOS->get_components())[i]->pEOS; 
        return log(EOS->reduce.p/p)+5.373*(1 + EOS->accentric)*(1-EOS->reduce.T/T);
    };

    void saturation_D_pure(HelmholtzEOSMixtureBackend *HEOS, long double rhomolar, saturation_D_pure_options &options);
    void saturation_T_pure(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_options &options);
    void saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_Akasaka_options &options);
    
    /**
    */
    struct saturation_PHSU_pure_options{
        enum specified_variable_options{IMPOSED_HL, IMPOSED_HV, IMPOSED_PL, IMPOSED_PV, IMPOSED_SL, IMPOSED_SV, IMPOSED_UL, IMPOSED_UV, IMPOSED_INVALID_INPUT};
        bool use_guesses, ///< True to start off at the values specified by rhoL, rhoV, T
             use_logdelta; ///< True to use partials with respect to log(delta) rather than delta
        int specified_variable;
        long double omega, rhoL, rhoV, pL, pV;
        saturation_PHSU_pure_options(){ specified_variable = IMPOSED_INVALID_INPUT; use_guesses = true; omega = 1.0; }
    };
    /**

    */
    void saturation_PHSU_pure(HelmholtzEOSMixtureBackend *HEOS, long double specified_value, saturation_PHSU_pure_options &options);

    long double successive_substitution(HelmholtzEOSMixtureBackend *HEOS, const long double beta, long double T, long double p, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &options);
    void x_and_y_from_K(long double beta, const std::vector<long double> &K, const std::vector<long double> &z, std::vector<long double> &x, std::vector<long double> &y);

    /*! A wrapper function around the residual to find the initial guess for the bubble point temperature
    \f[
    r = \sum_i \frac{z_i(K_i-1)}{1-beta+beta*K_i}
    \f]
    */
    class WilsonK_resid : public FuncWrapper1D
    {
    public:
        int input_type;
	    double T, p, beta;
	    const std::vector<long double> *z;
        std::vector<long double> *K;
	    HelmholtzEOSMixtureBackend *HEOS;

	    WilsonK_resid(HelmholtzEOSMixtureBackend *HEOS, double beta, double imposed_value, int input_type, const std::vector<long double> &z, std::vector<long double> &K){ 
            this->z = &z; this->K = &K; this->HEOS = HEOS; this->beta = beta; this->input_type = input_type;
            if (input_type == imposed_T){
                this->T = imposed_value;
            }
            else{
                this->p = imposed_value;
            }
	    };
	    double call(double input_value){
		    double summer = 0;
            if (input_type == imposed_T){
                p = input_value; // Iterate on pressure
            }
            else{
                T = input_value; // Iterate on temperature, pressure imposed
            }
		    for (unsigned int i = 0; i< (*z).size(); i++) {
			    (*K)[i] = exp(Wilson_lnK_factor(HEOS,T,p,i));
			    summer += (*z)[i]*((*K)[i]-1)/(1-beta+beta*(*K)[i]);
		    }
		    return summer;
	    };
    };
    inline double saturation_preconditioner(HelmholtzEOSMixtureBackend *HEOS, double input_value, int input_type, const std::vector<long double> &z)
    {
	    double ptriple = 0, pcrit = 0, Ttriple = 0, Tcrit = 0;
	    
	    for (unsigned int i = 0; i < z.size(); i++)
	    {
            EquationOfState *EOS = (HEOS->get_components())[i]->pEOS; 

		    ptriple += EOS->ptriple*z[i];
            pcrit += EOS->reduce.p*z[i];
		    Ttriple += EOS->Ttriple*z[i];
		    Tcrit += EOS->reduce.T*z[i];
	    }

        if (input_type == imposed_T)
        {
            return exp(log(pcrit/ptriple)/(Tcrit-Ttriple)*(input_value-Ttriple)+log(ptriple));
        }
        else if (input_type == imposed_p)
        {
            return 1/(1/Tcrit-(1/Ttriple-1/Tcrit)/log(pcrit/ptriple)*log(input_value/pcrit));
        }
        else{ throw ValueError();}
    }
    inline double saturation_Wilson(HelmholtzEOSMixtureBackend *HEOS, double beta, double input_value, int input_type, const std::vector<long double> &z, double guess)
    {
	    double T;

	    std::string errstr;

	    // Find first guess for T using Wilson K-factors
        WilsonK_resid Resid(HEOS, beta, input_value, input_type, z, HEOS->get_K());
	    T = Secant(Resid, guess, 0.001, 1e-10, 100, errstr);
	
	    if (!ValidNumber(T)){throw ValueError("saturation_p_Wilson failed to get good T");}
	    return T;
    }
    struct SuccessiveSubstitutionStep
    {
        long double T,p;
    };

    /*!
    A class to do newton raphson solver for VLE given guess values for vapor-liquid equilibria.  This class will then be included in the Mixture class

    A class is used rather than a function so that it is easier to store iteration histories, additional output values, etc.
    */
    class newton_raphson_VLE_GV
    {
    public:
	    long double error_rms, rhobar_liq, rhobar_vap, T, p, max_rel_change;
	    unsigned int N;
	    bool logging;
	    int Nsteps;
	    STLMatrix J;
        HelmholtzEOSMixtureBackend *SatL, *SatV;
	    std::vector<long double> K, x, y, phi_ij_liq, phi_ij_vap, dlnphi_drho_liq, dlnphi_drho_vap, r, negative_r, dXdS, neg_dFdS;
	    std::vector<SuccessiveSubstitutionStep> step_logger;

	    newton_raphson_VLE_GV(){};

	    void resize(unsigned int N);
	
	    // Reset the state of all the internal variables
	    void pre_call()
	    {
		    K.clear(); x.clear(); y.clear();  phi_ij_liq.clear(); 
            phi_ij_vap.clear(); dlnphi_drho_liq.clear(), dlnphi_drho_vap.clear(),
            step_logger.clear(); error_rms = 1e99; Nsteps = 0;
		    rhobar_liq = _HUGE; rhobar_vap = _HUGE; T = _HUGE; p = _HUGE;
	    };

	    /*! Call the Newton-Raphson VLE Solver

	    This solver must be passed reasonable guess values for the mole fractions, 
	    densities, etc.  You may want to take a few steps of successive substitution
	    before you start with Newton Raphson.

	    @param HEOS Temperature [K]
	    @param z Pressure [Pa]
	    @param z Bulk mole fractions [-]
	    @param K Array of K-factors [-]
	    */
	    void call(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO);

	    /*! Build the arrays for the Newton-Raphson solve

	    This method builds the Jacobian matrix, the sensitivity matrix, etc.

	    @param beta Void fraction [-] (0: bubble, 1: dew)
	    @param T Temperature [K]
	    @param p Pressure [Pa]
	    @param z Bulk mole fractions [-]
	    @param K Array of K-factors [-]
	    */
	    void build_arrays(HelmholtzEOSMixtureBackend *HEOS, long double beta, long double T, long double rhomolar_liq, const long double rho_vapor, const std::vector<long double> &z, std::vector<long double> & K);

        /** Check the derivatives in the Jacobian using numerical derivatives.
        */
        void check_Jacobian(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO);
    };
};

namespace PhaseEnvelope
{
    class PhaseEnvelopeLog
    {
    public:
	    std::vector< std::vector<long double> > K, lnK, x, y;
	    std::vector<long double> T, p, lnT, lnp, rhomolar_liq, rhomolar_vap, lnrhomolar_liq, lnrhomolar_vap;
        void resize(std::size_t N)
        {
            K.resize(N);
            lnK.resize(N);
            x.resize(N);
            y.resize(N);
        }
	    void store_variables(const long double T, 
		                     const long double p, 
						     const long double rhomolar_liq, 
						     const long double rhomolar_vap, 
						     const std::vector<long double> & K,
						     const std::vector<long double> & x, 
						     const std::vector<long double> & y)
	    {
            std::size_t N = K.size();
		    this->p.push_back(p);
		    this->T.push_back(T);
		    this->lnT.push_back(log(T));
		    this->lnp.push_back(log(p));
		    this->rhomolar_liq.push_back(rhomolar_liq);
		    this->rhomolar_vap.push_back(rhomolar_vap);
		    this->lnrhomolar_liq.push_back(log(rhomolar_liq));
		    this->lnrhomolar_vap.push_back(log(rhomolar_vap));
		    for (unsigned int i = 0; i < N; i++)
		    {
			    this->K[i].push_back(K[i]);
			    this->lnK[i].push_back(log(K[i]));
                this->x[i].push_back(x[i]);
			    this->y[i].push_back(y[i]);
		    }
	    };
    };
    class PhaseEnvelope_GV
    {
    public:
	    PhaseEnvelopeLog bubble, dew;

        void build(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, SaturationSolvers::mixture_VLE_IO &IO);
    };
};

} /* namespace CoolProp*/

#endif