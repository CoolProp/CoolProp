
#ifndef VLEROUTINES_H
#define VLEROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

namespace SaturationSolvers
{
    struct saturation_T_pure_Akasaka_options{
        bool use_guesses; ///< true to start off at the values specified by rhoL, rhoV
        CoolPropDbl omega, rhoL, rhoV, pL, pV;
        saturation_T_pure_Akasaka_options(){omega = _HUGE; rhoV = _HUGE; rhoL = _HUGE; pV = _HUGE, pL = _HUGE;}
    };
    struct saturation_T_pure_options{
        bool use_guesses; ///< true to start off at the values specified by rhoL, rhoV
        CoolPropDbl omega, rhoL, rhoV, pL, pV, p, T;
        saturation_T_pure_options(){omega = _HUGE; rhoV = _HUGE; rhoL = _HUGE; rhoL = _HUGE; pV = _HUGE, pL = _HUGE; T = _HUGE;}
    };
    
    struct saturation_D_pure_options{
        enum imposed_rho_options{IMPOSED_RHOL, IMPOSED_RHOV};
        bool use_guesses, ///< True to start off at the values specified by rhoL, rhoV, T
             use_logdelta; ///< True to use partials with respect to log(delta) rather than delta
        CoolPropDbl omega, rhoL, rhoV, pL, pV;
        int imposed_rho;
        saturation_D_pure_options(){ use_logdelta = true; omega = 1.0;} // Defaults
    };

    enum sstype_enum {imposed_T, imposed_p};
    struct mixture_VLE_IO
    {
        int sstype, Nstep_max;
        CoolPropDbl rhomolar_liq, rhomolar_vap, p, T, beta;
        std::vector<CoolPropDbl> x, y, K;
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
    static CoolPropDbl Wilson_lnK_factor(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, CoolPropDbl p, std::size_t i){ 
        EquationOfState *EOS = (HEOS.get_components())[i]->pEOS; 
        return log(EOS->reduce.p/p)+5.373*(1 + EOS->acentric)*(1-EOS->reduce.T/T);
    };

    void saturation_D_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl rhomolar, saturation_D_pure_options &options);
    void saturation_T_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_options &options);
    void saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options &options);
    void saturation_T_pure_Maxwell(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options &options);
    
    /**
    */
    struct saturation_PHSU_pure_options{
        enum specified_variable_options{IMPOSED_HL, IMPOSED_HV, IMPOSED_PL, IMPOSED_PV, IMPOSED_SL, IMPOSED_SV, IMPOSED_UL, IMPOSED_UV, IMPOSED_INVALID_INPUT};
        bool use_guesses, ///< True to start off at the values specified by rhoL, rhoV, T
             use_logdelta; ///< True to use partials with respect to log(delta) rather than delta
        specified_variable_options specified_variable;
        CoolPropDbl omega, rhoL, rhoV, pL, pV, T, p;
        saturation_PHSU_pure_options(){ specified_variable = IMPOSED_INVALID_INPUT; use_guesses = true; omega = 1.0; }
    };
    /**

    */
    void saturation_PHSU_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl specified_value, saturation_PHSU_pure_options &options);

    /* \brief This is a backup saturation_p solver for the case where the Newton solver cannot approach closely enough the solution
     *
     * This is especially a problem at low pressures where catastrophic truncation error occurs, especially in the saturated vapor side
     * 
     * @param HEOS The Helmholtz EOS backend instance to be used
     * @param p Imposed pressure in kPa
     * @param options Options to be passed to the function (at least T, rhoL and rhoV must be provided)
     */
    void saturation_P_pure_1D_T(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl p, saturation_PHSU_pure_options &options);
    
    /* \brief This is a backup saturation_T solver for the case where the Newton solver cannot approach closely enough the solution
     *
     * This is especially a problem at low pressures where catastrophic truncation error occurs, especially in the saturated vapor side
     * 
     * @param HEOS The Helmholtz EOS backend instance to be used
     * @param T Imposed temperature in K
     * @param options Options to be passed to the function (at least p, rhoL and rhoV must be provided)
     */
    void saturation_T_pure_1D_P(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_options &options);

    /* \brief A robust but slow solver in the very-near-critical region
     * 
     * This solver operates in the following fashion:
     * 1. Using a bounded interval for rho'':[rhoc, rhoc-??], guess a value for rho''
     * 2. For guessed value of rho'' and given value of T, calculate p
     * 3. Using a Brent solver on the other co-existing phase (rho'), calculate the (bounded) value of rho' that yields the same pressure
     * 4. Use another outer Brent solver on rho'' to enforce the same Gibbs function between liquid and vapor
     * 5. Fin.
     * 
     * @param HEOS The Helmholtz EOS backend instance to be used
     * @param ykey The CoolProp::parameters key to be imposed - one of iT or iP
     * @param y The value for the imposed variable
     */
    void saturation_critical(HelmholtzEOSMixtureBackend &HEOS, CoolProp::parameters ykey, CoolPropDbl y);
        
    void successive_substitution(HelmholtzEOSMixtureBackend &HEOS,
                                        const CoolPropDbl beta,
                                        CoolPropDbl T,
                                        CoolPropDbl p,
                                        const std::vector<CoolPropDbl> &z,
                                        std::vector<CoolPropDbl> &K,
                                        mixture_VLE_IO &options);
    /** \brief Extract the mole fractions of liquid (x) and vapor (y) given the bulk composition (z), vapor mole fraction and K-factors
     * @param beta Vapor molar fraction [-]
     * @param K K-factors for the components [-]
     * @param z Bulk molar composition [-]
     * @param x Liquid molar composition [-]
     * @param y Vapor molar composition [-]
     */
    void x_and_y_from_K(CoolPropDbl beta, const std::vector<CoolPropDbl> &K, const std::vector<CoolPropDbl> &z, std::vector<CoolPropDbl> &x, std::vector<CoolPropDbl> &y);

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
        const std::vector<CoolPropDbl> *z;
        std::vector<CoolPropDbl> *K;
        HelmholtzEOSMixtureBackend *HEOS;

        WilsonK_resid(HelmholtzEOSMixtureBackend &HEOS, double beta, double imposed_value, int input_type, const std::vector<CoolPropDbl> &z, std::vector<CoolPropDbl> &K){ 
            this->z = &z; this->K = &K; this->HEOS = &HEOS; this->beta = beta; this->input_type = input_type;
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
                (*K)[i] = exp(Wilson_lnK_factor(*HEOS,T,p,i));
                summer += (*z)[i]*((*K)[i]-1)/(1-beta+beta*(*K)[i]);
            }
            return summer;
        };
    };
    inline double saturation_preconditioner(HelmholtzEOSMixtureBackend &HEOS, double input_value, int input_type, const std::vector<CoolPropDbl> &z)
    {
        double ptriple = 0, pcrit = 0, Ttriple = 0, Tcrit = 0;
        
        for (unsigned int i = 0; i < z.size(); i++)
        {
            EquationOfState *EOS = (HEOS.get_components())[i]->pEOS; 

            ptriple += EOS->sat_min_liquid.p*z[i];
            pcrit += EOS->reduce.p*z[i];
            Ttriple += EOS->sat_min_liquid.T*z[i];
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
    inline double saturation_Wilson(HelmholtzEOSMixtureBackend &HEOS, double beta, double input_value, int input_type, const std::vector<CoolPropDbl> &z, double guess)
    {
        double T;

        std::string errstr;

        // Find first guess for T using Wilson K-factors
        WilsonK_resid Resid(HEOS, beta, input_value, input_type, z, HEOS.get_K());
        T = Secant(Resid, guess, 0.001, 1e-10, 100, errstr);
    
        if (!ValidNumber(T)){throw ValueError("saturation_p_Wilson failed to get good T");}
        return T;
    }
    struct SuccessiveSubstitutionStep
    {
        CoolPropDbl T,p;
    };
    
    struct newton_raphson_twophase_options{
        enum imposed_variable_options {NO_VARIABLE_IMPOSED = 0, P_IMPOSED, T_IMPOSED};
        int Nstep_max;
        std::size_t Nsteps;
        CoolPropDbl beta, omega, rhomolar_liq, rhomolar_vap, pL, pV, p, T, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap;
        imposed_variable_options imposed_variable;
        std::vector<CoolPropDbl> x, y, z;
        newton_raphson_twophase_options(){ Nstep_max = 30; Nsteps = 0; beta = -1; omega =1;} // Defaults
    };

    /** \brief A class to do newton raphson solver for mixture VLE for p,Q or T,Q
     * 
     * A class is used rather than a function so that it is easier to store iteration histories, additional output values, etc.
     * 
     * As in Gernert, FPE, 2014, except that only one of T and P are known
     * 
     * The independent variables are \f$N-1\f$ mole fractions in liquid, \f$N-1\f$ mole fractions in vapor, and the non-specified variable in p or T, for a total of \f$2N-1\f$ independent variables
     * 
     * First N residuals are from
     * 
     * \f$F_k = \ln f_i(T,p,\mathbf{x}) - \ln f_i(T,p,\mathbf{y})\f$ for \f$i = 1, ... N\f$ and \f$k=i\f$
     * 
     * Derivatives are the same as for the saturation solver \ref newton_raphson_saturation
     * 
     * Second N-1 residuals are from
     * 
     * \f$F_k = \dfrac{z_i-x_i}{y_i-x_i} - \beta_{spec}\f$ for \f$ i = 1, ... N-2\f$ and \f$k = i+N\f$
     * 
     * Gernert eq. 35
     * 
     * \f$\dfrac{\partial F_k}{\partial x_i} = \dfrac{z_i-y_i}{(y_i-x_i)^2}\f$
     * 
     * Gernert eq. 36
     * 
     * \f$\dfrac{\partial F_k}{\partial y_i} = -\dfrac{z_i-x_i}{(y_i-x_i)^2}\f$
     * 
     * \f$\dfrac{\partial F_k}{\partial T} = 0\f$ Because x, y and T are independent by definition of the formulation
     * 
     * \f$\dfrac{\partial F_k}{\partial p} = 0\f$ Because x, y and p are independent by definition of the formulation
     */
    class newton_raphson_twophase
    {
        public:
        newton_raphson_twophase_options::imposed_variable_options imposed_variable;
        CoolPropDbl error_rms, rhomolar_liq, rhomolar_vap, T, p, min_rel_change, beta;
        std::size_t N;
        bool logging;
        int Nsteps;
        STLMatrix J;
        HelmholtzEOSMixtureBackend *HEOS;
        std::vector<CoolPropDbl> K, x, y, z, r, negative_r, err_rel;
        std::vector<SuccessiveSubstitutionStep> step_logger;

        newton_raphson_twophase(){};

        void resize(unsigned int N);
    
        // Reset the state of all the internal variables
        void pre_call()
        {
            K.clear(); x.clear(); y.clear();  step_logger.clear(); error_rms = 1e99; Nsteps = 0;
            rhomolar_liq = _HUGE; rhomolar_vap = _HUGE; T = _HUGE; p = _HUGE;
        };

        /** \brief Call the Newton-Raphson VLE Solver
         *
         * This solver must be passed reasonable guess values for the mole fractions, 
         * densities, etc.  You may want to take a few steps of successive substitution
         * before you start with Newton Raphson.
         * 
         * @param HEOS HelmholtzEOSMixtureBackend instance
         * @param IO The input/output data structure
         */
        void call(HelmholtzEOSMixtureBackend &HEOS, newton_raphson_twophase_options &IO);

        /* \brief Build the arrays for the Newton-Raphson solve
         * 
         */
        void build_arrays();
    };
    

    struct newton_raphson_saturation_options{
        enum imposed_variable_options {NO_VARIABLE_IMPOSED = 0, P_IMPOSED, RHOV_IMPOSED, T_IMPOSED};
        int Nstep_max;
        bool bubble_point;
        std::size_t Nsteps;
        CoolPropDbl omega, rhomolar_liq, rhomolar_vap, pL, pV, p, T, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap;
        imposed_variable_options imposed_variable;
        std::vector<CoolPropDbl> x, y;
        newton_raphson_saturation_options(){ Nstep_max = 30;  Nsteps = 0;} // Defaults
    };

    /** \brief A class to do newton raphson solver mixture bubble point and dew point calculations
     * 
     * A class is used rather than a function so that it is easier to store iteration histories, additional output 
     * values, etc.  This class is used in \ref PhaseEnvelopeRoutines for the construction of the phase envelope
     * 
     * This class only handles bubble and dew lines.  The independent variables are the first N-1 mole fractions 
     * in the incipient phase along with one of T, p, or \f$\rho''\f$.
     * 
     * These methods are based on the work of Gernert, FPE, 2014, the thesis of Gernert, as well as much 
     * help from Andreas Jaeger of Uni. Bochum.
     * 
     * There are N residuals from
     * 
     * \f$F_i = \ln f_i(T, p, \mathbf{x}) - \ln f_i(T, p, \mathbf{y})\f$ for \f$i = 1, ... N\f$
     * 
     * if either T or p are imposed. In this case a solver is used to find \f$\rho\f$ given \f$T\f$ and \f$p\f$.  
     * 
     * In the case that \f$\rho''\f$ is imposed, the case is much nicer and the first \f$N\f$ residuals are
     * 
     * \f$F_i = \ln f_i(T, \rho', \mathbf{x}) - \ln f_i(T, \rho'', \mathbf{y})\f$ for \f$i = 1, ... N\f$ 
     * 
     * which requires no iteration. A final residual is set up to ensure the same pressures in the phases:
     * 
     * \f$p(T,\rho',\mathbf{x})-p(T,\rho'',\mathbf{y}) = 0\f$
     * 
     * Documentation of the derivatives needed can be found in the work of Gernert, FPE, 2014
     */
    class newton_raphson_saturation
    {
        public:
        newton_raphson_saturation_options::imposed_variable_options imposed_variable;
        CoolPropDbl error_rms, rhomolar_liq, rhomolar_vap, T, p, min_rel_change;
        std::size_t N;
        bool logging;
        bool bubble_point;
        int Nsteps;
        STLMatrix J;
        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl dTsat_dPsat, dPsat_dTsat;
        std::vector<CoolPropDbl> K, x, y, r, negative_r, err_rel;
        std::vector<SuccessiveSubstitutionStep> step_logger;

        newton_raphson_saturation(){};

        void resize(std::size_t N);
    
        // Reset the state of all the internal variables
        void pre_call()
        {
            K.clear(); x.clear(); y.clear();  
            step_logger.clear(); error_rms = 1e99; Nsteps = 0;
            rhomolar_liq = _HUGE; rhomolar_vap = _HUGE; T = _HUGE; p = _HUGE;
        };

        /** \brief Call the Newton-Raphson VLE Solver
         *
         * This solver must be passed reasonable guess values for the mole fractions, 
         * densities, etc.  You may want to take a few steps of successive substitution
         * before you start with Newton Raphson.
         * 
         * @param HEOS HelmholtzEOSMixtureBackend instance
         * @param z Bulk mole fractions [-]
         * @param z_incipient Initial guesses for the mole fractions of the incipient phase [-]
         * @param IO The input/output data structure
         */
        void call(HelmholtzEOSMixtureBackend &HEOS, const std::vector<CoolPropDbl> &z, std::vector<CoolPropDbl> &z_incipient, newton_raphson_saturation_options &IO);

        /** \brief Build the arrays for the Newton-Raphson solve
         * 
         * This method builds the Jacobian matrix, the sensitivity matrix, etc.
         * 
         */
        void build_arrays();

        /** \brief Check the derivatives in the Jacobian using numerical derivatives.
         */
        void check_Jacobian();
    };
};
    
} /* namespace CoolProp*/

#endif
