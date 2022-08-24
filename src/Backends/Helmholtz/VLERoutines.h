
#ifndef VLEROUTINES_H
#define VLEROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

namespace SaturationSolvers {
struct saturation_T_pure_Akasaka_options
{
    bool use_guesses;  ///< true to start off at the values specified by rhoL, rhoV
    CoolPropDbl omega, rhoL, rhoV, pL, pV;
    saturation_T_pure_Akasaka_options(bool use_guesses = false)
      : use_guesses(use_guesses), omega(_HUGE), rhoL(_HUGE), rhoV(_HUGE), pL(_HUGE), pV(_HUGE) {}
};
struct saturation_T_pure_options
{
    CoolPropDbl omega, rhoL, rhoV, pL, pV, p, T;
    saturation_T_pure_options() : omega(_HUGE), rhoL(_HUGE), rhoV(_HUGE), pL(_HUGE), pV(_HUGE), p(_HUGE), T(_HUGE) {}
};

struct saturation_D_pure_options
{
    enum imposed_rho_options
    {
        IMPOSED_RHOL,
        IMPOSED_RHOV
    };
    bool use_guesses,  ///< True to start off at the values specified by rhoL, rhoV, T
      use_logdelta;    ///< True to use partials with respect to log(delta) rather than delta
    CoolPropDbl omega, rhoL, rhoV, pL, pV;
    int imposed_rho;
    saturation_D_pure_options() : use_guesses(false), rhoL(_HUGE), rhoV(_HUGE), pL(_HUGE), pV(_HUGE), imposed_rho(0) {
        use_logdelta = true;
        omega = 1.0;
    }  // Defaults
};

enum sstype_enum
{
    imposed_T,
    imposed_p
};
struct mixture_VLE_IO
{
    sstype_enum sstype;
    int Nstep_max;
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
static CoolPropDbl Wilson_lnK_factor(const HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, CoolPropDbl p, std::size_t i) {
    const double pci = HEOS.get_fluid_constant(i, iP_critical);
    const double Tci = HEOS.get_fluid_constant(i, iT_critical);
    const double omegai = HEOS.get_fluid_constant(i, iacentric_factor);
    return log(pci / p) + 5.373 * (1 + omegai) * (1 - Tci / T);
};

void saturation_D_pure(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl rhomolar, saturation_D_pure_options& options);
void saturation_T_pure(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, saturation_T_pure_options& options);
void saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options& options);
void saturation_T_pure_Maxwell(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options& options);

/**
    */
struct saturation_PHSU_pure_options
{
    enum specified_variable_options
    {
        IMPOSED_HL,
        IMPOSED_HV,
        IMPOSED_PL,
        IMPOSED_PV,
        IMPOSED_SL,
        IMPOSED_SV,
        IMPOSED_UL,
        IMPOSED_UV,
        IMPOSED_INVALID_INPUT
    };
    bool use_guesses,  ///< True to start off at the values specified by rhoL, rhoV, T
      use_logdelta;    ///< True to use partials with respect to log(delta) rather than delta
    specified_variable_options specified_variable;
    CoolPropDbl omega, rhoL, rhoV, pL, pV, T, p;
    saturation_PHSU_pure_options() : use_logdelta(true), rhoL(_HUGE), rhoV(_HUGE), pL(_HUGE), pV(_HUGE), T(_HUGE), p(_HUGE) {
        specified_variable = IMPOSED_INVALID_INPUT;
        use_guesses = true;
        omega = 1.0;
    }
};
/**

    */
void saturation_PHSU_pure(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl specified_value, saturation_PHSU_pure_options& options);

/* \brief This is a backup saturation_p solver for the case where the Newton solver cannot approach closely enough the solution
     *
     * This is especially a problem at low pressures where catastrophic truncation error occurs, especially in the saturated vapor side
     *
     * @param HEOS The Helmholtz EOS backend instance to be used
     * @param p Imposed pressure in kPa
     * @param options Options to be passed to the function (at least T, rhoL and rhoV must be provided)
     */
void saturation_P_pure_1D_T(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl p, saturation_PHSU_pure_options& options);

/* \brief This is a backup saturation_T solver for the case where the Newton solver cannot approach closely enough the solution
     *
     * This is especially a problem at low pressures where catastrophic truncation error occurs, especially in the saturated vapor side
     *
     * @param HEOS The Helmholtz EOS backend instance to be used
     * @param T Imposed temperature in K
     * @param options Options to be passed to the function (at least p, rhoL and rhoV must be provided)
     */
void saturation_T_pure_1D_P(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, saturation_T_pure_options& options);

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
void saturation_critical(HelmholtzEOSMixtureBackend& HEOS, CoolProp::parameters ykey, CoolPropDbl y);

void successive_substitution(HelmholtzEOSMixtureBackend& HEOS, const CoolPropDbl beta, CoolPropDbl T, CoolPropDbl p,
                             const std::vector<CoolPropDbl>& z, std::vector<CoolPropDbl>& K, mixture_VLE_IO& options);
/** \brief Extract the mole fractions of liquid (x) and vapor (y) given the bulk composition (z), vapor mole fraction and K-factors
     * @param beta Vapor molar fraction [-]
     * @param K K-factors for the components [-]
     * @param z Bulk molar composition [-]
     * @param x Liquid molar composition [-]
     * @param y Vapor molar composition [-]
     */
void x_and_y_from_K(CoolPropDbl beta, const std::vector<CoolPropDbl>& K, const std::vector<CoolPropDbl>& z, std::vector<CoolPropDbl>& x,
                    std::vector<CoolPropDbl>& y);

/*! A wrapper function around the residual to find the initial guess for the bubble point temperature
    \f[
    r = \sum_i \frac{z_i(K_i-1)}{1-beta+beta*K_i}
    \f]
    */
class WilsonK_resid : public FuncWrapper1D
{
   public:
    sstype_enum input_type;
    double T, p, beta;
    const std::vector<CoolPropDbl>& z;
    std::vector<CoolPropDbl>& K;
    const HelmholtzEOSMixtureBackend& HEOS;

    WilsonK_resid(const HelmholtzEOSMixtureBackend& HEOS, double beta, double imposed_value, sstype_enum input_type,
                  const std::vector<CoolPropDbl>& z, std::vector<CoolPropDbl>& K)
      : input_type(input_type),
        T(imposed_value),
        p(imposed_value),
        beta(beta),
        z(z),
        K(K),
        HEOS(HEOS) {}  // if input_type == imposed_T -> use T, else use p; init both
    double call(double input_value) {
        double summer = 0;
        if (input_type == imposed_T) {
            p = input_value;  // Iterate on pressure
        } else {
            T = input_value;  // Iterate on temperature, pressure imposed
        }
        for (unsigned int i = 0; i < z.size(); i++) {
            K[i] = exp(Wilson_lnK_factor(HEOS, T, p, i));
            summer += z[i] * (K[i] - 1) / (1 - beta + beta * K[i]);
        }
        return summer;
    }
};
inline double saturation_preconditioner(HelmholtzEOSMixtureBackend& HEOS, double input_value, sstype_enum input_type,
                                        const std::vector<CoolPropDbl>& z) {
    double ptriple = 0, pcrit = 0, Ttriple = 0, Tcrit = 0;

    if (HEOS.get_components().empty()) {
        return -1;
    }

    for (unsigned int i = 0; i < z.size(); i++) {
        Tcrit += HEOS.get_fluid_constant(i, iT_critical) * z[i];
        pcrit += HEOS.get_fluid_constant(i, iP_critical) * z[i];
        Ttriple += HEOS.get_fluid_constant(i, iT_triple) * z[i];
        ptriple += HEOS.get_fluid_constant(i, iP_triple) * z[i];
    }
    // Return an invalid number if either triple point temperature or pressure are not available
    if (!ValidNumber(Ttriple) || !ValidNumber(ptriple)) {
        return _HUGE;
    }

    if (input_type == imposed_T) {
        return exp(log(pcrit / ptriple) / (Tcrit - Ttriple) * (input_value - Ttriple) + log(ptriple));
    } else if (input_type == imposed_p) {
        return 1 / (1 / Tcrit - (1 / Ttriple - 1 / Tcrit) / log(pcrit / ptriple) * log(input_value / pcrit));
    } else {
        throw ValueError();
    }
}
/**
     * Wilson gives the K-factor as
     * \f[
     * \ln K_i = \ln\left(\frac{p_{c,i}}{p}\right)+5.373(1+\omega_i)\left(1-\frac{T_{c,i}}{T}\right)
     * \f]
     *
     * From Rachford-Rice:
     * \f[
     * \sum_i \frac{x_i(K_i-1)}{1 - \beta + \beta K_i} = 0
     * \f]
     * When \f$T\f$ is known for \f$\beta=0$, \f$p\f$can be obtained from
     * \f[
     * -1+\sum_i K_ix_i=0,
     * \f]
     * or
     * \f[
     * p = \sum_i x_ip_{c,i}\exp(5.373(1+\omega_i)(1-T_{c,i}/T).
     * \f]
     * Or when \f$T\f$ is known for \f$\beta=1\f$, \f$p\f$can be obtained from
     * \f[
     * -1+\sum_ix_i=0,
     * \f]
     * or
     * \f[
     * p = \left[ \sum_i \frac{y_i}{p_{c,i}\exp(5.373(1+\omega_i)(1-T_{c,i}/T)} \right]^{-1}
     * \f]
     */
inline double saturation_Wilson(HelmholtzEOSMixtureBackend& HEOS, double beta, double input_value, sstype_enum input_type,
                                const std::vector<CoolPropDbl>& z, double guess) {
    double out = 0;
    std::string errstr;

    // If T is input and beta = 0 or beta = 1, explicit solution for p is possible
    if (input_type == imposed_T && (std::abs(beta) < 1e-12 || std::abs(beta - 1) < 1e-12)) {
        const std::vector<double> z = HEOS.get_mole_fractions_ref();
        bool beta0 = std::abs(beta) < 1e-12;  // True is beta is approx. zero
        for (int i = 0; i < static_cast<int>(z.size()); ++i) {
            double pci = HEOS.get_fluid_constant(i, iP_critical);
            double Tci = HEOS.get_fluid_constant(i, iT_critical);
            double omegai = HEOS.get_fluid_constant(i, iacentric_factor);
            if (beta0) {
                out += z[i] * pci * exp(5.373 * (1 + omegai) * (1 - Tci / input_value));
            } else {
                out += z[i] / (pci * exp(5.373 * (1 + omegai) * (1 - Tci / input_value)));
            }
        }
        if (!beta0) {       // beta = 1
            out = 1 / out;  // summation is for 1/p, take reciprocal to get p
        }
        std::vector<CoolPropDbl>& K = HEOS.get_K();
        for (int i = 0; i < static_cast<int>(z.size()); ++i) {
            double pci = HEOS.get_fluid_constant(i, iP_critical);
            double Tci = HEOS.get_fluid_constant(i, iT_critical);
            double omegai = HEOS.get_fluid_constant(i, iacentric_factor);
            K[i] = pci / out * exp(5.373 * (1 + omegai) * (1 - Tci / input_value));
        }
    } else {
        // Find first guess for output variable using Wilson K-factors
        WilsonK_resid Resid(HEOS, beta, input_value, input_type, z, HEOS.get_K());
        if (guess < 0 || !ValidNumber(guess))
            out = Brent(Resid, 50, 10000, 1e-10, 1e-10, 100);
        else
            out = Secant(Resid, guess, 0.001, 1e-10, 100);
        if (!ValidNumber(out)) {
            throw ValueError("saturation_p_Wilson failed to get good output value");
        }
    }
    return out;
}
struct SuccessiveSubstitutionStep
{
    CoolPropDbl T, p;
};

struct newton_raphson_twophase_options
{
    enum imposed_variable_options
    {
        NO_VARIABLE_IMPOSED = 0,
        P_IMPOSED,
        T_IMPOSED
    };
    int Nstep_max;
    std::size_t Nsteps;
    CoolPropDbl beta, omega, rhomolar_liq, rhomolar_vap, pL, pV, p, T, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap;
    imposed_variable_options imposed_variable;
    std::vector<CoolPropDbl> x, y, z;
    newton_raphson_twophase_options()
      : Nstep_max(30),
        Nsteps(0),
        beta(-1),
        omega(1),
        rhomolar_liq(_HUGE),
        rhomolar_vap(_HUGE),
        pL(_HUGE),
        pV(_HUGE),
        p(_HUGE),
        T(_HUGE),
        hmolar_liq(_HUGE),
        hmolar_vap(_HUGE),
        smolar_liq(_HUGE),
        smolar_vap(_HUGE),
        imposed_variable(NO_VARIABLE_IMPOSED) {}  // Defaults
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
    HelmholtzEOSMixtureBackend* HEOS;
    newton_raphson_twophase_options::imposed_variable_options imposed_variable;
    CoolPropDbl error_rms, rhomolar_liq, rhomolar_vap, T, p, min_rel_change, beta;
    std::size_t N;
    bool logging;
    int Nsteps;
    Eigen::MatrixXd J;
    Eigen::Vector2d r, err_rel;
    std::vector<CoolPropDbl> K, x, y, z;
    std::vector<SuccessiveSubstitutionStep> step_logger;

    newton_raphson_twophase()
      : HEOS(NULL),
        imposed_variable(newton_raphson_twophase_options::NO_VARIABLE_IMPOSED),
        error_rms(_HUGE),
        rhomolar_liq(_HUGE),
        rhomolar_vap(_HUGE),
        T(_HUGE),
        p(_HUGE),
        min_rel_change(_HUGE),
        beta(_HUGE),
        N(0),
        logging(false),
        Nsteps(0){};

    void resize(unsigned int N);

    // Reset the state of all the internal variables
    void pre_call() {
        K.clear();
        x.clear();
        y.clear();
        step_logger.clear();
        error_rms = 1e99;
        Nsteps = 0;
        rhomolar_liq = _HUGE;
        rhomolar_vap = _HUGE;
        T = _HUGE;
        p = _HUGE;
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
    void call(HelmholtzEOSMixtureBackend& HEOS, newton_raphson_twophase_options& IO);

    /* \brief Build the arrays for the Newton-Raphson solve
         *
         */
    void build_arrays();
};

struct newton_raphson_saturation_options
{
    enum imposed_variable_options
    {
        NO_VARIABLE_IMPOSED = 0,
        P_IMPOSED,
        RHOV_IMPOSED,
        T_IMPOSED
    };
    int Nstep_max;
    bool bubble_point;
    std::size_t Nsteps;
    CoolPropDbl omega, rhomolar_liq, rhomolar_vap, pL, pV, p, T, hmolar_liq, hmolar_vap, smolar_liq, smolar_vap;
    imposed_variable_options imposed_variable;
    std::vector<CoolPropDbl> x, y;
    newton_raphson_saturation_options()
      : bubble_point(false),
        omega(_HUGE),
        rhomolar_liq(_HUGE),
        rhomolar_vap(_HUGE),
        pL(_HUGE),
        pV(_HUGE),
        p(_HUGE),
        T(_HUGE),
        hmolar_liq(_HUGE),
        hmolar_vap(_HUGE),
        smolar_liq(_HUGE),
        smolar_vap(_HUGE),
        imposed_variable(NO_VARIABLE_IMPOSED) {
        Nstep_max = 30;
        Nsteps = 0;
    }  // Defaults
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
    Eigen::MatrixXd J;
    HelmholtzEOSMixtureBackend* HEOS;
    CoolPropDbl dTsat_dPsat, dPsat_dTsat;
    std::vector<CoolPropDbl> K, x, y;
    Eigen::VectorXd r, err_rel;
    std::vector<SuccessiveSubstitutionStep> step_logger;

    newton_raphson_saturation(){};

    void resize(std::size_t N);

    // Reset the state of all the internal variables
    void pre_call() {
        step_logger.clear();
        error_rms = 1e99;
        Nsteps = 0;
        rhomolar_liq = _HUGE;
        rhomolar_vap = _HUGE;
        T = _HUGE;
        p = _HUGE;
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
    void call(HelmholtzEOSMixtureBackend& HEOS, const std::vector<CoolPropDbl>& z, std::vector<CoolPropDbl>& z_incipient,
              newton_raphson_saturation_options& IO);

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

struct PTflash_twophase_options
{
    int Nstep_max;
    std::size_t Nsteps;
    CoolPropDbl omega, rhomolar_liq, rhomolar_vap, pL, pV, p, T;
    std::vector<CoolPropDbl> x,  ///< Liquid mole fractions
      y,                         ///< Vapor mole fractions
      z;                         ///< Bulk mole fractions
    PTflash_twophase_options() : omega(_HUGE), rhomolar_liq(_HUGE), rhomolar_vap(_HUGE), pL(_HUGE), pV(_HUGE), p(_HUGE), T(_HUGE) {
        Nstep_max = 30;
        Nsteps = 0;  // Defaults
    }
};

class PTflash_twophase
{
   public:
    double error_rms;
    bool logging;
    int Nsteps;
    Eigen::MatrixXd J;
    Eigen::VectorXd r, err_rel;
    HelmholtzEOSMixtureBackend& HEOS;
    PTflash_twophase_options& IO;
    std::vector<SuccessiveSubstitutionStep> step_logger;

    PTflash_twophase(HelmholtzEOSMixtureBackend& HEOS, PTflash_twophase_options& IO) : HEOS(HEOS), IO(IO){};

    /** \brief Call the Newton-Raphson Solver to solve the equilibrium conditions
         *
         * This solver must be passed reasonable guess values for the mole fractions,
         * densities, etc.  You may want to take a few steps of successive substitution
         * before you start with Newton Raphson.
         *
         */
    void solve();

    /** \brief Build the arrays for the Newton-Raphson solve
         *
         * This method builds the Jacobian matrix, the sensitivity matrix, etc.
         *
         */
    void build_arrays();
};
};  // namespace SaturationSolvers

namespace StabilityRoutines {

/** \brief Evaluate phase stability
     * Based on the work of Gernert et al., J. Chem. Thermodyn., 2014 http://dx.doi.org/10.1016/j.fluid.2014.05.012
     */
class StabilityEvaluationClass
{
   protected:
    HelmholtzEOSMixtureBackend& HEOS;
    std::vector<double> lnK, K, K0, x, y, xL, xH;
    const std::vector<double>& z;
    double rhomolar_liq, rhomolar_vap, beta, tpd_liq, tpd_vap, DELTAG_nRT;
    double m_T,  ///< The temperature to be used (if specified, otherwise that from HEOS)
      m_p;       ///< The pressure to be used (if specified, otherwise that from HEOS)
   private:
    bool _stable;
    bool debug;

   public:
    StabilityEvaluationClass(HelmholtzEOSMixtureBackend& HEOS)
      : HEOS(HEOS),
        z(HEOS.get_mole_fractions_doubleref()),
        rhomolar_liq(-1),
        rhomolar_vap(-1),
        beta(-1),
        tpd_liq(10000),
        tpd_vap(100000),
        DELTAG_nRT(10000),
        m_T(-1),
        m_p(-1),
        _stable(false),
        debug(false){};
    /** \brief Specify T&P, otherwise they are loaded the HEOS instance
         */
    void set_TP(double T, double p) {
        m_T = T;
        m_p = p;
    };
    /** \brief Calculate the liquid and vapor phase densities based on the guess values
         */
    void rho_TP_w_guesses();
    /** \brief Calculate the liquid and vapor phase densities using the global analysis
         */
    void rho_TP_global();
    /** \brief Calculate the liquid and vapor phase densities based on SRK, with Peneloux volume translation afterwards
         */
    void rho_TP_SRK_translated();

    /** \brief Calculate trial compositions
         */
    void trial_compositions();
    /** \brief Successive substitution
         */
    void successive_substitution(int num_steps);
    /** \brief Check stability
         * 1. Check stability by looking at tpd', tpd'' and \f$ \Delta G/(nRT)\f$
         * 2. Do a full TPD analysis
         */
    void check_stability();
    /** \brief Return best estimate for the stability of the point
         */
    bool is_stable() {
        trial_compositions();
        successive_substitution(3);
        check_stability();
        return _stable;
    }
    /// Accessor for liquid-phase composition and density
    void get_liq(std::vector<double>& x, double& rhomolar) {
        x = this->x;
        rhomolar = rhomolar_liq;
    }
    /// Accessor for vapor-phase composition and density
    void get_vap(std::vector<double>& y, double& rhomolar) {
        y = this->y;
        rhomolar = rhomolar_vap;
    }
};

}; /* namespace StabilityRoutines*/

}; /* namespace CoolProp*/

#endif
