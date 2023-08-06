/**
This file contains flash routines in which the state is unknown,
and a solver of some kind must be used to obtain temperature and
density, the two state variables upon which the equation of
state is based.
*/

// ***************************************************************
// *******************  FLASH ROUTINES  **************************
// ***************************************************************

#ifndef FLASHROUTINES_H
#define FLASHROUTINES_H

#include "HelmholtzEOSMixtureBackend.h"
#include "Solvers.h"

namespace CoolProp {

/**
This class is a friend class of HelmholtzEOSMixtureBackend, therefore the
static methods contained in it have access to the private and
protected variables in the HelmholtzEOSMixtureBackend instance.

In this way the Flash routines can be kept in their own separate file
and not pollute the HelmholtzEOSMixtureBackend namespace
*/
class FlashRoutines
{
   public:
    template <class T>
    T static g_RachfordRice(const std::vector<T>& z, const std::vector<T>& lnK, T beta) {
        // g function from Rachford-Rice
        T summer = 0;
        for (std::size_t i = 0; i < z.size(); i++) {
            T Ki = exp(lnK[i]);
            summer += z[i] * (Ki - 1) / (1 - beta + beta * Ki);
        }
        return summer;
    }
    template <class T>
    T static dgdbeta_RachfordRice(const std::vector<T>& z, const std::vector<T>& lnK, T beta) {
        // derivative of g function from Rachford-Rice with respect to beta
        T summer = 0;
        for (std::size_t i = 0; i < z.size(); i++) {
            T Ki = exp(lnK[i]);
            summer += -z[i] * pow((Ki - 1) / (1 - beta + beta * Ki), 2);
        }
        return summer;
    }

    /// Flash for given pressure and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void PQ_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Flash for given pressure and (molar) quality with guess values provided
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param guess The GuessesStructure to be used
    static void PQ_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess);

    /// Flash for given temperature and (molar) quality with guess values provided
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param guess The GuessesStructure to be used
    static void QT_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess);

    /// Flash for given pressure and temperature with guess values provided for molar density
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param guess The GuessesStructure to be used
    static void PT_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess);

    /// Flash for given temperature and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void QT_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Flash for given molar entropy and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void QS_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Flash for given molar density and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void DQ_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Flash for given molar enthalpy and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param Tguess (optional) The guess temperature in K to start from, ignored if < 0
    static void HQ_flash(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl Tguess = -1);

    /// Flash for mixture given temperature or pressure and (molar) quality
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The parameter that is imposed, either iT or iP
    /// @param value The value for the imposed parameter
    static void PT_Q_flash_mixtures(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl value);

    /// Flash for given pressure and temperature
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void PT_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Flash for given pressure and temperature for mixtures
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void PT_flash_mixtures(HelmholtzEOSMixtureBackend& HEOS);

    /// Use Peng-Robinson to get guess for temperature for given density and pressure
    static double T_DP_PengRobinson(HelmholtzEOSMixtureBackend& HEOS, double rhomolar, double p);

    /// Flash for given density and pressure
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void DP_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// The flash routine for T given and one of H,S,U
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param T The temperature
    /// @param value The value for the other variable
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iDmolar, iHmolar, iSmolar, iUmolar
    static void solver_for_rho_given_T_oneof_HSU(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, CoolPropDbl value, parameters other);

    /// A generic flash routine for the pairs (T,D), (T,H), (T,S), and (T,U).  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iDmolar, iHmolar, iSmolar, iUmolar
    static void DHSU_T_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other);

    /// A generic flash routine for the pairs (P,H), (P,S), and (P,U).  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iHmolar, iSmolar, iUmolar
    static void HSU_P_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other);

    /// The single-phase flash routine for the pairs (P,H), (P,S), and (P,U).  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iHmolar, iSmolar, iUmolar
    /// @param T0 The initial guess value for the temperature [K]
    /// @param rhomolar0 The initial guess value for the density [mol/m^3]
    static void HSU_P_flash_singlephase_Newton(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl T0, CoolPropDbl rhomolar0);

    /// The single-phase flash routine for the pairs (P,H), (P,S), and (P,U).  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iHmolar, iSmolar, iUmolar
    /// @param value The value of the other input
    /// @param Tmin The lower temperature limit [K]
    /// @param Tmax The higher temperature limit [K]
    /// @param phase The phase of the fluid that we should get back
    static void HSU_P_flash_singlephase_Brent(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl value, CoolPropDbl Tmin,
                                              CoolPropDbl Tmax, phases phase);

    /// A generic flash routine for the pairs (D,H), (D,S), and (D,U) for twophase state.  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iP, iHmolar, iSmolar, iUmolar
    static void HSU_D_flash_twophase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value);

    /// A generic flash routine for the pairs (D,P), (D,H), (D,S), and (D,U).  Similar analysis is needed
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param other The index for the other input from CoolProp::parameters; allowed values are iP, iHmolar, iSmolar, iUmolar
    static void HSU_D_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other);

    /// A flash routine for (H,S)
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    static void HS_flash(HelmholtzEOSMixtureBackend& HEOS);

    /// Randomly generate a single phase set of inputs for T and p - searches entire single-phase region
    /// @param HEOS The HelmholtzEOSMixtureBackend to be used
    /// @param T The temperature in K
    /// @param p The pressure in Pa
    static void HS_flash_generate_TP_singlephase_guess(HelmholtzEOSMixtureBackend& HEOS, double& T, double& p);

    struct HS_flash_singlephaseOptions
    {
        double omega;
        HS_flash_singlephaseOptions() {
            omega = 1.0;
        }
    };
    static void HS_flash_singlephase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec,
                                     HS_flash_singlephaseOptions& options);

    struct HS_flash_twophaseOptions
    {
        double omega;
        HS_flash_twophaseOptions() {
            omega = 1.0;
        }
    };
    static void HS_flash_twophase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec,
                                  HS_flash_twophaseOptions& options);
};

/** A residual function for the rho(T,P) solver
 */
class solver_TP_resid : public FuncWrapper1DWithDeriv
{
   public:
    HelmholtzEOSMixtureBackend* HEOS;
    CoolPropDbl T, p, rhor, tau, R_u, delta;

    solver_TP_resid(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, CoolPropDbl p)
      : HEOS(&HEOS),
        T(T),
        p(p),
        rhor(HEOS.get_reducing_state().rhomolar),
        tau(HEOS.get_reducing_state().T / T),
        R_u(HEOS.gas_constant()),
        delta(-_HUGE) {}
    double call(double rhomolar) {
        delta = rhomolar / rhor;  // needed for derivative
        HEOS->update_DmolarT_direct(rhomolar, T);
        CoolPropDbl peos = HEOS->p();
        return (peos - p) / p;
    };
    double deriv(double rhomolar) {
        // dp/drho|T / pspecified
        return R_u * T * (1 + 2 * delta * HEOS->dalphar_dDelta() + pow(delta, 2) * HEOS->d2alphar_dDelta2()) / p;
    };
};

/** A residual function for the f(P, Y) solver
 */
class PY_singlephase_flash_resid : public FuncWrapper1D
{
   public:
    HelmholtzEOSMixtureBackend* HEOS;
    CoolPropDbl p;
    parameters other;
    CoolPropDbl value;
    PY_singlephase_flash_resid(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl p, parameters other, CoolPropDbl value)
      : HEOS(&HEOS), p(p), other(other), value(value) {
        // Specify the state to avoid saturation calls, but only if phase is subcritical
        if (HEOS.phase() == iphase_liquid || HEOS.phase() == iphase_gas) {
            HEOS.specify_phase(HEOS.phase());
        }
    };
    double call(double T) {

        // Run the solver with T,P as inputs;
        HEOS->update(PT_INPUTS, p, T);

        CoolPropDbl rhomolar = HEOS->rhomolar();
        HEOS->update(DmolarT_INPUTS, rhomolar, T);
        // Get the value of the desired variable
        CoolPropDbl eos = HEOS->keyed_output(other);

        // Difference between the two is to be driven to zero
        return eos - value;
    };
};

} /* namespace CoolProp */
#endif /* FLASHROUTINES_H */
