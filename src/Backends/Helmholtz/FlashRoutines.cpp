#include "VLERoutines.h"
#include "FlashRoutines.h"

#include <cmath>
#include "HelmholtzEOSMixtureBackend.h"
#include "HelmholtzEOSBackend.h"
#include "PhaseEnvelopeRoutines.h"
#include "Configuration.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>
#    include "Backends/Cubics/CubicBackend.h"
#endif

#include "boost/math/tools/toms748_solve.hpp"

namespace CoolProp {

void FlashRoutines::PT_flash_mixtures(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.PhaseEnvelope.built) {
        // Use the phase envelope if already constructed to determine phase boundary
        // Determine whether you are inside (two-phase) or outside (single-phase)
        SimpleState closest_state;
        std::size_t i = 0;
        bool twophase = PhaseEnvelopeRoutines::is_inside(HEOS.PhaseEnvelope, iP, HEOS._p, iT, HEOS._T, i, closest_state);
        if (!twophase && HEOS._T > closest_state.T) {
            // Gas solution - bounded between phase envelope temperature and very high temperature
            //
            // Start with a guess value from SRK
            CoolPropDbl rhomolar_guess = HEOS.solver_rho_Tp_SRK(HEOS._T, HEOS._p, iphase_gas);

            solver_TP_resid resid(HEOS, HEOS._T, HEOS._p);
            std::string errstr;
            HEOS.specify_phase(iphase_gas);
            try {
                // Try using Newton's method
                CoolPropDbl rhomolar = Newton(resid, rhomolar_guess, 1e-10, 100);
                // Make sure the solution is within the bounds
                if (!is_in_closed_range(static_cast<CoolPropDbl>(closest_state.rhomolar), static_cast<CoolPropDbl>(0.0), rhomolar)) {
                    throw ValueError("out of range");
                }
                HEOS.update_DmolarT_direct(rhomolar, HEOS._T);
            } catch (...) {
                // If that fails, try a bounded solver
                CoolPropDbl rhomolar = Brent(resid, closest_state.rhomolar, 1e-10, DBL_EPSILON, 1e-10, 100);
                // Make sure the solution is within the bounds
                if (!is_in_closed_range(static_cast<CoolPropDbl>(closest_state.rhomolar), static_cast<CoolPropDbl>(0.0), rhomolar)) {
                    throw ValueError("out of range");
                }
            }
            HEOS.unspecify_phase();
            HEOS._Q = -1;
        } else {
            // Liquid solution
            throw ValueError();
        }
    } else {
        if (HEOS.imposed_phase_index == iphase_not_imposed) {
            // Blind flash call
            // Following the strategy of Gernert, 2014
            StabilityRoutines::StabilityEvaluationClass stability_tester(HEOS);
            if (!stability_tester.is_stable()) {
                // There is a phase split and liquid and vapor phases are formed
                CoolProp::SaturationSolvers::PTflash_twophase_options o;
                stability_tester.get_liq(o.x, o.rhomolar_liq);
                stability_tester.get_vap(o.y, o.rhomolar_vap);
                o.z = HEOS.get_mole_fractions();
                o.T = HEOS.T();
                o.p = HEOS.p();
                o.omega = 1.0;
                CoolProp::SaturationSolvers::PTflash_twophase solver(HEOS, o);
                solver.solve();
                HEOS._phase = iphase_twophase;
                HEOS._Q = (o.z[0] - o.x[0]) / (o.y[0] - o.x[0]);  // All vapor qualities are the same (these are the residuals in the solver)
                HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
            } else {
                // It's single-phase
                double rho = HEOS.solver_rho_Tp_global(HEOS.T(), HEOS.p(), 20000);
                HEOS.update_DmolarT_direct(rho, HEOS.T());
                HEOS._Q = -1;
                HEOS._phase = iphase_liquid;
            }
        } else {
            // It's single-phase, and phase is imposed
            double rho = HEOS.solver_rho_Tp(HEOS.T(), HEOS.p());
            HEOS.update_DmolarT_direct(rho, HEOS.T());
            HEOS._Q = -1;
            HEOS._phase = HEOS.imposed_phase_index;
        }
    }
}
void FlashRoutines::PT_flash(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {
        // At the critical point dP/drho -> 0, so solver_rho_Tp is ill-conditioned: a tight
        // pressure residual does not imply a tight density. Short-circuit when (T, p)
        // matches the critical point within 1e-10 relative on both axes. See issue #2738.
        CoolPropDbl Tc = HEOS.T_critical();
        CoolPropDbl pc = HEOS.p_critical();
        if (is_in_closed_range(Tc * (1 - 1e-10), Tc * (1 + 1e-10), static_cast<CoolPropDbl>(HEOS._T))
            && is_in_closed_range(pc * (1 - 1e-10), pc * (1 + 1e-10), static_cast<CoolPropDbl>(HEOS._p))) {
            HEOS._phase = iphase_critical_point;
            HEOS._T = Tc;
            HEOS._p = pc;
            HEOS._rhomolar = HEOS.rhomolar_critical();
            HEOS._Q = -1;
            return;
        }
        if (HEOS.imposed_phase_index == iphase_not_imposed)  // If no phase index is imposed (see set_components function)
        {
            // At very low temperature (near the triple point temp), the isotherms are VERY steep
            // Thus it can be very difficult to determine state based on ps = f(T)
            // So in this case, we do a phase determination based on p, generally it will be useful enough
            if (HEOS._T < 0.9 * HEOS.Ttriple() + 0.1 * HEOS.calc_Tmax_sat()) {
                // Find the phase, while updating all internal variables possible using the pressure
                bool saturation_called = false;
                HEOS.p_phase_determination_pure_or_pseudopure(iT, HEOS._T, saturation_called);
            } else {
                // Find the phase, while updating all internal variables possible using the temperature
                HEOS.T_phase_determination_pure_or_pseudopure(iP, HEOS._p);
            }
            // Check if twophase solution
            if (!HEOS.isHomogeneousPhase()) {
                throw ValueError("twophase not implemented yet");
            }
        } else {
            // Phase is imposed.  Update _phase in case it was reset elsewhere by another call
            HEOS._phase = HEOS.imposed_phase_index;
        }
        // Find density
        HEOS._rhomolar = HEOS.solver_rho_Tp(HEOS._T, HEOS._p);
        HEOS._Q = -1;
    } else {
        PT_flash_mixtures(HEOS);
    }
}

// Define the residual to be driven to zero
class solver_DP_resid : public FuncWrapper1DWithTwoDerivs
{
   public:
    HelmholtzEOSMixtureBackend* HEOS;
    CoolPropDbl rhomolar, p;
    solver_DP_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl rhomolar, CoolPropDbl p) : HEOS(HEOS), rhomolar(rhomolar), p(p) {}
    double call(double T) {
        HEOS->update_DmolarT_direct(rhomolar, T);
        CoolPropDbl peos = HEOS->p();
        CoolPropDbl r = (peos - p) / p;
        return r;
    };
    double deriv(double T) {
        // dp/dT|rho / pspecified
        return HEOS->first_partial_deriv(iP, iT, iDmolar) / p;
    };
    double second_deriv(double T) {
        // d2p/dT2|rho / pspecified
        return HEOS->second_partial_deriv(iP, iT, iDmolar, iT, iDmolar) / p;
    };
};

/***
\f[
\begin{array}{l}
p = \frac{{RT}}{{v - b}} - \frac{{a\alpha }}{{v\left( {v + b} \right)}}\\
\alpha  = \left( {1 + \kappa \left( {1 - \sqrt {{T_r}} } \right)} \right)\left( {1 + \kappa \left( {1 - \sqrt {{T_r}} } \right)} \right) = 1 + 2\kappa \left( {1 - \sqrt {{T_r}} } \right) + {\kappa ^2}{\left( {1 - \sqrt {{T_r}} } \right)^2}\\
\alpha  = 1 + 2\kappa \left( {1 - \sqrt {{T_r}} } \right) + {\kappa ^2}{\left( {1 - \sqrt {{T_r}} } \right)^2}\\
\alpha  = 1 + 2\kappa  - 2\kappa \sqrt {{T_r}}  + {\kappa ^2}\left[ {1 - 2\sqrt {{T_r}}  + {T_r}} \right]\\
T = {T_r}{T_c}\\
p = \frac{{R{T_r}{T_c}}}{{v - b}} - \frac{{a\left( {1 + 2\kappa  - 2\kappa \sqrt {{T_r}}  + {\kappa ^2}\left[ {1 - 2\sqrt {{T_r}}  + {T_r}} \right]} \right)}}{{v\left( {v + b} \right)}}\\
\\
{\rm{Factor in terms of }}\sqrt {{T_r}} \\
\\
p = \frac{{R{T_r}{T_c}}}{{v - b}} - \frac{{a\left( {1 + 2\kappa  + {\kappa ^2} - 2\kappa \sqrt {{T_r}}  + {\kappa ^2}\left[ { - 2\sqrt {{T_r}}  + {T_r}} \right]} \right)}}{{v\left( {v + b} \right)}}\\
p = \frac{{R{T_r}{T_c}}}{{v - b}} - \frac{{a\left( {1 + 2\kappa  + {\kappa ^2} - 2\kappa (1 + \kappa )\sqrt {{T_r}}  + {\kappa ^2}{T_r}} \right)}}{{v\left( {v + b} \right)}}\\
p = \frac{{R{T_r}{T_c}}}{{v - b}} - \frac{{a\left( {1 + 2\kappa  + {\kappa ^2}} \right)}}{{v\left( {v + b} \right)}} + \frac{{2a\kappa (1 + \kappa )}}{{v\left( {v + b} \right)}}\sqrt {{T_r}}  - \frac{{a{\kappa ^2}}}{{v\left( {v + b} \right)}}{T_r}\\
0 = \left[ {\frac{{R{T_c}}}{{v - b}} - \frac{{a{\kappa ^2}}}{{v\left( {v + b} \right)}}} \right]{T_r} + \frac{{2a\kappa (1 + \kappa )}}{{v\left( {v + b} \right)}}\sqrt {{T_r}}  - \frac{{a\left( {1 + 2\kappa  + {\kappa ^2}} \right)}}{{v\left( {v + b} \right)}} - p
\end{array}
\f]
 */
double FlashRoutines::T_DP_PengRobinson(HelmholtzEOSMixtureBackend& HEOS, double rhomolar, double p) {
    double omega = NAN, R = NAN, kappa = NAN, a = NAN, b = NAN, A = NAN, B = NAN, C = NAN, Tc = NAN, pc = NAN, V = 1 / rhomolar;
    omega = HEOS.acentric_factor();
    Tc = HEOS.T_critical();
    pc = HEOS.p_critical();
    R = HEOS.gas_constant();

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
    a = 0.457235 * R * R * Tc * Tc / pc;
    b = 0.077796 * R * Tc / pc;
    double den = V * V + 2 * b * V - b * b;

    // A sqrt(Tr)^2 + B sqrt(Tr) + C = 0
    A = R * Tc / (V - b) - a * kappa * kappa / (den);
    B = +2 * a * kappa * (1 + kappa) / (den);
    C = -a * (1 + 2 * kappa + kappa * kappa) / (den)-p;

    //D = B*B-4*A*C;

    double sqrt_Tr1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
    //double sqrt_Tr2 = (-B-sqrt(B*B-4*A*C))/(2*A);
    return sqrt_Tr1 * sqrt_Tr1 * Tc;
};

void FlashRoutines::DP_flash(HelmholtzEOSMixtureBackend& HEOS) {
    // Comment out the check for an imposed phase.  There's no code to handle if it is!
    // Solver below and flash calculations (if two phase) have to be called anyway.
    //
    //  if (HEOS.imposed_phase_index == iphase_not_imposed) // If no phase index is imposed (see set_components function)
    //  {
    if (HEOS.is_pure_or_pseudopure) {
        // Find the phase, while updating all internal variables possible using the pressure
        bool saturation_called = false;
        HEOS.p_phase_determination_pure_or_pseudopure(iDmolar, HEOS._rhomolar, saturation_called);

        if (HEOS.isHomogeneousPhase()) {
            CoolPropDbl T0 = NAN;
            if (HEOS._phase == iphase_liquid) {
                // If it is a liquid, start off at the ancillary value
                if (saturation_called) {
                    T0 = HEOS.SatL->T();
                } else {
                    T0 = HEOS._TLanc.pt();
                }
            } else if (HEOS._phase == iphase_supercritical_liquid) {
                // If it is a supercritical
                T0 = 1.1 * HEOS.T_critical();
            } else if (HEOS._phase == iphase_gas || HEOS._phase == iphase_supercritical_gas || HEOS._phase == iphase_supercritical) {
                // First, get a guess for density from Peng-Robinson
                T0 = T_DP_PengRobinson(HEOS, HEOS.rhomolar(), HEOS.p());
            } else {
                throw ValueError("I should never get here");
            }
            if (!std::isfinite(T0)) {
                throw ValueError("Starting value of T0 is not valid in DP_flash");
            }
            // Then, do the solver using the full EOS
            solver_DP_resid resid(&HEOS, HEOS.rhomolar(), HEOS.p());
            std::string errstr;
            Halley(resid, T0, 1e-10, 100);
            HEOS._Q = -1;
            // Update the state for conditions where the state was guessed
            HEOS.recalculate_singlephase_phase();
            if (!get_config_bool(DONT_CHECK_PROPERTY_LIMITS) && HEOS._T > 1.5 * HEOS.Tmax()) {
                throw CoolProp::OutOfRangeError(format("DP yielded T > 1.5Tmax w/ T (%g) K").c_str());
            }
        } else {
            // Nothing to do here; phase determination has handled this already
        }
    } else {
        throw NotImplementedError("DP_flash not ready for mixtures");
    }
    //  }
    //  TO DO:  Put the imposed phase check back in
    //          and provide the else code here if it is imposed.
}

class DQ_flash_residual : public FuncWrapper1DWithTwoDerivs
{
   public:
    HelmholtzEOSMixtureBackend& HEOS;
    double rhomolar, Q_target;
    DQ_flash_residual(HelmholtzEOSMixtureBackend& HEOS, double rhomolar, double Q_target) : HEOS(HEOS), rhomolar(rhomolar), Q_target(Q_target) {};
    double call(double T) {
        HEOS.update(QT_INPUTS, 0, T);  // Doesn't matter whether liquid or vapor, we are just doing a full VLE call for given T
        double rhoL = HEOS.saturated_liquid_keyed_output(iDmolar);
        double rhoV = HEOS.saturated_vapor_keyed_output(iDmolar);
        /// Error between calculated and target vapor quality based on densities
        return (1 / rhomolar - 1 / rhoL) / (1 / rhoV - 1 / rhoL) - Q_target;
    }
    double deriv(double T) {
        return _HUGE;
    }
    double second_deriv(double T) {
        return _HUGE;
    }
};

// Tolerance on Q to be treated as exactly 0 or exactly 1 for the superancillary
// strict-mode path. Centralized so the various callers stay in lockstep (#2773).
static constexpr double Q_BOUNDARY_TOL = 1e-10;

// Tolerance for deduplicating near-identical T-roots returned by get_x_for_y in
// adjacent monotonic intervals at an extremum (#2773 review I1).
static constexpr double T_ROOT_DEDUP_TOL = 1e-6;

// Sum the two IdealHelmholtzEnthalpyEntropyOffset contributions
// (EnthalpyEntropyOffsetCore from JSON parse and the user-mutable
// EnthalpyEntropyOffset) and apply the alpha0 prefactor. Disabled offsets
// canonicalize to (0, 0). The returned (a1, a2) is what goes into the
// SuperAncillary stamp and into the shift formula (Δh = R·T_red·Δa2,
// Δs = −R·Δa1). See #2773.
std::pair<double, double> FlashRoutines::alpha0_offset_total(HelmholtzEOSMixtureBackend& HEOS) {
    const auto& alpha0 = HEOS.get_components()[0].EOS().alpha0;
    const auto& core = alpha0.EnthalpyEntropyOffsetCore;
    const auto& user = alpha0.EnthalpyEntropyOffset;
    const double prefactor = alpha0.get_prefactor();
    const double a1_sum = (core.is_enabled() ? static_cast<double>(core.get_a1()) : 0.0) + (user.is_enabled() ? static_cast<double>(user.get_a1()) : 0.0);
    const double a2_sum = (core.is_enabled() ? static_cast<double>(core.get_a2()) : 0.0) + (user.is_enabled() ? static_cast<double>(user.get_a2()) : 0.0);
    return {prefactor * a1_sum, prefactor * a2_sum};
}

// Has the fluid got a superancillary AND is the requested Q exactly 0 or 1?
// Used by the no-guess default flashes to decide whether to route through the
// strict-mode superancillary path. Pseudo-pure fluids are rejected here too,
// because their caloric superancillaries are not built (see #2773 review I3).
bool FlashRoutines::sat_superanc_path_applies(HelmholtzEOSMixtureBackend& HEOS) {
    if (!HEOS.is_pure()) return false;
    auto p = HEOS.get_superanc();
    if (!p) return false;
    const double Q = HEOS._Q;
    return std::abs(Q) < Q_BOUNDARY_TOL || std::abs(Q - 1) < Q_BOUNDARY_TOL;
}

void FlashRoutines::DQ_flash(HelmholtzEOSMixtureBackend& HEOS) {
    if (!HEOS.is_pure_or_pseudopure) {
        throw NotImplementedError("DQ_flash not ready for mixtures");
    }
    HEOS.specify_phase(iphase_twophase);
    const double rhomolar = HEOS._rhomolar;
    const double Q = HEOS._Q;
    // Strict-mode superancillary path for Q ∈ {0, 1} on pure fluids: enumerate
    // every T-root and refuse to silently pick when there is more than one.
    // See GitHub #2773 / #2834.
    if (sat_superanc_path_applies(HEOS)) {
        const double T_super = resolve_T_via_superancillary(HEOS, iDmolar, rhomolar, std::nullopt, "DQ_flash");
        HEOS._T = T_super;
        QT_flash(HEOS);
        HEOS._rhomolar = rhomolar;
        HEOS._Q = Q;
        HEOS._phase = iphase_twophase;
        return;
    }
    // Fallback: Brent over [Tmin, Tmax] for fractional Q or fluids without
    // a superancillary. This was the original DQ_flash.
    SaturationSolvers::saturation_PHSU_pure_options options;
    options.use_logdelta = false;
    double Tmax = HEOS.T_critical() - 0.1;
    double Tmin = HEOS.Tmin() + 0.1;
    const double eps = 1e-12;
    if (rhomolar >= (HEOS.rhomolar_critical() + eps) && Q > (0 + eps)) {
        throw CoolProp::OutOfRangeError(
          format("DQ inputs are not defined for density (%g) above critical density (%g) and Q>0", rhomolar, HEOS.rhomolar_critical()).c_str());
    }
    DQ_flash_residual resid(HEOS, rhomolar, Q);
    Brent(resid, Tmin, Tmax, DBL_EPSILON, 1e-10, 100);
    HEOS._p = HEOS.SatV->p();
    HEOS._T = HEOS.SatV->T();
    HEOS._rhomolar = rhomolar;
    HEOS._Q = Q;
    HEOS._phase = iphase_twophase;
}
void FlashRoutines::HQ_flash(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl Tguess) {
    if (!HEOS.is_pure_or_pseudopure) {
        throw NotImplementedError("HQ_flash not ready for mixtures");
    }
    if (std::abs(HEOS.Q() - 1) > 1e-10) {
        throw ValueError(format("non-unity quality not currently allowed for HQ_flash"));
    }
    HEOS.specify_phase(iphase_twophase);
    const double hmolar = HEOS._hmolar;
    // Strict-mode superancillary path. See GitHub #2773 / #2834. Use the same
    // gate as DQ_flash and QS_flash so all three default flashes route
    // consistently (review I3 follow-up).
    if (sat_superanc_path_applies(HEOS)) {
        const double T_super = resolve_T_via_superancillary(HEOS, iHmolar, hmolar, std::nullopt, "HQ_flash");
        HEOS._T = T_super;
        QT_flash(HEOS);
        HEOS._hmolar = hmolar;
        HEOS._phase = iphase_twophase;
        return;
    }
    // Fallback: ancillary-seeded saturation_PHSU_pure for fluids without a
    // superancillary. (Note: ignores Tguess for branch selection — that
    // capability is exposed via update_with_guesses + HQ_flash_with_guesses.)
    (void)Tguess;
    SaturationSolvers::saturation_PHSU_pure_options options;
    options.use_logdelta = false;
    options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_HV;
    SaturationSolvers::saturation_PHSU_pure(HEOS, hmolar, options);
    HEOS._p = HEOS.SatV->p();
    HEOS._T = HEOS.SatV->T();
    HEOS._rhomolar = HEOS.SatV->rhomolar();
    HEOS._phase = iphase_twophase;
}
void FlashRoutines::QS_flash(HelmholtzEOSMixtureBackend& HEOS) {
    if (!HEOS.is_pure_or_pseudopure) {
        throw NotImplementedError("QS_flash not ready for mixtures");
    }
    if (std::abs(HEOS.smolar() - HEOS.get_state("reducing").smolar) < 0.001) {
        HEOS._p = HEOS.p_critical();
        HEOS._T = HEOS.T_critical();
        HEOS._rhomolar = HEOS.rhomolar_critical();
        HEOS._phase = iphase_critical_point;
        return;
    }
    HEOS.specify_phase(iphase_twophase);
    const double smolar = HEOS._smolar;
    // Strict-mode superancillary path. See GitHub #2773 / #2834.
    if (sat_superanc_path_applies(HEOS)) {
        const double T_super = resolve_T_via_superancillary(HEOS, iSmolar, smolar, std::nullopt, "QS_flash");
        HEOS._T = T_super;
        QT_flash(HEOS);
        HEOS._smolar = smolar;
        HEOS._phase = iphase_twophase;
        return;
    }
    // Fallback: ancillary-seeded saturation_PHSU_pure path.
    if (std::abs(HEOS.Q()) < 1e-10) {
        SaturationSolvers::saturation_PHSU_pure_options options;
        options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_SL;
        options.use_logdelta = false;
        SaturationSolvers::saturation_PHSU_pure(HEOS, smolar, options);
        HEOS._p = HEOS.SatL->p();
        HEOS._T = HEOS.SatL->T();
        HEOS._rhomolar = HEOS.SatL->rhomolar();
    } else if (std::abs(HEOS.Q() - 1) < 1e-10) {
        SaturationSolvers::saturation_PHSU_pure_options options;
        options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_SV;
        options.use_logdelta = false;
        SaturationSolvers::saturation_PHSU_pure(HEOS, smolar, options);
        HEOS._p = HEOS.SatV->p();
        HEOS._T = HEOS.SatV->T();
        HEOS._rhomolar = HEOS.SatV->rhomolar();
    } else {
        throw ValueError(format("non-zero or 1 quality not currently allowed for QS_flash"));
    }
    HEOS._phase = iphase_twophase;
}

// Translate a CoolProp parameters enum to the SuperAncillary's char key.
// Supports the three caloric saturation properties used by the #2773 paths.
static char param_to_superanc_key(parameters key) {
    switch (key) {
        case iDmolar:
            return 'D';
        case iHmolar:
            return 'H';
        case iSmolar:
            return 'S';
        default:
            throw ValueError(format("Unsupported parameters key for saturation superancillary: %d", static_cast<int>(key)));
    }
}

// Unified helper for both the no-guess default flashes and the
// *_flash_with_guesses paths (#2773). Looks up the saturation superancillary
// for property `key`, handles caloric lazy-build, and returns the
// disambiguated T-root.
//
// If guess_T is set, picks the monotonic sub-interval whose temperature range
// contains guess_T (with tie-break by midpoint distance for boundary cases)
// and TOMS748-solves there. If guess_T is std::nullopt, enumerates every root,
// dedups near-identical roots produced at extrema by adjacent intervals, and
// either returns the single root, throws MultipleSolutionsError, or throws
// OutOfRangeError. The TOMS748 result is accepted as-is — superancillaries are
// accurate to ~1e-12 relative to the multi-precision EOS reference, well below
// the precision of any downstream EOS evaluation.
double FlashRoutines::resolve_T_via_superancillary(HelmholtzEOSMixtureBackend& HEOS, parameters key, double target_value,
                                                   std::optional<double> guess_T, const char* fn_name) {
    if (guess_T.has_value() && (!ValidNumber(*guess_T) || *guess_T <= 0)) {
        throw ValueError(format("%s requires a positive guess.T", fn_name));
    }
    if (!HEOS.is_pure_or_pseudopure) {
        throw NotImplementedError(format("%s not ready for mixtures", fn_name));
    }
    if (!HEOS.is_pure()) {
        throw NotImplementedError(format("%s: superancillary path is not supported for pseudo-pure fluids", fn_name));
    }
    auto superanc_ptr = HEOS.get_superanc();
    if (!superanc_ptr) {
        throw NotImplementedError(format("%s requires a superancillary; this fluid has none", fn_name));
    }
    const double Q = HEOS._Q;
    short Q_key = -1;
    if (std::abs(Q) < Q_BOUNDARY_TOL) {
        Q_key = 0;
    } else if (std::abs(Q - 1) < Q_BOUNDARY_TOL) {
        Q_key = 1;
    } else {
        throw ValueError(format("%s currently requires Q=0 or Q=1; got Q=%g", fn_name, Q));
    }
    const char k = param_to_superanc_key(key);
    if (k == 'H' || k == 'S') {
        HEOS.ensure_caloric_superancillaries();
    }
    auto& superanc = *superanc_ptr;
    if (!superanc.has_variable(k)) {
        throw NotImplementedError(format("%s: %c-superancillary unavailable", fn_name, k));
    }
    const auto& approx = superanc.get_approx1d(k, Q_key);
    const char* phase_name = (Q_key == 0) ? "liquid" : "vapor";

    // Shift the user's target value into the cache's reference frame for H
    // and S. The IdealHelmholtzEnthalpyEntropyOffset's effect on h(T) and
    // s(T) along the saturation curve is exactly a constant — Δh = R·T_red·Δa2
    // and Δs = −R·Δa1, with (a1, a2) summed across both EnthalpyEntropyOffset
    // and EnthalpyEntropyOffsetCore and scaled by the alpha0 prefactor. So
    // we only need to translate the target, never rebuild. See
    // ensure_HS_under_lock comment and #2773. For D the cache is
    // reference-state-independent intrinsically (no shift needed).
    double target_in_cache_frame = target_value;
    if (k == 'H' || k == 'S') {
        const auto stamp = superanc.get_caloric_alpha0_stamp();
        if (stamp.has_value()) {
            const auto [caller_a1, caller_a2] = alpha0_offset_total(HEOS);
            const double cache_a1 = stamp->first;
            const double cache_a2 = stamp->second;
            if (k == 'H') {
                // h_cache(T) = h_caller(T) + R·T_red·(a2_cache − a2_caller)
                const double R = HEOS.gas_constant();
                const double T_red = HEOS.get_reducing_state().T;
                target_in_cache_frame = target_value + R * T_red * (cache_a2 - caller_a2);
            } else {  // 'S'
                // s_cache(T) = s_caller(T) + (−R)·(a1_cache − a1_caller)
                const double R = HEOS.gas_constant();
                target_in_cache_frame = target_value - R * (cache_a1 - caller_a1);
            }
        }
    }

    if (guess_T.has_value()) {
        // Pick the monotonic sub-interval whose temperature range contains
        // guess_T (and which can reach target_value in its y-range), then
        // TOMS748 the unique root inside it. If guess_T sits exactly on an
        // interval boundary (an extremum), tie-break by midpoint distance
        // toward the interval the user is more plausibly aiming at.
        const double gT = *guess_T;
        const auto& intervals = approx.get_monotonic_intervals();
        const auto& expansions = approx.get_expansions();
        const auto* chosen = [&]() -> const auto* {
            const auto* best_match = static_cast<decltype(intervals.data())>(nullptr);
            double best_midpoint_dist = std::numeric_limits<double>::infinity();
            for (const auto& iv : intervals) {
                if (!iv.contains_y(target_in_cache_frame)) continue;
                if (gT >= iv.xmin && gT <= iv.xmax) {
                    const double mid_dist = std::abs(0.5 * (iv.xmin + iv.xmax) - gT);
                    if (mid_dist < best_midpoint_dist) {
                        best_midpoint_dist = mid_dist;
                        best_match = &iv;
                    }
                }
            }
            return best_match;
        }();
        if (chosen != nullptr) {
            for (const auto& ei : chosen->expansioninfo) {
                if (ei.contains_y(target_in_cache_frame)) {
                    const auto& e = expansions[ei.idx];
                    auto [xvalue, num_steps] = e.solve_for_x_count(target_in_cache_frame, ei.xmin, ei.xmax, 64, 100U, 1e-10);
                    (void)num_steps;
                    return xvalue;
                }
            }
        }
        throw SolutionError(format("%s: no T-root on saturated %s in the monotonic sub-interval containing guess.T=%g K for %c=%g; "
                                   "superancillary range [%g, %g] K",
                                   fn_name, phase_name, gT, k, target_value, approx.xmin(), approx.xmax()));
    }

    // No guess: enumerate every root, dedup near-identical roots from adjacent
    // intervals at an extremum, then strict-mode-or-throw.
    auto solns = approx.get_x_for_y(target_in_cache_frame, 64, 100U, 1e-10);
    if (solns.empty()) {
        throw OutOfRangeError(format("%s: no T-root on saturated %s for %c=%g; superancillary range [%g, %g] K", fn_name, phase_name, k, target_value,
                                     approx.xmin(), approx.xmax()));
    }
    std::sort(solns.begin(), solns.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<double> Ts;
    Ts.reserve(solns.size());
    for (const auto& s : solns) {
        if (Ts.empty() || std::abs(s.first - Ts.back()) > T_ROOT_DEDUP_TOL) {
            Ts.push_back(s.first);
        }
    }
    if (Ts.size() > 1) {
        std::string Ts_str;
        for (std::size_t i = 0; i < Ts.size(); ++i) {
            Ts_str += format("%g K", Ts[i]);
            if (i + 1 < Ts.size()) Ts_str += ", ";
        }
        throw MultipleSolutionsError(
          format("%s: %c=%g on saturated %s has %zu T-roots (%s); use update_with_guesses with guess.T to pick a branch (see GitHub #2773)", fn_name,
                 k, target_value, phase_name, Ts.size(), Ts_str.c_str()));
    }
    return Ts[0];
}

// Branch-disambiguating variant of DQ_flash (see GitHub #2773).
// rho_L(T) on water and D2O is non-monotonic near 4 °C / 11 °C, so the DQ flash
// can have two T-solutions for the same density. Use the rho_sat superancillary
// to pick the monotonic sub-interval whose x-range contains guess.T and TOMS748-solve there.
void FlashRoutines::DQ_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    const double rhomolar = HEOS._rhomolar;
    const double T_super = resolve_T_via_superancillary(HEOS, iDmolar, rhomolar, guess.T, "DQ_flash_with_guesses");
    HEOS._T = T_super;
    HEOS.specify_phase(iphase_twophase);
    // Despite the name, for a pure fluid with superancillaries this is a fast
    // path: QT_flash dispatches to the rho_sat / p_sat superancillary at T —
    // no iterative VLE solve. Sub-microsecond per call.
    QT_flash(HEOS);
    HEOS._rhomolar = rhomolar;
    HEOS._phase = iphase_twophase;
}

// Branch-disambiguating variant of HQ_flash (see GitHub #2773).
// The default HQ_flash routes through saturation_PHSU_pure with IMPOSED_HV,
// which silently picks one of the two T-roots when h_g(T) is non-monotonic
// (water saturated-vapor enthalpy peaks near 540 K). This variant uses the
// h_sat superancillary (built lazily on first use) to disambiguate.
void FlashRoutines::HQ_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    const double hmolar = HEOS._hmolar;
    const double T_super = resolve_T_via_superancillary(HEOS, iHmolar, hmolar, guess.T, "HQ_flash_with_guesses");
    HEOS._T = T_super;
    HEOS.specify_phase(iphase_twophase);
    // Fast path: superancillary eval at T, not an iterative VLE solve. See DQ_flash_with_guesses for details.
    QT_flash(HEOS);
    HEOS._hmolar = hmolar;
    HEOS._phase = iphase_twophase;
}

// Branch-disambiguating variant of QS_flash. Mirrors HQ_flash_with_guesses.
void FlashRoutines::QS_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    const double smolar = HEOS._smolar;
    const double T_super = resolve_T_via_superancillary(HEOS, iSmolar, smolar, guess.T, "QS_flash_with_guesses");
    HEOS._T = T_super;
    HEOS.specify_phase(iphase_twophase);
    // Fast path: superancillary eval at T, not an iterative VLE solve. See DQ_flash_with_guesses for details.
    QT_flash(HEOS);
    HEOS._smolar = smolar;
    HEOS._phase = iphase_twophase;
}

void FlashRoutines::QT_flash(HelmholtzEOSMixtureBackend& HEOS) {
    CoolPropDbl T = HEOS._T;
    CoolPropDbl Q = HEOS._Q;
    if (HEOS.is_pure_or_pseudopure) {

        if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
            auto superanc_ptr = HEOS.get_superanc();
            if (superanc_ptr) {
                auto& superanc = *superanc_ptr;

                CoolPropDbl Tcrit_num = superanc.get_Tcrit_num();
                if (T > Tcrit_num) {
                    throw ValueError(
                      format("Temperature to QT_flash [%0.8Lg K] may not be above the numerical critical point of %0.15Lg K", T, Tcrit_num));
                }
                auto rhoL = superanc.eval_sat(T, 'D', 0);
                auto rhoV = superanc.eval_sat(T, 'D', 1);
                auto p = superanc.eval_sat(T, 'P', 1);
                HEOS.SatL->update_TDmolarP_unchecked(T, rhoL, p);
                HEOS.SatV->update_TDmolarP_unchecked(T, rhoV, p);
                HEOS._p = p;
                HEOS._rhomolar = 1 / (Q / rhoV + (1 - Q) / rhoL);
                HEOS._phase = iphase_twophase;
                return;
            }
        }

        // The maximum possible saturation temperature
        // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
        CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() + 1e-13;

        // Check what the minimum limits for the equation of state are
        CoolPropDbl Tmin_satL = NAN, Tmin_satV = NAN, Tmin_sat = NAN;
        HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
        Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;

        // Get a reference to keep the code a bit cleaner
        const CriticalRegionSplines& splines = HEOS.components[0].EOS().critical_region_splines;

        if ((get_config_bool(CRITICAL_WITHIN_1UK) && std::abs(T - Tmax_sat) < 1e-6) || std::abs(T - Tmax_sat) < 1e-12) {
            // If exactly(ish) at the critical temperature, liquid and vapor have the critial density
            HEOS.SatL->update(DmolarT_INPUTS, HEOS.rhomolar_critical(), HEOS._T);
            HEOS.SatV->update(DmolarT_INPUTS, HEOS.rhomolar_critical(), HEOS._T);
            HEOS._rhomolar = HEOS.rhomolar_critical();
            HEOS._p = 0.5 * HEOS.SatV->p() + 0.5 * HEOS.SatL->p();
        } else if (!is_in_closed_range(Tmin_sat - 0.1, Tmax_sat, T) && (CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS) == false)) {
            throw ValueError(format("Temperature to QT_flash [%0.8Lg K] must be in range [%0.8Lg K, %0.8Lg K]", T, Tmin_sat - 0.1, Tmax_sat));
        } else if (get_config_bool(CRITICAL_SPLINES_ENABLED) && splines.enabled && HEOS._T > splines.T_min) {
            double rhoL = _HUGE, rhoV = _HUGE;
            // Use critical region spline if it has it and temperature is in its range
            splines.get_densities(T, splines.rhomolar_min, HEOS.rhomolar_critical(), splines.rhomolar_max, rhoL, rhoV);
            HEOS.SatL->update(DmolarT_INPUTS, rhoL, HEOS._T);
            HEOS.SatV->update(DmolarT_INPUTS, rhoV, HEOS._T);
            HEOS._p = 0.5 * HEOS.SatV->p() + 0.5 * HEOS.SatL->p();
            HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
        } else if (!(HEOS.components[0].EOS().pseudo_pure)) {
            // Set some input options
            SaturationSolvers::saturation_T_pure_Akasaka_options options(false);

            // Actually call the solver
            SaturationSolvers::saturation_T_pure_Maxwell(HEOS, HEOS._T, options);

            HEOS._p = 0.5 * HEOS.SatV->p() + 0.5 * HEOS.SatL->p();
            HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
        } else {
            // Pseudo-pure fluid
            CoolPropDbl rhoLanc = _HUGE, rhoVanc = _HUGE, rhoLsat = _HUGE, rhoVsat = _HUGE;
            if (std::abs(HEOS._Q) < DBL_EPSILON) {
                HEOS._p = HEOS.components[0].ancillaries.pL.evaluate(HEOS._T);  // These ancillaries are used explicitly
                rhoLanc = HEOS.components[0].ancillaries.rhoL.evaluate(HEOS._T);
                HEOS.SatL->update_TP_guessrho(HEOS._T, HEOS._p, rhoLanc);
                HEOS._rhomolar = HEOS.SatL->rhomolar();
            } else if (std::abs(HEOS._Q - 1) < DBL_EPSILON) {
                HEOS._p = HEOS.components[0].ancillaries.pV.evaluate(HEOS._T);  // These ancillaries are used explicitly
                rhoVanc = HEOS.components[0].ancillaries.rhoV.evaluate(HEOS._T);
                HEOS.SatV->update_TP_guessrho(HEOS._T, HEOS._p, rhoVanc);
                HEOS._rhomolar = HEOS.SatV->rhomolar();
            } else {
                throw CoolProp::ValueError(format("For pseudo-pure fluid, quality must be equal to 0 or 1.  Two-phase quality is not defined"));
            }

            try {
            } catch (...) {
                // Near the critical point, the behavior is not very nice, so we will just use the ancillary
                rhoLsat = rhoLanc;
                rhoVsat = rhoVanc;
            }
        }
        // Load the outputs
        HEOS._phase = iphase_twophase;
    } else {
        if (HEOS.PhaseEnvelope.built) {
            PT_Q_flash_mixtures(HEOS, iT, HEOS._T);
        } else {
            // Set some input options
            SaturationSolvers::mixture_VLE_IO options;
            options.sstype = SaturationSolvers::imposed_T;
            options.Nstep_max = 20;

            // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
            CoolPropDbl pguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._T, SaturationSolvers::imposed_T, HEOS.mole_fractions);

            // Use Wilson iteration to obtain updated guess for pressure
            pguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._T, SaturationSolvers::imposed_T, HEOS.mole_fractions, pguess);

            // Actually call the successive substitution solver
            SaturationSolvers::successive_substitution(HEOS, HEOS._Q, HEOS._T, pguess, HEOS.mole_fractions, HEOS.K, options);

            // -----
            // Newton-Raphson
            // -----

            SaturationSolvers::newton_raphson_saturation NR;
            SaturationSolvers::newton_raphson_saturation_options IO;

            IO.bubble_point = (HEOS._Q < 0.5);

            IO.x = options.x;
            IO.y = options.y;
            IO.rhomolar_liq = options.rhomolar_liq;
            IO.rhomolar_vap = options.rhomolar_vap;
            IO.T = options.T;
            IO.p = options.p;
            IO.Nstep_max = 30;

            IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::T_IMPOSED;

            if (IO.bubble_point) {
                // Compositions are z, z_incipient
                NR.call(HEOS, IO.x, IO.y, IO);
            } else {
                // Compositions are z, z_incipient
                NR.call(HEOS, IO.y, IO.x, IO);
            }

            HEOS._p = IO.p;
            HEOS._rhomolar = 1 / (HEOS._Q / IO.rhomolar_vap + (1 - HEOS._Q) / IO.rhomolar_liq);
        }
        // Load the outputs
        HEOS._phase = iphase_twophase;
        HEOS._p = HEOS.SatV->p();
        HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
        HEOS._T = HEOS.SatL->T();
    }
}

void get_Henrys_coeffs_FP(const std::string& CAS, double& A, double& B, double& C, double& Tmin, double& Tmax) {
    // Coeffs from Fernandez-Prini JPCRD 2003 DOI: 10.1063/1.1564818
    if (CAS == "7440-59-7")  //Helium
    {
        A = -3.52839;
        B = 7.12983;
        C = 4.47770;
        Tmin = 273.21;
        Tmax = 553.18;
    } else if (CAS == "7440-01-9")  // Ne
    {
        A = -3.18301;
        B = 5.31448;
        C = 5.43774;
        Tmin = 273.20;
        Tmax = 543.36;
    } else if (CAS == "7440-37-1")  // Ar
    {
        A = -8.40954;
        B = 4.29587;
        C = 10.52779;
        Tmin = 273.19;
        Tmax = 568.36;
    } else if (CAS == "7439-90-9")  // Kr
    {
        A = -8.97358;
        B = 3.61508;
        C = 11.29963;
        Tmin = 273.19;
        Tmax = 525.56;
    } else if (CAS == "7440-63-3")  // Xe
    {
        A = -14.21635;
        B = 4.00041;
        C = 15.60999;
        Tmin = 273.22;
        Tmax = 574.85;
    } else if (CAS == "1333-74-0")  // H2
    {
        A = -4.73284;
        B = 6.08954;
        C = 6.06066;
        Tmin = 273.15;
        Tmax = 636.09;
    } else if (CAS == "7727-37-9")  // N2
    {
        A = -9.67578;
        B = 4.72162;
        C = 11.70585;
        Tmin = 278.12;
        Tmax = 636.46;
    } else if (CAS == "7782-44-7")  // O2
    {
        A = -9.44833;
        B = 4.43822;
        C = 11.42005;
        Tmin = 274.15;
        Tmax = 616.52;
    } else if (CAS == "630-08-0")  // CO
    {
        A = -10.52862;
        B = 5.13259;
        C = 12.01421;
        Tmin = 278.15;
        Tmax = 588.67;
    } else if (CAS == "124-38-9")  // CO2
    {
        A = -8.55445;
        B = 4.01195;
        C = 9.52345;
        Tmin = 274.19;
        Tmax = 642.66;
    } else if (CAS == "7783-06-4")  // H2S
    {
        A = -4.51499;
        B = 5.23538;
        C = 4.42126;
        Tmin = 273.15;
        Tmax = 533.09;
    } else if (CAS == "74-82-8")  // CH4
    {
        A = -10.44708;
        B = 4.66491;
        C = 12.12986;
        Tmin = 275.46;
        Tmax = 633.11;
    } else if (CAS == "74-84-0")  // C2H6
    {
        A = -19.67563;
        B = 4.51222;
        C = 20.62567;
        Tmin = 275.44;
        Tmax = 473.46;
    } else if (CAS == "2551-62-4")  // SF6
    {
        A = -16.56118;
        B = 2.15289;
        C = 20.35440;
        Tmin = 283.14;
        Tmax = 505.55;
    } else {
        throw ValueError("Bad component in Henry's law constants");
    }
}
void FlashRoutines::PQ_flash(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.is_pure_or_pseudopure) {

        if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
            auto superanc_ptr = HEOS.get_superanc();
            if (superanc_ptr) {
                auto& superanc = *superanc_ptr;
                CoolPropDbl pmax_num = superanc.get_pmax();
                if (HEOS._p > pmax_num) {
                    throw ValueError(
                      format("Pressure to PQ_flash [%0.8Lg Pa] may not be above the numerical critical point of %0.15Lg Pa", HEOS._p, pmax_num));
                }
                auto T = superanc.get_T_from_p(HEOS._p);
                auto rhoL = superanc.eval_sat(T, 'D', 0);
                auto rhoV = superanc.eval_sat(T, 'D', 1);
                auto p = HEOS._p;
                HEOS.SatL->update_TDmolarP_unchecked(T, rhoL, p);
                HEOS.SatV->update_TDmolarP_unchecked(T, rhoV, p);
                HEOS._T = T;
                HEOS._p = p;
                HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
                HEOS._phase = iphase_twophase;
                return;
            }
        }

        if (HEOS.components[0].EOS().pseudo_pure) {
            // It is a pseudo-pure mixture

            HEOS._TLanc = HEOS.components[0].ancillaries.pL.invert(HEOS._p);
            HEOS._TVanc = HEOS.components[0].ancillaries.pV.invert(HEOS._p);
            // Get guesses for the ancillaries for density
            CoolPropDbl rhoL = HEOS.components[0].ancillaries.rhoL.evaluate(HEOS._TLanc);
            CoolPropDbl rhoV = HEOS.components[0].ancillaries.rhoV.evaluate(HEOS._TVanc);
            // Solve for the density
            HEOS.SatL->update_TP_guessrho(HEOS._TLanc, HEOS._p, rhoL);
            HEOS.SatV->update_TP_guessrho(HEOS._TVanc, HEOS._p, rhoV);

            // Load the outputs
            HEOS._phase = iphase_twophase;
            HEOS._p = HEOS._Q * HEOS.SatV->p() + (1 - HEOS._Q) * HEOS.SatL->p();
            HEOS._T = HEOS._Q * HEOS.SatV->T() + (1 - HEOS._Q) * HEOS.SatL->T();
            HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
        } else {
            // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
            CoolPropDbl pmax_sat = HEOS.calc_pmax_sat();

            // Check what the minimum limits for the equation of state are
            CoolPropDbl pmin_satL = NAN, pmin_satV = NAN, pmin_sat = NAN;
            HEOS.calc_pmin_sat(pmin_satL, pmin_satV);
            pmin_sat = std::max(pmin_satL, pmin_satV);

            // Check for being AT the critical point
            if (is_in_closed_range(pmax_sat * (1 - 1e-10), pmax_sat * (1 + 1e-10), static_cast<CoolPropDbl>(HEOS._p))) {
                // Load the outputs
                HEOS._phase = iphase_critical_point;
                HEOS._p = HEOS.p_critical();
                HEOS._rhomolar = HEOS.rhomolar_critical();
                HEOS._T = HEOS.T_critical();
                return;
            }

            // Check limits
            if (CoolProp::get_config_bool(DONT_CHECK_PROPERTY_LIMITS) == false) {
                if (!is_in_closed_range(pmin_sat * 0.999999, pmax_sat * 1.000001, static_cast<CoolPropDbl>(HEOS._p))) {
                    throw ValueError(format("Pressure to PQ_flash [%6g Pa] must be in range [%8Lg Pa, %8Lg Pa]", HEOS._p, pmin_sat, pmax_sat));
                }
            }
            // ------------------
            // It is a pure fluid
            // ------------------

            // Set some input options
            SaturationSolvers::saturation_PHSU_pure_options options;
            // Specified variable is pressure
            options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_PL;
            // Use logarithm of delta as independent variables
            options.use_logdelta = false;

            double increment = 0.4;

            try {
                for (double omega = 1.0; omega > 0; omega -= increment) {
                    try {
                        options.omega = omega;

                        // Actually call the solver
                        SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS._p, options);

                        // If you get here, there was no error, all is well
                        break;
                    } catch (...) {
                        if (omega < 1.1 * increment) {
                            throw;
                        }
                        // else we are going to try again with a smaller omega
                    }
                }
            } catch (...) {
                // We may need to polish the solution at low pressure
                SaturationSolvers::saturation_P_pure_1D_T(HEOS, HEOS._p, options);
            }

            // Load the outputs
            HEOS._phase = iphase_twophase;
            HEOS._p = HEOS._Q * HEOS.SatV->p() + (1 - HEOS._Q) * HEOS.SatL->p();
            HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
            HEOS._T = HEOS.SatL->T();
        }
    } else {
        if (HEOS.PhaseEnvelope.built) {
            PT_Q_flash_mixtures(HEOS, iP, HEOS._p);
        } else {

            // Set some input options
            SaturationSolvers::mixture_VLE_IO io;
            io.sstype = SaturationSolvers::imposed_p;
            io.Nstep_max = 10;

            // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
            CoolPropDbl Tguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

            // Use Wilson iteration to obtain updated guess for temperature
            Tguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);

            std::vector<CoolPropDbl> K = HEOS.K;

            if (get_config_bool(HENRYS_LAW_TO_GENERATE_VLE_GUESSES) && std::abs(HEOS._Q - 1) < 1e-10) {
                const std::vector<CoolPropFluid>& components = HEOS.get_components();
                std::size_t iWater = 0;
                double p1star = PropsSI("P", "T", Tguess, "Q", 1, "Water");
                const std::vector<CoolPropDbl> y = HEOS.mole_fractions;
                std::vector<CoolPropDbl> x(y.size());
                for (std::size_t i = 0; i < components.size(); ++i) {
                    if (components[i].CAS == "7732-18-5") {
                        iWater = i;
                        continue;
                    } else {
                        double A = NAN, B = NAN, C = NAN, Tmin = NAN, Tmax = NAN;
                        get_Henrys_coeffs_FP(components[i].CAS, A, B, C, Tmin, Tmax);
                        double T_R = Tguess / 647.096, tau = 1 - T_R;
                        double k_H = p1star * exp(A / T_R + B * pow(tau, 0.355) / T_R + C * pow(T_R, -0.41) * exp(tau));
                        x[i] = y[i] * HEOS._p / k_H;
                        //
                        K[i] = y[i] / x[i];
                    }
                }
                // Update water K factor
                double summer = 0;
                for (std::size_t i = 0; i < y.size(); ++i) {
                    if (i != iWater) {
                        summer += x[i];
                    }
                }
                x[iWater] = summer;
                K[iWater] = y[iWater] / x[iWater];
            }

            // Actually call the successive substitution solver
            SaturationSolvers::successive_substitution(HEOS, HEOS._Q, Tguess, HEOS._p, HEOS.mole_fractions, K, io);

            // -----
            // Newton-Raphson
            // -----

            SaturationSolvers::newton_raphson_saturation NR;
            SaturationSolvers::newton_raphson_saturation_options IO;

            IO.bubble_point = (HEOS._Q < 0.5);
            IO.x = io.x;
            IO.y = io.y;
            IO.rhomolar_liq = io.rhomolar_liq;
            IO.rhomolar_vap = io.rhomolar_vap;
            IO.T = io.T;
            IO.p = io.p;
            IO.Nstep_max = 30;
            IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;

            if (IO.bubble_point) {
                // Compositions are z, z_incipient
                NR.call(HEOS, IO.x, IO.y, IO);
            } else {
                // Compositions are z, z_incipient
                NR.call(HEOS, IO.y, IO.x, IO);
            }
        }

        // Load the outputs
        HEOS._phase = iphase_twophase;
        HEOS._p = HEOS.SatV->p();
        HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
        HEOS._T = HEOS.SatL->T();
    }
}

void FlashRoutines::PQ_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    IO.rhomolar_liq = guess.rhomolar_liq;
    IO.rhomolar_vap = guess.rhomolar_vap;
    IO.x = std::vector<CoolPropDbl>(guess.x.begin(), guess.x.end());
    IO.y = std::vector<CoolPropDbl>(guess.y.begin(), guess.y.end());
    IO.T = guess.T;
    IO.p = HEOS._p;
    IO.bubble_point = false;
    IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;

    if (std::abs(HEOS.Q()) < 1e-10) {
        IO.bubble_point = true;
        NR.call(HEOS, IO.x, IO.y, IO);
    } else if (std::abs(HEOS.Q() - 1) < 1e-10) {
        IO.bubble_point = false;
        NR.call(HEOS, IO.y, IO.x, IO);
    } else {
        throw ValueError(format("Quality must be 0 or 1"));
    }

    // Load the other outputs
    HEOS._phase = iphase_twophase;
    HEOS._rhomolar = 1 / (HEOS._Q / IO.rhomolar_vap + (1 - HEOS._Q) / IO.rhomolar_liq);
    HEOS._T = IO.T;
}
void FlashRoutines::QT_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    IO.rhomolar_liq = guess.rhomolar_liq;
    IO.rhomolar_vap = guess.rhomolar_vap;
    IO.x = std::vector<CoolPropDbl>(guess.x.begin(), guess.x.end());
    IO.y = std::vector<CoolPropDbl>(guess.y.begin(), guess.y.end());
    IO.T = HEOS._T;
    IO.p = guess.p;
    IO.bubble_point = false;
    IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::T_IMPOSED;

    if (get_debug_level() > 9) {
        std::cout << format(" QT w/ guess  p %g T %g dl %g dv %g x %s y %s\n", IO.p, IO.T, IO.rhomolar_liq, IO.rhomolar_vap,
                            vec_to_string(IO.x, "%g").c_str(), vec_to_string(IO.y, "%g").c_str());
    }

    if (std::abs(HEOS.Q()) < 1e-10) {
        IO.bubble_point = true;
        NR.call(HEOS, IO.x, IO.y, IO);
    } else if (std::abs(HEOS.Q() - 1) < 1e-10) {
        IO.bubble_point = false;
        NR.call(HEOS, IO.y, IO.x, IO);
    } else {
        throw ValueError(format("Quality must be 0 or 1"));
    }

    // Load the other outputs
    HEOS._p = IO.p;
    HEOS._phase = iphase_twophase;
    HEOS._rhomolar = 1 / (HEOS._Q / IO.rhomolar_vap + (1 - HEOS._Q) / IO.rhomolar_liq);
}

void FlashRoutines::PT_flash_with_guesses(HelmholtzEOSMixtureBackend& HEOS, const GuessesStructure& guess) {
    HEOS.solver_rho_Tp(HEOS.T(), HEOS.p(), guess.rhomolar);
    // Load the other outputs
    HEOS._phase = iphase_gas;  // Guessed for mixtures
    if (HEOS.is_pure_or_pseudopure) {
        if (HEOS._p > HEOS.p_critical()) {
            if (HEOS._T > HEOS.T_critical()) {
                HEOS._phase = iphase_supercritical;
            } else {
                HEOS._phase = iphase_supercritical_liquid;
            }
        } else {
            if (HEOS._T > HEOS.T_critical()) {
                HEOS._phase = iphase_supercritical_gas;
            } else if (HEOS._rhomolar > HEOS.rhomolar_critical()) {
                HEOS._phase = iphase_liquid;
            } else {
                HEOS._phase = iphase_gas;
            }
        }
    }

    HEOS._Q = -1;
}

void FlashRoutines::PT_Q_flash_mixtures(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl value) {

    // Find the intersections in the phase envelope
    std::vector<std::pair<std::size_t, std::size_t>> intersections =
      PhaseEnvelopeRoutines::find_intersections(HEOS.get_phase_envelope_data(), other, value);

    PhaseEnvelopeData& env = HEOS.PhaseEnvelope;

    enum quality_options
    {
        SATURATED_LIQUID,
        SATURATED_VAPOR,
        TWO_PHASE
    };
    quality_options quality;
    if (std::abs(HEOS._Q) < 100 * DBL_EPSILON) {
        quality = SATURATED_LIQUID;
    } else if (std::abs(HEOS._Q - 1) < 100 * DBL_EPSILON) {
        quality = SATURATED_VAPOR;
    } else if (HEOS._Q > 0 && HEOS._Q < 1) {
        quality = TWO_PHASE;
    } else {
        throw ValueError("Quality is not within 0 and 1");
    }

    if (quality == SATURATED_LIQUID || quality == SATURATED_VAPOR) {
        // *********************************************************
        //            Bubble- or dew-point calculation
        // *********************************************************
        // Find the correct solution
        std::vector<std::size_t> solutions;
        for (std::vector<std::pair<std::size_t, std::size_t>>::const_iterator it = intersections.begin(); it != intersections.end(); ++it) {
            if (std::abs(env.Q[it->first] - HEOS._Q) < 10 * DBL_EPSILON && std::abs(env.Q[it->second] - HEOS._Q) < 10 * DBL_EPSILON) {
                solutions.push_back(it->first);
            }
        }

        if (solutions.size() == 1) {

            std::size_t& imax = solutions[0];

            // Shift the solution if needed to ensure that imax+2 and imax-1 are both in range
            if (imax + 2 >= env.T.size()) {
                imax--;
            } else if (imax == 0) {
                imax++;
            }
            // Here imax+2 or imax-1 is still possibly out of range:
            // 1. If imax initially is 1, and env.T.size() <= 3, then imax will become 0.
            // 2. If imax initially is 0, and env.T.size() <= 2, then imax will become MAX_UINT.
            // 3. If imax+2 initially is more than env.T.size(), then single decrement will not bring it to range

            SaturationSolvers::newton_raphson_saturation NR;
            SaturationSolvers::newton_raphson_saturation_options IO;

            if (other == iP) {
                IO.p = HEOS._p;
                IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;
                // p -> rhomolar_vap
                IO.rhomolar_vap = CubicInterp(env.p, env.rhomolar_vap, imax - 1, imax, imax + 1, imax + 2, static_cast<CoolPropDbl>(IO.p));
                IO.T = CubicInterp(env.rhomolar_vap, env.T, imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
            } else if (other == iT) {
                IO.T = HEOS._T;
                IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::T_IMPOSED;
                // T -> rhomolar_vap
                IO.rhomolar_vap = CubicInterp(env.T, env.rhomolar_vap, imax - 1, imax, imax + 1, imax + 2, static_cast<CoolPropDbl>(IO.T));
                IO.p = CubicInterp(env.rhomolar_vap, env.p, imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
            } else {
                throw ValueError();
            }
            IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);

            if (quality == SATURATED_VAPOR) {
                IO.bubble_point = false;
                IO.y = HEOS.get_mole_fractions();  // Because Q = 1
                IO.x.resize(IO.y.size());
                for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
                {
                    IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
                }
                IO.x[IO.x.size() - 1] = 1 - std::accumulate(IO.x.begin(), IO.x.end() - 1, 0.0);
                NR.call(HEOS, IO.y, IO.x, IO);
            } else {
                IO.bubble_point = true;
                IO.x = HEOS.get_mole_fractions();  // Because Q = 0
                IO.y.resize(IO.x.size());
                // Phases are inverted, so "liquid" is actually the lighter phase
                std::swap(IO.rhomolar_liq, IO.rhomolar_vap);
                for (std::size_t i = 0; i < IO.y.size() - 1; ++i)  // First N-1 elements
                {
                    // Phases are inverted, so liquid mole fraction (x) of phase envelope is actually the vapor phase mole fraction
                    // Use the liquid density as well
                    IO.y[i] = CubicInterp(env.rhomolar_vap, env.x[i], imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_liq);
                }
                IO.y[IO.y.size() - 1] = 1 - std::accumulate(IO.y.begin(), IO.y.end() - 1, 0.0);
                NR.call(HEOS, IO.x, IO.y, IO);
            }
        } else if (solutions.size() == 0) {
            throw ValueError("No solution was found in PQ_flash");
        } else {
            throw ValueError("More than 1 solution was found in PQ_flash");
        }
    } else {
        // *********************************************************
        //      Two-phase calculation for given vapor quality
        // *********************************************************

        // Find the correct solution
        std::vector<std::size_t> liquid_solutions, vapor_solutions;
        for (std::vector<std::pair<std::size_t, std::size_t>>::const_iterator it = intersections.begin(); it != intersections.end(); ++it) {
            if (std::abs(env.Q[it->first] - 0) < 10 * DBL_EPSILON && std::abs(env.Q[it->second] - 0) < 10 * DBL_EPSILON) {
                liquid_solutions.push_back(it->first);
            }
            if (std::abs(env.Q[it->first] - 1) < 10 * DBL_EPSILON && std::abs(env.Q[it->second] - 1) < 10 * DBL_EPSILON) {
                vapor_solutions.push_back(it->first);
            }
        }

        if (liquid_solutions.size() != 1 || vapor_solutions.size() != 1) {
            throw ValueError(format("Number liquid solutions [%d] or vapor solutions [%d] != 1", liquid_solutions.size(), vapor_solutions.size()));
        }
        std::size_t iliq = liquid_solutions[0], ivap = vapor_solutions[0];

        SaturationSolvers::newton_raphson_twophase NR;
        SaturationSolvers::newton_raphson_twophase_options IO;
        IO.beta = HEOS._Q;

        CoolPropDbl rhomolar_vap_sat_vap = NAN, T_sat_vap = NAN, rhomolar_liq_sat_vap = NAN, rhomolar_liq_sat_liq = NAN, T_sat_liq = NAN,
                    rhomolar_vap_sat_liq = NAN, p_sat_liq = NAN, p_sat_vap = NAN;

        if (other == iP) {
            IO.p = HEOS._p;
            p_sat_liq = IO.p;
            p_sat_vap = IO.p;
            IO.imposed_variable = SaturationSolvers::newton_raphson_twophase_options::P_IMPOSED;

            // Calculate the interpolated values for beta = 0 and beta = 1
            rhomolar_vap_sat_vap = CubicInterp(env.p, env.rhomolar_vap, ivap - 1, ivap, ivap + 1, ivap + 2, static_cast<CoolPropDbl>(IO.p));
            T_sat_vap = CubicInterp(env.rhomolar_vap, env.T, ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);
            rhomolar_liq_sat_vap = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);

            // Phase inversion for liquid solution (liquid is vapor and vice versa)
            rhomolar_liq_sat_liq = CubicInterp(env.p, env.rhomolar_vap, iliq - 1, iliq, iliq + 1, iliq + 2, static_cast<CoolPropDbl>(IO.p));
            T_sat_liq = CubicInterp(env.rhomolar_vap, env.T, iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);
            rhomolar_vap_sat_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);
        } else if (other == iT) {
            IO.T = HEOS._T;
            T_sat_liq = IO.T;
            T_sat_vap = IO.T;
            IO.imposed_variable = SaturationSolvers::newton_raphson_twophase_options::T_IMPOSED;

            // Calculate the interpolated values for beta = 0 and beta = 1
            rhomolar_vap_sat_vap = CubicInterp(env.T, env.rhomolar_vap, ivap - 1, ivap, ivap + 1, ivap + 2, static_cast<CoolPropDbl>(IO.T));
            p_sat_vap = CubicInterp(env.rhomolar_vap, env.p, ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);
            rhomolar_liq_sat_vap = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);

            // Phase inversion for liquid solution (liquid is vapor and vice versa)
            rhomolar_liq_sat_liq = CubicInterp(env.T, env.rhomolar_vap, iliq - 1, iliq, iliq + 1, iliq + 2, static_cast<CoolPropDbl>(IO.T));
            p_sat_liq = CubicInterp(env.rhomolar_vap, env.p, iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);
            rhomolar_vap_sat_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);
        } else {
            throw ValueError();
        }

        // Weight the guesses by the vapor mole fraction
        IO.rhomolar_vap = IO.beta * rhomolar_vap_sat_vap + (1 - IO.beta) * rhomolar_vap_sat_liq;
        IO.rhomolar_liq = IO.beta * rhomolar_liq_sat_vap + (1 - IO.beta) * rhomolar_liq_sat_liq;
        IO.T = IO.beta * T_sat_vap + (1 - IO.beta) * T_sat_liq;
        IO.p = IO.beta * p_sat_vap + (1 - IO.beta) * p_sat_liq;

        IO.z = HEOS.get_mole_fractions();
        IO.x.resize(IO.z.size());
        IO.y.resize(IO.z.size());

        for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
        {
            CoolPropDbl x_sat_vap = CubicInterp(env.rhomolar_vap, env.x[i], ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);
            CoolPropDbl y_sat_vap = CubicInterp(env.rhomolar_vap, env.y[i], ivap - 1, ivap, ivap + 1, ivap + 2, rhomolar_vap_sat_vap);

            CoolPropDbl x_sat_liq = CubicInterp(env.rhomolar_vap, env.y[i], iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);
            CoolPropDbl y_sat_liq = CubicInterp(env.rhomolar_vap, env.x[i], iliq - 1, iliq, iliq + 1, iliq + 2, rhomolar_liq_sat_liq);

            IO.x[i] = IO.beta * x_sat_vap + (1 - IO.beta) * x_sat_liq;
            IO.y[i] = IO.beta * y_sat_vap + (1 - IO.beta) * y_sat_liq;
        }
        IO.x[IO.x.size() - 1] = 1 - std::accumulate(IO.x.begin(), IO.x.end() - 1, 0.0);
        IO.y[IO.y.size() - 1] = 1 - std::accumulate(IO.y.begin(), IO.y.end() - 1, 0.0);
        NR.call(HEOS, IO);
    }
}
void FlashRoutines::HSU_D_flash_twophase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value) {
    class Residual : public FuncWrapper1D
    {

       public:
        HelmholtzEOSMixtureBackend& HEOS;
        CoolPropDbl rhomolar_spec;  // Specified value for density
        parameters other;           // Key for other value
        CoolPropDbl value;          // value for S,H,U
        CoolPropDbl Qd;             // Quality from density
        Residual(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value)
          : HEOS(HEOS), rhomolar_spec(rhomolar_spec), other(other), value(value) {
            Qd = _HUGE;
        };
        double call(double T) {
            HEOS.update(QT_INPUTS, 0, T);
            HelmholtzEOSMixtureBackend &SatL = HEOS.get_SatL(), &SatV = HEOS.get_SatV();
            // Quality from density
            Qd = (1 / rhomolar_spec - 1 / SatL.rhomolar()) / (1 / SatV.rhomolar() - 1 / SatL.rhomolar());
            // Quality from other parameter (H,S,U)
            CoolPropDbl Qo = (value - SatL.keyed_output(other)) / (SatV.keyed_output(other) - SatL.keyed_output(other));
            // Residual is the difference between the two
            return Qo - Qd;
        }
    } resid(HEOS, rhomolar_spec, other, value);

    // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
    CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() - 1e-13;

    // Check what the minimum limits for the equation of state are
    CoolPropDbl Tmin_satL = NAN, Tmin_satV = NAN, Tmin_sat = NAN;
    HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
    Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;

    Brent(resid, Tmin_sat, Tmax_sat - 0.01, DBL_EPSILON, 1e-12, 20);
    // Solve once more with the final vapor quality
    HEOS.update(QT_INPUTS, resid.Qd, HEOS.T());
}
// D given and one of P,H,S,U
void FlashRoutines::HSU_D_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other) {
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1DWithTwoDerivs
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl rhomolar, value;
        parameters other;
        CoolPropDbl Tmin, Tmax;

        solver_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl rhomolar, CoolPropDbl value, parameters other, CoolPropDbl Tmin, CoolPropDbl Tmax)
          : HEOS(HEOS), rhomolar(rhomolar), value(value), other(other), Tmin(Tmin), Tmax(Tmax) {
            /// Something homogeneous to avoid flash calls
            HEOS->specify_phase(iphase_gas);
        };
        double call(double T) {
            HEOS->update_DmolarT_direct(rhomolar, T);
            double eos = HEOS->keyed_output(other);
            if (other == iP) {
                // For p, should use fractional error
                return (eos - value) / value;
            } else {
                // For everything else, use absolute error
                return eos - value;
            }
        };
        double deriv(double T) {
            if (other == iP) {
                return HEOS->first_partial_deriv(other, iT, iDmolar) / value;
            }
            return HEOS->first_partial_deriv(other, iT, iDmolar);
        };
        double second_deriv(double T) {
            if (other == iP) {
                return HEOS->second_partial_deriv(other, iT, iDmolar, iT, iDmolar) / value;
            }
            return HEOS->second_partial_deriv(other, iT, iDmolar, iT, iDmolar);
        };
        bool input_not_in_range(double T) {
            return (T < Tmin || T > Tmax);
        }
    };

    if (HEOS.is_pure_or_pseudopure) {
        CoolPropFluid& component = HEOS.components[0];

        shared_ptr<HelmholtzEOSMixtureBackend> Sat;
        CoolPropDbl rhoLtriple = component.triple_liquid.rhomolar;
        CoolPropDbl rhoVtriple = component.triple_vapor.rhomolar;
        // Check if in the "normal" region
        if (HEOS._rhomolar >= rhoVtriple && HEOS._rhomolar <= rhoLtriple) {
            CoolPropDbl yL = NAN, yV = NAN, value = NAN, y_solid = NAN;
            CoolPropDbl TLtriple = component.triple_liquid.T;  ///TODO: separate TL and TV for ppure
            CoolPropDbl TVtriple = component.triple_vapor.T;

            // First check if solid (below the line connecting the triple point values) - this is an error for now
            switch (other) {
                case iSmolar:
                    yL = HEOS.calc_smolar_nocache(TLtriple, rhoLtriple);
                    yV = HEOS.calc_smolar_nocache(TVtriple, rhoVtriple);
                    value = HEOS._smolar;
                    break;
                case iHmolar:
                    yL = HEOS.calc_hmolar_nocache(TLtriple, rhoLtriple);
                    yV = HEOS.calc_hmolar_nocache(TVtriple, rhoVtriple);
                    value = HEOS._hmolar;
                    break;
                case iUmolar:
                    yL = HEOS.calc_umolar_nocache(TLtriple, rhoLtriple);
                    yV = HEOS.calc_umolar_nocache(TVtriple, rhoVtriple);
                    value = HEOS._umolar;
                    break;
                case iP:
                    yL = HEOS.calc_pressure_nocache(TLtriple, rhoLtriple);
                    yV = HEOS.calc_pressure_nocache(TVtriple, rhoVtriple);
                    value = HEOS._p;
                    break;
                default:
                    throw ValueError(format("Input is invalid"));
            }
            y_solid = (yV - yL) / (1 / rhoVtriple - 1 / rhoLtriple) * (1 / HEOS._rhomolar - 1 / rhoLtriple) + yL;

            if (value < y_solid) {
                throw ValueError(format("Other input [%d:%g] is solid", other, value));
            }

            // Check if other is above the saturation value.
            SaturationSolvers::saturation_D_pure_options optionsD;
            optionsD.omega = 1;
            optionsD.use_logdelta = false;
            optionsD.max_iterations = 200;
            for (int i_try = 0; i_try < 7; i_try++) {
                try {
                    if (HEOS._rhomolar > HEOS.rhomolar_critical()) {
                        optionsD.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOL;
                        SaturationSolvers::saturation_D_pure(HEOS, HEOS._rhomolar, optionsD);
                        // SatL and SatV have the saturation values
                        Sat = HEOS.SatL;
                    } else {
                        optionsD.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOV;
                        SaturationSolvers::saturation_D_pure(HEOS, HEOS._rhomolar, optionsD);
                        // SatL and SatV have the saturation values
                        Sat = HEOS.SatV;
                    }
                    break;  // good solve
                } catch (CoolPropBaseError) {
                    optionsD.omega /= 2;
                    optionsD.max_iterations *= 2;
                    if (i_try >= 6) {
                        throw;
                    }
                }
            }

            // If it is above, it is not two-phase and either liquid, vapor or supercritical
            if (value > Sat->keyed_output(other)) {
                solver_resid resid(&HEOS, HEOS._rhomolar, value, other, Sat->keyed_output(iT), HEOS.Tmax() * 1.5);
                CoolPropDbl T_converged;
                try {
                    T_converged = Halley(resid, 0.5 * (Sat->keyed_output(iT) + HEOS.Tmax() * 1.5), 1e-10, 100);
                } catch (...) {
                    T_converged = Brent(resid, Sat->keyed_output(iT), HEOS.Tmax() * 1.5, DBL_EPSILON, 1e-12, 100);
                }
                HEOS.unspecify_phase();
                // Re-evaluate state at the converged (rho, T) so the alphar
                // cache matches the final temperature, not the last solver
                // iterate. Without this the next h/s/p query returns values
                // computed from stale alphar derivatives (#1907).
                HEOS.update_DmolarT_direct(HEOS._rhomolar, T_converged);
                HEOS._Q = 10000;
                // Update the phase flag
                HEOS.recalculate_singlephase_phase();
            } else {
                // Now we know that temperature is between Tsat(D) +- tolerance and the minimum temperature for the fluid
                if (other == iP) {
                    // Iterate to find T(p), its just a saturation call

                    // Set some input options
                    SaturationSolvers::saturation_PHSU_pure_options optionsPHSU;
                    // Specified variable is pressure
                    optionsPHSU.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_PL;
                    // Use logarithm of delta as independent variables
                    optionsPHSU.use_logdelta = false;

                    // Actually call the solver
                    SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS._p, optionsPHSU);

                    // Load the outputs
                    HEOS._phase = iphase_twophase;
                    HEOS._Q = (1 / HEOS._rhomolar - 1 / HEOS.SatL->rhomolar()) / (1 / HEOS.SatV->rhomolar() - 1 / HEOS.SatL->rhomolar());
                    HEOS._T = HEOS.SatL->T();
                } else {
                    // Residual is difference in quality calculated from density and quality calculated from the other parameter
                    // Iterate to find T
                    HSU_D_flash_twophase(HEOS, HEOS._rhomolar, other, value);
                    HEOS._phase = iphase_twophase;
                }
            }
        }
        // Check if vapor/solid region below triple point vapor density
        else if (HEOS._rhomolar < component.triple_vapor.rhomolar) {
            CoolPropDbl y = NAN, value = NAN;
            CoolPropDbl TVtriple = component.triple_vapor.T;  //TODO: separate TL and TV for ppure

            // If value is above the value calculated from X(Ttriple, _rhomolar), it is vapor
            switch (other) {
                case iSmolar:
                    y = HEOS.calc_smolar_nocache(TVtriple, HEOS._rhomolar);
                    value = HEOS._smolar;
                    break;
                case iHmolar:
                    y = HEOS.calc_hmolar_nocache(TVtriple, HEOS._rhomolar);
                    value = HEOS._hmolar;
                    break;
                case iUmolar:
                    y = HEOS.calc_umolar_nocache(TVtriple, HEOS._rhomolar);
                    value = HEOS._umolar;
                    break;
                case iP:
                    y = HEOS.calc_pressure_nocache(TVtriple, HEOS._rhomolar);
                    value = HEOS._p;
                    break;
                default:
                    throw ValueError(format("Input is invalid"));
            }
            if (value > y) {
                solver_resid resid(&HEOS, HEOS._rhomolar, value, other, TVtriple, HEOS.Tmax() * 1.5);
                HEOS._phase = iphase_gas;
                CoolPropDbl T_converged;
                try {
                    T_converged = Halley(resid, 0.5 * (TVtriple + HEOS.Tmax() * 1.5), DBL_EPSILON, 100);
                } catch (...) {
                    T_converged = Brent(resid, TVtriple, HEOS.Tmax() * 1.5, DBL_EPSILON, 1e-12, 100);
                }
                // Re-evaluate at converged (rho, T) so alphar cache is
                // consistent with the final state (#1907).
                HEOS.unspecify_phase();
                HEOS.update_DmolarT_direct(HEOS._rhomolar, T_converged);
                HEOS._phase = iphase_gas;
                HEOS._Q = 10000;
            } else {
                throw ValueError(format("D < DLtriple %g %g", value, y));
            }

        }
        // Check in the liquid/solid region above the triple point density
        else {
            CoolPropDbl y = NAN, value = NAN;
            CoolPropDbl TLtriple = component.EOS().Ttriple;

            // If value is above the value calculated from X(Ttriple, _rhomolar), it is vapor
            switch (other) {
                case iSmolar:
                    y = HEOS.calc_smolar_nocache(TLtriple, HEOS._rhomolar);
                    value = HEOS._smolar;
                    break;
                case iHmolar:
                    y = HEOS.calc_hmolar_nocache(TLtriple, HEOS._rhomolar);
                    value = HEOS._hmolar;
                    break;
                case iUmolar:
                    y = HEOS.calc_umolar_nocache(TLtriple, HEOS._rhomolar);
                    value = HEOS._umolar;
                    break;
                case iP:
                    y = HEOS.calc_pressure_nocache(TLtriple, HEOS._rhomolar);
                    value = HEOS._p;
                    break;
                default:
                    throw ValueError(format("Input is invalid"));
            }
            if (value > y) {
                solver_resid resid(&HEOS, HEOS._rhomolar, value, other, TLtriple, HEOS.Tmax() * 1.5);
                HEOS._phase = iphase_liquid;
                CoolPropDbl T_converged;
                try {
                    T_converged = Halley(resid, 0.5 * (TLtriple + HEOS.Tmax() * 1.5), DBL_EPSILON, 100);
                } catch (...) {
                    T_converged = Brent(resid, TLtriple, HEOS.Tmax() * 1.5, DBL_EPSILON, 1e-12, 100);
                }
                // Re-evaluate at converged (rho, T) so alphar cache is
                // consistent with the final state (#1907).
                HEOS.unspecify_phase();
                HEOS.update_DmolarT_direct(HEOS._rhomolar, T_converged);
                HEOS._phase = iphase_liquid;
                HEOS._Q = 10000;
            } else {
                throw ValueError(format("D < DLtriple %g %g", value, y));
            }
        }
        // Update the state for conditions where the state was guessed
        if (HEOS.phase() != iphase_twophase) {
            HEOS.recalculate_singlephase_phase();
        }
    } else
        throw NotImplementedError("PHSU_D_flash not ready for mixtures");
}

void FlashRoutines::HSU_P_flash_singlephase_Newton(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl T0, CoolPropDbl rhomolar0) {
    double A[2][2], B[2][2];
    CoolPropDbl y = _HUGE;
    HelmholtzEOSMixtureBackend _HEOS(HEOS.get_components());
    _HEOS.update(DmolarT_INPUTS, rhomolar0, T0);
    CoolPropDbl Tc = HEOS.calc_T_critical();
    CoolPropDbl rhoc = HEOS.calc_rhomolar_critical();
    CoolPropDbl R = HEOS.gas_constant();
    CoolPropDbl p = HEOS.p();
    switch (other) {
        case iHmolar:
            y = HEOS.hmolar();
            break;
        case iSmolar:
            y = HEOS.smolar();
            break;
        default:
            throw ValueError("other is invalid in HSU_P_flash_singlephase_Newton");
    }

    CoolPropDbl worst_error = 999;
    int iter = 0;
    bool failed = false;
    CoolPropDbl omega = 1.0, f2 = NAN, df2_dtau = NAN, df2_ddelta = NAN;
    CoolPropDbl tau = _HEOS.tau(), delta = _HEOS.delta();
    while (worst_error > 1e-6 && failed == false) {

        // All the required partial derivatives
        CoolPropDbl a0 = _HEOS.calc_alpha0_deriv_nocache(0, 0, HEOS.mole_fractions, tau, delta, Tc, rhoc);
        CoolPropDbl da0_ddelta = _HEOS.calc_alpha0_deriv_nocache(0, 1, HEOS.mole_fractions, tau, delta, Tc, rhoc);
        CoolPropDbl da0_dtau = _HEOS.calc_alpha0_deriv_nocache(1, 0, HEOS.mole_fractions, tau, delta, Tc, rhoc);
        CoolPropDbl d2a0_dtau2 = _HEOS.calc_alpha0_deriv_nocache(2, 0, HEOS.mole_fractions, tau, delta, Tc, rhoc);
        CoolPropDbl d2a0_ddelta_dtau = 0.0;

        CoolPropDbl ar = _HEOS.calc_alphar_deriv_nocache(0, 0, HEOS.mole_fractions, tau, delta);
        CoolPropDbl dar_dtau = _HEOS.calc_alphar_deriv_nocache(1, 0, HEOS.mole_fractions, tau, delta);
        CoolPropDbl dar_ddelta = _HEOS.calc_alphar_deriv_nocache(0, 1, HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_ddelta_dtau = _HEOS.calc_alphar_deriv_nocache(1, 1, HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_ddelta2 = _HEOS.calc_alphar_deriv_nocache(0, 2, HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_dtau2 = _HEOS.calc_alphar_deriv_nocache(2, 0, HEOS.mole_fractions, tau, delta);

        CoolPropDbl f1 = delta / tau * (1 + delta * dar_ddelta) - p / (rhoc * R * Tc);
        CoolPropDbl df1_dtau = (1 + delta * dar_ddelta) * (-delta / tau / tau) + delta / tau * (delta * d2ar_ddelta_dtau);
        CoolPropDbl df1_ddelta = (1.0 / tau) * (1 + 2.0 * delta * dar_ddelta + delta * delta * d2ar_ddelta2);
        switch (other) {
            case iHmolar: {
                f2 = (1 + delta * dar_ddelta) + tau * (da0_dtau + dar_dtau) - tau * y / (R * Tc);
                df2_dtau = delta * d2ar_ddelta_dtau + da0_dtau + dar_dtau + tau * (d2a0_dtau2 + d2ar_dtau2) - y / (R * Tc);
                df2_ddelta = (dar_ddelta + delta * d2ar_ddelta2) + tau * (d2a0_ddelta_dtau + d2ar_ddelta_dtau);
                break;
            }
            case iSmolar: {
                f2 = tau * (da0_dtau + dar_dtau) - ar - a0 - y / R;
                df2_dtau = tau * (d2a0_dtau2 + d2ar_dtau2) + (da0_dtau + dar_dtau) - dar_dtau - da0_dtau;
                df2_ddelta = tau * (d2a0_ddelta_dtau + d2ar_ddelta_dtau) - dar_ddelta - da0_ddelta;
                break;
            }
            default:
                throw ValueError("other variable in HSU_P_flash_singlephase_Newton is invalid");
        }

        //First index is the row, second index is the column
        A[0][0] = df1_dtau;
        A[0][1] = df1_ddelta;
        A[1][0] = df2_dtau;
        A[1][1] = df2_ddelta;

        //double det = A[0][0]*A[1][1]-A[1][0]*A[0][1];

        MatInv_2(A, B);
        tau -= omega * (B[0][0] * f1 + B[0][1] * f2);
        delta -= omega * (B[1][0] * f1 + B[1][1] * f2);

        if (std::abs(f1) > std::abs(f2))
            worst_error = std::abs(f1);
        else
            worst_error = std::abs(f2);

        if (!ValidNumber(f1) || !ValidNumber(f2)) {
            throw SolutionError(format("Invalid values for inputs p=%g y=%g for fluid %s in HSU_P_flash_singlephase", p, y, _HEOS.name().c_str()));
        }

        iter += 1;
        if (iter > 100) {
            throw SolutionError(format("HSU_P_flash_singlephase did not converge with inputs p=%g h=%g for fluid %s", p, y, _HEOS.name().c_str()));
        }
    }

    HEOS.update(DmolarT_INPUTS, rhoc * delta, Tc / tau);
}
void FlashRoutines::HSU_P_flash_singlephase_Brent(HelmholtzEOSMixtureBackend& HEOS, parameters other, CoolPropDbl value, CoolPropDbl Tmin,
                                                  CoolPropDbl Tmax, phases phase) {
    if (!ValidNumber(HEOS._p)) {
        throw ValueError("value for p in HSU_P_flash_singlephase_Brent is invalid");
    };
    if (!ValidNumber(value)) {
        throw ValueError("value for other in HSU_P_flash_singlephase_Brent is invalid");
    };
    class solver_resid : public FuncWrapper1DWithTwoDerivs
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl p, value;
        parameters other;
        int iter;
        CoolPropDbl eos0, eos1, rhomolar, rhomolar0, rhomolar1;
        CoolPropDbl Tmin, Tmax;

        solver_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl p, CoolPropDbl value, parameters other, double Tmin, double Tmax)
          : HEOS(HEOS),
            p(p),
            value(value),
            other(other),
            iter(0),
            eos0(-_HUGE),
            eos1(-_HUGE),
            rhomolar(_HUGE),
            rhomolar0(_HUGE),
            rhomolar1(_HUGE),
            Tmin(Tmin),
            Tmax(Tmax) {
            // Specify the state to avoid saturation calls, but only if phase is subcritical
            switch (CoolProp::phases phase = HEOS->phase()) {
                case iphase_liquid:
                case iphase_gas:
                    HEOS->specify_phase(phase);
                default:
                    // Otherwise don't do anything (this is to make compiler happy)
                    {}
            }
        }
        double call(double T) {

            if (iter < 2 || std::abs(rhomolar1 / rhomolar0 - 1) > 0.05) {
                // Run the solver with T,P as inputs; but only if the last change in density was greater than a few percent
                HEOS->update(PT_INPUTS, p, T);
            } else {
                // Run the solver with T,P as inputs; but use the guess value for density from before
                HEOS->update_TP_guessrho(T, p, rhomolar);
            }

            // Get the value of the desired variable
            CoolPropDbl eos = HEOS->keyed_output(other);

            // Store the value of density
            rhomolar = HEOS->rhomolar();

            if (verbosity > 0 && iter == 0) {
                std::cout << format("T: %0.15g rho: %0.15g eos: %0.15g", T, rhomolar, eos);
            }

            // Difference between the two is to be driven to zero
            CoolPropDbl r = eos - value;

            // Store values for later use if there are errors
            if (iter == 0) {
                eos0 = eos;
                rhomolar0 = rhomolar;
            } else if (iter == 1) {
                eos1 = eos;
                rhomolar1 = rhomolar;
            } else {
                eos0 = eos1;
                eos1 = eos;
                rhomolar0 = rhomolar1;
                rhomolar1 = rhomolar;
            }

            iter++;
            return r;
        };
        double deriv(double T) {
            return HEOS->first_partial_deriv(other, iT, iP);
        }
        double second_deriv(double T) {
            return HEOS->second_partial_deriv(other, iT, iP, iT, iP);
        }
        bool input_not_in_range(double x) {
            return (x < Tmin || x > Tmax);
        }
    };
    solver_resid resid(&HEOS, HEOS._p, value, other, Tmin, Tmax);

    class resid_2D : public FuncWrapperND
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl p, value;
        parameters other;
        int iter;
        std::vector<std::vector<double>> J = {{-1.0, -1.0}, {-1.0, -1.0}};

        resid_2D(HelmholtzEOSMixtureBackend* HEOS, double p, CoolPropDbl value, parameters other, double Tmin, double Tmax)
          : HEOS(HEOS), p(p), value(value), other(other), iter(0) {}
        std::vector<double> call(const std::vector<double>& x) override {
            double T = x[0], rhomolar = x[1];
            HEOS->update_DmolarT_direct(rhomolar, T);
            J[0][0] = HEOS->first_partial_deriv(iP, iT, iDmolar) / p;
            J[0][1] = HEOS->first_partial_deriv(iP, iDmolar, iT) / p;
            J[1][0] = HEOS->first_partial_deriv(other, iT, iDmolar);
            J[1][1] = HEOS->first_partial_deriv(other, iDmolar, iT);

            return {(HEOS->p() - p) / p, HEOS->keyed_output(other) - value};
        }
        std::vector<std::vector<double>> Jacobian(const std::vector<double>& /*x*/) override {
            // pre-calculated in call function
            return J;
        }
    };
    resid_2D solver_resid2d(&HEOS, HEOS._p, value, other, Tmin, Tmax);

    // Get residual values at the bounds
    double resid_Tmin = resid.call(Tmin);
    double rhomolar_Tmin = HEOS.rhomolar();

    double resid_Tmax = resid.call(Tmax);
    double rhomolar_Tmax = HEOS.rhomolar();

    // For the derivative-based methods, figure out which point to start from
    bool use_min = std::abs(resid_Tmin) < std::abs(resid_Tmax);
    double Tstart = use_min ? Tmin : Tmax;
    double rhomolarstart = use_min ? rhomolar_Tmin : rhomolar_Tmax;

    try {
        if (get_debug_level() > 0) {
            resid.verbosity = 1;
        }
        if (resid_Tmin * resid_Tmax < 0) {
            // The residual values bound zero, use the TOMS748 method (no derivatives)
            // See: https://www.boost.org/doc/libs/1_58_0/libs/math/doc/html/math_toolkit/internals1/roots2.html#math_toolkit.internals1.roots2.algorithm_toms_748_alefeld_potra
            //
            // It is like a supercharged version of Brent's method, which is practically guaranteed
            // to converge for any continuous function, and take the optimal step among bisection
            // and higher-order methods
            resid.iter = 0;
            boost::math::uintmax_t max_iter = 100;

            auto f = [&resid](const double T) { return resid.call(T); };
            // Want 44 bits to be correct, tolerance is 2^(1-bits) ::
            // >>> 2**(1-44)
            // 1.1368683772161603e-13
            auto [l, r] = toms748_solve(f, Tmin, Tmax, resid_Tmin, resid_Tmax, boost::math::tools::eps_tolerance<double>(44), max_iter);
            if (!is_in_closed_range(Tmin, Tmax, static_cast<CoolPropDbl>(resid.HEOS->T()))) {
                throw ValueError(format("TOMS748 method yielded out of bound T of %g", static_cast<CoolPropDbl>(resid.HEOS->T())));
            }

            // Un-specify the phase of the fluid
            HEOS.unspecify_phase();
            HEOS.recalculate_singlephase_phase();
        } else {
            resid.iter = 0;
            Halley(resid, Tstart, 1e-12, 100);
            if (!is_in_closed_range(Tmin, Tmax, static_cast<CoolPropDbl>(resid.HEOS->T())) || resid.HEOS->phase() != phase) {
                throw ValueError("Halley's method was unable to find a solution in HSU_P_flash_singlephase_Brent");
            }
            // Un-specify the phase of the fluid
            HEOS.unspecify_phase();
            HEOS.recalculate_singlephase_phase();
        }
    } catch (...) {
        try {
            double p_critical_ = HEOS.p_critical();
            if (0.95 * p_critical_ > HEOS._p || HEOS._p > p_critical_) {
                throw;
            }
            if (get_debug_level() > 0) {
                std::cout << resid.errstring << std::endl;
            }
            std::vector<double> x0 = {Tstart, rhomolarstart};
            NDNewtonRaphson_Jacobian(&solver_resid2d, x0, 1e-12, 20, 1.0);
            if (!is_in_closed_range(Tmin, Tmax, static_cast<CoolPropDbl>(solver_resid2d.HEOS->T())) || solver_resid2d.HEOS->phase() != phase) {
                throw ValueError("2D Newton method was unable to find a solution in HSU_P_flash_singlephase_Brent");
            }
            // Un-specify the phase of the fluid
            HEOS.unspecify_phase();
            HEOS.recalculate_singlephase_phase();
        } catch (...) {
            if (get_debug_level() > 0) {
                std::cout << resid.errstring << std::endl;
            }
            // Un-specify the phase of the fluid
            HEOS.unspecify_phase();

            // Determine why you were out of range if you can
            //
            CoolPropDbl eos0 = resid.eos0, eos1 = resid.eos1;
            std::string name = get_parameter_information(other, "short");
            std::string units = get_parameter_information(other, "units");
            if (eos1 > eos0 && value > eos1) {
                throw ValueError(
                  format("HSU_P_flash_singlephase_Brent could not find a solution because %s [%Lg %s] is above the maximum value of %0.12Lg %s",
                         name.c_str(), value, units.c_str(), eos1, units.c_str()));
            }
            if (eos1 > eos0 && value < eos0) {
                throw ValueError(
                  format("HSU_P_flash_singlephase_Brent could not find a solution because %s [%Lg %s] is below the minimum value of %0.12Lg %s",
                         name.c_str(), value, units.c_str(), eos0, units.c_str()));
            }
            throw;
        }
    }
}

// P given and one of H, S, or U
void FlashRoutines::HSU_P_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other) {
    bool saturation_called = false;
    CoolPropDbl value = NAN;

    // Find the phase, while updating all internal variables possible
    switch (other) {
        case iSmolar:
            value = HEOS.smolar();
            break;
        case iHmolar:
            value = HEOS.hmolar();
            break;
        case iUmolar:
            value = HEOS.umolar();
            break;
        default:
            throw ValueError(format("Input for other [%s] is invalid", get_parameter_information(other, "long").c_str()));
    }
    if (HEOS.is_pure_or_pseudopure) {

        // Find the phase, while updating all internal variables possible
        HEOS.p_phase_determination_pure_or_pseudopure(other, value, saturation_called);

        if (HEOS.isHomogeneousPhase()) {
            // Now we use the single-phase solver to find T,rho given P,Y using a
            // bounded 1D solver by adjusting T and using given value of p
            CoolPropDbl Tmin = NAN, Tmax = NAN;
            switch (HEOS._phase) {
                case iphase_gas: {
                    Tmax = 1.5 * HEOS.Tmax();
                    if (HEOS._p < HEOS.p_triple()) {
                        Tmin = std::max(HEOS.Tmin(), HEOS.Ttriple());
                    } else {

                        if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
                            auto superanc_ptr = HEOS.get_superanc();
                            if (superanc_ptr) {
                                auto& superanc = *superanc_ptr;
                                CoolPropDbl pmax_num = superanc.get_pmax();
                                if (HEOS._p > pmax_num) {
                                    throw ValueError(
                                      format("Pressure to PQ_flash [%0.8Lg Pa] may not be above the numerical critical point of %0.15Lg Pa", HEOS._p,
                                             pmax_num));
                                }
                                Tmin = superanc.get_T_from_p(HEOS._p);
                                break;
                            }
                        }
                        if (saturation_called) {
                            Tmin = HEOS.SatV->T();
                        } else {
                            Tmin = HEOS._TVanc.pt() - 0.01;
                        }
                    }
                    break;
                }
                case iphase_liquid: {

                    // Sometimes the minimum pressure for the melting line is a bit above the triple point pressure
                    if (HEOS.has_melting_line() && HEOS._p > HEOS.calc_melting_line(iP_min, -1, -1)) {
                        Tmin = HEOS.calc_melting_line(iT, iP, HEOS._p) - 1e-3;
                    } else {
                        Tmin = HEOS.Tmin() - 1e-3;
                    }

                    if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
                        auto superanc_ptr = HEOS.get_superanc();
                        if (superanc_ptr) {
                            auto& superanc = *superanc_ptr;
                            CoolPropDbl pmax_num = superanc.get_pmax();
                            if (HEOS._p > pmax_num) {
                                throw ValueError(format(
                                  "Pressure to PQ_flash [%0.8Lg Pa] may not be above the numerical critical point of %0.15Lg Pa", HEOS._p, pmax_num));
                            }
                            Tmax = superanc.get_T_from_p(HEOS._p);
                            break;
                        }
                    }

                    if (saturation_called) {
                        Tmax = HEOS.SatL->T();
                    } else {
                        Tmax = HEOS._TLanc.pt() + 0.01;
                    }

                    break;
                }
                case iphase_supercritical_liquid:
                case iphase_supercritical_gas:
                case iphase_supercritical: {
                    Tmax = 1.5 * HEOS.Tmax();
                    // Sometimes the minimum pressure for the melting line is a bit above the triple point pressure
                    if (HEOS.has_melting_line() && HEOS._p > HEOS.calc_melting_line(iP_min, -1, -1)) {
                        Tmin = HEOS.calc_melting_line(iT, iP, HEOS._p) - 1e-3;
                    } else {
                        Tmin = HEOS.Tmin() - 1e-3;
                    }
                    break;
                }
                default: {
                    throw ValueError(format("Not a valid homogeneous state"));
                }
            }
            try {
                HSU_P_flash_singlephase_Brent(HEOS, other, value, Tmin, Tmax, HEOS._phase);
            } catch (std::exception& e) {
                throw ValueError(format("unable to solve 1phase PY flash with Tmin=%Lg, Tmax=%Lg due to error: %s", Tmin, Tmax, e.what()));
            }
            HEOS._Q = -1;
            // Update the state for conditions where the state was guessed
            HEOS.recalculate_singlephase_phase();
        }
    } else {
        if (HEOS.PhaseEnvelope.built) {
            // Determine whether you are inside or outside
            SimpleState closest_state;
            std::size_t iclosest = 0;
            bool twophase = PhaseEnvelopeRoutines::is_inside(HEOS.PhaseEnvelope, iP, HEOS._p, other, value, iclosest, closest_state);

            if (!twophase) {
                PY_singlephase_flash_resid resid(HEOS, HEOS._p, other, value);
                // If that fails, try a bounded solver
                Brent(resid, closest_state.T + 10, 1000, DBL_EPSILON, 1e-10, 100);
                HEOS.unspecify_phase();
            } else {
                throw ValueError("two-phase solution for Y");
            }

        } else {
            throw ValueError("phase envelope must be built to carry out HSU_P_flash for mixture");
        }
    }
}
void FlashRoutines::solver_for_rho_given_T_oneof_HSU(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl T, CoolPropDbl value, parameters other) {
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1DWithTwoDerivs
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl T, value;
        parameters other;

        solver_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl T, CoolPropDbl value, parameters other)
          : HEOS(HEOS), T(T), value(value), other(other) {}
        double call(double rhomolar) {
            HEOS->update_DmolarT_direct(rhomolar, T);
            double eos = HEOS->keyed_output(other);
            return eos - value;
        };
        double deriv(double rhomolar) {
            return HEOS->first_partial_deriv(other, iDmolar, iT);
        }
        double second_deriv(double rhomolar) {
            return HEOS->second_partial_deriv(other, iDmolar, iT, iDmolar, iT);
        }
    };
    solver_resid resid(&HEOS, T, value, other);

    double T_critical_ = (HEOS.is_pure_or_pseudopure) ? HEOS.T_critical() : HEOS._crit.T;

    // Supercritical temperature
    if (HEOS._T > T_critical_) {
        CoolPropDbl yc = NAN, ymin = NAN, y = NAN;
        CoolPropDbl rhoc = (HEOS.is_pure_or_pseudopure) ? HEOS.rhomolar_critical() : HEOS._crit.rhomolar;
        CoolPropDbl rhomin = 1e-10;

        // Determine limits for the other variable
        switch (other) {
            case iSmolar: {
                yc = HEOS.calc_smolar_nocache(HEOS._T, rhoc);
                ymin = HEOS.calc_smolar_nocache(HEOS._T, rhomin);
                y = HEOS._smolar;
                break;
            }
            case iHmolar: {
                yc = HEOS.calc_hmolar_nocache(HEOS._T, rhoc);
                ymin = HEOS.calc_hmolar_nocache(HEOS._T, rhomin);
                y = HEOS._hmolar;
                break;
            }
            case iUmolar: {
                yc = HEOS.calc_umolar_nocache(HEOS._T, rhoc);
                ymin = HEOS.calc_umolar_nocache(HEOS._T, rhomin);
                y = HEOS._umolar;
                break;
            }
            default:
                throw ValueError();
        }
        if (is_in_closed_range(yc, ymin, y)) {
            Brent(resid, rhoc, rhomin, LDBL_EPSILON, 1e-9, 100);
        } else if (y < yc) {
            // Increase rhomelt until it bounds the solution
            int step_count = 0;
            while (!is_in_closed_range(ymin, yc, y)) {
                rhoc *= 1.1;  // Increase density by a few percent
                switch (other) {
                    case iSmolar:
                        yc = HEOS.calc_smolar_nocache(HEOS._T, rhoc);
                        break;
                    case iHmolar:
                        yc = HEOS.calc_hmolar_nocache(HEOS._T, rhoc);
                        break;
                    case iUmolar:
                        yc = HEOS.calc_umolar_nocache(HEOS._T, rhoc);
                        break;
                    default:
                        throw ValueError(format("Input is invalid"));
                }
                if (step_count > 30) {
                    throw ValueError(format("Even by increasing rhoc, not able to bound input; input %Lg is not in range %Lg,%Lg", y, yc, ymin));
                }
                step_count++;
            }
            Brent(resid, rhomin, rhoc, LDBL_EPSILON, 1e-9, 100);
        } else {
            throw ValueError(format("input %Lg is not in range %Lg,%Lg,%Lg", y, yc, ymin));
        }
        // Update the state (T > Tc)
        if (HEOS._p < HEOS.p_critical()) {
            HEOS._phase = iphase_supercritical_gas;
        } else {
            HEOS._phase = iphase_supercritical;
        }
    }
    // Subcritical temperature liquid
    else if ((HEOS._phase == iphase_liquid) || (HEOS._phase == iphase_supercritical_liquid)) {
        CoolPropDbl ymelt = NAN, yL = NAN, y = NAN;
        CoolPropDbl rhomelt = HEOS.components[0].triple_liquid.rhomolar;
        CoolPropDbl rhoL = static_cast<double>(HEOS._rhoLanc);

        switch (other) {
            case iSmolar: {
                ymelt = HEOS.calc_smolar_nocache(HEOS._T, rhomelt);
                yL = HEOS.calc_smolar_nocache(HEOS._T, rhoL);
                y = HEOS._smolar;
                break;
            }
            case iHmolar: {
                ymelt = HEOS.calc_hmolar_nocache(HEOS._T, rhomelt);
                yL = HEOS.calc_hmolar_nocache(HEOS._T, rhoL);
                y = HEOS._hmolar;
                break;
            }
            case iUmolar: {
                ymelt = HEOS.calc_umolar_nocache(HEOS._T, rhomelt);
                yL = HEOS.calc_umolar_nocache(HEOS._T, rhoL);
                y = HEOS._umolar;
                break;
            }
            default:
                throw ValueError();
        }

        CoolPropDbl rhomolar_guess = (rhomelt - rhoL) / (ymelt - yL) * (y - yL) + rhoL;

        try {
            Halley(resid, rhomolar_guess, 1e-8, 100);
        } catch (...) {
            Secant(resid, rhomolar_guess, 0.0001 * rhomolar_guess, 1e-12, 100);
        }
    }
    // Subcritical temperature gas
    else if (HEOS._phase == iphase_gas) {
        CoolPropDbl rhomin = 1e-14;
        CoolPropDbl rhoV = static_cast<double>(HEOS._rhoVanc);

        try {
            Halley(resid, 0.5 * (rhomin + rhoV), 1e-8, 100);
        } catch (...) {
            try {
                Brent(resid, rhomin, rhoV, LDBL_EPSILON, 1e-12, 100);
            } catch (...) {
                throw ValueError();
            }
        }
    } else {
        throw ValueError(format("phase to solver_for_rho_given_T_oneof_HSU is invalid"));
    }
};

void FlashRoutines::DHSU_T_flash(HelmholtzEOSMixtureBackend& HEOS, parameters other) {
    if (HEOS.imposed_phase_index != iphase_not_imposed) {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
        double T_critical_ = (HEOS.is_pure_or_pseudopure) ? HEOS.T_critical() : HEOS._crit.T;
        // The remaining code in this branch was added to set some needed parameters if phase is imposed,
        // since HEOS.T_phase_determination_pure_or_pseudopure() is not being called.
        if (HEOS._T < T_critical_)  //
        {
            // TODO: is it a bug that this branch can be accessed for mixtures?
            HEOS._rhoVanc = HEOS.components[0].ancillaries.rhoV.evaluate(HEOS._T);
            HEOS._rhoLanc = HEOS.components[0].ancillaries.rhoL.evaluate(HEOS._T);
            if (HEOS._phase == iphase_liquid) {
                HEOS._Q = -1000;
            } else if (HEOS._phase == iphase_gas) {
                HEOS._Q = 1000;
            } else if (HEOS._phase == iphase_twophase) {
                // Actually have to use saturation information sadly
                // For the given temperature, find the saturation state
                // Run the saturation routines to determine the saturation densities and pressures
                HelmholtzEOSMixtureBackend HEOS1(HEOS.components);
                SaturationSolvers::saturation_T_pure_options options;
                SaturationSolvers::saturation_T_pure(HEOS1, HEOS._T, options);

                if (other != iDmolar) {
                    // Update the states
                    if (HEOS.SatL) HEOS.SatL->update(DmolarT_INPUTS, HEOS._rhoLanc, HEOS._T);
                    if (HEOS.SatV) HEOS.SatV->update(DmolarT_INPUTS, HEOS._rhoVanc, HEOS._T);
                    // Update the two-Phase variables
                    HEOS._rhoLmolar = HEOS.SatL->rhomolar();
                    HEOS._rhoVmolar = HEOS.SatV->rhomolar();
                }

                CoolPropDbl Q = NAN;

                switch (other) {
                    case iDmolar:
                        Q = (1 / HEOS.rhomolar() - 1 / HEOS1.SatL->rhomolar()) / (1 / HEOS1.SatV->rhomolar() - 1 / HEOS1.SatL->rhomolar());
                        break;
                    case iSmolar:
                        Q = (HEOS.smolar() - HEOS1.SatL->smolar()) / (HEOS1.SatV->smolar() - HEOS1.SatL->smolar());
                        break;
                    case iHmolar:
                        Q = (HEOS.hmolar() - HEOS1.SatL->hmolar()) / (HEOS1.SatV->hmolar() - HEOS1.SatL->hmolar());
                        break;
                    case iUmolar:
                        Q = (HEOS.umolar() - HEOS1.SatL->umolar()) / (HEOS1.SatV->umolar() - HEOS1.SatL->umolar());
                        break;
                    default:
                        throw ValueError(format("bad input for other"));
                }
                if (Q < 0) {
                    HEOS._Q = -1;
                } else if (Q > 1) {
                    HEOS._Q = 1;
                } else {
                    HEOS._Q = Q;
                    // Load the outputs
                    HEOS._p = HEOS._Q * HEOS1.SatV->p() + (1 - HEOS._Q) * HEOS1.SatL->p();
                    HEOS._rhomolar = 1 / (HEOS._Q / HEOS.SatV->rhomolar() + (1 - HEOS._Q) / HEOS.SatL->rhomolar());
                }
            } else if (HEOS._phase == iphase_supercritical_liquid) {
                HEOS._Q = -1000;
            } else
                throw ValueError(format("Temperature specified is not the imposed phase region."));
        } else if (HEOS._T > T_critical_ && HEOS._T > HEOS.components[0].EOS().Ttriple) {
            HEOS._Q = 1e9;
        }
    } else {
        if (HEOS.is_pure_or_pseudopure) {
            // Find the phase, while updating all internal variables possible
            switch (other) {
                case iDmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iDmolar, HEOS._rhomolar);
                    break;
                case iSmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iSmolar, HEOS._smolar);
                    break;
                case iHmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iHmolar, HEOS._hmolar);
                    break;
                case iUmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iUmolar, HEOS._umolar);
                    break;
                default:
                    throw ValueError(format("Input is invalid"));
            }
        } else {
            HEOS._phase = iphase_gas;
            throw NotImplementedError("DHSU_T_flash does not support mixtures (yet)");
        }
    }

    //if (HEOS.isHomogeneousPhase() && !ValidNumber(HEOS._p)) // original, pre 1352
    // only the solver requires single phase
    if (((other == iDmolar) || HEOS.isHomogeneousPhase()) && !ValidNumber(HEOS._p))  // post 1352
    {
        switch (other) {
            case iDmolar:
                break;
            case iHmolar:
                solver_for_rho_given_T_oneof_HSU(HEOS, HEOS._T, HEOS._hmolar, iHmolar);
                break;
            case iSmolar:
                solver_for_rho_given_T_oneof_HSU(HEOS, HEOS._T, HEOS._smolar, iSmolar);
                break;
            case iUmolar:
                solver_for_rho_given_T_oneof_HSU(HEOS, HEOS._T, HEOS._umolar, iUmolar);
                break;
            default:
                break;
        }
        HEOS.calc_pressure();
        HEOS._Q = -1;
    }
    if (HEOS.is_pure_or_pseudopure && HEOS._phase != iphase_twophase) {
        // Update the state for conditions where the state was guessed
        HEOS.recalculate_singlephase_phase();
    }
}
void FlashRoutines::HS_flash_twophase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec,
                                      HS_flash_twophaseOptions& options) {
    class Residual : public FuncWrapper1D
    {

       public:
        HelmholtzEOSMixtureBackend& HEOS;
        CoolPropDbl hmolar, smolar, Qs;
        Residual(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec)
          : HEOS(HEOS), hmolar(hmolar_spec), smolar(smolar_spec), Qs(_HUGE) {};
        double call(double T) {
            HEOS.update(QT_INPUTS, 0, T);
            HelmholtzEOSMixtureBackend &SatL = HEOS.get_SatL(), &SatV = HEOS.get_SatV();
            // Quality from entropy
            Qs = (smolar - SatL.smolar()) / (SatV.smolar() - SatL.smolar());
            // Quality from enthalpy
            CoolPropDbl Qh = (hmolar - SatL.hmolar()) / (SatV.hmolar() - SatL.hmolar());
            // Residual is the difference between the two
            return Qh - Qs;
        }
    } resid(HEOS, hmolar_spec, smolar_spec);

    // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
    CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() - 1e-13;

    // Check what the minimum limits for the equation of state are
    CoolPropDbl Tmin_satL = NAN, Tmin_satV = NAN, Tmin_sat = NAN;
    HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
    Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;

    Brent(resid, Tmin_sat, Tmax_sat - 0.01, DBL_EPSILON, 1e-12, 20);
    // Run once more with the final vapor quality
    HEOS.update(QT_INPUTS, resid.Qs, HEOS.T());
}
void FlashRoutines::HS_flash_singlephase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec,
                                         HS_flash_singlephaseOptions& options) {
    int iter = 0;
    double resid = 9e30, resid_old = 9e30;
    CoolProp::SimpleState reducing = HEOS.get_state("reducing");
    do {
        // Independent variables are T0 and rhomolar0, residuals are matching h and s
        Eigen::Vector2d r;
        Eigen::Matrix2d J;
        r(0) = HEOS.hmolar() - hmolar_spec;
        r(1) = HEOS.smolar() - smolar_spec;
        J(0, 0) = HEOS.first_partial_deriv(iHmolar, iTau, iDelta);
        J(0, 1) = HEOS.first_partial_deriv(iHmolar, iDelta, iTau);
        J(1, 0) = HEOS.first_partial_deriv(iSmolar, iTau, iDelta);
        J(1, 1) = HEOS.first_partial_deriv(iSmolar, iDelta, iTau);
        // Step in v obtained from Jv = -r
        Eigen::Vector2d v = J.colPivHouseholderQr().solve(-r);
        bool good_solution = false;
        double tau0 = HEOS.tau(), delta0 = HEOS.delta();
        // Calculate the old residual after the last step
        resid_old = sqrt(POW2(HEOS.hmolar() - hmolar_spec) + POW2(HEOS.smolar() - smolar_spec));
        for (double frac = 1.0; frac > 0.001; frac /= 2) {
            try {
                // Calculate new values
                double tau_new = tau0 + options.omega * frac * v(0);
                double delta_new = delta0 + options.omega * frac * v(1);
                double T_new = reducing.T / tau_new;
                double rhomolar_new = delta_new * reducing.rhomolar;
                // Update state with step
                HEOS.update(DmolarT_INPUTS, rhomolar_new, T_new);
                resid = sqrt(POW2(HEOS.hmolar() - hmolar_spec) + POW2(HEOS.smolar() - smolar_spec));
                if (resid > resid_old) {
                    throw ValueError(format("residual not decreasing; frac: %g, resid: %g, resid_old: %g", frac, resid, resid_old));
                }
                good_solution = true;
                break;
            } catch (...) {
                HEOS.clear();
                continue;
            }
        }
        if (!good_solution) {
            throw ValueError(format("Not able to get a solution"));
        }
        iter++;
        if (iter > 50) {
            throw ValueError(format("HS_flash_singlephase took too many iterations; residual is %g; prior was %g", resid, resid_old));
        }
    } while (std::abs(resid) > 1e-9);
}
void FlashRoutines::HS_flash_generate_TP_singlephase_guess(HelmholtzEOSMixtureBackend& HEOS, double& T, double& p) {
    // Randomly obtain a starting value that is single-phase
    double logp = ((double)rand() / (double)RAND_MAX) * (log(HEOS.pmax()) - log(HEOS.p_triple())) + log(HEOS.p_triple());
    T = ((double)rand() / (double)RAND_MAX) * (HEOS.Tmax() - HEOS.Ttriple()) + HEOS.Ttriple();
    p = exp(logp);
}
void FlashRoutines::HS_flash(HelmholtzEOSMixtureBackend& HEOS) {
    // Use TS flash and iterate on T (known to be between Tmin and Tmax)
    // in order to find H
    double hmolar = HEOS.hmolar(), smolar = HEOS.smolar();
    class Residual : public FuncWrapper1D
    {
       public:
        HelmholtzEOSMixtureBackend& HEOS;
        CoolPropDbl hmolar, smolar;
        Residual(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec)
          : HEOS(HEOS), hmolar(hmolar_spec), smolar(smolar_spec) {};
        double call(double T) {
            HEOS.update(SmolarT_INPUTS, smolar, T);
            double r = HEOS.hmolar() - hmolar;
            return r;
        }
    } resid(HEOS, hmolar, smolar);
    std::string errstr;
    // Find minimum temperature
    bool good_Tmin = false;
    double Tmin = HEOS.Ttriple();
    double rmin = NAN;
    do {
        try {
            rmin = resid.call(Tmin);
            good_Tmin = true;
        } catch (...) {
            Tmin += 0.5;
        }
        if (Tmin > HEOS.Tmax()) {
            throw ValueError("Cannot find good Tmin");
        }
    } while (!good_Tmin);

    // Find maximum temperature
    bool good_Tmax = false;
    double Tmax = HEOS.Tmax() * 1.01;  // Just a little above, so if we use Tmax as input, it should still work
    double rmax = NAN;
    do {
        try {
            rmax = resid.call(Tmax);
            good_Tmax = true;
        } catch (...) {
            Tmax -= 0.1;
        }
        if (Tmax < Tmin) {
            throw ValueError("Cannot find good Tmax");
        }
    } while (!good_Tmax);
    if (rmin * rmax > 0 && std::abs(rmax) < std::abs(rmin)) {
        throw CoolProp::ValueError(format("HS inputs correspond to temperature above maximum temperature of EOS [%g K]", HEOS.Tmax()));
    }
    Brent(resid, Tmin, Tmax, DBL_EPSILON, 1e-10, 100);
}

#if defined(ENABLE_CATCH)

TEST_CASE("PD with T very large should yield error", "[PDflash]") {
    shared_ptr<HelmholtzEOSBackend> HEOS(new HelmholtzEOSBackend("R134a"));
    double Tc = HEOS->T_critical();
    HEOS->update(DmassT_INPUTS, 1.1, 1.5 * Tc);
    CHECK_THROWS(HEOS->update(DmassP_INPUTS, 2, 5 * HEOS->p()));
}

TEST_CASE("Stability testing", "[stability]") {
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS(new HelmholtzEOSMixtureBackend(strsplit("n-Propane&n-Butane&n-Pentane&n-Hexane", '&')));
    std::vector<double> z(4);
    z[0] = 0.1;
    z[1] = 0.2;
    z[2] = 0.3;
    z[3] = 0.4;
    HEOS->set_mole_fractions(z);

    HEOS->update(PQ_INPUTS, 101325, 0);
    double TL = HEOS->T();

    HEOS->update(PQ_INPUTS, 101325, 1);
    double TV = HEOS->T();

    SECTION("Liquid (feed is stable)") {
        StabilityRoutines::StabilityEvaluationClass stability_tester(*HEOS);
        for (double T = TL - 1; T >= 100; T -= 1) {
            stability_tester.set_TP(T, 101325);
            CAPTURE(T);
            CHECK_NOTHROW(stability_tester.is_stable());
        }
    }
    SECTION("Vapor (feed is stable)") {
        StabilityRoutines::StabilityEvaluationClass stability_tester(*HEOS);
        for (double T = TV + 1; T <= 500; T += 1) {
            stability_tester.set_TP(T, 101325);
            CAPTURE(T);
            CHECK_NOTHROW(stability_tester.is_stable());
        }
    }
    SECTION("Two-phase (feed is unstable)") {
        StabilityRoutines::StabilityEvaluationClass stability_tester(*HEOS);
        stability_tester.set_TP((TV + TL) / 2.0, 101325);
        CHECK(stability_tester.is_stable() == false);
    }
}

TEST_CASE("Test critical points for methane + H2S", "[critical_points]") {
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS(new HelmholtzEOSMixtureBackend(strsplit("Methane&H2S", '&')));

    double zz[] = {0.998, 0.97, 0.9475, 0.94, 0.93, 0.86, 0.85, 0.84, 0.75, 0.53, 0.52, 0.51, 0.49, 0.36, 0.24, 0.229, 0.09};
    int Npts[] = {2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 1, 1, 1, 1};  // Number of critical points that should be found
    int imax = sizeof(zz) / sizeof(double);

    for (int i = 0; i < imax; ++i) {
        double z0 = zz[i];
        std::vector<double> z(2);
        z[0] = z0;
        z[1] = 1 - z0;
        HEOS->set_mole_fractions(z);
        CAPTURE(z0);
        std::vector<CriticalState> pts = HEOS->all_critical_points();
        CHECK(pts.size() == Npts[i]);
    }
}

TEST_CASE("Test critical points for nitrogen + ethane with HEOS", "[critical_points]") {
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS(new HelmholtzEOSMixtureBackend(strsplit("Nitrogen&Ethane", '&')));
    std::vector<double> zz = linspace(0.001, 0.999, 21);
    int failure_count = 0;
    for (int i = 0; i < static_cast<std::size_t>(zz.size()); ++i) {
        double z0 = zz[i];
        std::vector<double> z(2);
        z[0] = z0;
        z[1] = 1 - z0;
        HEOS->set_mole_fractions(z);
        CAPTURE(z0);
        std::vector<CriticalState> pts;
        try {
            pts = HEOS->all_critical_points();
        } catch (std::exception& e) {
            CAPTURE(e.what());
            failure_count++;
        }
    }
    // Only an error if more than half fail;
    CHECK(failure_count < 10);
}

TEST_CASE("Test critical points for nitrogen + ethane with PR", "[critical_points]") {
    shared_ptr<PengRobinsonBackend> HEOS(new PengRobinsonBackend(strsplit("Nitrogen&Ethane", '&')));
    HEOS->set_binary_interaction_double(0, 1, "kij", 0.0407);  // Ramırez-Jimenez et al.
    std::vector<double> zz = linspace(0.001, 0.999, 21);
    for (int i = 0; i < static_cast<std::size_t>(zz.size()); ++i) {
        double z0 = zz[i];
        std::vector<double> z(2);
        z[0] = z0;
        z[1] = 1 - z0;
        HEOS->set_mole_fractions(z);
        CAPTURE(z0);
        std::vector<CriticalState> pts;
        CHECK_NOTHROW(pts = HEOS->all_critical_points());
    }
}

#endif

} /* namespace CoolProp */
