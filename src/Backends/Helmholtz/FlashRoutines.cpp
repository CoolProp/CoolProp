#include "VLERoutines.h"
#include "FlashRoutines.h"

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "CoolProp/CoolProp.h"
#include "HelmholtzEOSMixtureBackend.h"
#include "HelmholtzEOSBackend.h"
#include "PhaseEnvelopeRoutines.h"
#include "CoolProp/Configuration.h"
#include "CoolProp/fluids/MeltingCaloric.h"

#if defined(ENABLE_CATCH)
#    include <catch2/catch_all.hpp>
#    include "Backends/Cubics/CubicBackend.h"
#endif

#include "boost/math/tools/toms748_solve.hpp"

namespace CoolProp {

void FlashRoutines::PT_flash_mixtures(HelmholtzEOSMixtureBackend& HEOS) {
    if (HEOS.PhaseEnvelope.built && HEOS.PhaseEnvelope.closed) {
        // Use the closed phase envelope to classify the (T, p) point
        try {
            SimpleState closest_state;
            std::size_t iclosest = 0;
            const PhaseEnvelopeData& env = HEOS.PhaseEnvelope;
            bool twophase = PhaseEnvelopeRoutines::is_inside(env, iP, HEOS._p, iT, HEOS._T, iclosest, closest_state);

            if (!twophase) {
                // Single-phase: determine gas vs liquid from temperature relative to closest envelope point
                // Save T and p before solver calls — solver_rho_Tp may corrupt
                // HEOS._T/_p on failure (Householder sets them to -inf).
                const CoolPropDbl T_saved = HEOS._T;
                const CoolPropDbl p_saved = HEOS._p;
                phases phase_guess = (T_saved > closest_state.T) ? iphase_gas : iphase_liquid;
                HEOS.specify_phase(phase_guess);
                double rho;
                try {
                    rho = HEOS.solver_rho_Tp(T_saved, p_saved);
                } catch (...) {
                    // If guessed phase fails, try the other branch
                    phases alt = (phase_guess == iphase_gas) ? iphase_liquid : iphase_gas;
                    HEOS.specify_phase(alt);
                    rho = HEOS.solver_rho_Tp(T_saved, p_saved);
                }
                HEOS.unspecify_phase();
                HEOS.update_DmolarT_direct(rho, T_saved);
                HEOS._Q = -1;
                HEOS._phase = (rho < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
            } else {
                // Two-phase: seed PT flash from the closest envelope point
                std::size_t N = env.K.size();
                CoolProp::SaturationSolvers::PTflash_twophase_options o;
                o.x.resize(N);
                o.y.resize(N);
                for (std::size_t k = 0; k < N; ++k) {
                    o.x[k] = env.x[k][iclosest];
                    o.y[k] = env.y[k][iclosest];
                }
                o.rhomolar_liq = env.rhomolar_liq[iclosest];
                o.rhomolar_vap = env.rhomolar_vap[iclosest];
                o.z = HEOS.get_mole_fractions();
                o.T = HEOS._T;
                o.p = HEOS._p;
                o.omega = 1.0;
                CoolProp::SaturationSolvers::PTflash_twophase solver(HEOS, o);
                solver.solve();
                // Beta-fallback: collapse trivial splits to single-phase
                if (o.beta < 1e-10) {
                    HEOS.update_DmolarT_direct(o.rhomolar_liq, HEOS._T);
                    HEOS._phase = (o.rhomolar_liq < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
                    HEOS._Q = -1;
                } else if (o.beta > 1.0 - 1e-10) {
                    HEOS.update_DmolarT_direct(o.rhomolar_vap, HEOS._T);
                    HEOS._phase = (o.rhomolar_vap < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
                    HEOS._Q = -1;
                } else {
                    HEOS._phase = iphase_twophase;
                    HEOS._Q = o.beta;
                    HEOS._rhomolar = 1.0 / (o.beta / o.rhomolar_vap + (1.0 - o.beta) / o.rhomolar_liq);
                }
            }
            return;
        } catch (...) {
            HEOS.unspecify_phase();
            // Envelope-guided flash failed; fall through to blind flash below
        }
    }

    if (HEOS.imposed_phase_index == iphase_not_imposed) {
        // Blind flash call
        // Instantiates stability tester (dispatches to Michelsen or Gernert based on config)
        StabilityRoutines::StabilityEvaluationClass stability_tester(HEOS);
        bool do_twophase = !stability_tester.is_stable();
        CoolProp::SaturationSolvers::PTflash_twophase_options o;
        o.z = HEOS.get_mole_fractions();
        bool wilson_seeded = false;

        // CoolProp-zgpy: the Michelsen stability trials can converge to the trivial
        // (feed) solution near a phase boundary -- notably for cubic mixtures at high
        // vapor fraction -- and report a false "stable" verdict, so a genuinely
        // two-phase state gets mislabeled single-phase liquid.  Cross-check a "stable"
        // verdict against a cheap ideal (Wilson K-factor) estimate: if it places the
        // feed strictly between its bubble and dew points, attempt a split seeded from
        // that estimate.  A genuinely single-phase point either fails the Wilson test
        // here or collapses back to single phase via the beta~0/1 fallback below.
        // (Lever from jakobreichert's PR #2720.)
        if (!do_twophase) {
            try {
                if (CoolProp::SaturationSolvers::guess_split_from_wilson(HEOS, o.x, o.y, o.rhomolar_liq, o.rhomolar_vap, o.z, HEOS.T(), HEOS.p(),
                                                                         10)) {
                    do_twophase = true;
                    wilson_seeded = true;
                }
            } catch (const CoolProp::CoolPropBaseError&) {
                // Speculative seed failed (e.g. no liquid density root at this state);
                // trust the stability verdict and stay on the single-phase path.
                do_twophase = false;
                wilson_seeded = false;
            }
        }

        if (do_twophase) {
            // There is a phase split and liquid and vapor phases are formed
            if (!wilson_seeded) {
                stability_tester.get_liq(o.x, o.rhomolar_liq);
                stability_tester.get_vap(o.y, o.rhomolar_vap);
            }
            o.T = HEOS.T();
            o.p = HEOS.p();
            o.omega = 1.0;
            CoolProp::SaturationSolvers::PTflash_twophase solver(HEOS, o);
            if (wilson_seeded) {
                // The Wilson split is speculative (the stability test said "stable"): a
                // failure to converge, or a trivial (x == y) result, must fall back to the
                // single-phase path rather than abort the flash or publish a bogus split.
                try {
                    solver.solve();
                } catch (const CoolProp::CoolPropBaseError&) {
                    do_twophase = false;
                }
                if (do_twophase) {
                    CoolPropDbl spread = 0;
                    for (std::size_t i = 0; i < o.z.size(); ++i)
                        spread = std::max(spread, std::abs(o.x[i] - o.y[i]));
                    if (spread < 1e-6) do_twophase = false;  // trivial split -> single phase
                }
            } else {
                solver.solve();
            }
        }

        if (do_twophase) {
            // Fallback block: Catches mathematically trivial splits (beta ~ 0 or 1) caused by boundary
            // floating-point noise in the Michelsen TPD solver, or false positives if using the legacy solver.
            if (o.beta < 1e-10) {
                HEOS.update_DmolarT_direct(o.rhomolar_liq, HEOS.T());
                HEOS._phase = (o.rhomolar_liq < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
                HEOS._Q = -1;
            } else if (o.beta > 1.0 - 1e-10) {
                HEOS.update_DmolarT_direct(o.rhomolar_vap, HEOS.T());
                HEOS._phase = (o.rhomolar_vap < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
                HEOS._Q = -1;
            } else {
                HEOS._phase = iphase_twophase;
                HEOS._Q = o.beta;
                HEOS._rhomolar = 1.0 / (o.beta / o.rhomolar_vap + (1.0 - o.beta) / o.rhomolar_liq);
            }
        } else {
            // It's single-phase -- find the density.
            // Save T and p before any solver calls — solver_rho_Tp may corrupt
            // HEOS._T/_p on failure (Householder sets them to -inf).
            const CoolPropDbl T_saved = HEOS.T();
            const CoolPropDbl p_saved = HEOS.p();

            // Solve SRK cubic for both gas and liquid roots to decide which
            // HEOS branch(es) to solve.  When both roots are valid, solve
            // both HEOS densities and pick the phase with lower Gibbs energy.
            // This mirrors the pattern in solver_rho_Tp_global().
            double rho;
            double rho_srk_gas = HEOS.solver_rho_Tp_SRK(T_saved, p_saved, iphase_gas);
            double rho_srk_liq = HEOS.solver_rho_Tp_SRK(T_saved, p_saved, iphase_liquid);
            bool gas_ok = rho_srk_gas > 0 && ValidNumber(rho_srk_gas);
            bool liq_ok = rho_srk_liq > 0 && ValidNumber(rho_srk_liq);

            if (gas_ok && liq_ok) {
                // Both SRK roots valid — solve both HEOS roots, pick lower Gibbs
                double rho_gas = -1, rho_liq = -1;
                HEOS.specify_phase(iphase_gas);
                try {
                    rho_gas = HEOS.solver_rho_Tp(T_saved, p_saved);
                } catch (...) {
                }
                HEOS.specify_phase(iphase_liquid);
                try {
                    rho_liq = HEOS.solver_rho_Tp(T_saved, p_saved);
                } catch (...) {
                }
                HEOS.unspecify_phase();

                if (rho_gas > 0 && rho_liq > 0) {
                    double G_gas = HEOS.calc_gibbsmolar_nocache(T_saved, rho_gas);
                    double G_liq = HEOS.calc_gibbsmolar_nocache(T_saved, rho_liq);
                    rho = (G_liq <= G_gas) ? rho_liq : rho_gas;
                } else if (rho_gas > 0 || rho_liq > 0) {
                    rho = (rho_gas > 0) ? rho_gas : rho_liq;
                } else {
                    throw ValueError("Unable to obtain either HEOS density root in PT_flash_mixtures");
                }
            } else {
                // Only one SRK root valid (or neither) — solve that branch,
                // fall back to the other if HEOS throws.
                phases primary = gas_ok ? iphase_gas : (liq_ok ? iphase_liquid : iphase_gas);
                phases fallback = (primary == iphase_gas) ? iphase_liquid : iphase_gas;
                HEOS.specify_phase(primary);
                try {
                    rho = HEOS.solver_rho_Tp(T_saved, p_saved);
                } catch (...) {
                    HEOS.specify_phase(fallback);
                    try {
                        rho = HEOS.solver_rho_Tp(T_saved, p_saved);
                    } catch (...) {
                        HEOS.unspecify_phase();
                        throw;
                    }
                }
                HEOS.unspecify_phase();
            }
            HEOS.update_DmolarT_direct(rho, T_saved);
            HEOS._Q = -1;
            HEOS._phase = (rho < HEOS.rhomolar_reducing()) ? iphase_gas : iphase_liquid;
        }
    } else {
        // It's single-phase, and phase is imposed
        const CoolPropDbl T_saved = HEOS.T();
        const CoolPropDbl p_saved = HEOS.p();
        double rho = HEOS.solver_rho_Tp(T_saved, p_saved);
        HEOS.update_DmolarT_direct(rho, T_saved);
        HEOS._Q = -1;
        HEOS._phase = HEOS.imposed_phase_index;
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
    double call(double T) override {
        HEOS->update_DmolarT_direct(rhomolar, T);
        CoolPropDbl peos = HEOS->p();
        CoolPropDbl r = (peos - p) / p;
        return r;
    };
    double deriv(double T) override {
        // dp/dT|rho / pspecified
        return HEOS->first_partial_deriv(iP, iT, iDmolar) / p;
    };
    double second_deriv(double T) override {
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
    double call(double T) override {
        HEOS.update(QT_INPUTS, 0, T);  // Doesn't matter whether liquid or vapor, we are just doing a full VLE call for given T
        double rhoL = HEOS.saturated_liquid_keyed_output(iDmolar);
        double rhoV = HEOS.saturated_vapor_keyed_output(iDmolar);
        /// Error between calculated and target vapor quality based on densities
        return (1 / rhomolar - 1 / rhoL) / (1 / rhoV - 1 / rhoL) - Q_target;
    }
    double deriv(double T) override {
        return _HUGE;
    }
    double second_deriv(double T) override {
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
    const double a1_sum =
      (core.is_enabled() ? static_cast<double>(core.get_a1()) : 0.0) + (user.is_enabled() ? static_cast<double>(user.get_a1()) : 0.0);
    const double a2_sum =
      (core.is_enabled() ? static_cast<double>(core.get_a2()) : 0.0) + (user.is_enabled() ? static_cast<double>(user.get_a2()) : 0.0);
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
        case iUmolar:
            return 'U';
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
    if (k == 'H' || k == 'S' || k == 'U') {
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
    if (k == 'H' || k == 'S' || k == 'U') {
        const auto stamp = superanc.get_caloric_alpha0_stamp();
        if (stamp.has_value()) {
            const auto [caller_a1, caller_a2] = alpha0_offset_total(HEOS);
            const double cache_a1 = stamp->first;
            const double cache_a2 = stamp->second;
            if (k == 'H' || k == 'U') {
                // h_cache(T) = h_caller(T) + R·T_red·(a2_cache − a2_caller); U
                // inherits the H shift since u = h − p/rho and p/rho is offset-free.
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
                // Newton-step damping: 3 iters (omega = 1.0, 0.6, 0.2)
                // — bounded loop, no accumulation issue worth rewriting.
                for (double omega = 1.0; omega > 0; omega -= increment) {  // NOLINT(cert-flp30-c)
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
        for (const auto& intersection : intersections) {
            if (std::abs(env.Q[intersection.first] - HEOS._Q) < 10 * DBL_EPSILON
                && std::abs(env.Q[intersection.second] - HEOS._Q) < 10 * DBL_EPSILON) {
                solutions.push_back(intersection.first);
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
        for (const auto& intersection : intersections) {
            if (std::abs(env.Q[intersection.first] - 0) < 10 * DBL_EPSILON && std::abs(env.Q[intersection.second] - 0) < 10 * DBL_EPSILON) {
                liquid_solutions.push_back(intersection.first);
            }
            if (std::abs(env.Q[intersection.first] - 1) < 10 * DBL_EPSILON && std::abs(env.Q[intersection.second] - 1) < 10 * DBL_EPSILON) {
                vapor_solutions.push_back(intersection.first);
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
          : HEOS(HEOS), rhomolar_spec(rhomolar_spec), other(other), value(value), Qd(_HUGE) {

            };
        double call(double T) override {
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
    // Two-phase residual for the superancillary "happy path": for a trial
    // saturation temperature T it reads the saturated densities (and pressure)
    // straight off the superancillary, sets the SatL/SatV sub-states, and
    // returns the difference between the vapor quality implied by the density
    // and the vapor quality implied by the caloric input (H/S/U).  Driving
    // this to zero locates the two-phase saturation temperature.  See
    // CoolProp-j3n; ported from the saD_HSU branch.
    class solver_resid_2phase : public FuncWrapper1D
    {
       public:
        // Reference members mirror the sibling residual classes in this file
        // (e.g. HSU_D_flash_twophase::Residual); the object is a short-lived
        // local that never outlives HEOS/superanc.
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
        HelmholtzEOSMixtureBackend& HEOS;
        // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
        EquationOfState::SuperAncillary_t& superanc;
        CoolPropDbl rhomolar_spec;  // Specified density
        parameters other;           // Key for the caloric input (EOS-mode)
        CoolPropDbl value;          // Value of H/S/U in the CALLER reference frame
        // Fast caloric-superancillary mode: read the saturated caloric values
        // straight off the H/S/U caloric superancillary (Chebyshev, ~35 ns)
        // instead of a full EOS property evaluation (~1.3 us).  ca_key is 'H',
        // 'S' or 'U'.  value_cache is `value` shifted into the superancillary's
        // reference frame.  When use_ca is false the residual
        // uses the full EOS (fallback for fluids whose caloric SA is absent).
        bool use_ca;
        char ca_key;
        CoolPropDbl value_cache;
        CoolPropDbl Qd;  // Quality implied by density (output)
        solver_resid_2phase(HelmholtzEOSMixtureBackend& HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value,
                            EquationOfState::SuperAncillary_t& superanc, bool use_ca, char ca_key, CoolPropDbl value_cache)
          : HEOS(HEOS),
            superanc(superanc),
            rhomolar_spec(rhomolar_spec),
            other(other),
            value(value),
            use_ca(use_ca),
            ca_key(ca_key),
            value_cache(value_cache),
            Qd(_HUGE) {};
        double call(double T) override {
            const double rhoL = superanc.eval_sat(T, 'D', 0);
            const double rhoV = superanc.eval_sat(T, 'D', 1);
            // Quality from the specified density
            Qd = (1 / rhomolar_spec - 1 / rhoL) / (1 / rhoV - 1 / rhoL);
            CoolPropDbl Qo;
            if (use_ca) {
                // Saturated caloric values straight off the superancillary (H, S or U),
                // in its (stamped) reference frame; value_cache is the target in the
                // same frame.
                const double yL = superanc.eval_sat(T, ca_key, 0);
                const double yV = superanc.eval_sat(T, ca_key, 1);
                Qo = (value_cache - yL) / (yV - yL);
            } else {
                // Full-EOS fallback.
                const double p = superanc.eval_sat(T, 'P', 1);
                HEOS.SatL->update_TDmolarP_unchecked(T, rhoL, p);
                HEOS.SatV->update_TDmolarP_unchecked(T, rhoV, p);
                const double yL = HEOS.SatL->keyed_output(other);
                const double yV = HEOS.SatV->keyed_output(other);
                Qo = (value - yL) / (yV - yL);
            }
            const double resid = Qo - Qd;
            if (!std::isfinite(resid)) {
                throw ValueError(format("HSU_D superancillary resid not finite @ T=%g K; Qo=%g; Qd=%g", T, Qo, Qd));
            }
            return resid;
        }
    };

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
        // RAII: ensure the imposed gas phase is cleared on every exit
        // path, including exceptions thrown from Halley/Brent.  Without
        // this, the imposed phase leaks into subsequent state queries,
        // which take fast paths that assume gas-phase ideal-gas
        // guesses and converge to metastable roots in the compressed-
        // liquid region (e.g. CO2 LIQUID).  See solver_rho_Tp; #2738.
        ~solver_resid() {
            // This RAII cleanup runs from a destructor (implicitly noexcept):
            // a throw escaping here would call std::terminate.  unspecify_phase()
            // only resets a flag and shouldn't throw, but guard it so a future
            // change can't silently turn this into a crash.
            try {
                HEOS->unspecify_phase();
            } catch (...) {  // NOLINT(bugprone-empty-catch) — a dtor must not let an exception escape
            }
        }
        double call(double T) override {
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
        double deriv(double T) override {
            if (other == iP) {
                return HEOS->first_partial_deriv(other, iT, iDmolar) / value;
            }
            return HEOS->first_partial_deriv(other, iT, iDmolar);
        };
        double second_deriv(double T) override {
            if (other == iP) {
                return HEOS->second_partial_deriv(other, iT, iDmolar, iT, iDmolar) / value;
            }
            return HEOS->second_partial_deriv(other, iT, iDmolar, iT, iDmolar);
        };
        bool input_not_in_range(double T) override {
            return (T < Tmin || T > Tmax);
        }
    };

    if (HEOS.is_pure_or_pseudopure) {

        // ===================== superancillary "happy path" =====================
        // For a pure fluid with a built superancillary and the caloric inputs
        // H/S/U, use the superancillary saturation-density roots to robustly
        // decide single- vs two-phase and to bracket the temperature solve.
        // This avoids the saturation_D_pure instability that makes the legacy
        // routine fail near the critical density, and clamps the quality at the
        // saturation boundary.  Any failure here falls through to the legacy
        // ancillary "sad path" below.  See CoolProp-j3n (ported from saD_HSU).
        //
        // Escape hatch: COOLPROP_DISABLE_SUPERANC_HSU_D (read once) forces the
        // legacy path -- a field kill-switch and an A/B regression-test handle.
        static const bool superanc_hsu_d_disabled = (std::getenv("COOLPROP_DISABLE_SUPERANC_HSU_D") != nullptr);
        if (!superanc_hsu_d_disabled && get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()
            && (other == iHmolar || other == iSmolar || other == iUmolar)) {
            std::shared_ptr<EquationOfState::SuperAncillary_t> sa_ptr = HEOS.get_superanc();
            if (sa_ptr) {
                // Capture the inputs OUTSIDE the try so the catch can restore them.
                // Only read the two members the flash actually consumes -- the density
                // and the one relevant caloric input -- never the others (which are
                // uninitialized for a D+X flash). The solves overwrite these members at
                // trial states, and both the legacy path and later candidate intervals
                // re-read them, so restore them on every non-committing exit.
                const double rho_in = HEOS._rhomolar;
                const CoolPropDbl value = (other == iHmolar) ? HEOS._hmolar : ((other == iSmolar) ? HEOS._smolar : HEOS._umolar);
                auto restore_inputs = [&]() {
                    HEOS._rhomolar = rho_in;
                    if (other == iHmolar) {
                        HEOS._hmolar = value;
                    } else if (other == iSmolar) {
                        HEOS._smolar = value;
                    } else {
                        HEOS._umolar = value;
                    }
                    HEOS.unspecify_phase();
                };
                try {
                    EquationOfState::SuperAncillary_t& superanc = *sa_ptr;

                    // Decide once whether the two-phase residual can read the saturated
                    // caloric values straight off the (fast) caloric superancillary
                    // rather than the EOS.  H, S and U all have superancillaries (built
                    // together; U shares H's reference frame since u = h - p/rho).  The
                    // cached caloric SA is expressed in a stamped reference frame, so
                    // shift the target `value` into that frame by the constant offset
                    // (see resolve_T_via_superancillary / #2773).
                    const char ca_key = (other == iHmolar) ? 'H' : ((other == iSmolar) ? 'S' : 'U');
                    bool use_ca = false;
                    CoolPropDbl value_cache = value;
                    HEOS.ensure_caloric_superancillaries();
                    if (superanc.has_variable(ca_key)) {
                        double delta = 0.0;  // caller -> cache frame shift (constant in T)
                        const auto stamp = superanc.get_caloric_alpha0_stamp();
                        if (stamp.has_value()) {
                            const auto [caller_a1, caller_a2] = alpha0_offset_total(HEOS);
                            const double R = HEOS.gas_constant();
                            if (ca_key == 'S') {
                                delta = -R * (stamp->first - caller_a1);
                            } else {  // 'H' or 'U': u inherits the H enthalpy shift (p/rho is offset-free)
                                const double T_red = HEOS.get_reducing_state().T;
                                delta = R * T_red * (stamp->second - caller_a2);
                            }
                        }
                        value_cache = value + delta;
                        use_ca = true;
                    }

                    const double Tmin_sa = superanc.get_Tmin();
                    const double Tmax_1phase = HEOS.Tmax() * 1.5;
                    const double tol = 1e-12;
                    // 44 correct bits -> tolerance 2^(1-44) ~ 1.14e-13
                    auto eps44 = boost::math::tools::eps_tolerance<double>(44);

                    // Finalize a converged single-phase state at the current (rho, T).
                    auto finalize_1phase = [](HelmholtzEOSMixtureBackend& H) {
                        H._Q = 10000;
                        H._p = H.calc_pressure_nocache(H.T(), H.rhomolar());
                        H.unspecify_phase();
                        H.recalculate_singlephase_phase();
                    };
                    // Final safety gate: the committed state must reproduce BOTH inputs
                    // (density and the caloric value).  A solve can occasionally report
                    // success on a spurious root -- e.g. where X(rho,T) is non-monotonic
                    // at extreme compression and a bracket collapses onto a boundary
                    // temperature.  Refusing such a state here turns a silent wrong
                    // answer into a clean fall-through to the legacy path.
                    auto committed_ok = [&]() -> bool {
                        const double dout = HEOS.rhomolar();
                        if (!std::isfinite(dout) || std::abs(dout / rho_in - 1) > 1e-7) return false;
                        const double xout = HEOS.keyed_output(other);
                        return std::isfinite(xout) && std::abs(xout - value) <= 1e-6 * std::abs(value) + 1e-3;
                    };
                    // True if (rho, T) lies strictly inside the two-phase dome,
                    // i.e. a single-phase root there would be metastable.  Only
                    // meaningful below the critical temperature.
                    auto inside_dome = [&](double T) -> bool {
                        if (T >= superanc.get_Tcrit_num()) return false;
                        const double rhoV = superanc.eval_sat(T, 'D', 1);
                        const double rhoL = superanc.eval_sat(T, 'D', 0);
                        return HEOS._rhomolar > rhoV && HEOS._rhomolar < rhoL;
                    };
                    // Bracketed single-phase solve on [Ta, Tb]; returns true if it
                    // converged to a genuine (non-metastable) single-phase root.
                    auto solve_1phase = [&](double Ta, double Tb) -> bool {
                        try {
                            solver_resid resid(&HEOS, HEOS._rhomolar, value, other, Ta, Tb);
                            const double fa = resid.call(Ta), fb = resid.call(Tb);
                            double Tconv;
                            if (std::abs(fa) < tol) {
                                Tconv = Ta;
                            } else if (std::abs(fb) < tol) {
                                Tconv = Tb;
                            } else if (fa * fb < 0) {
                                boost::math::uintmax_t max_iter = 100;
                                auto f = [&resid](double T) { return resid.call(T); };
                                auto [l, r] = toms748_solve(f, Ta, Tb, fa, fb, eps44, max_iter);
                                Tconv = 0.5 * (l + r);
                            } else {
                                return false;
                            }
                            // Belt-and-suspenders: reject a converged root that lies
                            // inside the dome (metastable).  Interval classification
                            // below should never hand us such a bracket, but the
                            // superancillary/EOS saturation curves differ at the ~1e-8
                            // level, so guard the boundary anyway.
                            if (inside_dome(Tconv)) return false;
                            HEOS.update_DmolarT_direct(HEOS._rhomolar, Tconv);
                            finalize_1phase(HEOS);
                            return true;
                        } catch (const std::exception&) {
                            return false;
                        }
                    };
                    // Bracketed two-phase solve on [Ta, Tb]; returns true if it converged
                    // to a genuine two-phase state.
                    auto solve_2phase = [&](double Ta, double Tb) -> bool {
                        try {
                            solver_resid_2phase resid(HEOS, HEOS._rhomolar, other, value, superanc, use_ca, ca_key, value_cache);
                            const double fa = resid.call(Ta), fb = resid.call(Tb);
                            double Tsol;
                            if (std::abs(fa) < tol) {
                                Tsol = Ta;
                            } else if (std::abs(fb) < tol) {
                                Tsol = Tb;
                            } else if (fa * fb < 0) {
                                boost::math::uintmax_t max_iter = 100;
                                auto f = [&resid](double T) { return resid.call(T); };
                                auto [l, r] = toms748_solve(f, Ta, Tb, fa, fb, eps44, max_iter);
                                Tsol = 0.5 * (l + r);
                            } else {
                                return false;
                            }
                            // EOS polish: the fast path read the saturated caloric
                            // values off the superancillary, which carries a ~1e-8
                            // deviation from the EOS.  T_sol is an excellent seed, so a
                            // couple of secant steps on the EOS-caloric residual (full
                            // EOS h/s/u of each phase) clean the solution up to EOS
                            // precision at a few EOS evaluations -- vs the ~12 the fast
                            // path saved.  U already uses the EOS, so nothing to polish.
                            // The polish is on by default (EOS-exact) but can be turned
                            // off via HSU_D_TWOPHASE_EOS_POLISH for callers that prefer
                            // the faster raw superancillary solution (~1e-8 deviation).
                            double Qd_final;
                            if (use_ca && get_config_bool(HSU_D_TWOPHASE_EOS_POLISH)) {
                                solver_resid_2phase resid_eos(HEOS, HEOS._rhomolar, other, value, superanc, /*use_ca=*/false, ca_key, value_cache);
                                double t0 = Tsol, r0 = resid_eos.call(t0);
                                // Probe just above Tsol, but keep the second secant seed
                                // strictly inside the two-phase bracket [Ta, Tb]: when the
                                // raw solve lands on a boundary root, Tsol*(1+1e-7) can step
                                // past Tb and make resid_eos.call() fail, needlessly kicking
                                // a valid happy-path state back to the legacy path.
                                double t1 = Tsol * (1.0 + 1e-7);
                                if (t1 <= Ta || t1 >= Tb) {
                                    t1 = 0.5 * (Tsol + ((Tsol - Ta) > (Tb - Tsol) ? Ta : Tb));
                                }
                                double r1 = resid_eos.call(t1);
                                for (int it = 0; it < 4; ++it) {
                                    if (!std::isfinite(r1) || std::abs(r1) < tol || std::abs(t1 - t0) <= 1e-13 * t1) break;
                                    // Secant step; a vanishing denominator (r1 == r0) makes tn
                                    // non-finite, which the guard below turns into a clean break
                                    // -- avoids the exact float compare CodeQL flags.
                                    const double tn = t1 - r1 * (t1 - t0) / (r1 - r0);
                                    if (!std::isfinite(tn)) break;
                                    t0 = t1;
                                    r0 = r1;
                                    t1 = tn;
                                    r1 = resid_eos.call(t1);
                                }
                                Tsol = t1;
                                resid_eos.call(Tsol);
                                Qd_final = resid_eos.Qd;
                            } else {
                                resid.call(Tsol);  // refresh Qd at the solution temperature
                                Qd_final = resid.Qd;
                            }
                            // Reject a spurious crossing: if the implied quality is outside
                            // [0, 1] the specified density is NOT within the two-phase band
                            // at Tsol, so this is a single-phase point whose quality residual
                            // merely happened to cross zero.  Clamping Qd here would return a
                            // saturated-boundary density (~rho_crit) instead of the input
                            // density.  Reject so the loop proceeds to the single-phase
                            // interval.  (Dual of the inside_dome guard in solve_1phase.)
                            constexpr double Qeps = 1e-8;
                            if (Qd_final < -Qeps || Qd_final > 1.0 + Qeps) return false;
                            HEOS.update(QT_INPUTS, std::clamp(Qd_final, 0.0, 1.0), Tsol);
                            return true;
                        } catch (const std::exception&) {
                            return false;
                        }
                    };

                    // Reliable saturation temperatures at the specified density.
                    // These are the ONLY temperatures at which the density crosses
                    // the saturation curve, so dome membership is constant between
                    // consecutive roots.  Build candidate intervals
                    // [Tmin_sa .. roots .. Tcrit .. Tmax], classify each by testing
                    // its midpoint with inside_dome(), and run the matching solve.
                    // This never hands the single-phase solver a dome-interior
                    // bracket, and it handles uniformly: the normal topology (one
                    // root), the rho ~ rho_crit collapse (roots pile up at Tcrit, so
                    // the whole sub-critical span is one two-phase interval), and the
                    // water/heavy-water liquid density anomaly (a non-monotonic
                    // liquid branch yields multiple liquid-side roots).
                    auto Tsats = superanc.get_all_intersections('D', HEOS._rhomolar, 48, 100, 1e-13);
                    std::sort(Tsats.begin(), Tsats.end());

                    const double Tcrit = superanc.get_Tcrit_num();
                    // rhoV and rhoL coalesce at Tcrit and the quality ratio degenerates
                    // to 0/0, so keep two-phase brackets just shy of it.
                    const double Tcrit_2phase = Tcrit - std::max(1e-6, 1e-9 * Tcrit);

                    std::vector<double> edges{Tmin_sa};
                    for (const auto& pr : Tsats) {
                        if (pr.first > Tmin_sa && pr.first < Tcrit) edges.push_back(pr.first);
                    }
                    edges.push_back(Tcrit);
                    edges.push_back(Tmax_1phase);

                    for (std::size_t i = 0; i + 1 < edges.size(); ++i) {
                        const double a = edges[i], b = edges[i + 1];
                        if (b - a < 1e-10) continue;
                        const double mid = 0.5 * (a + b);
                        // Clamp the two-phase upper bound just shy of Tcrit, but if that
                        // collapses the bracket (ub <= a, i.e. the last saturation
                        // intersection lands inside the critical guard band) fall back to
                        // the single-phase solve rather than hand solve_2phase an inverted
                        // interval -- which would fail immediately, exactly the near-critical
                        // case this path is meant to recover.
                        bool solved;
                        if (mid < Tcrit && inside_dome(mid)) {
                            const double ub = std::min(b, Tcrit_2phase);
                            solved = (ub > a) ? solve_2phase(a, ub) : solve_1phase(a, b);
                        } else {
                            solved = solve_1phase(a, b);
                        }
                        if (solved) {
                            if (committed_ok()) return;
                            // Converged but did not reproduce the inputs (a spurious root);
                            // undo the state mutation before trying the next interval.
                            restore_inputs();
                        }
                    }
                    throw ValueError("HSU_D superancillary: no candidate interval reproduced the inputs");
                } catch (const std::exception& e) {
                    if (get_debug_level() > 0) {
                        std::cout << "HSU_D_flash superancillary path failed (" << e.what() << "); using legacy path\n";
                    }
                    // Restore the inputs (the happy path may have mutated HEOS), then
                    // fall through to the legacy ancillary path below.
                    restore_inputs();
                }
            }
        }
        // ===================== legacy ancillary "sad path" =====================

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
                } catch (const CoolPropBaseError&) {
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
                // value <= y means the requested state sits below the triple-point
                // isochore (T < Ttriple).  Compressed liquid still exists there
                // above the melting curve (e.g. water down to ~251 K near the
                // 209 MPa melting minimum; IAPWS-95 stays valid), so extend the
                // single-phase solve below the triple point to the coldest point
                // on the melting line rather than rejecting outright.  Only fluids
                // with a melting line have a defined sub-triple liquid region
                // (CoolProp-lgk).
                if (!HEOS.has_melting_line()) {
                    throw ValueError(format("D < DLtriple %g %g", value, y));
                }
                // Coldest liquid temperature = the global minimum of the melting
                // curve.  The curve can be non-monotonic in p (water dips to
                // ~251 K near 209 MPa) and its minimum can sit at a non-smooth
                // segment corner, so neither the endpoint-based
                // MeltingLineVariables::Tmin nor a single coarse scan pinpoints
                // it: sweep coarsely in log(p), then refine linearly around the
                // best point.  For normal fluids whose melting temperature rises
                // monotonically with pressure this stays at ~Ttriple, so the sign
                // check below still rejects (correctly -- they have no sub-triple
                // liquid region).
                MeltingLineVariables& ml = component.ancillaries.melting_line;
                auto Tmelt_at = [&](CoolPropDbl p) -> CoolPropDbl {
                    try {
                        return HEOS.melting_line(iT, iP, p);
                    } catch (...) {
                        return _HUGE;  // p outside this curve's valid range; ignore this sample
                    }
                };
                CoolPropDbl T_lo = ml.Tmin;
                if (ml.pmin > 0 && ml.pmax > ml.pmin) {
                    const int Ncoarse = 100;
                    CoolPropDbl p_best = ml.pmin;
                    for (int i = 0; i <= Ncoarse; ++i) {
                        const CoolPropDbl p_scan = ml.pmin * std::pow(ml.pmax / ml.pmin, static_cast<CoolPropDbl>(i) / Ncoarse);
                        const CoolPropDbl Tm = Tmelt_at(p_scan);
                        if (Tm < T_lo) {
                            T_lo = Tm;
                            p_best = p_scan;
                        }
                    }
                    // Refine within +/- one coarse log-step of the best pressure.
                    const CoolPropDbl step = std::pow(ml.pmax / ml.pmin, 1.0 / Ncoarse);
                    const CoolPropDbl p_a = (std::max)(ml.pmin, p_best / step);
                    const CoolPropDbl p_b = (std::min)(ml.pmax, p_best * step);
                    const int Nfine = 200;
                    for (int i = 0; i <= Nfine; ++i) {
                        const CoolPropDbl p_scan = p_a + (p_b - p_a) * static_cast<CoolPropDbl>(i) / Nfine;
                        const CoolPropDbl Tm = Tmelt_at(p_scan);
                        if (Tm < T_lo) {
                            T_lo = Tm;
                        }
                    }
                }
                // y_lo is the caloric value on the requested isochore at the
                // coldest liquid temperature.  Using the curve's global minimum
                // (rather than the melting point at this exact density) makes the
                // floor conservative on rejection: states below it are genuinely
                // solid, while a thin metastable sliver just below the per-density
                // melting point may still be accepted -- it round-trips correctly.
                CoolPropDbl y_lo = NAN;
                switch (other) {
                    case iSmolar:
                        y_lo = HEOS.calc_smolar_nocache(T_lo, HEOS._rhomolar);
                        break;
                    case iHmolar:
                        y_lo = HEOS.calc_hmolar_nocache(T_lo, HEOS._rhomolar);
                        break;
                    case iUmolar:
                        y_lo = HEOS.calc_umolar_nocache(T_lo, HEOS._rhomolar);
                        break;
                    default:  // iP
                        y_lo = HEOS.calc_pressure_nocache(T_lo, HEOS._rhomolar);
                        break;
                }
                if (!ValidNumber(y_lo) || value < y_lo) {
                    // Colder than the coldest liquid attainable at this density:
                    // genuinely below the melting line (solid / out of range).
                    throw ValueError(format("D+X below melting line [other=%d]: %g < %g at Tmelt_min %g K", other, value, y_lo, T_lo));
                }
                solver_resid resid(&HEOS, HEOS._rhomolar, value, other, T_lo, TLtriple);
                HEOS._phase = iphase_liquid;
                CoolPropDbl T_converged;
                try {
                    T_converged = Halley(resid, 0.5 * (T_lo + TLtriple), DBL_EPSILON, 100);
                } catch (...) {
                    T_converged = Brent(resid, T_lo, TLtriple, DBL_EPSILON, 1e-12, 100);
                }
                // Re-evaluate at converged (rho, T) so alphar cache is
                // consistent with the final state (#1907).
                HEOS.unspecify_phase();
                HEOS.update_DmolarT_direct(HEOS._rhomolar, T_converged);
                HEOS._phase = iphase_liquid;
                HEOS._Q = 10000;
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
    // Save the target pressure at function entry.  HEOS._p can be CORRUPTED by
    // the inner density solve (update(PT_INPUTS, p, T)) if it throws part-way
    // through its Newton iteration -- it can be left holding a garbage value
    // (e.g. a negative pressure).  The exception-fallback guard below must test
    // the *intended* pressure, not this possibly-corrupted live value, otherwise
    // it branches on garbage and re-throws instead of running the 2-D Newton
    // fallback.  Always read p_target (not HEOS._p) in the fallback guard.
    const CoolPropDbl p_target = HEOS._p;

    // Env toggles for the direct-EOS cache-bypass path (CoolProp-cdj).
    //
    // Default: ON for pure fluids — bit-equivalent within ULP to master, ~17%
    // faster end-to-end on the PXcdj_sweep, validated against 362,004 grid
    // points with zero new failures.  Opt out by setting PXFLASH_DIRECT_EOS=0.
    //
    //   PXFLASH_DIRECT_EOS=0   - disable the direct path; uses master's full
    //                            cached path.  All other values (incl. unset)
    //                            keep the default-on direct path.
    //                            (Cold probes always keep master's
    //                            update(PT_INPUTS) so the SRK + Householder4
    //                            basin classification is preserved exactly.)
    //   PXFLASH_INNER_NEWTON   - OFF (default) = Householder3 inner solver;
    //                            ON           = plain Newton (no 2nd/3rd derivs).
    //
    // is_pure() gate excludes mixtures and pseudo-pure here -- the property
    // reconstruction formulas below are written for a single residual block.
    static const bool pxflash_direct_eos = []() {
        const char* env = std::getenv("PXFLASH_DIRECT_EOS");
        if (env == nullptr) return true;
        return std::string(env) != "0";
    }();
    static const bool pxflash_inner_newton = (std::getenv("PXFLASH_INNER_NEWTON") != nullptr);
    const bool direct_path_enabled = pxflash_direct_eos && HEOS.is_pure();

    class solver_resid : public FuncWrapper1DWithTwoDerivs
    {
       public:
        HelmholtzEOSMixtureBackend* HEOS;
        CoolPropDbl p, value;
        parameters other;
        int iter;
        CoolPropDbl eos0, eos1, rhomolar, rhomolar0, rhomolar1;
        CoolPropDbl Tmin, Tmax;
        bool direct_path_enabled;

        solver_resid(HelmholtzEOSMixtureBackend* HEOS, CoolPropDbl p, CoolPropDbl value, parameters other, double Tmin, double Tmax,
                     bool direct_path_enabled)
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
            Tmax(Tmax),
            direct_path_enabled(direct_path_enabled) {
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
        // Inner-density Newton: solve p(T, rho) = p_target using only
        // residual_helmholtz->all (bypasses CachedElement).  Derivatives
        // match SolverTPResid::deriv/second_deriv/third_deriv exactly.
        // Throws on failure so the caller can fall back to the master path.
        //
        // ar_out / a0_out: populated with the HelmholtzDerivatives (residual
        // and ideal-gas) evaluated AT the converged rho.  Loop checks
        // convergence on (rho, ar) BEFORE taking the step so ar_out is never
        // stale; a0_out is computed once at the converged rho immediately
        // before return.  This makes direct_density_newton the sole producer
        // of derivative sets and lets direct_property_eval be a pure consumer
        // (no ->all(), no calc_all_alpha0_derivs_nocache).
        static double direct_density_newton(HelmholtzEOSMixtureBackend& H, const std::vector<CoolPropDbl>& x, double Tr, double rhor, double R,
                                            double T, double p_target, double rho_initial, bool use_plain_newton, HelmholtzDerivatives& ar_out,
                                            HelmholtzDerivatives& a0_out) {
            const double tau = Tr / T;
            double rho = rho_initial;
            for (int it = 0; it < 15; ++it) {
                const double delta = rho / rhor;
                ar_out = H.residual_helmholtz->all(H, x, tau, delta, false);
                const double alphar_d = ar_out.dalphar_ddelta;
                const double alphar_dd = ar_out.d2alphar_ddelta2;
                const double p = rho * R * T * (1.0 + delta * alphar_d);
                const double f = p - p_target;
                const double dpdrho = R * T * (1.0 + 2.0 * delta * alphar_d + delta * delta * alphar_dd);
                if (dpdrho <= 0.0) {
                    // Mechanically unstable region (inside spinodal); the cached path's
                    // solver_rho_Tp has phase-imposed basin-recovery for this case.
                    // Throw to fall back to master's cached path via the catch in call().
                    throw ValueError("direct-EOS Newton: dpdrho <= 0 (unstable region)");
                }
                // Residual convergence: ar_out matches the returned rho exactly.
                if (std::abs(f) < 1e-12 * std::abs(p_target)) {
                    a0_out = H.calc_all_alpha0_derivs_nocache(x, tau, delta, Tr, rhor);
                    return rho;
                }
                double dx;
                if (use_plain_newton) {
                    if (!std::isfinite(dpdrho) || std::abs(dpdrho) < 1e-300) {
                        throw ValueError("direct-EOS density Newton: singular dpdrho");
                    }
                    dx = -f / dpdrho;
                } else {
                    const double alphar_ddd = ar_out.d3alphar_ddelta3;
                    const double alphar_dddd = ar_out.d4alphar_ddelta4;
                    const double d2pdrho2 = (R * T / rhor) * (2.0 * alphar_d + 4.0 * delta * alphar_dd + delta * delta * alphar_ddd);
                    const double d3pdrho3 = (R * T / (rhor * rhor)) * (6.0 * alphar_dd + 6.0 * delta * alphar_ddd + delta * delta * alphar_dddd);
                    // Householder3 step (same form as solvers.cpp Householder4 body)
                    const double num = f * (dpdrho * dpdrho - f * d2pdrho2 / 2.0);
                    const double den = dpdrho * dpdrho * dpdrho - f * dpdrho * d2pdrho2 + d3pdrho3 * f * f / 6.0;
                    dx = (!std::isfinite(den) || std::abs(den) < 1e-300) ? (-f / dpdrho) : (-num / den);
                }
                if (!std::isfinite(dx)) throw ValueError("direct-EOS density Newton: non-finite step");
                double scale = 1.0;
                while (rho + scale * dx <= 0.0 && scale > 1e-6)
                    scale *= 0.5;
                const double rho_new = rho + scale * dx;
                if (!std::isfinite(rho_new) || rho_new <= 0.0) throw ValueError("direct-EOS density Newton: non-finite rho");
                // Step-size convergence: the step was taken (rho moved), so the
                // ar_out from BEFORE the step is now stale w.r.t. rho_new.
                // Refresh ar_out AND a0_out at rho_new before returning so the
                // caller can safely consume them.  This is the rare path (most
                // converge on |f| < tol above); one extra ->all() + one extra
                // alpha0 vs. master's ~9 ->all()s saved.
                if (std::abs(scale * dx) / rho_new < 1e-12) {
                    rho = rho_new;
                    const double delta_new = rho / rhor;
                    ar_out = H.residual_helmholtz->all(H, x, tau, delta_new, false);
                    a0_out = H.calc_all_alpha0_derivs_nocache(x, tau, delta_new, Tr, rhor);
                    return rho;
                }
                rho = rho_new;
            }
            throw ValueError("direct-EOS density Newton: max iters");
        }
        // Direct property reconstruction at (T, rho) using HelmholtzDerivatives
        // sets supplied by direct_density_newton at convergence.  Matches
        // calc_hmolar_nocache / calc_smolar_nocache / calc_umolar_nocache exactly.
        // Pure consumer: no ->all() and no calc_all_alpha0_derivs_nocache here;
        // direct_density_newton is the sole producer of derivative sets so the
        // cost-model accounting stays honest (one ar + one a0 per outer probe).
        static double direct_property_eval(HelmholtzEOSMixtureBackend& H, const std::vector<CoolPropDbl>& x, double Tr, double rhor, double R,
                                           double T, double rho, parameters other, const HelmholtzDerivatives& ar, const HelmholtzDerivatives& a0) {
            const double tau = Tr / T;
            const double delta = rho / rhor;
            const double atau_total = a0.dalphar_dtau + ar.dalphar_dtau;
            const double a_total = a0.alphar + ar.alphar;
            switch (other) {
                case iHmolar:
                    return R * T * (1.0 + tau * atau_total + delta * ar.dalphar_ddelta);
                case iSmolar:
                    return R * (tau * atau_total - a_total);
                case iUmolar:
                    return R * T * tau * atau_total;
                default:
                    throw ValueError("direct-EOS property: unsupported `other` key");
            }
        }
        double call(double T) override {
            // Direct-path WARM branch only: cold probes go through master's
            // update(PT_INPUTS) so the SRK seed + Householder4 + basin-recovery
            // (iphase_liquid/gas retries inside solver_rho_Tp) are byte-identical
            // to master.  iter < 2 cold iters are amortized; the win is from
            // bypassing the cached inner Householder4 on the warm probes that
            // dominate TOMS748's iteration count.
            const bool want_warm = (iter >= 2 && std::abs(rhomolar1 / rhomolar0 - 1) <= 0.05);
            if (direct_path_enabled && want_warm) {
                try {
                    // WARM: bypass the cached inner solve; reuse the carried
                    // rhomolar (identical to what update_TP_guessrho would seed).
                    // ar AND a0 are filled by direct_density_newton at the
                    // converged rho and reused by direct_property_eval -- saves
                    // one ->all() + one alpha0 call per outer probe at the same
                    // (tau, delta).  direct_density_newton is the sole producer;
                    // direct_property_eval is a pure consumer (symmetry).
                    HelmholtzDerivatives ar, a0;
                    const double rho_solved =
                      direct_density_newton(*HEOS, HEOS->get_mole_fractions_ref(), HEOS->get_reducing_state().T, HEOS->get_reducing_state().rhomolar,
                                            HEOS->gas_constant(), T, p, rhomolar, pxflash_inner_newton, ar, a0);
                    CoolPropDbl eos = direct_property_eval(*HEOS, HEOS->get_mole_fractions_ref(), HEOS->get_reducing_state().T,
                                                           HEOS->get_reducing_state().rhomolar, HEOS->gas_constant(), T, rho_solved, other, ar, a0);
                    // Commit (T, rho, p) to the HEOS cache so deriv/second_deriv
                    // (used by Halley) and any downstream readers see the
                    // converged state.  We use update_TDmolarP_unchecked
                    // instead of update_DmolarT_direct because the latter
                    // re-derives p from (T, rho) via calc_pressure() (one more
                    // ->dalphar_dDelta()), whereas here p == p_target by
                    // construction of the Newton convergence test.
                    HEOS->update_TDmolarP_unchecked(T, rho_solved, p);
                    rhomolar = rho_solved;

                    if (verbosity > 0 && iter == 0) {
                        std::cout << format("T: %0.15g rho: %0.15g eos: %0.15g", T, rhomolar, eos);
                    }
                    CoolPropDbl r = eos - value;
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
                } catch (...) {
                    // Defensive per-probe fallback to the master cached path.
                    // fall through
                }
            }

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
        double deriv(double T) override {
            return HEOS->first_partial_deriv(other, iT, iP);
        }
        double second_deriv(double T) override {
            return HEOS->second_partial_deriv(other, iT, iP, iT, iP);
        }
        bool input_not_in_range(double x) override {
            return (x < Tmin || x > Tmax);
        }
    };
    solver_resid resid(&HEOS, HEOS._p, value, other, Tmin, Tmax, direct_path_enabled);

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
            // TOMS748 correct-bits tolerance is 2^(1-bits).  44 bits (~1.1e-13)
            // is far tighter than thermophysical properties need; 30 bits
            // (~9.3e-10 on T) keeps the worst P+X round-trip well under the 1e-6
            // hard limit while cutting ~7% of the outer iterations.
            // >>> 2**(1-30) == 9.31e-10
            auto [l, r] = toms748_solve(f, Tmin, Tmax, resid_Tmin, resid_Tmax, boost::math::tools::eps_tolerance<double>(30), max_iter);
            // Re-evaluate at the midpoint of the final bracket so HEOS is left at
            // the converged state (TOMS748's last probe may be at a slightly
            // different point within the bracket).  Same discipline as the
            // supercritical-cold narrow-TOMS748 path below.
            resid.iter = 0;
            resid.call(0.5 * (l + r));
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
            // Decide whether to attempt the 2-D Newton fallback.  IMPORTANT: test
            // p_target (saved at entry) NOT HEOS._p -- the inner density solve can
            // leave HEOS._p corrupted when it throws mid-Newton, and branching on
            // that garbage pressure is what caused supercritical-cold P+X states to
            // re-throw instead of recovering (see the p_target comment above).
            //   - near-critical (0.95*pc < p <= pc): the original fallback case.
            //   - supercritical (p > pc): compressed-liquid corner where TOMS748
            //     can't bracket zero across the wide T range; engage the fallback
            //     here too (this is the Nitrogen/PS fix).
            const bool in_near_critical = (p_target > 0.95 * p_critical_ && p_target <= p_critical_);
            const bool in_supercritical = (p_target > p_critical_);
            if (!in_near_critical && !in_supercritical) {
                throw;
            }
            if (get_debug_level() > 0) {
                std::cout << resid.errstring << '\n';
            }
            if (in_supercritical) {
                // Supercritical-cold strategy.  TOMS748 threw because the EOS has
                // numerical issues at high T inside the wide supercritical bracket;
                // the solution is a compressed-liquid state at T << T_crit.
                // Step 1: try a narrowed TOMS748 on [Tmin, 1.05*T_crit], which
                //   avoids the high-T instability and the non-monotone region.
                // Step 2: if that bracket doesn't sign-change (or throws), fall
                //   back to the 2-D Newton from the high-density Tmin seed.
                const double T_crit = HEOS.T_critical();
                const double Tnarrow = T_crit * 1.05;
                bool narrow_toms748_ok = false;
                if (Tnarrow < Tmax) {
                    resid.iter = 0;
                    double resid_Tnarrow = 0.0;
                    bool have_narrow_bracket = false;
                    try {
                        resid_Tnarrow = resid.call(Tnarrow);
                        have_narrow_bracket = (resid_Tnarrow * resid_Tmin < 0);
                    } catch (const std::exception& e) {
                        if (get_debug_level() > 0) {
                            std::cout << "narrow TOMS748 probe threw: " << e.what() << '\n';
                        }
                        // narrow_toms748_ok stays false -> fall through to 2-D Newton backstop
                    }
                    if (have_narrow_bracket) {
                        resid.iter = 0;
                        boost::math::uintmax_t max_iter_narrow = 100;
                        auto fn = [&resid](const double T) { return resid.call(T); };
                        try {
                            auto [nl, nr] = toms748_solve(fn, Tmin, Tnarrow, resid_Tmin, resid_Tnarrow, boost::math::tools::eps_tolerance<double>(30),
                                                          max_iter_narrow);
                            // Re-evaluate at the midpoint of the final bracket so HEOS
                            // is left at the converged state (TOMS748's last call may
                            // be at a slightly different point).
                            const double T_conv = 0.5 * (nl + nr);
                            resid.iter = 0;
                            resid.call(T_conv);
                            if (is_in_closed_range(Tmin, Tmax, static_cast<CoolPropDbl>(resid.HEOS->T()))) {
                                HEOS.unspecify_phase();
                                HEOS.recalculate_singlephase_phase();
                                narrow_toms748_ok = true;
                            }
                        } catch (const std::exception& e) {
                            // narrow TOMS748 also threw -- intentionally swallow and fall
                            // through to the 2-D Newton fallback below (narrow_toms748_ok
                            // stays false).  Surface the reason only at debug level.
                            if (get_debug_level() > 0) {
                                std::cout << "narrow TOMS748 fallback threw: " << e.what() << '\n';
                            }
                        }
                    }
                }
                if (!narrow_toms748_ok) {
                    // 2-D Newton from the high-density Tmin seed (liquid-like basin).
                    std::vector<double> x0 = {Tmin, rhomolar_Tmin};
                    NDNewtonRaphson_Jacobian(&solver_resid2d, x0, 1e-12, 20, 1.0);
                    // Acceptance check: T in range AND both residuals small.
                    // NDNewtonRaphson_Jacobian can early-exit on a tiny step while
                    // the residual is still large; verify the state actually reproduces
                    // the inputs rather than returning a silently-wrong root.
                    {
                        const double T_sol = static_cast<CoolPropDbl>(solver_resid2d.HEOS->T());
                        const double p_resid_rel = std::abs(HEOS.p() - p_target) / p_target;
                        const double other_scale = std::abs(value) > 1.0 ? std::abs(value) : 1.0;
                        const double other_resid_rel = std::abs(HEOS.keyed_output(other) - value) / other_scale;
                        if (!is_in_closed_range(Tmin, Tmax, T_sol) || p_resid_rel > 1e-6 || other_resid_rel > 1e-6) {
                            throw ValueError("2D Newton method was unable to find a solution in HSU_P_flash_singlephase_Brent");
                        }
                    }
                    HEOS.unspecify_phase();
                    HEOS.recalculate_singlephase_phase();
                }
            } else {
                std::vector<double> x0 = {Tstart, rhomolarstart};
                NDNewtonRaphson_Jacobian(&solver_resid2d, x0, 1e-12, 20, 1.0);
                const double T_sol = static_cast<CoolPropDbl>(solver_resid2d.HEOS->T());
                const double p_resid_rel = std::abs(HEOS.p() - p_target) / p_target;
                const double other_scale = std::abs(value) > 1.0 ? std::abs(value) : 1.0;
                const double other_resid_rel = std::abs(HEOS.keyed_output(other) - value) / other_scale;
                if (!is_in_closed_range(Tmin, Tmax, T_sol) || solver_resid2d.HEOS->phase() != phase || p_resid_rel > 1e-6 || other_resid_rel > 1e-6) {
                    throw ValueError("2D Newton method was unable to find a solution in HSU_P_flash_singlephase_Brent");
                }
                // Un-specify the phase of the fluid
                HEOS.unspecify_phase();
                HEOS.recalculate_singlephase_phase();
            }
        } catch (...) {
            if (get_debug_level() > 0) {
                std::cout << resid.errstring << '\n';
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
        double call(double rhomolar) override {
            HEOS->update_DmolarT_direct(rhomolar, T);
            double eos = HEOS->keyed_output(other);
            return eos - value;
        };
        double deriv(double rhomolar) override {
            return HEOS->first_partial_deriv(other, iDmolar, iT);
        }
        double second_deriv(double rhomolar) override {
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
        // Update the state (T > Tc). Honor an imposed phase here: a caller who used
        // specify_phase(iphase_supercritical_gas/liquid) would otherwise have it
        // silently rewritten by this density-based classifier. (#2718)
        if (HEOS.imposed_phase_index == iphase_not_imposed) {
            if (HEOS._p < HEOS.p_critical()) {
                HEOS._phase = iphase_supercritical_gas;
            } else {
                HEOS._phase = iphase_supercritical;
            }
        }
    }
    // Subcritical temperature liquid
    else if ((HEOS._phase == iphase_liquid) || (HEOS._phase == iphase_supercritical_liquid)) {
        CoolPropDbl ymelt = NAN, yL = NAN, y = NAN;
        CoolPropDbl rhomelt = HEOS.components[0].triple_liquid.rhomolar;
        // Self-contained seed evaluation: prefer the cached ancillary set by
        // upstream phase-determination, otherwise fall through to superancillary
        // (numerical-critical-point accurate), then to the polynomial ancillary.
        // Lets DHSU_T_flash skip the unconditional ancillary evaluation on the
        // imposed-phase fast path without breaking this consumer (#2718).
        CoolPropDbl rhoL = NAN;
        if (HEOS._rhoLanc) {
            rhoL = static_cast<double>(HEOS._rhoLanc);
        } else if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
            auto superanc_ptr = HEOS.get_superanc();
            if (superanc_ptr) {
                // get_superanc() only reports cache existence, not domain coverage —
                // eval_sat() can still throw near the numerical critical point or a
                // domain edge. Swallow the throw so the polynomial fallback below stays
                // reachable for previously valid subcritical flashes.
                try {
                    rhoL = superanc_ptr->eval_sat(HEOS._T, 'D', 0);
                } catch (...) {
                    rhoL = NAN;
                }
            }
        }
        if (!ValidNumber(rhoL)) {
            rhoL = HEOS.components[0].ancillaries.rhoL.evaluate(HEOS._T);
        }

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
        // See companion block in the liquid branch above (#2718).
        CoolPropDbl rhoV = NAN;
        if (HEOS._rhoVanc) {
            rhoV = static_cast<double>(HEOS._rhoVanc);
        } else if (get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
            auto superanc_ptr = HEOS.get_superanc();
            if (superanc_ptr) {
                // See companion try/catch in the liquid branch above.
                try {
                    rhoV = superanc_ptr->eval_sat(HEOS._T, 'D', 1);
                } catch (...) {
                    rhoV = NAN;
                }
            }
        }
        if (!ValidNumber(rhoV)) {
            rhoV = HEOS.components[0].ancillaries.rhoV.evaluate(HEOS._T);
        }

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
            if (HEOS._phase == iphase_liquid) {
                HEOS._Q = -1000;
            } else if (HEOS._phase == iphase_gas) {
                HEOS._Q = 1000;
            } else if (HEOS._phase == iphase_twophase) {
                // The ancillary rhoL/rhoV densities are only needed by the twophase
                // sub-branch (and only when other != iDmolar — they seed SatL/SatV
                // updates). Keep the evaluation localized so the cheaper single-phase
                // imposed-phase paths (liquid, gas, supercritical_liquid) avoid the
                // two polynomial evaluations on every update() (#2718).
                // TODO: is it a bug that this branch can be accessed for mixtures?
                HEOS._rhoVanc = HEOS.components[0].ancillaries.rhoV.evaluate(HEOS._T);
                HEOS._rhoLanc = HEOS.components[0].ancillaries.rhoL.evaluate(HEOS._T);
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
    if (HEOS.is_pure_or_pseudopure && HEOS._phase != iphase_twophase && HEOS.imposed_phase_index == iphase_not_imposed) {
        // Update the state for conditions where the state was guessed. When the user
        // imposed a phase via specify_phase(), honor it: recalculate_singlephase_phase()
        // would silently overwrite the caller's choice (e.g. flip iphase_gas to
        // iphase_liquid above rhomolar_critical), and the overwrite is also wasted
        // work on the imposed-phase fast path (#2718).
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
        double call(double T) override {
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
        // Geometric halving of the Newton step length (10 iters from
        // 1.0 down to ~1/1024) — geometric, not summation; no FP
        // counter-accumulation issue.
        for (double frac = 1.0; frac > 0.001; frac /= 2) {  // NOLINT(cert-flp30-c)
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
namespace {
// ===========================================================================
// Single-phase H,S superancillary "happy path" (the cascade).  Ported from the
// validated prototype in src/Tests/CoolProp-Tests-HS.cpp (bd CoolProp-j3n.4).
//
// Three legs, each a homotopy that anchors where a state of the target entropy
// provably exists and then marches in (T, ln rho) holding (h,s) on a chord:
//   1) saturation-anchor   : anchor = saturated state whose entropy == s_t
//   2) supercritical isentrope : anchor on the T=Tmax isotherm (dome-free)
//   3) ideal-gas departure : anchor at lambda=0 (no dome at all)
// A candidate is accepted only if it reproduces (h,s) AND is fully intrinsically
// stable (dp/drho|_T > 0 AND cv > 0) within [Tmin,Tmax]; that makes the stable
// (h,s)->(T,rho) root unique, so no knowledge of the true (T,rho) is needed.
// See the docs notebook (Web/coolprop/HSFlash.ipynb) for the full derivation.
// ===========================================================================
using HS_SA_t = EquationOfState::SuperAncillary_t;

// RAII: impose single phase so update_DmolarT_direct + h/s/derivative queries
// take the single-phase EOS path everywhere (including metastable iterates).
struct HSGasGuard
{
    HelmholtzEOSMixtureBackend& H;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)
    explicit HSGasGuard(HelmholtzEOSMixtureBackend& H) : H(H) {
        H.specify_phase(iphase_gas);
    }
    // unspecify_phase() only resets an enum flag and does not throw.
    // NOLINTNEXTLINE(bugprone-exception-escape)
    ~HSGasGuard() {
        H.unspecify_phase();
    }
    HSGasGuard(const HSGasGuard&) = delete;
    HSGasGuard& operator=(const HSGasGuard&) = delete;
    HSGasGuard(HSGasGuard&&) = delete;
    HSGasGuard& operator=(HSGasGuard&&) = delete;
};

// Shared (T, w=ln rho) homotopy corrector: homotope (h,s) linearly from the
// anchor's values to the target with adaptive subdivision; the dp/drho|_T>0 guard
// keeps the corrector on the mechanically stable branch.  Caller imposes the gas
// phase.  Returns true and fills T_out/rho_out on success.
bool hs_corrector(HelmholtzEOSMixtureBackend& H, double T0, double rho0, double h_t, double s_t, double& T_out, double& rho_out,
                  double Tlo_override = -1.0) {
    const double Rgas = H.gas_constant(), Tsc = H.T_critical();
    const double hscale = Rgas * Tsc, sscale = Rgas;  // dimensional scales (h passes through 0 at the ref state)
    const double Tlo = (Tlo_override > 0.0) ? Tlo_override : H.Tmin() * (1.0 - 2e-2);
    const double Thi = 1.5 * H.Tmax();  // 2% sub-Tmin slack default; melting leg passes a lower floor
    auto eval = [&](double T, double rho) { H.update_DmolarT_direct(rho, T); };
    eval(T0, rho0);
    const double h0 = H.hmolar(), s0 = H.smolar();
    for (int N = 1; N <= 128; N *= 2) {
        double T = T0, w = std::log(rho0);
        bool failed = false;
        for (int k = 1; k <= N && !failed; ++k) {
            const double lam = static_cast<double>(k) / N;
            const double ht = h0 + lam * (h_t - h0), st = s0 + lam * (s_t - s0);
            bool conv = false;
            double best_norm = 1e300, best_T = T, best_w = w;
            for (int iter = 0; iter < 40; ++iter) {
                const double rho = std::exp(w);
                eval(T, rho);
                const double rh = H.hmolar() - ht, rs = H.smolar() - st;
                const double norm = std::abs(rh) / hscale + std::abs(rs) / sscale;
                if (norm < best_norm) {
                    best_norm = norm;
                    best_T = T;
                    best_w = w;
                }
                if (norm < 1e-11) {
                    conv = true;
                    break;
                }
                const double hT = H.first_partial_deriv(iHmolar, iT, iDmolar), hr = H.first_partial_deriv(iHmolar, iDmolar, iT);
                const double sT = H.first_partial_deriv(iSmolar, iT, iDmolar), sr = H.first_partial_deriv(iSmolar, iDmolar, iT);
                const double a11 = hT, a12 = rho * hr, a21 = sT, a22 = rho * sr;
                const double det = a11 * a22 - a12 * a21;
                if (!std::isfinite(det) || std::abs(det) < 1e-300) break;
                const double dT = -(a22 * rh - a12 * rs) / det, dw = -(-a21 * rh + a11 * rs) / det;
                double f = 1.0;
                while ((T + f * dT < Tlo || T + f * dT > Thi || std::abs(f * dw) > 2.0) && f > 1e-6)
                    f *= 0.5;
                for (int g = 0; g < 8 && f > 1e-3; ++g) {  // stay on the mechanically stable branch
                    eval(T + f * dT, std::exp(w + f * dw));
                    if (H.first_partial_deriv(iP, iDmolar, iT) > 0) break;
                    f *= 0.5;
                }
                T += f * dT;
                w += f * dw;
            }
            if (!conv && best_norm < 1e-8) {
                T = best_T;
                w = best_w;
                conv = true;
            }
            if (!conv) failed = true;
        }
        if (!failed) {
            T_out = T;
            rho_out = std::exp(w);
            return true;
        }
    }
    return false;
}

// Target entropy shifted into the caloric superancillary's stamped frame (#2773).
double hs_s_to_cache(HelmholtzEOSMixtureBackend& H, HS_SA_t& sa, double s_t) {
    const auto stamp = sa.get_caloric_alpha0_stamp();
    if (!stamp.has_value()) return s_t;
    const double a1 = FlashRoutines::alpha0_offset_total(H).first;
    return s_t - H.gas_constant() * (stamp->first - a1);
}
// Target enthalpy shifted into the caloric superancillary's stamped frame (#2773).
double hs_h_to_cache(HelmholtzEOSMixtureBackend& H, HS_SA_t& sa, double h_t) {
    const auto stamp = sa.get_caloric_alpha0_stamp();
    if (!stamp.has_value()) return h_t;
    const double a2 = FlashRoutines::alpha0_offset_total(H).second;
    return h_t + H.gas_constant() * H.get_reducing_state().T * (stamp->second - a2);
}

// Leg 1: saturation anchor.  Anchor at the saturated state whose entropy equals
// the target's (found by the superancillary's native S-inversion), then march.
bool hs_leg_saturation(HelmholtzEOSMixtureBackend& H, HS_SA_t& sa, double h_t, double s_t, double& T_out, double& rho_out) {
    HSGasGuard guard(H);
    const double Tc = sa.get_Tcrit_num(), Tmin = sa.get_Tmin();
    const double Tcap = Tc - std::max(1e-4, 1e-6 * Tc);
    const double a = Tmin + 1e-3, b = Tcap;
    H.ensure_caloric_superancillaries();
    const double s_cache = hs_s_to_cache(H, sa, s_t);
    // Disambiguate the saturation branch by enthalpy closeness to the target.  When
    // the caloric H superancillary is present, read the saturated enthalpy straight
    // off it (Chebyshev, no EOS call) instead of evaluating the EOS at each
    // candidate -- this is the dominant leg, so it trims EOS evaluations off the
    // common case (and several for retrograde fluids with multiple roots).
    const bool use_ca_h = sa.has_variable('H');
    const double h_cache = use_ca_h ? hs_h_to_cache(H, sa, h_t) : h_t;
    double T0 = 0, rho0 = 0, best_hgap = 1e300;
    bool found = false, near_critical = false;
    auto consider = [&](double Tcand, short Q) {
        const double rho = sa.eval_sat(Tcand, 'D', Q);
        double hgap;
        if (use_ca_h) {
            hgap = std::abs(sa.eval_sat(Tcand, 'H', Q) - h_cache);
        } else {
            H.update_DmolarT_direct(rho, Tcand);
            hgap = std::abs(H.hmolar() - h_t);
        }
        if (std::isfinite(hgap) && hgap < best_hgap) {
            best_hgap = hgap;
            T0 = Tcand;
            rho0 = rho;
            found = true;
        }
    };
    for (const auto& pr : sa.get_all_intersections('S', s_cache, 48, 100, 1e-12)) {
        const double Tcand = pr.first;
        if (Tcand >= b) near_critical = true;  // anchor pinned against Tc => effectively critical/supercritical entropy
        if (Tcand <= a || Tcand >= b) continue;
        const short Qroot = (std::abs(sa.eval_sat(Tcand, 'S', 0) - s_cache) <= std::abs(sa.eval_sat(Tcand, 'S', 1) - s_cache)) ? 0 : 1;
        consider(Tcand, Qroot);
    }
    if (!found) {
        if (near_critical) return false;  // hand to the supercritical isentrope leg
        const double s_crit = [&] {
            H.update_DmolarT_direct(sa.eval_sat(Tcap, 'D', 0), Tcap);
            return H.smolar();
        }();
        if (s_t > s_crit) {  // dilute gas: T from the monotone 1D enthalpy inversion
            rho0 = sa.get_rhocrit_num() * 1e-4;
            auto hres = [&](double T) {
                H.update_DmolarT_direct(rho0, T);
                return H.hmolar() - h_t;
            };
            double Ta = Tmin + 1e-3, Tb = 1.5 * H.Tmax(), fa = hres(Ta), fb = hres(Tb);
            if (fa * fb <= 0) {
                boost::math::uintmax_t it = 60;
                auto [l, r] = toms748_solve(hres, Ta, Tb, fa, fb, boost::math::tools::eps_tolerance<double>(30), it);
                T0 = 0.5 * (l + r);
            } else {
                T0 = (std::abs(fa) <= std::abs(fb)) ? Ta : Tb;
            }
        } else {  // cold compressed liquid: triple-point saturated liquid
            T0 = a;
            rho0 = sa.eval_sat(a, 'D', 0);
        }
    }
    return hs_corrector(H, T0, rho0, h_t, s_t, T_out, rho_out);
}

// Leg 2: supercritical isentrope.  Anchor on the dome-free T=Tmax isotherm at the
// density whose entropy equals s_t, then march (holding s == s_t).
bool hs_leg_isentrope(HelmholtzEOSMixtureBackend& H, double h_t, double s_t, double& T_out, double& rho_out) {
    HSGasGuard guard(H);
    const double Ta = H.Tmax(), rhoc = H.rhomolar_critical();
    auto sres = [&](double rho) {
        H.update_DmolarT_direct(rho, Ta);
        return H.smolar() - s_t;
    };
    const double ra = 1e-6 * rhoc, ga = sres(ra);
    double rb = 6.0 * rhoc, gb = sres(rb);
    for (int e = 0; e < 12 && ga * gb > 0; ++e) {
        const double rb_next = rb * 2.0;
        double gb_next;
        try {
            gb_next = sres(rb_next);
        } catch (...) {
            break;
        }
        if (!std::isfinite(gb_next)) break;
        rb = rb_next;
        gb = gb_next;
    }
    double rho0;
    if (ga * gb > 0) {
        rho0 = (std::abs(ga) <= std::abs(gb)) ? ra : rb;
    } else {
        boost::math::uintmax_t it = 80;
        auto [l, r] = toms748_solve(sres, ra, rb, ga, gb, boost::math::tools::eps_tolerance<double>(40), it);
        rho0 = 0.5 * (l + r);
    }
    return hs_corrector(H, Ta, rho0, h_t, s_t, T_out, rho_out);
}

// Leg 3: ideal-gas departure.  Scale residual Helmholtz by lambda; anchor at
// lambda=0 (ideal gas, no dome) and continue lambda:0->1.  Self-contained.
bool hs_leg_departure(HelmholtzEOSMixtureBackend& H, double h_t, double s_t, double& T_out, double& rho_out) {
    HSGasGuard guard(H);
    const double Rg = H.gas_constant();
    const auto& red = H.get_reducing_state();
    const double Tr = red.T, rhor = red.rhomolar;
    const auto& x = H.get_mole_fractions_ref();
    const double Tmin = H.Tmin(), Tmax = H.Tmax();
    // Properties (and the Jacobian entries the Newton step needs) of the
    // lambda-scaled model alpha = alpha0 + lam*alphar: h, s, the stability
    // marker prho (= dp/drho|_T, must stay > 0), and the partials
    // hT=dh/dT, hr=dh/drho, sT=ds/dT, sr=ds/drho.
    struct LP
    {
        double h, s, prho, hT, hr, sT, sr;
    };
    // Evaluate the above at (T, rho) for continuation parameter lam in [0, 1].
    auto lprops = [&](double T, double rho, double lam) -> LP {
        const double tau = Tr / T, delta = rho / rhor;
        // Each calc_*_deriv_nocache computes the FULL HelmholtzDerivatives set and
        // returns one component, so the per-term calls recompute everything 6-9x.
        // Evaluate the ideal-gas and residual blocks once each instead.
        const HelmholtzDerivatives a0d = H.calc_all_alpha0_derivs_nocache(x, tau, delta, Tr, rhor);
        const HelmholtzDerivatives ard_all = H.residual_helmholtz->all(H, x, tau, delta, false);
        const double a0 = a0d.alphar, a0t = a0d.dalphar_dtau, a0tt = a0d.d2alphar_dtau2;
        const double ar = ard_all.alphar, art = ard_all.dalphar_dtau, ard = ard_all.dalphar_ddelta;
        const double artt = ard_all.d2alphar_dtau2, ardd = ard_all.d2alphar_ddelta2, artd = ard_all.d2alphar_ddelta_dtau;
        const double at = a0t + lam * art, att = a0tt + lam * artt;
        const double Hh = 1.0 + tau * at + delta * lam * ard;
        LP L{};
        L.h = Rg * T * Hh;
        L.s = Rg * (tau * at - (a0 + lam * ar));
        L.prho = Rg * T * (1.0 + 2.0 * delta * lam * ard + delta * delta * lam * ardd);
        const double dHdtau = at + tau * att + delta * lam * artd;
        const double dHddelta = lam * ard + tau * lam * artd + delta * lam * ardd;
        L.hT = Rg * (Hh - tau * dHdtau);
        L.hr = Rg * T * dHddelta / rhor;
        L.sT = -Rg * tau * tau * att / T;
        L.sr = Rg * (tau * lam * artd - 1.0 / delta - lam * ard) / rhor;
        // A non-finite property/derivative (e.g. a degenerate EOS evaluation) must
        // not leak a NaN into the continuation/Newton loop; throw so the cascade
        // defers this leg cleanly (the per-term calc_*_deriv_nocache previously
        // threw on a non-finite alpha0 -- the batched accessor does not).
        if (!std::isfinite(L.h) || !std::isfinite(L.s) || !std::isfinite(L.prho) || !std::isfinite(L.hT) || !std::isfinite(L.hr)
            || !std::isfinite(L.sT) || !std::isfinite(L.sr)) {
            throw ValueError("hs_leg_departure: non-finite lambda-model property/derivative");
        }
        return L;
    };
    // lambda=0 anchor: the ideal-gas limit has no dome, so (h,s)->(T,rho) is a
    // pair of decoupled 1-D solves -- T from h(T) at the reducing density, then
    // rho from s(T0, rho).  This seeds the continuation below.
    double T0 = 0, rho0 = 0;
    {
        auto hig = [&](double T) { return lprops(T, rhor, 0.0).h - h_t; };
        double Ta = 0.3 * Tmin, Tb = 3.0 * Tmax, fa = hig(Ta), fb = hig(Tb);
        if (fa * fb > 0) return false;
        boost::math::uintmax_t it = 80;
        auto [l, r] = toms748_solve(hig, Ta, Tb, fa, fb, boost::math::tools::eps_tolerance<double>(40), it);
        T0 = 0.5 * (l + r);
        auto sig = [&](double rho) { return lprops(T0, rho, 0.0).s - s_t; };
        double ra = 1e-8 * rhor, rb = 50.0 * rhor, ga = sig(ra), gb = sig(rb);
        if (ga * gb > 0) {
            rho0 = (std::abs(ga) <= std::abs(gb)) ? ra : rb;
        } else {
            boost::math::uintmax_t it2 = 80;
            auto [l2, r2] = toms748_solve(sig, ra, rb, ga, gb, boost::math::tools::eps_tolerance<double>(40), it2);
            rho0 = 0.5 * (l2 + r2);
        }
    }
    // Continuation lam: 0 -> 1 in N equal steps, solving in (T, w=log rho) so rho
    // stays positive.  Retry with a finer schedule (N doubling) if a coarse one
    // loses the path; first N that reaches lam=1 wins.
    const double hscale = Rg * Tr, sscale = Rg;
    for (int N = 1; N <= 256; N *= 2) {
        double T = T0, w = std::log(rho0);
        bool failed = false;
        for (int k = 1; k <= N && !failed; ++k) {
            const double lam = static_cast<double>(k) / N;
            bool conv = false;
            double best_norm = 1e300, best_T = T, best_w = w;
            for (int iter = 0; iter < 40; ++iter) {
                const double rho = std::exp(w);
                LP L = lprops(T, rho, lam);
                const double rh = L.h - h_t, rs = L.s - s_t;
                const double norm = std::abs(rh) / hscale + std::abs(rs) / sscale;
                if (norm < best_norm) {
                    best_norm = norm;
                    best_T = T;
                    best_w = w;
                }
                if (norm < 1e-11) {
                    conv = true;
                    break;
                }
                // 2x2 Newton step on (rh, rs); columns scaled by rho since we
                // step in w = log rho (drho = rho*dw).
                const double a11 = L.hT, a12 = rho * L.hr, a21 = L.sT, a22 = rho * L.sr;
                const double det = a11 * a22 - a12 * a21;
                if (!std::isfinite(det) || std::abs(det) < 1e-300) break;
                const double dT = -(a22 * rh - a12 * rs) / det, dw = -(-a21 * rh + a11 * rs) / det;
                // Damp the step: first keep T in range and cap the log-rho move,
                // then back off until the trial state is mechanically stable (prho>0).
                double f = 1.0;
                while ((T + f * dT <= 0 || T + f * dT > 3.0 * Tmax || std::abs(f * dw) > 2.0) && f > 1e-6)
                    f *= 0.5;
                for (int g = 0; g < 8 && f > 1e-3; ++g) {
                    if (lprops(T + f * dT, std::exp(w + f * dw), lam).prho > 0) break;
                    f *= 0.5;
                }
                T += f * dT;
                w += f * dw;
            }
            if (!conv && best_norm < 1e-8) {
                T = best_T;
                w = best_w;
                conv = true;
            }
            if (!conv) failed = true;
        }
        if (!failed) {
            T_out = T;
            rho_out = std::exp(w);
            return true;
        }
    }
    return false;
}

// Accept only a faithful (h,s) reproduction that is in-range and fully stable
// (dp/drho|_T>0 AND cv>0): over that region the (h,s)->(T,rho) map is injective,
// so a stable+matching root IS the unique physical root.
// Tmin_override: if > 0, replaces H.Tmin() as the lower bound (used by the
// melting-line leg to allow sub-triple compressed-liquid states whose EOS IS
// valid below the triple-point temperature).
bool hs_accept(HelmholtzEOSMixtureBackend& H, double T, double rho, double h_t, double s_t, double Tmin_override = -1.0) {
    if (!std::isfinite(T) || !std::isfinite(rho) || rho <= 0) return false;
    const double Rg = H.gas_constant(), Tc = H.T_critical();
    const double Tmin_eff = (Tmin_override > 0) ? Tmin_override : H.Tmin();
    if (T < Tmin_eff * (1 - 1e-6) || T > H.Tmax() * (1 + 1e-6)) return false;
    HSGasGuard guard(H);
    H.update_DmolarT_direct(rho, T);
    if (std::abs(H.hmolar() - h_t) > 1e-6 * Rg * Tc || std::abs(H.smolar() - s_t) > 1e-6 * Rg) return false;
    double cv = -1;
    try {
        cv = H.cvmolar();
    } catch (...) {
        return false;
    }
    if (!std::isfinite(cv) || cv <= 0) return false;
    return H.first_partial_deriv(iP, iDmolar, iT) > 0;
}

// True if (T,rho) lies strictly inside the two-phase dome (a metastable single-
// phase extension).  Such a root reproduces a two-phase (h,s) but is NOT the
// physical answer, so the cascade must reject it and defer to the two-phase solve.
bool hs_inside_dome(HS_SA_t& sa, double T, double rho) {
    if (T >= sa.get_Tcrit_num()) return false;
    return rho > sa.eval_sat(T, 'D', 1) && rho < sa.eval_sat(T, 'D', 0);
}

// Shift a target h/s into a MeltingCaloric's stamped build frame (mirror of
// hs_h_to_cache/hs_s_to_cache, reading the melting caloric's own stamp).
double meltcal_s_to_cache(HelmholtzEOSMixtureBackend& H, const MeltingCaloric& mc, double s_t) {
    const auto stamp = mc.stamp();
    if (!stamp.has_value()) return s_t;
    const double a1 = FlashRoutines::alpha0_offset_total(H).first;
    return s_t - H.gas_constant() * (stamp->first - a1);
}
double meltcal_h_to_cache(HelmholtzEOSMixtureBackend& H, const MeltingCaloric& mc, double h_t) {
    const auto stamp = mc.stamp();
    if (!stamp.has_value()) return h_t;
    const double a2 = FlashRoutines::alpha0_offset_total(H).second;
    return h_t + H.gas_constant() * H.get_reducing_state().T * (stamp->second - a2);
}

// Leg 4: melting-line anchor for the cold compressed-liquid corner (incl. sub-triple T).
bool hs_leg_melting(HelmholtzEOSMixtureBackend& H, const MeltingCaloric& mc, double h_t, double s_t, double& T_out, double& rho_out) {
    HSGasGuard guard(H);  // corrector requires the caller to impose the gas phase (mirror hs_leg_saturation)
    const double s_cache = meltcal_s_to_cache(H, mc, s_t);
    const double h_cache = meltcal_h_to_cache(H, mc, h_t);
    double T0 = 0, rho0 = 0;
    if (!mc.seed_for_hs(s_cache, h_cache, T0, rho0)) return false;
    // Floor at the melting-curve minimum temperature (water folds to ~251 K),
    // minus a small slack, so the corrector can reach the sub-triple region.
    double Tlo = -1.0;
    const double cTmin = mc.curve_Tmin();
    if (std::isfinite(cTmin) && cTmin > 0) Tlo = cTmin * (1.0 - 1e-3);
    return hs_corrector(H, T0, rho0, h_t, s_t, T_out, rho_out, Tlo);
}

// Try the three legs in order; accept the first result that is stable, matching,
// AND outside the dome (a genuine single-phase state).  An in-dome (metastable)
// root means the input is actually two-phase, so reject it -- the dispatcher then
// falls through to the two-phase solve.
bool hs_cascade(HelmholtzEOSMixtureBackend& H, HS_SA_t* sa, double h_t, double s_t, double& T_out, double& rho_out) {
    auto good = [&](double T, double rho) { return hs_accept(H, T, rho, h_t, s_t) && !(sa != nullptr && hs_inside_dome(*sa, T, rho)); };
    // Run one leg; on success commit T_out/rho_out.  A leg whose rootfind wanders
    // out of the EOS domain can throw -- that just defers to the next leg.
    auto run = [&](auto&& leg) -> bool {
        try {
            double T = 0, rho = 0;
            if (leg(T, rho) && good(T, rho)) {
                T_out = T;
                rho_out = rho;
                return true;
            }
        } catch (...) {  // NOLINT(bugprone-empty-catch): a throwing leg defers to the next
        }
        return false;
    };
    if (sa != nullptr && run([&](double& T, double& rho) { return hs_leg_saturation(H, *sa, h_t, s_t, T, rho); })) return true;
    if (run([&](double& T, double& rho) { return hs_leg_isentrope(H, h_t, s_t, T, rho); })) return true;
    if (run([&](double& T, double& rho) { return hs_leg_departure(H, h_t, s_t, T, rho); })) return true;
    // Leg 4: melting-line anchor for the cold compressed-liquid corner. Gated by
    // config + env kill-switch; built lazily and cached.
    // Uses a dedicated run_melt that accepts a lower Tmin bound (the melting-curve
    // minimum temperature, ~251 K for water) instead of H.Tmin() (~273 K), so
    // sub-triple compressed-liquid states -- which are physically valid EOS states
    // but lie below the triple-point temperature -- pass the acceptance gate.
    // All other acceptance criteria (h/s residual, dp/drho > 0, cv > 0, outside
    // dome) are identical to the standard run()/good() path.
    static const bool meltcal_disabled = (std::getenv("COOLPROP_DISABLE_MELTING_CALORIC_HS") != nullptr);
    if (!meltcal_disabled && get_config_bool(ENABLE_MELTING_CALORIC_HS)) {
        if (auto mc = get_melting_caloric_cached(H)) {
            const double melt_Tmin = mc->curve_Tmin();
            auto good_melt = [&](double T, double rho) {
                return hs_accept(H, T, rho, h_t, s_t, melt_Tmin) && !(sa != nullptr && hs_inside_dome(*sa, T, rho));
            };
            auto run_melt = [&](auto&& leg) -> bool {
                try {
                    double T = 0, rho = 0;
                    if (leg(T, rho) && good_melt(T, rho)) {
                        T_out = T;
                        rho_out = rho;
                        return true;
                    }
                } catch (...) {  // NOLINT(bugprone-empty-catch): a throwing leg defers to the next
                }
                return false;
            };
            if (run_melt([&](double& T, double& rho) { return hs_leg_melting(H, *mc, h_t, s_t, T, rho); })) return true;
        }
    }
    return false;
}

// Cheap superancillary two-phase detector (no EOS): scan the Qh==Qs residual on
// the caloric superancillary for a sign change with quality strictly inside (0,1).
// Used only as an OPTIMIZATION -- a true routes the point straight to the
// EOS-exact two-phase solve, skipping the (doomed) single-phase cascade.  It is
// deliberately conservative: a false negative (or a near-critical miss) is safe
// because the cascade rejects the resulting in-dome metastable root and the
// dispatcher's two-phase fallback then runs; a false positive is safe because the
// two-phase solve fails to reproduce (h,s) and we fall through to the cascade.
bool hs_two_phase_likely(HelmholtzEOSMixtureBackend& H, HS_SA_t& sa, double h_t, double s_t) {
    H.ensure_caloric_superancillaries();
    if (!sa.has_variable('H') || !sa.has_variable('S')) return false;
    const double hc = hs_h_to_cache(H, sa, h_t), sc = hs_s_to_cache(H, sa, s_t);
    const double Tmin = sa.get_Tmin(), Tc = sa.get_Tcrit_num();
    auto Qh_minus_Qs = [&](double T) -> double {
        const double sL = sa.eval_sat(T, 'S', 0), sV = sa.eval_sat(T, 'S', 1);
        const double hL = sa.eval_sat(T, 'H', 0), hV = sa.eval_sat(T, 'H', 1);
        return (hc - hL) / (hV - hL) - (sc - sL) / (sV - sL);
    };
    // Coarse scan for a sign change: the endpoint near Tc is ill-conditioned
    // (rhoL->rhoV, the quality ratios degenerate), so scan rather than bracket the
    // ends -- this is what a naive two-endpoint test missed at low quality.
    const int M = 40;
    const double Tlo = Tmin + 1e-3, Thi = Tc - std::max(0.5, 1e-3 * Tc);
    double Tprev = Tlo, fprev = Qh_minus_Qs(Tlo);
    for (int i = 1; i <= M; ++i) {
        const double T = Tlo + (Thi - Tlo) * i / M;
        const double f = Qh_minus_Qs(T);
        if (std::isfinite(fprev) && std::isfinite(f) && fprev * f <= 0) {
            boost::math::uintmax_t it = 60;
            auto [l, r] = toms748_solve(Qh_minus_Qs, Tprev, T, fprev, f, boost::math::tools::eps_tolerance<double>(40), it);
            const double Tsol = 0.5 * (l + r);
            const double sL = sa.eval_sat(Tsol, 'S', 0), sV = sa.eval_sat(Tsol, 'S', 1);
            const double Q = (sc - sL) / (sV - sL);
            if (Q > 1e-6 && Q < 1.0 - 1e-6) return true;  // strictly interior => two-phase
        }
        Tprev = T;
        fprev = f;
    }
    return false;
}

}  // namespace

bool FlashRoutines::hs_corrector_probe(HelmholtzEOSMixtureBackend& H, double T0, double rho0, double h_t, double s_t, double& T_out, double& rho_out,
                                       double Tlo_override) {
    HSGasGuard guard(H);  // hs_corrector requires caller to impose single-phase
    return hs_corrector(H, T0, rho0, h_t, s_t, T_out, rho_out, Tlo_override);
}

void FlashRoutines::HS_flash(HelmholtzEOSMixtureBackend& HEOS) {
    // ===================== superancillary "happy path" =====================
    // For a pure fluid with a built superancillary, robustly classify two-phase
    // vs single-phase from (h,s) and solve each fast: two-phase by the Qh==Qs
    // saturation solve (mirror of the D+X happy path), single-phase by the
    // cascade above.  Any failure falls through to the legacy "sad path" below.
    // Kill-switch COOLPROP_DISABLE_SUPERANC_HS forces the legacy path.
    {
        static const bool hs_disabled = (std::getenv("COOLPROP_DISABLE_SUPERANC_HS") != nullptr);
        if (!hs_disabled && get_config_bool(ENABLE_SUPERANCILLARIES) && HEOS.is_pure()) {
            std::shared_ptr<EquationOfState::SuperAncillary_t> sa_ptr = HEOS.get_superanc();
            if (sa_ptr) {
                const double h_t = HEOS._hmolar, s_t = HEOS._smolar;
                auto restore_inputs = [&]() {
                    HEOS._hmolar = h_t;
                    HEOS._smolar = s_t;
                    HEOS.unspecify_phase();
                };
                auto reproduces = [&]() -> bool {
                    const double hh = HEOS.hmolar(), ss = HEOS.smolar();
                    return std::isfinite(hh) && std::isfinite(ss) && std::abs(hh - h_t) <= 1e-6 * std::abs(h_t) + 1e-3
                           && std::abs(ss - s_t) <= 1e-6 * std::abs(s_t) + 1e-5;
                };
                try {
                    EquationOfState::SuperAncillary_t& sa = *sa_ptr;
                    bool tried_2ph = false;
                    auto try_twophase = [&]() -> bool {
                        HS_flash_twophaseOptions opt;
                        HS_flash_twophase(HEOS, h_t, s_t, opt);
                        tried_2ph = true;
                        return reproduces();
                    };
                    // (0) Fast two-phase screen (no EOS): if (h,s) is clearly inside the
                    // dome, solve two-phase directly and skip the doomed single-phase
                    // cascade (~50-100 wasted evals for two-phase inputs).
                    if (hs_two_phase_likely(HEOS, sa, h_t, s_t)) {
                        if (try_twophase()) return;
                    }
                    // (1) Single-phase cascade (the common case, ~16 evals).  It accepts
                    // only a stable root OUTSIDE the dome, so for a two-phase (h,s) every
                    // leg lands on an in-dome metastable root and is rejected -> the
                    // cascade fails and we fall through to (2).  This dome-rejection is
                    // the robust safety net behind the conservative screen in (0).
                    double T = 0, rho = 0;
                    if (hs_cascade(HEOS, sa_ptr.get(), h_t, s_t, T, rho)) {
                        HEOS.update_DmolarT_direct(rho, T);
                        HEOS._Q = 10000;
                        HEOS._p = HEOS.calc_pressure_nocache(HEOS.T(), HEOS.rhomolar());
                        HEOS.unspecify_phase();
                        HEOS.recalculate_singlephase_phase();
                        if (reproduces()) return;
                    }
                    // (2) Two-phase fallback (if the screen did not already try it): the
                    // cascade found no stable single-phase root, so solve the EOS-exact
                    // Qh==Qs problem and accept only if it reproduces (h,s).
                    if (!tried_2ph) {
                        if (try_twophase()) return;
                    }
                    // Every fall-through to the legacy path below MUST see the original
                    // inputs: the cascade / two-phase attempts mutate HEOS._hmolar,
                    // _smolar (via update_DmolarT_direct -> clear()), so restore them
                    // unconditionally here.  (Skipping this on the screen-fired +
                    // cascade-failed path corrupted the legacy solve's target.)
                    restore_inputs();
                } catch (const std::exception& e) {
                    if (get_debug_level() > 0) {
                        std::cout << "HS_flash superancillary path failed (" << e.what() << "); using legacy path\n";
                    }
                    restore_inputs();
                }
            }
        }
    }
    // ===================== legacy "sad path" =====================
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
        double call(double T) override {
            HEOS.update(SmolarT_INPUTS, smolar, T);
            double r = HEOS.hmolar() - hmolar;
            return r;
        }
    } resid(HEOS, hmolar, smolar);
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
    shared_ptr<HelmholtzEOSBackend> HEOS = std::make_shared<HelmholtzEOSBackend>("R134a");
    double Tc = HEOS->T_critical();
    HEOS->update(DmassT_INPUTS, 1.1, 1.5 * Tc);
    CHECK_THROWS(HEOS->update(DmassP_INPUTS, 2, 5 * HEOS->p()));
}

TEST_CASE("Stability testing", "[stability]") {
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS =
      std::make_shared<HelmholtzEOSMixtureBackend>(strsplit("n-Propane&n-Butane&n-Pentane&n-Hexane", '&'));
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
        // Unit step over an integer range — T values are bit-exact doubles, no accumulation.
        for (double T = TL - 1; T >= 100; T -= 1) {  // NOLINT(cert-flp30-c)
            stability_tester.set_TP(T, 101325);
            CAPTURE(T);
            CHECK_NOTHROW(stability_tester.is_stable());
        }
    }
    SECTION("Vapor (feed is stable)") {
        StabilityRoutines::StabilityEvaluationClass stability_tester(*HEOS);
        // Unit step over an integer range — see liquid section above.
        for (double T = TV + 1; T <= 500; T += 1) {  // NOLINT(cert-flp30-c)
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
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS = std::make_shared<HelmholtzEOSMixtureBackend>(strsplit("Methane&H2S", '&'));

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
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS = std::make_shared<HelmholtzEOSMixtureBackend>(strsplit("Nitrogen&Ethane", '&'));
    std::vector<double> zz = linspace(0.001, 0.999, 21);
    int failure_count = 0;
    for (double z0 : zz) {
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
    shared_ptr<PengRobinsonBackend> HEOS = std::make_shared<PengRobinsonBackend>(strsplit("Nitrogen&Ethane", '&'));
    HEOS->set_binary_interaction_double(0, 1, "kij", 0.0407);  // Ramırez-Jimenez et al.
    std::vector<double> zz = linspace(0.001, 0.999, 21);
    for (double z0 : zz) {
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
