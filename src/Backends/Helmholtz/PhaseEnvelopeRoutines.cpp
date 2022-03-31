#ifndef PHASEENVELOPE_H
#define PHASEENVELOPE_H

#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "PhaseEnvelopeRoutines.h"
#include "PhaseEnvelope.h"
#include "CoolPropTools.h"
#include "Configuration.h"
#include "CPnumerics.h"

namespace CoolProp {

void PhaseEnvelopeRoutines::build(HelmholtzEOSMixtureBackend& HEOS, const std::string& level) {
    if (HEOS.get_mole_fractions_ref().empty()) {
        throw ValueError("Mole fractions have not been set yet.");
    }
    bool debug = get_debug_level() > 0 || false;
    if (HEOS.get_mole_fractions_ref().size() == 1) {
        // It's a pure fluid
        PhaseEnvelopeData& env = HEOS.PhaseEnvelope;
        env.resize(HEOS.mole_fractions.size());

        // Breakpoints in the phase envelope
        std::vector<CoolPropDbl> Tbp, Qbp;
        std::vector<std::size_t> Nbp;
        // Triple point vapor up to Tmax_sat
        Tbp.push_back(HEOS.Ttriple());
        Qbp.push_back(1.0);
        Nbp.push_back(40);

        if (HEOS.is_pure()) {
            // Up to critical point, back to triple point on the liquid side
            Tbp.push_back(HEOS.T_critical() - 1e-3);
            Qbp.push_back(0.0);
            Tbp.push_back(HEOS.Ttriple());
            Nbp.push_back(40);
        } else {
            SimpleState max_sat_T = HEOS.get_state("max_sat_T"), max_sat_p = HEOS.get_state("max_sat_p"), crit = HEOS.get_state("critical");
            if (max_sat_T.rhomolar < crit.rhomolar && max_sat_T.rhomolar < max_sat_p.rhomolar) {
                Tbp.push_back(HEOS.calc_Tmax_sat());
                if (max_sat_p.rhomolar < crit.rhomolar) {
                    // psat_max density less than critical density
                    Qbp.push_back(1.0);
                    Qbp.push_back(1.0);
                    Tbp.push_back(max_sat_p.T);
                    Tbp.push_back(crit.T);
                } else {
                    // Vapor line density less than critical density
                    Qbp.push_back(1.0);
                    Qbp.push_back(0.0);
                    Nbp.push_back(10);
                    Tbp.push_back(crit.T);
                    Tbp.push_back(max_sat_p.T);
                }
                Nbp.push_back(10);
                Nbp.push_back(10);
                Qbp.push_back(0.0);
                Nbp.push_back(40);
                Tbp.push_back(HEOS.Ttriple());
            } else {
                throw ValueError(format(""));
            }
        }

        for (std::size_t i = 0; i < Tbp.size() - 1; ++i) {
            CoolPropDbl Tmin = Tbp[i], Tmax = Tbp[i + 1];
            std::size_t N = Nbp[i];
            for (CoolPropDbl T = Tmin; is_in_closed_range(Tmin, Tmax, T); T += (Tmax - Tmin) / (N - 1)) {
                try {
                    HEOS.update(QT_INPUTS, Qbp[i], T);
                } catch (...) {
                    continue;
                }
                if (Qbp[i] > 0.5) {
                    env.store_variables(HEOS.T(), HEOS.p(), HEOS.saturated_liquid_keyed_output(iDmolar), HEOS.saturated_vapor_keyed_output(iDmolar),
                                        HEOS.saturated_liquid_keyed_output(iHmolar), HEOS.saturated_vapor_keyed_output(iHmolar),
                                        HEOS.saturated_liquid_keyed_output(iSmolar), HEOS.saturated_vapor_keyed_output(iSmolar),
                                        std::vector<CoolPropDbl>(1, 1.0), std::vector<CoolPropDbl>(1, 1.0));
                } else {
                    env.store_variables(HEOS.T(), HEOS.p(), HEOS.saturated_vapor_keyed_output(iDmolar), HEOS.saturated_liquid_keyed_output(iDmolar),
                                        HEOS.saturated_vapor_keyed_output(iHmolar), HEOS.saturated_liquid_keyed_output(iHmolar),
                                        HEOS.saturated_vapor_keyed_output(iSmolar), HEOS.saturated_liquid_keyed_output(iSmolar),
                                        std::vector<CoolPropDbl>(1, 1.0), std::vector<CoolPropDbl>(1, 1.0));
                }
            }
        }
    } else {
        // It's a mixture
        // --------------

        // First we try to generate all the critical points.  This
        // is very useful
        std::vector<CriticalState> critpts;
        //        try{
        //             critpts = HEOS.all_critical_points();
        //            //throw CoolProp::ValueError("critical points disabled");
        //        }
        //        catch(std::exception &e)
        //        {
        //            if (debug){ std::cout << e.what() << std::endl; }
        //        };

        std::size_t failure_count = 0;
        // Set some input options
        SaturationSolvers::mixture_VLE_IO io;
        io.sstype = SaturationSolvers::imposed_p;
        io.Nstep_max = 20;

        // Set the pressure to a low pressure
        HEOS._p = get_config_double(PHASE_ENVELOPE_STARTING_PRESSURE_PA);  //[Pa]
        HEOS._Q = 1;

        // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
        CoolPropDbl Tguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

        // Use Wilson iteration to obtain updated guess for temperature
        Tguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);

        // Actually call the successive substitution solver
        io.beta = 1;
        SaturationSolvers::successive_substitution(HEOS, HEOS._Q, Tguess, HEOS._p, HEOS.mole_fractions, HEOS.K, io);

        // Use the residual function based on x_i, T and rho' as independent variables.  rho'' is specified
        SaturationSolvers::newton_raphson_saturation NR;
        SaturationSolvers::newton_raphson_saturation_options IO;

        IO.bubble_point = false;  // Do a "dewpoint" calculation all the way around
        IO.x = io.x;
        IO.y = HEOS.mole_fractions;
        IO.rhomolar_liq = io.rhomolar_liq;
        IO.rhomolar_vap = io.rhomolar_vap;
        IO.T = io.T;
        IO.p = io.p;
        IO.Nstep_max = 30;

        /*
        IO.p = 1e5;
        IO.rhomolar_liq = 17257.17130;
        IO.rhomolar_vap = 56.80022884;
        IO.T = 219.5200523;
        IO.x[0] = 0.6689704673;
        IO.x[1] = 0.3310295327;
        */

        //IO.rhomolar_liq *= 1.2;

        IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;

        NR.call(HEOS, IO.y, IO.x, IO);

        // Switch to density imposed
        IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::RHOV_IMPOSED;

        bool dont_extrapolate = false;

        PhaseEnvelopeData& env = HEOS.PhaseEnvelope;
        env.resize(HEOS.mole_fractions.size());

        std::size_t iter = 0,  //< The iteration counter
          iter0 = 0;           //< A reference point for the counter, can be increased to go back to linear interpolation
        CoolPropDbl factor = 1.05;

        for (;;) {
        top_of_loop:;  // A goto label so that nested loops can break out to the top of this loop

            if (failure_count > 5) {
                // Stop since we are stuck at a bad point
                //throw SolutionError("stuck");
                return;
            }

            if (iter - iter0 > 0) {
                IO.rhomolar_vap *= factor;
            }
            if (dont_extrapolate) {
                // Reset the step to a reasonably small size
                factor = 1.0001;
            } else if (iter - iter0 == 2) {
                IO.T = LinearInterp(env.rhomolar_vap, env.T, iter - 2, iter - 1, IO.rhomolar_vap);
                IO.rhomolar_liq = LinearInterp(env.rhomolar_vap, env.rhomolar_liq, iter - 2, iter - 1, IO.rhomolar_vap);
                for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
                {
                    IO.x[i] = LinearInterp(env.rhomolar_vap, env.x[i], iter - 2, iter - 1, IO.rhomolar_vap);
                }
            } else if (iter - iter0 == 3) {
                IO.T = QuadInterp(env.rhomolar_vap, env.T, iter - 3, iter - 2, iter - 1, IO.rhomolar_vap);
                IO.rhomolar_liq = QuadInterp(env.rhomolar_vap, env.rhomolar_liq, iter - 3, iter - 2, iter - 1, IO.rhomolar_vap);
                for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
                {
                    IO.x[i] = QuadInterp(env.rhomolar_vap, env.x[i], iter - 3, iter - 2, iter - 1, IO.rhomolar_vap);
                }
            } else if (iter - iter0 > 3) {
                // Use the spline interpolation class of Devin Lane: http://shiftedbits.org/2011/01/30/cubic-spline-interpolation/
                Spline<double, double> spl_T(env.rhomolar_vap, env.T);
                IO.T = spl_T.interpolate(IO.rhomolar_vap);
                Spline<double, double> spl_rho(env.rhomolar_vap, env.rhomolar_liq);
                IO.rhomolar_liq = spl_rho.interpolate(IO.rhomolar_vap);

                // Check if there is a large deviation from linear interpolation - this suggests a step size that is so large that a minima or maxima of the interpolation function is crossed
                CoolPropDbl T_linear = LinearInterp(env.rhomolar_vap, env.T, iter - 2, iter - 1, IO.rhomolar_vap);
                if (std::abs((T_linear - IO.T) / IO.T) > 0.1) {
                    // Try again, but with a smaller step
                    IO.rhomolar_vap /= factor;
                    factor = 1 + (factor - 1) / 2;
                    failure_count++;
                    continue;
                }
                for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
                {
                    // Use the spline interpolation class of Devin Lane: http://shiftedbits.org/2011/01/30/cubic-spline-interpolation/
                    Spline<double, double> spl(env.rhomolar_vap, env.x[i]);
                    IO.x[i] = spl.interpolate(IO.rhomolar_vap);

                    if (IO.x[i] < 0 || IO.x[i] > 1) {
                        // Try again, but with a smaller step
                        IO.rhomolar_vap /= factor;
                        factor = 1 + (factor - 1) / 2;
                        failure_count++;
                        goto top_of_loop;
                    }
                }
            }

            // The last mole fraction is sum of N-1 first elements
            IO.x[IO.x.size() - 1] = 1 - std::accumulate(IO.x.begin(), IO.x.end() - 1, 0.0);

            // Uncomment to check guess values for Newton-Raphson
            //std::cout << "\t\tdv " << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " T " << IO.T << " x " << vec_to_string(IO.x, "%0.10Lg") << std::endl;

            // Dewpoint calculation, liquid (x) is incipient phase
            try {
                NR.call(HEOS, IO.y, IO.x, IO);
                if (!ValidNumber(IO.rhomolar_liq) || !ValidNumber(IO.p) || !ValidNumber(IO.T)) {
                    throw ValueError("Invalid number");
                }
                // Reject trivial solution
                if (std::abs(IO.rhomolar_liq - IO.rhomolar_vap) < 1e-3) {
                    throw ValueError("Trivial solution");
                }
                // Reject negative presssure
                if (IO.p < 0) {
                    throw ValueError("negative pressure");
                }
                // Reject steps with enormous steps in temperature
                if (!env.T.empty() && std::abs(env.T[env.T.size() - 1] - IO.T) > 100) {
                    throw ValueError("Change in temperature too large");
                }
            } catch (std::exception& e) {
                if (debug) {
                    std::cout << e.what() << std::endl;
                }
                //std::cout << IO.T << " " << IO.p << std::endl;
                // Try again, but with a smaller step
                IO.rhomolar_vap /= factor;
                if (iter < 4) {
                    throw ValueError(format("Unable to calculate at least 4 points in phase envelope; quitting"));
                }
                IO.rhomolar_liq = QuadInterp(env.rhomolar_vap, env.rhomolar_liq, iter - 3, iter - 2, iter - 1, IO.rhomolar_vap);
                factor = 1 + (factor - 1) / 2;
                failure_count++;
                continue;
            }

            if (debug) {
                std::cout << "dv " << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " T " << IO.T << " p " << IO.p << " hl " << IO.hmolar_liq
                          << " hv " << IO.hmolar_vap << " sl " << IO.smolar_liq << " sv " << IO.smolar_vap << " x " << vec_to_string(IO.x, "%0.10Lg")
                          << " Ns " << IO.Nsteps << " factor " << factor << std::endl;
            }
            env.store_variables(IO.T, IO.p, IO.rhomolar_liq, IO.rhomolar_vap, IO.hmolar_liq, IO.hmolar_vap, IO.smolar_liq, IO.smolar_vap, IO.x, IO.y);

            iter++;

            //            CoolPropDbl abs_rho_difference = std::abs((IO.rhomolar_liq - IO.rhomolar_vap)/IO.rhomolar_liq);

            //            bool next_crosses_crit = false;
            //            if (it_critpts != critpts.end() ){
            //                // Density at the next critical point
            //                double rhoc = (*it_critpts).rhomolar;
            //                // Next vapor density that will be used
            //                double rho_next = IO.rhomolar_vap*factor;
            //                // If the signs of the differences are different, you have crossed
            //                // the critical point density and have a phase inversion
            //                // on your hands
            //                next_crosses_crit = ((IO.rhomolar_vap-rhoc)*(rho_next-rhoc) < 0);
            //            }

            //            // Critical point jump
            //            if (next_crosses_crit || (abs_rho_difference < 0.01 && IO.rhomolar_liq  > IO.rhomolar_vap)){
            //                //std::cout << "dv" << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
            //                CoolPropDbl rhoc_approx = 0.5*IO.rhomolar_liq + 0.5*IO.rhomolar_vap;
            //                if (it_critpts != critpts.end() ){
            //                    // We actually know what the critical point is to numerical precision
            //                    rhoc_approx = (*it_critpts).rhomolar;
            //                }
            //                CoolPropDbl rho_vap_new = 1.05*rhoc_approx;
            //                // Linearly interpolate to get new guess for T
            //                IO.T = LinearInterp(env.rhomolar_vap,env.T,iter-2,iter-1,rho_vap_new);
            //                IO.rhomolar_liq = LinearInterp(env.rhomolar_vap, env.rhomolar_liq, iter-2, iter-1, rho_vap_new);
            //                for (std::size_t i = 0; i < IO.x.size()-1; ++i){
            //                    IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], iter-4, iter-3, iter-2, iter-1, rho_vap_new);
            //                }
            //                IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
            //                factor = rho_vap_new/IO.rhomolar_vap;
            //                dont_extrapolate = true; // So that we use the mole fractions we calculated here instead of the extrapolated values
            //                if (debug) std::cout << "[CRIT jump] new values: dv " << rho_vap_new << " dl " << IO.rhomolar_liq << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
            //                iter0 = iter - 1; // Back to linear interpolation again
            //                continue;
            //            }

            dont_extrapolate = false;
            if (iter < 5) {
                continue;
            }
            if (IO.Nsteps > 10) {
                factor = 1 + (factor - 1) / 10;
            } else if (IO.Nsteps > 5) {
                factor = 1 + (factor - 1) / 3;
            } else if (IO.Nsteps <= 4) {
                factor = 1 + (factor - 1) * 2;
            }
            // Min step is 1.01
            factor = std::max(factor, static_cast<CoolPropDbl>(1.01));
            // As we approach the critical point, control step size
            if (std::abs(IO.rhomolar_liq / IO.rhomolar_vap - 1) < 4) {
                // Max step is 1.1
                factor = std::min(factor, static_cast<CoolPropDbl>(1.1));
            }

            // Stop if the pressure is below the starting pressure
            // or if the composition of one of the phases becomes almost pure
            CoolPropDbl max_fraction = *std::max_element(IO.x.begin(), IO.x.end());
            if (iter > 4 && (IO.p < env.p[0] || std::abs(1.0 - max_fraction) < 1e-9)) {
                env.built = true;
                if (debug) {
                    std::cout << format("envelope built.\n");
                    std::cout << format("closest fraction to 1.0: distance %g\n", 1 - max_fraction);
                }

                // Now we refine the phase envelope to add some points in places that are still pretty rough
                refine(HEOS, level);

                return;
            }

            // Reset the failure counter
            failure_count = 0;
        }
    }
}

void PhaseEnvelopeRoutines::refine(HelmholtzEOSMixtureBackend& HEOS, const std::string& level) {
    bool debug = (get_debug_level() > 0 || false);
    PhaseEnvelopeData& env = HEOS.PhaseEnvelope;
    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::RHOV_IMPOSED;
    IO.bubble_point = false;
    IO.y = HEOS.get_mole_fractions();

    double acceptable_pdiff = 0.5;
    double acceptable_rhodiff = 0.25;
    int N = 5;  // Number of steps of refining
    if (level == "veryfine") {
        acceptable_pdiff = 0.1;
        acceptable_rhodiff = 0.1;
    }
    if (level == "none") {
        return;
    }
    std::size_t i = 0;
    do {

        // Don't do anything if change in density and pressure is small enough
        if ((std::abs(env.rhomolar_vap[i] / env.rhomolar_vap[i + 1] - 1) < acceptable_rhodiff)
            && (std::abs(env.p[i] / env.p[i + 1] - 1) < acceptable_pdiff)) {
            i++;
            continue;
        }

        // Ok, now we are going to do some more refining in this step

        // Vapor densities for this step, vapor density monotonically increasing
        const double rhomolar_vap_start = env.rhomolar_vap[i], rhomolar_vap_end = env.rhomolar_vap[i + 1];

        double factor = pow(rhomolar_vap_end / rhomolar_vap_start, 1.0 / N);

        int failure_count = 0;
        for (double rhomolar_vap = rhomolar_vap_start * factor; rhomolar_vap < rhomolar_vap_end; rhomolar_vap *= factor) {
            IO.rhomolar_vap = rhomolar_vap;
            IO.x.resize(IO.y.size());
            if (i < env.T.size() - 3) {
                IO.T = CubicInterp(env.rhomolar_vap, env.T, i, i + 1, i + 2, i + 3, IO.rhomolar_vap);
                IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, i, i + 1, i + 2, i + 3, IO.rhomolar_vap);
                for (std::size_t j = 0; j < IO.x.size() - 1; ++j) {  // First N-1 elements
                    IO.x[j] = CubicInterp(env.rhomolar_vap, env.x[j], i, i + 1, i + 2, i + 3, IO.rhomolar_vap);
                }
            } else {
                IO.T = CubicInterp(env.rhomolar_vap, env.T, i, i - 1, i - 2, i - 3, IO.rhomolar_vap);
                IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, i, i - 1, i - 2, i - 3, IO.rhomolar_vap);
                for (std::size_t j = 0; j < IO.x.size() - 1; ++j) {  // First N-1 elements
                    IO.x[j] = CubicInterp(env.rhomolar_vap, env.x[j], i, i - 1, i - 2, i - 3, IO.rhomolar_vap);
                }
            }
            IO.x[IO.x.size() - 1] = 1 - std::accumulate(IO.x.begin(), IO.x.end() - 1, 0.0);
            try {
                NR.call(HEOS, IO.y, IO.x, IO);
                if (!ValidNumber(IO.rhomolar_liq) || !ValidNumber(IO.p)) {
                    throw ValueError("invalid numbers");
                }
                env.insert_variables(IO.T, IO.p, IO.rhomolar_liq, IO.rhomolar_vap, IO.hmolar_liq, IO.hmolar_vap, IO.smolar_liq, IO.smolar_vap, IO.x,
                                     IO.y, i + 1);
                if (debug) {
                    std::cout << "dv " << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " T " << IO.T << " p " << IO.p << " hl " << IO.hmolar_liq
                              << " hv " << IO.hmolar_vap << " sl " << IO.smolar_liq << " sv " << IO.smolar_vap << " x "
                              << vec_to_string(IO.x, "%0.10Lg") << " Ns " << IO.Nsteps << std::endl;
                }
            } catch (...) {
                failure_count++;
                continue;
            }
            i++;
        }
        // If we had a failure, we don't want to get stuck on this value of i,
        // so we bump up one and keep moving
        if (failure_count > 0) {
            i++;
        }
    } while (i < env.T.size() - 1);
}
double PhaseEnvelopeRoutines::evaluate(const PhaseEnvelopeData& env, parameters output, parameters iInput1, double value1, std::size_t& i) {
    int _i = static_cast<int>(i);
    std::vector<double> const *x, *y;

    switch (output) {
        case iT:
            y = &(env.T);
            break;
        case iP:
            y = &(env.p);
            break;
        case iDmolar:
            y = &(env.rhomolar_vap);
            break;
        case iHmolar:
            y = &(env.hmolar_vap);
            break;
        case iSmolar:
            y = &(env.smolar_vap);
            break;
        case iCpmolar:
            y = &(env.cpmolar_vap);
            break;
        case iCvmolar:
            y = &(env.cvmolar_vap);
            break;
        case iviscosity:
            y = &(env.viscosity_vap);
            break;
        case iconductivity:
            y = &(env.conductivity_vap);
            break;
        case ispeed_sound:
            y = &(env.speed_sound_vap);
            break;
        default:
            throw ValueError("Pointer to vector y is unset in is_inside");
    }

    double inval = value1;
    switch (iInput1) {
        case iT:
            x = &(env.T);
            break;
        case iP:
            x = &(env.lnp);
            inval = log(value1);
            break;
        case iDmolar:
            x = &(env.rhomolar_vap);
            break;
        case iHmolar:
            x = &(env.hmolar_vap);
            break;
        case iSmolar:
            x = &(env.smolar_vap);
            break;
        default:
            throw ValueError("Pointer to vector x is unset in is_inside");
    }
    if (_i + 2 >= static_cast<int>(y->size())) {
        _i--;
    }
    if (_i + 1 >= static_cast<int>(y->size())) {
        _i--;
    }
    if (_i - 1 < 0) {
        _i++;
    }

    double outval = CubicInterp(*x, *y, _i - 1, _i, _i + 1, _i + 2, inval);
    i = static_cast<std::size_t>(_i);
    return outval;
}
void PhaseEnvelopeRoutines::finalize(HelmholtzEOSMixtureBackend& HEOS) {
    // No finalization for pure or pseudo-pure fluids
    if (HEOS.get_mole_fractions_ref().size() == 1) {
        return;
    }

    enum maxima_points
    {
        PMAX_SAT = 0,
        TMAX_SAT = 1
    };
    std::size_t imax;  // Index of the maximal temperature or pressure

    PhaseEnvelopeData& env = HEOS.PhaseEnvelope;

    // Find the index of the point with the highest temperature
    std::size_t iTmax = std::distance(env.T.begin(), std::max_element(env.T.begin(), env.T.end()));

    // Find the index of the point with the highest pressure
    std::size_t ipmax = std::distance(env.p.begin(), std::max_element(env.p.begin(), env.p.end()));

    // Determine if the phase envelope corresponds to a Type I mixture
    // For now we consider a mixture to be Type I if the pressure at the
    // end of the envelope is lower than max pressure pressure
    env.TypeI = env.p[env.p.size() - 1] < env.p[ipmax];

    // Approximate solutions for the maxima of the phase envelope
    // See method in Gernert.  We use our spline class to find the coefficients
    if (env.TypeI) {
        for (int imaxima = 0; imaxima <= 1; ++imaxima) {
            maxima_points maxima;
            if (imaxima == PMAX_SAT) {
                maxima = PMAX_SAT;
            } else if (imaxima == TMAX_SAT) {
                maxima = TMAX_SAT;
            } else {
                throw ValueError("I don't understand your maxima index");
            }

            // Spline using the points around it
            SplineClass spline;
            if (maxima == TMAX_SAT) {
                imax = iTmax;
                if (iTmax > env.T.size() - 3) {
                    iTmax -= 2;
                }
                spline.add_4value_constraints(env.rhomolar_vap[iTmax - 1], env.rhomolar_vap[iTmax], env.rhomolar_vap[iTmax + 1],
                                              env.rhomolar_vap[iTmax + 2], env.T[iTmax - 1], env.T[iTmax], env.T[iTmax + 1], env.T[iTmax + 2]);
            } else {
                imax = ipmax;
                if (ipmax > env.p.size() - 3) {
                    ipmax -= 2;
                }
                spline.add_4value_constraints(env.rhomolar_vap[ipmax - 1], env.rhomolar_vap[ipmax], env.rhomolar_vap[ipmax + 1],
                                              env.rhomolar_vap[ipmax + 2], env.p[ipmax - 1], env.p[ipmax], env.p[ipmax + 1], env.p[ipmax + 2]);
            }
            spline.build();  // y = a*rho^3 + b*rho^2 + c*rho + d

            // Take derivative
            // dy/drho = 3*a*rho^2 + 2*b*rho + c
            // Solve quadratic for derivative to find rho
            int Nsoln = 0;
            double rho0 = _HUGE, rho1 = _HUGE, rho2 = _HUGE;
            solve_cubic(0, 3 * spline.a, 2 * spline.b, spline.c, Nsoln, rho0, rho1, rho2);

            SaturationSolvers::newton_raphson_saturation_options IO;
            IO.rhomolar_vap = _HUGE;
            // Find the correct solution
            if (Nsoln == 1) {
                IO.rhomolar_vap = rho0;
            } else if (Nsoln == 2) {
                if (is_in_closed_range(env.rhomolar_vap[imax - 1], env.rhomolar_vap[imax + 1], rho0)) {
                    IO.rhomolar_vap = rho0;
                }
                if (is_in_closed_range(env.rhomolar_vap[imax - 1], env.rhomolar_vap[imax + 1], rho1)) {
                    IO.rhomolar_vap = rho1;
                }
            } else {
                throw ValueError("More than 2 solutions found");
            }

            class solver_resid : public FuncWrapper1D
            {
               public:
                std::size_t imax;
                maxima_points maxima;
                HelmholtzEOSMixtureBackend* HEOS;
                SaturationSolvers::newton_raphson_saturation NR;
                SaturationSolvers::newton_raphson_saturation_options IO;
                solver_resid(HelmholtzEOSMixtureBackend& HEOS, std::size_t imax, maxima_points maxima) {
                    this->HEOS = &HEOS, this->imax = imax;
                    this->maxima = maxima;
                };
                double call(double rhomolar_vap) {
                    PhaseEnvelopeData& env = HEOS->PhaseEnvelope;
                    IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::RHOV_IMPOSED;
                    IO.bubble_point = false;
                    IO.rhomolar_vap = rhomolar_vap;
                    IO.y = HEOS->get_mole_fractions();
                    IO.x = IO.y;  // Just to give it good size
                    if (imax >= env.T.size() - 2) {
                        imax -= 2;
                    }
                    IO.T = CubicInterp(env.rhomolar_vap, env.T, imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
                    IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
                    for (std::size_t i = 0; i < IO.x.size() - 1; ++i)  // First N-1 elements
                    {
                        IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], imax - 1, imax, imax + 1, imax + 2, IO.rhomolar_vap);
                    }
                    IO.x[IO.x.size() - 1] = 1 - std::accumulate(IO.x.begin(), IO.x.end() - 1, 0.0);
                    NR.call(*HEOS, IO.y, IO.x, IO);
                    if (maxima == TMAX_SAT) {
                        return NR.dTsat_dPsat;
                    } else {
                        return NR.dPsat_dTsat;
                    }
                };
            };

            solver_resid resid(HEOS, imax, maxima);
            try {
                double rho = Brent(resid, IO.rhomolar_vap * 0.95, IO.rhomolar_vap * 1.05, DBL_EPSILON, 1e-12, 100);

                // If maxima point is greater than density at point from the phase envelope, increase index by 1 so that the
                // insertion will happen *after* the point in the envelope since density is monotonically increasing.
                if (rho > env.rhomolar_vap[imax]) {
                    imax++;
                }

                env.insert_variables(resid.IO.T, resid.IO.p, resid.IO.rhomolar_liq, resid.IO.rhomolar_vap, resid.IO.hmolar_liq, resid.IO.hmolar_vap,
                                     resid.IO.smolar_liq, resid.IO.smolar_vap, resid.IO.x, resid.IO.y, imax);
            } catch (...) {
                // Don't do the insertion
            }
        }
    }

    // Find the index of the point with the highest temperature
    env.iTsat_max = std::distance(env.T.begin(), std::max_element(env.T.begin(), env.T.end()));

    // Find the index of the point with the highest pressure
    env.ipsat_max = std::distance(env.p.begin(), std::max_element(env.p.begin(), env.p.end()));
}

std::vector<std::pair<std::size_t, std::size_t>> PhaseEnvelopeRoutines::find_intersections(const PhaseEnvelopeData& env, parameters iInput,
                                                                                           double value) {
    std::vector<std::pair<std::size_t, std::size_t>> intersections;

    for (std::size_t i = 0; i < env.p.size() - 1; ++i) {
        bool matched = false;
        switch (iInput) {
            case iP:
                if (is_in_closed_range(env.p[i], env.p[i + 1], value)) {
                    matched = true;
                }
                break;
            case iT:
                if (is_in_closed_range(env.T[i], env.T[i + 1], value)) {
                    matched = true;
                }
                break;
            case iHmolar:
                if (is_in_closed_range(env.hmolar_vap[i], env.hmolar_vap[i + 1], value)) {
                    matched = true;
                }
                break;
            case iSmolar:
                if (is_in_closed_range(env.smolar_vap[i], env.smolar_vap[i + 1], value)) {
                    matched = true;
                }
                break;
            default:
                throw ValueError(format("bad index to find_intersections"));
        }

        if (matched) {
            intersections.push_back(std::pair<std::size_t, std::size_t>(i, i + 1));
        }
    }
    return intersections;
}
bool PhaseEnvelopeRoutines::is_inside(const PhaseEnvelopeData& env, parameters iInput1, CoolPropDbl value1, parameters iInput2, CoolPropDbl value2,
                                      std::size_t& iclosest, SimpleState& closest_state) {
    // Find the indices that bound the solution(s)
    std::vector<std::pair<std::size_t, std::size_t>> intersections = find_intersections(env, iInput1, value1);

    if (get_debug_level() > 5) {
        std::cout << format("is_inside(%Lg,%Lg); iTsat_max=%d; ipsat_max=%d\n", value1, value2, env.iTsat_max, env.ipsat_max);
    }
    // Check whether input is above max value
    if (iInput1 == iT && 0 < env.iTsat_max && env.iTsat_max < env.T.size() && value1 > env.T[env.iTsat_max]) {
        return false;
    }
    if (iInput1 == iP && 0 < env.ipsat_max && env.ipsat_max < env.p.size() && value1 > env.p[env.ipsat_max]) {
        return false;
    }

    // If number of intersections is 0, input is out of range, quit
    if (intersections.size() == 0) {
        throw ValueError(format("Input is out of range for primary value [%Lg], inputs were (%s,%Lg,%s,%Lg); no intersections found", value1,
                                get_parameter_information(iInput1, "short").c_str(), value1, get_parameter_information(iInput2, "short").c_str(),
                                value2));
    }

    // If number of intersections is 1, input will be determined based on the single intersection
    // Need to know if values increase or decrease to the right of the intersection point
    if (intersections.size() % 2 != 0) {
        throw ValueError("Input is weird; odd number of intersections found");
    }

    // If number of intersections is even, might be a bound
    if (intersections.size() % 2 == 0) {
        if (intersections.size() != 2) {
            throw ValueError("for now only even value accepted is 2");
        }
        std::vector<std::size_t> other_indices(4, 0);
        std::vector<double> const* y;
        std::vector<double> other_values(4, 0);
        other_indices[0] = intersections[0].first;
        other_indices[1] = intersections[0].second;
        other_indices[2] = intersections[1].first;
        other_indices[3] = intersections[1].second;

        switch (iInput2) {
            case iT:
                y = &(env.T);
                break;
            case iP:
                y = &(env.p);
                break;
            case iDmolar:
                y = &(env.rhomolar_vap);
                break;
            case iHmolar:
                y = &(env.hmolar_vap);
                break;
            case iSmolar:
                y = &(env.smolar_vap);
                break;
            default:
                throw ValueError("Pointer to vector y is unset in is_inside");
        }

        other_values[0] = (*y)[other_indices[0]];
        other_values[1] = (*y)[other_indices[1]];
        other_values[2] = (*y)[other_indices[2]];
        other_values[3] = (*y)[other_indices[3]];

        CoolPropDbl min_other = *(std::min_element(other_values.begin(), other_values.end()));
        CoolPropDbl max_other = *(std::max_element(other_values.begin(), other_values.end()));

        if (get_debug_level() > 5) {
            std::cout << format("is_inside: min: %Lg max: %Lg val: %Lg\n", min_other, max_other, value2);
        }

        // If by using the outer bounds of the second variable, we are outside the range,
        // then the value is definitely not inside the phase envelope and we don't need to
        // do any more analysis.
        if (!is_in_closed_range(min_other, max_other, value2)) {
            std::vector<CoolPropDbl> d(4, 0);
            d[0] = std::abs(other_values[0] - value2);
            d[1] = std::abs(other_values[1] - value2);
            d[2] = std::abs(other_values[2] - value2);
            d[3] = std::abs(other_values[3] - value2);

            // Index of minimum distance in the other_values vector
            std::size_t idist = std::distance(d.begin(), std::min_element(d.begin(), d.end()));
            // Index of closest point in the phase envelope
            iclosest = other_indices[idist];

            // Get the state for the point which is closest to the desired value - this
            // can be used as a bounding value in the outer single-phase flash routine
            // since you know (100%) that it is a good bound
            closest_state.T = env.T[iclosest];
            closest_state.p = env.p[iclosest];
            closest_state.rhomolar = env.rhomolar_vap[iclosest];
            closest_state.hmolar = env.hmolar_vap[iclosest];
            closest_state.smolar = env.smolar_vap[iclosest];
            closest_state.Q = env.Q[iclosest];

            if (get_debug_level() > 5) {
                std::cout << format("is_inside: it is not inside") << std::endl;
            }
            return false;
        } else {
            // Now we have to do a saturation flash call in order to determine whether or not we are inside the phase envelope or not

            // First we can interpolate using the phase envelope to get good guesses for the necessary values
            CoolPropDbl y1 = evaluate(env, iInput2, iInput1, value1, intersections[0].first);
            CoolPropDbl y2 = evaluate(env, iInput2, iInput1, value1, intersections[1].first);
            if (is_in_closed_range(y1, y2, value2)) {
                if (std::abs(y1 - value2) < std::abs(y2 - value2)) {
                    iclosest = intersections[0].first;
                } else {
                    iclosest = intersections[1].first;
                }
                // Get the state for the point which is closest to the desired value - this
                // can be used as a bounding value in the outer single-phase flash routine
                // since you know (100%) that it is a good bound
                closest_state.T = env.T[iclosest];
                closest_state.p = env.p[iclosest];
                closest_state.rhomolar = env.rhomolar_vap[iclosest];
                closest_state.hmolar = env.hmolar_vap[iclosest];
                closest_state.smolar = env.smolar_vap[iclosest];
                closest_state.Q = env.Q[iclosest];
                return true;
            } else {
                return false;
            }
        }
    } else {
        throw ValueError("You have a funny number of intersections in is_inside");
    }
}

} /* namespace CoolProp */

#endif
