#include "VLERoutines.h"
#include "FlashRoutines.h"
#include "HelmholtzEOSMixtureBackend.h"
#include "HelmholtzEOSBackend.h"
#include "PhaseEnvelopeRoutines.h"
#include "Configuration.h"

#if defined(ENABLE_CATCH)
#include "catch.hpp"
#endif

namespace CoolProp{

template<class T> T g_RachfordRice(const std::vector<T> &z, const std::vector<T> &lnK, T beta)
{
    // g function from Rachford-Rice
    T summer = 0;
    for (std::size_t i = 0; i < z.size(); i++)
    {
        T Ki = exp(lnK[i]);
        summer += z[i]*(Ki-1)/(1-beta+beta*Ki);
    }
    return summer;
}
template<class T> T dgdbeta_RachfordRice(const std::vector<T> &z, const std::vector<T> &lnK, T beta)
{
    // derivative of g function from Rachford-Rice with respect to beta
    T summer = 0;
    for (std::size_t i = 0; i < z.size(); i++)
    {
        T Ki = exp(lnK[i]);
        summer += -z[i]*pow((Ki-1)/(1-beta+beta*Ki),2);
    }
    return summer;
}

void FlashRoutines::PT_flash_mixtures(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.PhaseEnvelope.built){
        // Use the phase envelope if already constructed to determine phase boundary
        // Determine whether you are inside (two-phase) or outside (single-phase)
        SimpleState closest_state;
        std::size_t i;
        bool twophase = PhaseEnvelopeRoutines::is_inside(HEOS, iP, HEOS._p, iT, HEOS._T, i, closest_state);
        if (!twophase && HEOS._T > closest_state.T){
            // Gas solution - bounded between phase envelope temperature and very high temperature
            //
            // Start with a guess value from SRK
            CoolPropDbl rhomolar_guess = HEOS.solver_rho_Tp_SRK(HEOS._T, HEOS._p, iphase_gas);
            
            solver_TP_resid resid(HEOS, HEOS._T, HEOS._p);
            std::string errstr;
            HEOS.specify_phase(iphase_gas);
            try{
                // Try using Newton's method
                CoolPropDbl rhomolar = Newton(resid, rhomolar_guess, 1e-10, 100, errstr);
                // Make sure the solution is within the bounds
                if (!is_in_closed_range(static_cast<CoolPropDbl>(closest_state.rhomolar), 0.0L, rhomolar)){
                    throw ValueError("out of range");
                }
                HEOS.update_DmolarT_direct(rhomolar, HEOS._T);
            }
            catch(...){
                // If that fails, try a bounded solver
                CoolPropDbl rhomolar = Brent(resid, closest_state.rhomolar, 1e-10, DBL_EPSILON, 1e-10, 100, errstr);
                // Make sure the solution is within the bounds
                if (!is_in_closed_range(static_cast<CoolPropDbl>(closest_state.rhomolar), 0.0L, rhomolar)){
                    throw ValueError("out of range");
                }
            }
            HEOS.unspecify_phase();
        }
        else{
            // Liquid solution
            throw ValueError();
        }
    }
    else{
        // Following the strategy of Gernert, 2014
        
        // Step 1 a) Get lnK factors using Wilson
        std::vector<CoolPropDbl> lnK(HEOS.get_mole_fractions().size());
        for (std::size_t i = 0; i < lnK.size(); ++i){
            lnK[i] = SaturationSolvers::Wilson_lnK_factor(HEOS, HEOS._T, HEOS._p, i);
        }
        
        // Use Rachford-Rice to check whether you are in a homogeneous phase
        CoolPropDbl g_RR_0 = g_RachfordRice(HEOS.get_mole_fractions(), lnK, 0.0L);
        if (g_RR_0 < 0){
            // Subcooled liquid - done
            CoolPropDbl rhomolar_guess = HEOS.solver_rho_Tp_SRK(HEOS._T, HEOS._p, iphase_liquid);
            HEOS.specify_phase(iphase_liquid);
            HEOS.update_TP_guessrho(HEOS._T, HEOS._p, rhomolar_guess);
            HEOS.unspecify_phase();
            return;
        }
        else{
            CoolPropDbl g_RR_1 = g_RachfordRice(HEOS.get_mole_fractions(), lnK, 1.0L);
            if (g_RR_1 > 0){
                // Superheated vapor - done
                CoolPropDbl rhomolar_guess = HEOS.solver_rho_Tp_SRK(HEOS._T, HEOS._p, iphase_gas);
                HEOS.specify_phase(iphase_gas);
                HEOS.update_TP_guessrho(HEOS._T, HEOS._p, rhomolar_guess);
                HEOS.unspecify_phase();
                return;
            }
        }
    }
}
void FlashRoutines::PT_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.imposed_phase_index == iphase_not_imposed) // If no phase index is imposed (see set_components function)
    {
        if (HEOS.is_pure_or_pseudopure)
        {
            // At very low temperature (near the triple point temp), the isotherms are VERY steep
            // Thus it can be very difficult to determine state based on ps = f(T)
            // So in this case, we do a phase determination based on p, generally it will be useful enough
            if (HEOS._T < 0.9*HEOS.Ttriple() + 0.1*HEOS.calc_Tmax_sat())
            {
                // Find the phase, while updating all internal variables possible using the pressure
                bool saturation_called = false;
                HEOS.p_phase_determination_pure_or_pseudopure(iT, HEOS._T, saturation_called);
            }
            else{
                // Find the phase, while updating all internal variables possible using the temperature
                HEOS.T_phase_determination_pure_or_pseudopure(iP, HEOS._p);
            }
            // Check if twophase solution
            if (!HEOS.isHomogeneousPhase())
            {
                throw ValueError("twophase not implemented yet");
            }
        }
        else{
            PT_flash_mixtures(HEOS);
        }
    }
    
    // Find density
    HEOS._rhomolar = HEOS.solver_rho_Tp(HEOS._T, HEOS._p);
    HEOS._Q = -1;
}
void FlashRoutines::HQ_flash(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl Tguess)
{
    SaturationSolvers::saturation_PHSU_pure_options options;
    options.use_logdelta = false;
    HEOS.specify_phase(iphase_twophase);
    if (Tguess < 0){
        options.use_guesses = true;
        options.T = Tguess;
        CoolProp::SaturationAncillaryFunction &rhoL = HEOS.get_components()[0]->ancillaries.rhoL;
        CoolProp::SaturationAncillaryFunction &rhoV = HEOS.get_components()[0]->ancillaries.rhoV;
        options.rhoL = rhoL.evaluate(Tguess);
        options.rhoV = rhoV.evaluate(Tguess);
    }
    if (HEOS.is_pure_or_pseudopure){
        if (std::abs(HEOS.Q()-1) > 1e-10){throw ValueError(format("non-unity quality not currently allowed for HQ_flash"));}
        // Do a saturation call for given h for vapor, first with ancillaries, then with full saturation call
        options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_HV;
        SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS.hmolar(), options);
        HEOS._p = HEOS.SatV->p();
        HEOS._T = HEOS.SatV->T();
        HEOS._rhomolar = HEOS.SatV->rhomolar();
        HEOS._phase = iphase_twophase;
    }
    else{
        throw NotImplementedError("QS_flash not ready for mixtures");
    }
}
void FlashRoutines::QS_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure){
        
        if (std::abs(HEOS.smolar() - HEOS.get_state("reducing").smolar) < 0.001)
        {
            HEOS._p = HEOS.p_critical();
            HEOS._T = HEOS.T_critical();
            HEOS._rhomolar = HEOS.rhomolar_critical();
            HEOS._phase = iphase_critical_point;
        }
        else if (std::abs(HEOS.Q()) < 1e-10){
            // Do a saturation call for given s for liquid, first with ancillaries, then with full saturation call
            SaturationSolvers::saturation_PHSU_pure_options options;
            options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_SL;
            options.use_logdelta = false;
            HEOS.specify_phase(iphase_twophase);
            SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS.smolar(), options);
            HEOS._p = HEOS.SatL->p();
            HEOS._T = HEOS.SatL->T();
            HEOS._rhomolar = HEOS.SatL->rhomolar();
            HEOS._phase = iphase_twophase;
        }
        else if (std::abs(HEOS.Q()-1) < 1e-10)
        {
            // Do a saturation call for given s for vapor, first with ancillaries, then with full saturation call
            SaturationSolvers::saturation_PHSU_pure_options options;
            options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_SV;
            options.use_logdelta = false;
            HEOS.specify_phase(iphase_twophase);
            SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS.smolar(), options);
            HEOS._p = HEOS.SatV->p();
            HEOS._T = HEOS.SatV->T();
            HEOS._rhomolar = HEOS.SatV->rhomolar();
            HEOS._phase = iphase_twophase;
        }
        else{
            throw ValueError(format("non-zero or 1 quality not currently allowed for QS_flash"));
        }
    }
    else{
        throw NotImplementedError("QS_flash not ready for mixtures");
    }
}
void FlashRoutines::QT_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    CoolPropDbl T = HEOS._T;
    if (HEOS.is_pure_or_pseudopure)
    {
        // The maximum possible saturation temperature
        // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
        CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() + 1e-13;
        
        // Check what the minimum limits for the equation of state are
        CoolPropDbl Tmin_satL, Tmin_satV, Tmin_sat;
        HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
        Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;
        
        // Get a reference to keep the code a bit cleaner
        CriticalRegionSplines &splines = HEOS.components[0]->pEOS->critical_region_splines;
        
        // If exactly(ish) at the critical temperature, liquid and vapor have the critial density
        if ((get_config_bool(CRITICAL_WITHIN_1UK) && std::abs(T-Tmax_sat)< 1e-6) || std::abs(T-Tmax_sat)< 1e-12){
             HEOS.SatL->update(DmolarT_INPUTS, HEOS.rhomolar_critical(), HEOS._T);
             HEOS.SatV->update(DmolarT_INPUTS, HEOS.rhomolar_critical(), HEOS._T);
             HEOS._rhomolar = HEOS.rhomolar_critical();
             HEOS._p = 0.5*HEOS.SatV->p() + 0.5*HEOS.SatL->p();
        }
        else if (!is_in_closed_range(Tmin_sat, Tmax_sat, T)){
            throw ValueError(format("Temperature to QT_flash [%0.8Lg K] must be in range [%0.8Lg K, %0.8Lg K]", T, Tmin_sat, Tmax_sat));
        }
        else if (get_config_bool(CRITICAL_SPLINES_ENABLED) && splines.enabled && HEOS._T > splines.T_min){
            double rhoL = _HUGE, rhoV = _HUGE;
            // Use critical region spline if it has it and temperature is in its range
            splines.get_densities(T, splines.rhomolar_min, HEOS.rhomolar_critical(), splines.rhomolar_max, rhoL, rhoV);
            HEOS.SatL->update(DmolarT_INPUTS, rhoL, HEOS._T);
            HEOS.SatV->update(DmolarT_INPUTS, rhoV, HEOS._T);
            HEOS._p = 0.5*HEOS.SatV->p() + 0.5*HEOS.SatL->p();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        }
        else if (!(HEOS.components[0]->pEOS->pseudo_pure))
        {
            // Set some imput options
            SaturationSolvers::saturation_T_pure_Akasaka_options options;
            options.use_guesses = false;

            // Actually call the solver
            SaturationSolvers::saturation_T_pure_Maxwell(HEOS, HEOS._T, options);
            
            HEOS._p = 0.5*HEOS.SatV->p() + 0.5*HEOS.SatL->p();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        }
        else{
            // Pseudo-pure fluid
            CoolPropDbl rhoLanc = _HUGE, rhoVanc = _HUGE, rhoLsat = _HUGE, rhoVsat = _HUGE;
            CoolPropDbl psatLanc = HEOS.components[0]->ancillaries.pL.evaluate(HEOS._T); // These ancillaries are used explicitly
            CoolPropDbl psatVanc = HEOS.components[0]->ancillaries.pV.evaluate(HEOS._T); // These ancillaries are used explicitly
            try{
                rhoLanc = HEOS.components[0]->ancillaries.rhoL.evaluate(HEOS._T);
                rhoVanc = HEOS.components[0]->ancillaries.rhoV.evaluate(HEOS._T);

                if (!ValidNumber(rhoLanc) || !ValidNumber(rhoVanc))
                {
                    throw ValueError("pseudo-pure failed");
                }

                HEOS.SatL->update_TP_guessrho(HEOS._T, psatLanc, rhoLanc);
                HEOS.SatV->update_TP_guessrho(HEOS._T, psatVanc, rhoVanc);
                if (!ValidNumber(rhoLsat) || !ValidNumber(rhoVsat) ||
                     std::abs(rhoLsat/rhoLanc-1) > 0.5 || std::abs(rhoVanc/rhoVsat-1) > 0.5)
                {
                    throw ValueError("pseudo-pure failed");
                }
            }
            catch (...){
                // Near the critical point, the behavior is not very nice, so we will just use the ancillary
                rhoLsat = rhoLanc;
                rhoVsat = rhoVanc;
            }
            HEOS._p = HEOS._Q*psatVanc + (1-HEOS._Q)*psatLanc;
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
            HEOS.SatL->update(DmolarT_INPUTS, rhoLsat, HEOS._T);
            HEOS.SatV->update(DmolarT_INPUTS, rhoVsat, HEOS._T);
        }
        // Load the outputs
        HEOS._phase = iphase_twophase;
    }
    else
    {
        if(HEOS.PhaseEnvelope.built){
            PT_Q_flash_mixtures(HEOS, iT, HEOS._T);
        }
        else{
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

            HEOS._p = options.p;
            HEOS._rhomolar = 1/(HEOS._Q/options.rhomolar_vap+(1-HEOS._Q)/options.rhomolar_liq);
        }
        // Load the outputs
        HEOS._phase = iphase_twophase;
        HEOS._p = HEOS.SatV->p();
        HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        HEOS._T = HEOS.SatL->T();
    }
}
void FlashRoutines::PQ_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        if (HEOS.components[0]->pEOS->pseudo_pure){
            // It is a pseudo-pure mixture
            
            HEOS._TLanc = HEOS.components[0]->ancillaries.pL.invert(HEOS._p);
            HEOS._TVanc = HEOS.components[0]->ancillaries.pV.invert(HEOS._p);
            // Get guesses for the ancillaries for density
            CoolPropDbl rhoL = HEOS.components[0]->ancillaries.rhoL.evaluate(HEOS._TLanc);
            CoolPropDbl rhoV = HEOS.components[0]->ancillaries.rhoV.evaluate(HEOS._TVanc);
            // Solve for the density
            HEOS.SatL->update_TP_guessrho(HEOS._TLanc, HEOS._p, rhoL);
            HEOS.SatV->update_TP_guessrho(HEOS._TVanc, HEOS._p, rhoV);
            
            // Load the outputs
            HEOS._phase = iphase_twophase;
            HEOS._p = HEOS._Q*HEOS.SatV->p() + (1 - HEOS._Q)*HEOS.SatL->p();
            HEOS._T = HEOS._Q*HEOS.SatV->T() + (1 - HEOS._Q)*HEOS.SatL->T();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        }
        else{
            // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
            CoolPropDbl pmax_sat = HEOS.calc_pmax_sat();
            
            // Check what the minimum limits for the equation of state are
            CoolPropDbl pmin_satL, pmin_satV, pmin_sat;
            HEOS.calc_pmin_sat(pmin_satL, pmin_satV);
            pmin_sat = std::max(pmin_satL, pmin_satV);
            
            // Check for being AT the critical point
            if (is_in_closed_range(pmax_sat*(1-1e-10), pmax_sat*(1+1e-10), static_cast<CoolPropDbl>(HEOS._p))){
                // Load the outputs
                HEOS._phase = iphase_critical_point;
                HEOS._p = HEOS.p_critical();
                HEOS._rhomolar = HEOS.rhomolar_critical();
                HEOS._T = HEOS.T_critical();
                return;
            }
            
            // Check limits
            if (!is_in_closed_range(pmin_sat*0.999999, pmax_sat*1.000001, static_cast<CoolPropDbl>(HEOS._p))){
                throw ValueError(format("Pressure to PQ_flash [%6g Pa] must be in range [%8Lg Pa, %8Lg Pa]",HEOS._p, pmin_sat, pmax_sat));
            }
            // ------------------
            // It is a pure fluid
            // ------------------

            // Set some imput options
            SaturationSolvers::saturation_PHSU_pure_options options;
            // Specified variable is pressure
            options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_PL;
            // Use logarithm of delta as independent variables
            options.use_logdelta = false;
            
            double increment = 0.4;

            try{
                for (double omega = 1.0; omega > 0; omega -= increment){
                    try{
                        options.omega = omega;
                        
                        // Actually call the solver
                        SaturationSolvers::saturation_PHSU_pure(HEOS, HEOS._p, options);
                        
                        // If you get here, there was no error, all is well
                        break;
                    }
                    catch(...){
                        if (omega < 1.1*increment){
                            throw;
                        }
                        // else we are going to try again with a smaller omega
                    }
                }
            }
            catch(...){
                // We may need to polish the solution at low pressure
                SaturationSolvers::saturation_P_pure_1D_T(HEOS, HEOS._p, options);
            }

            // Load the outputs
            HEOS._phase = iphase_twophase;
            HEOS._p = HEOS._Q*HEOS.SatV->p() + (1 - HEOS._Q)*HEOS.SatL->p();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
            HEOS._T = HEOS.SatL->T();
        }
    }
    else
    {
        if (HEOS.PhaseEnvelope.built){
            PT_Q_flash_mixtures(HEOS, iP, HEOS._p);
        }
        else{
            // Set some imput options
            SaturationSolvers::mixture_VLE_IO io;
            io.sstype = SaturationSolvers::imposed_p;
            io.Nstep_max = 10;

            // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
            CoolPropDbl Tguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

            // Use Wilson iteration to obtain updated guess for temperature
            Tguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);

            // Actually call the successive substitution solver
            SaturationSolvers::successive_substitution(HEOS, HEOS._Q, Tguess, HEOS._p, HEOS.mole_fractions, HEOS.K, io);
        }
                    
        // Load the outputs
        HEOS._phase = iphase_twophase;
        HEOS._p = HEOS.SatV->p();
        HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        HEOS._T = HEOS.SatL->T();
    }
}

void FlashRoutines::PT_Q_flash_mixtures(HelmholtzEOSMixtureBackend &HEOS, parameters other, CoolPropDbl value)
{
    
    // Find the intersections in the phase envelope
    std::vector< std::pair<std::size_t, std::size_t> > intersections = PhaseEnvelopeRoutines::find_intersections(HEOS, other, value);
    
    PhaseEnvelopeData &env = HEOS.PhaseEnvelope;
    
    enum quality_options{SATURATED_LIQUID, SATURATED_VAPOR, TWO_PHASE};
    quality_options quality;
    if (std::abs(HEOS._Q) < 100*DBL_EPSILON){
        quality = SATURATED_LIQUID;
    }
    else if (std::abs(HEOS._Q - 1) < 100*DBL_EPSILON){
        quality = SATURATED_VAPOR;
    }
    else if (HEOS._Q > 0 && HEOS._Q < 1){
        quality = TWO_PHASE;
    }
    else{
        throw ValueError("Quality is not within 0 and 1");
    }
    
    if (quality == SATURATED_LIQUID || quality == SATURATED_VAPOR)
    {
        // *********************************************************
        //            Bubble- or dew-point calculation
        // *********************************************************
        // Find the correct solution
        std::vector<std::size_t> solutions;
        for (std::vector< std::pair<std::size_t, std::size_t> >::const_iterator it = intersections.begin(); it != intersections.end(); ++it){
            if (std::abs(env.Q[it->first] - HEOS._Q) < 10*DBL_EPSILON && std::abs(env.Q[it->second] - HEOS._Q) < 10*DBL_EPSILON ){
                solutions.push_back(it->first);
            }
        }
            
        if (solutions.size() == 1){
            
            std::size_t &imax = solutions[0];
            
            // Shift the solution if needed to ensure that imax+2 and imax-1 are both in range
            if (imax+2 >= env.T.size()){ imax--; }
            else if (imax-1 < 0){ imax++; }
            
            SaturationSolvers::newton_raphson_saturation NR;
            SaturationSolvers::newton_raphson_saturation_options IO;
            
            if (other == iP){
                IO.p = HEOS._p;
                IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::P_IMPOSED;
                // p -> rhomolar_vap
                IO.rhomolar_vap = CubicInterp(env.p, env.rhomolar_vap, imax-1, imax, imax+1, imax+2, static_cast<CoolPropDbl>(IO.p));
                IO.T = CubicInterp(env.rhomolar_vap, env.T, imax-1, imax, imax+1, imax+2, IO.rhomolar_vap);
            }
            else if (other == iT){
                IO.T = HEOS._T;
                IO.imposed_variable = SaturationSolvers::newton_raphson_saturation_options::T_IMPOSED;
                // T -> rhomolar_vap
                IO.rhomolar_vap = CubicInterp(env.T, env.rhomolar_vap, imax-1, imax, imax+1, imax+2, static_cast<CoolPropDbl>(IO.T));
                IO.p = CubicInterp(env.rhomolar_vap, env.p, imax-1, imax, imax+1, imax+2, IO.rhomolar_vap);
            }
            IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, imax-1, imax, imax+1, imax+2, IO.rhomolar_vap);
            
            if (quality == SATURATED_VAPOR){
                IO.bubble_point = false;
                IO.y = HEOS.get_mole_fractions(); // Because Q = 1
                IO.x.resize(IO.y.size());
                for (std::size_t i = 0; i < IO.x.size()-1; ++i) // First N-1 elements
                {
                    IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], imax-1, imax, imax+1, imax+2, IO.rhomolar_vap);
                }
                IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
                NR.call(HEOS, IO.y, IO.x, IO);
            }
            else{
                IO.bubble_point = true;
                IO.x = HEOS.get_mole_fractions(); // Because Q = 0
                IO.y.resize(IO.x.size());
                // Phases are inverted, so "liquid" is actually the lighter phase
                std::swap(IO.rhomolar_liq, IO.rhomolar_vap);
                for (std::size_t i = 0; i < IO.y.size()-1; ++i) // First N-1 elements
                {
                    // Phases are inverted, so liquid mole fraction (x) of phase envelope is actually the vapor phase mole fraction
                    // Use the liquid density as well
                    IO.y[i] = CubicInterp(env.rhomolar_vap, env.x[i], imax-1, imax, imax+1, imax+2, IO.rhomolar_liq);
                }
                IO.y[IO.y.size()-1] = 1 - std::accumulate(IO.y.begin(), IO.y.end()-1, 0.0);
                NR.call(HEOS, IO.x, IO.y, IO);
            }
        }
        else if (solutions.size() == 0){
            throw ValueError("No solution was found in PQ_flash");
        }
        else{
            throw ValueError("More than 1 solution was found in PQ_flash");
        }
    }
    else{
        // *********************************************************
        //      Two-phase calculation for given vapor quality
        // *********************************************************
        
         // Find the correct solution
        std::vector<std::size_t> liquid_solutions, vapor_solutions;
        for (std::vector< std::pair<std::size_t, std::size_t> >::const_iterator it = intersections.begin(); it != intersections.end(); ++it){
            if (std::abs(env.Q[it->first] - 0) < 10*DBL_EPSILON && std::abs(env.Q[it->second] - 0) < 10*DBL_EPSILON ){
                liquid_solutions.push_back(it->first);
            }
            if (std::abs(env.Q[it->first] - 1) < 10*DBL_EPSILON && std::abs(env.Q[it->second] - 1) < 10*DBL_EPSILON ){
                vapor_solutions.push_back(it->first);
            }
        }
        
        if (liquid_solutions.size() != 1 || vapor_solutions.size() != 1){
            throw ValueError(format("Number liquid solutions [%d] or vapor solutions [%d] != 1", liquid_solutions.size(), vapor_solutions.size() ));
        }
        std::size_t iliq = liquid_solutions[0], ivap = vapor_solutions[0];
        
        SaturationSolvers::newton_raphson_twophase NR;
        SaturationSolvers::newton_raphson_twophase_options IO;
        IO.beta = HEOS._Q;
        
        CoolPropDbl rhomolar_vap_sat_vap, T_sat_vap, rhomolar_liq_sat_vap, rhomolar_liq_sat_liq, T_sat_liq, rhomolar_vap_sat_liq, p_sat_liq, p_sat_vap;

        if (other == iP){
            IO.p = HEOS._p;
            p_sat_liq = IO.p; p_sat_vap = IO.p;
            IO.imposed_variable = SaturationSolvers::newton_raphson_twophase_options::P_IMPOSED;
            
            // Calculate the interpolated values for beta = 0 and beta = 1
            rhomolar_vap_sat_vap = CubicInterp(env.p, env.rhomolar_vap, ivap-1, ivap, ivap+1, ivap+2, static_cast<CoolPropDbl>(IO.p));
            T_sat_vap = CubicInterp(env.rhomolar_vap, env.T, ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            rhomolar_liq_sat_vap = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            
            // Phase inversion for liquid solution (liquid is vapor and vice versa)
            rhomolar_liq_sat_liq = CubicInterp(env.p, env.rhomolar_vap, iliq-1, iliq, iliq+1, iliq+2, static_cast<CoolPropDbl>(IO.p)); 
            T_sat_liq = CubicInterp(env.rhomolar_vap, env.T, iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
            rhomolar_vap_sat_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
        }
        else if (other == iT){
            IO.T = HEOS._T;
            T_sat_liq = IO.T; T_sat_vap = IO.T;
            IO.imposed_variable = SaturationSolvers::newton_raphson_twophase_options::T_IMPOSED;
            
            // Calculate the interpolated values for beta = 0 and beta = 1
            rhomolar_vap_sat_vap = CubicInterp(env.T, env.rhomolar_vap, ivap-1, ivap, ivap+1, ivap+2, static_cast<CoolPropDbl>(IO.T));
            p_sat_vap = CubicInterp(env.rhomolar_vap, env.p, ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            rhomolar_liq_sat_vap = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            
            // Phase inversion for liquid solution (liquid is vapor and vice versa)
            rhomolar_liq_sat_liq = CubicInterp(env.T, env.rhomolar_vap, iliq-1, iliq, iliq+1, iliq+2, static_cast<CoolPropDbl>(IO.T)); 
            p_sat_liq = CubicInterp(env.rhomolar_vap, env.p, iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
            rhomolar_vap_sat_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
        }
        else{
            throw ValueError();
        }
        
        // Weight the guesses by the vapor mole fraction
        IO.rhomolar_vap = IO.beta*rhomolar_vap_sat_vap + (1-IO.beta)*rhomolar_vap_sat_liq;
        IO.rhomolar_liq = IO.beta*rhomolar_liq_sat_vap + (1-IO.beta)*rhomolar_liq_sat_liq;
        IO.T = IO.beta*T_sat_vap + (1-IO.beta)*T_sat_liq;
        IO.p = IO.beta*p_sat_vap + (1-IO.beta)*p_sat_liq;
        
        IO.z = HEOS.get_mole_fractions();
        IO.x.resize(IO.z.size());
        IO.y.resize(IO.z.size());
        
        for (std::size_t i = 0; i < IO.x.size()-1; ++i) // First N-1 elements
        {
            CoolPropDbl x_sat_vap = CubicInterp(env.rhomolar_vap, env.x[i], ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            CoolPropDbl y_sat_vap = CubicInterp(env.rhomolar_vap, env.y[i], ivap-1, ivap, ivap+1, ivap+2, rhomolar_vap_sat_vap);
            
            CoolPropDbl x_sat_liq = CubicInterp(env.rhomolar_vap, env.y[i], iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
            CoolPropDbl y_sat_liq = CubicInterp(env.rhomolar_vap, env.x[i], iliq-1, iliq, iliq+1, iliq+2, rhomolar_liq_sat_liq);
            
            IO.x[i] = IO.beta*x_sat_vap + (1-IO.beta)*x_sat_liq;
            IO.y[i] = IO.beta*y_sat_vap + (1-IO.beta)*y_sat_liq;
        }
        IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
        IO.y[IO.y.size()-1] = 1 - std::accumulate(IO.y.begin(), IO.y.end()-1, 0.0);
        std::vector<CoolPropDbl> &XX = IO.x;
        std::vector<CoolPropDbl> &YY = IO.y;
        NR.call(HEOS, IO);
    }
}
void FlashRoutines::HSU_D_flash_twophase(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value)
{
    class Residual : public FuncWrapper1D
    {
        
    public:
        HelmholtzEOSMixtureBackend &HEOS;
        CoolPropDbl rhomolar_spec, Qd, Qo;
		parameters other;
		CoolPropDbl value;
        Residual(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl rhomolar_spec, parameters other, CoolPropDbl value) : HEOS(HEOS), rhomolar_spec(rhomolar_spec), other(other), value(value){};
        double call(double T){
            HEOS.update(QT_INPUTS, 0, T);
            HelmholtzEOSMixtureBackend &SatL = HEOS.get_SatL(),
                                       &SatV = HEOS.get_SatV();
            // Quality from density
            Qd = (1/HEOS.rhomolar()-1/SatL.rhomolar())/(1/SatV.rhomolar()-1/SatL.rhomolar());
            // Quality from other parameter (H,S,U)
            Qo = (value-SatL.keyed_output(other))/(SatV.keyed_output(other)-SatL.keyed_output(other));
            // Residual is the difference between the two
            return Qo-Qd;
        }
    } resid(HEOS, rhomolar_spec, other, value);
    
    std::string errstr;
    // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
    CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() - 1e-13;
    
    // Check what the minimum limits for the equation of state are
    CoolPropDbl Tmin_satL, Tmin_satV, Tmin_sat;
    HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
    Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;
        
    Brent(resid, Tmin_sat, Tmax_sat-0.01, DBL_EPSILON, 1e-12, 20, errstr);
}
// D given and one of P,H,S,U
void FlashRoutines::PHSU_D_flash(HelmholtzEOSMixtureBackend &HEOS, parameters other)
{
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl r, eos, rhomolar, value, T;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend *HEOS, CoolPropDbl rhomolar, CoolPropDbl value, int other) : HEOS(HEOS), rhomolar(rhomolar), value(value), other(other){};
        double call(double T){
            this->T = T;
            switch(other)
            {
            case iP:
                eos = HEOS->calc_pressure_nocache(T, rhomolar); break;
            case iSmolar:
                eos = HEOS->calc_smolar_nocache(T, rhomolar); break;
            case iHmolar:
                eos = HEOS->calc_hmolar_nocache(T, rhomolar); break;
            case iUmolar:
                eos = HEOS->calc_umolar_nocache(T, rhomolar); break;
            default:
                throw ValueError(format("Input not supported"));
            }

            r = eos - value;
            return r;
        };
    };

    std::string errstring;

    if (HEOS.imposed_phase_index != iphase_not_imposed)
    {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
    }
    else
    {
        if (HEOS.is_pure_or_pseudopure)
        {
            CoolPropFluid * component = HEOS.components[0];

            shared_ptr<HelmholtzEOSMixtureBackend> Sat;
            CoolPropDbl rhoLtriple = component->triple_liquid.rhomolar;
            CoolPropDbl rhoVtriple = component->triple_vapor.rhomolar;
            // Check if in the "normal" region
            if (HEOS._rhomolar >= rhoVtriple && HEOS._rhomolar <= rhoLtriple)
            {
                CoolPropDbl yL, yV, value, y_solid;
                CoolPropDbl TLtriple = component->triple_liquid.T; ///TODO: separate TL and TV for ppure
                CoolPropDbl TVtriple = component->triple_vapor.T;

                // First check if solid (below the line connecting the triple point values) - this is an error for now
                switch (other)
                {
                    case iSmolar:
                        yL = HEOS.calc_smolar_nocache(TLtriple, rhoLtriple); yV = HEOS.calc_smolar_nocache(TVtriple, rhoVtriple); value = HEOS._smolar; break;
                    case iHmolar:
                        yL = HEOS.calc_hmolar_nocache(TLtriple, rhoLtriple); yV = HEOS.calc_hmolar_nocache(TVtriple, rhoVtriple); value = HEOS._hmolar; break;
                    case iUmolar:
                        yL = HEOS.calc_umolar_nocache(TLtriple, rhoLtriple); yV = HEOS.calc_umolar_nocache(TVtriple, rhoVtriple); value = HEOS._umolar; break;
                    case iP:
                        yL = HEOS.calc_pressure_nocache(TLtriple, rhoLtriple); yV = HEOS.calc_pressure_nocache(TVtriple, rhoVtriple); value = HEOS._p; break;
                    default:
                        throw ValueError(format("Input is invalid"));
                }
                y_solid = (yV-yL)/(1/rhoVtriple-1/rhoLtriple)*(1/HEOS._rhomolar-1/rhoLtriple) + yL;

                if (value < y_solid){ throw ValueError(format("Other input [%d:%g] is solid", other, value));}

                // Check if other is above the saturation value.
                SaturationSolvers::saturation_D_pure_options optionsD;
                optionsD.omega = 1;
                optionsD.use_logdelta = false;
                if (HEOS._rhomolar > HEOS._crit.rhomolar)
                {
                    optionsD.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOL;
                    SaturationSolvers::saturation_D_pure(HEOS, HEOS._rhomolar, optionsD);
                    // SatL and SatV have the saturation values
                    Sat = HEOS.SatL;
                }
                else
                {
                    optionsD.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOV;
                    SaturationSolvers::saturation_D_pure(HEOS, HEOS._rhomolar, optionsD);
                    // SatL and SatV have the saturation values
                    Sat = HEOS.SatV;
                }

                // If it is above, it is not two-phase and either liquid, vapor or supercritical
                if (value > Sat->keyed_output(other))
                {
                    solver_resid resid(&HEOS, HEOS._rhomolar, value, other);
                    HEOS._T = Brent(resid, Sat->keyed_output(iT), HEOS.Tmax()*1.5, DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
					// Update the phase flag
					HEOS.recalculate_singlephase_phase();
                }
                else
                {
					// Now we know that temperature is between Tsat(D) +- tolerance and the minimum temperature for the fluid
					if (other == iP){
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
						HEOS._Q = (1/HEOS._rhomolar-1/HEOS.SatL->rhomolar())/(1/HEOS.SatV->rhomolar()-1/HEOS.SatL->rhomolar());
						HEOS._T = HEOS.SatL->T();
					}
					else{
						// Residual is difference in quality calculated from density and quality calculated from the other parameter
						// Iterate to find T
						HSU_D_flash_twophase(HEOS, HEOS._rhomolar, other, value);
					}
                }
            }
            // Check if vapor/solid region below triple point vapor density
            else if (HEOS._rhomolar < component->triple_vapor.rhomolar)
            {
                CoolPropDbl y, value;
                CoolPropDbl TVtriple = component->triple_vapor.T; //TODO: separate TL and TV for ppure

                // If value is above the value calculated from X(Ttriple, _rhomolar), it is vapor
                switch (other)
                {
                    case iSmolar:
                        y = HEOS.calc_smolar_nocache(TVtriple, HEOS._rhomolar); value = HEOS._smolar; break;
                    case iHmolar:
                        y = HEOS.calc_hmolar_nocache(TVtriple, HEOS._rhomolar); value = HEOS._hmolar; break;
                    case iUmolar:
                        y = HEOS.calc_umolar_nocache(TVtriple, HEOS._rhomolar); value = HEOS._umolar; break;
                    case iP:
                        y = HEOS.calc_pressure_nocache(TVtriple, HEOS._rhomolar); value = HEOS._p; break;
                    default:
                        throw ValueError(format("Input is invalid"));
                }
                if (value > y)
                {
                    solver_resid resid(&HEOS, HEOS._rhomolar, value, other);
                    HEOS._phase = iphase_gas;
                    HEOS._T = Brent(resid, TVtriple, HEOS.Tmax()*1.5, DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
                }
                else
                {
                    throw ValueError(format("D < DLtriple %g %g", value, y));
                }

            }
            // Check in the liquid/solid region above the triple point density
            else
            {
                CoolPropDbl y, value;
                CoolPropDbl TLtriple = component->pEOS->Ttriple;

                // If value is above the value calculated from X(Ttriple, _rhomolar), it is vapor
                switch (other)
                {
                    case iSmolar:
                        y = HEOS.calc_smolar_nocache(TLtriple, HEOS._rhomolar); value = HEOS._smolar; break;
                    case iHmolar:
                        y = HEOS.calc_hmolar_nocache(TLtriple, HEOS._rhomolar); value = HEOS._hmolar; break;
                    case iUmolar:
                        y = HEOS.calc_umolar_nocache(TLtriple, HEOS._rhomolar); value = HEOS._umolar; break;
                    case iP:
                        y = HEOS.calc_pressure_nocache(TLtriple, HEOS._rhomolar); value = HEOS._p; break;
                    default:
                        throw ValueError(format("Input is invalid"));
                }
                if (value > y)
                {
                    solver_resid resid(&HEOS, HEOS._rhomolar, value, other);
                    HEOS._phase = iphase_liquid;
                    HEOS._T = Brent(resid, TLtriple, HEOS.Tmax()*1.5, DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
                }
                else
                {
                    throw ValueError(format("D < DLtriple %g %g", value, y));
                }
            }
        }
        else
            throw NotImplementedError("PHSU_D_flash not ready for mixtures");
    }
}

void FlashRoutines::HSU_P_flash_singlephase_Newton(HelmholtzEOSMixtureBackend &HEOS, parameters other, CoolPropDbl T0, CoolPropDbl rhomolar0)
{
    double A[2][2], B[2][2];
    CoolPropDbl y = _HUGE;
    HelmholtzEOSMixtureBackend _HEOS(HEOS.get_components());
    _HEOS.update(DmolarT_INPUTS, rhomolar0, T0);
    CoolPropDbl Tc = HEOS.calc_T_critical();
    CoolPropDbl rhoc = HEOS.calc_rhomolar_critical();
    CoolPropDbl R = HEOS.gas_constant();
    CoolPropDbl p = HEOS.p();
    switch (other)
    {
        case iHmolar: y = HEOS.hmolar(); break;
        case iSmolar: y = HEOS.smolar(); break;
        default: throw ValueError("other is invalid in HSU_P_flash_singlephase_Newton");
    }
    
    CoolPropDbl worst_error = 999;
    int iter = 0;
    bool failed = false;
    CoolPropDbl omega = 1.0, f2, df2_dtau, df2_ddelta;
    CoolPropDbl tau = _HEOS.tau(), delta = _HEOS.delta();
    while (worst_error>1e-6 && failed == false)
    {
        
        // All the required partial derivatives
        CoolPropDbl a0 = _HEOS.calc_alpha0_deriv_nocache(0,0,HEOS.mole_fractions, tau, delta,Tc,rhoc);
        CoolPropDbl da0_ddelta = _HEOS.calc_alpha0_deriv_nocache(0,1,HEOS.mole_fractions, tau, delta,Tc,rhoc);
        CoolPropDbl da0_dtau = _HEOS.calc_alpha0_deriv_nocache(1,0,HEOS.mole_fractions, tau, delta,Tc,rhoc);
        CoolPropDbl d2a0_dtau2 = _HEOS.calc_alpha0_deriv_nocache(2,0,HEOS.mole_fractions, tau, delta,Tc,rhoc);
        CoolPropDbl d2a0_ddelta_dtau = 0.0;
        
        CoolPropDbl ar = _HEOS.calc_alphar_deriv_nocache(0,0,HEOS.mole_fractions, tau, delta);
        CoolPropDbl dar_dtau = _HEOS.calc_alphar_deriv_nocache(1,0,HEOS.mole_fractions, tau, delta);
        CoolPropDbl dar_ddelta = _HEOS.calc_alphar_deriv_nocache(0,1,HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_ddelta_dtau = _HEOS.calc_alphar_deriv_nocache(1,1,HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_ddelta2 = _HEOS.calc_alphar_deriv_nocache(0,2,HEOS.mole_fractions, tau, delta);
        CoolPropDbl d2ar_dtau2 = _HEOS.calc_alphar_deriv_nocache(2,0,HEOS.mole_fractions, tau, delta);

        CoolPropDbl f1 = delta/tau*(1+delta*dar_ddelta)-p/(rhoc*R*Tc);
        CoolPropDbl df1_dtau = (1+delta*dar_ddelta)*(-delta/tau/tau)+delta/tau*(delta*d2ar_ddelta_dtau);
        CoolPropDbl df1_ddelta = (1.0/tau)*(1+2.0*delta*dar_ddelta+delta*delta*d2ar_ddelta2);
        switch (other)
        {
            case iHmolar:
            {
                f2 = (1+delta*dar_ddelta)+tau*(da0_dtau+dar_dtau) - tau*y/(R*Tc);
                df2_dtau = delta*d2ar_ddelta_dtau+da0_dtau+dar_dtau+tau*(d2a0_dtau2+d2ar_dtau2) - y/(R*Tc);
                df2_ddelta = (dar_ddelta+delta*d2ar_ddelta2)+tau*(d2a0_ddelta_dtau+d2ar_ddelta_dtau);    
                break;
            }
            case iSmolar:
            {
                f2 = tau*(da0_dtau+dar_dtau)-ar-a0-y/R;
                df2_dtau = tau*(d2a0_dtau2 + d2ar_dtau2)+(da0_dtau+dar_dtau)-dar_dtau-da0_dtau;
                df2_ddelta = tau*(d2a0_ddelta_dtau+d2ar_ddelta_dtau)-dar_ddelta-da0_ddelta;
                break;
            }
            default:
                throw ValueError("other variable in HSU_P_flash_singlephase_Newton is invalid");
        }

        //First index is the row, second index is the column
        A[0][0]=df1_dtau;
        A[0][1]=df1_ddelta;
        A[1][0]=df2_dtau;
        A[1][1]=df2_ddelta;

        //double det = A[0][0]*A[1][1]-A[1][0]*A[0][1];

        MatInv_2(A,B);
        tau -= omega*(B[0][0]*f1+B[0][1]*f2);
        delta -= omega*(B[1][0]*f1+B[1][1]*f2);

        if (std::abs(f1) > std::abs(f2))
            worst_error = std::abs(f1);
        else
            worst_error = std::abs(f2);

        if (!ValidNumber(f1) || !ValidNumber(f2))
        {
            throw SolutionError(format("Invalid values for inputs p=%g y=%g for fluid %s in HSU_P_flash_singlephase", p, y, _HEOS.name().c_str()));
        }

        iter += 1;
        if (iter>100)
        {
            throw SolutionError(format("HSU_P_flash_singlephase did not converge with inputs p=%g h=%g for fluid %s", p, y,_HEOS.name().c_str()));
        }
    }
    
    HEOS.update(DmolarT_INPUTS, rhoc*delta, Tc/tau);
}
void FlashRoutines::HSU_P_flash_singlephase_Brent(HelmholtzEOSMixtureBackend &HEOS, parameters other, CoolPropDbl value, CoolPropDbl Tmin, CoolPropDbl Tmax)
{
    if (!ValidNumber(HEOS._p)){throw ValueError("value for p in HSU_P_flash_singlephase_Brent is invalid");};
    if (!ValidNumber(value)){throw ValueError("value for other in HSU_P_flash_singlephase_Brent is invalid");};
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl r, eos, p, value, T, rhomolar;
        int other;
        int iter;
        CoolPropDbl r0, r1, T1, T0, eos0, eos1, pp;
        solver_resid(HelmholtzEOSMixtureBackend *HEOS, CoolPropDbl p, CoolPropDbl value, int other) : 
                HEOS(HEOS), p(p), value(value), other(other)
                {
                    iter = 0;
                    // Specify the state to avoid saturation calls, but only if phase is subcritical
                    if (HEOS->phase() == iphase_liquid || HEOS->phase() == iphase_gas ){
                        HEOS->specify_phase(HEOS->phase());
                    }
                };
        double call(double T){

            this->T = T;

            // Run the solver with T,P as inputs;
            HEOS->update(PT_INPUTS, p, T);
            
            rhomolar = HEOS->rhomolar();
            HEOS->update(DmolarT_INPUTS, rhomolar, T);
            // Get the value of the desired variable
            eos = HEOS->keyed_output(other);
            pp = HEOS->p();

            // Difference between the two is to be driven to zero
            r = eos - value;
            
            // Store values for later use if there are errors
            if (iter == 0){ 
                r0 = r; T0 = T; eos0 = eos;
            }
            else if (iter == 1){
                r1 = r; T1 = T; eos1 = eos; 
            }
            else{
                r0 = r1; T0 = T1; eos0 = eos1;
                r1 = r;  T1 = T; eos1 = eos;
            }

            iter++;
            return r;
        };
    };
    solver_resid resid(&HEOS, HEOS._p, value, other);
    
    std::string errstr;
    try{
        Brent(resid, Tmin, Tmax, DBL_EPSILON, 1e-12, 100, errstr);
        // Un-specify the phase of the fluid
        HEOS.unspecify_phase();
    }
    catch(...){
        // Un-specify the phase of the fluid
        HEOS.unspecify_phase();
        
        // Determine why you were out of range if you can
        // 
        CoolPropDbl eos0 = resid.eos0, eos1 = resid.eos1;
        std::string name = get_parameter_information(other,"short");
        std::string units = get_parameter_information(other,"units");
        if (eos1 > eos0 && value > eos1){
            throw ValueError(format("HSU_P_flash_singlephase_Brent could not find a solution because %s [%Lg %s] is above the maximum value of %0.12Lg %s", name.c_str(), value, units.c_str(), eos1, units.c_str()));
        }
        if (eos1 > eos0 && value < eos0){
            throw ValueError(format("HSU_P_flash_singlephase_Brent could not find a solution because %s [%Lg %s] is below the minimum value of %0.12Lg %s", name.c_str(), value, units.c_str(), eos0, units.c_str()));
        }
        throw;
    }
}

// P given and one of H, S, or U
void FlashRoutines::HSU_P_flash(HelmholtzEOSMixtureBackend &HEOS, parameters other)
{
    bool saturation_called = false;
    CoolPropDbl value;
    if (HEOS.imposed_phase_index != iphase_not_imposed)
    {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
    }
    else
    {
        // Find the phase, while updating all internal variables possible
        switch (other)
        {
            case iSmolar:
                value = HEOS.smolar(); break;
            case iHmolar:
                value = HEOS.hmolar(); break;
            case iUmolar:
                value = HEOS.umolar(); break;
            default:
                throw ValueError(format("Input for other [%s] is invalid", get_parameter_information(other, "long").c_str()));
        }
        if (HEOS.is_pure_or_pseudopure)
        {

            // Find the phase, while updating all internal variables possible
            HEOS.p_phase_determination_pure_or_pseudopure(other, value, saturation_called);
            
            if (HEOS.isHomogeneousPhase())
            {
                // Now we use the single-phase solver to find T,rho given P,Y using a 
                // bounded 1D solver by adjusting T and using given value of p
                CoolPropDbl Tmin, Tmax;
                switch(HEOS._phase)
                {
                    case iphase_gas:
                    {
                        Tmax = 1.5*HEOS.Tmax();
                        if (saturation_called){ Tmin = HEOS.SatV->T();}else{Tmin = HEOS._TVanc.pt();}
                        break;
                    }
                    case iphase_liquid:
                    {
                        if (saturation_called){ Tmax = HEOS.SatL->T();}else{Tmax = HEOS._TLanc.pt();}
                        
                        // Sometimes the minimum pressure for the melting line is a bit above the triple point pressure
                        if (HEOS.has_melting_line() && HEOS._p > HEOS.calc_melting_line(iP_min, -1, -1)){
                            Tmin = HEOS.calc_melting_line(iT, iP, HEOS._p)-1e-3;
                        }
                        else{
                            Tmin = HEOS.Tmin()-1e-3;
                        }
                        break;
                    }
                    case iphase_supercritical_liquid:
                    case iphase_supercritical_gas:
                    case iphase_supercritical:
                    {
                        Tmax = 1.5*HEOS.Tmax();
                        // Sometimes the minimum pressure for the melting line is a bit above the triple point pressure
                        if (HEOS.has_melting_line() && HEOS._p > HEOS.calc_melting_line(iP_min, -1, -1)){
                            Tmin = HEOS.calc_melting_line(iT, iP, HEOS._p)-1e-3;
                        }
                        else{
                            Tmin = HEOS.Tmin()-1e-3;
                        }
                        break;
                    }
                    default:
                    { throw ValueError(format("Not a valid homogeneous state")); }
                }
                try{
                    HSU_P_flash_singlephase_Brent(HEOS, other, value, Tmin, Tmax);
                }
                catch(std::exception &e){
                    throw ValueError(format("unable to solve 1phase PY flash with Tmin=%Lg, Tmax=%Lg due to error: %s",Tmin, Tmax, e.what()));
                }
                HEOS._Q = -1;
            }
        }
        else
        {
            std::cout << format("PHSU flash for mixture\n");
            if (HEOS.PhaseEnvelope.built){
                // Determine whether you are inside or outside
                SimpleState closest_state;
                std::size_t iclosest;
                std::cout << format("pre is inside\n");
                bool twophase = PhaseEnvelopeRoutines::is_inside(HEOS, iP, HEOS._p, other, value, iclosest, closest_state);
                std::cout << format("post is inside\n");
                
                std::string errstr;
                if (!twophase){
                    PY_singlephase_flash_resid resid(HEOS, HEOS._p, other, value);
                    // If that fails, try a bounded solver
                    CoolPropDbl rhomolar = Brent(resid, closest_state.T+10, 1000, DBL_EPSILON, 1e-10, 100, errstr);
                    HEOS.unspecify_phase();
                }
                else{
                    throw ValueError("two-phase solution for Y");
                }
                
            }
            else{
                throw ValueError("phase envelope must be built");
            }
        }
    }    
}
void FlashRoutines::DHSU_T_flash(HelmholtzEOSMixtureBackend &HEOS, parameters other)
{
    if (HEOS.imposed_phase_index != iphase_not_imposed)
    {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
    }
    else
    {
        if (HEOS.is_pure_or_pseudopure)
        {
            // Find the phase, while updating all internal variables possible
            switch (other)
            {
                case iDmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iDmolar, HEOS._rhomolar); break;
                case iSmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iSmolar, HEOS._smolar); break;
                case iHmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iHmolar, HEOS._hmolar); break;
                case iUmolar:
                    HEOS.T_phase_determination_pure_or_pseudopure(iUmolar, HEOS._umolar); break;
                default:
                    throw ValueError(format("Input is invalid"));
            }
        }
        else
        {
            HEOS._phase = iphase_gas;
            throw NotImplementedError("DHSU_T_flash does not support mixtures (yet)");
            // Find the phase, while updating all internal variables possible
        }
    }

    if (HEOS.isHomogeneousPhase() && !ValidNumber(HEOS._p))
    {
        switch (other)
        {
            case iDmolar:
                break;
            case iHmolar:
                HEOS._rhomolar = HEOS.solver_for_rho_given_T_oneof_HSU(HEOS._T, HEOS._hmolar, iHmolar); break;
            case iSmolar:
                HEOS._rhomolar = HEOS.solver_for_rho_given_T_oneof_HSU(HEOS._T, HEOS._smolar, iSmolar); break;
            case iUmolar:
                HEOS._rhomolar = HEOS.solver_for_rho_given_T_oneof_HSU(HEOS._T, HEOS._umolar, iUmolar); break;
            default:
                break;
        }
        HEOS.calc_pressure();
        HEOS._Q = -1;
    }
    // Update the state for conditions where the state was guessed
    
}
void FlashRoutines::HS_flash_twophase(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec, HS_flash_twophaseOptions &options)
{
    class Residual : public FuncWrapper1D
    {
        
    public:
        HelmholtzEOSMixtureBackend &HEOS;
        CoolPropDbl hmolar, smolar, Qs, Qh;
        Residual(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec) : HEOS(HEOS), hmolar(hmolar_spec), smolar(smolar_spec){};
        double call(double T){
            HEOS.update(QT_INPUTS, 0, T);
            HelmholtzEOSMixtureBackend &SatL = HEOS.get_SatL(),
                                       &SatV = HEOS.get_SatV();
            // Quality from entropy
            Qs = (smolar-SatL.smolar())/(SatV.smolar()-SatL.smolar());
            // Quality from enthalpy
            Qh = (hmolar-SatL.hmolar())/(SatV.hmolar()-SatL.hmolar());
            // Residual is the difference between the two
            return Qh-Qs;
        }
    } resid(HEOS, hmolar_spec, smolar_spec);
    
    std::string errstr;
    // Critical point for pure fluids, slightly different for pseudo-pure, very different for mixtures
    CoolPropDbl Tmax_sat = HEOS.calc_Tmax_sat() - 1e-13;
    
    // Check what the minimum limits for the equation of state are
    CoolPropDbl Tmin_satL, Tmin_satV, Tmin_sat;
    HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
    Tmin_sat = std::max(Tmin_satL, Tmin_satV) - 1e-13;
        
    Brent(resid, Tmin_sat, Tmax_sat-0.01, DBL_EPSILON, 1e-12, 20, errstr);
}
void FlashRoutines::HS_flash_singlephase(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl hmolar_spec, CoolPropDbl smolar_spec, HS_flash_singlephaseOptions &options)
{
    int iter = 0;
    double resid = 9e30, resid_old = 9e30;
    CoolProp::SimpleState reducing = HEOS.get_state("reducing");
    do{
        // Independent variables are T0 and rhomolar0, residuals are matching h and s
        Eigen::Vector2d r;
        Eigen::Matrix2d J;
        r(0) = HEOS.hmolar() - hmolar_spec;
        r(1) = HEOS.smolar() - smolar_spec;
        J(0,0) = HEOS.first_partial_deriv(iHmolar, iTau, iDelta);
        J(0,1) = HEOS.first_partial_deriv(iHmolar, iDelta, iTau);
        J(1,0) = HEOS.first_partial_deriv(iSmolar, iTau, iDelta);
        J(1,1) = HEOS.first_partial_deriv(iSmolar, iDelta, iTau);
        // Step in v obtained from Jv = -r
        Eigen::Vector2d v = J.colPivHouseholderQr().solve(-r);
        bool good_solution = false;
        double tau0 = HEOS.tau(), delta0 = HEOS.delta();
        // Calculate the old residual after the last step
        resid_old = sqrt(POW2(HEOS.hmolar() - hmolar_spec) + POW2(HEOS.smolar() - smolar_spec));
        for (double frac = 1.0; frac > 0.001; frac /= 2)
        {
            try{
                // Calculate new values
                double tau_new = tau0 + options.omega*frac*v(0);
                double delta_new = delta0 + options.omega*frac*v(1);
                double T_new = reducing.T/tau_new;
                double rhomolar_new = delta_new*reducing.rhomolar;
                // Update state with step
                HEOS.update(DmolarT_INPUTS, rhomolar_new, T_new);
                resid = sqrt(POW2(HEOS.hmolar() - hmolar_spec) + POW2(HEOS.smolar() - smolar_spec));
                if (resid > resid_old){
                    throw ValueError(format("residual not decreasing; frac: %g, resid: %g, resid_old: %g", frac, resid, resid_old));
                }
                good_solution = true;
                break;
            }
            catch(...){
                HEOS.clear();
                continue;
            }            
        }
        if (!good_solution){
            throw ValueError(format("Not able to get a solution" ));
        }
        iter++;
        if (iter > 50){
            throw ValueError(format("HS_flash_singlephase took too many iterations; residual is %g; prior was %g", resid, resid_old));
        }
    }
    while(std::abs(resid) > 1e-9);
}
void FlashRoutines::HS_flash_generate_TP_singlephase_guess(HelmholtzEOSMixtureBackend &HEOS, double &T, double &p)
{
    // Randomly obtain a starting value that is single-phase
    double logp = ((double)rand()/(double)RAND_MAX)*(log(HEOS.pmax())-log(HEOS.p_triple()))+log(HEOS.p_triple());
    T = ((double)rand()/(double)RAND_MAX)*(HEOS.Tmax()-HEOS.Ttriple())+HEOS.Ttriple();
    p = exp(logp);
}
void FlashRoutines::HS_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.imposed_phase_index != iphase_not_imposed)
    {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
    }
    else
    {
        enum solution_type_enum{not_specified = 0, single_phase_solution, two_phase_solution};
        solution_type_enum solution;
        
        shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> HEOS_copy(new CoolProp::HelmholtzEOSMixtureBackend(HEOS.components));
        
        // Find maxima states if needed
        // Cache the maximum enthalpy saturation state;
        HEOS.calc_hsat_max();
        // For weird fluids like the siloxanes, there can also be a maximum 
        // entropy along the vapor saturation line. Try to find it if it has one
        HEOS.calc_ssat_max();
        
        CoolProp::SimpleState crit = HEOS.get_state("reducing");
        CoolProp::SimpleState &tripleL = HEOS.components[0]->triple_liquid,
                              &tripleV = HEOS.components[0]->triple_vapor;
                              
        double first_maxima_in_saturation_entropy;
        if (HEOS.ssat_max.exists == SsatSimpleState::SSAT_MAX_DOES_EXIST){
            first_maxima_in_saturation_entropy = HEOS.ssat_max.smolar;
        }
        else{
            first_maxima_in_saturation_entropy = tripleV.smolar;
        }
        
        // Enthalpy at solid line for given entropy
        double hsolid = (tripleV.hmolar-tripleL.hmolar)/(tripleV.smolar-tripleL.smolar)*(HEOS.smolar()-tripleL.smolar) + tripleL.hmolar;
        // Part A - first check if HS is below triple line formed by connecting the triple point states
        // If so, it is solid, and not supported
        if (HEOS.hmolar() < hsolid-0.1){ // -0.1 is for a buffer
            throw ValueError(format("Enthalpy [%g J/mol] is below solid enthalpy [%g J/mol] for entropy [%g J/mol/K]", HEOS.hmolar(), hsolid-0.1, HEOS.smolar()));
        }
        /* Now check up to the first maxima in saturated vapor entropy.
         * For nicely behaved fluids, this means all the way up to the triple point vapor
         * For not-nicely behaved fluids with a local maxima on the saturated vapor entropy curve, 
         * the entropy at the local maxima in entropy
         * 
         * In this section, there can only be one solution for the saturation temperature, which is why the solver
         * is divided up in this way.
         */
        else if (HEOS.smolar() < first_maxima_in_saturation_entropy){
            double Q;
            if (HEOS.smolar() < crit.smolar){ 
                Q = 0; // Liquid part
            }
            else{
                Q = 1; // Vapor part
            }
            
            // Update the temporary instance with saturated entropy 
            HEOS_copy->update(QSmolar_INPUTS, Q, HEOS.smolar());
                
            // Check if above the saturation enthalpy for given entropy
            // If it is, the inputs are definitely single-phase.  We are done here
            if (HEOS.hmolar() > HEOS_copy->hmolar()){
                solution = single_phase_solution;
            }
            else{
                // C2: It is below hsat(ssat)
                // Either two-phase, or vapor (for funky saturation curves like the siloxanes)
                
                /* If enthalpy is between enthalpy at maxima in h and triple
                 * point enthalpy, search in enthalpy to find possible solutions
                 * that yield the correct enthalpy
                 * 
                 * There should only be one solution since we have already bypassed the local maxima in enthalpy
                 */
                
                HEOS_copy->update_HmolarQ_with_guessT(HEOS.hmolar(), 1, HEOS.hsat_max.T);
                
                if (HEOS.smolar() > HEOS_copy->smolar()){
                    solution = single_phase_solution;
                }
                else{
                    // C2a: It is below ssatV(hsatV) --> two-phase
                    solution = two_phase_solution;
                }
            }
        }
        // Part D - Check higher limit
        else if (HEOS.smolar() > tripleV.smolar){
            solution = single_phase_solution;
            HEOS_copy->update(PT_INPUTS, HEOS_copy->p_triple(), 0.5*HEOS_copy->Tmin() + 0.5*HEOS_copy->Tmax());
        }
        // Part E - HEOS.smolar() > crit.hmolar > tripleV.smolar
        else{
            // Now branch depending on the saturated vapor curve
            // If maximum 
            throw ValueError(format("partE HEOS.smolar() = %g tripleV.smolar = %g", HEOS.smolar(), tripleV.smolar));
        }
        
        switch (solution){
            case single_phase_solution:
            {
                // Fixing it to be gas is probably sufficient
                HEOS_copy->specify_phase(iphase_gas);
                HS_flash_singlephaseOptions options;
                options.omega = 1.0;
                try{
                    // Do the flash calcs starting from the guess value
                    HS_flash_singlephase(*HEOS_copy, HEOS.hmolar(), HEOS.smolar(), options);
                    // Copy the results
                    HEOS.update(DmolarT_INPUTS, HEOS_copy->rhomolar(), HEOS_copy->T());
                    break;
                }
                catch(...){
                    try{
                        // Trying again with another guessed value
                        HEOS_copy->update(DmolarT_INPUTS, HEOS.rhomolar_critical()*1.3, HEOS.Tmax());
                        HS_flash_singlephase(*HEOS_copy, HEOS.hmolar(), HEOS.smolar(), options);
                        // Copy the results
                        HEOS.update(DmolarT_INPUTS, HEOS_copy->rhomolar(), HEOS_copy->T());
                        break;
                    }
                    catch (...){
                        // Trying again with another guessed value
                        HEOS_copy->update(DmolarT_INPUTS, HEOS.rhomolar_critical(), 0.5*HEOS.Tmax() + 0.5*HEOS.T_critical());
                        HS_flash_singlephase(*HEOS_copy, HEOS.hmolar(), HEOS.smolar(), options);
                        // Copy the results
                        HEOS.update(DmolarT_INPUTS, HEOS_copy->rhomolar(), HEOS_copy->T());
                        break;
                    }
                }
            }
            case two_phase_solution:
            {
                HS_flash_twophaseOptions options;
                HS_flash_twophase(*HEOS_copy, HEOS.hmolar(), HEOS.smolar(), options);
                HEOS.update_internal(*HEOS_copy);
                break;
            }
            default:
                throw ValueError("solution not set");
        }
    }
}

#if defined(ENABLE_CATCH)

TEST_CASE("PD with T very large should yield error","[PDflash]")
{
	shared_ptr<HelmholtzEOSBackend> HEOS(new HelmholtzEOSBackend("R134a"));
	double Tc = HEOS->T_critical();
	HEOS->update(DmassT_INPUTS, 1.1, 1.5*Tc);
	CHECK_THROWS(HEOS->update(DmassP_INPUTS, 2, 5*HEOS->p()));
}
#endif

} /* namespace CoolProp */
