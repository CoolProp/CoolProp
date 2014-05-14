#include "VLERoutines.h"
#include "FlashRoutines.h"
#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

void FlashRoutines::PT_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    // Find the phase, while updating all internal variables possible
    HEOS.T_phase_determination_pure_or_pseudopure(iP, HEOS._p);

    if (!HEOS.isHomogeneousPhase())
    {
        throw ValueError("twophase not implemented yet");
    }
    else
    {
        // Find density
        HEOS._rhomolar = HEOS.solver_rho_Tp(HEOS._T, HEOS._p);
    }
}

void FlashRoutines::QT_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        if (!(HEOS.components[0]->pEOS->pseudo_pure))
        {
            // Set some imput options
            SaturationSolvers::saturation_T_pure_options options;
            options.omega = 1.0;
            options.use_guesses = false;
            // Actually call the solver
            SaturationSolvers::saturation_T_pure(&HEOS, HEOS._T, options);
            // Load the outputs
            HEOS._p = HEOS._Q*HEOS.SatV->p() + (1- HEOS._Q)*HEOS.SatL->p();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
        }
        else{
            // Pseudo-pure fluid
            long double rhoLanc, rhoVanc, rhoLsat, rhoVsat;
            long double psatLanc = HEOS.components[0]->ancillaries.pL.evaluate(HEOS._T); // These ancillaries are used explicitly
            long double psatVanc = HEOS.components[0]->ancillaries.pV.evaluate(HEOS._T); // These ancillaries are used explicitly
            try{
                rhoLanc = HEOS.components[0]->ancillaries.rhoL.evaluate(HEOS._T);
                rhoVanc = HEOS.components[0]->ancillaries.rhoV.evaluate(HEOS._T);

                if (!ValidNumber(rhoLanc) || !ValidNumber(rhoVanc))
                {
                    throw ValueError("pseudo-pure failed");
                }

                rhoLsat = HEOS.solver_rho_Tp(HEOS._T, psatLanc, rhoLanc);
                rhoVsat = HEOS.solver_rho_Tp(HEOS._T, psatVanc, rhoLanc);
                if (!ValidNumber(rhoLsat) || !ValidNumber(rhoVsat) || 
                     fabs(rhoLsat/rhoLanc-1) > 0.1 || fabs(rhoVanc/rhoVsat-1) > 0.1)
                {
                    throw ValueError("pseudo-pure failed");
                }
            }
            catch (std::exception &){
                // Near the critical point, the behavior is not very nice, so we will just use the ancillary near the critical point
                rhoLsat = rhoLanc;
                rhoVsat = rhoVanc;
            }
            HEOS._p = HEOS._Q*psatVanc + (1-HEOS._Q)*psatLanc;
            HEOS._rhomolar = 1/(HEOS._Q/rhoVsat + (1-HEOS._Q)/rhoLsat);
        }
    }
    else
    {
        // Set some imput options
        SaturationSolvers::mixture_VLE_IO options;
        options.sstype = SaturationSolvers::imposed_T;
        options.Nstep_max = 5;

        // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
        long double pguess = SaturationSolvers::saturation_preconditioner(&HEOS, HEOS._T, SaturationSolvers::imposed_T, HEOS.mole_fractions);

        // Use Wilson iteration to obtain updated guess for pressure
        pguess = SaturationSolvers::saturation_Wilson(&HEOS, HEOS._Q, HEOS._T, SaturationSolvers::imposed_T, HEOS.mole_fractions, pguess);
        
        // Actually call the successive substitution solver
        SaturationSolvers::successive_substitution(&HEOS, HEOS._Q, HEOS._T, pguess, HEOS.mole_fractions, HEOS.K, options);
    }
}
void FlashRoutines::PQ_flash(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        double pc = HEOS.components[0]->pEOS->reduce.p;
        double Tc = HEOS.components[0]->pEOS->reduce.T;
        double Tt = HEOS.components[0]->pEOS->Ttriple;
        double pt = HEOS.components[0]->pEOS->ptriple;
        double Tsat_guess = 1/(1/Tc-(1/Tt-1/Tc)/log(pc/pt)*log(HEOS._p/pc));

        if (HEOS.components[0]->pEOS->pseudo_pure){
            // It is a psedo-pure mixture
            throw NotImplementedError("PQ_flash not implemented for pseudo-pure fluids");
        }
        else{
            // It is a pure fluid

            // Set some imput options
            SaturationSolvers::saturation_PHSU_pure_options options;
            // Specified variable is pressure
            options.specified_variable = SaturationSolvers::saturation_PHSU_pure_options::IMPOSED_PV;
            // Use logarithm of delta as independent variables
            options.use_logdelta = false;
            // Actually call the solver
            SaturationSolvers::saturation_PHSU_pure(&HEOS, HEOS._p, options);
            
            // Load the outputs
            HEOS._p = HEOS._Q*HEOS.SatV->p() + (1 - HEOS._Q)*HEOS.SatL->p();
            HEOS._rhomolar = 1/(HEOS._Q/HEOS.SatV->rhomolar() + (1 - HEOS._Q)/HEOS.SatL->rhomolar());
            HEOS._T = HEOS.SatL->T();
        }
    }
    else
    {
        // Set some imput options
        SaturationSolvers::mixture_VLE_IO io;
        io.sstype = SaturationSolvers::imposed_p;
        io.Nstep_max = 20;

        // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
        long double Tguess = SaturationSolvers::saturation_preconditioner(&HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

        // Use Wilson iteration to obtain updated guess for temperature
        Tguess = SaturationSolvers::saturation_Wilson(&HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);
        
        // Actually call the successive substitution solver
        SaturationSolvers::successive_substitution(&HEOS, HEOS._Q, Tguess, HEOS._p, HEOS.mole_fractions, HEOS.K, io);
        
        PhaseEnvelope::PhaseEnvelope_GV ENV_GV;
        ENV_GV.build(&HEOS, HEOS.mole_fractions, HEOS.K, io);
    }
}
// D given and one of P,H,S,U
void FlashRoutines::PHSU_D_flash(HelmholtzEOSMixtureBackend &HEOS, int other)
{
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:
        
        HelmholtzEOSMixtureBackend *HEOS;
        long double r, eos, rhomolar, value, T;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend *HEOS, long double rhomolar, long double value, int other) : HEOS(HEOS), rhomolar(rhomolar), value(value), other(other){};
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

    if (HEOS.imposed_phase_index > -1) 
    {
        // Use the phase defined by the imposed phase
        HEOS._phase = HEOS.imposed_phase_index;
    }
    else
    {
        if (HEOS.is_pure_or_pseudopure)
        {
            CoolPropFluid * component = HEOS.components[0];
            HelmholtzEOSMixtureBackend *Sat;
            long double rhoLtriple = component->pEOS->rhoLtriple;
            long double rhoVtriple = component->pEOS->rhoVtriple;
            // Check if in the "normal" region
            if (HEOS._rhomolar >= rhoVtriple && HEOS._rhomolar <= rhoLtriple)
            {
                long double yL, yV, value, y_solid;
                long double TLtriple = component->pEOS->Ttriple; //TODO: separate TL and TV for ppure
                long double TVtriple = component->pEOS->Ttriple; //TODO: separate TL and TV for ppure
                
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
                SaturationSolvers::saturation_D_pure_options options;
                options.omega = 1;
                options.use_logdelta = false;
                if (HEOS._rhomolar > HEOS._crit.rhomolar)
                {
                    options.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOL;
                    SaturationSolvers::saturation_D_pure(&HEOS, HEOS._rhomolar, options);
                    // SatL and SatV have the saturation values
                    Sat = HEOS.SatL;
                }
                else
                {
                    options.imposed_rho = SaturationSolvers::saturation_D_pure_options::IMPOSED_RHOV;
                    SaturationSolvers::saturation_D_pure(&HEOS, HEOS._rhomolar, options);
                    // SatL and SatV have the saturation values
                    Sat = HEOS.SatV;
                }

                // If it is above, it is not two-phase and either liquid, vapor or supercritical
                if (value > Sat->keyed_output(other))
                {
                    solver_resid resid(&HEOS, HEOS._rhomolar, value, other);

                    HEOS._T = Brent(resid, Sat->keyed_output(iT), HEOS.Tmax(), DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
                }
                else
                {
                    throw NotImplementedError("Two-phase for PHSU_D_flash not supported yet");
                }

            }
            // Check if vapor/solid region below triple point vapor density
            else if (HEOS._rhomolar < component->pEOS->rhoVtriple)
            {
                long double y, value;
                long double TVtriple = component->pEOS->Ttriple; //TODO: separate TL and TV for ppure

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

                    HEOS._T = Brent(resid, TVtriple, HEOS.Tmax(), DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
                }
                else
                {
                    throw ValueError(format("D < DLtriple"));
                }

            }
            // Check in the liquid/solid region above the triple point density
            else 
            {
                long double y, value;
                long double TLtriple = component->pEOS->Ttriple; //TODO: separate TL and TV for ppure

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

                    HEOS._T = Brent(resid, TLtriple, HEOS.Tmax(), DBL_EPSILON, 1e-12, 100, errstring);
                    HEOS._Q = 10000;
                    HEOS.calc_pressure();
                }
                else
                {
                    throw ValueError(format("D < DLtriple"));
                }
            }
        }
        else 
            throw NotImplementedError("PHSU_D_flash not ready for mixtures");
    }
}
void FlashRoutines::HSU_P_flash(HelmholtzEOSMixtureBackend &HEOS, int other)
{
    throw NotImplementedError("HSU_P_flash Not implemented yet");
}
void FlashRoutines::DHSU_T_flash(HelmholtzEOSMixtureBackend &HEOS, int other)
{
    if (HEOS.imposed_phase_index > -1) 
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
}

} /* namespace CoolProp */