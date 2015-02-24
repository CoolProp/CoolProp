
#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "MatrixMath.h"
#include "MixtureDerivatives.h"

namespace CoolProp {
    
void SaturationSolvers::saturation_critical(HelmholtzEOSMixtureBackend &HEOS, parameters ykey, CoolPropDbl y){
    
    class inner_resid : public FuncWrapper1D{
        public:
        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl T, desired_p, rhomolar_liq, calc_p, rhomolar_crit;
        
        inner_resid(HelmholtzEOSMixtureBackend *HEOS, CoolPropDbl T, CoolPropDbl desired_p)
            : HEOS(HEOS), T(T), desired_p(desired_p){};
        double call(double rhomolar_liq){
            this->rhomolar_liq = rhomolar_liq;
            HEOS->SatL->update(DmolarT_INPUTS, rhomolar_liq, T);
            calc_p = HEOS->SatL->p();
            std::cout << format("inner p: %0.16Lg; res: %0.16Lg", calc_p, calc_p - desired_p) << std::endl;
            return calc_p - desired_p;
        }
    };
    
    // Define the outer residual to be driven to zero - this is the equality of 
    // Gibbs function for both co-existing phases
    class outer_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        parameters ykey;
        CoolPropDbl y;
        CoolPropDbl r, T, rhomolar_liq, rhomolar_vap, value, p, gL, gV, rhomolar_crit;
        int other;

        outer_resid(HelmholtzEOSMixtureBackend &HEOS, CoolProp::parameters ykey, CoolPropDbl y) 
            : HEOS(&HEOS), ykey(ykey), y(y){
                rhomolar_crit = HEOS.rhomolar_critical();
            };
        double call(double rhomolar_vap){
            this->y = y;
            
            // Calculate the other variable (T->p or p->T) for given vapor density
            if (ykey == iT){
                T = y;
                HEOS->SatV->update(DmolarT_INPUTS, rhomolar_vap, y);
                this->p = HEOS->SatV->p();
                std::cout << format("outer p: %0.16Lg",this->p) << std::endl;
                inner_resid inner(HEOS, T, p);
                std::string errstr2;
                rhomolar_liq = Brent(inner, rhomolar_crit*1.5, rhomolar_crit*(1+1e-8), LDBL_EPSILON, 1e-10, 100, errstr2);
            }
            HEOS->SatL->update(DmolarT_INPUTS, rhomolar_liq, T);
            HEOS->SatV->update(DmolarT_INPUTS, rhomolar_vap, T);
            
            // Calculate the Gibbs functions for liquid and vapor
            gL = HEOS->SatL->gibbsmolar();
            gV = HEOS->SatV->gibbsmolar();
            
            // Residual is difference in Gibbs function
            r = gL - gV;
            
            return this->p;
        };
    };
    outer_resid resid(HEOS, iT, y);
    
    double rhomolar_crit = HEOS.rhomolar_critical();
    
    std::string errstr;
    Brent(&resid, rhomolar_crit*(1-1e-8), rhomolar_crit*0.5, DBL_EPSILON, 1e-9, 20, errstr);
}

void SaturationSolvers::saturation_T_pure_1D_P(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_options &options)
{
    
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl r, T, rhomolar_liq, rhomolar_vap, value, p, gL, gV;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, CoolPropDbl rhomolar_liq_guess, CoolPropDbl rhomolar_vap_guess) 
            : HEOS(&HEOS), T(T), rhomolar_liq(rhomolar_liq_guess), rhomolar_vap(rhomolar_vap_guess){};
        double call(double p){
            this->p = p;
            // Recalculate the densities using the current guess values
            HEOS->SatL->update_TP_guessrho(T, p, rhomolar_liq);
            HEOS->SatV->update_TP_guessrho(T, p, rhomolar_vap);
            
            // Calculate the Gibbs functions for liquid and vapor
            gL = HEOS->SatL->gibbsmolar();
            gV = HEOS->SatV->gibbsmolar();
            
            // Residual is difference in Gibbs function
            r = gL - gV;
            
            return r;
        };
    };
    solver_resid resid(HEOS, T, options.rhoL, options.rhoV);
    
    if (!ValidNumber(options.p)){throw ValueError(format("options.p is not valid in saturation_T_pure_1D_P for T = %Lg",T));};
    if (!ValidNumber(options.rhoL)){throw ValueError(format("options.rhoL is not valid in saturation_T_pure_1D_P for T = %Lg",T));};
    if (!ValidNumber(options.rhoV)){throw ValueError(format("options.rhoV is not valid in saturation_T_pure_1D_P for T = %Lg",T));};
    
    std::string errstr;
    try{
        Secant(resid, options.p, options.p*1.1, 1e-10, 100, errstr);
    }
    catch(std::exception &){
        CoolPropDbl pmax = std::min(options.p*1.03, static_cast<CoolPropDbl>(HEOS.p_critical()+1e-6));
        CoolPropDbl pmin = std::max(options.p*0.97, static_cast<CoolPropDbl>(HEOS.p_triple()-1e-6));
        Brent(resid, pmin, pmax, LDBL_EPSILON, 1e-8, 100, errstr);
    }
}

void SaturationSolvers::saturation_P_pure_1D_T(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl p, saturation_PHSU_pure_options &options){
    
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        CoolPropDbl r, p, rhomolar_liq, rhomolar_vap, value, T, gL, gV;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl p, CoolPropDbl rhomolar_liq_guess, CoolPropDbl rhomolar_vap_guess) 
            : HEOS(&HEOS), p(p), rhomolar_liq(rhomolar_liq_guess), rhomolar_vap(rhomolar_vap_guess){};
        double call(double T){
            this->T = T;
            // Recalculate the densities using the current guess values
            HEOS->SatL->update_TP_guessrho(T, p, rhomolar_liq);
            HEOS->SatV->update_TP_guessrho(T, p, rhomolar_vap);
            
            // Calculate the Gibbs functions for liquid and vapor
            gL = HEOS->SatL->gibbsmolar();
            gV = HEOS->SatV->gibbsmolar();
            
            // Residual is difference in Gibbs function
            r = gL - gV;
            
            return r;
        };
    };
    solver_resid resid(HEOS, p, options.rhoL, options.rhoV);
    
    if (!ValidNumber(options.T)){throw ValueError("options.T is not valid in saturation_P_pure_1D_T");};
    if (!ValidNumber(options.rhoL)){throw ValueError("options.rhoL is not valid in saturation_P_pure_1D_T");};
    if (!ValidNumber(options.rhoV)){throw ValueError("options.rhoV is not valid in saturation_P_pure_1D_T");};
    
    std::string errstr;
    CoolPropDbl Tmax = std::min(options.T + 2, static_cast<CoolPropDbl>(HEOS.T_critical()-1e-6));
    CoolPropDbl Tmin = std::max(options.T - 2, static_cast<CoolPropDbl>(HEOS.Ttriple()+1e-6));
    Brent(resid, Tmin, Tmax, LDBL_EPSILON, 1e-11, 100, errstr);
}
    
void SaturationSolvers::saturation_PHSU_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl specified_value, saturation_PHSU_pure_options &options)
{
    /*
    This function is inspired by the method of Akasaka:

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from
    Helmholtz Energy Equations of State",
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */
    std::vector<CoolPropDbl> negativer(3,_HUGE), v;
    std::vector<std::vector<CoolPropDbl> > J(3, std::vector<CoolPropDbl>(3,_HUGE));

    HEOS.calc_reducing_state();
    const SimpleState & reduce = HEOS.get_reducing_state();
    CoolProp::SimpleState crit = HEOS.get_state("reducing");
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    CoolPropDbl T, rhoL, rhoV, pL, pV, hL, sL, hV, sV;
    CoolPropDbl deltaL=0, deltaV=0, tau=0, error;
    int iter=0, specified_parameter;

    // Use the density ancillary function as the starting point for the solver
    try
    {
        if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_PL || options.specified_variable == saturation_PHSU_pure_options::IMPOSED_PV)
        {
            // Invert liquid density ancillary to get temperature
            // TODO: fit inverse ancillaries too
            try{
                T = HEOS.get_components()[0]->ancillaries.pL.invert(specified_value);
            }
            catch(std::exception &)
            {
                throw ValueError("Unable to invert ancillary equation");
            }
        }
        else if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_HL)
        {
            CoolProp::SimpleState hs_anchor = HEOS.get_state("hs_anchor");
            // Ancillary is deltah = h - hs_anchor.h
            try{ T = HEOS.get_components()[0]->ancillaries.hL.invert(specified_value - hs_anchor.hmolar); }
            catch(std::exception &){
                throw ValueError("Unable to invert ancillary equation for hL");
            }
        }
        else if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_HV)
        {
            class Residual : public FuncWrapper1D
            {
                public:
                CoolPropFluid *component;
                double h;
                Residual(CoolPropFluid &component, double h){
                    this->component = &component;
                    this->h = h;
                }
                double call(double T){
                    CoolPropDbl h_liq = component->ancillaries.hL.evaluate(T) + component->pEOS->hs_anchor.hmolar;
                    return h_liq + component->ancillaries.hLV.evaluate(T) - h;
                };
            };
            Residual resid(*(HEOS.get_components()[0]), HEOS.hmolar());
            
            // Ancillary is deltah = h - hs_anchor.h
            std::string errstr;
            CoolPropDbl Tmin_satL, Tmin_satV;
            HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
            double Tmin = Tmin_satL;
            double Tmax = HEOS.calc_Tmax_sat();
            try{ T = Brent(resid, Tmin-3, Tmax + 1, DBL_EPSILON, 1e-10, 50, errstr); }
            catch(std::exception &){
                shared_ptr<HelmholtzEOSMixtureBackend> HEOS_copy(new HelmholtzEOSMixtureBackend(HEOS.get_components()));
                HEOS_copy->update(QT_INPUTS, 1, Tmin); double hTmin = HEOS_copy->hmolar();
                HEOS_copy->update(QT_INPUTS, 1, Tmax); double hTmax = HEOS_copy->hmolar();
                T = (Tmax-Tmin)/(hTmax-hTmin)*(HEOS.hmolar()-hTmin) + Tmin;
            }
        }
        else if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_SL)
        {
            CoolPropFluid &component = *(HEOS.get_components()[0]);
            CoolProp::SaturationAncillaryFunction &anc = component.ancillaries.sL;
            CoolProp::SimpleState hs_anchor = HEOS.get_state("hs_anchor");
            // If near the critical point, use a near critical guess value for T
            if (std::abs(HEOS.smolar() - crit.smolar) < std::abs(component.ancillaries.sL.get_max_abs_error()))
            {
                T = std::max(0.99*crit.T, crit.T-0.1);
            }
            else{
                CoolPropDbl Tmin, Tmax, Tmin_satV;
                HEOS.calc_Tmin_sat(Tmin, Tmin_satV);
                Tmax = HEOS.calc_Tmax_sat();
                // Ancillary is deltas = s - hs_anchor.s
                // First try a conventional call
                try{
                    T = anc.invert(specified_value - hs_anchor.smolar, Tmin, Tmax);
                }
                catch(std::exception &){
                    try{ 
                        T = anc.invert(specified_value - hs_anchor.smolar, Tmin - 3, Tmax + 3); 
                    }
                    catch(std::exception &){
                        double vmin = anc.evaluate(Tmin);
                        double vmax = anc.evaluate(Tmax);
                        if (std::abs(specified_value - hs_anchor.smolar) < std::abs(vmax)){
                            T = Tmax - 0.1;
                        }
                        else{
                            throw ValueError(format("Unable to invert ancillary equation for sL for value %Lg with Tminval %g and Tmaxval %g ", specified_value - hs_anchor.smolar, vmin, vmax));
                        }
                    }
                }
            }
        }
        else if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_SV)
        {
            CoolPropFluid &component = *(HEOS.get_components()[0]);
            CoolProp::SimpleState hs_anchor = HEOS.get_state("hs_anchor");
            class Residual : public FuncWrapper1D
            {
                public:
                CoolPropFluid *component;
                double s;
                Residual(CoolPropFluid &component, double s){
                    this->component = &component;
                    this->s = s;
                }
                double call(double T){
                    CoolPropDbl s_liq = component->ancillaries.sL.evaluate(T) + component->pEOS->hs_anchor.smolar;
                    CoolPropDbl resid = s_liq + component->ancillaries.sLV.evaluate(T) - s;
                    
                    return resid;
                };
            };
            Residual resid(component, HEOS.smolar());
            
            // Ancillary is deltas = s - hs_anchor.s
            std::string errstr;
            CoolPropDbl Tmin_satL, Tmin_satV;
            HEOS.calc_Tmin_sat(Tmin_satL, Tmin_satV);
            double Tmin = Tmin_satL;
            double Tmax = HEOS.calc_Tmax_sat();
            try{
                T = Brent(resid, Tmin-3, Tmax, DBL_EPSILON, 1e-10, 50, errstr);
            }
            catch(std::exception &){
                CoolPropDbl vmax = resid.call(Tmax);
                // If near the critical point, use a near critical guess value for T
                if (std::abs(specified_value - hs_anchor.smolar) < std::abs(vmax)){
                    T = std::max(0.99*crit.T, crit.T-0.1);
                }
                else{
                    shared_ptr<HelmholtzEOSMixtureBackend> HEOS_copy(new HelmholtzEOSMixtureBackend(HEOS.get_components()));
                    HEOS_copy->update(QT_INPUTS, 1, Tmin); double sTmin = HEOS_copy->smolar();
                    HEOS_copy->update(QT_INPUTS, 1, Tmax); double sTmax = HEOS_copy->smolar();
                    T = (Tmax-Tmin)/(sTmax-sTmin)*(HEOS.smolar()-sTmin) + Tmin;
                }
            }
        }
        else
        {
            throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }
        // If T from the ancillaries is above the critical temp, this will cause failure 
        // in ancillaries for rhoV and rhoL, decrease if needed
        T = std::min(T, static_cast<CoolPropDbl>(HEOS.T_critical()-0.1));

        // Evaluate densities from the ancillary equations
        rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T);
        rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T);

        // Apply a single step of Newton's method to improve guess value for liquid
        // based on the error between the gas pressure (which is usually very close already)
        // and the liquid pressure, which can sometimes (especially at low pressure),
        // be way off, and often times negative
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        double rhoL_updated = rhoL -(SatL->p()-SatV->p())/SatL->first_partial_deriv(iP, iDmolar, iT);
        
        // Accept the update if the liquid density is greater than the vapor density
        if (rhoL_updated > rhoV){ rhoL = rhoL_updated; }
        
        // Update the state again with the better guess for the liquid density
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);

        deltaL = rhoL/reduce.rhomolar;
        deltaV = rhoV/reduce.rhomolar;
        tau = reduce.T/T;
    }
    catch(NotImplementedError &)
    {
        throw;
    }

    do{
        /*if (get_debug_level()>8){
            std::cout << format("%s:%d: right before the derivs with deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
        }*/

        pL = SatL->p();
        hL = SatL->hmolar();
        sL = SatL->smolar();
        pV = SatV->p();
        hV = SatV->hmolar();
        sV = SatV->smolar();

        // These derivatives are needed for both cases
        CoolPropDbl alpharL = SatL->alphar();
        CoolPropDbl alpharV = SatV->alphar();
        CoolPropDbl dalphar_dtauL = SatL->dalphar_dTau();
        CoolPropDbl dalphar_dtauV = SatV->dalphar_dTau();
        CoolPropDbl d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        CoolPropDbl d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        CoolPropDbl dalphar_ddeltaL = SatL->dalphar_dDelta();
        CoolPropDbl dalphar_ddeltaV = SatV->dalphar_dDelta();
        CoolPropDbl d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        CoolPropDbl d2alphar_ddelta2V = SatV->d2alphar_dDelta2();

        // -r_1 (equate the pressures)
        negativer[0] = -(deltaV*(1+deltaV*dalphar_ddeltaV)-deltaL*(1+deltaL*dalphar_ddeltaL));
        // -r_2 (equate the gibbs energy)
        negativer[1] = -(deltaV*dalphar_ddeltaV+alpharV+log(deltaV)-deltaL*dalphar_ddeltaL-alpharL-log(deltaL));
        switch (options.specified_variable){
            case saturation_PHSU_pure_options::IMPOSED_PL:
                // -r_3 (equate calculated pressure and specified liquid pressure)
                negativer[2] = -(pL/specified_value - 1); break;
            case saturation_PHSU_pure_options::IMPOSED_PV:
                // -r_3 (equate calculated pressure and specified vapor pressure)
                negativer[2] = -(pV/specified_value - 1); break;
            case saturation_PHSU_pure_options::IMPOSED_HL:
                // -r_3 (equate calculated liquid enthalpy and specified liquid enthalpy)
                negativer[2] = -(hL - specified_value); break;
            case saturation_PHSU_pure_options::IMPOSED_HV:
                // -r_3 (equate calculated vapor enthalpy and specified vapor enthalpy)
                negativer[2] = -(hV - specified_value); break;
            case saturation_PHSU_pure_options::IMPOSED_SL:
                // -r_3 (equate calculated liquid entropy and specified liquid entropy)
                negativer[2] = -(sL - specified_value); break;
            case saturation_PHSU_pure_options::IMPOSED_SV:
                // -r_3 (equate calculated vapor entropy and specified vapor entropy)
                negativer[2] = -(sV - specified_value); break;
            default:
                throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }

        // dr1_dtau
        J[0][0] = pow(deltaV,2)*d2alphar_ddelta_dtauV-pow(deltaL,2)*d2alphar_ddelta_dtauL;
        // dr2_dtau
        J[1][0] = deltaV*d2alphar_ddelta_dtauV+dalphar_dtauV-deltaL*d2alphar_ddelta_dtauL-dalphar_dtauL;

        if (options.use_logdelta){
            // dr_1/d_log(delta'')
            J[0][1] = -deltaL-2*pow(deltaL,2)*dalphar_ddeltaL-pow(deltaL,3)*d2alphar_ddelta2L;
            // dr_2/d_log(delta'')
            J[1][1] = -pow(deltaL,2)*d2alphar_ddelta2L-2*deltaL*dalphar_ddeltaL-1;
        }
        else{
            // dr_1/ddelta''
            J[0][1] = -1-2*deltaL*dalphar_ddeltaL-pow(deltaL,2)*d2alphar_ddelta2L;
            // dr_2/ddelta''
            J[1][1] = -1/deltaL-2*dalphar_ddeltaL-deltaL*d2alphar_ddelta2L;
        }

        if (options.use_logdelta){
            // dr_1/d_log(delta'')
            J[0][2] = deltaV+2*pow(deltaV,2)*dalphar_ddeltaV+pow(deltaV,3)*d2alphar_ddelta2V;
            // dr_2/d_log(delta'')
            J[1][2] = 1+2*deltaV*dalphar_ddeltaV+1+pow(deltaV,2)*d2alphar_ddelta2V;
        }
        else{
            // dr_1/ddelta''
            J[0][2] = 1+2*deltaV*dalphar_ddeltaV+pow(deltaV,2)*d2alphar_ddelta2V;
            // dr_2/ddelta''
            J[1][2] = deltaV*d2alphar_ddelta2V+2*dalphar_ddeltaV+1/deltaV;
        }

        // Derivatives of the specification equation
        switch (options.specified_variable){
            case saturation_PHSU_pure_options::IMPOSED_PL:
                // dr_3/dtau
                J[2][0] = SatL->first_partial_deriv(iP,iTau,iDelta)/specified_value;
                if (options.use_logdelta){
                    // dr_3/d(log(delta'))
                    J[2][1] = deltaL*SatL->first_partial_deriv(iP,iDelta,iTau)/specified_value;
                }
                else{
                    // dr_3/ddelta'
                    J[2][1] = SatL->first_partial_deriv(iP,iDelta,iTau)/specified_value;
                }
                // dr_3/ddelta'' (liquid pressure not a function of vapor density)
                J[2][2] = 0;
                specified_parameter = CoolProp::iP;
                break;
            case saturation_PHSU_pure_options::IMPOSED_PV:
                // dr_3/dtau
                J[2][0] = SatV->first_partial_deriv(iP,iTau,iDelta)/specified_value;
                // dr_3/ddelta' (vapor pressure not a function of liquid density)
                J[2][1] = 0;
                if (options.use_logdelta){
                    // dr_3/d(log(delta'')
                    J[2][2] = deltaV*SatV->first_partial_deriv(iP,iDelta,iTau)/specified_value;
                }
                else{
                    // dr_3/ddelta''
                    J[2][2] = SatV->first_partial_deriv(iP,iDelta,iTau)/specified_value;
                }
                specified_parameter = CoolProp::iP;
                break;
            case saturation_PHSU_pure_options::IMPOSED_HL:
                // dr_3/dtau
                J[2][0] = SatL->first_partial_deriv(iHmolar,iTau,iDelta);
                // dr_3/ddelta'
                J[2][1] = SatL->first_partial_deriv(iHmolar,iDelta,iTau);
                if (options.use_logdelta){ J[2][1]*=deltaL;}
                // dr_3/ddelta''
                J[2][2] = 0; //(liquid enthalpy not a function of vapor density)
                specified_parameter = CoolProp::iHmolar;
                break;
            case saturation_PHSU_pure_options::IMPOSED_HV:
                // dr_3/dtau
                J[2][0] = SatV->first_partial_deriv(iHmolar,iTau,iDelta);
                // dr_3/ddelta'
                J[2][1] = 0; //(vapor enthalpy not a function of liquid density)
                // dr_3/ddelta''
                J[2][2] = SatV->first_partial_deriv(iHmolar,iDelta,iTau);
                if (options.use_logdelta){ J[2][2]*=deltaV;}
                specified_parameter = CoolProp::iHmolar;
                break;
            case saturation_PHSU_pure_options::IMPOSED_SL:
                // dr_3/dtau
                J[2][0] = SatL->first_partial_deriv(iSmolar,iTau,iDelta);
                // dr_3/ddelta'
                J[2][1] = SatL->first_partial_deriv(iSmolar,iDelta,iTau);
                if (options.use_logdelta){ J[2][1] *= deltaL; }
                // dr_3/ddelta''
                J[2][2] = 0; //(liquid entropy not a function of vapor density)
                specified_parameter = CoolProp::iSmolar;
                break;
            case saturation_PHSU_pure_options::IMPOSED_SV:
                // dr_3/dtau
                J[2][0] = SatV->first_partial_deriv(iSmolar,iTau,iDelta);
                // dr_3/ddelta'
                J[2][1] = 0; //(vapor enthalpy not a function of liquid density)
                // dr_3/ddelta''
                J[2][2] = SatV->first_partial_deriv(iSmolar,iDelta,iTau);
                if (options.use_logdelta){ J[2][2]*=deltaV;}
                specified_parameter = CoolProp::iSmolar;
                break;
            default:
                throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }

        v = linsolve(J, negativer);
        
        // Conditions for an acceptable step are:
        // a) tau > 1
        // b) rhoL > rhoV or deltaL > deltaV
        double tau0 = tau, deltaL0 = deltaL, deltaV0 = deltaV;
        for (double omega_local = 1.0; omega_local > 0.1; omega_local /= 1.1)
        {
            tau = tau0 + omega_local*options.omega*v[0];
            if (options.use_logdelta){
                deltaL = exp(log(deltaL0)+omega_local*options.omega*v[1]);
                deltaV = exp(log(deltaV0)+omega_local*options.omega*v[2]);
            }
            else{
                deltaL = deltaL0 + omega_local*options.omega*v[1];
                deltaV = deltaV0 + omega_local*options.omega*v[2];
            }
            if (tau > 1 && deltaL > deltaV){
                break;
            }
        }

        rhoL = deltaL*reduce.rhomolar;
        rhoV = deltaV*reduce.rhomolar;
        T = reduce.T/tau;
        
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);

        error = sqrt(pow(negativer[0], 2)+pow(negativer[1], 2)+pow(negativer[2], 2));
        iter++;
        if (T < 0)
        {
            throw SolutionError(format("saturation_PHSU_pure solver T < 0"));
        }
        if (iter > 50){
            // Set values back into the options structure for use in next solver
            options.rhoL = rhoL; options.rhoV = rhoV; options.T = T;
            // Error out
            std::string info = get_parameter_information(specified_parameter, "short");
            throw SolutionError(format("saturation_PHSU_pure solver did not converge after 50 iterations for %s=%Lg current error is %Lg", info.c_str(), specified_value, error));
        }
    }
    while (error > 1e-9);
}
void SaturationSolvers::saturation_D_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl rhomolar, saturation_D_pure_options &options)
{
    /*
    This function is inspired by the method of Akasaka:

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from
    Helmholtz Energy Equations of State",
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */
    std::vector<CoolPropDbl> r(2,_HUGE), v;
    std::vector<std::vector<CoolPropDbl> > J(2, std::vector<CoolPropDbl>(2,_HUGE));

    HEOS.calc_reducing_state();
    const SimpleState & reduce = HEOS.get_reducing_state();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    CoolPropDbl T, rhoL,rhoV;
    CoolPropDbl deltaL=0, deltaV=0, tau=0, error, p_error;
    int iter=0;

    // Use the density ancillary function as the starting point for the solver
    try
    {
        if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOL)
        {
            // Invert liquid density ancillary to get temperature
            // TODO: fit inverse ancillaries too
            T = HEOS.get_components()[0]->ancillaries.rhoL.invert(rhomolar);
            rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T);
            rhoL = rhomolar;
        }
        else if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOV)
        {
            // Invert vapor density ancillary to get temperature
            // TODO: fit inverse ancillaries too
            T = HEOS.get_components()[0]->ancillaries.rhoV.invert(rhomolar);
            rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T);
            rhoV = rhomolar;
        }
        else
        {
            throw ValueError(format("imposed rho to saturation_D_pure [%d%] is invalid",options.imposed_rho));
        }

        deltaL = rhoL/reduce.rhomolar;
        deltaV = rhoV/reduce.rhomolar;
        tau = reduce.T/T;
    }
    catch(NotImplementedError &e)
    {
        throw e;
    }

    do{
        /*if (get_debug_level()>8){
            std::cout << format("%s:%d: right before the derivs with deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
        }*/

        // Calculate once to save on calls to EOS
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);

        CoolPropDbl pL = SatL->p();
        CoolPropDbl pV = SatV->p();

        // These derivatives are needed for both cases
        CoolPropDbl dalphar_dtauL = SatL->dalphar_dTau();
        CoolPropDbl dalphar_dtauV = SatV->dalphar_dTau();
        CoolPropDbl d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        CoolPropDbl d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        CoolPropDbl alpharL = SatL->alphar();
        CoolPropDbl alpharV = SatV->alphar();
        CoolPropDbl dalphar_ddeltaL = SatL->dalphar_dDelta();
        CoolPropDbl dalphar_ddeltaV = SatV->dalphar_dDelta();

        // -r_1
        r[0] = -(deltaV*(1+deltaV*dalphar_ddeltaV)-deltaL*(1+deltaL*dalphar_ddeltaL));
        // -r_2
        r[1] =  -(deltaV*dalphar_ddeltaV+alpharV+log(deltaV)-deltaL*dalphar_ddeltaL-alpharL-log(deltaL));

        // dr1_dtau
        J[0][0] = pow(deltaV,2)*d2alphar_ddelta_dtauV-pow(deltaL,2)*d2alphar_ddelta_dtauL;
        // dr2_dtau
        J[1][0] = deltaV*d2alphar_ddelta_dtauV+dalphar_dtauV-deltaL*d2alphar_ddelta_dtauL-dalphar_dtauL;

        if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOL)
        {
            CoolPropDbl d2alphar_ddelta2V = SatV->d2alphar_dDelta2();
            if (options.use_logdelta)
            {
                J[0][1] = deltaV+2*pow(deltaV,2)*dalphar_ddeltaV+pow(deltaV,3)*d2alphar_ddelta2V;
                J[1][1] = pow(deltaV,2)*d2alphar_ddelta2V+2*deltaV*dalphar_ddeltaV+1;
            }
            else
            {
                J[0][1] = 1+2*deltaV*dalphar_ddeltaV+pow(deltaV,2)*d2alphar_ddelta2V;
                J[1][1] = deltaV*d2alphar_ddelta2V+2*dalphar_ddeltaV+1/deltaV;
            }
        }
        else if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOV)
        {
            CoolPropDbl d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
            if (options.use_logdelta)
            {
                J[0][1] = -deltaL-2*pow(deltaL,2)*dalphar_ddeltaL-pow(deltaL,3)*d2alphar_ddelta2L;
                J[1][1] = -pow(deltaL,2)*d2alphar_ddelta2L-2*deltaL*dalphar_ddeltaL-1;
            }
            else
            {
                J[0][1] = -1-2*deltaL*dalphar_ddeltaL-pow(deltaL,2)*d2alphar_ddelta2L;
                J[1][1] = -deltaL*d2alphar_ddelta2L-2*dalphar_ddeltaL-1/deltaL;
            }
        }

        //double DET = J[0][0]*J[1][1]-J[0][1]*J[1][0];

        v = linsolve(J, r);

        tau += options.omega*v[0];

        if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOL)
        {
            if (options.use_logdelta)
                deltaV = exp(log(deltaV)+options.omega*v[1]);
            else
                deltaV += v[1];
        }
        else
        {
            if (options.use_logdelta)
                deltaL = exp(log(deltaL)+options.omega*v[1]);
            else
                deltaL += v[1];
        }

        rhoL = deltaL*reduce.rhomolar;
        rhoV = deltaV*reduce.rhomolar;
        T = reduce.T/tau;
        
        p_error = (pL-pV)/pL;

        error = sqrt(pow(r[0], 2)+pow(r[1], 2));
        iter++;
        if (T < 0)
        {
            throw SolutionError(format("saturation_D_pure solver T < 0"));
        }
        if (iter > 200){
            throw SolutionError(format("saturation_D_pure solver did not converge after 100 iterations with rho: %g mol/m^3",rhomolar));
        }
    }
    while (error > 1e-9);
    CoolPropDbl p_error_limit = 1e-3;
    if (std::abs(p_error) > p_error_limit){
        throw SolutionError(format("saturation_D_pure solver abs error on p [%Lg] > limit [%Lg]", p_error, p_error_limit));
    }
}
void SaturationSolvers::saturation_T_pure(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_options &options)
{
    // Set some imput options
    SaturationSolvers::saturation_T_pure_Akasaka_options _options;
    _options.omega = 1.0;
    _options.use_guesses = false;
    try{
        // Actually call the solver
        SaturationSolvers::saturation_T_pure_Maxwell(HEOS, T, _options);
    }
    catch(std::exception &){
        try{
            // Actually call the solver
            SaturationSolvers::saturation_T_pure_Akasaka(HEOS, T, _options);
        }
        catch(std::exception &){
            // If there was an error, store values for use in later solvers
            options.pL = _options.pL;
            options.pV = _options.pV;
            options.rhoL = _options.rhoL;
            options.rhoV = _options.rhoV;
            options.p = _options.pL;
            SaturationSolvers::saturation_T_pure_1D_P(HEOS, T, options);
        }
    }
}
void SaturationSolvers::saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options &options)
{
    // Start with the method of Akasaka

    /*
    This function implements the method of Akasaka

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from
    Helmholtz Energy Equations of State",
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */

    HEOS.calc_reducing_state();
    const SimpleState & reduce = HEOS.get_reducing_state();
    CoolPropDbl R_u = HEOS.gas_constant();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    CoolPropDbl rhoL = _HUGE, rhoV = _HUGE,JL,JV,KL,KV,dJL,dJV,dKL,dKV;
    CoolPropDbl DELTA, deltaL=0, deltaV=0, error, PL, PV, stepL, stepV;
    int iter=0;
    
    try
    {
        if (options.use_guesses)
        {
            // Use the guesses provided in the options structure
            rhoL = options.rhoL;
            rhoV = options.rhoV;
        }
        else
        {
            // Use the density ancillary function as the starting point for the solver
            
            // If very close to the critical temp, evaluate the ancillaries for a slightly lower temperature
            if (T > 0.99*HEOS.get_reducing_state().T){
                rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T-0.1);
                rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T-0.1);
            }
            else{
                rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T);
                rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T);
                
                // Apply a single step of Newton's method to improve guess value for liquid
                // based on the error between the gas pressure (which is usually very close already)
                // and the liquid pressure, which can sometimes (especially at low pressure),
                // be way off, and often times negative
                SatL->update(DmolarT_INPUTS, rhoL, T);
                SatV->update(DmolarT_INPUTS, rhoV, T);
                
                // Update the guess for liquid density using density solver with vapor pressure 
                // and liquid density guess from ancillaries
                HEOS.specify_phase(iphase_liquid);
                rhoL = HEOS.solver_rho_Tp(T, SatV->p(), rhoL);
                HEOS.unspecify_phase();
            }
        }

        deltaL = rhoL/reduce.rhomolar;
        deltaV = rhoV/reduce.rhomolar;
    }
    catch(NotImplementedError &)
    {
        /*double Tc = crit.T;
        double pc = crit.p.Pa;
        double w = 6.67228479e-09*Tc*Tc*Tc-7.20464352e-06*Tc*Tc+3.16947758e-03*Tc-2.88760012e-01;
        double q = -6.08930221451*w -5.42477887222;
        double pt = exp(q*(Tc/T-1))*pc;*/

        //double rhoL = density_Tp_Soave(T, pt, 0), rhoV = density_Tp_Soave(T, pt, 1);

        //deltaL = rhoL/reduce.rhomolar;
        //deltaV = rhoV/reduce.rhomolar;
        //tau = reduce.T/T;
    }
    //if (get_debug_level()>5){
    //        std::cout << format("%s:%d: Akasaka guess values deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
    //    }

    do{
        /*if (get_debug_level()>8){
            std::cout << format("%s:%d: right before the derivs with deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
        }*/

        // Calculate once to save on calls to EOS
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        CoolPropDbl alpharL = SatL->alphar();
        CoolPropDbl alpharV = SatV->alphar();
        CoolPropDbl dalphar_ddeltaL = SatL->dalphar_dDelta();
        CoolPropDbl dalphar_ddeltaV = SatV->dalphar_dDelta();
        CoolPropDbl d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        CoolPropDbl d2alphar_ddelta2V = SatV->d2alphar_dDelta2();

        JL = deltaL * (1 + deltaL*dalphar_ddeltaL);
        JV = deltaV * (1 + deltaV*dalphar_ddeltaV);
        KL = deltaL*dalphar_ddeltaL + alpharL + log(deltaL);
        KV = deltaV*dalphar_ddeltaV + alpharV + log(deltaV);

        PL = R_u*reduce.rhomolar*T*JL;
        PV = R_u*reduce.rhomolar*T*JV;

        // At low pressure, the magnitude of d2alphar_ddelta2L and d2alphar_ddelta2V are enormous, truncation problems arise for all the partials
        dJL = 1 + 2*deltaL*dalphar_ddeltaL + deltaL*deltaL*d2alphar_ddelta2L;
        dJV = 1 + 2*deltaV*dalphar_ddeltaV + deltaV*deltaV*d2alphar_ddelta2V;
        dKL = 2*dalphar_ddeltaL + deltaL*d2alphar_ddelta2L + 1/deltaL;
        dKV = 2*dalphar_ddeltaV + deltaV*d2alphar_ddelta2V + 1/deltaV;

        DELTA = dJV*dKL-dJL*dKV;

        error = sqrt((KL-KV)*(KL-KV)+(JL-JV)*(JL-JV));

        //  Get the predicted step
        stepL = options.omega/DELTA*( (KV-KL)*dJV-(JV-JL)*dKV);
        stepV = options.omega/DELTA*( (KV-KL)*dJL-(JV-JL)*dKL);

        CoolPropDbl deltaL0 = deltaL, deltaV0 = deltaV;
        // Conditions for an acceptable step are:
        // a) rhoL > rhoV or deltaL > deltaV
        for (double omega_local = 1.0; omega_local > 0.1; omega_local /= 1.1)
        {
            deltaL = deltaL0 + omega_local*stepL;
            deltaV = deltaV0 + omega_local*stepV;

            if (deltaL > 1 && deltaV < 1 && deltaV > 0){
                break;
            }
        }
        
        rhoL = deltaL*reduce.rhomolar;
        rhoV = deltaV*reduce.rhomolar;
        iter++;
        if (iter > 100){
            throw SolutionError(format("Akasaka solver did not converge after 100 iterations"));
        }
    }
    while (error > 1e-10 && std::abs(stepL) > 10*DBL_EPSILON*std::abs(stepL) && std::abs(stepV) > 10*DBL_EPSILON*std::abs(stepV));
    
    CoolPropDbl p_error_limit = 1e-3;
    CoolPropDbl p_error = (PL - PV)/PL;
    if (std::abs(p_error) > p_error_limit){
        options.pL = PL;
        options.pV = PV;
        options.rhoL = rhoL;
        options.rhoV = rhoV;
        throw SolutionError(format("saturation_T_pure_Akasaka solver abs error on p [%g] > limit [%g]", std::abs(p_error), p_error_limit));
    }
}

CoolPropDbl sign(CoolPropDbl x)
{
    if (x > 0){
        return 1;
    }
    else{
        return -1;
    }
}

void SaturationSolvers::saturation_T_pure_Maxwell(HelmholtzEOSMixtureBackend &HEOS, CoolPropDbl T, saturation_T_pure_Akasaka_options &options)
{

    /*
    This function implements the method of 

    Ancillary equations are used to get a sensible starting point
    */

    HEOS.calc_reducing_state();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;
    CoolProp::SimpleState &crit = HEOS.get_components()[0]->crit;
    CoolPropDbl rhoL = _HUGE, rhoV = _HUGE, error = 999, DeltavL, DeltavV, pL, pV, p;
    int iter=0, small_step_count = 0;
    
    try
    {
        if (options.use_guesses)
        {
            // Use the guesses provided in the options structure
            rhoL = options.rhoL;
            rhoV = options.rhoV;
        }
        else
        {
            // Use the density ancillary function as the starting point for the solver
            
            // If very close to the critical temp, evaluate the ancillaries for a slightly lower temperature
            if (T > 0.9999*HEOS.get_reducing_state().T){
                rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T-0.1);
                rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T-0.1);
            }
            else{
                rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T);
                rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T);
                p = HEOS.get_components()[0]->ancillaries.pV.evaluate(T);
                
                CoolProp::SimpleState &tripleL = HEOS.get_components()[0]->triple_liquid;
                CoolProp::SimpleState &tripleV = HEOS.get_components()[0]->triple_vapor;
                
                // If the guesses are terrible, apply a simple correction
                if (rhoL < crit.rhomolar*0.8 || rhoL > tripleL.rhomolar*1.2 || 
                    rhoV > crit.rhomolar*1.2 || rhoV < tripleV.rhomolar*0.8)
                {
                    // Lets assume that liquid density is more or less linear with T
                    rhoL = (crit.rhomolar - tripleL.rhomolar)/(crit.T - tripleL.T)*(T-tripleL.T)+tripleL.rhomolar;
                    // Then we calculate pressure from this density
                    SatL->update_DmolarT_direct(rhoL, T);
                    // Then we assume vapor to be ideal gas
                    rhoV = SatL->p()/(SatL->gas_constant()*T);
                    // Update the vapor state
                    SatV->update_DmolarT_direct(rhoV, T);
                }
                else{
                    SatL->update_DmolarT_direct(rhoL, T);
                    SatV->update_DmolarT_direct(rhoV, T);
                }
                if (get_debug_level() > 0){ std::cout << format("[Maxwell] ancillaries T: %0.16Lg rhoL: %0.16Lg rhoV: %0.16Lg pL: %g pV: %g\n", T, rhoL, rhoV, SatL->p(), SatV->p());}
                
                // Update the guess for liquid density using density solver with vapor pressure 
                // and liquid density guess from ancillaries, but only if the pressures are not 
                // close to each other
                if (std::abs((SatL->p()-p)/p) > 0.1){
                    for (int iii = 0; iii < 6; ++iii){
                        // Use Halley's method to update the liquid density (http://en.wikipedia.org/wiki/Halley%27s_method)
                        CoolPropDbl f = SatL->p()-SatV->p();
                        CoolPropDbl dpdrho = SatL->first_partial_deriv(iP,iDmolar,iT);
                        CoolPropDbl d2pdrho2 = SatL->second_partial_deriv(iP,iDmolar,iT,iDmolar,iT);
                        CoolPropDbl deltarhoLHalley = -(2*f*dpdrho)/(2*POW2(dpdrho)-f*d2pdrho2);
                        rhoL += deltarhoLHalley;
                        if (std::abs(deltarhoLHalley/rhoL) < DBL_EPSILON){
                            break;
                        }
                        SatL->update_DmolarT_direct(rhoL, T);
                    }
                }

                SatL->update_DmolarT_direct(rhoL, T);
                SatV->update_DmolarT_direct(rhoV, T);
                    
                // Update the guess for vapor density using density solver with vapor pressure 
                // and density guess from ancillaries, but only if the pressures are not 
                // close to each other
                if (std::abs((SatV->p()-p)/p) > 0.1){
                    HEOS.specify_phase(iphase_gas);
                    rhoV = SatV->solver_rho_Tp(T, p, rhoV);
                    HEOS.unspecify_phase();
                }
            }
        }
    }
    catch(NotImplementedError &)
    {}
    
    if (rhoL < crit.rhomolar){
        rhoL = 1.01*crit.rhomolar;
    }
    if (rhoV > crit.rhomolar){
        rhoV = 0.99*crit.rhomolar;
    }

    SatL->update_DmolarT_direct(rhoL, T);
    SatV->update_DmolarT_direct(rhoV, T);
    if (get_debug_level() > 0){ std::cout << format("[Maxwell] starting T: %0.16Lg rhoL: %Lg rhoV: %Lg pL: %Lg pV: %g\n", T, rhoL, rhoV, SatL->p(), SatV->p());}
    do{
        pL = SatL->p(); pV = SatV->p();
        CoolPropDbl vL = 1/SatL->rhomolar(), vV = 1/SatV->rhomolar();
        // Get alpha, the pressure derivative with volume at constant T
        // Given by (dp/drho|T)*drhodv
        CoolPropDbl alphaL = SatL->first_partial_deriv(iP, iDmolar, iT)*(-POW2(SatL->rhomolar()));
        CoolPropDbl alphaV = SatV->first_partial_deriv(iP, iDmolar, iT)*(-POW2(SatV->rhomolar()));
        
        // Total helmholtz energy for liquid and vapor
        CoolPropDbl RT = SatL->gas_constant()*T;
        CoolPropDbl helmholtzL = (SatL->calc_alpha0() + SatL->calc_alphar())*RT;
        CoolPropDbl helmholtzV = (SatV->calc_alpha0() + SatV->calc_alphar())*RT;
        
        // Calculate the mean pressure
        CoolPropDbl pM = (helmholtzL - helmholtzV)/(vV - vL);
        
        // Coefficients for the quadratic in the step
        CoolPropDbl A = 0.5*alphaL*(alphaL - alphaV);
        CoolPropDbl B = alphaL*(pL - pV - alphaV*(vL - vV));
        CoolPropDbl C = alphaV*(vL - vV)*(pM - pL) + 0.5*POW2(pL - pV);
        
        // Argument to square root
        CoolPropDbl sqrt_arg = std::abs(B*B/(4*A*A)-C/A);
        
        // If the argument to sqrt is very small, we multiply it by a large factor to make it 
        // larger, and then also divide the sqrt by the sqrt of the factor
        if (std::abs(sqrt_arg) > 1e-10){
            DeltavL = -0.5*B/A + sign((alphaL - alphaV)/alphaV)*sqrt(sqrt_arg);
        }
        else{
            // Scale the argument to sqrt() function to make it about 1.0, and divide by sqrt(factor) to yield a factor of 1
            CoolPropDbl powerfactor = -log10(sqrt_arg);
            DeltavL = -0.5*B/A + sign((alphaL - alphaV)/alphaV)*sqrt(sqrt_arg*powerfactor)/sqrt(powerfactor);
        }
        DeltavV = (pL - pV + alphaL*DeltavL)/alphaV;
        
        // Update the densities of liquid and vapor
        rhoL = 1/(vL + DeltavL);
        rhoV = 1/(vV + DeltavV);
        if (B*B/(4*A*A)-C/A < -10*DBL_EPSILON){
            rhoL *= 1.01;
            rhoV /= 1.01;
        }
        
        // Update the states again 
        SatL->update_DmolarT_direct(rhoL, T);
        SatV->update_DmolarT_direct(rhoV, T);
        
        // Calculate the error (here the relative error in pressure)
        error = std::abs((SatL->p() - SatV->p())/SatL->p());

        if (get_debug_level() > 0){ std::cout << format("[Maxwell] rhoL: %0.16Lg rhoV: %0.16Lg error: %Lg dvL/vL: %Lg dvV/vV: %Lg pL: %Lg pV: %Lg\n", rhoL, rhoV, error, DeltavL/vL, DeltavV/vV, pL, pV);}

        // If the step size is small, start a counter to allow the other density 
        // to be corrected a few times
        if (std::abs(DeltavL*rhoL) < 10*DBL_EPSILON || std::abs(DeltavV*rhoV) < 10*DBL_EPSILON){
            small_step_count++;
        }
        
        iter++;
        if (iter > 30){
            throw SolutionError(format("Maxwell solver did not converge after 30 iterations;  rhoL: %0.16Lg rhoV: %0.16Lg error: %Lg dvL/vL: %Lg dvV/vV: %Lg pL: %Lg pV: %Lg\n", rhoL, rhoV, error, DeltavL/vL, DeltavV/vV, pL, pV));
        }
    }
    while ((SatL->p() < 0) || (error > 1e-5 && small_step_count < 4));
    if (get_debug_level() > 0){ std::cout << format("[Maxwell] pL: %g pV: %g\n", SatL->p(), SatV->p());}
}

void SaturationSolvers::x_and_y_from_K(CoolPropDbl beta, const std::vector<CoolPropDbl> &K, const std::vector<CoolPropDbl> &z, std::vector<CoolPropDbl> &x, std::vector<CoolPropDbl> &y)
{
    for (unsigned int i=0; i < K.size(); i++)
    {
        double denominator = (1-beta+beta*K[i]); // Common denominator
        x[i] = z[i]/denominator;
        y[i] = K[i]*z[i]/denominator;
    }
}

void SaturationSolvers::successive_substitution(HelmholtzEOSMixtureBackend &HEOS, const CoolPropDbl beta, CoolPropDbl T, CoolPropDbl p, const std::vector<CoolPropDbl> &z,
                                                       std::vector<CoolPropDbl> &K, mixture_VLE_IO &options)
{
    int iter = 1;
    CoolPropDbl change, f, df, deriv_liq, deriv_vap;
    std::size_t N = z.size();
    std::vector<CoolPropDbl> ln_phi_liq, ln_phi_vap;
    ln_phi_liq.resize(N); ln_phi_vap.resize(N);

    std::vector<CoolPropDbl> &x = HEOS.SatL->get_mole_fractions_ref(), &y = HEOS.SatV->get_mole_fractions_ref();
    x_and_y_from_K(beta, K, z, x, y);

    HEOS.SatL->specify_phase(iphase_liquid);
    HEOS.SatV->specify_phase(iphase_gas);

    normalize_vector(x);
    normalize_vector(y);

    HEOS.SatL->set_mole_fractions(x);
    HEOS.SatV->set_mole_fractions(y);
    HEOS.SatL->calc_reducing_state();
    HEOS.SatV->calc_reducing_state();
    CoolPropDbl rhomolar_liq = HEOS.SatL->solver_rho_Tp_SRK(T, p, iphase_liquid); // [mol/m^3]
    CoolPropDbl rhomolar_vap = HEOS.SatV->solver_rho_Tp_SRK(T, p, iphase_gas); // [mol/m^3]

    HEOS.SatL->update_TP_guessrho(T, p, rhomolar_liq*1.5);
    HEOS.SatV->update_TP_guessrho(T, p, rhomolar_vap);

    do
    {
        HEOS.SatL->update_TP_guessrho(T, p, HEOS.SatL->rhomolar());
        HEOS.SatV->update_TP_guessrho(T, p, HEOS.SatV->rhomolar());

        f = 0;
        df = 0;

        x_N_dependency_flag xN_flag = XN_INDEPENDENT;
        for (std::size_t i = 0; i < N; ++i)
        {
            ln_phi_liq[i] = MixtureDerivatives::ln_fugacity_coefficient(*(HEOS.SatL.get()), i, xN_flag);
            ln_phi_vap[i] = MixtureDerivatives::ln_fugacity_coefficient(*(HEOS.SatV.get()), i, xN_flag);

            if (options.sstype == imposed_p){
                deriv_liq = MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(*(HEOS.SatL.get()), i, xN_flag);
                deriv_vap = MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(*(HEOS.SatV.get()), i, xN_flag);
            }
            else if (options.sstype == imposed_T){
                deriv_liq = MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(*(HEOS.SatL.get()), i, xN_flag);
                deriv_vap = MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(*(HEOS.SatV.get()), i, xN_flag);
            }
            else {throw ValueError();}

            K[i] = exp(ln_phi_liq[i]-ln_phi_vap[i]);

            f += z[i]*(K[i]-1)/(1-beta+beta*K[i]);

            double dfdK = K[i]*z[i]/pow(1-beta+beta*K[i],(int)2);
            df += dfdK*(deriv_liq-deriv_vap);
        }

        change = -f/df;

        if (options.sstype == imposed_p){
            T += change;
        }
        else if (options.sstype == imposed_T){
            p += change;
        }

        x_and_y_from_K(beta, K, z, x, y);
        normalize_vector(x);
        normalize_vector(y);
        HEOS.SatL->set_mole_fractions(x);
        HEOS.SatV->set_mole_fractions(y);

        iter += 1;
        if (iter > 50)
        {
            throw ValueError(format("saturation_p was unable to reach a solution within 50 iterations"));
        }
    }
    while(std::abs(f) > 1e-12 && iter < options.Nstep_max);

    HEOS.SatL->update_TP_guessrho(T, p, HEOS.SatL->rhomolar());
    HEOS.SatV->update_TP_guessrho(T, p, HEOS.SatV->rhomolar());

    options.p = HEOS.SatL->p();
    options.T = HEOS.SatL->T();
    options.rhomolar_liq = HEOS.SatL->rhomolar();
    options.rhomolar_vap = HEOS.SatV->rhomolar();
    options.x = x;
    options.y = y;
}
void SaturationSolvers::newton_raphson_saturation::resize(std::size_t N)
{
    this->N = N;
    x.resize(N);
    y.resize(N);

    if (imposed_variable == newton_raphson_saturation_options::RHOV_IMPOSED){
        r.resize(N+1);
        negative_r.resize(N+1);
        J.resize(N+1);
        for (std::size_t i = 0; i < N+1; ++i){
            J[i].resize(N+1);
        }
        err_rel.resize(N+1);
    }
    else if (imposed_variable == newton_raphson_saturation_options::P_IMPOSED || imposed_variable == newton_raphson_saturation_options::T_IMPOSED){
        r.resize(N);
        negative_r.resize(N);
        J.resize(N, std::vector<CoolPropDbl>(N, 0));
        err_rel.resize(N);
    }
    else{
        throw ValueError();
    }
}
void SaturationSolvers::newton_raphson_saturation::check_Jacobian()
{
    // References to the classes for concision
    HelmholtzEOSMixtureBackend &rSatL = *(HEOS->SatL.get()), &rSatV = *(HEOS->SatV.get());
    
    // Build the Jacobian and residual vectors
    build_arrays();

    // Make copies of the base
    CoolPropDbl T0 = T;
    std::vector<CoolPropDbl> r0 = r, x0 = x;
    STLMatrix J0 = J;
    CoolPropDbl rhomolar_liq0 = rSatL.rhomolar();
    CoolPropDbl rhomolar_vap0 = rSatV.rhomolar();
    
    {
        // Derivatives with respect to T
        double dT = 1e-3, T1 = T+dT, T2 = T-dT;
        this->T = T1;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<CoolPropDbl> r1 = r;
        this->T = T2;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<CoolPropDbl> r2 = r;
        
        std::vector<CoolPropDbl> diffn(N+1, _HUGE);
        for (std::size_t i = 0; i < N+1; ++i){
            diffn[i] = (r1[i]-r2[i])/(2*dT);
        }
        std::cout << format("For T\n");
        std::cout << "numerical: " << vec_to_string(diffn, "%0.11Lg") << std::endl;
        std::cout << "analytic: " << vec_to_string(get_col(J0, N-1), "%0.11Lg") << std::endl;
    }
    {
        // Derivatives with respect to rho'
        double drho = 1;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0+drho;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<CoolPropDbl> rr1 = r;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0-drho;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<CoolPropDbl> rr2 = r;
        
        std::vector<CoolPropDbl> difffn(N+1, _HUGE);
        for (std::size_t i = 0; i < N+1; ++i){
            difffn[i] = (rr1[i]-rr2[i])/(2*drho);
        }
        std::cout << format("For rho\n");
        std::cout << "numerical: " << vec_to_string(difffn, "%0.11Lg") << std::endl;
        std::cout << "analytic: " << vec_to_string(get_col(J0, N), "%0.11Lg") << std::endl;
    }
    for (std::size_t i = 0; i < x.size()-1;  ++i)
    {
        // Derivatives with respect to x[i]
        double dx = 1e-5;
        this->x = x0; 
        this->x[i] += dx; this->x[x.size()-1] -= dx;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays(); 
        std::vector<CoolPropDbl> r1 = this->r;
    
        this->x = x0; 
        this->x[i] -= dx; this->x[x.size()-1] += dx;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays(); 
        std::vector<CoolPropDbl> r2 = this->r;
    
        std::vector<CoolPropDbl> diffn(N+1, _HUGE);
        for (std::size_t j = 0; j < N+1; ++j){
            diffn[j] = (r1[j]-r2[j])/(2*dx);
        }
        std::cout << format("For x%d N %d\n", i, N);
        std::cout << "numerical: " << vec_to_string(diffn, "%0.11Lg") << std::endl;
        std::cout << "analytic: " << vec_to_string(get_col(J0, i), "%0.11Lg") << std::endl;
    }
}
void SaturationSolvers::newton_raphson_saturation::call(HelmholtzEOSMixtureBackend &HEOS, const std::vector<CoolPropDbl> &z, std::vector<CoolPropDbl> &z_incipient, newton_raphson_saturation_options &IO)
{
    int iter = 0;
    
    if (get_debug_level() > 9){std::cout << " NRsat::call:  p" << IO.p << " T" << IO.T << " dl" << IO.rhomolar_liq << " dv" << IO.rhomolar_vap << std::endl;}
    
    // Reset all the variables and resize
    pre_call();
    
    this->bubble_point = IO.bubble_point;
    rhomolar_liq = IO.rhomolar_liq;
    rhomolar_vap = IO.rhomolar_vap;
    T = IO.T;
    p = IO.p;
    imposed_variable = IO.imposed_variable;
    
    resize(z.size());
    
    if (bubble_point){
        // Bubblepoint, vapor (y) is the incipient phase
        x = z;
        y = z_incipient;
    }
    else{
        // Dewpoint, liquid (x) is the incipient phase
        y = z;
        x = z_incipient;
    }
    
    // Hold a pointer to the backend
    this->HEOS = &HEOS;

    //check_Jacobian();

    do
    {
        // Build the Jacobian and residual vectors
        build_arrays();

        // Solve for the step; v is the step with the contents
        // [delta(x_0), delta(x_1), ..., delta(x_{N-2}), delta(spec)]
        std::vector<CoolPropDbl> v = linsolve(J, negative_r);
        
        if (bubble_point){
            for (unsigned int i = 0; i < N-1; ++i){
                err_rel[i] = v[i]/y[i];
                y[i] += v[i];
            }        
            y[N-1] = 1 - std::accumulate(y.begin(), y.end()-1, 0.0);
        }
        else{
            for (unsigned int i = 0; i < N-1; ++i){
                err_rel[i] = v[i]/x[i];
                x[i] += v[i];
            }        
            x[N-1] = 1 - std::accumulate(x.begin(), x.end()-1, 0.0);
        }
        if (imposed_variable == newton_raphson_saturation_options::P_IMPOSED){
            T += v[N-1]; err_rel[N-1] = v[N-1]/T;
        }
        else if (imposed_variable == newton_raphson_saturation_options::T_IMPOSED){
            p += v[N-1]; err_rel[N-1] = v[N-1]/p;
        }
        else if (imposed_variable == newton_raphson_saturation_options::RHOV_IMPOSED){
            T += v[N-1]; err_rel[N-1] = v[N-1]/T;
            rhomolar_liq += v[N]; err_rel[N] = v[N]/rhomolar_liq;
        }
        else{
            throw ValueError("invalid imposed_variable");
        }
        //std::cout << format("\t%Lg ", this->error_rms) << T << " " << rhomolar_liq << " " << rhomolar_vap << " v " << vec_to_string(v, "%0.10Lg")  << " x " << vec_to_string(x, "%0.10Lg") << " r " << vec_to_string(r, "%0.10Lg") << std::endl;
        
        min_rel_change = min_abs_value(err_rel);
        iter++;
        
        if (iter == IO.Nstep_max){
            throw ValueError(format("newton_raphson_saturation::call reached max number of iterations [%d]",IO.Nstep_max));
        }
    }
    while(this->error_rms > 1e-9 && min_rel_change > 1000*DBL_EPSILON && iter < IO.Nstep_max);

    IO.Nsteps = iter;
    IO.p = p;
    IO.x = x; // Mole fractions in liquid
    IO.y = y; // Mole fractions in vapor
    IO.T = T;
    IO.rhomolar_liq = rhomolar_liq;
    IO.rhomolar_vap = rhomolar_vap;
    IO.hmolar_liq = HEOS.SatL.get()->hmolar();
    IO.hmolar_vap = HEOS.SatV.get()->hmolar();
    IO.smolar_liq = HEOS.SatL.get()->smolar();
    IO.smolar_vap = HEOS.SatV.get()->smolar();
}

void SaturationSolvers::newton_raphson_saturation::build_arrays()
{
    // References to the classes for concision
    HelmholtzEOSMixtureBackend &rSatL = *(HEOS->SatL.get()), &rSatV = *(HEOS->SatV.get());
    
    // Step 0:
    // -------
    // Set mole fractions for the incipient phase
    if (bubble_point){
        // Vapor is incipient phase, set its composition
        rSatV.set_mole_fractions(y);
        rSatL.set_mole_fractions(x);
    }
    else{
        // Liquid is incipient phase, set its composition
        rSatL.set_mole_fractions(x);
        rSatV.set_mole_fractions(y);
    }

    if (imposed_variable == newton_raphson_saturation_options::RHOV_IMPOSED){
        rSatL.update(DmolarT_INPUTS, rhomolar_liq, T);
        rSatV.update(DmolarT_INPUTS, rhomolar_vap, T);
    }
    else if (imposed_variable == newton_raphson_saturation_options::P_IMPOSED || imposed_variable == newton_raphson_saturation_options::T_IMPOSED){
        rSatL.update_TP_guessrho(T, p, rhomolar_liq); rhomolar_liq = rSatL.rhomolar();
        rSatV.update_TP_guessrho(T, p, rhomolar_vap); rhomolar_vap = rSatV.rhomolar();
    }
    else{
        throw ValueError("imposed variable not set for NR VLE");
    }
    
    // For diagnostic purposes calculate the pressures (no derivatives are evaluated)
    CoolPropDbl p_liq = rSatL.p();
    CoolPropDbl p_vap = rSatV.p();
    p = 0.5*(p_liq + p_vap);
    
    // Step 2:
    // -------
    // Build the residual vector and the Jacobian matrix

    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    
    if (imposed_variable == newton_raphson_saturation_options::RHOV_IMPOSED){
        // For the residuals F_i (equality of fugacities)
        for (std::size_t i = 0; i < N; ++i)
        {
            // Equate the liquid and vapor fugacities
            CoolPropDbl ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            CoolPropDbl ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
            r[i] = ln_f_liq - ln_f_vap;
            
            for (std::size_t j = 0; j < N-1; ++j){ // j from 0 to N-2
                if (bubble_point){
                    J[i][j] = -MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(rSatV, i, j, xN_flag);
                }
                else{ 
                    J[i][j] = MixtureDerivatives::dln_fugacity_dxj__constT_rho_xi(rSatL, i, j, xN_flag);
                }
            }
            J[i][N-1] = MixtureDerivatives::dln_fugacity_i_dT__constrho_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_i_dT__constrho_n(rSatV, i, xN_flag);
            J[i][N] = MixtureDerivatives::dln_fugacity_i_drho__constT_n(rSatL, i, xN_flag);
        }
        // ---------------------------------------------------------------
        // Derivatives of pL(T,rho',x)-p(T,rho'',y) with respect to inputs
        // ---------------------------------------------------------------
        r[N] = p_liq - p_vap;
        for (std::size_t j = 0; j < N-1; ++j){ // j from 0 to N-2
            J[N][j] = MixtureDerivatives::dpdxj__constT_V_xi(rSatL, j, xN_flag); // p'' not a function of x0
        }
        // Fixed composition derivatives
        J[N][N-1] = rSatL.first_partial_deriv(iP, iT, iDmolar)-rSatV.first_partial_deriv(iP, iT, iDmolar);
        J[N][N] = rSatL.first_partial_deriv(iP, iDmolar, iT);
    }
    else if (imposed_variable == newton_raphson_saturation_options::P_IMPOSED){
        // Independent variables are N-1 mole fractions of incipient phase and T
        
        // For the residuals F_i (equality of fugacities)
        for (std::size_t i = 0; i < N; ++i)
        {
            // Equate the liquid and vapor fugacities
            CoolPropDbl ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            CoolPropDbl ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
            r[i] = ln_f_liq - ln_f_vap;
            
            for (std::size_t j = 0; j < N-1; ++j){ // j from 0 to N-2
                if (bubble_point){
                    J[i][j] = -MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatV, i, j, xN_flag);
                }
                else{ 
                    J[i][j] = MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatL, i, j, xN_flag);
                }
            }
            J[i][N-1] = MixtureDerivatives::dln_fugacity_i_dT__constp_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_i_dT__constp_n(rSatV, i, xN_flag);
        }
    }
    else if (imposed_variable == newton_raphson_saturation_options::T_IMPOSED){
        // Independent variables are N-1 mole fractions of incipient phase and p
        
        // For the residuals F_i (equality of fugacities)
        for (std::size_t i = 0; i < N; ++i)
        {
            // Equate the liquid and vapor fugacities
            CoolPropDbl ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            CoolPropDbl ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
            r[i] = ln_f_liq - ln_f_vap;
            
            for (std::size_t j = 0; j < N-1; ++j){ // j from 0 to N-2
                if (bubble_point){
                    J[i][j] = -MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatV, i, j, xN_flag);
                }
                else{
                    J[i][j] = MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatL, i, j, xN_flag);
                }
            }
            J[i][N-1] = MixtureDerivatives::dln_fugacity_i_dp__constT_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_i_dp__constT_n(rSatV, i, xN_flag);
        }
    }
    else{
        throw ValueError();
    }

    // Flip all the signs of the entries in the residual vector since we are solving Jv = -r, not Jv=r
    // Also calculate the rms error of the residual vector at this step
    error_rms = 0;
    for (unsigned int i = 0; i < J.size(); ++i)
    {
        negative_r[i] = -r[i];
        error_rms += r[i]*r[i]; // Sum the squares
    }
    error_rms = sqrt(error_rms); // Square-root (The R in RMS)
    
    // Calculate derivatives along phase boundary;
    // Gernert thesis 3.96 and 3.97
    CoolPropDbl dQ_dPsat = 0, dQ_dTsat = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        dQ_dPsat += x[i]*(MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rSatV, i, xN_flag));
        dQ_dTsat += x[i]*(MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(rSatV, i, xN_flag));
    }
    dTsat_dPsat = -dQ_dPsat/dQ_dTsat;
    dPsat_dTsat = -dQ_dTsat/dQ_dPsat;
}

void SaturationSolvers::newton_raphson_twophase::call(HelmholtzEOSMixtureBackend &HEOS, newton_raphson_twophase_options &IO)
{
    int iter = 0;
    
    if (get_debug_level() > 9){std::cout << " NRsat::call:  p" << IO.p << " T" << IO.T << " dl" << IO.rhomolar_liq << " dv" << IO.rhomolar_vap << std::endl;}
    
    // Reset all the variables and resize
    pre_call();
    
    rhomolar_liq = IO.rhomolar_liq;
    rhomolar_vap = IO.rhomolar_vap;
    T = IO.T;
    p = IO.p;
    imposed_variable = IO.imposed_variable;
    x = IO.x;
    y = IO.y;
    z = IO.z;
    beta = IO.beta;
    
    this->N = z.size();
    x.resize(N);
    y.resize(N);
    r.resize(2*N-1);
    negative_r.resize(2*N-1);
    J.resize(2*N-1, std::vector<CoolPropDbl>(2*N-1, 0));
    err_rel.resize(2*N-1);
    
    // Hold a pointer to the backend
    this->HEOS = &HEOS;

    do
    {
        // Build the Jacobian and residual vectors
        build_arrays();

        // Solve for the step; v is the step with the contents
        // [delta(x_0), delta(x_1), ..., delta(x_{N-2}), delta(spec)]
        
        // Uncomment to see Jacobian and residual at every step
        // std::cout << vec_to_string(J, "%0.12Lg") << std::endl;
        // std::cout << vec_to_string(negative_r, "%0.12Lg") << std::endl;
        
        std::vector<CoolPropDbl> v = linsolve(J, negative_r);
        for (unsigned int i = 0; i < N-1; ++i){
            err_rel[i] = v[i]/x[i];
            x[i] += v[i];
            err_rel[i+(N-1)] = v[i+(N-1)]/y[i];
            y[i] += v[i+(N-1)];
        }        
        x[N-1] = 1 - std::accumulate(x.begin(), x.end()-1, 0.0);
        y[N-1] = 1 - std::accumulate(y.begin(), y.end()-1, 0.0);
            
        if (imposed_variable == newton_raphson_twophase_options::P_IMPOSED){
            T += v[2*N-2]; err_rel[2*N-2] = v[2*N-2]/T;
        }
        else if (imposed_variable == newton_raphson_twophase_options::T_IMPOSED){
            p += v[2*N-2]; err_rel[2*N-2] = v[2*N-2]/p;
        }
        else{
            throw ValueError("invalid imposed_variable");
        }
        //std::cout << format("\t%Lg ", this->error_rms) << T << " " << rhomolar_liq << " " << rhomolar_vap << " v " << vec_to_string(v, "%0.10Lg")  << " x " << vec_to_string(x, "%0.10Lg") << " r " << vec_to_string(r, "%0.10Lg") << std::endl;
        
        min_rel_change = min_abs_value(err_rel);
        iter++;
        
        if (iter == IO.Nstep_max){
            throw ValueError(format("newton_raphson_saturation::call reached max number of iterations [%d]",IO.Nstep_max));
        }
    }
    while(this->error_rms > 1e-9 && min_rel_change > 1000*DBL_EPSILON && iter < IO.Nstep_max);

    IO.Nsteps = iter;
    IO.p = p;
    IO.x = x; // Mole fractions in liquid
    IO.y = y; // Mole fractions in vapor
    IO.T = T;
    IO.rhomolar_liq = rhomolar_liq;
    IO.rhomolar_vap = rhomolar_vap;
    IO.hmolar_liq = HEOS.SatL.get()->hmolar();
    IO.hmolar_vap = HEOS.SatV.get()->hmolar();
    IO.smolar_liq = HEOS.SatL.get()->smolar();
    IO.smolar_vap = HEOS.SatV.get()->smolar();
}

void SaturationSolvers::newton_raphson_twophase::build_arrays()
{
    // References to the classes for concision
    HelmholtzEOSMixtureBackend &rSatL = *(HEOS->SatL.get()), &rSatV = *(HEOS->SatV.get());
    
    // Step 0:
    // -------
    // Set mole fractions 
    rSatL.set_mole_fractions(x);
    rSatV.set_mole_fractions(y);

    //std::vector<CoolPropDbl> &x = rSatL.get_mole_fractions();
    //std::vector<CoolPropDbl> &y = rSatV.get_mole_fractions();

    rSatL.update_TP_guessrho(T, p, rhomolar_liq); rhomolar_liq = rSatL.rhomolar();
    rSatV.update_TP_guessrho(T, p, rhomolar_vap); rhomolar_vap = rSatV.rhomolar();
    
    // For diagnostic purposes calculate the pressures (no derivatives are evaluated)
    CoolPropDbl p_liq = rSatL.p();
    CoolPropDbl p_vap = rSatV.p();
    p = 0.5*(p_liq + p_vap);
    
    // Step 2:
    // -------
    // Build the residual vector and the Jacobian matrix
    
    x_N_dependency_flag xN_flag = XN_DEPENDENT;
    
    // Form of residuals do not depend on which variable is imposed
    for (std::size_t i = 0; i < N; ++i)
    {
        // Equate the liquid and vapor fugacities
        CoolPropDbl ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
        CoolPropDbl ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
        r[i] = ln_f_liq - ln_f_vap; // N of these
    
        if (i != N-1){
            // Equate the specified vapor mole fraction and that given defined by the ith component
            r[i+N] = (z[i]-x[i])/(y[i]-x[i]) - beta; // N-1 of these
        }
    }

    // First part of derivatives with respect to ln f_i
    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = 0; j < N-1; ++j){
            J[i][j] = MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatL, i, j, xN_flag);
            J[i][j+N-1] = -MixtureDerivatives::dln_fugacity_dxj__constT_p_xi(rSatV, i, j, xN_flag);
        }
                
        // Last derivative with respect to either T or p depending on what is imposed
        if (imposed_variable == newton_raphson_twophase_options::P_IMPOSED){
            J[i][2*N-2] = MixtureDerivatives::dln_fugacity_i_dT__constp_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_i_dT__constp_n(rSatV, i, xN_flag);
        }
        else if (imposed_variable == newton_raphson_twophase_options::T_IMPOSED){
            J[i][2*N-2] = MixtureDerivatives::dln_fugacity_i_dp__constT_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_i_dp__constT_n(rSatV, i, xN_flag);
        }
        else{
            throw ValueError();
        }
    }
    // Derivatives with respect to the vapor mole fractions residual
    for (std::size_t i = 0; i < N-1; ++i)
    {
        std::size_t k = i + N; // N ln f_i residuals
        J[k][i] =  (z[i]-y[i])/pow(y[i]-x[i], 2);
        J[k][i+(N-1)] = -(z[i]-x[i])/pow(y[i]-x[i], 2);
    }

    // Flip all the signs of the entries in the residual vector since we are solving Jv = -r, not Jv=r
    // Also calculate the rms error of the residual vector at this step
    error_rms = 0;
    for (unsigned int i = 0; i < J.size(); ++i)
    {
        negative_r[i] = -r[i];
        error_rms += r[i]*r[i]; // Sum the squares
    }
    error_rms = sqrt(error_rms); // Square-root (The R in RMS)
}

} /* namespace CoolProp*/
