
#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "MatrixMath.h"
#include "MixtureDerivatives.h"

namespace CoolProp {
    
void SaturationSolvers::saturation_critical(HelmholtzEOSMixtureBackend &HEOS, parameters ykey, long double y){
    
    class inner_resid : public FuncWrapper1D{
        public:
        HelmholtzEOSMixtureBackend *HEOS;
        long double T, desired_p, rhomolar_liq, calc_p, rhomolar_crit;
        
        inner_resid(HelmholtzEOSMixtureBackend *HEOS, long double T, long double desired_p)
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
        long double y;
        long double r, T, rhomolar_liq, rhomolar_vap, value, p, gL, gV, rhomolar_crit;
        int other;

        outer_resid(HelmholtzEOSMixtureBackend &HEOS, CoolProp::parameters ykey, long double y) 
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

void SaturationSolvers::saturation_T_pure_1D_P(HelmholtzEOSMixtureBackend &HEOS, long double T, saturation_T_pure_options &options)
{
    
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        long double r, T, rhomolar_liq, rhomolar_vap, value, p, gL, gV;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend &HEOS, long double T, long double rhomolar_liq_guess, long double rhomolar_vap_guess) 
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
        long double pmax = std::min(options.p*1.03, static_cast<long double>(HEOS.p_critical()+1e-6));
        long double pmin = std::max(options.p*0.97, static_cast<long double>(HEOS.p_triple()-1e-6));
        Brent(resid, pmin, pmax, LDBL_EPSILON, 1e-8, 100, errstr);
    }
}

void SaturationSolvers::saturation_P_pure_1D_T(HelmholtzEOSMixtureBackend &HEOS, long double p, saturation_PHSU_pure_options &options){
    
    // Define the residual to be driven to zero
    class solver_resid : public FuncWrapper1D
    {
    public:

        HelmholtzEOSMixtureBackend *HEOS;
        long double r, p, rhomolar_liq, rhomolar_vap, value, T, gL, gV;
        int other;

        solver_resid(HelmholtzEOSMixtureBackend &HEOS, long double p, long double rhomolar_liq_guess, long double rhomolar_vap_guess) 
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
    long double Tmax = std::min(options.T + 2, static_cast<long double>(HEOS.T_critical()-1e-6));
    long double Tmin = std::max(options.T - 2, static_cast<long double>(HEOS.Ttriple()+1e-6));
    Brent(resid, Tmin, Tmax, LDBL_EPSILON, 1e-11, 100, errstr);
}
    
void SaturationSolvers::saturation_PHSU_pure(HelmholtzEOSMixtureBackend &HEOS, long double specified_value, saturation_PHSU_pure_options &options)
{
    /*
    This function is inspired by the method of Akasaka:

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from
    Helmholtz Energy Equations of State",
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */
    std::vector<long double> negativer(3,_HUGE), v;
    std::vector<std::vector<long double> > J(3, std::vector<long double>(3,_HUGE));

    HEOS.calc_reducing_state();
    const SimpleState & reduce = HEOS.get_reducing_state();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    long double T, rhoL, rhoV, pL, pV;
    long double deltaL=0, deltaV=0, tau=0, error;
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
			catch(std::exception &e)
			{
				throw ValueError("Unable to invert ancillary equation");
			}
        }
        else
        {
            throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }
		if (T > HEOS.T_critical()-1){ T -= 1; }

        // Evaluate densities from the ancillary equations
        rhoV = HEOS.get_components()[0]->ancillaries.rhoV.evaluate(T);
        rhoL = HEOS.get_components()[0]->ancillaries.rhoL.evaluate(T);

        // Apply a single step of Newton's method to improve guess value for liquid
        // based on the error between the gas pressure (which is usually very close already)
        // and the liquid pressure, which can sometimes (especially at low pressure),
        // be way off, and often times negative
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        rhoL += -(SatL->p()-SatV->p())/SatL->first_partial_deriv(iP, iDmolar, iT);
        
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
        pV = SatV->p();

        // These derivatives are needed for both cases
        long double alpharL = SatL->alphar();
        long double alpharV = SatV->alphar();
        long double dalphar_dtauL = SatL->dalphar_dTau();
        long double dalphar_dtauV = SatV->dalphar_dTau();
        long double d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        long double d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        long double dalphar_ddeltaL = SatL->dalphar_dDelta();
        long double dalphar_ddeltaV = SatV->dalphar_dDelta();
        long double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        long double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();

        // -r_1
        negativer[0] = -(deltaV*(1+deltaV*dalphar_ddeltaV)-deltaL*(1+deltaL*dalphar_ddeltaL));
        // -r_2
        negativer[1] = -(deltaV*dalphar_ddeltaV+alpharV+log(deltaV)-deltaL*dalphar_ddeltaL-alpharL-log(deltaL));
        switch (options.specified_variable){
            case saturation_PHSU_pure_options::IMPOSED_PL:
                // -r_3
                negativer[2] = -(pL/specified_value - 1); break;
            case saturation_PHSU_pure_options::IMPOSED_PV:
                // -r_3
                negativer[2] = -(pV/specified_value - 1); break;
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
            default:
                throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }

        v = linsolve(J, negativer);

        tau += options.omega*v[0];

        if (options.use_logdelta){
            deltaL = exp(log(deltaL)+options.omega*v[1]);
            deltaV = exp(log(deltaV)+options.omega*v[2]);
        }
        else{
            deltaL += options.omega*v[1];
            deltaV += options.omega*v[2];
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
        if (iter > 25){
            // Set values back into the options structure for use in next solver
            options.rhoL = rhoL; options.rhoV = rhoV; options.T = T;
            // Error out
			std::string info = get_parameter_information(specified_parameter, "short");
            throw SolutionError(format("saturation_PHSU_pure solver did not converge after 25 iterations for %s=%Lg current error is %Lg", info.c_str(), specified_value, error));
        }
    }
    while (error > 1e-9);
}
void SaturationSolvers::saturation_D_pure(HelmholtzEOSMixtureBackend &HEOS, long double rhomolar, saturation_D_pure_options &options)
{
    /*
    This function is inspired by the method of Akasaka:

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from
    Helmholtz Energy Equations of State",
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */
    std::vector<long double> r(2,_HUGE), v;
    std::vector<std::vector<long double> > J(2, std::vector<long double>(2,_HUGE));

    HEOS.calc_reducing_state();
    const SimpleState & reduce = HEOS.get_reducing_state();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    long double T, rhoL,rhoV;
    long double deltaL=0, deltaV=0, tau=0, error, p_error;
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

        long double pL = SatL->p();
        long double pV = SatV->p();

        // These derivatives are needed for both cases
        long double dalphar_dtauL = SatL->dalphar_dTau();
        long double dalphar_dtauV = SatV->dalphar_dTau();
        long double d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        long double d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        long double alpharL = SatL->alphar();
        long double alpharV = SatV->alphar();
        long double dalphar_ddeltaL = SatL->dalphar_dDelta();
        long double dalphar_ddeltaV = SatV->dalphar_dDelta();

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
            long double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();
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
            long double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
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
	long double p_error_limit = 1e-3;
	if (std::abs(p_error) > p_error_limit){
		throw SolutionError(format("saturation_D_pure solver abs error on p [%Lg] > limit [%Lg]", p_error, p_error_limit));
	}
}
void SaturationSolvers::saturation_T_pure(HelmholtzEOSMixtureBackend &HEOS, long double T, saturation_T_pure_options &options)
{
    // Set some imput options
    SaturationSolvers::saturation_T_pure_Akasaka_options _options;
    _options.omega = 1.0;
    _options.use_guesses = false;
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
void SaturationSolvers::saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend &HEOS, long double T, saturation_T_pure_Akasaka_options &options)
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
    long double R_u = HEOS.gas_constant();
    shared_ptr<HelmholtzEOSMixtureBackend> SatL = HEOS.SatL,
                                           SatV = HEOS.SatV;

    long double rhoL = _HUGE, rhoV = _HUGE,JL,JV,KL,KV,dJL,dJV,dKL,dKV;
    long double DELTA, deltaL=0, deltaV=0, error, PL, PV, stepL, stepV;
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
				rhoL += -(SatL->p()-SatV->p())/SatL->first_partial_deriv(iP, iDmolar, iT);
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
    //		std::cout << format("%s:%d: Akasaka guess values deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
    //	}

    do{
        /*if (get_debug_level()>8){
            std::cout << format("%s:%d: right before the derivs with deltaL = %g deltaV = %g tau = %g\n",__FILE__,__LINE__,deltaL, deltaV, tau).c_str();
        }*/

        // Calculate once to save on calls to EOS
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        long double alpharL = SatL->alphar();
        long double alpharV = SatV->alphar();
        long double dalphar_ddeltaL = SatL->dalphar_dDelta();
        long double dalphar_ddeltaV = SatV->dalphar_dDelta();
        long double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        long double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();

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

        if (deltaL+stepL > 1 && deltaV+stepV < 1 && deltaV+stepV > 0){
            deltaL += stepL;  deltaV += stepV;
            rhoL = deltaL*reduce.rhomolar; rhoV = deltaV*reduce.rhomolar;
        }
        else{
            throw ValueError(format("rhosatPure_Akasaka densities are crossed"));
        }
        iter++;
        if (iter > 100){
            throw SolutionError(format("Akasaka solver did not converge after 100 iterations"));
        }
    }
    while (error > 1e-10 && std::abs(stepL) > 10*DBL_EPSILON*std::abs(stepL) && std::abs(stepV) > 10*DBL_EPSILON*std::abs(stepV));
	
	long double p_error_limit = 1e-3;
	long double p_error = (PL - PV)/PL;
	if (std::abs(p_error) > p_error_limit){
        options.pL = PL;
        options.pV = PV;
        options.rhoL = rhoL;
        options.rhoV = rhoV;
		throw SolutionError(format("saturation_T_pure_Akasaka solver abs error on p [%g] > limit [%g]", std::abs(p_error), p_error_limit));
	}
}

void SaturationSolvers::x_and_y_from_K(long double beta, const std::vector<long double> &K, const std::vector<long double> &z, std::vector<long double> &x, std::vector<long double> &y)
{
    for (unsigned int i=0; i < K.size(); i++)
    {
        double denominator = (1-beta+beta*K[i]); // Common denominator
        x[i] = z[i]/denominator;
        y[i] = K[i]*z[i]/denominator;
    }
}

void SaturationSolvers::successive_substitution(HelmholtzEOSMixtureBackend &HEOS, const long double beta, long double T, long double p, const std::vector<long double> &z,
                                                       std::vector<long double> &K, mixture_VLE_IO &options)
{
    int iter = 1;
    long double change, f, df, deriv_liq, deriv_vap;
    std::size_t N = z.size();
    std::vector<long double> ln_phi_liq, ln_phi_vap;
    ln_phi_liq.resize(N); ln_phi_vap.resize(N);

    std::vector<long double> &x = HEOS.SatL->get_mole_fractions(), &y = HEOS.SatV->get_mole_fractions();
    x_and_y_from_K(beta, K, z, x, y);

    HEOS.SatL->specify_phase(iphase_liquid);
    HEOS.SatV->specify_phase(iphase_gas);

    normalize_vector(x);
    normalize_vector(y);

    HEOS.SatL->set_mole_fractions(x);
    HEOS.SatV->set_mole_fractions(y);
    HEOS.SatL->calc_reducing_state();
    HEOS.SatV->calc_reducing_state();
    long double rhomolar_liq = HEOS.SatL->solver_rho_Tp_SRK(T, p, iphase_liquid); // [mol/m^3]
    long double rhomolar_vap = HEOS.SatV->solver_rho_Tp_SRK(T, p, iphase_gas); // [mol/m^3]

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

    HEOS.SatL->update_TP_guessrho(T, p, rhomolar_liq);
    HEOS.SatV->update_TP_guessrho(T, p, rhomolar_vap);

    options.p = HEOS.SatL->p();
    options.T = HEOS.SatL->T();
    options.rhomolar_liq = HEOS.SatL->rhomolar();
    options.rhomolar_vap = HEOS.SatV->rhomolar();
    options.x = x;
    options.y = y;
}
void SaturationSolvers::newton_raphson_saturation::resize(unsigned int N)
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
        J.resize(N, std::vector<long double>(N, 0));
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
    long double T0 = T;
    std::vector<long double> r0 = r, x0 = x;
    STLMatrix J0 = J;
    long double rhomolar_liq0 = rSatL.rhomolar();
    long double rhomolar_vap0 = rSatV.rhomolar();
    
    {
        // Derivatives with respect to T
        double dT = 1e-3, T1 = T+dT, T2 = T-dT;
        this->T = T1;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<long double> r1 = r;
        this->T = T2;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<long double> r2 = r;
        
        std::vector<long double> diffn(N+1, _HUGE);
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
        std::vector<long double> rr1 = r;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0-drho;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays();
        std::vector<long double> rr2 = r;
        
        std::vector<long double> difffn(N+1, _HUGE);
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
        std::vector<long double> r1 = this->r;
    
        this->x = x0; 
        this->x[i] -= dx; this->x[x.size()-1] += dx;
        this->T = T0;
        this->rhomolar_liq = rhomolar_liq0;
        this->rhomolar_vap = rhomolar_vap0;
        build_arrays(); 
        std::vector<long double> r2 = this->r;
    
        std::vector<long double> diffn(N+1, _HUGE);
        for (std::size_t j = 0; j < N+1; ++j){
            diffn[j] = (r1[j]-r2[j])/(2*dx);
        }
        std::cout << format("For x%d N %d\n", i, N);
        std::cout << "numerical: " << vec_to_string(diffn, "%0.11Lg") << std::endl;
        std::cout << "analytic: " << vec_to_string(get_col(J0, i), "%0.11Lg") << std::endl;
    }
}
void SaturationSolvers::newton_raphson_saturation::call(HelmholtzEOSMixtureBackend &HEOS, const std::vector<long double> &z, std::vector<long double> &z_incipient, newton_raphson_saturation_options &IO)
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
        std::vector<long double> v = linsolve(J, negative_r);
        
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
    std::vector<long double> &x = rSatL.get_mole_fractions();
    std::vector<long double> &y = rSatV.get_mole_fractions();
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
    long double p_liq = rSatL.p();
    long double p_vap = rSatV.p();
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
            long double ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            long double ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
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
            long double ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            long double ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
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
            long double ln_f_liq = log(MixtureDerivatives::fugacity_i(rSatL, i, xN_flag));
            long double ln_f_vap = log(MixtureDerivatives::fugacity_i(rSatV, i, xN_flag));
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
    long double dQ_dPsat = 0, dQ_dTsat = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        dQ_dPsat += x[i]*(MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_coefficient_dp__constT_n(rSatV, i, xN_flag));
        dQ_dTsat += x[i]*(MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(rSatL, i, xN_flag) - MixtureDerivatives::dln_fugacity_coefficient_dT__constp_n(rSatV, i, xN_flag));
    }
    dTsat_dPsat = -dQ_dPsat/dQ_dTsat;
    dPsat_dTsat = -dQ_dTsat/dQ_dPsat;
}

} /* namespace CoolProp*/
