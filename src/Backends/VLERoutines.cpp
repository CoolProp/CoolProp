
#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "MatrixMath.h"

namespace CoolProp {

void SaturationSolvers::saturation_PHSU_pure(HelmholtzEOSMixtureBackend *HEOS, long double specified_value, saturation_PHSU_pure_options &options)
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
    
    HEOS->calc_reducing_state();
    const SimpleState & reduce = HEOS->get_reducing();
    HelmholtzEOSMixtureBackend *SatL = HEOS->SatL, *SatV = HEOS->SatV;
    const std::vector<long double> & mole_fractions = HEOS->get_mole_fractions();
    const long double R_u = HEOS->gas_constant();
    
    long double T, rhoL,rhoV;
    long double deltaL=0, deltaV=0, tau=0, error;
    int iter=0;

    // Use the density ancillary function as the starting point for the solver
    try
    {
        if (options.specified_variable == saturation_PHSU_pure_options::IMPOSED_PL || options.specified_variable == saturation_PHSU_pure_options::IMPOSED_PV)
        {
            // Invert liquid density ancillary to get temperature 
            // TODO: fit inverse ancillaries too
            T = HEOS->get_components()[0]->ancillaries.pL.invert(specified_value);
        }
        else
        {
            throw ValueError(format("options.specified_variable to saturation_PHSU_pure [%d] is invalid",options.specified_variable));
        }

        // Evaluate densities from the ancillary equations
        rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T);
        rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T);

        // Apply a single step of Newton's method to improve guess value for liquid
        // based on the error between the gas pressure (which is usually very close already)
        // and the liquid pressure, which can sometimes (especially at low pressure),
        // be way off, and often times negative
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        rhoL += -(SatL->p()-SatV->p())/SatL->first_partial_deriv(iP, iDmolar, iT);

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

        // Calculate once to save on calls to EOS
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);

        double pL = SatL->p();
        double pV = SatV->p();
        
        // These derivatives are needed for both cases
        double alpharL = SatL->alphar();
        double alpharV = SatV->alphar();
        double dalphar_dtauL = SatL->dalphar_dTau();
        double dalphar_dtauV = SatV->dalphar_dTau();
        double d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        double d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        double dalphar_ddeltaL = SatL->dalphar_dDelta();
        double dalphar_ddeltaV = SatV->dalphar_dDelta();
        double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();

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
        
        error = sqrt(pow(negativer[0], 2)+pow(negativer[1], 2)+pow(negativer[2], 2));
        iter++;
        if (T < 0)
        {
            throw SolutionError(format("saturation_PHSU_pure solver T < 0"));
        }
        if (iter > 200){
            throw SolutionError(format("saturation_PHSU_pure solver did not converge after 100 iterations with specified value: %g with index %d",specified_value, options.specified_variable));
        }
    }
    while (error > 1e-11);
}
void SaturationSolvers::saturation_D_pure(HelmholtzEOSMixtureBackend *HEOS, long double rhomolar, saturation_D_pure_options &options)
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
    
    HEOS->calc_reducing_state();
    const SimpleState & reduce = HEOS->get_reducing();
    HelmholtzEOSMixtureBackend *SatL = HEOS->SatL, *SatV = HEOS->SatV;
    const std::vector<long double> & mole_fractions = HEOS->get_mole_fractions();
    
    long double T, rhoL,rhoV;
    long double deltaL=0, deltaV=0, tau=0, error;
    int iter=0;

    // Use the density ancillary function as the starting point for the solver
    try
    {
        if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOL)
        {
            // Invert liquid density ancillary to get temperature 
            // TODO: fit inverse ancillaries too
            T = HEOS->get_components()[0]->ancillaries.rhoL.invert(rhomolar);
            rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T);
            rhoL = rhomolar;
        }
        else if (options.imposed_rho == saturation_D_pure_options::IMPOSED_RHOV)
        {
            // Invert vapor density ancillary to get temperature 
            // TODO: fit inverse ancillaries too
            T = HEOS->get_components()[0]->ancillaries.rhoV.invert(rhomolar);
            rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T);
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

        double pL = SatL->p();
        double pV = SatV->p();
        
        // These derivatives are needed for both cases
        double dalphar_dtauL = SatL->dalphar_dTau();
        double dalphar_dtauV = SatV->dalphar_dTau();
        double d2alphar_ddelta_dtauL = SatL->d2alphar_dDelta_dTau();
        double d2alphar_ddelta_dtauV = SatV->d2alphar_dDelta_dTau();
        double alpharL = SatL->alphar();
        double alpharV = SatV->alphar();
        double dalphar_ddeltaL = SatL->dalphar_dDelta();
        double dalphar_ddeltaV = SatV->dalphar_dDelta();
        
        
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
            double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();
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
            double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
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

        double DET = J[0][0]*J[1][1]-J[0][1]*J[1][0];

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
}
void SaturationSolvers::saturation_T_pure(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_options &options)
{
    // Set some imput options
    SaturationSolvers::saturation_T_pure_Akasaka_options _options;
    _options.omega = 1.0;
    _options.use_guesses = false;
    // Actually call the solver
    SaturationSolvers::saturation_T_pure_Akasaka(HEOS, T, _options);
}
void SaturationSolvers::saturation_T_pure_Akasaka(HelmholtzEOSMixtureBackend *HEOS, long double T, saturation_T_pure_Akasaka_options &options)
{
    // Start with the method of Akasaka

    /*
    This function implements the method of Akasaka 

    R. Akasaka,"A Reliable and Useful Method to Determine the Saturation State from 
    Helmholtz Energy Equations of State", 
    Journal of Thermal Science and Technology v3 n3,2008

    Ancillary equations are used to get a sensible starting point
    */
    
    HEOS->calc_reducing_state();
    const SimpleState & reduce = HEOS->get_reducing();
    long double R_u = HEOS->calc_gas_constant();
    HelmholtzEOSMixtureBackend *SatL = HEOS->SatL, *SatV = HEOS->SatV;

    long double rhoL,rhoV,JL,JV,KL,KV,dJL,dJV,dKL,dKV;
    long double DELTA, deltaL=0, deltaV=0, tau=0, error, PL, PV, stepL, stepV;
    int iter=0;
    // Use the density ancillary function as the starting point for the solver
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
            // If very close to the critical temp, evaluate the ancillaries for a slightly lower temperature
            if (T > 0.99*HEOS->get_reducing().T){
                rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T-1);
                rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T-1);
            }
            else{
                rhoL = HEOS->get_components()[0]->ancillaries.rhoL.evaluate(T);
                rhoV = HEOS->get_components()[0]->ancillaries.rhoV.evaluate(T);
            }
        }    

        // Apply a single step of Newton's method to improve guess value for liquid
        // based on the error between the gas pressure (which is usually very close already)
        // and the liquid pressure, which can sometimes (especially at low pressure),
        // be way off, and often times negative
        SatL->update(DmolarT_INPUTS, rhoL, T);
        SatV->update(DmolarT_INPUTS, rhoV, T);
        rhoL += -(SatL->p()-SatV->p())/SatL->first_partial_deriv(iP, iDmolar, iT);

        deltaL = rhoL/reduce.rhomolar;
        deltaV = rhoV/reduce.rhomolar;
        tau = reduce.T/T;
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
        double alpharL = SatL->alphar();
        double alpharV = SatV->alphar();
        double dalphar_ddeltaL = SatL->dalphar_dDelta();
        double dalphar_ddeltaV = SatV->dalphar_dDelta();
        double d2alphar_ddelta2L = SatL->d2alphar_dDelta2();
        double d2alphar_ddelta2V = SatV->d2alphar_dDelta2();
        
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

        error = fabs(KL-KV)+fabs(JL-JV);

        //  Get the predicted step
        stepL = options.omega/DELTA*( (KV-KL)*dJV-(JV-JL)*dKV);
        stepV = options.omega/DELTA*( (KV-KL)*dJL-(JV-JL)*dKL);
        
        if (deltaL+stepL > 1 && deltaV+stepV < 1 && deltaV+stepV > 0){
            deltaL += stepL;  deltaV += stepV; 
            rhoL = deltaL*reduce.rhomolar; rhoV = deltaV*reduce.rhomolar;
        }
        else{
            throw ValueError(format("rhosatPure_Akasaka failed"));
        }
        iter++;
        if (iter > 100){
            throw SolutionError(format("Akasaka solver did not converge after 100 iterations"));
        }
    }
    while (error > 1e-10 && fabs(stepL) > 10*DBL_EPSILON*fabs(stepL) && fabs(stepV) > 10*DBL_EPSILON*fabs(stepV));
}

void SaturationSolvers::x_and_y_from_K(long double beta, const std::vector<long double> &K, const std::vector<long double> &z, std::vector<long double> &x, std::vector<long double> &y)
{
    for (unsigned int i=0; i < K.size(); i++)
    {
        double denominator = (1-beta+beta*K[i]); // Common denominator
        x[i] = z[i]/denominator;
        y[i] = K[i]*z[i]/denominator;
    }
    //normalize_vector(x);
    //normalize_vector(y);
}

long double SaturationSolvers::successive_substitution(HelmholtzEOSMixtureBackend *HEOS, const long double beta, long double T, long double p, const std::vector<long double> &z, 
                                                       std::vector<long double> &K, mixture_VLE_IO &options)
{
    int iter = 1;
    long double change, f, df, deriv_liq, deriv_vap;
    std::size_t N = z.size();
    std::vector<long double> x, y, ln_phi_liq, ln_phi_vap;
    ln_phi_liq.resize(N); ln_phi_vap.resize(N); x.resize(N); y.resize(N);

    x_and_y_from_K(beta, K, z, x, y);
    HelmholtzEOSMixtureBackend *SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
                               *SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid);
    SatV->specify_phase(iphase_gas);

    SatL->set_mole_fractions(x);
    SatV->set_mole_fractions(y);
    long double rhomolar_liq = SatL->solver_rho_Tp_SRK(T, p, iphase_liquid); // [mol/m^3]
    long double rhomolar_vap = SatV->solver_rho_Tp_SRK(T, p, iphase_gas); // [mol/m^3]
    
    do
    {
        
        SatL->update_TP_guessrho(T, p, rhomolar_liq);
        SatV->update_TP_guessrho(T, p, rhomolar_vap);

        f = 0;
        df = 0;

        for (std::size_t i = 0; i < N; ++i)
        {
            ln_phi_liq[i] = SatL->mixderiv_ln_fugacity_coefficient(i);
            ln_phi_vap[i] = SatV->mixderiv_ln_fugacity_coefficient(i);

            if (options.sstype == imposed_p){
                deriv_liq = SatL->mixderiv_dln_fugacity_coefficient_dT__constp_n(i);
                deriv_vap = SatV->mixderiv_dln_fugacity_coefficient_dT__constp_n(i);
            }
            else if (options.sstype == imposed_T){
                deriv_liq = SatL->mixderiv_dln_fugacity_coefficient_dp__constT_n(i);
                deriv_vap = SatV->mixderiv_dln_fugacity_coefficient_dp__constT_n(i);
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
        SatL->set_mole_fractions(x);
        SatV->set_mole_fractions(y);

        iter += 1;
        if (iter > 50)
        {
            return _HUGE;
            //throw ValueError(format("saturation_p was unable to reach a solution within 50 iterations"));
        }
        rhomolar_liq = SatL->rhomolar();
        rhomolar_vap = SatV->rhomolar();
    }
    while(fabs(f) > 1e-12 || iter < options.Nstep_max);

    SatL->update_TP_guessrho(T, p, rhomolar_liq);
    SatV->update_TP_guessrho(T, p, rhomolar_vap);
    
    double pL = SatL->calc_pressure();
    double pV = SatV->calc_pressure();


    options.rhomolar_liq = SatL->rhomolar();
    options.rhomolar_vap = SatV->rhomolar();
    options.p = pL;
    options.T = T;
    options.x = SatL->get_mole_fractions();
    options.y = SatV->get_mole_fractions();

    delete SatL; delete SatV;
}

void SaturationSolvers::newton_raphson_VLE_GV::resize(unsigned int N)
{
    this->N = N;
    x.resize(N); 
    y.resize(N); 
    phi_ij_liq.resize(N); 
    phi_ij_vap.resize(N);
    dlnphi_drho_liq.resize(N);
    dlnphi_drho_vap.resize(N);

    r.resize(N+2);
    negative_r.resize(N+2);
    J.resize(N+2, std::vector<long double>(N+2, 0));

    neg_dFdS.resize(N+2);
    dXdS.resize(N+2);

    // Fill the vector -dFdS with zeros (Gerg Eqn. 7.132)
    std::fill(neg_dFdS.begin(), neg_dFdS.end(), (double)0.0);
    // Last entry is 1
    neg_dFdS[N+1] = 1.0;
}
void SaturationSolvers::newton_raphson_VLE_GV::check_Jacobian(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO)
{
    // Reset all the variables and resize
    pre_call();
    std::size_t N = K.size();
    resize(N);

    SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
    SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid);
    SatV->specify_phase(iphase_gas);

    long double rhomolar_liq0 = IO.rhomolar_liq;
    const long double rhomolar_vap0 = IO.rhomolar_vap;
    long double T0 = IO.T;
    long double beta0 = IO.beta;

    // Build the Jacobian and residual vectors for given inputs of K_i,T,p
    build_arrays(HEOS,beta0,T0,rhomolar_liq0,rhomolar_vap0,z,K);
    
    // Make copies of the base
    std::vector<long double> r0 = r;
    STLMatrix J0 = J;
    STLMatrix Jnum = J;

    for (std::size_t i = 0; i < N+2; ++i)
    {
        for (std::size_t j = 0; j < N+2; ++j)
        {
            Jnum[i][j] = _HUGE;
        }    
    }

    for (std::size_t j = 0; j < N; ++j)
    {
        std::vector<long double> KK = K;
        KK[j] += 1e-6;
        build_arrays(HEOS,beta0,T0,rhomolar_liq0,rhomolar_vap0,z,KK);
        std::vector<long double> r1 = r;
        for (std::size_t i = 0; i < N+2; ++i)
        {
            Jnum[i][j] = (r1[i]-r0[i])/(log(KK[j])-log(K[j]));
        }
        std::cout << vec_to_string(get_col(Jnum,j),"%12.11f") << std::endl;
        std::cout << vec_to_string(get_col(J,j),"%12.11f") << std::endl;
    }
    
    build_arrays(HEOS,beta0,T0+1e-6,rhomolar_liq0,rhomolar_vap0,z,K);
    std::vector<long double> r1 = r, JN = r;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        Jnum[i][N] = (r1[i]-r0[i])/(log(T0+1e-6)-log(T0));
    }
    std::cout << vec_to_string(get_col(Jnum,N),"%12.11f") << std::endl;
    std::cout << vec_to_string(get_col(J,N),"%12.11f") << std::endl;

    // Build the Jacobian and residual vectors for given inputs of K_i,T,p
    build_arrays(HEOS,beta0,T0,rhomolar_liq0+1e-3,rhomolar_vap0,z,K);
    std::vector<long double> r2 = r, JNp1 = r;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        Jnum[i][N+1] = (r2[i]-r0[i])/(log(rhomolar_liq0+1e-3)-log(rhomolar_liq0));
    }
    std::cout << vec_to_string(get_col(Jnum, N+1),"%12.11f") << std::endl;
    std::cout << vec_to_string(get_col(J,N+1),"%12.11f") << std::endl;

    delete SatL; delete SatV;
}
void SaturationSolvers::newton_raphson_VLE_GV::call(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, mixture_VLE_IO &IO)
{
    int iter = 0;

    // Reset all the variables and resize
    pre_call();
    resize(K.size());

    SatL = new HelmholtzEOSMixtureBackend(HEOS->get_components()), 
    SatV = new HelmholtzEOSMixtureBackend(HEOS->get_components());
    SatL->specify_phase(iphase_liquid); // So it will always just use single-phase solution 
    SatV->specify_phase(iphase_gas); // So it will always just use single-phase solution

    do
    {
        // Build the Jacobian and residual vectors for given inputs of K_i,T,p
        build_arrays(HEOS,IO.beta,IO.T,IO.rhomolar_liq,IO.rhomolar_vap,z,K);

        // Solve for the step; v is the step with the contents 
        // [delta(lnK0), delta(lnK1), ..., delta(lnT), delta(lnrho')]
        std::vector<long double> v = linsolve(J, negative_r);

        max_rel_change = max_abs_value(v);

        // Set the variables again, the same structure independent of the specified variable
        for (unsigned int i = 0; i < N; i++)
        {
            K[i] = exp(log(K[i]) + v[i]);
            if (!ValidNumber(K[i]))
            {
                throw ValueError(format("K[i] (%g) is invalid",K[i]).c_str());
            }
        }
        IO.T = exp(log(IO.T) + v[N]);
        IO.rhomolar_liq = exp(log(IO.rhomolar_liq) + v[N+1]);

        if (fabs(IO.T) > 1e6)
        {
            /*std::cout << "J = " << vec_to_string(J,"%16.15g");
            std::cout << "nr = " << vec_to_string(r,"%16.15g");*/
            throw ValueError("Temperature or p has bad value");
        }
        
        //std::cout << iter << " " << T << " " << p << " " << error_rms << std::endl;
        iter++;
    }
    while(this->error_rms > 1e-8 && max_rel_change > 1000*LDBL_EPSILON && iter < IO.Nstep_max);
    Nsteps = iter;
    IO.p = p;
    IO.x = x; // Mole fractions in liquid
    IO.y = y; // Mole fractions in vapor
}

void SaturationSolvers::newton_raphson_VLE_GV::build_arrays(HelmholtzEOSMixtureBackend *HEOS, long double beta, long double T, long double rhomolar_liq, const long double rhomolar_vap, const std::vector<long double> &z, std::vector<long double> &K)
{
    // Step 0:
    // --------
    // Calculate the mole fractions in liquid and vapor phases
    x_and_y_from_K(beta, K, z, x, y);

    // Set the mole fractions in the classes
    SatL->set_mole_fractions(x);
    SatV->set_mole_fractions(y);

    // Update the liquid and vapor classes
    SatL->update(DmolarT_INPUTS, rhomolar_liq, T);
    SatV->update(DmolarT_INPUTS, rhomolar_vap, T);

    // For diagnostic purposes calculate the pressures (no derivatives are evaluated)
    long double p_liq = SatL->p();
    long double p_vap = SatV->p();
    p = 0.5*(p_liq+p_vap);

    // Step 2:
    // -------
    // Build the residual vector and the Jacobian matrix

    // For the residuals F_i
    for (unsigned int i = 0; i < N; ++i)
    {
        long double ln_phi_liq = SatL->mixderiv_ln_fugacity_coefficient(i);
        long double phi_iT_liq = SatL->mixderiv_dln_fugacity_coefficient_dT__constrho_n(i);
        dlnphi_drho_liq[i] = SatL->mixderiv_dln_fugacity_coefficient_drho__constT_n(i);
        for (unsigned int j = 0; j < N; ++j)
        {
            // I think this is wrong.
            phi_ij_liq[j] = SatL->mixderiv_ndln_fugacity_coefficient_dnj__constT_p(i,j) + (SatL->mixderiv_partial_molar_volume(i)/(SatL->gas_constant()*T)-1/p)*SatL->mixderiv_ndpdni__constT_V_nj(i); // 7.126 from GERG monograph
        }

        long double ln_phi_vap = SatV->mixderiv_ln_fugacity_coefficient(i);
        long double phi_iT_vap = SatV->mixderiv_dln_fugacity_coefficient_dT__constrho_n(i);
        dlnphi_drho_vap[i] = SatV->mixderiv_dln_fugacity_coefficient_drho__constT_n(i);
        for (unsigned int j = 0; j < N; ++j)
        {
            // I think this is wrong.
            phi_ij_vap[j] = SatV->mixderiv_ndln_fugacity_coefficient_dnj__constT_p(i,j) + (SatV->mixderiv_partial_molar_volume(i)/(SatV->gas_constant()*T)-1/p)*SatV->mixderiv_ndpdni__constT_V_nj(i); ; // 7.126 from GERG monograph
        }
        
        r[i] = log(K[i]) + ln_phi_vap - ln_phi_liq;
        // dF_i/d(ln(K_j))
        for (unsigned int j = 0; j < N; ++j)
        {	
            J[i][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2)*((1-beta)*phi_ij_vap[j]+beta*phi_ij_liq[j])+Kronecker_delta(i,j);
        }
        // dF_{i}/d(ln(T))
        J[i][N] = T*(phi_iT_vap-phi_iT_liq);
        // dF_{i}/d(ln(rho'))
        J[i][N+1] = -rhomolar_liq*dlnphi_drho_liq[i];
    }

    double summer1 = 0;
    for (unsigned int i = 0; i < N; ++i)
    {
        // Although the definition of this term is given by 
        // y[i]-x[i], when x and y are normalized, you get 
        // the wrong values.  Why? No idea.
        summer1 += z[i]*(K[i]-1)/(1-beta+beta*K[i]); 
    }
    r[N] = summer1;

    // For the residual term F_{N}, only non-zero derivatives are with respect
    // to ln(K[i])
    for (unsigned int j = 0; j < N; ++j)
    {
        J[N][j] = K[j]*z[j]/pow(1-beta+beta*K[j],(int)2);
    }

    // For the residual term F_{N+1} = p'-p''
    r[N+1] = p_liq-p_vap;
    for (unsigned int j = 0; j < N; ++j)
    {	
        J[N+1][j] = HEOS->gas_constant()*T*K[j]*z[j]/pow(1-beta+beta*K[j],(int)2)*((1-beta)*dlnphi_drho_vap[j]+beta*dlnphi_drho_liq[j]);
    }
    // dF_{N+1}/d(ln(T))
    J[N+1][N] = T*(SatL->mixderiv_dpdT__constV_n() - SatV->mixderiv_dpdT__constV_n());
    // dF_{N+1}/d(ln(rho'))
    J[N+1][N+1] = rhomolar_liq*SatL->mixderiv_dpdrho__constT_n();

    // Flip all the signs of the entries in the residual vector since we are solving Jv = -r, not Jv=r
    // Also calculate the rms error of the residual vector at this step
    error_rms = 0;
    for (unsigned int i = 0; i < N+2; ++i)
    {
        negative_r[i] = -r[i];
        error_rms += r[i]*r[i]; // Sum the squares
    }
    error_rms = sqrt(error_rms); // Square-root (The R in RMS)
}

void PhaseEnvelope::PhaseEnvelope_GV::build(HelmholtzEOSMixtureBackend *HEOS, const std::vector<long double> &z, std::vector<long double> &K, SaturationSolvers::mixture_VLE_IO &IO)
{
    // Use the residual function based on ln(K_i), ln(T) and ln(rho') as independent variables.  rho'' is specified
    SaturationSolvers::newton_raphson_VLE_GV NRVLE;
    SaturationSolvers::mixture_VLE_IO IO_NRVLE = IO;
    bubble.resize(z.size());
    dew.resize(z.size());
        
    // HACK
    IO_NRVLE.beta = 1.0;
    IO_NRVLE.Nstep_max = 30;
    int iter = 0;
    long double factor = IO_NRVLE.rhomolar_vap*0.01;
    for (;;)
    {
        if (iter > 0){ IO_NRVLE.rhomolar_vap += factor;}
        if (iter == 2 || (factor > 2 && factor < 0.24))
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnT,iter-2,iter-1,x));
            IO_NRVLE.rhomolar_liq = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnrhomolar_liq,iter-2,iter-1,x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(LinearInterp(dew.lnrhomolar_vap,dew.lnK[i],iter-2,iter-1,x));
            }
        }
        else if (iter == 3)
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnT,iter-3,iter-2,iter-1,x));
            IO_NRVLE.rhomolar_liq = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnrhomolar_liq,iter-3,iter-2,iter-1,x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(QuadInterp(dew.lnrhomolar_vap,dew.lnK[i],iter-3,iter-2,iter-1,x));
            }
        }
        else if (iter > 3)
        {
            long double x = log(IO_NRVLE.rhomolar_vap);
            IO_NRVLE.T = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnT, iter-4, iter-3, iter-2, iter-1, x));
            IO_NRVLE.rhomolar_liq = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnrhomolar_liq, iter-4, iter-3, iter-2, iter-1, x));
            for (std::size_t i = 0; i < K.size(); ++i)
            {
                K[i] = exp(CubicInterp(dew.lnrhomolar_vap, dew.lnK[i], iter-4, iter-3, iter-2, iter-1, x));
            }
        }
        /*if (IO_NRVLE.T > 344)
        {
            NRVLE.check_Jacobian(HEOS,z,K,IO_NRVLE);
        }*/
        NRVLE.call(HEOS,z,K,IO_NRVLE);
        dew.store_variables(IO_NRVLE.T,IO_NRVLE.p,IO_NRVLE.rhomolar_liq,IO_NRVLE.rhomolar_vap,K,IO_NRVLE.x,IO_NRVLE.y);
        iter ++;
        std::cout << format("%g %g %g %g %g %d %g\n",IO_NRVLE.p,IO_NRVLE.rhomolar_liq,IO_NRVLE.rhomolar_vap,IO_NRVLE.T,K[0],NRVLE.Nsteps,factor);
        if (iter < 5){continue;}
        if (NRVLE.Nsteps > 10)
        {
            factor /= 5;
        }
        else if (NRVLE.Nsteps > 5)
        {
            factor /= 1.2;
        }
        else if (NRVLE.Nsteps <= 4)
        {
            factor *= 1.2;
        }
        // Min step is 0.1 mol/m^3
        factor = std::max(factor,static_cast<long double>(0.1));
    }
}

} /* namespace CoolProp*/