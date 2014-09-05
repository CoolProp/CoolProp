#ifndef PHASEENVELOPE_H
#define PHASEENVELOPE_H

#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "PhaseEnvelopeRoutines.h"
#include "PhaseEnvelope.h"

namespace CoolProp{

void PhaseEnvelopeRoutines::build(HelmholtzEOSMixtureBackend &HEOS, const std::vector<long double> &z, SaturationSolvers::newton_raphson_saturation_options &IO_global)
{
    // Set some imput options
    SaturationSolvers::mixture_VLE_IO io;
    io.sstype = SaturationSolvers::imposed_p;
    io.Nstep_max = 2;

    // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
    long double Tguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

    // Use Wilson iteration to obtain updated guess for temperature
    Tguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);
    
    /*
    double mult_factor = 1.01;
    for (;IO.rhomolar_vap < 50000; IO.rhomolar_vap *= mult_factor)
    {
        if (IO.bubble_point){
            // Bubblepoint, vapor (y) is incipient phase
            NR.call(HEOS, IO.x, IO.y, IO);
        }
        else{
            
        }
        
        if (std::abs((IO.rhomolar_liq - IO.rhomolar_vap)/IO.rhomolar_liq) < 0.2){
            mult_factor = 1.0001;
        }
        else{
            mult_factor = 1.1;
        }
        
        
        
        std::cout << IO.p << " " << IO.T << " " << IO.rhomolar_liq << " " << IO.rhomolar_vap << std::endl;
    }
    */
        
    // Use the residual function based on ln(K_i), ln(T) and ln(rho') as independent variables.  rho'' is specified
    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    
    IO.bubble_point = false; // Do a "dewpoint" calculation all the way around
    IO.x = io.x;
    IO.y = io.y;
    IO.rhomolar_liq = io.rhomolar_liq;
    IO.rhomolar_vap = io.rhomolar_vap;
    IO.T = io.T;
    IO.p = io.p;
    IO.Nstep_max = 100;
    
    PhaseEnvelopeData &env = HEOS.PhaseEnvelope;

    int iter = 0;
    long double factor = 1.05;
    for (;;)
    {
        long double x = IO.rhomolar_vap;
        if (iter > 0){ IO.rhomolar_vap *= factor;}
        if (iter == 2)
        {
            IO.T = LinearInterp(env.rhomolar_vap,env.T,iter-2,iter-1,x);
            IO.rhomolar_liq = LinearInterp(env.rhomolar_vap,env.rhomolar_liq,iter-2,iter-1,x);
            for (std::size_t i = 0; i < IO.x.size(); ++i)
            {
                IO.x[i] = LinearInterp(env.rhomolar_vap, env.x[i], iter-2, iter-1, x);
            }
        }
        else if (iter == 3)
        {
            IO.T = QuadInterp(env.rhomolar_vap, env.T, iter-3, iter-2, iter-1, x);
            IO.rhomolar_liq = QuadInterp(env.rhomolar_vap, env.rhomolar_liq, iter-3, iter-2, iter-1, x);
            for (std::size_t i = 0; i < IO.x.size(); ++i)
            {
                IO.x[i] = QuadInterp(env.rhomolar_vap, env.x[i], iter-3, iter-2, iter-1, x);
            }
        }
        else if (iter > 3)
        {
            IO.T = CubicInterp(env.rhomolar_vap, env.T, iter-4, iter-3, iter-2, iter-1, x);
            IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iter-4, iter-3, iter-2, iter-1, x);
            for (std::size_t i = 0; i < IO.x.size(); ++i)
            {
                IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], iter-4, iter-3, iter-2, iter-1, x);
            }
        }
        // The last mole fraction is sum of N-1 first elements
        IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
        
        // Dewpoint, liquid (x) is incipient phase
        NR.call(HEOS, IO.y, IO.x, IO);
        
        /*
        // Critical point jump
        if (std::abs((IO.rhomolar_liq - IO.rhomolar_vap)/IO.rhomolar_liq) < 0.05 && IO.rhomolar_liq  > IO.rhomolar_vap){
            long double rhoc_approx = 0.5*IO.rhomolar_liq + 0.5*IO.rhomolar_vap;
            long double rho_vap_new = 2*rhoc_approx - IO.rhomolar_vap;
            IO.rhomolar_liq = 2*rhoc_approx - IO.rhomolar_liq;
            mult_factor = rho_vap_new/IO.rhomolar_vap;
            //std::cout << "dv" << rho_vap_new << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
            for (std::size_t i = 0; i < IO.x.size()-1; ++i){
                IO.x[i] = IO.x[i] - 2*(IO.x[i] - IO.y[i]);
            }
            IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
            //std::cout << "dv" << rho_vap_new << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
        }
        */
        env.store_variables(IO.T, IO.p, IO.rhomolar_liq, IO.rhomolar_vap, IO.x, IO.y);
        iter ++;
        if (iter < 5){continue;}
        if (IO.Nsteps > 10)
        {
            factor /= 5;
        }
        else if (IO.Nsteps > 5)
        {
            factor /= 1.2;
        }
        else if (IO.Nsteps <= 4)
        {
            factor *= 1.2;
        }
        // Min step is 0.1 mol/m^3
        factor = std::max(factor,static_cast<long double>(0.1));
    }
}

} /* namespace CoolProp */

#endif