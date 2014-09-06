#ifndef PHASEENVELOPE_H
#define PHASEENVELOPE_H

#include "HelmholtzEOSMixtureBackend.h"
#include "VLERoutines.h"
#include "PhaseEnvelopeRoutines.h"
#include "PhaseEnvelope.h"

namespace CoolProp{

void PhaseEnvelopeRoutines::build(HelmholtzEOSMixtureBackend &HEOS)
{
    // Set some imput options
    SaturationSolvers::mixture_VLE_IO io;
    io.sstype = SaturationSolvers::imposed_p;
    io.Nstep_max = 2;
    
    // Set the pressure to a low pressure 
    HEOS._p = 100000;
    
    // Get an extremely rough guess by interpolation of ln(p) v. T curve where the limits are mole-fraction-weighted
    long double Tguess = SaturationSolvers::saturation_preconditioner(HEOS, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions);

    // Use Wilson iteration to obtain updated guess for temperature
    Tguess = SaturationSolvers::saturation_Wilson(HEOS, HEOS._Q, HEOS._p, SaturationSolvers::imposed_p, HEOS.mole_fractions, Tguess);
    
    // Actually call the successive substitution solver
    io.beta = 1;
    HEOS._Q = 1;
    SaturationSolvers::successive_substitution(HEOS, HEOS._Q, Tguess, HEOS._p, HEOS.mole_fractions, HEOS.K, io);
        
    // Use the residual function based on x_i, T and rho' as independent variables.  rho'' is specified
    SaturationSolvers::newton_raphson_saturation NR;
    SaturationSolvers::newton_raphson_saturation_options IO;
    
    IO.bubble_point = false; // Do a "dewpoint" calculation all the way around
    IO.x = io.x;
    IO.y = HEOS.mole_fractions;
    IO.rhomolar_liq = io.rhomolar_liq;
    IO.rhomolar_vap = io.rhomolar_vap;
    IO.T = io.T;
    IO.p = io.p;
    IO.Nstep_max = 100;
    
    bool dont_extrapolate = false;
    
    PhaseEnvelopeData &env = HEOS.PhaseEnvelope;
    env.resize(HEOS.mole_fractions.size());

    std::size_t iter = 0, //< The iteration counter
                iter0 = 0; //< A reference point for the counter, can be increased to go back to linear interpolation
    long double factor = 1.05;
    for (;;)
    {
        long double x = IO.rhomolar_vap;
        if (iter - iter0 > 0){ IO.rhomolar_vap *= factor;}
        if (dont_extrapolate)
        {
            // Reset the step to a reasonably small size
            factor = 1.0001;
        }
        else if (iter - iter0 == 2)
        {
            IO.T = LinearInterp(env.rhomolar_vap,env.T,iter-2,iter-1,x);
            IO.rhomolar_liq = LinearInterp(env.rhomolar_vap,env.rhomolar_liq,iter-2,iter-1,x);
            for (std::size_t i = 0; i < IO.x.size()-1; ++i) // First N-1 elements
            {            
                IO.x[i] = LinearInterp(env.rhomolar_vap, env.x[i], iter-2, iter-1, x);
            }
        }
        else if (iter - iter0 == 3)
        {
            IO.T = QuadInterp(env.rhomolar_vap, env.T, iter-3, iter-2, iter-1, x);
            IO.rhomolar_liq = QuadInterp(env.rhomolar_vap, env.rhomolar_liq, iter-3, iter-2, iter-1, x);
            for (std::size_t i = 0; i < IO.x.size()-1; ++i) // First N-1 elements
            {
                IO.x[i] = QuadInterp(env.rhomolar_vap, env.x[i], iter-3, iter-2, iter-1, x);
            }
        }
        else if (iter - iter0 > 3)
        {
            IO.T = CubicInterp(env.rhomolar_vap, env.T, iter-4, iter-3, iter-2, iter-1, x);
            IO.rhomolar_liq = CubicInterp(env.rhomolar_vap, env.rhomolar_liq, iter-4, iter-3, iter-2, iter-1, x);
            for (std::size_t i = 0; i < IO.x.size()-1; ++i) // First N-1 elements
            {
                IO.x[i] = CubicInterp(env.rhomolar_vap, env.x[i], iter-4, iter-3, iter-2, iter-1, x);
            }
        }
        // The last mole fraction is sum of N-1 first elements
        IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
        
        // Dewpoint calculation, liquid (x) is incipient phase
        NR.call(HEOS, IO.y, IO.x, IO);
        
        std::cout << "dv " << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " T " << IO.T << " x " << vec_to_string(IO.x, "%0.10Lg") << " y " << vec_to_string(IO.y, "%0.10Lg") << " Ns " << IO.Nsteps << std::endl;
        env.store_variables(IO.T, IO.p, IO.rhomolar_liq, IO.rhomolar_vap, IO.x, IO.y);
        
        iter ++;
        
        // Critical point jump
        if (std::abs((IO.rhomolar_liq - IO.rhomolar_vap)/IO.rhomolar_liq) < 0.05 && IO.rhomolar_liq  > IO.rhomolar_vap){
            //std::cout << "dv" << IO.rhomolar_vap << " dl " << IO.rhomolar_liq << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
            long double rhoc_approx = 0.5*IO.rhomolar_liq + 0.5*IO.rhomolar_vap;
            long double rho_vap_new = 2*rhoc_approx - IO.rhomolar_vap;
            // Linearly interpolate to get new guess for T
            IO.T = LinearInterp(env.rhomolar_vap,env.T,iter-2,iter-1,rho_vap_new);
            IO.rhomolar_liq = 2*rhoc_approx - IO.rhomolar_liq;
            for (std::size_t i = 0; i < IO.x.size()-1; ++i){
                IO.x[i] = 2*IO.y[i] - IO.x[i];
            }
            IO.x[IO.x.size()-1] = 1 - std::accumulate(IO.x.begin(), IO.x.end()-1, 0.0);
            factor = rho_vap_new/IO.rhomolar_vap;
            dont_extrapolate = true; // So that we use the mole fractions we calculated here instead of the extrapolated values
            //std::cout << "dv " << rho_vap_new << " dl " << IO.rhomolar_liq << " " << vec_to_string(IO.x, "%0.10Lg") << " " << vec_to_string(IO.y, "%0.10Lg") << std::endl;
            iter0 = iter - 1; // Back to linear interpolation again
            continue;
        }
        
        dont_extrapolate = false;
        if (iter < 5){continue;}
        if (IO.Nsteps > 10)
        {
            factor = 1 + (factor-1)/10;
        }
        else if (IO.Nsteps > 5)
        {
            factor = 1 + (factor-1)/3;
        }
        else if (IO.Nsteps <= 4)
        {
            factor = 1 + (factor-1)*2;
        }
        // Min step is 1.01
        factor = std::max(factor, static_cast<long double>(1.01));
    }
}

} /* namespace CoolProp */

#endif