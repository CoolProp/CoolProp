
#include "TransportRoutines.h"
#include "CoolPropFluid.h"

namespace CoolProp{

long double TransportRoutines::general_dilute_gas_viscosity(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        long double Tstar = HEOS.T()/HEOS.components[0]->transport.epsilon_over_k;
        long double sigma_nm = HEOS.components[0]->transport.sigma_eta*1e9; // 1e9 to convert from m to nm
        long double molar_mass_kgkmol = HEOS.molar_mass()*1000; // 1000 to convert from kg/mol to kg/kmol
    
        // The nondimensional empirical collision integral from Neufeld 
        // Neufeld, P. D.; Janzen, A. R.; Aziz, R. A. Empirical Equations to Calculate 16 of the Transport Collision Integrals (l,s)* 
        // for the Lennard-Jones (12-6) Potential. J. Chem. Phys. 1972, 57, 1100-1102
        long double OMEGA22 = 1.16145*pow(Tstar, static_cast<long double>(-0.14874))+0.52487*exp(-0.77320*Tstar)+2.16178*exp(-2.43787*Tstar);

        // The dilute gas component - 
        return 26.692e-9*sqrt(molar_mass_kgkmol*HEOS.T())/(pow(sigma_nm, 2)*OMEGA22); // Pa-s
    }
    else{
        throw NotImplementedError("TransportRoutines::general_dilute_gas_viscosity is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::dilute_gas_viscosity(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteGasCollisionIntegralData &data = HEOS.components[0]->transport.viscosity_dilute.collision_integral;
        const std::vector<long double> &a = data.a, &t = data.t;
        const long double C = data.C, molar_mass = data.molar_mass;

        long double S;
        // Unit conversions and variable definitions
        const long double Tstar = HEOS.T()/HEOS.components[0]->transport.epsilon_over_k;
        const long double sigma_nm = HEOS.components[0]->transport.sigma_eta*1e9; // 1e9 to convert from m to nm
        const long double molar_mass_kgkmol = molar_mass*1000; // 1000 to convert from kg/mol to kg/kmol

        /// Both the collision integral \f$\mathfrak{S}^*\f$ and effective cross section \f$\Omega^{(2,2)}\f$ have the same form, 
        /// in general we don't care which is used.  The are related through \f$\Omega^{(2,2)} = (5/4)\mathfrak{S}^*\f$ 
        /// see Vesovic(JPCRD, 1990) for CO\f$_2\f$ for further information
        long double summer = 0, lnTstar = log(Tstar);
        for (std::size_t i = 0; i < a.size(); ++i)
        {
            summer += a[i]*pow(lnTstar,t[i]);
        }
        S = exp(summer);
        
        // The dilute gas component
        return C*sqrt(molar_mass_kgkmol*HEOS.T())/(pow(sigma_nm, 2)*S); // Pa-s
    }
    else{
        throw NotImplementedError("TransportRoutines::dilute_gas_viscosity is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::modified_Batschinski_Hildebrand_viscosity_term(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        CoolProp::ViscosityModifiedBatschinskiHildebrandData &HO = HEOS.components[0]->transport.viscosity_higher_order.modified_Batschinski_Hildebrand;

        long double delta = HEOS.rhomolar()/HO.rhomolar_reduce, tau = HO.T_reduce/HEOS.T();

        // The first term that is formed of powers of tau and delta
        long double S = 0;
        for (unsigned int i = 0; i < HO.a.size(); ++i){
            S += HO.a[i]*pow(delta, HO.d1[i])*pow(tau, HO.t1[i])*exp(HO.gamma[i]*pow(delta, HO.l[i]));
        }

        // For the terms that multiplies the bracketed term with delta and delta0
        long double F = 0;
        for (unsigned int i = 0; i < HO.f.size(); ++i){
            F += HO.f[i]*pow(delta, HO.d2[i])*pow(tau, HO.t2[i]);
        }

        // for delta_0
        long double summer_numer = 0;
        for (unsigned int i = 0; i < HO.g.size(); ++i){
            summer_numer += HO.g[i]*pow(tau, HO.h[i]);
        }
        long double summer_denom = 0;
        for (unsigned int i = 0; i < HO.p.size(); ++i){
            summer_denom += HO.p[i]*pow(tau, HO.q[i]);
        }
        long double delta0 = summer_numer/summer_denom;

        double Tr = 1/tau;
        // The higher-order-term component
        return S + F*(1/(delta0-delta)-1/delta0); // Pa-s
    }
    else{
        throw NotImplementedError("TransportRoutines::modified_Batschinski_Hildebrand_viscosity_term is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::initial_density_dependence_viscosity_term(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        // Retrieve values from the state class
        CoolProp::ViscosityRainWaterFriendData &data = HEOS.components[0]->transport.viscosity_initial.rainwater_friend;
        const std::vector<long double> &b = data.b, &t = data.t;

        long double B_eta, B_eta_star;
        long double Tstar = HEOS.T()/HEOS.components[0]->transport.epsilon_over_k; // [no units]
        long double sigma = HEOS.components[0]->transport.sigma_eta; // [m]
        
        long double summer = 0;
        for (unsigned int i = 0; i < b.size(); ++i){
            summer += b[i]*pow(Tstar, t[i]);
        }
        B_eta_star = summer; // [no units]
        B_eta = 6.02214129e23*pow(sigma, 3)*B_eta_star; // [m^3/mol]
        return B_eta; // [m^3/mol]
    }
    else{
        throw NotImplementedError("TransportRoutines::initial_density_dependence_viscosity_term is only for pure and pseudo-pure");
    }
}


}; /* namespace CoolProp */