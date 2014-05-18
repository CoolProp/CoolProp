
#include "TransportRoutines.h"

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
        throw NotImplementedError("TransportRoutines::GeneralDiluteGasViscosity is only for pure and pseudo-pure");
    }
}

}; /* namespace CoolProp */