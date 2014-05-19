
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



long double TransportRoutines::dilute_gas_viscosity(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        // HARD CODED FOR PROPANE FOR TESTING PURPOSES
        double _a[]={0.25104574,-0.47271238,0,0.060836515}, _t[] = {0,1,2,3}, C = 0.141824e-6/sqrt(44.0956), D = sqrt(44.0956)*C;
        std::vector<long double> a(_a,_a+sizeof(_a)/sizeof(_a[0]));
        std::vector<long double> t(_t,_t+sizeof(_t)/sizeof(_t[0]));
        HEOS.components[0]->transport.epsilon_over_k = 263.88;
        HEOS.components[0]->transport.sigma_eta = 0.49748e-9;
        long double molar_mass = 44.0956/1000;

        long double S;
        long double Tstar = HEOS.T()/HEOS.components[0]->transport.epsilon_over_k;
        long double sigma_nm = HEOS.components[0]->transport.sigma_eta*1e9; // 1e9 to convert from m to nm
        long double molar_mass_kgkmol = molar_mass*1000; // 1000 to convert from kg/mol to kg/kmol

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
        //// HARD CODED FOR PROPANE FOR TESTING PURPOSES
        double eta_star, _a[]={0.25104574,-0.47271238,0,0.060836515,0}, _t[] = {0,1,2,3,4}, C = 0.021357e-6;
        std::vector<long double> a(_a,_a+sizeof(_a)/sizeof(_a[0]));
        std::vector<long double> t(_t,_t+sizeof(_t)/sizeof(_t[0]));
        HEOS.components[0]->transport.epsilon_over_k = 263.88;
        HEOS.components[0]->transport.sigma_eta = 0.49748e-9;
        std::vector< std::vector<long double> > e(6, std::vector<long double>(3,0));
	    e[2][0] =  35.9873030195*1e-6; e[2][1] = -180.512188564*1e-6; e[2][2] =  87.7124888223*1e-6;
	    e[3][0] = -105.773052525*1e-6; e[3][1] =  205.319740877*1e-6; e[3][2] = -129.210932610*1e-6;
	    e[4][0] =  58.9491587759*1e-6; e[4][1] = -129.740033100*1e-6; e[4][2] =  76.6280419971*1e-6;
	    e[5][0] = -9.59407868475*1e-6; e[5][1] =  21.0726986598*1e-6; e[5][2] = -14.3971968187*1e-6;
        double _f[] = {1616.88405374e-6},summer_num,summer_denom;
        std::vector<long double> f(_f,_f+sizeof(_f)/sizeof(_f[0]));
        std::vector<long double> g(2,1),h(2,0);
        g[0] = 2.50053938863; g[1] = 2.50053938863*0.860516059264; h[0] = 0; h[1] = -0.5;
        std::vector<long double> p(1,1),q(1,0);
        long double delta = HEOS.rhomolar()/5000.0, tau = 369.825/HEOS.T();
        //// END HARD CODED
        
        long double summer = 0;
        for (int i = 0; i < e.size(); ++i){
            for (int j = 0; j < e[0].size(); ++j){
                summer += e[i][j]*pow(delta,i)*pow(tau,j);
            }
        }
        long double S = summer;

        summer = 0;
        for (int i = 0; i < f.size(); ++i){
            summer += f[i]/pow(tau, i);
        }
        long double F = summer;

        // for delta_0
        summer_num = 0;
        for (int i = 0; i < g.size(); ++i){
            summer_num += g[i]*pow(tau, h[i]);
        }
        summer_denom = 0;
        for (int i = 0; i < p.size(); ++i){
            summer_denom += p[i]*pow(tau, q[i]);
        }
        long double delta0 = summer_num/summer_denom;

        // The higher-order-term component
        return S+F*delta*(1/(delta0-delta)-1/delta0); // Pa-s
    }
    else{
        throw NotImplementedError("TransportRoutines::GeneralDiluteGasViscosity is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::initial_density_dependence_viscosity_term(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        //// HARD CODED FOR PROPANE FOR TESTING PURPOSES
        double _b[]={-19.572881,219.73999,-1015.3226,2471.01251,-3375.1717,2491.6597,-787.26086,14.085455,-0.34664158};
        std::vector<long double> b(_b,_b+sizeof(_b)/sizeof(_b[0]));
        double _c[]={0,-0.25,-0.5,-0.75,-1.0,-1.25,-1.5,-2.5,-5.5};
        std::vector<long double> c(_c,_c+sizeof(_c)/sizeof(_c[0]));
        HEOS.components[0]->transport.epsilon_over_k = 263.88;
        HEOS.components[0]->transport.sigma_eta = 0.49748e-9;
        //// END HARD CODED

        long double B_eta, B_eta_star;
        long double Tstar = HEOS.T()/HEOS.components[0]->transport.epsilon_over_k; // [no units]
        long double sigma = HEOS.components[0]->transport.sigma_eta; // [m]
        
        long double summer = 0;
        for (int i = 0; i < b.size(); ++i){
            summer += b[i]*pow(Tstar, c[i]);
        }
        B_eta_star = summer; // [no units]
        B_eta = 6.02214129e23*pow(sigma,3)*B_eta_star; // [m^3/mol]
        return B_eta; // [m^3/mol]
    }
    else{
        throw NotImplementedError("TransportRoutines::GeneralDiluteGasViscosity is only for pure and pseudo-pure");
    }
}


}; /* namespace CoolProp */