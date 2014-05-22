
#include "TransportRoutines.h"
#include "CoolPropFluid.h"

namespace CoolProp{

long double TransportRoutines::viscosity_dilute_kinetic_theory(HelmholtzEOSMixtureBackend &HEOS)
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
        throw NotImplementedError("TransportRoutines::viscosity_dilute_kinetic_theory is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::viscosity_dilute_collision_integral(HelmholtzEOSMixtureBackend &HEOS)
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
        throw NotImplementedError("TransportRoutines::viscosity_dilute_collision_integral is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::viscosity_dilute_powers_of_T(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteGasPowersOfT &data = HEOS.components[0]->transport.viscosity_dilute.powers_of_T;
        const std::vector<long double> &a = data.a, &t = data.t;

        long double summer = 0, T = HEOS.T();
        for (std::size_t i = 0; i < a.size(); ++i)
        {
            summer += a[i]*pow(T, t[i]);
        }
        return summer;
    }
    else{
        throw NotImplementedError("TransportRoutines::viscosity_dilute_powers_of_T is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::viscosity_dilute_collision_integral_powers_of_T(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        // Retrieve values from the state class
        CoolProp::ViscosityDiluteCollisionIntegralPowersOfTstarData &data = HEOS.components[0]->transport.viscosity_dilute.collision_integral_powers_of_Tstar;
        const std::vector<long double> &a = data.a, &t = data.t;

        long double summer = 0, Tstar = HEOS.T()/data.T_reducing;
        for (std::size_t i = 0; i < a.size(); ++i)
        {
            summer += a[i]*pow(Tstar, t[i]);
        }
        return data.C*sqrt(HEOS.T())/summer;
    }
    else{
        throw NotImplementedError("TransportRoutines::viscosity_dilute_collision_integral_powers_of_T is only for pure and pseudo-pure");
    }

}
long double TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        CoolProp::ViscosityModifiedBatschinskiHildebrandData &HO = HEOS.components[0]->transport.viscosity_higher_order.modified_Batschinski_Hildebrand;

        long double delta = HEOS.rhomolar()/HO.rhomolar_reduce, tau = HO.T_reduce/HEOS.T();

        // The first term that is formed of powers of tau (Tc/T) and delta (rho/rhoc)
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

        // The higher-order-term component
        return S + F*(1/(delta0-delta)-1/delta0); // Pa-s
    }
    else{
        throw NotImplementedError("TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand is only for pure and pseudo-pure");
    }
}

long double TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend(HelmholtzEOSMixtureBackend &HEOS)
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
        throw NotImplementedError("TransportRoutines::viscosity_initial_density_dependence_Rainwater_Friend is only for pure and pseudo-pure");
    }
}

static void visc_Helper(double Tbar, double rhobar, double *mubar_0, double *mubar_1)
{
	std::vector<std::vector<long double> > H(6,std::vector<long double>(7,0));
    double sum;
	int i,j;

	// Dilute-gas component
	*mubar_0=100.0*sqrt(Tbar)/(1.67752+2.20462/Tbar+0.6366564/powInt(Tbar,2)-0.241605/powInt(Tbar,3));

	//Fill in zeros in H
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			H[i][j]=0;
		}
	}

	//Set non-zero parameters of H
	H[0][0]=5.20094e-1;
	H[1][0]=8.50895e-2;
	H[2][0]=-1.08374;
	H[3][0]=-2.89555e-1;

	H[0][1]=2.22531e-1;
	H[1][1]=9.99115e-1;
	H[2][1]=1.88797;
	H[3][1]=1.26613;
	H[5][1]=1.20573e-1;

	H[0][2]=-2.81378e-1;
	H[1][2]=-9.06851e-1;
	H[2][2]=-7.72479e-1;
	H[3][2]=-4.89837e-1;
	H[4][2]=-2.57040e-1;

	H[0][3]=1.61913e-1;
	H[1][3]=2.57399e-1;

	H[0][4]=-3.25372e-2;
	H[3][4]=6.98452e-2;

	H[4][5]=8.72102e-3;

	H[3][6]=-4.35673e-3;
	H[5][6]=-5.93264e-4;

	// Finite density component
	sum=0;
	for (i=0;i<=5;i++)
	{
		for (j=0;j<=6;j++)
		{
			sum+=powInt(1/Tbar-1,i)*(H[i][j]*powInt(rhobar-1,j));
		}
	}
	*mubar_1=exp(rhobar*sum);
}
long double TransportRoutines::viscosity_water_hardcoded(HelmholtzEOSMixtureBackend &HEOS)
{
	double x_mu=0.068,qc=1/1.9,qd=1/1.1,nu=0.630,gamma=1.239,zeta_0=0.13,LAMBDA_0=0.06,Tbar_R=1.5, pstar, Tstar, rhostar;
	double delta,tau,mubar_0,mubar_1,mubar_2,drhodp,drhodp_R,DeltaChibar,zeta,w,L,Y,psi_D,Tbar,rhobar;
	double drhobar_dpbar,drhobar_dpbar_R,R_Water;
	
    pstar = 22.064e6;
    Tstar = 647.096;
    rhostar = 322;
	Tbar = HEOS.T()/Tstar;
    rhobar = HEOS.keyed_output(CoolProp::iDmass)/rhostar;
    R_Water = HEOS.gas_constant()/HEOS.molar_mass();

	// Dilute and finite gas portions
	visc_Helper(Tbar, rhobar, &mubar_0, &mubar_1);

	///************************ Critical Enhancement ************************
	delta=rhobar;
	// "Normal" calculation
	tau=1/Tbar;
	drhodp=1/(R_Water*HEOS.T()*(1+2*delta*HEOS.dalphar_dDelta()+delta*delta*HEOS.d2alphar_dDelta2()));
	drhobar_dpbar = pstar/rhostar*drhodp;
	// "Reducing" calculation
	tau=1/Tbar_R;
	drhodp_R=1/(R_Water*Tbar_R*Tstar*(1+2*delta*HEOS.dalphar_dDelta()+delta*delta*HEOS.d2alphar_dDelta2()));
	drhobar_dpbar_R = pstar/rhostar*drhodp_R;
	
	DeltaChibar=rhobar*(drhobar_dpbar-drhobar_dpbar_R*Tbar_R/Tbar);
	if (DeltaChibar<0)
		DeltaChibar=0;
	zeta=zeta_0*pow(DeltaChibar/LAMBDA_0,nu/gamma);
	if (zeta<0.3817016416){
		Y=1.0/5.0*qc*zeta*powInt(qd*zeta,5)*(1-qc*zeta+powInt(qc*zeta,2)-765.0/504.0*powInt(qd*zeta,2));
	}
	else
	{
		psi_D=acos(pow(1+powInt(qd*zeta,2),-1.0/2.0));
		w=sqrt(fabs((qc*zeta-1)/(qc*zeta+1)))*tan(psi_D/2.0);
		if (qc*zeta>1){
			L=log((1+w)/(1-w));
		}
		else{
			L=2*atan(fabs(w));
		}
		Y=1.0/12.0*sin(3*psi_D)-1/(4*qc*zeta)*sin(2*psi_D)+1.0/powInt(qc*zeta,2)*(1-5.0/4.0*powInt(qc*zeta,2))*sin(psi_D)-1.0/powInt(qc*zeta,3)*((1-3.0/2.0*powInt(qc*zeta,2))*psi_D-pow(fabs(powInt(qc*zeta,2)-1),3.0/2.0)*L);
	}
	mubar_2=exp(x_mu*Y);

	return (mubar_0*mubar_1*mubar_2)/1e6;
}

long double TransportRoutines::viscosity_hydrogen_higher_order_hardcoded(HelmholtzEOSMixtureBackend &HEOS)
{
    long double Tr = HEOS.T()/33.145, rhor = HEOS.keyed_output(CoolProp::iDmass)/90.5;
    long double c[] = {0, 6.43449673e-6, 4.56334068e-2, 2.32797868e-1, 9.58326120e-1, 1.27941189e-1, 3.63576595e-1};
    return c[1]*pow(rhor,2)*exp(c[2]*Tr+c[3]/Tr+c[4]*pow(rhor,2)/(c[5]+Tr)+c[6]*pow(rhor,6));
}

long double TransportRoutines::viscosity_hexane_higher_order_hardcoded(HelmholtzEOSMixtureBackend &HEOS)
{
    long double Tr = HEOS.T()/507.82, rhor = HEOS.keyed_output(CoolProp::iDmass)/233.182;

    // Output is in Pa-s
    double c[] = {2.53402335/1e6, -9.724061002/1e6, 0.469437316, 158.5571631, 72.42916856/1e6, 10.60751253, 8.628373915, -6.61346441, -2.212724566};
	return pow(rhor,static_cast<long double>(2.0/3.0))*sqrt(Tr)*(c[0]/Tr+c[1]/(c[2]+Tr+c[3]*rhor*rhor)+c[4]*(1+rhor)/(c[5]+c[6]*Tr+c[7]*rhor+rhor*rhor+c[8]*rhor*Tr));
}

long double TransportRoutines::viscosity_higher_order_friction_theory(HelmholtzEOSMixtureBackend &HEOS)
{
    if (HEOS.is_pure_or_pseudopure)
    {
        CoolProp::ViscosityFrictionTheoryData &F = HEOS.components[0]->transport.viscosity_higher_order.friction_theory;

        long double tau = F.T_reduce/HEOS.T(), kii, krrr, kaaa, krr, kdrdr;

        double psi1 = exp(tau)-F.c1;
	    double psi2 = exp(pow(tau,2))-F.c2;

	    double ki = (F.Ai[0] + F.Ai[1]*psi1 + F.Ai[2]*psi2)*tau;

        double ka = (F.Aa[0] + F.Aa[1]*psi1 + F.Aa[2]*psi2)*pow(tau, F.Na);
	    double kr = (F.Ar[0] + F.Ar[1]*psi1 + F.Ar[2]*psi2)*pow(tau, F.Nr);
	    double kaa = (F.Aaa[0] + F.Aaa[1]*psi1 + F.Aaa[2]*psi2)*pow(tau, F.Naa);
        if (F.Arr.empty()){
            krr = 0;
            kdrdr = (F.Adrdr[0] + F.Adrdr[1]*psi1 + F.Adrdr[2]*psi2)*pow(tau, F.Nrr);
        }
        else{
	        krr = (F.Arr[0] + F.Arr[1]*psi1 + F.Arr[2]*psi2)*pow(tau, F.Nrr);
            kdrdr = 0;
        }

        if (!F.Aii.empty() && !F.Arrr.empty() && !F.Aaaa.empty()){
	        kii = (F.Aii[0] + F.Aii[1]*psi1 + F.Aii[2]*psi2)*pow(tau, F.Nii);
	        krrr = (F.Arrr[0] + F.Arrr[1]*psi1 + F.Arrr[2]*psi2)*pow(tau, F.Nrrr);
	        kaaa = (F.Aaaa[0] + F.Aaaa[1]*psi1 + F.Aaaa[2]*psi2)*pow(tau, F.Naaa);
        }
        else{
            kii = 0;
            krrr = 0;
            kaaa = 0;
        }

        double p = HEOS.p()/1e5; // [bar]; 1e5 for conversion from Pa -> bar
	    double pr = HEOS.T()*HEOS.first_partial_deriv(CoolProp::iP, CoolProp::iT, CoolProp::iDmolar)/1e5; // [bar/K]; 1e5 for conversion from Pa -> bar
	    double pa = p - pr; //[bar]
        double pid = HEOS.rhomolar() * HEOS.gas_constant() * HEOS.T() / 1e5; // [bar]; 1e5 for conversion from Pa -> bar
	    double deltapr = pr - pid;

	    double eta_f = ka*pa + kr*deltapr + ki*pid + kaa*pa*pa + kdrdr*deltapr*deltapr + krr*pr*pr + kii*pid*pid + krrr*pr*pr*pr + kaaa*pa*pa*pa;
        
        return eta_f; //[Pa-s]
    }
    else{
        throw NotImplementedError("TransportRoutines::viscosity_higher_order_friction_theory is only for pure and pseudo-pure");
    }


}

}; /* namespace CoolProp */