#include "CubicBackend.h"

void CoolProp::AbstractCubicBackend::setup(){
	// Reset the residual Helmholtz energy class
	residual_helmholtz.reset(new CubicResidualHelmholtz(this));
	// If pure, set the mole fractions to be unity
	if (is_pure_or_pseudopure){
		mole_fractions = std::vector<CoolPropDbl>(1, 1.0);
	}
	// Now set the reducing function for the mixture
    Reducing.reset(new ConstantReducingFunction(cubic->T_r, cubic->rho_r));
}

void CoolProp::AbstractCubicBackend::get_linear_reducing_parameters(double &rhomolar_r, double &T_r){
    // In the case of models where the reducing temperature is not a function of composition (SRK, PR, etc.), 
    // we need to use an appropriate value for T_r and v_r, here we use a linear weighting
    T_r = 0; 
    double v_r = 0;
    const std::vector<double> Tc = cubic->get_Tc(), pc = cubic->get_pc();
    for (std::size_t i = 0; i < mole_fractions.size(); ++i){   
        T_r += mole_fractions[i]*Tc[i];
        // Curve fit from all the pure fluids in CoolProp (thanks to recommendation of A. Kazakov)
        double v_c_Lmol = 2.14107171795*(Tc[i]/pc[i]*1000)+0.00773144012514; // [L/mol]
        v_r += mole_fractions[i]*v_c_Lmol/1000.0;
    }
    rhomolar_r = 1/v_r;
}

void CoolProp::AbstractCubicBackend::get_critical_point_search_radii(double &R_delta, double &R_tau)
{
    // Get the starting values from the base class
    CoolProp::HelmholtzEOSMixtureBackend::get_critical_point_search_radii(R_delta, R_tau);

    // Now we scale them to get the appropriate search radii
    double Tr_GERGlike, rhor_GERGlike;
    get_linear_reducing_parameters(rhor_GERGlike, Tr_GERGlike);
    R_delta *= rhor_GERGlike/rhomolar_reducing();
    R_tau *= T_reducing()/Tr_GERGlike;
}

bool CoolProp::AbstractCubicBackend::get_critical_is_terminated(double &delta, double &tau)
{
    // If the volume is less than the mixture covolume, stop.  The mixture covolume is the
    // smallest volume that is physically allowed for a cubic EOS
    double b = get_cubic()->bm_term(mole_fractions); // [m^3/mol]
    double v = 1/(delta*rhomolar_reducing()); //[m^3/mol]
    bool covolume_check = v < 1.2*b;
    
    return covolume_check;
}

void CoolProp::AbstractCubicBackend::get_critical_point_starting_values(double &delta0, double &tau0){
    
    // Get the starting values from the base class
    CoolProp::HelmholtzEOSMixtureBackend::get_critical_point_starting_values(delta0, tau0);

    // The starting tau and delta values can be thought of as being given with the 
    // GERG reducing values, or tau0 = Tr_GERG/T = 0.xxx and delta0 = rho/rhor_GERG = 0.xxx 
    // Then we need to multiply by 
    // 
    // Tr_cubic/Tr_GERGlike & rhor_GERGlike/rhor_cubic 
    // 
    // to get shifted values:
    //
    // delta0 = rho/rhor_GERG*(rhor_GERGlike/rhor_cubic)
    // tau0 = Tr_GERG/T*(Tr_cubic/Tr_GERGlike)
    // 
    double Tr_GERGlike, rhor_GERGlike;
    get_linear_reducing_parameters(rhor_GERGlike, Tr_GERGlike);
    delta0 *= rhor_GERGlike/rhomolar_reducing();
    tau0 *= T_reducing()/Tr_GERGlike;
}

CoolPropDbl CoolProp::AbstractCubicBackend::calc_alphar_deriv_nocache(const int nTau, const int nDelta, const std::vector<CoolPropDbl> & mole_fractions, const CoolPropDbl &tau, const CoolPropDbl &delta){
    bool cache_values = true;
    HelmholtzDerivatives derivs = residual_helmholtz->all(*this, get_mole_fractions_ref(), cache_values);
    switch (nTau){
        case 0:
        {
            switch (nDelta){
                case 0: return derivs.alphar;
                case 1: return derivs.dalphar_ddelta;
                case 2: return derivs.d2alphar_ddelta2;
                case 3: return derivs.d3alphar_ddelta3;
                case 4: return derivs.d4alphar_ddelta4;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
            break;
        }
        case 1:
        {
            switch (nDelta){
                case 0: return derivs.dalphar_dtau;
                case 1: return derivs.d2alphar_ddelta_dtau;
                case 2: return derivs.d3alphar_ddelta2_dtau;
                case 3: return derivs.d4alphar_ddelta3_dtau;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
            break;
        }
        case 2:
        {
            switch (nDelta){
                case 0: return derivs.d2alphar_dtau2;
                case 1: return derivs.d3alphar_ddelta_dtau2;
                case 2: return derivs.d4alphar_ddelta2_dtau2;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        case 3:
        {
            switch (nDelta){
                case 0: return derivs.d3alphar_dtau3;
                case 1: return derivs.d4alphar_ddelta_dtau3;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        case 4:
        {
            switch (nDelta){
                case 0: return derivs.d4alphar_dtau4;
                default: throw ValueError(format("nDelta (%d) is invalid",nDelta));
            }
        }
        default: throw ValueError(format("nTau (%d) is invalid",nTau));
    }
}