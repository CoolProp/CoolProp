#include "CoolProp.h"
#include "MixtureDerivatives.h"
#include <iostream>
using namespace CoolProp;
int main()
{
    // Ethane/Propane mixture, 25/75 molar
    std::vector<std::string> components(2,"Ethane"); components[1] = "Propane";
    std::vector<long double> z(2,0.25); z[1] = 0.75;
    
    shared_ptr<HelmholtzEOSMixtureBackend> HEOS(new HelmholtzEOSMixtureBackend(components));
	HelmholtzEOSMixtureBackend &rHEOS = *(HEOS.get());
    HEOS->set_mole_fractions(z);
	HEOS->specify_phase(iphase_gas); // So that we don't do a 
    HEOS->update(DmolarT_INPUTS, 300, 300);
    
    std::vector<std::string> terms;
    terms.push_back("p");
	terms.push_back("p2(deriv)");
	terms.push_back("rhor");
	terms.push_back("Tr");
	terms.push_back("dalphar_dDelta");
	terms.push_back("dTr_dxi");
	terms.push_back("drhor_dxi");	
	terms.push_back("ndpdV__constT_n");
	terms.push_back("dpdxj__constT_V_xi");

	terms.push_back("ln_fugacity_coefficient");
    terms.push_back("ndpdni__constT_V_nj");
	terms.push_back("dln_fugacity_coefficient_dxj__constT_p_xi");
	terms.push_back("d2nalphar_dxj_dni__constT_V");
	
	terms.push_back("d2alphar_dxi_dDelta");
    
    for (std::vector<std::string>::iterator it = terms.begin(); it != terms.end(); ++it)
    {
        if (!it->compare("p")){
            printf("p: %0.16Lg\n", HEOS->p());
        }
		else if (!it->compare("p2(deriv)")){
			printf("p calculated by rho*R*T*(1+deltadar_dDelta): %0.16Lg\n", HEOS->rhomolar()*HEOS->gas_constant()*HEOS->T()*(1+HEOS->delta()*HEOS->dalphar_dDelta()));
        }
		else if (!it->compare("dalphar_dDelta")){
            printf("dalphar_dDelta: %0.16Lg\n", HEOS->dalphar_dDelta());
        }
		else if (!it->compare("rhor")){
			printf("rhor: %0.16Lg\n", HEOS->get_reducing_state().rhomolar);
        }
		else if (!it->compare("Tr")){
			printf("Tr: %0.16Lg\n", HEOS->get_reducing_state().T);
        }
		else if (!it->compare("dTr_dxi")){
			printf("dTr_dxi: %0.16Lg\n", HEOS->Reducing.p->dTrdxi__constxj(rHEOS.get_mole_fractions(), 0, XN_DEPENDENT));
        }
		else if (!it->compare("drhor_dxi")){
			printf("drhor_dxi: %0.16Lg\n", HEOS->Reducing.p->drhormolardxi__constxj(rHEOS.get_mole_fractions(), 0, XN_DEPENDENT));
        }
		else if(!it->compare("ndpdV__constT_n")){
            printf("ndpdV__constT_n: %0.16Lg\n", MixtureDerivatives::ndpdV__constT_n(rHEOS));
        }
		else if(!it->compare("ln_fugacity_coefficient")){
			printf("ln_fugacity_coefficient(0): %0.16Lg\n", MixtureDerivatives::ln_fugacity_coefficient(rHEOS, 0, XN_DEPENDENT));
            printf("ln_fugacity_coefficient(1): %0.16Lg\n", MixtureDerivatives::ln_fugacity_coefficient(rHEOS, 1, XN_DEPENDENT));
        }
        else if(!it->compare("dln_fugacity_coefficient_dxj__constT_p_xi")){
            printf("dln_fugacity_coefficient_dxj__constT_p_xi(0,0): %0.16Lg\n", MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(rHEOS, 0, 0, XN_DEPENDENT));
            //printf("dln_fugacity_coefficient_dxj__constT_p_xi(0,1): %0.16Lg\n", MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(rHEOS, 0, 1, XN_DEPENDENT));
			printf("dln_fugacity_coefficient_dxj__constT_p_xi(1,0): %0.16Lg\n", MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(rHEOS, 1, 0, XN_DEPENDENT));
            //printf("dln_fugacity_coefficient_dxj__constT_p_xi(1,1): %0.16Lg\n", MixtureDerivatives::dln_fugacity_coefficient_dxj__constT_p_xi(rHEOS, 1, 1, XN_DEPENDENT));
        }
		else if(!it->compare("d2nalphar_dxj_dni__constT_V")){
            printf("d2nalphar_dxj_dni__constT_V(0,0): %0.16Lg\n", MixtureDerivatives::d2nalphar_dxj_dni__constT_V(rHEOS, 0, 0, XN_DEPENDENT));
            //printf("d2nalphar_dxj_dni__constT_V(0,1): %0.16Lg\n", MixtureDerivatives::d2nalphar_dxj_dni__constT_V(rHEOS, 0, 1, XN_DEPENDENT));
			printf("d2nalphar_dxj_dni__constT_V(1,0): %0.16Lg\n", MixtureDerivatives::d2nalphar_dxj_dni__constT_V(rHEOS, 1, 0, XN_DEPENDENT));
			//printf("d2nalphar_dxj_dni__constT_V(1,1): %0.16Lg\n", MixtureDerivatives::d2nalphar_dxj_dni__constT_V(rHEOS, 1, 1, XN_DEPENDENT));
        }
		else if(!it->compare("d2alphar_dxi_dDelta")){
			printf("d2alphar_dxi_dDelta(0): %0.16Lg\n", MixtureDerivatives::d2alphar_dxi_dDelta(rHEOS, 0, XN_DEPENDENT));
            //printf("d2alphar_dxi_dDelta(1): %0.16Lg\n", MixtureDerivatives::d2alphar_dxi_dDelta(rHEOS, 1, XN_DEPENDENT));
        }
		else if(!it->compare("ndpdni__constT_V_nj")){
            printf("ndpdni__constT_V_nj(0): %0.16Lg\n", MixtureDerivatives::ndpdni__constT_V_nj(rHEOS, 0, XN_DEPENDENT));
            //printf("ndpdni__constT_V_nj(1): %0.16Lg\n", MixtureDerivatives::ndpdni__constT_V_nj(rHEOS, 1, XN_DEPENDENT));
        }
		else if(!it->compare("dpdxj__constT_V_xi")){
            printf("dpdxj__constT_V_xi(0): %0.16Lg\n", MixtureDerivatives::dpdxj__constT_V_xi(rHEOS, 0, XN_DEPENDENT));
            //printf("dpdxj__constT_V_xi(1): %0.16Lg\n", MixtureDerivatives::dpdxj__constT_V_xi(rHEOS, 1, XN_DEPENDENT));
        }
    }    
    
    return EXIT_SUCCESS;
}