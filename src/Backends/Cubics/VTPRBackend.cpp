
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

#include "VTPRBackend.h"
#include "Configuration.h"
#include "Exceptions.h"

static UNIFAQLibrary::UNIFAQParameterLibrary lib;

void CoolProp::VTPRBackend::setup(const std::vector<std::string> &names, bool generate_SatL_and_SatV){

    R = get_config_double(R_U_CODATA);
    
    // Reset the residual Helmholtz energy class
    residual_helmholtz.reset(new CubicResidualHelmholtz(this));
    
    // If pure, set the mole fractions to be unity
    if (is_pure_or_pseudopure){
        mole_fractions = std::vector<CoolPropDbl>(1, 1.0);
        mole_fractions_double = std::vector<double>(1, 1.0);
    }
    
    // Now set the reducing function for the mixture
    Reducing.reset(new ConstantReducingFunction(cubic->T_r, cubic->rho_r));

    VTPRCubic * _cubic= static_cast<VTPRCubic *>(cubic.get());
    _cubic->get_unifaq().set_components("name", names);
    _cubic->get_unifaq().set_interaction_parameters();

    // Store the fluid names
    m_fluid_names = names;
    
    // Resize the vectors
    resize(names.size());
    
    // Top-level class can hold copies of the base saturation classes,
    // saturation classes cannot hold copies of the saturation classes
    if (generate_SatL_and_SatV)
    {
        bool SatLSatV = false;
        SatL.reset(this->get_copy(SatLSatV));
        SatL->specify_phase(iphase_liquid);
        linked_states.push_back(SatL);
        SatV.reset(this->get_copy(SatLSatV));
        SatV->specify_phase(iphase_gas);
        linked_states.push_back(SatV);

		if (is_pure_or_pseudopure) {
			std::vector<CoolPropDbl> z(1, 1.0);
			set_mole_fractions(z);
			SatL->set_mole_fractions(z);
			SatV->set_mole_fractions(z);
		}
    }
}

const UNIFAQLibrary::UNIFAQParameterLibrary & CoolProp::VTPRBackend::LoadLibrary(){
    if (!lib.is_populated()){
        std::string UNIFAQ_path = get_config_string(VTPR_UNIFAQ_PATH);
        if (UNIFAQ_path.empty()){
            throw ValueError("You must provide the path to the UNIFAQ library files as VTPR_UNIFAQ_PATH");
        }
        if (!(UNIFAQ_path[UNIFAQ_path.size()-1] == '\\' || UNIFAQ_path[UNIFAQ_path.size()-1] == '/')){
            throw ValueError("VTPR_UNIFAQ_PATH must end with / or \\ character");
        }
        std::string group_path = UNIFAQ_path + "/group_data.json";
        std::string groups = get_file_contents(group_path.c_str());
        std::string interaction_path = UNIFAQ_path + "/interaction_parameters.json";
        std::string interaction = get_file_contents(interaction_path.c_str());
        std::string decomps_path = UNIFAQ_path + "/decompositions.json";
        std::string decomps = get_file_contents(decomps_path.c_str());
        lib.populate(groups, interaction, decomps);
    }
    return lib;
}

#ifdef ENABLE_CATCH
#include "catch.hpp"

#include "Backends/Cubics/CubicBackend.h"

using namespace CoolProp;

TEST_CASE("VTPR test","[VTPR]")
{
    shared_ptr<VTPRBackend> VTPR(new VTPRBackend(strsplit("Ethane&n-Propane&n-Butane",'&')));
    std::vector<double> z(3); z[0] = 0.1; z[1] = 0.2; z[2] = 0.7;
    VTPR->set_mole_fractions(z);
    
    SECTION("dam_dxi"){
        shared_ptr<AbstractCubic> cubic = VTPR->get_cubic();
        double tau = 0.001, dz = 1e-6;
        std::vector<double> zp = z, zm = z;
        zp[0] += dz; zm[0] -= dz;
        
        double dam_dxi_num = (cubic->am_term(tau, zp, 0) - cubic->am_term(tau, zm, 0))/(2*dz);
        double dam_dxi_ana = cubic->d_am_term_dxi(tau, z, 0, 0, XN_INDEPENDENT);
        double diff = dam_dxi_num-dam_dxi_ana;
        CHECK(std::abs(diff)<1e-6);
    }
}

#endif
