
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

    // Set the pure fluid flag
    is_pure_or_pseudopure = (N == 1);
    
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
    
    // Set the alpha function for the backend
    set_alpha_from_components();
    
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

    // Resize the vectors (including linked states)
    resize(names.size());
}

void CoolProp::VTPRBackend::set_alpha_from_components(){
    
    VTPRCubic * _cubic= static_cast<VTPRCubic *>(cubic.get());
    const std::vector<UNIFAQLibrary::Component> &components = _cubic->get_unifaq().get_components();
    
    /// If components is not present, you are using a vanilla cubic, so don't do anything
    if (components.empty()){ return; }
    
    for (std::size_t i = 0; i < N; ++i){
        const std::string &alpha_type = components[i].alpha_type;
        if (alpha_type != "default"){
            const std::vector<double> &c = components[i].alpha_coeffs;
            shared_ptr<AbstractCubicAlphaFunction> acaf;
            if (alpha_type == "Twu"){
                acaf.reset(new TwuAlphaFunction(get_cubic()->a0_ii(i), c[0], c[1], c[2], get_cubic()->T_r/get_cubic()->get_Tc()[i]));
            }
            else if (alpha_type == "MathiasCopeman" || alpha_type == "Mathias-Copeman"){
                acaf.reset(new MathiasCopemanAlphaFunction(get_cubic()->a0_ii(i), c[0], c[1], c[2], get_cubic()->T_r / get_cubic()->get_Tc()[i]));
            }
            else{
                throw ValueError("alpha function is not understood");
            }
            cubic->set_alpha_function(i, acaf);
        }
    }
}

CoolPropDbl CoolProp::VTPRBackend::calc_molar_mass(void)
{
    double summer = 0;
    for (unsigned int i = 0; i < N; ++i){
        summer += mole_fractions[i]*molemass[i];
    }
    return summer;
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
        if (!XN_INDEPENDENT) {
            zp[2] -= dz; zm[2] += dz;
        }
        
        double dam_dxi_num = (cubic->am_term(tau, zp, 0) - cubic->am_term(tau, zm, 0))/(2*dz);
        double dam_dxi_ana = cubic->d_am_term_dxi(tau, z, 0, 0, XN_INDEPENDENT);
        double diff = dam_dxi_num-dam_dxi_ana;
        CHECK(std::abs(diff)<1e-6);
    }
}

#endif
