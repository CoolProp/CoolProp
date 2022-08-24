
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

#include "VTPRBackend.h"
#include "Configuration.h"
#include "Exceptions.h"

static UNIFACLibrary::UNIFACParameterLibrary lib;

void CoolProp::VTPRBackend::setup(const std::vector<std::string>& names, bool generate_SatL_and_SatV) {

    R = get_config_double(R_U_CODATA);

    // Set the pure fluid flag
    is_pure_or_pseudopure = (N == 1);

    // Reset the residual Helmholtz energy class
    residual_helmholtz.reset(new CubicResidualHelmholtz(this));

    // If pure, set the mole fractions to be unity
    if (is_pure_or_pseudopure) {
        mole_fractions = std::vector<CoolPropDbl>(1, 1.0);
        mole_fractions_double = std::vector<double>(1, 1.0);
    }

    // Now set the reducing function for the mixture
    Reducing.reset(new ConstantReducingFunction(cubic->get_Tr(), cubic->get_rhor()));

    VTPRCubic* _cubic = static_cast<VTPRCubic*>(cubic.get());
    _cubic->get_unifaq().set_components("name", names);
    _cubic->get_unifaq().set_interaction_parameters();

    // Store the fluid names
    m_fluid_names = names;

    // Set the alpha function for the backend
    set_alpha_from_components();

    // Set the ideal-gas helmholtz energy based on the components in use;
    set_alpha0_from_components();

    // Top-level class can hold copies of the base saturation classes,
    // saturation classes cannot hold copies of the saturation classes
    if (generate_SatL_and_SatV) {
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

void CoolProp::VTPRBackend::set_alpha_from_components() {

    VTPRCubic* _cubic = static_cast<VTPRCubic*>(cubic.get());
    const std::vector<UNIFACLibrary::Component>& components = _cubic->get_unifaq().get_components();

    /// If components is not present, you are using a vanilla cubic, so don't do anything
    if (components.empty()) {
        return;
    }

    for (std::size_t i = 0; i < N; ++i) {
        const std::string& alpha_type = components[i].alpha_type;
        if (alpha_type != "default") {
            const std::vector<double>& c = components[i].alpha_coeffs;
            shared_ptr<AbstractCubicAlphaFunction> acaf;
            if (alpha_type == "Twu") {
                acaf.reset(new TwuAlphaFunction(get_cubic()->a0_ii(i), c[0], c[1], c[2], get_cubic()->get_Tr() / get_cubic()->get_Tc()[i]));
            } else if (alpha_type == "MathiasCopeman" || alpha_type == "Mathias-Copeman") {
                acaf.reset(
                  new MathiasCopemanAlphaFunction(get_cubic()->a0_ii(i), c[0], c[1], c[2], get_cubic()->get_Tr() / get_cubic()->get_Tc()[i]));
            } else {
                throw ValueError("alpha function is not understood");
            }
            cubic->set_alpha_function(i, acaf);
        }
    }
}

CoolPropDbl CoolProp::VTPRBackend::calc_molar_mass(void) {
    double summer = 0;
    for (unsigned int i = 0; i < N; ++i) {
        summer += mole_fractions[i] * molemass[i];
    }
    return summer;
}

void CoolProp::VTPRBackend::set_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter,
                                                          const double value) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }    
    cubic->set_interaction_parameter(i, j, parameter, value);
    for (std::vector<shared_ptr<HelmholtzEOSMixtureBackend>>::iterator it = linked_states.begin(); it != linked_states.end(); ++it) {
        (*it)->set_binary_interaction_double(i, j, parameter, value);
    }
};

void CoolProp::VTPRBackend::set_Q_k(const size_t sgi, const double value) {
    cubic->set_Q_k(sgi, value);
};

double CoolProp::VTPRBackend::get_binary_interaction_double(const std::size_t i, const std::size_t j, const std::string& parameter) {
    // bound-check indices
    if (i < 0 || i >= N) {
        if (j < 0 || j >= N) {
            throw ValueError(format("Both indices i [%d] and j [%d] are out of bounds. Must be between 0 and %d.", i, j, N-1));
        } else {
            throw ValueError(format("Index i [%d] is out of bounds. Must be between 0 and %d.", i, N-1));
        }
    } else if (j < 0 || j >= N) {
        throw ValueError(format("Index j [%d] is out of bounds. Must be between 0 and %d.", j, N-1));
    }    
    return cubic->get_interaction_parameter(i, j, parameter);
};

const UNIFACLibrary::UNIFACParameterLibrary& CoolProp::VTPRBackend::LoadLibrary() {
    if (!lib.is_populated() || get_config_bool(VTPR_ALWAYS_RELOAD_LIBRARY)) {
        std::string UNIFAC_path = get_config_string(VTPR_UNIFAC_PATH);
        if (UNIFAC_path.empty()) {
            throw ValueError("You must provide the path to the UNIFAC library files as VTPR_UNIFAC_PATH");
        }
        if (!(UNIFAC_path[UNIFAC_path.size() - 1] == '\\' || UNIFAC_path[UNIFAC_path.size() - 1] == '/')) {
            throw ValueError("VTPR_UNIFAC_PATH must end with / or \\ character");
        }
        std::string group_path = UNIFAC_path + "group_data.json";
        std::string groups = get_file_contents(group_path.c_str());
        std::string interaction_path = UNIFAC_path + "interaction_parameters.json";
        std::string interaction = get_file_contents(interaction_path.c_str());
        std::string decomps_path = UNIFAC_path + "decompositions.json";
        std::string decomps = get_file_contents(decomps_path.c_str());
        lib.populate(groups, interaction, decomps);
    }
    return lib;
}

CoolPropDbl CoolProp::VTPRBackend::calc_fugacity_coefficient(std::size_t i) {
    //double slower = log(HelmholtzEOSMixtureBackend::calc_fugacity_coefficient(i));
    VTPRCubic* _cubic = static_cast<VTPRCubic*>(cubic.get());
    std::vector<double> here = _cubic->ln_fugacity_coefficient(mole_fractions, rhomolar(), p(), T());
    return exp(here[i]);
}

#ifdef ENABLE_CATCH
#    include <catch2/catch_all.hpp>

#    include "Backends/Cubics/CubicBackend.h"

using namespace CoolProp;

TEST_CASE("VTPR test", "[VTPR]") {
    shared_ptr<VTPRBackend> VTPR(new VTPRBackend(strsplit("Ethane&n-Propane&n-Butane", '&')));
    std::vector<double> z(3);
    z[0] = 0.1;
    z[1] = 0.2;
    z[2] = 0.7;
    VTPR->set_mole_fractions(z);

    SECTION("dam_dxi") {
        shared_ptr<AbstractCubic> cubic = VTPR->get_cubic();
        double tau = 0.001, dz = 1e-6;
        std::vector<double> zp = z, zm = z;
        zp[0] += dz;
        zm[0] -= dz;
        if (!XN_INDEPENDENT) {
            zp[2] -= dz;
            zm[2] += dz;
        }

        double dam_dxi_num = (cubic->am_term(tau, zp, 0) - cubic->am_term(tau, zm, 0)) / (2 * dz);
        double dam_dxi_ana = cubic->d_am_term_dxi(tau, z, 0, 0, XN_INDEPENDENT);
        double diff = dam_dxi_num - dam_dxi_ana;
        CHECK(std::abs(diff) < 1e-6);
    }
}

#endif
