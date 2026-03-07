/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSBACKEND_H_
#define HELMHOLTZEOSBACKEND_H_

#include <vector>
#include "HelmholtzEOSMixtureBackend.h"
#include "Fluids/FluidLibrary.h"
#include "MixtureParameters.h"
#include "DataStructures.h"

namespace CoolProp {

///Templates for turning vectors (1D-matrices) into strings
inline std::string vecstring_to_string(const std::vector<std::string>& a) {
    std::stringstream out;
    out << "[ " << format("%s", a[0].c_str());
    for (size_t j = 1; j < a.size(); j++) {
        out << ", " << format("%s", a[j].c_str());
    }
    out << " ]";
    return out.str();
};

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend
{
   public:
    HelmholtzEOSBackend() {};
    HelmholtzEOSBackend(CoolPropFluid Fluid) {
        set_components(std::vector<CoolPropFluid>(1, Fluid));
    };
    HelmholtzEOSBackend(const std::string& name) : HelmholtzEOSMixtureBackend() {
        Dictionary dict;
        std::vector<double> mole_fractions;
        std::vector<CoolPropFluid> components;
        CoolProp::JSONFluidLibrary& library = get_library();
        // Strip optional JSON config suffix (e.g. Water{"EOS":"Wagner-JPCRD-2002"})
        auto brace_pos = name.find('{');
        const std::string clean_name = (brace_pos == std::string::npos) ? name : name.substr(0, brace_pos);
        const std::string json_str = (brace_pos == std::string::npos) ? std::string() : name.substr(brace_pos);
        if (is_predefined_mixture(clean_name, dict)) {
            std::vector<std::string> fluids = dict.get_string_vector("fluids");
            mole_fractions = dict.get_double_vector("mole_fractions");
            if (get_debug_level() > 0) {
                std::cout << "Got the fluids" << vecstring_to_string(fluids) << std::endl;
                std::cout << "Got the fractions" << vec_to_string(mole_fractions, "%g") << std::endl;
            }
            if (!json_str.empty()) {
                throw ValueError("JSON EOS config is not supported for predefined mixtures like \"" + clean_name + "\"");
            }
            for (unsigned int i = 0; i < fluids.size(); ++i) {
                components.push_back(library.get(fluids[i]));
            }
        } else {
            components.push_back(library.get(clean_name));
            // Pre-set selected_EOS_index before set_components so the reducing
            // function is built from the correct EOS on the first construction pass.
            apply_json_to_fluid(components[0], json_str);
            mole_fractions.push_back(1.);
        }
        // Set the components
        set_components(components);
        // Set the mole fractions
        set_mole_fractions(std::vector<CoolPropDbl>(mole_fractions.begin(), mole_fractions.end()));
        if (get_debug_level() > 0) {
            std::cout << "successfully set up state" << std::endl;
        }
    };
    virtual ~HelmholtzEOSBackend() {};
    std::string backend_name(void) {
        return get_backend_string(HEOS_BACKEND_PURE);
    }
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
