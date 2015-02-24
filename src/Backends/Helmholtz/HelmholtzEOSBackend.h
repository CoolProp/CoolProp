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

namespace CoolProp {

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend  {
public:
    HelmholtzEOSBackend();
    HelmholtzEOSBackend(CoolPropFluid *pFluid){set_components(std::vector<CoolPropFluid*>(1,pFluid));};
    HelmholtzEOSBackend(const std::string &name){
        Dictionary dict;
        std::vector<double> mole_fractions;
        std::vector<CoolPropFluid*> components;
        if (is_predefined_mixture(name, dict)){
            std::vector<std::string> fluids = dict.get_string_vector("fluids");
            mole_fractions = dict.get_double_vector("mole_fractions");
            
            components.resize(fluids.size());
            for (unsigned int i = 0; i < components.size(); ++i){
                components[i] = &(get_library().get(fluids[i]));
            }
        }
        else{
            components.push_back(&(get_library().get(name))); // Until now it's empty
            mole_fractions.push_back(1.);
        }
        // Set the components
        set_components(components);
        // Set the mole fractions
        set_mole_fractions(std::vector<CoolPropDbl>(mole_fractions.begin(), mole_fractions.end()));
    };
    virtual ~HelmholtzEOSBackend(){};
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
