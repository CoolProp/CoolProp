/*
 * AbstractBackend.cpp
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#include "REFPROPBackend.h"

#include "CoolPropTools.h"

#if defined(_MSC_VER)
#    define _CRTDBG_MAP_ALLOC
#    ifndef _CRT_SECURE_NO_WARNINGS
#        define _CRT_SECURE_NO_WARNINGS
#    endif
#    include <crtdbg.h>
#    include <sys/stat.h>
#else
#    include <sys/stat.h>
#endif

#include <string>

namespace CoolProp {

REFPROPBackend::REFPROPBackend(const std::string& fluid_name) {

    // Do the REFPROP instantiation for this fluid
    std::vector<std::string> component_names(1, fluid_name);
    construct(component_names);

    // Set the mole fraction to 1 in the base class (we can't set the mole fraction in this class,
    // otherwise a NotImplementedError will be returned)
    if (get_mole_fractions().empty()) {
        std::vector<CoolPropDbl> x(1, 1.0);  // (one element with value of 1.0)
        REFPROPMixtureBackend::set_mole_fractions(x);
    }
}

REFPROPBackend::~REFPROPBackend() {
    // TODO Auto-generated destructor stub
}

} /* namespace CoolProp */
