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

namespace CoolProp {

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend  {
public:
    HelmholtzEOSBackend();
    HelmholtzEOSBackend(CoolPropFluid *pFluid){set_components(std::vector<CoolPropFluid*>(1,pFluid));};
    HelmholtzEOSBackend(const std::string &name){set_components(std::vector<CoolPropFluid*>(1,&(get_library().get(name))));};
    virtual ~HelmholtzEOSBackend(){};
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
