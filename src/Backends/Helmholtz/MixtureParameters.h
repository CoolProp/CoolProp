#ifndef MIXTURE_PARAMETERS_H
#define MIXTURE_PARAMETERS_H

#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

class MixtureParameters
{
public:
    static void set_mixture_parameters(HelmholtzEOSMixtureBackend &HEOS);
};

} /* namespace CoolProp */
#endif