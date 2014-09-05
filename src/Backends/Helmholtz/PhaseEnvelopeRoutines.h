#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

class PhaseEnvelopeRoutines{
public:
    static void build(HelmholtzEOSMixtureBackend &HEOS, const std::vector<long double> &z, SaturationSolvers::newton_raphson_saturation_options &IO);
};
    
} /* namespace CoolProp */