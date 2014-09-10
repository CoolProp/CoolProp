#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp{

class PhaseEnvelopeRoutines{
public:
    static void build(HelmholtzEOSMixtureBackend &HEOS);
    
    /** \brief Finalize the phase envelope and calculate maxima values, critical point, etc.
     */
    static void finalize(HelmholtzEOSMixtureBackend &HEOS);
};
    
} /* namespace CoolProp */